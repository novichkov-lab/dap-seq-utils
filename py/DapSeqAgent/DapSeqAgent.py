#!/usr/bin/python3
import os
import gzip
from subprocess import Popen, PIPE, CalledProcessError
from DapSeqAgent.ReferenceData.ReferenceLibrary import ReferenceLibrary
from DapSeqAgent.Task import Task
from DapSeqAgent.MacsParser.MacsParser import MacsParser
from DapSeqAgent.MemeParser.MemeParser import MemeParser
from DapSeqAgent.JSONUtil import export_task,import_task
from DapSeqAgent.Report import generate_text_report, generate_mapping_report

class DapSeqAgent:
    
    def __init__(self, filename, indir, outdir, sample=None):
        self.task = Task(filename)
        self.indir = indir
        self.outdir = outdir
        self.ref_libraries = self.load_libraries()
    
    def load_libraries(self):
        ret_val = {}
        for sample in self.task.samples.keys():
            ret_val[self.task.samples[sample]['ref_library']] = None
        for lib_id in ret_val:
            ret_val[lib_id] = ReferenceLibrary(lib_id)
        return ret_val
        
    def run_task(self):
        # Make directories
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)
        if not os.path.isdir(os.path.join(self.outdir, 'fastq')):
            os.mkdir(os.path.join(self.outdir, 'fastq'))
        if not os.path.isdir(os.path.join(self.outdir, 'focus')):
            os.mkdir(os.path.join(self.outdir, 'focus'))
        if not os.path.isdir(os.path.join(self.outdir, 'bam')):
            os.mkdir(os.path.join(self.outdir, 'bam'))
        if not os.path.isdir(os.path.join(self.outdir, 'peaks')):
            os.mkdir(os.path.join(self.outdir, 'peaks'))
        if not os.path.isdir(os.path.join(self.outdir, 'motifs')):
            os.mkdir(os.path.join(self.outdir, 'motifs'))
        # Try to load task JSON
        if os.path.exists(os.path.join(self.outdir, 'task.json')):
            self.task = import_task(os.path.join(self.outdir, 'task.json'))
            self.task.load_task()
        
        # Run QC and mapping
        for sample in self.task.samples.keys():
            filtered_fastq, trimmomatic_log = run_trimmomatic_se(self.task.samples[sample]['fastq'], self.outdir)
            if trimmomatic_log:
                self.task.samples[sample]['info']['trimmomatic_log'] = trimmomatic_log
            #trimmed_fastq = run_read_trimming(self.task.samples[sample]['fastq'], os.path.join(self.outdir, 'fastq'))
            #filtered_fastq = run_read_filtering(trimmed_fastq, os.path.join(self.outdir, 'focus'))
            if 'average_read_length' not in self.task.samples[sample]['info']:
                self.task.samples[sample]['info']['average_read_length'] = get_average_length(filtered_fastq)
            self.map_reads(sample, filtered_fastq)
            
        # Run peak calling
        for sample in self.task.samples.keys():
            if self.task.samples[sample]['reads aligned'] < 10:
                continue
            bam_file = os.path.join(self.outdir, 'bam', sample + '.sorted.bam')
            control_bam_file = os.path.join(self.outdir, 'bam', self.task.samples[sample]['control'] + '.sorted.bam')
            if bam_file == control_bam_file:
                continue
            outfile = os.path.join(self.outdir, 'peaks',sample + '_vs_' + self.task.samples[sample]['control'])
            if not os.path.exists(outfile + '_peaks.xls'):
                run_macs2(bam_file, control_bam_file, str(self.ref_libraries[self.task.samples[sample]['ref_library']].genome.size), outfile)
            macs2_outfile = outfile + '_peaks.xls'
            parser = MacsParser(macs2_outfile, protein=self.task.samples[sample]['protein'])
            self.task.samples[sample]['peaks'] = parser.peaks
            if len(self.task.samples[sample]['peaks']) > 1:
                self.task.samples[sample]['motifs'] = self.predict_motifs(sample, self.task.samples[sample]['peaks'])

        export_task(self.task, os.path.join(self.outdir, 'task.json'))
        self.generate_report()
    
    def generate_report(self):
        outfile = os.path.join(self.outdir, 'report.txt')
        for sample in self.task.samples.keys():
            if 'peaks' not in self.task.samples[sample]:
                continue
            for peak in self.task.samples[sample]['peaks']:
                peak.annotate_peak(self.ref_libraries[self.task.samples[sample]['ref_library']])
        generate_text_report(self.task, outfile)
        generate_mapping_report(self.task, os.path.join(self.outdir,'mapping_report.tsv'))
        
        
    def map_reads(self, sample, infile):
        if not os.path.exists(infile):
            return
        if os.path.exists(os.path.join(self.outdir, 'bam', sample + '.sorted.bam')):
            return
        if infile.endswith('.gz'):
            fastq = infile[:-3]
            with gzip.open(infile, 'rt') as f:
                with open(fastq, 'w') as of:
                    of.writelines(f)
        else:
            fastq = infile
        print('Starting bowtie', fastq)
        sam_file = os.path.join(self.outdir, 'bam', sample + '.sam')
        process_args = ['bowtie',
                        '-S',
                        self.ref_libraries[self.task.samples[sample]['ref_library']].bowtie_path,
                        '-m',
                        '1',
                        fastq,
                        sam_file
                        #~ ,
                        #~ '1>>',
                        #~ 'mapping_quality.txt',
                        #~ '2>&1'
                        ]
        bowtie_log = []
        with Popen(process_args, stdout=PIPE, stderr=PIPE, bufsize=1, universal_newlines=True) as p:
            with open (os.path.join(self.outdir, 'mapping_quality.txt'), 'a') as of:
                of.write(sample + ': Bowtie output\n')
                for line in p.stdout:
                    of.write(line)
                    bowtie_log.append(line)
                    print(line, end='') 
                for line in p.stderr:
                    of.write(line)
                    bowtie_log.append(line)
                    print(line, end='') 
                
        if p.returncode != 0:
            raise CalledProcessError(p.returncode, p.args)
        print ('Bowtie finished')
        for line in bowtie_log:
            if line.startswith('# reads processed'):
                k,v = line[2:].split(': ')
                self.task.samples[sample]['info'][k] = int(v)
            elif line.startswith('# reads with at least one reported alignment'):
                k,v = line[2:].split(': ')
                v1,v2 = v.split(' (')
                self.task.samples[sample]['info']['reads aligned'] = int(v1)
                

        if infile.endswith('.gz'):
            os.remove(fastq)
        bam_file = os.path.join(self.outdir, 'bam', sample + '.bam')
        process_args = ['samtools', 'view', '-bS', '-o',
                        bam_file,
                        sam_file
                        ]
        with Popen(process_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')
        if p.returncode != 0:
            raise CalledProcessError(p.returncode, p.args)

        os.remove(sam_file)
        
        process_args = ['samtools', 
                        'sort', 
                        '-o',
                        os.path.join(self.outdir, 'bam', sample + '.sorted.bam'),
                        bam_file
                        ]
        with Popen(process_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')
        if p.returncode != 0:
            raise CalledProcessError(p.returncode, p.args)
        
        os.remove(bam_file)
        print('Mapping finished:', sample)
        
    def predict_motifs(self, sample, peaks):
        temp_fastafile = os.path.join(self.outdir, 'temp.fasta')
        outfile = os.path.join(self.outdir, 'motifs', sample + '_meme.txt')
        # Generate fasta file with peak sequences
        peaks_dict = {}
        peak_index = 0
        with open(temp_fastafile, 'w') as of:
            for peak in peaks:
                peak_sequence = self.ref_libraries[self.task.samples[sample]['ref_library']].get_sequence(peak.sequence_id,peak.start, peak.end)
                peak_index += 1
                peaks_dict[str(peak_index)] = peak
                if len(peak_sequence) > 12:
                    of.write('>' + str(peak_index) + '\n')
                    of.write(peak_sequence + '\n')
        # Run MEME
        run_meme(temp_fastafile, outfile, self.outdir)
        # Parse MEME output
        meme_parser = MemeParser(outfile, peaks_dict)
        os.remove(temp_fastafile)
        return meme_parser.motifs
            
            
def run_trimmomatic_se(infile, outdir):
    print('Run TrimmomaticSE for ', infile)
    fastq_outdir = os.path.join(outdir, 'fastq')
    fastq_basename = os.path.basename(infile)
    if fastq_basename.endswith('.gz'):
        outfile = os.path.join(fastq_outdir, fastq_basename[:-3] + '.trimmed.gz')
    else:
        outfile = os.path.join(fastq_outdir, fastq_basename + '.trimmed')
    if os.path.exists(outfile):
#        run_focus(outfile, os.path.join(outdir, 'focus'))
        return outfile, None
    process_args = ['TrimmomaticSE',
                    '-threads',
                    '12',
                    '-phred33',
                    infile,
                    outfile,
                    'ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10',
                    'LEADING:3',
                    'TRAILING:3',
                    'SLIDINGWINDOW:4:14',
                    'MINLEN:50'
                    ]
    log = []
    with Popen(process_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
            log.append(line)
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)
    print ('TrimmomaticSE finished')
    run_focus(outfile, os.path.join(outdir, 'focus'))
    return outfile, log
    
    

def run_read_trimming(infile, outdir):
    if infile.endswith('.gz'):
        fastq = infile[:-3]
        with gzip.open(infile, 'rt') as f:
            with open(fastq, 'w') as of:
                of.writelines(f)
    else:
        fastq = infile
    print('Trimming reads', infile)
    process_args = ['SolexaQA++',
                    'dynamictrim',
                    fastq,
                    '-d',
                    outdir
                    ]
    with Popen(process_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)
    print ('SolexaQA++ finished')
    if infile.endswith('.gz'):
        os.remove(fastq)
    return os.path.join(outdir, os.path.basename(fastq) + '.trimmed')

def run_read_filtering(infile, focus_dir):
    min_length = "30"
    print('Filtering reads', infile)
    process_args = ['SolexaQA++',
                    'lengthsort',
                    infile,
                    '-l',
                    min_length
                    ]
    with Popen(process_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)
    print ('SolexaQA++ finished')
    os.remove(infile)
    with open(infile + '.single', 'rt') as f_in, gzip.open(infile + '.single.gz', 'wt') as f_out:
        f_out.writelines(f_in)
    with open(infile + '.discard', 'rt') as f_in, gzip.open(infile + '.discard.gz', 'wt') as f_out:
        f_out.writelines(f_in)
    os.remove(infile + '.discard')
    run_focus(infile + '.single', focus_dir)
    os.remove(infile + '.single')
    return infile + '.single.gz'

def run_focus(infile, outdir):
    print('Running FOCUS', infile)
    focus_infile = infile
    if infile.endswith('.gz'):
        focus_infile = infile[:-3]
        with gzip.open(infile, 'rt') as f_in, open(infile[:-3], 'w') as f_out:
            f_out.writelines(f_in)
    
    prefix = os.path.basename(infile)
    process_args = ['python3',
                    '/home/aekazakov/Soft/FOCUS/FOCUS/focus_app/focus_my_1file.py',
                    '-q',
                    focus_infile,
                    '-o',
                    outdir,
                    '-p',
                    prefix
                    ]
    with Popen(process_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)
    print ('FOCUS finished')
    if infile.endswith('.gz'):
        os.remove(focus_infile)
    
    
def run_macs2 (bam, control_bam, size, outfile):
    print('Starting MACS2')
    process_args = ['macs2',
                    'callpeak',
                    '-t',
                    bam,
                    '-c',
                    control_bam,
                    '-f',
                    'BAM',
                    '-g',
                    size,
                    '-n',
                    outfile,
                    '-q',
                    '0.0001',
                    '--extsize',
                    '250',
                    '--nomodel',
                    '--keep-dup=auto'
                    ]
    with Popen(process_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)
    print ('MACS2 finished')

def run_meme (temp_fastafile, outfile, output_dir):
    if os.path.exists(outfile):
        return
    print('Running MEME')
    ret_val = []
    process_args = ['meme',
                    temp_fastafile,
                    '-o',
                    output_dir,
                    '-text',
                    '-dna',
                    '-mod',
                    'anr',
                    '-nmotifs',
                    '1',
                    '-minw',
                    '12',
                    '-maxw',
                    '30',
                    '-revcomp',
                    '-pal',
                    '-maxsize',
                    '1000000'
                    ]
    with Popen(process_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            ret_val.append(line)
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)
    with open(outfile, 'w') as of:
        of.write(''.join(ret_val))
    print ('MEME finished')

def get_average_length(fastq):
    read_count = 0
    sequence_length = 0
    line_counter = 0
    if fastq.endswith('.gz'):
        fh = gzip.open(fastq, 'rt')
    else:
        fh = open(fastq, 'r')
    for line in fh:
        if line_counter == 1:
            read_count += 1
            sequence_length += len(line) - 1
        line_counter += 1
        if line_counter == 4:
            line_counter = 0
    fh.close()
    if read_count > 0:
        return sequence_length/read_count
    else:
        return 0.0
