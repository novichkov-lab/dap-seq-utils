import csv
from DapSeqAgent.ReferenceData.Genome import Genome
from DapSeqAgent.ReferenceData.Gene import Gene
from DapSeqAgent.ReferenceData.Contig import Contig
from DapSeqAgent.ReferenceData.Operon import Operon
from DapSeqAgent.ReferenceData.Utr import Utr
from DapSeqAgent.ReferenceData.filepaths import file_paths

class ReferenceLibrary:
    def __init__(self, genome_id):
        self.genome_id = genome_id
        self.genome = None
        self.bowtie_path = ''
        self.fitgenomics_path = ''
        self.operons_path = ''
        self.replicons_path = ''
        self.gene_utr_path = ''
        self.operon_utr_path = ''
        self.load_library()
        
    def load_library(self):
        self.bowtie_path = file_paths[self.genome_id]['bowtie_lib']
        self.fitgenomics_path = file_paths[self.genome_id]['fitgenomics_file']
        self.operons_path = file_paths[self.genome_id]['operons_file']
        self.replicons_path = file_paths[self.genome_id]['replicons_file']
        self.sequence_path = file_paths[self.genome_id]['sequence_file']
        self.utrs_path = file_paths[self.genome_id]['utrs_file']
        self.genome = self.load_genome()
        self.genome.size = file_paths[self.genome_id]['size']
        
    def load_genome(self):
        genome = Genome(self.genome_id, name=file_paths[self.genome_id]['name'])

        with open(self.replicons_path, 'r') as f:
            tsvin = csv.reader(f, delimiter='\t')
            for row in tsvin:
                contig = Contig(row[0], size=int(row[1]))
                genome.contigs[row[0]] = contig
        with open(self.fitgenomics_path, 'r') as f:
            tsvin = csv.reader(f, delimiter='\t')
            for row in tsvin:
                if row[0] == 'locusId':
                    continue
                gene = Gene(row[0], self.genome_id, locus_tag=row[1], \
                        sequence_id=row[2], start=int(row[3]), end=int(row[4]), strand=row[5])
                gene.description=row[7]
                genome.contigs[row[2]].add_gene(gene)
        
        with open(self.operons_path, 'r') as f:
            tsvin = csv.reader(f, delimiter='\t')
            for row in tsvin:
                if row[0] == 'Operon number':
                    continue
                operon = Operon(row[0], self.genome_id, row[1], \
                        start=int(row[2]), end=int(row[3]), strand=row[4])
                operon.gene_count = int(row[5])
                if row[4] == '+':
                    for gene in row[6:]:
                        if gene != '':
                            operon.gene_ids.append(gene)
                elif row[4] == '-':
                    for gene in list(reversed(row[6:])):
                        if gene != '':
                            operon.gene_ids.append(gene)
                genome.contigs[row[1]].add_operon(operon)

        with open(self.sequence_path, 'r') as f:
            current_id = ''
            sequence_lines = []
            for line in f:
                line = line.rstrip('\n\r')
                if line.startswith ('>'):
                    if current_id != '':
                        if current_id in genome.contigs:
                            genome.contigs[current_id].sequence = ''.join(sequence_lines)
                        else:
                            raise ValueError('Sequence import error: contig ' + current_id + ' not found')
                    current_id = line[1:].split(' ')[0]
                    if current_id.endswith('|'):
                        current_id = current_id[:-1]
                    sequence_lines = []
                else:
                    sequence_lines.append(line)
            if current_id in genome.contigs:
                genome.contigs[current_id].sequence = ''.join(sequence_lines)
            else:
                raise ValueError('Sequence import error: contig ' + current_id + ' not found')

        with open(self.utrs_path, 'r') as f:
            tsvin = csv.reader(f, delimiter='\t')
            for row in tsvin:
                if row[0] == '#Locus':
                    continue
                utr = Utr(utr_id = row[0],
                        genome_id = self.genome_id,
                        sequence_id = row[1],
                        start = int(row[2]),
                        end = int(row[3])
                        )
                if row[5] != '':
                    utr.rev_gene_id = row[5]
                if row[6] != '':
                    utr.fwd_gene_id = row[6]
                if row[7] != '':
                    utr.rev_operon_id = row[7]
                if row[8] != '':
                    utr.fwd_operon_id = row[7]
                genome.contigs[row[1]].add_utr(utr)
            
        for contig_id in genome.contigs.keys():
            #genome.contigs[contig_id].calculate_utrs()
            print('Contig',contig_id,':', len(genome.contigs[contig_id].genes),'genes found')
            print('Contig',contig_id,':', len(genome.contigs[contig_id].operons),'operons found')
            print('Contig',contig_id,':', len(genome.contigs[contig_id].utrs),'UTRs found')
        return genome
    
    def get_genes(self, sequence_id, start, end):
        ret_val = []
        for c in self.genome.contigs.keys():
            contig = self.genome.contigs[c]
            if contig.contig_id != sequence_id:
                continue
            for gene in contig.genes:
                if gene.start >= start and gene.start <= end:
                    ret_val.append(gene.gene_id)
                elif gene.end >= start and gene.end <= end:
                    ret_val.append(gene.gene_id)
                elif gene.start < start and gene.end > start:
                    ret_val.append(gene.gene_id)
        return ret_val

    def get_operons(self, sequence_id, start, end):
        ret_val = []
        for c in self.genome.contigs.keys():
            contig = self.genome.contigs[c]
            if contig.contig_id != sequence_id:
                continue
            for operon in contig.operons:
                if operon.start >= start and operon.start <= end:
                    ret_val.append(str(operon))
                elif operon.end >= start and operon.end <= end:
                    ret_val.append(str(operon))
                elif operon.start < start and operon.end > start:
                    ret_val.append(str(operon))
        return ret_val

    def get_utrs(self, sequence_id, start, end):
        ret_val = []
        for c in self.genome.contigs.keys():
            contig = self.genome.contigs[c]
            if contig.contig_id != sequence_id:
                continue
            for utr in contig.utrs:
                if utr.start >= start and utr.start <= end:
                    ret_val.append(utr.utr_id)
                elif utr.end >= start and utr.end <= end:
                    ret_val.append(utr.utr_id)
                elif utr.start < start and utr.end > start:
                    ret_val.append(utr.utr_id)
        return ret_val

    def get_fuction_by_gene_id(self, sequence_id, gene_id):
        ret_val = []
        for c in self.genome.contigs.keys():
            contig = self.genome.contigs[c]
            if contig.contig_id != sequence_id:
                continue
            for gene in contig.genes:
                if gene.gene_id == gene_id:
                    return gene.description
        return ''

    def get_fuction_by_locus_tag(self, sequence_id, locus_tag):
        ret_val = []
        for c in self.genome.contigs.keys():
            contig = self.genome.contigs[c]
            if contig.contig_id != sequence_id:
                continue
            for gene in contig.genes:
                if gene.locus_tag == locus_tag:
                    return gene.description
        return ''

    def get_utr_by_id(self, utr_id):
        for c in self.genome.contigs.keys():
            contig = self.genome.contigs[c]
            for utr in contig.utrs:
                if utr.utr_id == utr_id:
                    return utr
        return None

    def get_operon_id_by_first_gene(self,gene_id):
        for c in self.genome.contigs.keys():
            contig = self.genome.contigs[c]
            for operon in contig.operons:
                if operon.gene_ids[0] == gene_id:
                    return operon.operon_id
        return None
        
    #~ def is_noncoding(self, sequence_id, position):
        #~ ret_val = False
        #~ for c in self.genome.contigs.keys():
            #~ contig = self.genome.contigs[c]
            #~ if contig.contig_id != sequence_id:
                #~ continue
            #~ for utr in contig.utrs:
                #~ if utr.start <= position and utr.end >= position:
                    #~ ret_val = True
                    #~ break
        #~ return ret_val

    def get_sequence(self, sequence_id, start, end):
        end = end - 1 
        return self.genome.contigs[sequence_id].sequence[start:end]
