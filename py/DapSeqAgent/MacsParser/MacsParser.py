import os, csv, re
from DapSeqAgent.utils import autovivify
from DapSeqAgent.MacsParser.Peak import Peak

class MacsParser:

    def __init__(self, macsfile, protein=None):
        self.filename = macsfile
        self.protein = protein
        self.peaks = []
        self.band_width = None
        self.tag_size = None
        self.tags_treatment_total = None
        self.tags_treatment_filtered = None
        self.tags_treatment_redundant_rate = None
        self.tags_control_total = None
        self.tags_control_filtered = None
        self.tags_control_redundant_rate = None
        
        self.parse_file()

    def parse_file(self):
        with open(self.filename, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')
                if line.startswith('#'):
                    if line.startswith('# band width'):
                        line_tokens = line.split(' = ')
                        self.band_width = int(line_tokens[-1])
                    elif line.startswith('# tag size is determined as'):
                        line_tokens = line.split(' ')
                        self.tag_size = int(line_tokens[-2])
                    elif line.startswith('# total tags in treatment'):
                        line_tokens = line.split(': ')
                        self.tags_treatment_total = int(line_tokens[-1])
                    elif line.startswith('# tags after filtering in treatment:'):
                        line_tokens = line.split(': ')
                        self.tags_treatment_filtered = int(line_tokens[-1])
                    elif line.startswith('# Redundant rate in treatment:'):
                        line_tokens = line.split(': ')
                        self.tags_treatment_redundant_rate = float(line_tokens[-1])
                    elif line.startswith('# total tags in control:'):
                        line_tokens = line.split(': ')
                        self.tags_control_total = int(line_tokens[-1])
                    elif line.startswith('# tags after filtering in control:'):
                        line_tokens = line.split(': ')
                        self.tags_control_filtered = int(line_tokens[-1])
                    elif line.startswith('# Redundant rate in control:'):
                        line_tokens = line.split(': ')
                        self.tags_control_redundant_rate = float(line_tokens[-1])
                        
                elif line.startswith('chr'):
                    pass
                elif line == '':
                    pass
                else:
                    line_tokens = line.split('\t')
                    seq_id = line_tokens[0].split('|')[0]
                    peak = Peak(line_tokens[9].split('/')[-1], \
                                    seq_id, \
                                    int(line_tokens[1]), \
                                    int(line_tokens[2]), \
                                    abs_summit = int(line_tokens[4]), \
                                    length = int(line_tokens[3]), \
                                    pileup = float(line_tokens[5]), \
                                    log_pvalue = float(line_tokens[6]), \
                                    log_qvalue = float(line_tokens[8]), \
                                    fold_enrichment = float(line_tokens[7]),
                                    protein_id = self.protein
                                    )
                    self.peaks.append(peak)

    def filter_peaks(self, pvalue_cutoff=None, qvalue_cutoff=None, enrichment_cutoff=None, length_cutoff=None):
        for i,peak in enumerate(self.peaks):
            if pvalue_cutoff is not None and peak.log_pvalue > pvalue_cutoff:
                self.peaks[i].filtered = True
            if qvalue_cutoff is not None and peak.log_qvalue > qvalue_cutoff:
                self.peaks[i].filtered = True
            if enrichment_cutoff is not None and peak.fold_enrichment < enrichment_cutoff:
                self.peaks[i].filtered = True
            if length_cutoff is not None and peak.length > length_cutoff:
                self.peaks[i].filtered = True
            
