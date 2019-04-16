import re
from DapSeqAgent.MemeParser.Motif import Motif
from DapSeqAgent.MemeParser.Site import Site

class MemeParser:
    def __init__(self, memefile, peaks_dict):
        self.memefile = memefile
        self.motifs = []
        self.parse_meme_output(peaks_dict)
        
    def parse_meme_output(self, peaks_dict):
        sitesection_begins = re.compile('\s*Motif \d+ sites sorted by position p-value\s*')
        sitesection_ends = re.compile('\s*Motif \d+ block diagrams\s*')
        site_lines = []
        sites_flag = False
        motif_number = 0
        with open(self.memefile, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')
                if re.match(sitesection_begins, line):
                    sites_flag = True
                elif re.match(sitesection_ends, line):
                    motif_number += 1
                    self.motifs.append(parse_sites('motif_' + str(motif_number), site_lines, peaks_dict))
                    sites_flag = False
                    site_lines = []
                elif sites_flag:
                    site_lines.append(line)
    

def parse_sites(motif_id, site_lines, peaks_dict):
    motif = Motif(motif_id)
    for line in site_lines:
        line_tokens = line.split()
        if len(line_tokens) < 6:
            continue
        peak_index = line_tokens[0]
        if peak_index not in peaks_dict:
            continue
        peak = peaks_dict[peak_index]
        site_id = peak.peak_id
        strand = line_tokens[1]
        sequence = line_tokens[5]
        start = peak.start + int(line_tokens[2])
        end = start + len(sequence)
        site = Site(site_id=site_id, sequence=sequence, start=start, end=end, strand=strand)
        motif.sites.append(site)
    return motif
