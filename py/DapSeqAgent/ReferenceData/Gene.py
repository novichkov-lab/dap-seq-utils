
class Gene:
    def __init__(self, gene_id, genome_id, locus_tag=None, sequence_id=None, start=None, end=None, strand=None, utr_start=None, utr_end=None, operon_id = None):
        self.gene_id = gene_id
        self.genome_id = genome_id
        self.locus_tag = locus_tag
        self.sequence_id = sequence_id
        self.start = start
        self.end = end
        self.strand = strand
        self.utr_start = utr_start
        self.utr_end = utr_end
        self.operon_id = operon_id
        self.identifiers = {}
        self.description = ''

    def __str__(self):
        return "\t".join([self.gene_id, self.locus_tag, self.genome_id,\
                        str(self.start), str(self.end), self.strand, \
                        self.description])
    
