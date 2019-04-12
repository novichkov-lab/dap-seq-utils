
class Utr:
    def __init__(self, utr_id = None, genome_id = None, sequence_id=None, start=None, end=None, strand=None, gene_id=None, operon_id = None):
        self.utr_id = utr_id
        self.genome_id = genome_id
        self.sequence_id = sequence_id
        self.start = start
        self.end = end
        self.strand = strand
        self.gene_id = gene_id
        self.operon_id = operon_id

    def __str__(self):
        return "\t".join([self.utr_id, self.genome_id, self.sequence_id,\
                        str(self.start), str(self.end), self.strand, \
                        self.gene_id])
    
