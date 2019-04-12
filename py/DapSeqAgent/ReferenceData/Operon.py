
class Operon:
    def __init__(self, operon_id, genome_id, sequence_id, start=None, end=None, strand=None):
        self.operon_id = operon_id
        self.genome_id = genome_id
        self.sequence_id = sequence_id
        self.start = start
        self.end = end
        self.strand = strand
        self.gene_count = 0
        self.gene_ids = []

    def __str__(self):
        return self.operon_id + ':' + self.gene_ids[0] + '-' + self.gene_ids[-1]
