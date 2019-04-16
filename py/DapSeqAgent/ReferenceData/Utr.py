
class Utr:
    def __init__(self, utr_id = None, genome_id = None, sequence_id=None, start=None, end=None, strand=None, rev_gene_id = None, fwd_gene_id=None, rev_operon_id = None, fwd_operon_id = None):
        self.utr_id = utr_id
        self.genome_id = genome_id
        self.sequence_id = sequence_id
        self.start = start
        self.end = end
        self.strand = strand
        self.rev_gene_id = rev_gene_id
        self.fwd_gene_id = fwd_gene_id
        self.rev_operon_id = rev_operon_id
        self.fwd_operon_id = fwd_operon_id

    def __str__(self):
        if self.rev_gene_id is None:
            rev_gene = ''
        else:
            rev_gene = self.rev_gene_id
        if self.fwd_gene_id is None:
            fwd_gene = ''
        else:
            fwd_gene = self.fwd_gene_id
        if self.rev_operon_id is None:
            rev_operon = ''
        else:
            rev_operon = self.rev_operon_id
        if self.fwd_operon_id is None:
            fwd_operon = ''
        else:
            fwd_operon = self.fwd_operon_id
        return "\t".join([self.utr_id, self.genome_id, self.sequence_id,\
                        str(self.start), str(self.end), self.strand, \
                        rev_gene, fwd_gene, rev_operon, fwd_operon])
    
