from DapSeqAgent.ReferenceData.Utr import Utr

class Contig:
    def __init__(self, contig_id, size=0):
        self.contig_id = contig_id
        self.size = size
        self.operons = []
        self.genes = []
        self.utrs = []

    def add_gene(self,gene):
        self.genes.append(gene)
    
    def add_operon(self,operon):
        self.operons.append(operon)
        
    def calculate_utrs(self):
        if len(self.genes) == 0:
            print('No genes in contig', self.contig_id, 'UTR search skipped')
            return
        upstream_gene = None
        upstream_gene_end = 0 # Sequence positions start from 1
        utr_index = 1
        for gene_index,gene in enumerate(self.genes):
            if gene.strand == '+':
                utr_start = upstream_gene_end + 1
                utr_end = gene.start - 1
                if utr_end > utr_start + 1:
                    operon_id = None
                    for operon in self.operons:
                        if gene.gene_id in operon.gene_ids:
                            operon_id = operon.operon_id
                    utr = Utr(utr_id = 'utr_' + str(utr_index),
                            genome_id = gene.genome_id,
                            sequence_id = gene.sequence_id,
                            start = utr_start,
                            end = utr_end,
                            strand = gene.strand,
                            gene_id = gene.gene_id)
                    if operon_id is not None:
                        utr.operon_id = operon_id
                    self.utrs.append(utr)
                    
            elif gene.strand == '-':
                if gene_index == len(self.genes) - 1: # the last gene
                    utr_start = gene.end + 1
                    utr_end = self.size
                else:
                    utr_start = gene.end + 1
                    utr_end = self.genes[gene_index + 1].start - 1
                if utr_end > utr_start + 1:
                    operon_id = None
                    for operon in self.operons:
                        if gene.gene_id in operon.gene_ids:
                            operon_id = operon.operon_id
                    
                    utr = Utr(utr_id = 'utr_' + str(utr_index),
                            genome_id = gene.genome_id,
                            sequence_id = gene.sequence_id,
                            start = utr_start,
                            end = utr_end,
                            strand = gene.strand,
                            gene_id = gene.gene_id)
                    if operon_id is not None:
                        utr.operon_id = operon_id
                    self.utrs.append(utr)
            upstream_gene_end = gene.end
