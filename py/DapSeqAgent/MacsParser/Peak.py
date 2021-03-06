
class Peak:
    def __init__(self, peak_id=None, seq_id=None, start=None, end=None, abs_summit=None, length=None, pileup=None, log_pvalue=None, log_qvalue=None, fold_enrichment=None, protein_id = None):
        self.peak_id = peak_id
        self.sequence_id = seq_id
        self.start = start
        self.end = end
        self.abs_summit = abs_summit
        self.length = length
        self.pileup = pileup
        self.log_pvalue = log_pvalue
        self.log_qvalue = log_qvalue
        self.fold_enrichment = fold_enrichment
        self.protein_id = protein_id
        self.genes = []
        self.operons = []
        self.utrs = []
        self.status = ''
        self.noncoding = None
        self.filtered = False
        self.gene_functions = ''
        self.utr_functions = ''

    def __str__(self):
        return "\t".join([self.peak_id, self.sequence_id,\
                        str(self.start), str(self.end), str(self.length),\
                        str(self.abs_summit), str(self.pileup), str(self.log_pvalue),\
                        str(self.log_qvalue), str(self.fold_enrichment),\
                        self.protein_id])

    def annotate_peak(self, ref_library):
        self.genes = ref_library.get_genes(self.sequence_id, self.start, self.end)
        self.utrs = ref_library.get_utrs(self.sequence_id, self.start, self.end)
        operons = set()
        for utr_id in self.utrs:
            utr = ref_library.get_utr_by_id(utr_id)
            if utr.rev_operon_id:
                operons.add(utr.rev_operon_id)
            if utr.fwd_operon_id:
                operons.add(utr.fwd_operon_id)
        for gene_id in self.genes:
            operon_id = ref_library.get_operon_id_by_first_gene(gene_id)
            if operon_id is not None:
                operons.add(operon_id)
#        self.operons = ref_library.get_operons(self.sequence_id, self.start, self.end)
        self.operons = list(operons)
        self.status = ''
        if len(self.genes) != 0 or len(self.operons) != 0:
            self.status = 'Intragenic'
            if len(self.utrs) != 0:
                self.status = 'Upstream/Intragenic'
        elif len(self.utrs) != 0:
            self.status = 'Upstream'
#        self.status = ref_library.get_interval_status(self.sequence_id, self.start, self.end)#is_noncoding(self.sequence_id, self.abs_summit)
        if len(self.genes) > 0:
            self.gene_functions = ';'.join([ref_library.get_fuction_by_gene_id(self.sequence_id, x) for x in self.genes])
        if len(self.utrs) > 0:
            self.utr_functions = ';'.join([ref_library.get_fuction_by_gene_id(self.sequence_id, x) for x in self.utrs])
        pass
