
class Genome:
    def __init__(self, genome_id, name=None, ncbi_taxonomy=None):
        self.genome_id = genome_id
        self.name = name
        self.ncbi_taxonomy = ncbi_taxonomy
        self.size = 0
        self.contigs = {}
