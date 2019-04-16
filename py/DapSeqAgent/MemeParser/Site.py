class Site:
    def __init__(self, site_id=None, sequence=None, start=None, end=None, strand=None):
        self.site_id = site_id
        self.start = start
        self.end = end
        self.strand = strand
        self.sequence = sequence
