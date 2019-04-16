import csv
from collections import defaultdict

class Task:
    def __init__(self, infile=None):
        self.taskfile = infile
        self.samples=defaultdict(dict)
        if infile is not None:
            self.load_task()
        
    def load_task(self):
        with open(self.taskfile, 'r') as f:
            tsvin = csv.reader(f, delimiter='\t')
            for row in tsvin:
                self.samples[row[0]]['ref_library'] = row[1]
                self.samples[row[0]]['protein'] = row[2]
                self.samples[row[0]]['replicate'] = row[3]
                self.samples[row[0]]['treatment'] = row[4]
                self.samples[row[0]]['control'] = row[5]
                self.samples[row[0]]['fastq'] = row[6]
                self.samples[row[0]]['info'] = {}
            
    
