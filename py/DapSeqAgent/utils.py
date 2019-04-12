import operator
from collections import defaultdict

def autovivify(levels=1, final=dict):
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))

def sanitize_file_name(filename):
    filename = filename.replace(' ', '_')
    filename = filename.replace("'", "")
    filename = filename.replace('"', '')
    return filename

