#!/usr/bin/python3
import os,json
from DapSeqAgent.MacsParser.Peak import Peak
from DapSeqAgent.Task import Task

def decode_object(o):
    if '__Peak__' in o:
        a = Peak()
        a.__dict__.update(o['__Peak__'])
        return a
    elif '__Task__' in o:
        a = Task()
        a.__dict__.update(o['__Task__'])
        return a
    return o

def decode_task(o):
    if '__Peak__' in o:
        a = Peak()
        a.__dict__.update(o['__Peak__'])
        return a
    elif '__Site__' in o:
        a = Site()
        a.__dict__.update(o['__Site__'])
        return a
    elif '__Motif__' in o:
        a = Motif()
        a.__dict__.update(o['__Motif__'])
        return a
    elif '__Task__' in o:
        a = Task()
        a.__dict__.update(o['__Task__'])
        return a
    return o

def import_task(infile):
    deserialized = None
    with open (infile, 'r') as f:
        deserialized = json.load(f, object_hook=decode_object)
        f.closed
    return deserialized

def export_task(task, outfile):
    #outfile = os.path.join(task.outdir, 'task.json')
    #print pretty JSON: print(json.dumps(task,indent=4, cls=CustomEncoder))
    with open (outfile, 'w') as of:
        json.dump(task,of,cls=CustomEncoder)
        of.closed

class CustomEncoder(json.JSONEncoder):

     def default(self, o):

         return {'__{}__'.format(o.__class__.__name__): o.__dict__}

