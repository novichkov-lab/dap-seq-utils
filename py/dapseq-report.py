#!/usr/bin/python
import sys,argparse
from DapSeqAgent.DapSeqAgent import DapSeqAgent
from DapSeqAgent.JSONUtil import import_task

def get_args():
    desc = '''This program runs DAP-seq pipeline.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-t', dest='task', type=str, help='Path to task file')
    parser.add_argument('-i', dest='indir', type=str, help='Input directory')
    parser.add_argument('-o', dest='outdir', type=str, help='Output directory')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args

def main():
    args = get_args()
    agent = DapSeqAgent(args.task, args.indir, args.outdir)
    if os.path.exists(os.path.join(args.outdir, 'task.json')):
        agent.task = import_task(os.path.join(args.outdir, 'task.json'))
        agent.task.load_task()
        agent.generate_report()
    print('Done!')

if __name__=='__main__':
    main()

