#!/usr/bin/python
import sys,argparse
from DapSeqAgent.DapSeqAgent import DapSeqAgent

def get_args():
    desc = '''This program runs DAP-seq pipeline.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-t', dest='task', type=str, help='Path to task file')
    parser.add_argument('-i', dest='indir', type=str, help='Input directory')
    parser.add_argument('-o', dest='outdir', type=str, help='Output directory')
    parser.add_argument('-s', dest='sample', type=str, default=None,
                        help='Sample ID (optional)')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args

def main():
    args = get_args()
    if args.sample is None:
        agent = DapSeqAgent(args.task, args.indir, args.outdir)
    else:
        agent = DapSeqAgent(args.task, args.indir, args.outdir, sample = args.sample)
    agent.run_task()
    print('Done!')

if __name__=='__main__':
    main()

