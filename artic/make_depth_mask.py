#!/usr/bin/env python
from Bio import SeqIO
import sys
import vcf
import subprocess
from collections import defaultdict
import os.path
import operator
from .vcftagprimersites import read_bed_file
import itertools

def collect_depths(bamfile):
    if not os.path.exists(bamfile):
        raise SystemExit("bamfile %s doesn't exist" % (bamfile,))

    print(bamfile, file=sys.stderr)

    p = subprocess.Popen(['samtools', 'depth', bamfile],
                             stdout=subprocess.PIPE)
    out, err = p.communicate()
    depths = defaultdict(dict)
    for ln in out.decode('utf-8').split("\n"):
       if ln:
          contig, pos, depth = ln.split("\t")
          depths[contig][int(pos)] = int(depth)
    return depths

# from https://www.geeksforgeeks.org/python-make-a-list-of-intervals-with-sequential-numbers/
def intervals_extract(iterable): 
    iterable = sorted(set(iterable)) 
    for key, group in itertools.groupby(enumerate(iterable), 
    lambda t: t[1] - t[0]): 
        group = list(group) 
        yield [group[0][1], group[-1][1]] 

def go(args):
    depths = collect_depths(args.bamfile)

    mask_vector = []

    seq = list(SeqIO.parse(args.reference, "fasta"))[0]
    cons = list(seq.seq)

    for n, c in enumerate(cons):
        try:
            depth = depths[seq.id][n+1]
        except:
            depth = 0

        if depth < args.depth:
            mask_vector.append(n)

    intervals = list(intervals_extract(mask_vector))

    maskfh = open(args.outfile, 'w')
    for i in intervals:
        maskfh.write("%s\t%s\t%s\n" % (seq.id, i[0]+1, i[1]+1))
    maskfh.close()

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--depth', type=int, default=20)
    parser.add_argument('reference')
    parser.add_argument('bamfile')
    parser.add_argument('outfile')

    args = parser.parse_args()
    go(args)

if __name__ == "__main__":
    main()

