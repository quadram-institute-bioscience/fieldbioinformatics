#!/usr/bin/env python3
import sys
from Bio import SeqIO
import dnaio
import tempfile
import os
import glob
import gzip
import fnmatch
import shutil
# import pandas as pd
import numpy as np
from collections import defaultdict
from mimetypes import guess_type
from functools import partial
from math import log10
from random import random
from pathlib import Path
import click

import time

def get_mean_qscore(quality):
    phred = dict(zip("""!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~""",
                     range(0, 94)))
    phred_score = [phred[q] for q in quality]
    _numerator = [10**(score/-10) for score in phred_score]
    return -10 * log10(np.mean(_numerator))


def run(parser, args):
    # fastq_directory,min_length, max_length, phred_score, output
    start_time = time.time()
    fastq_files = [str(fq) for fq in Path(args.directory).iterdir() if fq.name.endswith('.fastq') or fq.name.endswith('.fastq.gz')]
    i = 0
    total_seqs = 0
    dups = set()

    if len(fastq_files) >=1 :
        if not args.output:
            fastq_outfn = "%s_%s.fastq.gz" % (args.prefix, os.path.basename(args.directory))
        else:
            fastq_outfn = args.output
        with dnaio.FastqWriter(fastq_outfn) as fo:
            for fastq_file in fastq_files:
                with dnaio.open(fastq_file) as fi:
                    for record in fi:
                        total_seqs += 1
                        seq_len = len(record.sequence)
                        if args.min_length <= seq_len <= args.max_length and get_mean_qscore(record.qualities) >= args.quality:
                            i += 1
                            # print(
                                # f"{record.name.split(' ')[0]}: {get_mean_qscore(record.qualities)}")
                            if record.name not in dups:
                                fo.write(record)
                                dups.add(record.name)
                        # if i == 10:
                            # break
    
    print(f"Finished in {time.time() - start_time} seconds", file=sys.stderr)
    print(f"Included reads: {i} out of {total_seqs}", file=sys.stderr)
    print(f"Duplicated reads {dups}", file=sys.stderr)

# if __name__ == "__main__":
#     run()
