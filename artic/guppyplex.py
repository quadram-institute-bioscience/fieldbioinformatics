import sys
from Bio import SeqIO
import tempfile
import os
import glob
import fnmatch
import shutil
import pandas as pd
from collections import defaultdict
from mimetypes import guess_type
from functools import partial
from math import log10
from random import random

# the ambition of this module is to merge some of the functionality from the "gather" and "demultiplex" tasks;
# we are assuming the input is already guppy-demultiplexed data (one-pot barcoding - like) - this workflow will
# also assume that Medaka rather than Nanopolish will be used and will skip the sequencing_summary step - basic
# QC metrics will be derived during the fastq parsing step ...

# This method will also allow for the analysis of gzipped fastq files = makes sense of space

# Have tried to be minimally invasive to the existing code - maintain FieldBioinformatics style


def get_read_mean_quality(record):
    return -10 * log10((10 ** (pd.Series(record.letter_annotations["phred_quality"]) / -10)).mean())


def run(parser, args):
    fastq = defaultdict(list)
    files = os.listdir(args.directory)
    fastq[args.directory].extend([os.path.join(args.directory, f) for f in files if fnmatch.fnmatch(f, '*.fastq*')])

    for barcode_directory, fastq in list(fastq.items()):
        if len(fastq):
            fastq_outfn = "%s_%s.fastq" % (args.prefix, os.path.basename(barcode_directory))
            outfh = open(fastq_outfn, "w")
            print("Processing %s files in %s" % (len(fastq), barcode_directory), file=sys.stderr)

            dups = set()

            for file in fastq[0:len(fastq)-2]:
                encoding = guess_type(file)[1]
                _open = open
                # only accommodating gzip compression at present
                if encoding == "gzip":
                    _open = partial(gzip.open, mode="rt")
                with _open(file) as f:
                    try:
                        for rec in SeqIO.parse(f, "fastq"):
                            if args.max_length and len(rec) > args.max_length:
                                continue
                            if args.min_length and len(rec) < args.min_length:
                                continue
                            if not args.skip_quality_check and get_read_mean_quality(rec) < args.quality:
                                continue
                            if args.sample < 1:
                                r = random()
                                if r >= args.sample:
                                    continue

                            if rec.id not in dups:
                                SeqIO.write([rec], outfh, "fastq")
                                dups.add(rec.id)
                    except ValueError:
                       pass

            outfh.close()
            print(f"{fastq_outfn}\t{len(dups)}")
