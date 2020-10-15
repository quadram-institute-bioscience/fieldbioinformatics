#!/usr/bin/env python3
# Script for filtering out short reads for coronaHiT method
# Thanh Le Viet - Quadram Institute Bioscience, 2020

from Bio import SeqIO
import pysam
import logging
from logging import DEBUG, INFO
import pathlib


def mismatch_rate(segment, ref_seq):
    """
    Parameters
    ----------
        ref_seq : SeqIO
            Reference sequence 
        segment: AlignmentFile object
            Aligment 
    Returns
    ------- 
        float
            Mismatched rate
    """
    alg_pos = segment.get_aligned_pairs(matches_only=True)
    mm = 0
    for read_pos, ref_pos in alg_pos:
        if segment.query_sequence[read_pos] != ref_seq[ref_pos]:
            mm +=1
    mismatch_rate = (mm/segment.query_alignment_length)*100
    return mismatch_rate


def filter_chit_reads(reference, infile, output_name=None, mismatch_rate_threshold=25, mapq=1):
    """
    Parameters
    ----------
        reference : Fasta file
            Reference file 
        infile: Bam file
            Alignment bam file
        output_name: Value
            Value for output bamfile name, if not defined, infile name will be used
        mismatch_rate_threshold: Value
            Threshold of mismatch rate per aligned read
        mapq: Value
            MAPQ value
    Returns
    ------- 
        bamfile
            Two bam files: The bam file with less mismatched reads and high MAPQ will be used for variant calling and the other file is for \
            storing excluded reads
    """
    record = list(SeqIO.parse(reference, "fasta"))[0]
    ref_seq = record.seq
    
    bamFile = pysam.AlignmentFile(infile, 'rb')
    bam_header = bamFile.header.copy().to_dict()
    bamFile_path = pathlib.Path(infile).parent
    logging.debug(f"In bamfile: {infile}")
    logging.debug(f"Mismatch rate: {mismatch_rate_threshold}\tMAPQ:{mapq}")
    if output_name is None:
        output_name = pathlib.Path(infile).name
    logging.debug(f"Output name: {output_name}")
    inc_bam = f"{output_name}.primertrimmed.rg.sorted.bam"
    logging.debug(f"Passed reads will be written to: {inc_bam}")
    exc_bam = f"{output_name}.chit.excluded.bam"
    logging.debug(f"Passed reads will be written to: {exc_bam}")
    inc_bam_file = pysam.AlignmentFile(inc_bam, "wb", header=bam_header)
    exc_bam_file = pysam.AlignmentFile(exc_bam, "wb", header=bam_header)
    for segment in bamFile:
        if segment.get_tag('RG') != 'unmatched':
            inc_bam_file.write(segment)

        elif segment.get_tag('RG') == 'unmatched' and segment.mapq >= mapq:
            mmr = mismatch_rate(segment, ref_seq)
            if mmr >= mismatch_rate_threshold:
                 exc_bam_file.write(segment)
            else:
                inc_bam_file.write(segment)
        else:
            exc_bam_file.write(segment)
    
    exc_bam_file.close()
    inc_bam_file.close()
    logging.debug(f"Finished filtering!")
    try:
        pysam.index(inc_bam)
        logging.debug(f"Created index for {inc_bam}")
        pysam.index(exc_bam)
        logging.debug(f"Created index for {exc_bam}")
    except:
        logging.debug(f"Building index for output bam files failed.")

    bamFile.close()

def go(args):
    if args.verbose:
        logging.basicConfig(format='%(name)s-%(levelname)s-%(message)s', datefmt='%d-%b-%y %H:%M:%S', level=DEBUG)

    assert pathlib.Path(args.bamfile).exists(), f"{args.bamfile} is not existed."
    assert pathlib.Path(args.reference).exists(), f"{args.bamfile} is not existed."
    filter_chit_reads(reference=args.reference, infile=args.bamfile, output_name=args.output_name, mismatch_rate_threshold=args.mismatch_rate, mapq=args.mapq)


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--mismatch-rate', type=float, default=25.0, help=r"Mismatch rate (default: %(default)s)")
    parser.add_argument('--mapq', type=int, default=1, help="MAPQ threshold (default: %(default)s)")
    parser.add_argument('--output-name', type=str, default=None, help="Name for output bam files (default: %(default)s)", required=False)
    parser.add_argument('-v','--verbose', action='store_true', default=True, help="Debug mode")
    parser.add_argument('reference')
    parser.add_argument('bamfile')
    args = parser.parse_args()
    go(args)

if __name__ == "__main__":
    main()