#!/usr/bin/env python3
# Script for trimming primers for coronaHiT data
# Some chunks of code were copied from the original align_trim.py
# Thanh Le Viet - Quadram Institute Bioscience, 2020

from Bio import SeqIO
import pysam
import logging
from logging import DEBUG, INFO
import subprocess
import pathlib
import sys
import operator
from collections import defaultdict

from .align_trim import find_primer
from .vcftagprimersites import read_bed_file


def trim(infile, bedfile, reference, output_name, singularity_samtool=None, threads=2):
    if output_name is None:
        output_name="amplicon_clipped"
    
    samtools_ampliconclip_cmd = f"samtools ampliconclip -@ {threads} -b {bedfile} {infile} --reference {reference} --strand -O bam --both-end --original -f {output_name}.ampliconclip.stats"
    pipe_cmd = "|"
    samtools_sort = f"samtools sort -@ {threads} -o {output_name}.sorted.trimmed.bam -"
    SINGULARITY_SAMTOOLS = ""
    if singularity_samtool is not None:
        SINGULARITY_SAMTOOLS = f"singularity exec {singularity_samtool}"
    cmd = f"{SINGULARITY_SAMTOOLS} {samtools_ampliconclip_cmd} {pipe_cmd} {SINGULARITY_SAMTOOLS} {samtools_sort}".strip()
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.check_returncode()
    

def go(args):
    if args.verbose:
        logging.basicConfig(format='%(name)s-%(levelname)s-%(message)s',
                            datefmt='%d-%b-%y %H:%M:%S', level=DEBUG)
    
    scheme_file = f"{args.scheme_directory}/nCoV-2019.scheme.bed"
    bedfile = f"{args.scheme_directory}/nCoV-2019.bed"
    reference = f"{args.scheme_directory}/nCoV-2019.reference.fasta"

    assert pathlib.Path(args.bamfile).exists(
        ), f"{args.bamfile} is not existed."

    assert pathlib.Path(args.scheme_directory).exists(
    ), f"{args.scheme_directory} is not existed."
    
    assert pathlib.Path(scheme_file).exists(
    ), f"{scheme_file} is not existed."

    trim(args.bamfile, bedfile=bedfile, reference=reference,
         output_name=args.output_name, singularity_samtool=args.singularity_samtools, threads = args.threads)
    
    # prepare the report outfile
    if args.report:
        reportfh = open(args.report, "w")
        print("QueryName\tReferenceStart\tReferenceEnd\tPrimerPair\tPrimer1\tPrimer1Start\tPrimer2\tPrimer2Start\tIsSecondary\tIsSupplementary\tStart\tEnd\tCorrectlyPaired\tPool", file=reportfh)

    # open the primer scheme and get the pools
    bed = read_bed_file(scheme_file)
    primer_pos = defaultdict(dict)
    for b in bed:
        primer_pos[b['Primer_ID']]['start'] = b['start']
        primer_pos[b['Primer_ID']]['end'] = b['end']
    
    # print(primer_pos)

    pools = set([row['PoolName'] for row in bed])
    pools.add('unmatched')

    # open the input SAM file and process read groups
    infile = pysam.AlignmentFile(f"{args.output_name}.sorted.trimmed.bam", "rb")
    bam_header = infile.header.copy().to_dict()
    bam_header['RG'] = []
    for pool in pools:
        read_group = {}
        read_group['ID'] = pool
        bam_header['RG'].append(read_group)
    # prepare the alignment outfile
    outfile = pysam.AlignmentFile(f"{args.output_name}.primertrimmed.rg.sorted.bam", "wb", header=bam_header)
    # iterate over the alignment segments in the input SAM file
    # set up a counter to track amplicon abundance
    counter = defaultdict(int)
    for segment in infile:
        # filter out unmapped and supplementary alignment segments
        if segment.is_unmapped:
            print("%s skipped as unmapped" %
                  (segment.query_name), file=sys.stderr)
            continue
        if segment.is_supplementary:
            print("%s skipped as supplementary" %
                  (segment.query_name), file=sys.stderr)
            continue

        # locate the nearest primers to this alignment segment
        p1 = find_primer(bed, segment.reference_start, '+')
        p2 = find_primer(bed, segment.reference_end, '-')
        # logging.debug(f"Primer 1: {p1}")
        # logging.debug(f"Primer 1: {p2}")
        correctly_paired = p1[2]['Primer_ID'].replace(
            '_LEFT', '') == p2[2]['Primer_ID'].replace('_RIGHT', '')

        pool_name = 'unmatched'
        #Tag primer into bam file
        tags = segment.get_tags()
        tags += [("P1", p1[2]['Primer_ID']), ("P2", p2[2]['Primer_ID'])]
        # segment.set_tags(tags)
        if correctly_paired:
            if not segment.is_reverse: #Forward
                if int(segment.reference_start) >= int(p1[2]['end']) and int(segment.reference_end) <= int(p2[2]['start']):
                    pool_name = p1[2]['PoolName']
            else: #Reverse
                if int(segment.reference_end) <= int(p2[2]['start']) and int(segment.reference_start) >= int(p1[2]['end']):
                    pool_name = p1[2]['PoolName']
            tags += [("PL", p1[2]['Primer_ID']), ("PR", p2[2]['Primer_ID'])] #Set calibrate primer tags
            segment.set_tag('RG', pool_name)
        else:
            if not segment.is_reverse: #Forward
                if int(segment.reference_start) < int(p1[2]['start']):
                    if int(segment.reference_end <= int(p2[2]['start'])):
                        left_primer = p2[2]['Primer_ID'].replace('RIGHT', 'LEFT')
                        if int(segment.reference_start) >= primer_pos[left_primer]['end']:
                            pool_name = p2[2]['PoolName']
                            tags += [("PL", left_primer),
                                     ("PR", p2[2]['Primer_ID'])]
                    #     else:
                    #         pool_name = 'unmatched'
                    # else:
                    #     pool_name = 'unmatched'
                else:
                    right_primer = p1[2]['Primer_ID'].replace('LEFT', 'RIGHT')
                    if int(segment.reference_end) <= primer_pos[right_primer]['start']:
                        pool_name = p1[2]['PoolName']
                        tags += [("PL", p1[2]['Primer_ID']),
                                 ("PR", right_primer)]
                    # pool_name = p1[2]['PoolName']
            else: #Reverse
                if int(segment.reference_end) > int(p2[2]['start']):
                    right_primer = p1[2]['Primer_ID'].replace('LEFT', 'RIGHT')
                    if int(segment.reference_start) > int(p1[2]['start']):
                        if int(segment.reference_end) <= primer_pos[right_primer]['start']:
                            pool_name = p1[2]['PoolName']
                            tags += [("PL", p1[2]['Primer_ID']),
                                     ("PR", right_primer)]
                    # else:
                    #     pool_name = 'unmatched'
                else:
                    left_primer = p2[2]['Primer_ID'].replace('RIGHT', 'LEFT')
                    if int(segment.reference_start) >= primer_pos[left_primer]['end']:
                        pool_name = p2[2]['PoolName']
                        tags += [("PL", left_primer),("PR", p2[2]['Primer_ID'])]
                # segment.set_tag('RG',pool_name)
            # update the report with this alignment segment + primer details
            # normalise if requested
        tags = [('RG', pool_name)] + tags
        segment.set_tags(tags)
        if args.normalise:
            # print(tags)
            left_p = tags[:-1][1]
            right_p = tags[:-3][1]
            pair = "%s-%s-%d" % (left_p, right_p, segment.is_reverse)
            counter[pair] += 1
            if counter[pair] > args.normalise:
                print("%s dropped as abundance theshold reached" %
                      (segment.query_name), file=sys.stderr)
                continue
        report = f"{segment.query_name}\t{segment.reference_start}\t{segment.reference_end}\t{p1[2]['Primer_ID']}_{p2[2]['Primer_ID']}\t{p1[2]['Primer_ID']}\t{abs(p1[1])}\t{p2[2]['Primer_ID']}\t{abs(p2[1])}\t{segment.is_secondary}\t{segment.is_supplementary}\t{p1[2]['start']}\t{p2[2]['end']}\t{correctly_paired}\t{segment.get_tag('RG')}"
        if args.report:
            print(report, file=reportfh)
            
        outfile.write(segment)
    infile.close()
    outfile.close()


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-name', type=str, default=None, help="Name for output bam files (default: %(default)s)", required=False)
    parser.add_argument('--report', type=str, help='Output report to file')
    parser.add_argument('--normalise', type=int,
                        help='Subsample to n coverage per strand')
    parser.add_argument('--singularity-samtools', type=str, default=None,
                        help="Full absolute path to singularity samtools >= 1.11 (default: %(default)s)", required=False)
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of threads (default: %(default)d)')
    parser.add_argument('-v','--verbose', action='store_true', default=True, help="Debug mode")
    parser.add_argument('bamfile')
    parser.add_argument('scheme_directory')
    # parser.add_argument('reference')
    
    args = parser.parse_args()
    go(args)

if __name__ == "__main__":
    main()
