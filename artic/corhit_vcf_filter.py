#!/usr/bin/env python3
# Script for filtering out VCF for coronaHiT method
# Filter duplicated snp/indel in unmatched
# Code were modified from vcf_filter.py by Nick Loman
# Thanh Le Viet - Quadram Institute Bioscience, 2020

import vcf
from operator import attrgetter
from collections import defaultdict
import logging

def go(args):
    if args.verbose:
        logging.basicConfig(format='%(asctime)s-%(levelname)s-%(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.DEBUG)
    
    vcf_reader = vcf.Reader(filename=args.input_vcf)
    # vcf for longshot
    vcf_keep_writer = vcf.Writer(open(args.keep_vcf, 'w'), template=vcf_reader)
    # vcf for low quality SNPs/indel in unmatched
    vcf_remove_writer = vcf.Writer(
        open(args.remove_vcf, 'w'), template=vcf_reader)

    variants = [v for v in vcf_reader]

    group_variants = defaultdict(list)
    for v in variants:
        indx = "%s-%s" % (v.CHROM, v.POS)
        group_variants[indx].append(v)

    for v in variants:
        if v.POS==4850:
            logging.debug(f"FUCK IT POS: {v.POS} REF: {v.REF[0]}:{len(v.REF[0])} ALT: {v.ALT[0]}:{len(v.ALT[0])} INFO: {v.INFO['Pool']} {v.FORMAT}")
        indx = "%s-%s" % (v.CHROM, v.POS)
        if len(group_variants[indx]) > 1:
            if v.INFO['Pool'] != 'unmatched':
                vcf_keep_writer.write_record(v)
                logging.debug(f"Duplicated {indx} POS: {v.POS} ALT: {v.ALT} INFO: {v.INFO['Pool']} {v.var_type}")
            else:
                vcf_remove_writer.write_record(v)
        else:
            if v.INFO['Pool'] != 'unmatched':
                vcf_keep_writer.write_record(v)                
            else:
                if v.var_type == "snp":
                    vcf_keep_writer.write_record(v)
                    logging.debug(
                        f"Keep<=={indx} POS: {v.POS} REF: {v.REF[0]}:{len(v.REF[0])} ALT: {v.ALT[0]}:{len(v.ALT[0])} INFO: {v.INFO['Pool']} {v.FORMAT}")
                else:
                    vcf_remove_writer.write_record(v)
                    logging.debug(
                    f"Remove==>{indx} POS: {v.POS} REF: {v.REF[0]}:{len(v.REF[0])} ALT: {v.ALT[0]}:{len(v.ALT[0])} INFO: {v.INFO['Pool']} {v.FORMAT}")

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('input_vcf')
    parser.add_argument('keep_vcf')
    parser.add_argument('remove_vcf')
    args = parser.parse_args()

    go(args)

if __name__ == "__main__":
    main()


