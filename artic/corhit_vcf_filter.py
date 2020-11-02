#!/usr/bin/env python3
# Script for filtering out duplicated snps with lowest QUAL scores in VCF for coronaHiT method 
# Code were modified from vcf_filter.py by Nick Loman
# Thanh Le Viet - Quadram Institute Bioscience, 2020

import vcf
from operator import attrgetter
from collections import defaultdict
import logging
import operator

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
        indx = "%s-%s" % (v.CHROM, v.POS)
        #Find max QUAL in the duplicated snps
        max_qual_snp = max(group_variants[indx], key=operator.attrgetter('QUAL'))
        if len(group_variants[indx]) > 1:
            if (v.INFO['Pool'] == max_qual_snp.INFO['Pool']):
                vcf_keep_writer.write_record(v)
                logging.debug(f"Keep puplicated {indx} POS: {v.POS} ALT: {v.ALT} INFO: {v.INFO['Pool']} {v.var_type} QUAL:{v.QUAL}")
            else:
                vcf_remove_writer.write_record(v)
                logging.debug(
                    f"Remove duplicated {indx} POS: {v.POS} ALT: {v.ALT} INFO: {v.INFO['Pool']} {v.var_type} QUAL:{v.QUAL}")
        else:
            vcf_keep_writer.write_record(v)                
        
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


