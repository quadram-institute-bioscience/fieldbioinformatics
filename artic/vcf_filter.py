import vcf
import sys
from operator import attrgetter
from collections import defaultdict
from .vcftagprimersites import read_bed_file

class NanoporeFilter:
    def __init__(self, *args, **kwargs):
        pass

    def check_filter(self, v):
        total_reads = float(v.INFO['TotalReads'])
        qual = v.QUAL
        strandbias = float(v.INFO['StrandFisherTest'])

        if qual / total_reads < 3:
            return False

        if len(v.ALT) > 1:
            print ("This code does not support multiple genotypes!")
            raise SystemExit

        ref = v.REF
        alt = v.ALT[0]

        if (len(alt) - len(ref) % 3):
            return False

        if v.is_indel:
            strand_fraction_by_strand = v.INFO['SupportFractionByStrand']
            if float(strand_fraction_by_strand[0]) < 0.5: 
                return False

            if float(strand_fraction_by_strand[1]) < 0.5:
                return False

        if total_reads < 20:
            return False

        return True

class MedakaFilter:
    def __init__(self, *args, **kwargs):
        pass

    def check_filter(self, v):
        if v.num_het:
            return False
        return True

class LongshotFilter:
    def __init__(self, *args, **kwargs):
        pass

    def check_filter(self, v):
        depth = v.INFO['DP']
        if depth < 20:
            return False

        if v.num_het:
            return False
        return True

def go(args):
    vcf_reader = vcf.Reader(filename=args.inputvcf)
    vcf_writer = vcf.Writer(open(args.output_pass_vcf, 'w'), vcf_reader)
    vcf_writer_filtered = vcf.Writer(open(args.output_fail_vcf, 'w'), vcf_reader)
    if args.nanopolish:
        filter = NanoporeFilter()
    elif args.medaka:
        filter = MedakaFilter()
    elif args.longshot:
        filter = LongshotFilter()
    else:
        print("Please specify a VCF type, i.e. --nanopolish or --medaka\n")
        raise SystemExit

    for v in vcf_reader:
        if filter.check_filter(v):
            vcf_writer.write_record(v)
        else:
            vcf_writer_filtered.write_record(v)

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--nanopolish', action='store_true')
    parser.add_argument('--medaka', action='store_true')
    parser.add_argument('--longshot', action='store_true')
    parser.add_argument('inputvcf')
    parser.add_argument('output_pass_vcf')
    parser.add_argument('output_fail_vcf')

    args = parser.parse_args()

    go(args)

if __name__ == "__main__":
    main()


