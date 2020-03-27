#Written by Nick Loman (@pathogenomenick)

import os
import sys
from Bio import SeqIO
from clint.textui import colored, puts, indent
from .vcftagprimersites import read_bed_file

def get_nanopolish_header(ref):
    recs = list(SeqIO.parse(open(ref), "fasta"))
    if len (recs) != 1:
        print("FASTA has more than one sequence", file=sys.stderr)
        raise SystemExit(1)

    return  "%s:%d-%d" % (recs[0].id, 1, len(recs[0])+1)

def run(parser, args):
    log = "%s.minion.log.txt" % (args.sample)
    logfh = open(log, 'w')

    if args.scheme.find('/') != -1:
        scheme_name, scheme_version = args.scheme.split('/')
    else:
        scheme_name = args.scheme
        scheme_version = "V1"

    ref = "%s/%s/%s/%s.reference.fasta" % (args.scheme_directory, scheme_name, scheme_version, scheme_name)
    bed = "%s/%s/%s/%s.scheme.bed" % (args.scheme_directory, scheme_name, scheme_version, scheme_name)

    if args.read_file:
        read_file = args.read_file
    else:
        read_file = "%s.fasta" % (args.sample)

    if not os.path.exists(ref):
        print(colored.red('Scheme reference file not found: ') + ref)
        raise SystemExit(1)
    if not os.path.exists(bed):
        print(colored.red('Scheme BED file not found: ') + bed)
        raise SystemExit(1)

    pools = set([row['PoolName'] for row in read_bed_file(bed)])

    cmds = []

    nanopolish_header = get_nanopolish_header(ref)

    # 3) index the ref & align with bwa"
    if not args.bwa:
        cmds.append("minimap2 -a -x map-ont -t %s %s %s | samtools sort -o %s.sorted.bam -" % (args.threads, ref, read_file, args.sample))
    else:
        cmds.append("bwa index %s" % (ref,))
        cmds.append("bwa mem -t %s -x ont2d %s %s | samtools view -bS - | samtools sort -o %s.sorted.bam -" % (args.threads, ref, read_file, args.sample))
    cmds.append("samtools index %s.sorted.bam" % (args.sample,))

    # 4) trim the alignments to the primer start sites and normalise the coverage to save time
    if args.normalise:
        normalise_string = '--normalise %d' % (args.normalise)
    else:
        normalise_string = ''
#    if args.medaka:
#        cmds.append("align_trim --no-read-groups --start %s %s --report %s.alignreport.txt < %s.sorted.bam 2> %s.alignreport.er | samtools sort -T %s - -o %s.trimmed.sorted.bam" % (normalise_string, bed, args.sample, args.sample, args.sample, args.sample, args.sample))
#        cmds.append("align_trim %s %s --no-read-groups --report %s.alignreport.txt < %s.sorted.bam 2> %s.alignreport.er | samtools sort -T %s - -o %s.primertrimmed.sorted.bam" % (normalise_string, bed, args.sample, args.sample, args.sample, args.sample, args.sample))
#        cmds.append("samtools index %s.trimmed.sorted.bam" % (args.sample))
#        cmds.append("samtools index %s.primertrimmed.sorted.bam" % (args.sample))
#    else:
    cmds.append("align_trim --start %s %s --report %s.alignreport.txt < %s.sorted.bam 2> %s.alignreport.er | samtools sort -T %s - -o %s.trimmed.rg.sorted.bam" % (normalise_string, bed, args.sample, args.sample, args.sample, args.sample, args.sample))
    cmds.append("align_trim %s %s --remove-incorrect-pairs --report %s.alignreport.txt < %s.sorted.bam 2> %s.alignreport.er | samtools sort -T %s - -o %s.primertrimmed.rg.sorted.bam" % (normalise_string, bed, args.sample, args.sample, args.sample, args.sample, args.sample))
    cmds.append("samtools index %s.trimmed.rg.sorted.bam" % (args.sample))
    cmds.append("samtools index %s.primertrimmed.rg.sorted.bam" % (args.sample))

    if args.medaka:
       for p in pools:
          cmds.append("samtools view -b -r \"%s\" %s.primertrimmed.rg.sorted.bam > %s.primertrimmed.%s.sorted.bam" % (p, args.sample, args.sample, p))
          cmds.append("samtools index %s.primertrimmed.%s.sorted.bam" % (args.sample, p))

    # 6) do variant calling using the raw signal alignment
    if args.medaka:
        for p in pools:
            if os.path.exists("%s.%s.hdf" % (args.sample, p)):
                os.remove("%s.%s.hdf" % (args.sample, p))
            cmds.append("medaka consensus %s.primertrimmed.%s.sorted.bam %s.%s.hdf" % (args.sample, p, args.sample, p))
            cmds.append("medaka variant %s %s.%s.hdf %s.%s.vcf" % (ref, args.sample, p, args.sample, p))

        #cmds.append("margin_cons_medaka --depth 20 --quality 10 %s %s.primertrimmed.medaka.vcf %s.primertrimmed.sorted.bam > %s.consensus.fasta 2> %s.report.txt" % (ref, args.sample, args.sample, args.sample, args.sample))
    else:
        if not args.skip_nanopolish:
            if args.nanopolish_read_file:
                    indexed_nanopolish_file = args.nanopolish_read_file
            else:
                    indexed_nanopolish_file = read_file

            for p in pools:
               cmds.append("nanopolish variants -x %s --progress -t %s --reads %s -o %s.%s.vcf -b %s.trimmed.rg.sorted.bam -g %s -w \"%s\" --ploidy 1 -m 0.15 --read-group %s" % (args.max_haplotypes, args.threads, indexed_nanopolish_file, args.sample, p, args.sample, ref, nanopolish_header, p))

    merge_vcf_cmd = "artic_vcf_merge %s %s" % (args.sample, bed)
    for p in pools:
        merge_vcf_cmd += " %s:%s.%s.vcf" % (p, args.sample, p)
    cmds.append(merge_vcf_cmd)

    if args.medaka:
        cmds.append("longshot -P 0.001 -F -A --no_haps --bam %s.primertrimmed.rg.sorted.bam --ref %s --out %s.longshot.vcf --potential_variants %s.merged.vcf" % (args.sample, ref, args.sample, args.sample))
        cmds.append("artic_vcf_filter --longshot %s.longshot.vcf %s.filtered.vcf" % (args.sample, args.sample))
    else:
        cmds.append("artic_vcf_filter --nanopolish %s.merged.vcf %s.filtered.vcf" % (args.sample, args.sample))

    cmds.append("artic_make_depth_mask %s %s.primertrimmed.rg.sorted.bam %s.coverage_mask.txt" % (ref, args.sample, args.sample))

    vcf_file = "%s.filtered.vcf" % (args.sample,)
    cmds.append("bgzip -f %s" % (vcf_file))
    cmds.append("tabix -p vcf %s.gz" % (vcf_file))

    cmds.append("bcftools consensus -f %s %s.gz -m %s.coverage_mask.txt -I -o %s.consensus.fasta" % (ref, vcf_file, args.sample, args.sample))

    if args.medaka:
        method = 'medaka'
    else:
        method = 'nanopolish'
    fasta_header = "%s/ARTIC/%s" % (args.sample, method)

    cmds.append("artic_fasta_header %s.consensus.fasta \"%s\"" % (args.sample, fasta_header))

            #python nanopore-scripts/expand-cigar.py --bam "$sample".primertrimmed.sorted.bam --fasta $ref | python nanopore-scripts/count-errors.py /dev/stdin > "$sample".errors.txt

            # 7) do phasing
            #nanopolish phase-reads --reads $sample.fasta --bam $sample.trimmed.sorted.bam --genome $ref $sample.vcf

            # 8) variant frequency plot
            #cmds.append("vcfextract %s > %s.variants.tab" % (args.sample, args.sample))

            # 8) filter the variants and produce a consensus
            # here we use the vcf file without primer binding site trimming (to keep nanopolish happy with flanks)
            # but we use the primertrimmed sorted bam file in order that primer binding sites do not count
            # for the depth calculation to determine any low coverage sites that need masking
            #cmds.append("margin_cons %s %s.vcf %s.primertrimmed.sorted.bam a > %s.consensus.fasta" % (ref, args.sample, args.sample, args.sample))

    for cmd in cmds:
# print(colored.green("Running: ") + cmd, file=sys.stderr)
        print(cmd)
        print(cmd, file=logfh)
        if not args.dry_run:
            retval = os.system(cmd)
            if retval != 0:
                print(colored.red('Command failed:' ) + cmd, file=sys.stderr)
                raise SystemExit(20)

    logfh.close()

