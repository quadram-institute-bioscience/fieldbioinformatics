#!/usr/bin/env python

import vcf
import sys
import subprocess
import csv
from collections import defaultdict


def merge_sites(canonical, alt):
    """Merges a canonical primer site with an alt site, producing an interval that encompasses both

    Parameters
    ----------
    canonical : dict
        The canonical primer site, provided as a dictionary of the bed file row
    alt : dict
        The alt primer site, provided as a dictionary of the bed file row

    Returns
    -------
    dict
        A dictionary of the merged site, where the dict represents a bed file row
    """
    # set up a new merged site to return
    mergedSite = {}
    mergedSite['Primer_ID'] = canonical['Primer_ID']
    mergedSite['start'] = canonical['start']
    mergedSite['end'] = canonical['end']

    # check the both the canonical and alt are the same direction
    if canonical['direction'] != alt['direction']:
        print(
            "could not merge alt with different orientation to canonical", file=sys.stderr)
        raise SystemExit(1)

    # merge the start/ends of the alt with the canonical to get the largest window possible
    if canonical['direction'] == '+':
        if alt['start'] < canonical['start']:
            mergedSite['start'] = alt['start']
        if alt['end'] > canonical['end']:
            mergedSite['end'] = alt['end']
    else:
        if alt['start'] > canonical['start']:
            mergedSite['start'] = alt['start']
        if alt['end'] < canonical['end']:
            mergedSite['end'] = alt['end']
    return mergedSite


def read_bed_file(fn):
    """Parses a bed file and collapses alts into canonical primer sites

    Parameters
    ----------
    fn : str
        The bedfile to parse

    Returns
    -------
    list
        A list of dictionaries, where each dictionary contains a row of the parsed bedfile.
        The available dictionary keys are - Primer_ID, direction, start, end
    """
    # use a dictionary of dictionaries to hold all the rows in the bed file by primer ID
    # this can then be used to collapse the alts into canonical primer sites
    bedFile = {}

    with open(fn) as csvfile:
        reader = csv.reader(csvfile, dialect='excel-tab')
        for row in reader:

            # remove empty fields
            row = list(filter(None, row))

            # read row of bed file into a dictionary
            bedrow = {}
            bedrow['Primer_ID'] = row[3]
            bedrow['PoolName'] = row[4]

            # check the bed format
            if len(row) >= 6:
                # new style bed
                bedrow['direction'] = row[5]
            elif len(row) == 5:
                # old style without directory
                if 'LEFT' in row[3]:
                    bedrow['direction'] = '+'
                elif 'RIGHT' in row[3]:
                    bedrow['direction'] = '-'
                else:
                    print("Malformed BED file!", file=sys.stderr)
                    raise SystemExit(1)
            else:
                print("Malformed BED file!", file=sys.stderr)
                raise SystemExit(1)

            # grab the direction and set the start and end of the site
            if bedrow['direction'] == '+':
                bedrow['end'] = int(row[2])
                bedrow['start'] = int(row[1])
            else:
                bedrow['end'] = int(row[1])
                bedrow['start'] = int(row[2])

            # get the primer base name (removes the alt tag)
            # NOTE: alts are assumed to have an ID ending in "_alt*"
            primerID = row[3].split('_alt')[0]

            # check if this primer ID is already in the dict
            if primerID not in bedFile:

                # add to the bed file and continue
                bedFile[primerID] = bedrow
                continue

            # otherwise, we've got a primer ID we've already seen so merge the alt
            mergedSite = merge_sites(bedFile[primerID], bedrow)

            # update the bedFile
            bedFile[primerID] = mergedSite

    # return the bedFile as a list
    return list(bedFile.values())


def overlaps(coords, pos):
    for v in coords:
        if pos >= v['start'] and pos <= v['end']:
            return v
    return False


if __name__ == "__main__":
    if sys.argv[1] not in sets:
        print("Invalid set")
        raise SystemExit(1)

    bedfile = read_bed_file(sys.argv[1])

    vcf_reader = vcf.Reader(filename=sys.argv[2])
    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
    for record in vcf_reader:
        v = overlaps(bedfile, record.POS)
        if v:
            record.INFO['PRIMER'] = v["Sequence_(5-3')"]

#	PP = list(record.INFO)
#	record.INFO = {}
#	record.INFO['PP'] = PP
#	record.INFO['DEPTH'] = depths[record.CHROM][record.POS]

        vcf_writer.write_record(record)
