# align_trim_test.py is the unit tests for alignment trimming
import pytest
import os
import pysam
from copy import copy
from cigar import Cigar

from artic import align_trim
from artic import vcftagprimersites

# pytest needs something like this to resolve where test data is kept:
TEST_DIR = os.path.dirname(os.path.abspath(__file__))

# dummy primers (using min required fields)
p1 = {
        "start": 0,
        "end": 10,
        "direction": "+",
        "primerID": "primer1_LEFT"
}
p2 = {
        "start": 30,
        "end": 40,
        "direction": "-",
        "primerID": "primer1_RIGHT"
}
p3 = {
        "start": 10,
        "end": 20,
        "direction": "+",
        "primerID": "primer2_LEFT"
}
p4 = {
        "start": 40,
        "end": 50,
        "direction": "-",
        "primerID": "primer2_RIGHT"
}

# primer scheme to hold dummy primers
dummyPrimerScheme = [p1, p2, p3, p4]

# actual the primer scheme for nCov
primerScheme = vcftagprimersites.read_bed_file(
        TEST_DIR + "/../test-data/primer-schemes/nCoV-2019/V3/nCoV-2019.scheme.bed")

# dummy alignment segment
seg = pysam.AlignedSegment()
seg.query_name = "0be29940-97ae-440e-b02c-07748edeceec"
seg.flag = 0
seg.reference_id = 0
seg.reference_start = 4294
seg.mapping_quality = 60
seg.cigarstring = "40S9M1D8M2D55M1D4M2I12M1I14M1D101M1D53M2D78M1I60M52S"
seg.query_sequence = "CAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTACAGTGCTTAAAAGTGTAAGTGCCTTTTACATTCTACCATCTATTATCTCTAATGAGAAGCAAGAAATTCTTGAACCTTCATACTTGGAATTTTGCGAGAAATGCTGCACATGCAGAAGAAACACGCAAATTAATGCCTGTCTGTGTGGAAACTAAAGCCATAGTTTCAACTATACAGCGTAAATATAAGGGTATTAAAATACAAGGGGTGTGGTTGATTATGGTGCTAGATTTTACTTTTACACCAGTAAAACAACTGGCGTCACTTATCAACACACTTAACGATCTAAATGAAACTCTTGTTACAATGCCACTTGGCTATGTAACACATGGCTTAGAATTTGGAAGAAGCTGCTCGGTATATGAGATCTCTCAAAGTGCCAGCTACAGTTTCTGTTGCGATTGCTGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACCCATG"
seg.query_qualities = [30] * 490

# test for the find_primer function
def test_find_primer():

    # test the primer finder on the primers themselves
    for primer in dummyPrimerScheme:
        if primer["direction"] == "+":
            result = align_trim.find_primer(
                dummyPrimerScheme, primer["start"], primer["direction"])
        else:
            result = align_trim.find_primer(
                dummyPrimerScheme, primer["end"], primer["direction"])
        assert result[2]["primerID"] == primer["primerID"], "find_primer did not produce the query primer, which should be nearest"

    # test against other ref positions
    result = align_trim.find_primer(
        dummyPrimerScheme, 8, "+")
    assert result[2]["primerID"] == "primer2_LEFT", "find_primer returned incorrect primer"
    result = align_trim.find_primer(
        dummyPrimerScheme, 25, "-")
    assert result[2]["primerID"] == "primer1_RIGHT", "find_primer returned incorrect primer"



# test for the soft_mask function
def test_mask():

    # get the nearest primers to the alignment segment
    p1 = align_trim.find_primer(primerScheme, seg.reference_start, '+')
    p2 = align_trim.find_primer(primerScheme, seg.reference_end, '-')

    # get the end of the first primer
    primer_position = p1[2]['end']

    # this should need a mask
    numToMask = primer_position - seg.reference_start
    assert numToMask > 0, "missed a soft masking opportunity"

    # but it is already softmasked in this region, so only the start pos should have changed
    oldStart = seg.reference_start
    masked = align_trim.soft_mask(seg, primer_position, False, True)
    assert masked == False and seg.reference_start == (oldStart + numToMask), "soft masking shouldn't be needed, only a reference_start increment"

    # get the end of the second primer
    primer_position = p2[2]['start']

    # this should need a mask
    numToMask = seg.reference_end - primer_position
    assert numToMask > 0, "missed a soft masking opportunity"
    oldStart = seg.reference_start
    masked = align_trim.soft_mask(seg, primer_position, True, True)

    # again, should not need softmasking
    assert masked == False and seg.reference_start == oldStart, "soft masking or reference_start increment shouldn't be needed"
