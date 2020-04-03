# minion_medaka_test.py is a funcitonal test for the minion pipeline command when run using the medaka option
import argparse
from Bio import SeqIO
import os
import pytest
import sys
import vcf

from . import pipeline


# help pytest resolve where test data is kept
TEST_DIR = os.path.dirname(os.path.abspath(__file__))

# run the test with a subsampled read set known to contain a deletion to the reference
testSample = "MT007544"

# set up the medaka command
medakaCMD = [
    "minion",
    "--medaka",
    "--read-file",
    TEST_DIR + "/../test-data/reads/MT007544.deletion.fastq",
    "--scheme-directory",
    TEST_DIR + "/../test-data/primer-schemes",
    "nCoV-2019/V1",
    testSample
]

# expectedDel is the deleted reference bases expected from these test reads
expectedDel = "ACGATCGAGTG"

# expectedDelPos is the start position of the deletion in the reference
expectedDelPos = 29749


def test_minion_medaka():
    """Test for the minion run function with the medaka pipeline.
    """

    # setup a parser
    parser = pipeline.init_pipeline_parser()

    # parse the arguments
    try:
        args = parser.parse_args(medakaCMD)
    except SystemExit:
        print("failed to parse valid command for `artic minion --medaka`")
        assert False

    # run
    try:
        args.func(parser, args)
    except SystemExit:
        print("artic minion exited early with an error")
        assert False

    # check longshot has reported the expected variant
    vcffn = "%s.longshot.vcf" % testSample
    assert os.path.exists(
        vcffn) == True, "test did not produce longshot vcf file"
    records = list(vcf.Reader(open(vcffn, 'r')))
    assert len(records) == 1, "incorrect number of variants reported by longshot"
    assert records[0].POS == expectedDelPos, "incorrect variant POS reported by longshot"
    assert records[0].REF == expectedDel, "incorrect variant REF reported by longshot"

    # check that the deletion to the reference has been added to the consensus
    fastafn = "%s.consensus.fasta" % testSample
    assert os.path.exists(
        fastafn) == True, "test did not produce consensus file"
    for record in SeqIO.parse(open(fastafn, 'r'), 'fasta'):
        assert expectedDel not in record.seq, "consensus not updated with known INDEL"
