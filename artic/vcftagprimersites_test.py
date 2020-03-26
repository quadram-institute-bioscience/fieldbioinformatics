# test_vcftagprimersites.py is the unit tests for the vcf primer site tagging
import pytest
import os

from artic import vcftagprimersites

# pytest needs something like this to resolve where test data is kept:
TEST_DIR = os.path.dirname(os.path.abspath(__file__))


def test_read_bed_file():

    # process the nCoV-2019 primer scheme
    bedfile = vcftagprimersites.read_bed_file(
        TEST_DIR + "/../test-data/primer-schemes/nCoV-2019/V3/nCoV-2019.scheme.bed")

    # check the the alts have been collapsed into a canonical primer site
    assert len(bedfile) == 196, "alts where not collapsed"

    # TODO: this is a starting point for the unit tests - more to come...
