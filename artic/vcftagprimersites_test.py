# test_vcftagprimersites.py is the unit tests for the vcf primer site tagging
import pytest
import os

from artic import vcftagprimersites

# pytest needs something like this to resolve where test data is kept:
TEST_DIR = os.path.dirname(os.path.abspath(__file__))


def test_read_bed_file():

    # process the nCoV-2019 V3 primer scheme
    primerScheme = vcftagprimersites.read_bed_file(
        TEST_DIR + "/../test-data/primer-schemes/nCoV-2019/V3/nCoV-2019.scheme.bed")

    # check the the alts have been collapsed into a canonical primer site
    assert len(primerScheme) == 196, "alts were not collapsed"

    # check access to primer scheme fields
    for row in primerScheme:
        assert 'Primer_ID' in row, "failed to parse primer scheme for Primer_ID"
        assert 'PoolName' in row, "failed to parse primer scheme for PoolName"

    # process the nCoV-2019 V2 primer scheme
    # this scheme has a single alt that has replaced the original primer so no alts should be collapsed
    primerScheme2 = vcftagprimersites.read_bed_file(
        TEST_DIR + "/../test-data/primer-schemes/nCoV-2019/V2/nCoV-2019.bed")
    assert len(primerScheme2) == 196, "no alts should have been collapsed"

    # TODO: this is a starting point for the unit tests - more to come...
