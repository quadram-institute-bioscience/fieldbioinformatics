# align_trim_test.py is the unit tests for alignment trimming
import pytest
import os

from artic import align_trim


def test_find_primer():

    # create some dummy primers
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

    # create a dummy bedfile of primers
    bedFile = [p1, p2, p3, p4]

    # test the primer finder on the primers themselves
    for primer in bedFile:
        if primer["direction"] == "+":
            result = align_trim.find_primer(
                bedFile, primer["start"], primer["direction"])
        else:
            result = align_trim.find_primer(
                bedFile, primer["end"], primer["direction"])
        assert result[2]["primerID"] == primer["primerID"], "find_primer did not produce the query primer, which should be nearest"

    # test against other ref positions
    result = align_trim.find_primer(
        bedFile, 8, "+")
    assert result[2]["primerID"] == "primer2_LEFT", "find_primer returned incorrect primer"
    result = align_trim.find_primer(
        bedFile, 25, "-")
    assert result[2]["primerID"] == "primer1_RIGHT", "find_primer returned incorrect primer"
