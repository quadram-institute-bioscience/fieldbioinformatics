---
title: tests
summary: Description of the available tests
authors:
  - Will Rowe
  - Nick Loman
date: 2020-03-30
---

All of the ARTIC tests are run as part of the Travis Continuous Integration testing. To run any of the tests yourself, it is assumed that you have downloaded the codebase and are using an appropiate environment:

```
git clone https://github.com/artic-network/fieldbioinformatics
cd fieldbioinformatics
conda env create -f environment.yml
conda activate artic
```

## Unit tests

We have begun writing unit tests for the ARTIC Python modules. These currently includes tests for the `align_trim` module, which performs amplicon soft clipping and alignment filtering, and `vcftagprimersites` which processes the ARTIC primer scheme. To run all available unit tests:

```
pytest -s artic/*_unit_test.py
```

## Pipeline tests

To test the complete workflow, as described in the [nCov SOP](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html), you can use the `test-runner.sh` bash script. You can test both the medaka and nanopolish workflows:

```
./test-runner.sh medaka
./test-runner.sh nanopolish
```

Both tests use a small subset of an Ebola virus amplicon sequencing run (flongle) which is included in the GitHub repository. You can download the full dataset [here](http://artic.s3.climb.ac.uk/run-folders/EBOV_Amplicons_flongle.tar.gz).

## Variant validation tests

Finally, we have also included some validation tests that will download several reference datasets, run the nanopolish/medaka workflows and then validate the reported variants. To run these:

```
pytest -s artic/minion_medaka_validator.py --remote-data
pytest -s artic/minion_nanopolish_validator.py --remote-data
```