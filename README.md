<div align="center">
    <h1>ARTIC pipeline</h1>
    <h3>a bioinformatics pipeline for working with virus sequencing data sequenced with nanopore</h3>
    <hr>
    <a href="https://travis-ci.org/artic-network/fieldbioinformatics"><img src="https://travis-ci.org/artic-network/fieldbioinformatics.svg?branch=master" alt="travis"></a>
    <a href='http://artic.readthedocs.io/en/latest/?badge=latest'><img src='https://readthedocs.org/projects/artic/badge/?version=latest' alt='Documentation Status'></a>
    <a href="https://bioconda.github.io/recipes/artic/README.html"><img src="https://anaconda.org/bioconda/artic/badges/downloads.svg" alt="bioconda"></a>
    <a href="https://github.com/artic-network/fieldbioinformatics/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-orange.svg" alt="License"></a>
</div>

---

## Overview

The `artic pipeline` is designed to help run the artic bioinformatics protocols; for example the [nCoV-2019 novel coronavirus protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html).

Features include:

* read filtering
* primer trimming
* amplicon coverage normalisation
* variant calling
* consensus building

There are **2 variants** of the pipeline, one which uses signal data (via [nanopolish](https://github.com/jts/nanopolish)) and one that does not (via [medaka](https://github.com/nanoporetech/medaka)).

## Installation

### conda

> coming soon

### source

#### installing the pipeline:

```sh
git clone https://github.com/artic-network/fieldbioinformatics
cd fieldbioinformatics
python setup.py install
```

#### installing dependencies:

The `artic pipeline` has several software dependencies. These dependencies are different for the `nanopolish` and `medaka` variants of the pipeline. We have provided some minimal conda environments to manage these dependencies.

Unfortunately, the required `nanopolish` and `medaka` versions are not compatible with each other, so **you will need 2 environments** if you want to try both the `nanopolish` and `medaka` variants of the pipeline.

```
# nanopolish:
conda env create -f test-data/env-nanopolish.yml && conda activate artic-nanopolish

# medaka:
conda env create -f test-data/env-medaka.yml && conda activate artic-medaka
```

>**Note:** currently, the latest conda nanopolish recipe does not support OSX, so OSX users must manage this dependency another way

#### testing:

First check the pipeline can be called.

```
artic -v
```

Now try the unit tests

```
pytest -s artic/*_unit_test.py
```

You can also try the functional tests

```
./test-runner.sh nanopolish

or

./test-runner.sh medaka
```


## Documentation

For documentation and getting started, see the [SOP on the artic website](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html).

Documentation via [read the docs](http://artic.readthedocs.io/en/latest/?badge=latest) is being written and will be available soon.
