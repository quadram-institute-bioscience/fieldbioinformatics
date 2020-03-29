<div align="center">
    <h1>ARTIC bioinformatics pipeline</h1>
    <h3>a pipeline for working with virus sequencing data sequenced with nanopore</h3>
    <hr>
    <a href="https://travis-ci.org/artic-network/fieldbioinformatics"><img src="https://travis-ci.org/artic-network/fieldbioinformatics.svg?branch=master" alt="travis"></a>
    <a href='http://artic.readthedocs.io/en/latest/?badge=latest'><img src='https://readthedocs.org/projects/artic/badge/?version=latest' alt='Documentation Status'></a>
    <a href="https://bioconda.github.io/recipes/artic/README.html"><img src="https://anaconda.org/bioconda/artic/badges/downloads.svg" alt="bioconda"></a>
    <a href="https://github.com/artic-network/fieldbioinformatics/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-orange.svg" alt="License"></a>
</div>

---

## Overview

> coming soon

## Installation

### conda

> coming soon

### source

```sh
git clone https://github.com/artic-network/fieldbioinformatics
cd fieldbioinformatics
python setup.py install
```

> note: there are a few dependencies required - there are a couple of minimal conda environments provided in this repository:

- for running the nanopolish version of the pipeline:

```sh
  conda env create -f environment.yml
  conda activate fieldbioinformatics
```

- for running the medaka version of the pipeline:

```sh
  conda env create -f environment-medaka.yml
  conda activate fieldbioinformatics-medaka
```

## Documentation

For documentation and getting started, see the [SOP on the artic website](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop-beta.html).

Documentation via [read the docs](http://artic.readthedocs.io/en/latest/?badge=latest) is being written and will be available soon.
