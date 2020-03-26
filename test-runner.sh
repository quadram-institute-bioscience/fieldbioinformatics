#!/usr/bin/env bash
set -e
#
# test-runner.sh runs a the entire ARTIC field bioinformatics pipeline using a minimal set
# of data (the Mayinga barcode from an Ebola amplicon library sequenced on a flongle).
#
# full data available: http://artic.s3.climb.ac.uk/run-folders/EBOV_Amplicons_flongle.tar.gz
#
# usage:
#       ./test-runner.sh [medaka|nanopolish]
#
#   specify either medaka or nanopolish to run the respective branch of the pipeline
#
###########################################################################################
# Setup the data, commands and the testing function.

# data
inputData="../test-data/ebov-flongle/"
primerSchemes="../test-data/primer-schemes"
primerScheme="IturiEBOV/V1"
prefix="ebov"
barcode="03"
threads=2

# pipeline commands
gatherCmd="artic gather \
        --min-length 400 \
        --max-length 800 \
        --prefix ${prefix} \
        --directory ${inputData}"

demultiplexCmd="artic demultiplex \
            --threads ${threads} \
            ${prefix}_fastq_pass.fastq"

nanopolishCmd="nanopolish index \
            -d ${inputData} \
            -s ${prefix}_sequencing_summary.txt \
            ${prefix}_fastq_pass.fastq"

minionNanopolishCmd="artic minion \
                --normalise 200 \
                --threads ${threads} \
                --scheme-directory ${primerSchemes} \
                --read-file ${prefix}_fastq_pass-NB${barcode}.fastq \
                --nanopolish-read-file ${prefix}_fastq_pass.fastq \
                ${primerScheme} \
                ${prefix}"

minionMedakaCmd="artic minion \
            --normalise 200 \
            --threads ${threads} \
            --scheme-directory ${primerSchemes} \
            --read-file ${prefix}_fastq_pass-NB${barcode}.fastq \
            --medaka \
            ${primerScheme} \
            ${prefix}"

# colours
NC='\033[0m'
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'

# cmdTester is function to run a command and check for failure
function cmdTester {
    echo "###########################################################################################"
    echo -e "${BLUE}Running:${NC} $*"
    echo
    "$@"
    echo
    local status=$?
    if [ $status -ne 0 ]
    then
        echo -e "${RED}FAIL${NC}" >&2
    else
        echo -e "${GREEN}PASS${NC}" >&2
    fi
    echo
    return $status
}

###########################################################################################
# Run the tests.

# setup a tmp directory to work in
mkdir tmp && cd tmp || exit

# check that nanopolish or medaka is specified
if [ "$1" == "nanopolish" ] || [ "$1" == "medaka" ]; then
    echo -e "${BLUE}Starting tests...${NC}"
    echo -e "${BLUE} - using the $1 branch${NC}"
    echo
else
    echo "please specify medaka or nanopolish"
    echo "./test-runner.sh [medaka|nanopolish]"
    exit 1
fi

# collect the reads
cmdTester $gatherCmd

# demultiplex
cmdTester $demultiplexCmd

# run nanopolish or medaka branch
if [ "$1" == "nanopolish" ]
then
    cmdTester $nanopolishCmd
    cmdTester $minionNanopolishCmd
else
    cmdTester $minionMedakaCmd
fi

###########################################################################################
# Check the output and clean up.

# check output created
echo -e "${BLUE}Checking output...${NC}"
if test -f "${prefix}.consensus.fasta"
then
    echo -e "${GREEN} - consensus found${NC}"
else
    echo -e "${RED} - no consensus found${NC}"
    exit 1
fi

# TODO: add more checks....

# cleanup
cd .. && rm -r tmp
echo -e "${BLUE}Done.${NC}"