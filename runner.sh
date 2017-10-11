#!/bin/bash

FASTA_FILE=${1}
RESULTS_DIR=${2}

mkdir -p ${RESULTS_DIR}
cp ${FASTA_FILE} ${RESULTS_DIR}/input_seqs.fa
cd /data/projects/rishi/vdjones
source ./env/bin/activate
git rev-parse --short HEAD >> algorithm.info
python vdjones.py ${RESULTS_DIR}/input_seqs.fa ${RESULTS_DIR}
