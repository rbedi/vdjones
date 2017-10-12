#!/bin/bash

FASTA_FILE=${1}
RESULTS_DIR=${2}

mkdir -p ${RESULTS_DIR}
cp ${FASTA_FILE} ${RESULTS_DIR}/input_seqs.fa

if [ ! -z "$3" ]
  then
    IMGT_ANNOT_FILE=${3}
    cp ${IMGT_ANNOT_FILE} ${RESULTS_DIR}/input_seqs_imgt_annot.fa
fi

# source /data/projects/rishi/sw/anaconda2/bin/activate

cd /data/projects/rishi/vdjones
git rev-parse --short HEAD >> ${RESULTS_DIR}/algorithm.info
/data/projects/rishi/sw/anaconda2/bin/python vdjones.py ${RESULTS_DIR}/input_seqs.fa ${RESULTS_DIR}
