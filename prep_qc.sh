#!/bin/bash
RESULTS_DIR=${1}

cd ${RESULTS_DIR}
mkdir -p qc

# cp full-alg-v-segs.fa qc/full-alg-v-segs.fa
# cp most-freq-v-segs.fa qc/most-freq-v-segs.fa

cd qc
# fasta-vdj-pipeline-unified.pl --file=full-alg-v-segs.fa
# fasta-vdj-pipeline-unified.pl --file=most-freq-v-segs.fa

source /data/projects/rishi/sw/anaconda2/bin/activate
cd /data/projects/rishi/vdjones
python qc.py ${RESULTS_DIR} > ${RESULTS_DIR}/qc/stats.txt