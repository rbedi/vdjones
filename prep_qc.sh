cp dominant.fa qc/full-alg-v-segs.fa
cp most_freq.fa qc/most-freq-v-segs.fa

cd qc
fasta-vdj-pipeline-unified.pl --file=full-alg-v-segs.fa
fasta-vdj-pipeline-unified.pl --file=most-freq-v-segs.fa


