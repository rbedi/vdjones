HUMANS=(091 104 130 158 182 219 231 272 275 276)

for i in "${HUMANS[@]}"
do
    OUTDIR=/data/projects/rishi/final-paper/human/${i}
    INFILE=/data/projects/rishi/humans-separated/original-fasta-files/${i}-all.fa
    bsub -n 8 -e ${OUTDIR}/err -o ${OUTDIR}/out "/data/projects/rishi/vdjones/runner.sh ${INFILE} ${OUTDIR}"
done


