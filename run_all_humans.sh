# !/bin/bash
HUMANS=(091 104 130 158 182 219 231 272 275 276)

for i in "${HUMANS[@]}"
do
    OUTDIR=/data/projects/rishi/final-paper/human/${i}
    INFILE=/data/projects/rishi/humans-separated/cluster/${i}/${i}-all.fa
    IMGT_ANNOT_FILE=/data/projects/rishi/humans-separated/original-fasta-files/${i}-all.fa
    bsub -n 8 -e ${OUTDIR}/err -o ${OUTDIR}/out "/data/projects/rishi/vdjones/runner.sh ${INFILE} ${OUTDIR} ${IMGT_ANNOT_FILE}"
done



######## FOR QC; doesn't work properly

# HUMANS=(091)

# for i in "${HUMANS[@]}"
# do
#     OUTDIR=/data/projects/rishi/final-paper/human/${i}
#     bsub -n 1 -L /bin/bash -e ${OUTDIR}/qc/err -o ${OUTDIR}/qc/out "/bin/bash -l -c /data/projects/rishi/vdjones/prep_qc.sh ${OUTDIR}"
# done


