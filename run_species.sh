#!/bin/bash

SPECIES=("dragonfish" "rockcod" "killifish" "llama-1" "llama-2" "alpaca-1" "alpaca-2")


INFILE=/data/projects/rishi/final-paper/alpaca-1/PFMID15.vd.clean.VDJ.H3.L3.CH1.dnaH3.fa
OUTDIR=/data/projects/rishi/final-paper/species/alpaca-1
bsub -n 4 -e ${OUTDIR}/err -o ${OUTDIR}/out "/data/projects/rishi/vdjones/runner.sh ${INFILE} ${OUTDIR}"

INFILE=/data/projects/rishi/final-paper/alpaca-2/PFMID16.vd.clean.VDJ.H3.L3.CH1.dnaH3.fa
OUTDIR=/data/projects/rishi/final-paper/species/alpaca-2
bsub -n 4 -e ${OUTDIR}/err -o ${OUTDIR}/out "/data/projects/rishi/vdjones/runner.sh ${INFILE} ${OUTDIR}"

INFILE=/data/projects/rishi/final-paper/llama-2/PFMID14.vd.clean.VDJ.H3.L3.CH1.dnaH3.fa
OUTDIR=/data/projects/rishi/final-paper/species/llama-2
bsub -n 4 -e ${OUTDIR}/err -o ${OUTDIR}/out "/data/projects/rishi/vdjones/runner.sh ${INFILE} ${OUTDIR}"

INFILE=/data/projects/rishi/final-paper/llama-1/PFMID13.vd.clean.VDJ.H3.L3.CH1.dnaH3.fa
OUTDIR=/data/projects/rishi/final-paper/species/llama-1
bsub -n 4 -e ${OUTDIR}/err -o ${OUTDIR}/out "/data/projects/rishi/vdjones/runner.sh ${INFILE} ${OUTDIR}"

INFILE=/data/projects/rishi/final-paper/dragonfish/dragonfish-raw.clean.VDJ.H3.L3.CH1.dnaH3.fa
OUTDIR=/data/projects/rishi/final-paper/species/dragonfish
bsub -n 4 -e ${OUTDIR}/err -o ${OUTDIR}/out "/data/projects/rishi/vdjones/runner.sh ${INFILE} ${OUTDIR}"

INFILE=/data/projects/rishi/final-paper/killifish/killifish-vdj.dna.clean.VDJ.H3.L3.CH1.dnaH3.fa
OUTDIR=/data/projects/rishi/final-paper/species/killifish
bsub -n 4 -e ${OUTDIR}/err -o ${OUTDIR}/out "/data/projects/rishi/vdjones/runner.sh ${INFILE} ${OUTDIR}"

INFILE=/data/projects/rishi/final-paper/rockcod/rockcod-raw.clean.VDJ.H3.L3.CH1.dnaH3.fa
OUTDIR=/data/projects/rishi/final-paper/species/rockcod
bsub -n 4 -e ${OUTDIR}/err -o ${OUTDIR}/out "/data/projects/rishi/vdjones/runner.sh ${INFILE} ${OUTDIR}"




