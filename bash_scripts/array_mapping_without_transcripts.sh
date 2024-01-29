#!/bin/bash

#SBATCH --array=1-12
#SBATCH --time=13:00:00
#SBATCH --mem=10000M
#SBATCH --cpus-per-task=4
#SBATCH --job-name=hisat2_mapping
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pall

# define variables
WORKDIR="/data/users/harribas/RNA_seq_course/breast_cancer_project/mapping"
OUTDIR="$WORKDIR/results"
SAMPLELIST="$WORKDIR/samplelist.tsv"

SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

OUTFILE="$OUTDIR/${SAMPLE}.sam"

mkdir -p $OUTDIR

module load UHTS/Aligner/hisat/2.2.1
hisat2 -p 8 -x /data/users/harribas/RNA_seq_course/breast_cancer_project/reference_genome/genome -1 $READ1 -2 $READ2 -S $OUTFILE
