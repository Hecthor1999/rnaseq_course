#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4000M
#SBATCH --cpus-per-task=1
#SBATCH --job-name=array_index_bam
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pall

# define variables
WORKDIR="/data/users/harribas/RNA_seq_course/breast_cancer_project"
OUTDIR="$WORKDIR/index_bam"
SAMPLELIST="$WORKDIR/scripts/samplelist.tsv"

SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
BAM=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $5; exit}' $SAMPLELIST`

OUTFILE="$OUTDIR/${SAMPLE}_indexed.bam"

############################

module load UHTS/Analysis/samtools/1.10
samtools index $BAM
