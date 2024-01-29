#!/bin/bash
#SBATCH --job-name=fastQC_samples
#SBATCH --output=/data/users/harribas/RNA_seq_course/breast_cancer_project/fastqc/fastQC_samples.out
#SBATCH --error=/data/users/harribas/RNA_seq_course/breast_cancer_project/fastqc/fastQC_samples.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --time=01:30:00
#SBATCH --partition=pall

module add UHTS/Quality_control/fastqc/0.11.9

mkdir -p /data/users/harribas/RNA_seq_course/breast_cancer_project/fastqc
srun fastqc --extract /data/users/harribas/RNA_seq_course/breast_cancer_project/reads/* --threads 1 --outdir /data/users/harribas/RNA_seq_course/breast_cancer_project/fastqc
