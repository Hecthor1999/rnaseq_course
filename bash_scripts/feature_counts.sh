#!/bin/bash

#SBATCH --job-name="feature_counts"
#SBATCH --output=featureCounts.out
#SBATCH --error=featureCounts.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=20:00:00
#SBATCH --mem=10G
#SBATCH --partition=pall

module load UHTS/Analysis/subread/2.0.1;

bam_path="/data/users/harribas/RNA_seq_course/breast_cancer_project/sam-bam"
ref_path="/data/users/harribas/RNA_seq_course/breast_cancer_project/reference_genome"

featureCounts -p -t exon -g gene_id -a ${ref_path}/genome.gtf -o featurecounts.txt ${bam_path}/*sorted.bam 
