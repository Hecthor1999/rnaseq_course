#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=06:00:00
#SBATCH --job-name=hisat2_index_without_transcripts
#SBATCH --output=hisat2_index.out
#SBATCH --error=hisat2_index.err
#SBATCH --mail-user=hector.arribasarias@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --partition=pall

module load UHTS/Aligner/hisat/2.2.1;

hisat2-build -p 16 genome.fa genome
