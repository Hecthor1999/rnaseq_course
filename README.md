# RNA_seq_course_breast_cancer
Repository for the documented code used in the analysis performed during the RNA sequencing course

Detection of differential gene expression from bulk RNA sequencing data of three different breast cancer subtypes
# 1.Getting started
(Step carried out in the IBU cluster)

creation of working spaces:
```
$ cd /data/Users/harribas
$ mkdir RNA_seq_course
$ mkdir breast_cancer_project
$ cd RNA_seq_course/breast_cancer_project/
$ mkdir reference_genome
$ mkdir mapping
$ mkdir index_bam
$ mkdir sam-bam
$ mkdir feature_counts
$ mkdir scripts
```
create softlink to reads folder where the data is stored
```
$ ln -s /data/courses/rnaseq_course/breastcancer_de/reads ./
```

# 2.Quality check 
(Step carried out in the IBU cluster)

create script
```
$ vim run_fastqc.sh
```
run script through job submission
```
$ sbatch run_fastqc.sh
```
# 3.Mapping reads to the reference genome
(Step carried out in the IBU cluster)
```
$ cd ../reference_genome/
```
Upload the latest reference genome sequence and associated annotation, Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz (2023-04-21) and Homo_sapiens.GRCh38.110.gtf
Unzip reference genome
```
$ gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
change the names of Homo_sapiens.GRCh38.dna.primary_assembly.fa to genome.fa and Homo_sapiens.GRCh38.110.gtf to genome.gtf (for easier use) 
```
$ mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
$ mv Homo_sapiens.GRCh38.110.gtf genome.gtf
```
create script to index reference genome and submit it
```
$ vim index_without_transcripts.sh
$ sbatch index_without_transcripts.sh
```
---------------------------------------------
For mapping

go to directory, create file samplelist.tsv with paths to reads for slurm-array
```
$ cd ../mapping/
$ touch samplelist.tsv
```
create and run mapping script
```
$ vim array_mapping_without_transcripts
$ sbatch array_mapping_without_transcripts
```
--------------------------------------------------
sam to bam 

update samplelist.tsv with sam files location and move it to scripts folder, go to sam-bam directory
```
$ cd ../sam-bam/
```
create and run array_sam_to_bam.sh which converts sam files to bam files and sorts bam files.
```
$ vim array_sam_to_bam.sh
$ sbatch array_sam_to_bam.sh
```
---------------------------------------------------
index bam files

go to index_bam directory
```
$ cd ../index_bam
```
create and run array_index_bam.sh
```
$ vim array_index_bam.sh
$ sbatch array_index_bam.sh
```
download bam files with scp command to inspect with Integrative Genomic Viewer (IGV)

# 4.Count the number of reads per gene
(Step carried out in the IBU cluster)

go to feature_counts directory, create and run feature_counts.sh
```
$ cd ../feature_counts
$ vim feature_counts.sh
$ sbatch feature_counts.sh
```
Download featurecounts.txt to local machine with scp command for further downstream analysis

# 5. Exploratory analysis, differential expression analysis, Overrepresentation analysis
(Steps carried out on local machine)

create metadata.txt
create and run R script DeSeq2_pipeline.R















