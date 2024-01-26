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
$ ln -s /data/courses/rnaseq_course/breastcancer_de/reads ./
$ mkdir reference_genome
$ mkdir mapping
$ mkdir index_bam
$ mkdir sam-bam
$ mkdir feature_counts
$ mkdir scripts
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
$ cd reference_genome/
```
Upload the latest reference genome sequence and associated annotation, Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz (2023-04-21) and Homo_sapiens.GRCh38.110.gtf
Unzip reference genome
```
$ gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
change the name of Homo_sapiens.GRCh38.dna.primary_assembly.fa to genome.fa (for easier use)
```
$ mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
```
create script to index reference genome and submit it
```
$ vim index_without_transcripts.sh
$ sbatch index_without_transcripts.sh
```
---------------------------------------------
mapping
go to directory, create file samplelist.tsv with paths to reads for slurm-array
```
$ cd mapping/
$ touch samplelist.tsv
```
create and run mapping script
```
$ 
```



# 4.Count the number of reads per gene
(Step carried out in the IBU cluster)



















