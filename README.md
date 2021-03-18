# RNA-seq
### RNA-seq alignment pipeline

### Table of Contents

1. [Installation](https://github.com/ben-laufer/RNA-seq#installation)
2. [Indexing](https://github.com/ben-laufer/RNA-seq#indexing)
3. [Download FASTQ Files](https://github.com/ben-laufer/RNA-seq#download-fastq-files)
4. [Alignment Pipeline](https://github.com/ben-laufer/RNA-seq#main-pipeline)
   1. [Cluster](https://github.com/ben-laufer/RNA-seq#cluster)
   2. [Node](https://github.com/ben-laufer/RNA-seq#node)
5. [Quality Control](https://github.com/ben-laufer/RNA-seq#quality-control)
5. [Differential Gene Expression Analysis](https://github.com/ben-laufer/RNA-seq#differential-gene-expression-analysis)

## Installation

This workflow utilizes the following packages, which need to be installed and in your path:
1. [Trim Galore!](https://github.com/FelixKrueger/TrimGalore)
2. [STAR](https://github.com/alexdobin/STAR)
3. [Samtools](http://www.htslib.org)
4. [MultiQC](http://multiqc.info)

I recommend using [Bioconda](https://bioconda.github.io) to install and manage the package updates, which can be accomplished by:

`conda install -c trim-galore star samtools multiqc`

The node version of the alignment script also requires GNU parellel, which can be installed by:

`conda install -c parallel`

## Indexing

Both the genome (FASTA) and gene annotations (GTF) are required to create the STAR indexes. An example of how to download the required files and create the index for 150 bp paired end reads for mm10 is provided in the `star_index_mm10.sh` script. The indexes are specific to both species and read length. You will have to edit `export mainPath` to be the main directory on your cluster.

The complete genome folder structure should appear as:

```
├── genomes
│   ├── mm10
│   │   ├── star_150
```

## Download FASTQ Files

An example of how to download and check FASTQ files from SLIMS is provided in `01-getFASTQ.sh`.

Line 4 allows you to specify the project directory name, while lines 6 and 7 will require your SLIMS information (replace the examples).

This script assumes the unique sample name is in the first string before the underscore delimiter. If this isn't the case then change lines 38 and 49 accordingly. For example the following will take the first three strings based on the underscore delimiter:

`awk -F '_' '{print $1"_"$2"_"$3}'`

This script will place everything in the current working directory and then submit the results to `02-align_cluster.sh`, so make sure that script is placed in the project directory or remove lines 70-77 if you want to use the node version.

## Alignment Pipeline

The primary pipeline is `02-align.sh`, where there are both cluster (`02-align_cluster.sh`) and node (`02-align.sh`) versions. This script will [trim](https://github.com/alexdobin/STAR/issues/455#issuecomment-407539412), align, and generate a generate a [gene count matrix](https://www.biostars.org/p/218995/) for each sample. You will have to edit `genomes` to be the genome directory on your cluster as well as either edit or remove `export PYTHON_EGG_CACHE`. Both scripts assume 2color chemistry (NovaSeq and NextSeq). If your data is from the HiSeq and MiSeq, the you should replace `--2colour` with `--quality`.

This pipeline requires a specific directory structure, which is made by `01-getFASTQ.sh`; however you can also create it yourself if you did not use that script:

1.	Create a parent directory for the project
2.	Within that parent project directory, add a text file called “task_samples.txt”, where each new line contains the unique part of sample name exactly as it appears on the fastq read pair files, aside from the file extension (“.fastq.gz”). Only list a sample once.

Overall, the directory tree structure should be the following:

```
├── Project
│   ├── raw_sequences
│   │   ├── sample1_R1_001.fastq.gz
│   │   ├── sample1_R2_001.fastq.gz
│   │   ├── sample2_R1_001.fastq.gz
│   │   ├── sample2_R2_001.fastq.gz
│   ├── task_samples.txt
```

### Cluster 

The cluster version of script can be executed using the following command from the working directory: 

`sbatch --array=1-96 02-align_cluster.sh`.

The `--array` argument represents which samples in the list to run. Here we are running samples 1 to 96. You could run select samples using the following format --array=2,4-12.

### Node

The node version of the script can be run from the working directory, where you specify the number of parallel runs through the `-j` argument on lines 107 and 108 in the parallel calls. Each parallel process requires 8 cores and 32 GB of ram. 

## Quality Control

The `03-QC.sh` script is the final step in the process, which a generates a html report and places the gene count matrices into a single folder. It requires `03-multiqc_config.yaml` to be placed in the same folder and you can edit this file to contain your project information.

## Differential Gene Expression Analysis

The `04-limma-voom.R` script provides an example of a differential gene expression analysis pipeline in R.
