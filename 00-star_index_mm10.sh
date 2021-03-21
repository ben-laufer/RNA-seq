#!/bin/bash
#
#SBATCH --job-name=star_index
#SBATCH --workdir /share/lasallelab/genomes
#SBATCH --ntasks=30 # Number of cores/threads
#SBATCH --mem=64000 # Ram in Mb
#SBATCH --partition=production 
#SBATCH --time=0-01:00:00

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
##########################################################################################

###################
# Run Information #
###################

start=`date +%s`

hostname

THREADS=${SLURM_NTASKS}
MEM=$((${SLURM_MEM_PER_CPU}/1024))

echo "Allocated threads: ${THREADS}"
echo "Allocated memory:  ${MEM}"

################
# Load Modules #
################

module load star/2.7.3a

###########################
# Download Genome (FASTA) #
###########################
# Must have un-placed/un-localized scaffolds, so use .dna.primary.assembly
# Section 2.2.1 https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

export mainPath="/share/lasallelab"

# UCSC
#mkdir -p ${mainPath}/genomes/mm10
#cd ${mainPath}/genomes/mm10
#rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz .
#gunzip mm10.fa.gz

# Ensembl
mkdir -p ${mainPath}/genomes/GRCm38
cd ${mainPath}/genomes/GRCm38
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz .
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

########################
# Download Genes (GTF) #
########################

# UCSC
#rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz .
#gunzip mm10.refGene.gtf.gz

# Ensembl
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz .
gunzip Mus_musculus.GRCm38.102.gtf.gz

####################
# Build Star Index #
####################
# the splice-junction-data-base-overhang parameter should have a value of read length – 1

mkdir -p star_150/

call="STAR \
--runThreadN 30 \
--runMode genomeGenerate \
--genomeDir star_150/ \
--genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa \
--sjdbGTFfile Mus_musculus.GRCm38.102.gtf \
--sjdbOverhang 149"

echo $call
eval $call

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
