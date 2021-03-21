align(){

	######################
	# Set Up Environment #
	######################

	directory=${PWD}/
	sample=$1
	rawpath=${directory}raw_sequences/
	mappath=${directory}${sample}
	fastq1=${rawpath}${sample}_*_R1_001.fastq.gz
	fastq2=${rawpath}${sample}_*_R2_001.fastq.gz
	trim1=${sample}_*_val_1.fq.gz
	trim2=${sample}_*_val_2.fq.gz
	BAM=${sample}_Aligned.sortedByCoord.out.bam

	########
	# Trim #
	########
	# Use 2color for NovaSeq and NextSeq, replace with quality for HiSeq and MiSeq
	# Should trim for STAR: https://github.com/alexdobin/STAR/issues/455#issuecomment-407539412

	mkdir ${mappath}

	call="trim_galore \
	--paired \
	--cores 2 \
	--2colour 20 \
	--fastqc \
	--output_dir ${mappath} \
	${fastq1} \
	${fastq2}" 

	echo $call
	eval $call

	#########
	# Align #
	#########
	# adjust threads and genome directory
	# Use zcat command for fastq.gz https://www.biostars.org/p/243683/
	# ENCODE options from section 3.3.2 of STAR manual
	# Use qauntMode to get GeneCounts for R https://www.biostars.org/p/218995/

	cd ${mappath}

	call="STAR \
	--runThreadN 8 \
	--genomeDir /share/lasallelab/genomes/GRCm38/star_150/
	--readFilesIn ${trim1} ${trim2} \
	--readFilesCommand zcat \
	--outFilterType BySJout \
	--outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 
	--outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ${sample}_ \
	--quantMode GeneCounts"

	echo $call
	eval $call

	#########
	# Index #
	#########

	call="samtools \
	index \
	-@ 7 \
	${BAM}"

	echo $call
	eval $call

}
export -f align

################
# Load Modules #
################

module load star/2.7.3a
module load samtools/1.11
export PYTHON_EGG_CACHE="${mainPath}/programs/CpG_Me"
module load trim_galore/0.6.6
source activate cutadapt-2.10
export mainPath="/share/lasallelab"

#######
# Run #
#######

cd /share/lasallelab/Ben/PEBBLES/RNA/

mkdir alignLogs

parallel --dry-run --will-cite --results alignLogs -j 4 "align {}" :::: task_samples.txt
parallel --verbose --will-cite --results alignLogs -j 4 "align {}" :::: task_samples.txt



