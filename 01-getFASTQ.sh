#!/bin/bash

# Arguments
project=RNA-seq
#https://bioinformatics.ucdavis.edu/research-computing/documentation/archiving-slims-data/
SLIMSstring=gz3bl6zibe
SLIMSdir=Un_DTDB126/Project_BL_Nova150P_Laufer

echo "Creating directory for $project"
mkdir ${project}

echo "Downloading fastq files for $SLIMSstring"
call="rsync \
-avL \
slimsdata.genomecenter.ucdavis.edu::slims/${SLIMSstring}/ \
${project}"
	
echo ${call}
eval ${call}

echo "md5sum check"
cd ${project}/${SLIMSdir}
if md5sum -c \@md5Sum.md5
then
    echo
    echo "All files have the correct md5sum"
else
    echo "ERROR: Some files are corrupt or missing"
    exit 1
fi

echo "Moving undetermined files"
mkdir Other
mv Undetermined* Other

echo "Checking for the right number of unique sample IDs for both R1 and R2"
countFASTQ(){
	awk -F '_' '{print $1"_"$2}' | \
	sort -u | \
	wc -l
}
export -f countFASTQ

R1=`ls -1 *R1*.gz | countFASTQ`
R2=`ls -1 *R2*.gz | countFASTQ`

echo "Creating a file of unique IDs based on first two strings from underscore delimiter"
ls -1 *fastq.gz | \
awk -F '_' '{print $1"_"$2}' | \
sort -u > \
task_samples.txt

echo "Creating final directory"
mkdir ../../raw_sequences
mv task_samples.txt ../..
mv *.fastq.gz ../../raw_sequences

echo "Checking read pairs"
cd ../../raw_sequences
pairedReads=$(($(ls | wc -l)/2))
if [ ${pairedReads} = ${R1} ]
then
        echo "${pairedReads} will be aligned using STAR"
else
        echo "ERROR: Incorrect number of read pairs in raw_sequences"
        echo "There are ${pairedReads} but there should be ${R1}"
        exit 1
fi

echo "Submitting $genome alignment command to cluster for ${project}"
cd ..
call="sbatch \
--array=1-${pairedReads} \
02-align.sh"
	
echo ${call}
eval ${call}

echo "Done"
exit 0
