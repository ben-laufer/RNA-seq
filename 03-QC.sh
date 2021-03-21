#!/bin/bash
#
#SBATCH --job-name=RNA_QC
#SBATCH --ntasks=1 # Number of cores/threads
#SBATCH --mem=2000 # Ram in Mb
#SBATCH --partition=production 
#SBATCH --time=0-00:30:00

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
##########################################################################################

###################
# Run Information #
###################

start=`date +%s`

hostname

################
# Load Modules #
################

module load multiqc/1.9

###########
# MultiQC #
###########

call="multiqc
. \
--ignore alignLogs/ \
--ignore raw_sequences/ \
 --config 03-multiqc_config.yaml"

echo $call
eval $call

########
# Copy #
########

mkdir GeneCounts
"$(find `.` -name '*ReadsPerGene.out.tab' -print0 | xargs -0 cp -t GeneCounts)"

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
