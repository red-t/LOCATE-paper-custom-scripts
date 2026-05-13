#!/bin/bash

#########################################
### Print a little info for debugging ###
#########################################
echo "HOSTNAME: " $(hostname)
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
date
echo ""

################################
### Parameter initialization ###
################################
[ -z $FILELIST ] && echo "ERROR: FILELIST not set" && exit 1

# Get contigs list from reference genome
[ -z $REF_TEMPLATE ] && echo "ERROR: REF_TEMPLATE not set" && exit 1
CONTIGS=(`cut -f 1 ${REF_TEMPLATE}.fai`)

###############################
### Get files from FILELIST ###
###############################
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
toks=($str)
SAMPLE=${toks[0]}

echo "Sample: $SAMPLE"

###########################################
### Create temporary directory for work ###
###########################################
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t merge-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)
echo ""

###############################
### Copy inputs to temp dir ###
###############################
echo "Copying input files from ${SAMPLE}..."
for contig in ${CONTIGS[*]}
do
    mkdir ${contig} && cp ${SAMPLE}/${contig}/TGS.fasta ${contig}
done
echo "Done."
echo ""

###############
### Process ###
###############
echo "Merging..."
cat */TGS.fasta >> TGS.fasta
echo "Done." && echo ""

#################################
### Copy files to output dir ###
################################
echo "Copying results to destination..."
ls -lh
cp TGS.fasta $SAMPLE
echo "Done."
echo ""

################
### Clean up ###
################
echo "Cleaning up..."
cd /tmp
echo "Deleting temp dir: " $TMPDIR
rm -rd $TMPDIR
echo ""
echo "Script complete."
date