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
[ -z $CONDA_PATH ] && CONDA_PATH=/zata/zippy/zhongrenhu/Software/mambaforge/etc/profile.d/conda.sh
[ -z $CONDA_ENV ] && CONDA_ENV=TEMP3

[ -z $THREADS ] && THREADS=${SLURM_THREADS:-30}

###############################
### Get files from FILELIST ###
###############################
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
toks=($str)
REF=${toks[0]}
QUERY=${toks[1]}
PRESET=${toks[2]}

echo "Parameters:"
echo "  REF: $REF"
echo "  QUERY: $QUERY"
echo "  PRESET: $PRESET"
echo ""

###########################################
### Create temporary directory for work ###
###########################################
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t minimap2-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)
echo ""

###############################
### Copy inputs to temp dir ###
###############################
echo "Copying input file ${QUERY} to temp directory..."
cp $REF .
cp $QUERY .
echo "Done."
echo ""

####################################################
### Get filenames in temp by cutting of the path ###
####################################################
PREFIX=`basename ${QUERY%.f*a*}` && QUERY_FA=`basename $QUERY` && REF_FA=`basename $REF`
OUTPUTDIR=`dirname $QUERY`

##################################
### Activate conda environment ###
##################################
source $CONDA_PATH
conda activate $CONDA_ENV

########################
### Process the data ###
########################
echo "Running minimap2 for ${PREFIX}"
minimap2 -aYx $PRESET --MD -t $THREADS $REF_FA $QUERY_FA | samtools view -@ $THREADS -bhS - | samtools sort -@ $THREADS -o $PREFIX.bam -
samtools index -@ $THREADS $PREFIX.bam
echo ""

#################################
### Copy files to output dir ###
################################
echo "Copying results to destination..."
ls -lh
mkdir ${OUTPUTDIR}/${PRESET}
cp -r $PREFIX.bam ${OUTPUTDIR}/${PRESET}
cp -r $PREFIX.bam.bai ${OUTPUTDIR}/${PRESET}
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