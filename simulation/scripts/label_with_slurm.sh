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

# Default parameters
[ -z $NPROCESS ] && NPROCESS=7
[ -z $MINL ] && MINL=100
[ -z $SUBSIZE ] && SUBSIZE=5
[ -z $TEMP3_DIR ] && TEMP3_DIR=/zata/zippy/boxu/for_hzr/Software/LOCATE-paper-custom-scripts/TEMP3
[ -z $BLACKLIST ] && BLACKLIST=""
[ -z $GERM_MODEL ] && GERM_MODEL=""
[ -z $SOMA_MODEL ] && SOMA_MODEL=""

# Get contigs list from reference genome
[ -z $REF_TEMPLATE ] && echo "ERROR: REF_TEMPLATE not set" && exit 1
CONTIGS=(`cut -f 1 ${REF_TEMPLATE}.fai`)

# Program path
if [ -z $LABEL_DIR ]; then
    LABEL_DIR="/data/tusers/zhongrenhu/for_SMS/bin/demo3/test_with_simulation"
fi

###############################
### Get files from FILELIST ###
###############################
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
toks=($str)
BAM=${toks[0]}
RMSK=${toks[1]}
GAP=${toks[2]}
TE=${toks[3]}
WORKDIR=`dirname $BAM`
OUTPUTDIR=${WORKDIR}/result_no_secondary
WORKDIR=`dirname $WORKDIR`

echo "Parameters:"
echo "  BAM: $BAM"
echo "  RMSK: $RMSK"
echo "  GAP: $GAP"
echo "  TE: $TE"
echo "  WORKDIR: $WORKDIR"
echo "  OUTPUTDIR: $OUTPUTDIR"
echo ""

###########################################
### Create temporary directory for work ###
###########################################
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t label-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)

### Copy inputs to temp dir ###
echo "Copying input file ${BAM} to temp directory..."
for contig in ${CONTIGS[*]}; do
    cp -r $WORKDIR/$contig .
done
mkdir -p $OUTPUTDIR
cp -r $OUTPUTDIR . 2>/dev/null || mkdir -p result_no_secondary
cd result_no_secondary && echo "Changing wd: " $(pwd)
echo "Done."
echo ""

##################################
### Activate conda environment ###
##################################
source $CONDA_PATH
conda activate $CONDA_ENV

########################
### Process the data ###
########################
echo "Labeling for ${BAM}"
bash $LABEL_DIR/label_protocol_local.sh \
    -b $BAM \
    -d $TMPDIR \
    -r $RMSK \
    -g $GAP \
    -P $NPROCESS \
    -T $NPROCESS \
    --refTe $TE \
    --minLen $MINL \
    --subSize $SUBSIZE \
    --TEMP3 $TEMP3_DIR \
    --blacklist $BLACKLIST \
    --germModel $GERM_MODEL \
    --somaModel $SOMA_MODEL
echo ""

################################
### Copy files to output dir ###
################################
echo "Copying results to destination..."
ls -lh
cp AllIns* $OUTPUTDIR
cp *bed $OUTPUTDIR
cp *txt $OUTPUTDIR
cp merge* $OUTPUTDIR
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