#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=30G
#SBATCH -c 24
#SBATCH --array=1-900%20
#SBATCH --partition=12hours
#SBATCH --output=./logs/label-log-%A-%a.out


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
[ -z $FILELIST ] && FILELIST=your_path_to/simulation/simulate_training_data/dm6/filelist_examples/label_filelist_dm6
[ -z $CONTIGS ] && CONTIGS=(`cut -f 1 your_path_to/Dm6_no_alt.fa.fai`)

REPEAT='your_path_to/Dm6.rmsk.bed'
GAP='your_path_to/Dm6.gap.bed'
TE='your_path_to/Dm6.transposon.fa'
TEMP3='your_path_to/TEMP3'
BLACKLIST=''
germModel=''
somaModel=''

### Get parameters from FILELIST ###
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
toks=($str)
BAM=${toks[0]}
WORKDIR=`dirname $BAM`
OUTPUTDIR=${WORKDIR}/result_no_secondary
WORKDIR=`dirname $WORKDIR`


###########################################
### Create temporary directory for work ###
###########################################
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t zhongrenhu-tmp-XXXXXXXX)
cd $TMPDIR

### Copy inputs to temp dir ###
echo "Copying input file ${BAM} to temp directory..."
for contig in ${CONTIGS[*]}; do
    cp -r $WORKDIR/$contig .
done
cp -r $OUTPUTDIR .
cd result_no_secondary && echo "Changing wd: " $(pwd)
echo "Done."
echo ""


##################################
### Activate conda environment ###
##################################
source your_path_to/conda.sh
conda activate TEMP3


########################
### Process the data ###
########################
echo "Labeling for ${BAM}"
bash your_path_to/simulation/label/label_protocol.sh \
    -b $BAM \
    -d $WORKDIR \
    -r $REPEAT \
    -g $GAP \
    -P 10 \
    -T 5 \
    --refTe $TE \
    --blacklist $BLACKLIST \
    --minLen 100 \
    --TEMP3 $TEMP3 \
    --subSize 5 \
    --germModel $germModel \
    --somaModel $somaModel
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
cd $HOME
echo "Deleting temp dir: " $TMPDIR
rm -rd $TMPDIR
echo ""
echo "Script complete."
date
