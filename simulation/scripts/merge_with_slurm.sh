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
    mkdir ${contig}
    if [ -f ${SAMPLE}/${contig}/TGS.fastq ]; then
        cp ${SAMPLE}/${contig}/TGS.fastq ${contig}
    elif [ -f ${SAMPLE}/${contig}/TGS.fasta ]; then
        cp ${SAMPLE}/${contig}/TGS.fasta ${contig}
    fi
    if [ -f ${SAMPLE}/${contig}/NGS_1.fastq ]; then
        cp ${SAMPLE}/${contig}/NGS_1.fastq ${SAMPLE}/${contig}/NGS_2.fastq ${contig}/
    fi
done
echo "Done."
echo ""

###############
### Process ###
###############
echo "Merging..."
if ls */TGS.fastq 1>/dev/null 2>&1; then
    cat */TGS.fastq > TGS.fastq
    echo "Compressing TGS.fastq..."
    gzip -1 TGS.fastq
elif ls */TGS.fasta 1>/dev/null 2>&1; then
    cat */TGS.fasta > TGS.fasta
fi
if ls */NGS_1.fastq 1>/dev/null 2>&1; then
    cat */NGS_1.fastq > NGS_1.fastq
    cat */NGS_2.fastq > NGS_2.fastq
    echo "Compressing NGS PE files..."
    gzip -1 NGS_1.fastq NGS_2.fastq
fi
echo "Done." && echo ""

#################################
### Copy files to output dir ###
################################
echo "Copying results to destination..."
ls -lh
cp TGS.fastq.gz TGS.fasta $SAMPLE 2>/dev/null
cp NGS_1.fastq.gz NGS_2.fastq.gz $SAMPLE 2>/dev/null
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