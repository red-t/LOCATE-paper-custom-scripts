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

# Track overall start time
SCRIPT_START=$(date +%s)

################################
### Parameter initialization ###
################################
[ -z $FILELIST ] && echo "ERROR: FILELIST not set" && exit 1
[ -z $CONDA_PATH ] && CONDA_PATH=/zata/zippy/boxu/for_hzr/Software/miniforge3/etc/profile.d/conda.sh
[ -z $CONDA_ENV ] && CONDA_ENV=simulation

# Python unbuffered output for real-time log visibility
export PYTHONUNBUFFERED=1

# Default parameters
[ -z $NGS_LEN ] && NGS_LEN=150
[ -z $NGS_INNER ] && NGS_INNER=200
[ -z $NGS_STD ] && NGS_STD=20
[ -z $NGS_ERR ] && NGS_ERR=0.0005
[ -z $TGS_MAXL ] && TGS_MAXL=300000
[ -z $TGS_MINL ] && TGS_MINL=100
[ -z $SIMULATOR ] && SIMULATOR=builtin
[ -z $PBSIM_MODEL ] && PBSIM_MODEL=""

# Program path (auto-detect or use SIMULATION_DIR)
if [ -z $SIMULATION_DIR ]; then
    SIMULATION_DIR="/zata/zippy/boxu/for_hzr/Software/LOCATE-paper-custom-scripts/simulation"
fi

### Get parameters from FILELIST ###
echo "Reading filelist..."
str=`sed "${SLURM_ARRAY_TASK_ID}q;d" $FILELIST`
toks=($str)
WORKDIR=${toks[0]} && CONTIG=`basename ${WORKDIR}`
TE_FA=${toks[1]}
GENOME_SIZE=${toks[2]}
CONTIG_SIZE=${toks[3]}
POP_SIZE=${toks[4]}
SUB_POP_SIZE=${toks[5]}
DEPTH=${toks[6]}
PROTOCOL=${toks[7]}
[ $PROTOCOL == "ccs" ] && TGS_MEANL=13490
[ $PROTOCOL == "clr" ] && TGS_MEANL=7896
[ $PROTOCOL == "ont" ] && TGS_MEANL=7170

echo "Parameters:"
echo "  WORKDIR: $WORKDIR"
echo "  CONTIG: $CONTIG"
echo "  TE_FA: $TE_FA"
echo "  GENOME_SIZE: $GENOME_SIZE"
echo "  CONTIG_SIZE: $CONTIG_SIZE"
echo "  POP_SIZE: $POP_SIZE"
echo "  SUB_POP_SIZE: $SUB_POP_SIZE"
echo "  DEPTH: $DEPTH"
echo "  PROTOCOL: $PROTOCOL"
echo ""

###########################################
### Create temporary directory for work ###
###########################################
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t simulation-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)

### Copy inputs to temp dir ###
echo "Copying ${WORKDIR} to temp directory..."
cp -r $WORKDIR .
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
cd $CONTIG
N_SUB=`awk -v sub_pop=$SUB_POP_SIZE -v pop=$POP_SIZE 'BEGIN{n_sub=int(pop/sub_pop); print n_sub}'`
NGS_READS=`awk -v depth=$DEPTH -v g_l=$GENOME_SIZE -v c_l=$CONTIG_SIZE -v n_sub=$N_SUB -v r_l=150 -v inner=$NGS_INNER 'BEGIN{ratio=c_l/g_l; ngs_reads=int(ratio*(depth*g_l/(2*r_l+inner))/n_sub); print ngs_reads}'`
TGS_READS=`awk -v depth=$DEPTH -v g_l=$GENOME_SIZE -v c_l=$CONTIG_SIZE -v n_sub=$N_SUB -v r_l=$TGS_MEANL 'BEGIN{ratio=c_l/g_l; tgs_reads=int(ratio*(depth*g_l/r_l)/n_sub); print tgs_reads}'`

echo "Processing config: N_SUB=$N_SUB TGS_READS=$TGS_READS PROTOCOL=$PROTOCOL"
date

if [ -f $CONTIG.0.pgd ]; then
    for((j=0; j<$N_SUB; j++))
    do
        ITER_START=$(date +%s)
        echo ""
        echo "=========================================="
        echo "[ $(date) ] Iteration $j / $((N_SUB-1)) — Build sub population genome of $CONTIG"
        echo "=========================================="

        # build population genome
        if [ ! -f $CONTIG.$j.fa ]; then
            python $SIMULATION_DIR/build-population-genome.py \
                --pgd $CONTIG.$j.pgd \
                --te-seqs $TE_FA \
                --chassis $CONTIG.tmp.chasis.fasta \
                --output $CONTIG.$j.fa \
                --ins-seq "" \
                --sub_idx $j \
                --sub_size $SUB_POP_SIZE
        fi

        # generate TGS
        echo "[ $(date) ] Iteration $j — Generate TGS data ($PROTOCOL, simulator=$SIMULATOR)"
        if [ "$SIMULATOR" == "pbsim" ]; then
            if [ ! -f $CONTIG.${j}_tgs.fastq ]; then
                if [ -z "$PBSIM_MODEL" ]; then
                    echo "ERROR: SIMULATOR=pbsim but PBSIM_MODEL not set" && exit 1
                fi
                PBSIM_DEPTH=$(awk -v depth=$DEPTH -v pop=$POP_SIZE 'BEGIN{printf "%.2f", depth / pop}')
                pbsim --strategy wgs \
                      --method errhmm \
                      --errhmm $PBSIM_MODEL \
                      --genome $CONTIG.$j.fa \
                      --depth $PBSIM_DEPTH \
                      --prefix $CONTIG.${j}_tgs
                # errhmm 输出 gzipped fq, 解压合并后重命名 read 避免多 contig/j 间冲突
                zcat $CONTIG.${j}_tgs_*.fq.gz | \
                    awk -v prefix="${CONTIG}_${j}_" 'NR%4==1{print "@" prefix substr($0,2); next} {print}' \
                    > $CONTIG.${j}_tgs.fastq
                rm -f $CONTIG.${j}_tgs_*.fq.gz $CONTIG.${j}_tgs_*.ref $CONTIG.${j}_tgs_*.maf
            fi
        else
            if [ ! -f $CONTIG.${j}_tgs.fasta ]; then
                python $SIMULATION_DIR/generate_TGS.py \
                    --pg $CONTIG.$j.fa \
                    --reads $TGS_READS \
                    --fasta $CONTIG.${j}_tgs.fasta \
                    --tgs-maxl $TGS_MAXL \
                    --tgs-minl $TGS_MINL \
                    --protocol $PROTOCOL
            fi
        fi

        # remove intermediate sub-population genome
        rm $CONTIG.$j.fa

        ITER_END=$(date +%s)
        echo "[ $(date) ] Iteration $j done (took $((ITER_END - ITER_START)) seconds)"
    done
else
    echo "WARNING: $CONTIG.0.pgd not found — skipping all processing"
fi

cd ..
if [ -f $CONTIG/$CONTIG.0_tgs.fastq ]; then
    cat */*_tgs.fastq > TGS.fastq
elif [ -f $CONTIG/$CONTIG.0_tgs.fasta ]; then
    cat */*_tgs.fasta > TGS.fasta
fi

################################
### Copy files to output dir ###
################################
echo "Copying results to destination..."
ls -lh
cp -r TGS* $WORKDIR
echo "Done."
echo ""

################
### Clean up ###
################
echo "Cleaning up..."
cd /tmp
echo "Deleting temp dir: " $TMPDIR
rm -rd $TMPDIR
SCRIPT_END=$(date +%s)
echo ""
echo "Script complete."
echo "Total elapsed: $((SCRIPT_END - SCRIPT_START)) seconds"
date