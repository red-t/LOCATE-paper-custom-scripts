#!/bin/bash

path1="${path1:-path1}"
path2="${path2:-path2}"

start_time=$(date +%s)

DATA_PATH="${path2}/rawdata/simulation/IRGSPv1/new_ngs"
OUT_PATH_dir="${path2}/result/simulation/IRGSPv1_new"

ANNO_PATH="${path1}/" #"${path1}"
GENOME="IRGSPv1"
CPU="10"

TE="IRGSPv1.transposon" # ALUL1SVA

BWA_INDEX="${ANNO_PATH}/${GENOME}/BWAIndex"


GENOME_FA="${ANNO_PATH}/${GENOME}/${GENOME}.fa"
TE_ANNO_FA="${ANNO_PATH}/${GENOME}/${GENOME}.transposon.fa"
REPEATMASKER_FILE="${ANNO_PATH}/${GENOME}/${GENOME}.rmsk.bed"
REPEAT_CLASS="${ANNO_PATH}/${GENOME}/${GENOME}.transposon.class"





[ ! -d ${OUT_PATH_dir} ] && mkdir ${OUT_PATH_dir}

for SAMPLE in simulation_germ  # _clr simulation_germ_ccs simulation_soma_ont simulation_soma_clr  simulation_soma_ccs # simulation_soma # HG001 # HG002 # SRR9685183.sra
do
    for dep in 50 # 1 2 3 4 5 10 20 30 40 50
    do

        # break
        start_time=$(date +%s)
        OUT_PATH=${OUT_PATH_dir}/${dep}
        
        [ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}

        [ ! -d ${OUT_PATH}/melt ] && mkdir ${OUT_PATH}/melt
        [ ! -d ${OUT_PATH}/melt/${SAMPLE} ] && mkdir ${OUT_PATH}/melt/${SAMPLE}
        [ ! -d ${OUT_PATH}/TEMP2 ] && mkdir ${OUT_PATH}/TEMP2
        [ ! -d ${OUT_PATH}/TEMP2/${SAMPLE} ] && mkdir ${OUT_PATH}/TEMP2/${SAMPLE}

        # NGS mapping
        start_map=$(date +%s)
        
        LEFT=${DATA_PATH}/${SAMPLE}.${dep}X_R1.fastq
        RIGHT=${DATA_PATH}/${SAMPLE}.${dep}X_R2.fastq

        # if [ ! -f ${DATA_PATH}/${SAMPLE}_ngs.${dep}X.sam ];then
        #     bwa mem -t ${CPU} ${BWA_INDEX}/genome ${LEFT} ${RIGHT} > ${DATA_PATH}/${SAMPLE}_ngs.${dep}X.sam 2>/dev/null
        # fi

        # if [ ! -f ${DATA_PATH}/${SAMPLE}_ngs.${dep}X.sort.bam ];then
        #     samtools view -@ ${CPU} -h ${DATA_PATH}/${SAMPLE}_ngs.${dep}X.sam | sed '/^@HD/d' | samtools sort -@ ${CPU} -O bam -o ${DATA_PATH}/${SAMPLE}_ngs.${dep}X.sort.bam - 
        #     samtools index -@ ${CPU} ${DATA_PATH}/${SAMPLE}_ngs.${dep}X.sort.bam
        # fi
        
        end_map=$(date +%s)
        time_map=$(( $end_map - $start_map ))
        echo "time_map>>>>${time_map}"

        ngs_bam_file=${DATA_PATH}/${SAMPLE}_ngs.${dep}X.sort.bam

        # start_temp2=$(date +%s)

        # if [ ! -f ${OUT_PATH}/TEMP2/${SAMPLE}_${dep}X.insertion.bed ];then
        #     source ${path2}/miniconda3/envs/O_O/etc/profile.d/conda.sh
        #     conda activate py2 && TEMP2 insertion -l ${LEFT} -r ${RIGHT} -i ${ngs_bam_file} -g ${GENOME_FA} -R ${TE_ANNO_FA} -t ${REPEATMASKER_FILE} -o ${OUT_PATH}/TEMP2/${SAMPLE} -p ${SAMPLE}_${dep}X -d -c ${CPU} && conda deactivate
        # fi
        # end_temp2=$(date +%s)
        # time_temp2=$(( $end_temp2 - $start_temp2 ))


        ## map for MELT
        # if [ ! -f ${DATA_PATH}/for_melt/${SAMPLE}_ngs.${dep}X.sam ];then
        #     bwa mem -t ${CPU} ${BWA_INDEX}/genome_melt ${LEFT} ${RIGHT} > ${DATA_PATH}/for_melt/${SAMPLE}_ngs.${dep}X.sam 2>/dev/null
        # fi

        # if [ ! -f ${DATA_PATH}/for_melt/${SAMPLE}_ngs.${dep}X.sort.bam ];then
        #     samtools view -@ ${CPU} -h ${DATA_PATH}/for_melt/${SAMPLE}_ngs.${dep}X.sam | sed '/^@HD/d' | samtools sort -@ ${CPU} -O bam -o ${DATA_PATH}/for_melt/${SAMPLE}_ngs.${dep}X.sort.bam - 
        #     samtools index -@ ${CPU} ${DATA_PATH}/for_melt/${SAMPLE}_ngs.${dep}X.sort.bam
        # fi


        # MELT
        start_melt=$(date +%s)
        MELT_PATH="${path2}/tools/MELTv2.2.2"
        rm -r ${OUT_PATH}/melt/${SAMPLE}/*
        java -jar ${MELT_PATH}/MELT.jar Single \
                -a \
                -h ${ANNO_PATH}/${GENOME}/${GENOME}.fa \
                -bamfile ${DATA_PATH}/${SAMPLE}_ngs.${dep}X.sort.bam \
                -n ${ANNO_PATH}/${GENOME}/IRGSPv1.genes.bed \
                -t ${path1}/IRGSPv1/meltTEZip/mei_list.txt \
                -w ${OUT_PATH}/melt/${SAMPLE}/ \
                -d 100 \
                -nocleanup \
                -z 6000

        while read info; do
            te_info=(${info})
            melt_te=${te_info[0]%%_*}
            if [ -f ${OUT_PATH}/melt/${SAMPLE}/${melt_te}.master.bed ];then
                awk -v te=${melt_te} 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,te,$5,$6}' ${OUT_PATH}/melt/${SAMPLE}/${melt_te}.master.bed >> ${OUT_PATH}/melt/${SAMPLE}/${SAMPLE}_melt.bed
            fi
        done <<< "$( cat ${TRANSPOSON_SIZE} )"

        end_melt=$(date +%s)
        time_melt=$(( $end_melt - $start_melt ))

        echo -e "\n-----------\n" >> ${OUT_PATH_dir}/time_simulation_ngs_rice.melt.txt
        echo -e $(date) >> ${OUT_PATH_dir}/time_simulation_ngs_rice.melt.txt
        echo -e "map time\t"${time_map} >> ${OUT_PATH_dir}/time_simulation_ngs_rice.melt.txt
        echo -e "sample\tdepth\tTEMP2\tmelt" >> ${OUT_PATH_dir}/time_simulation_ngs_rice.melt.txt
        echo -e "${SAMPLE}\t${dep}\t${time_temp2}\t${time_melt}" >> ${OUT_PATH_dir}/time_simulation_ngs_rice.melt.txt


    done
done


