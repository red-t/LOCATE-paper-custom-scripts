#!/bin/bash

path1="${path1:-path1}"
path2="${path2:-path2}"

start_time=$(date +%s)

DATA_PATH="${path2}/rawdata/simulation"
OUT_PATH_dir="${path2}/result/simulation/dm6_new_ont_long"



ANNO_PATH="${path1}" #"${path1}"
GENOME="dm6"
CPU="10"

BWA_INDEX="${ANNO_PATH}/${GENOME}/BWAIndex"

TE="dm6.transposon" # ALUL1SVA

TRANSPOSON_SIZE="${ANNO_PATH}/${GENOME}/${GENOME}.transposon.size"
GENOME_FA="${ANNO_PATH}/${GENOME}/benchmark/${GENOME}_template.fa"
TE_ANNO_FA="${ANNO_PATH}/${GENOME}/benchmark/${GENOME}_clean.transposon.fa"
TE_ANNO_FA2="${ANNO_PATH}/${GENOME}/benchmark/${GENOME}_clean.transposon.2.fa"

REPEATMASKER_FILE="${ANNO_PATH}/${GENOME}/benchmark/${GENOME}.repeat.bed"
REPEATMASKER_FILE="${path1}/dm6_clean/dm6_clean.rmsk.bed"
REPEAT_CLASS="${ANNO_PATH}/${GENOME}/benchmark/${GENOME}_clean.transposon.class"
CHROM_LIST="${ANNO_PATH}/run_dm6.chrom.list"

library=${TE_ANNO_FA}



GraffiTE_PATH="${path2}/tools/GraffiTE"


for SAMPLE in simulation_germ_ont # simulation_germ_ont simulation_germ_clr # simulation_germ_clr simulation_germ_ont # simulation_germ_ccs simulation_germ_clr simulation_germ_ont  # simulation_germ_clr ## simulation_soma_ont simulation_soma_clr simulation_soma_ccs #   # simulation_soma # HG001 # HG002 # SRR9685183.sra
do
    for depth in 1 2 3 4 5 10 20 30 40 50 # 1 2 3 4 5 10 20 30 40 50
    do
        OUT_PATH=${OUT_PATH_dir}/${depth}/
        [ ! -d ${OUT_PATH}/GraffiTE/${SAMPLE} ] && mkdir -p ${OUT_PATH}/GraffiTE/${SAMPLE}

        if [[ ${SAMPLE}  =~ "clr" ]];then 
            template_map_mode="pb" 
        fi
        if [[ ${SAMPLE}  =~ "ont" ]];then 
            template_map_mode="ont" 
        fi
        if [[ ${SAMPLE}  =~ "ccs" ]];then 
            template_map_mode="hifi" 
        fi
        
        # if [ ! -f ${DATA_PATH}/${GENOME}/tgs/${SAMPLE}.${depth}X.fa ];then
        #     cp ${path2}/rawdata/simulation/dm6/tgs/${SAMPLE}.${depth}X.fa ${DATA_PATH}/${GENOME}/tgs/
        # fi

        echo "path,sample,type" > ${OUT_PATH}/GraffiTE/${SAMPLE}/longreads.csv
        echo "${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${depth}X.fa,${SAMPLE},${template_map_mode}" >> ${OUT_PATH}/GraffiTE/${SAMPLE}/longreads.csv
        cp ${OUT_PATH}/GraffiTE/${SAMPLE}/longreads.csv ${OUT_PATH}/GraffiTE/${SAMPLE}/reads.csv

        nextflow run ${GraffiTE_PATH}/main.nf \
            --cores 20 \
            --longreads ${OUT_PATH}/GraffiTE/${SAMPLE}/longreads.csv \
            --TE_library ${library} \
            --reference ${GENOME_FA} \
            --reads ${OUT_PATH}/GraffiTE/${SAMPLE}/reads.csv \
            --out ${OUT_PATH}/GraffiTE/${SAMPLE}/ \
            --tsd_time 4h \
            --map_asm_time 4h \
            --svim_asm_time 4h \
            --map_longreads_time 4h \
            --repeatmasker_time 4h \
            --genotype true [-with-singularity ${GraffiTE_PATH}/graffite_latest.sif]

        # tar -zcf ${OUT_PATH}/GraffiTE/${SAMPLE}/work.tar.gz -C ${OUT_PATH}/GraffiTE/${SAMPLE}/work . --remove-files
        # rm ${DATA_PATH}/${GENOME}/tgs/${SAMPLE}.${depth}X.fa

    done
done


