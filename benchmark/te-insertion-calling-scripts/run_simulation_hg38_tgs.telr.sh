#!/bin/bash

path1="${path1:-path1}"
path2="${path2:-path2}"

start_time=$(date +%s)

DATA_PATH="${path2}/rawdata/simulation/"
OUT_PATH_dir="${path2}/result/simulation/new_ont_long" # "${path2}/simulation/result/"

ANNO_PATH="${path1}" #"${path1}"
GENOME="hg38"
CPU="10"

TE="ALSE" # ALUL1SVA

TRANSPOSON_INDEX="${ANNO_PATH}/${GENOME}/Minimap2Index/${TE}.mmi"
TRANSPOSON_SIZE="${ANNO_PATH}/${GENOME}/${GENOME}.transposon.size"
GENOME_FA="${path1}/from_temp3/GRCh38_no_alt.fa" # "${ANNO_PATH}/${GENOME}/${GENOME}.fa"
GENOME_INDEX="${ANNO_PATH}/${GENOME}/Minimap2Index/${GENOME}.mmi"
TE_ANNO_FA="${ANNO_PATH}/from_temp3/GRCh38.transposon.fa" # "${ANNO_PATH}/${GENOME}/${TE}.fa"
TE_ANNO_FA2="${ANNO_PATH}/${GENOME}/${TE}.2.fa" # "${ANNO_PATH}/${GENOME}/${TE}.fa"
REPEATMASKER_FILE="${ANNO_PATH}/${GENOME}/${TE}.bed"
REPEAT_CLASS="${path1}/${GENOME}/hg38.repeat.class"

CHROM_LIST="${ANNO_PATH}/${GENOME}/run_hg38.chrom.list"




[ ! -d ${OUT_PATH_dir} ] && mkdir ${OUT_PATH_dir}

for SAMPLE in simulation_germ_ont  # simulation_germ_ont_0.02 simulation_soma_ont_0.02 # simulation_germ_clr simulation_germ_ccs simulation_germ_ont simulation_soma_ont simulation_soma_clr simulation_soma_ccs  #   # simulation_soma # HG001 # HG002 # SRR9685183.sra
do
    
    for dep in 1 2 3 4 5 10 20 30 40 50 # dm6
    do
        
        # break
        OUT_PATH=${OUT_PATH_dir}/${dep}
        [ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}

        [ ! -d ${OUT_PATH}/tldr ] && mkdir ${OUT_PATH}/tldr

        [ ! -d ${OUT_PATH}/TELR ] && mkdir ${OUT_PATH}/TELR
        [ ! -d ${OUT_PATH}/TELR/${SAMPLE} ] && mkdir ${OUT_PATH}/TELR/${SAMPLE}

        [ ! -d ${OUT_PATH}/TrEMOLO ] && mkdir ${OUT_PATH}/TrEMOLO
        [ ! -d ${OUT_PATH}/TrEMOLO/${SAMPLE} ] && mkdir ${OUT_PATH}/TrEMOLO/${SAMPLE}

        [ ! -d ${OUT_PATH}/LOCATE ] && mkdir ${OUT_PATH}/LOCATE
        [ ! -d ${OUT_PATH}/LOCATE/${SAMPLE} ] && mkdir ${OUT_PATH}/LOCATE/${SAMPLE}



        # TGS mapping
        echo ${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.bam

        bam_file=${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.bam # ${DATA_PATH}/${dep}/temp_bam/${SAMPLE}_tgs.${dep}X.bam
        
        # if [ ! -f ${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.fa ];then
        #     samtools fasta ${bam_file} -@ ${CPU} > ${DATA_PATH}/${GENOME}/tgs/${SAMPLE}.${dep}X.fa

        #     sed -i "s/:/_/g" ${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.fa
        #     sed -i "s/;/_/g" ${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.fa

        #     samtools faidx ${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.fa
        # fi



        # TELR
       
        if [ ! -f ${OUT_PATH}/TELR/${SAMPLE}/${SAMPLE}.${dep}X.telr.vcf ];then 
            start_telr=$(date +%s)
            source ${path2}/miniconda3/envs/O_O/etc/profile.d/conda.sh
            conda activate TELR && telr -i ${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.fa -r ${GENOME_FA} -l ${TE_ANNO_FA} -t 10 --aligner minimap2 -o ${OUT_PATH}/TELR/${SAMPLE} && conda deactivate
            end_telr=$(date +%s)
            time_telr=$(( $end_telr - $start_telr ))
        fi
        


            
        echo -e "\n-----------\n" >> ${OUT_PATH_dir}/time_simulation_tgs_new_ont_long.txt
        echo -e $(date) >> ${OUT_PATH_dir}/time_simulation_tgs_new_ont_long.txt
        echo -e "sample\tdepth\ttelr" >> ${OUT_PATH_dir}/time_simulation_tgs_new_ont_long.txt
        echo -e "${SAMPLE}\t${dep}\t${time_telr}" >> ${OUT_PATH_dir}/time_simulation_tgs_new_ont_long.txt

    done
    
done

# ${path2}/bin/run_intersection_simu_dm6.sh

