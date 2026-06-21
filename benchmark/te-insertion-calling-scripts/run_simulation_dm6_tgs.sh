#!/bin/bash

path1="${path1:-path1}"
path2="${path2:-path2}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

start_time=$(date +%s)

DATA_PATH="${path2}/rawdata/simulation/"
OUT_PATH_dir="${path2}/result/simulation/dm6_new_ont_long"

BWA_INDEX="${ANNO_PATH}/${GENOME}/BWAIndex"


ANNO_PATH="${path1}" #"${path1}"
GENOME="dm6"
CPU="10"


TE="dm6.transposon" # ALUL1SVA

TRANSPOSON_SIZE="${ANNO_PATH}/${GENOME}/${GENOME}.transposon.size"
GENOME_FA="${ANNO_PATH}/${GENOME}/benchmark/${GENOME}_template.fa"
TE_ANNO_FA="${ANNO_PATH}/${GENOME}/benchmark/${GENOME}_clean.transposon.fa"
TE_ANNO_FA2="${ANNO_PATH}/${GENOME}/benchmark/${GENOME}_clean.transposon.2.fa"

REPEATMASKER_FILE="${ANNO_PATH}/${GENOME}/benchmark/${GENOME}.repeat.bed"
REPEATMASKER_FILE="${path1}/dm6_clean/dm6_clean.rmsk.bed"
REPEAT_CLASS="${ANNO_PATH}/${GENOME}/benchmark/${GENOME}_clean.transposon.class"
CHROM_LIST="${ANNO_PATH}/run_dm6.chrom.list"



[ ! -d ${OUT_PATH_dir} ] && mkdir -p ${OUT_PATH_dir}

for SAMPLE in simulation_germ_ont  # simulation_germ_ont_0.02 simulation_soma_ont_0.02 # simulation_germ_clr simulation_germ_ccs simulation_germ_ont simulation_soma_ont simulation_soma_clr simulation_soma_ccs  #   # simulation_soma # HG001 # HG002 # SRR9685183.sra
do
    
    for dep in 50 # 1 2 3 4 5 10 20 30 40 50 # dm6
    do
        
        
        # break
        OUT_PATH=${OUT_PATH_dir}/${dep}

        # mv ${OUT_PATH}/TrEMOLO ${OUT_PATH}/TrEMOLO_rice_template
        
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
        
        if [ ! -f ${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.fa ];then
            samtools fasta ${bam_file} -@ ${CPU} > ${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.fa

            sed -i "s/:/_/g" ${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.fa
            sed -i "s/;/_/g" ${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.fa

            samtools faidx ${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.fa
        fi


        # tldr
        # l_tldr=`cat ${OUT_PATH}/tldr/${SAMPLE}.table.txt | wc -l`
        # if [ $l_tldr -eq 0 ];then
        # if [ ! -f ${OUT_PATH}/tldr/${SAMPLE}.table.txt ];then
        #     start_tldr=$(date +%s)
        #     source ${path2}/miniconda3/envs/O_O/etc/profile.d/conda.sh
        #     conda activate tldr && tldr -m 1 -p 10 --max_cluster_size 600 --flanksize 250 --max_te_len 12000 --min_te_len 150 --embed_minreads 0 -b ${bam_file} -e ${TE_ANNO_FA2} -r ${GENOME_FA} -c ${CHROM_LIST} --color_consensus -o ${OUT_PATH}/tldr/${SAMPLE} && conda deactivate
        #     end_tldr=$(date +%s)
        #     time_tldr=$(( $end_tldr - $start_tldr ))
        #     echo "time_tldr>>>>${time_tldr}"
        # fi
        


        # TELR
       
        # if [ ! -f ${OUT_PATH}/TELR/${SAMPLE}/${SAMPLE}.${dep}X.telr.vcf ];then 
        #     start_telr=$(date +%s)
        #     source ${path2}/miniconda3/envs/O_O/etc/profile.d/conda.sh
        #     conda activate TELR && telr -i ${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.fa -r ${GENOME_FA} -l ${TE_ANNO_FA} -t 10 --aligner minimap2 -o ${OUT_PATH}/TELR/${SAMPLE}  && conda deactivate
        #     end_telr=$(date +%s)
        #     time_telr=$(( $end_telr - $start_telr ))
        # fi
        


        # TrEMOLO

        if [ ! -f ${OUT_PATH}/TrEMOLO/${SAMPLE}/TE_INFOS.bed.1 ];then
            start_TrEMOLO=$(date +%s)
            cp "${SCRIPT_DIR}/config_template_dm6.yaml"  ${OUT_PATH}/TrEMOLO/${SAMPLE}/${SAMPLE}.yaml
            template_sample_fa=${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${dep}X.fa
            template_sample_fa_sed=` echo $template_sample_fa | sed 's#\/#\\\/#g'`

            template_outpath=${OUT_PATH}/TrEMOLO/${SAMPLE}
            template_outpath_sed=` echo $template_outpath | sed 's#\/#\\\/#g' `
            template_path1_sed=` echo $path1 | sed 's#\/#\\\/#g' `

            if [[ ${SAMPLE}  =~ "clr" ]];then 
                template_map_mode="map-pb" 
            fi
            if [[ ${SAMPLE}  =~ "ont" ]];then 
                template_map_mode="map-ont" 
            fi
            if [[ ${SAMPLE}  =~ "ccs" ]];then 
                template_map_mode="map-hifi" 
            fi

            echo ${template_map_mode}

            sed -i 's/template_sample_fa/'"${template_sample_fa_sed}"'/g' ${OUT_PATH}/TrEMOLO/${SAMPLE}/${SAMPLE}.yaml
            sed -i 's/template_outpath/'"${template_outpath_sed}"'/g' ${OUT_PATH}/TrEMOLO/${SAMPLE}/${SAMPLE}.yaml
            sed -i 's/template_map_mode/'"${template_map_mode}"'/g' ${OUT_PATH}/TrEMOLO/${SAMPLE}/${SAMPLE}.yaml
            sed -i 's/template_path1/'"${template_path1_sed}"'/g' ${OUT_PATH}/TrEMOLO/${SAMPLE}/${SAMPLE}.yaml

            source ${path2}/miniconda3/envs/O_O/etc/profile.d/conda.sh
            conda activate TrEMOLO && snakemake --snakefile ${path2}/tools/TrEMOLO/run.snk --configfile ${OUT_PATH}/TrEMOLO/${SAMPLE}/${SAMPLE}.yaml --cores 20 && conda deactivate
            end_TrEMOLO=$(date +%s)
            time_TrEMOLO=$(( $end_TrEMOLO - $start_TrEMOLO ))
        fi
            
        echo -e "\n-----------\n" >> ${OUT_PATH_dir}/time_simulation_tgs_rice.ont.txt
        echo -e $(date) >> ${OUT_PATH_dir}/time_simulation_tgs_rice.ont.txt
        echo -e "sample\tdepth\ttldr\ttelr\tTrEMOLO\tLOCATE" >> ${OUT_PATH_dir}/time_simulation_tgs_rice.ont.txt
        echo -e "${SAMPLE}\t${dep}\t${time_tldr}\t${time_telr}\t${time_TrEMOLO}\t${time_TEMP3}" >> ${OUT_PATH_dir}/time_simulation_tgs_rice.ont.txt

    done
    
done

# ${path2}/bin/run_intersection_simu_dm6.sh
