#!/bin/bash

path1="${path1:-path1}"
path2="${path2:-path2}"

start_time=$(date +%s)

DATA_PATH="${path2}/rawdata/simulation"
OUT_PATH_dir="${path2}/result/simulation/IRGSPv1_new_ont_long_time"

ANNO_PATH="${path1}/" #"${path1}"
GENOME="IRGSPv1"
CPU="10"

TE="IRGSPv1.transposon" # ALUL1SVA

BWA_INDEX="${ANNO_PATH}/${GENOME}/BWAIndex"


GENOME_FA="${ANNO_PATH}/IRGSPv1_v2/${GENOME}.fa"
TE_ANNO_FA="${ANNO_PATH}/IRGSPv1_v2/${GENOME}.transposon.fa"
TE_ANNO_FA2="${ANNO_PATH}/IRGSPv1_v2/${GENOME}.transposon.2.fa"
TRANSPOSON_SIZE="${ANNO_PATH}/IRGSPv1_v2/${GENOME}.transposon.size"


REPEATMASKER_FILE="${ANNO_PATH}/IRGSPv1_v2/${GENOME}.rmsk.bed"
REPEAT_CLASS="${ANNO_PATH}/IRGSPv1_v2/${GENOME}.transposon.class"
CHROM_LIST="${ANNO_PATH}/IRGSPv1_v2/${GENOME}.chrom.list"



for SAMPLE in simulation_germ_ont # simulation_germ_ccs simulation_germ_clr # simulation_germ_clr simulation_germ_ont # simulation_germ_ccs simulation_germ_clr simulation_germ_ont  # simulation_germ_clr ## simulation_soma_ont simulation_soma_clr simulation_soma_ccs #   # simulation_soma # HG001 # HG002 # SRR9685183.sra
do
    for depth in 1 2 3 4 5 10 20 30 40 50 # 1 2 3 4 5 10 20 30 40 50
    do
        ### copy data file between storage locations
        # if [ ! -f ${DATA_PATH}/${GENOME}/tgs/${SAMPLE}.${depth}X.bam ];then
        #     cp ${RAW_DATA_PATH}/${GENOME}/tgs/${SAMPLE}.${depth}X.bam ${DATA_PATH}/${GENOME}/tgs/${SAMPLE}.${depth}X.bam
        # fi

        # if [ ! -f ${DATA_PATH}/${GENOME}/tgs/${SAMPLE}.${depth}X.bam.bai ];then
        #     cp ${RAW_DATA_PATH}/${GENOME}/tgs/${SAMPLE}.${depth}X.bam.bai ${DATA_PATH}/${GENOME}/tgs/${SAMPLE}.${depth}X.bam.bai
        # fi
        
        bam_file=${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${depth}X.bam
        
        
        LOCATE_version="LOCATE" # "LOCATE.reftest" # "LOCATE.ra.2"
        # LOCATE_version="v2.8_2"
        OUT_PATH=${OUT_PATH_dir}/${depth}
        

        [ ! -d ${OUT_PATH}/${LOCATE_version} ] && mkdir ${OUT_PATH}/${LOCATE_version}
        # [ -d ${OUT_PATH}/${LOCATE_version}/${SAMPLE} ] && rm -r ${OUT_PATH}/${LOCATE_version}/${SAMPLE}
        [ ! -d ${OUT_PATH}/${LOCATE_version}/${SAMPLE} ] && mkdir ${OUT_PATH}/${LOCATE_version}/${SAMPLE}

        # [ -d ${OUT_PATH}/TEMP3/${SAMPLE} ] && rm -r ${OUT_PATH}/TEMP3/${SAMPLE}

        # [ ! -d ${OUT_PATH}/TEMP3/${SAMPLE} ] && mkdir ${OUT_PATH}/TEMP3/${SAMPLE}

        if [ ! -f ${OUT_PATH}/${LOCATE_version}/${SAMPLE}/${SAMPLE}.txt ];then

            start_TEMP3=$(date +%s)

            # rm -r ${path2}/for_SMS/tmp_*

            # rm -r ${OUT_PATH}/${LOCATE_version}/${SAMPLE}/*
            # cp -r ${OUT_PATH}/LOCATE.v2.7.old/${SAMPLE}/tmp_assm ${OUT_PATH}/${LOCATE_version}/${SAMPLE}/
            # rm ${OUT_PATH}/${LOCATE_version}/${SAMPLE}/tmp_assm/*assembled*


            BAM=${bam_file} # "${path2}/rawdata/GIAB/tgs/HG002_2/bam_hs37d5/HG002_2.sort.bam" # "${path2}/temp/giab/bench/hs37d5/test.md.bam" # "${OUT_PATH}/${SAMPLE}_tgs.${dep}X.region.sorted.bam"       # BAM       Alignment 文件路径
 
            REPEAT=${REPEATMASKER_FILE} # "${ANNO_PATH}/${GENOME}/from_TEMP3/dm3.rmsk.bed"
            
            GAP="${ANNO_PATH}/${GENOME}/${GENOME}.gap.bed"
            BlackList="${ANNO_PATH}/IRGSPv1_v2/IRGSPv1.blacklist.bed"
            
            TE_TEMP3="${TE_ANNO_FA}" # "${ANNO_PATH}/${GENOME}/from_TEMP3/dm3.transposon_for_simulaTE.fa"

            GERM_MODEL="${ANNO_PATH}/IRGSPv1_v2/IRGSPv1_HighFreq"
            SOMA_MODEL="${ANNO_PATH}/IRGSPv1_v2/IRGSPv1_LowFreq"
            NPROCESS=10 # NPROCESS  最大进程数
            NTHREADS=30 # NTHREADS  每个进程最大线程数
            MINL=100    # MINL      Minimum fragment length

            echo -e "locate -b ${BAM}  -r ${REPEAT}  -g ${GAP}    -B ${BlackList} -T ${TE_TEMP3}  -R ${GENOME_FA} -H ${GERM_MODEL} -L ${SOMA_MODEL} -t ${NTHREADS}  -l ${MINL} -o ${OUT_PATH}/${LOCATE_version}/${SAMPLE} -C ${REPEAT_CLASS} "

            source ${path2}/miniconda3/envs/O_O/etc/profile.d/conda.sh
            # conda activate locate_pub2 && locate -b ${BAM} \
            #                                 -r ${REPEAT} \
            #                                 -g ${GAP} \
            #                                 -B ${BlackList} \
            #                                 -T ${TE_TEMP3} \
            #                                 -R ${GENOME_FA} \
            #                                 -H ${GERM_MODEL} \
            #                                 -L ${SOMA_MODEL} \
            #                                 -t ${NTHREADS} \
            #                                 -l ${MINL} \
            #                                 -o ${OUT_PATH}/${LOCATE_version}/${SAMPLE} \
            #                                 -C ${REPEAT_CLASS} \
            #                         && conda deactivate
            conda activate locate_pub2 && locate -b ${BAM} \
                                            -r ${REPEAT} \
                                            -g ${GAP} \
                                            -B ${BlackList} \
                                            -T ${TE_TEMP3} \
                                            -R ${GENOME_FA} \
                                            -H ${GERM_MODEL} \
                                            -L ${SOMA_MODEL} \
                                            -l ${MINL} \
                                            -t ${NTHREADS} \
                                            -o ${OUT_PATH}/${LOCATE_version}/${SAMPLE} \
                                            -C ${REPEAT_CLASS} \
                                            -n 256 \
                                    && conda deactivate
            end_TEMP3=$(date +%s)
            time_TEMP3=$(( $end_TEMP3 - $start_TEMP3 ))

            echo -e "\n-----------\n" >> ${OUT_PATH_dir}/time_simulation_tgs_rice.ccs.txt
            echo -e $(date) >> ${OUT_PATH_dir}/time_simulation_tgs_rice.ccs.txt
            echo -e "sample\tdepth\ttldr\ttelr\tTrEMOLO\tLOCATE" >> ${OUT_PATH_dir}/time_simulation_tgs_rice.ccs.txt
            echo -e "${SAMPLE}\t${dep}\t${time_tldr}\t${time_telr}\t${time_TrEMOLO}\t${time_TEMP3}" >> ${OUT_PATH_dir}/time_simulation_tgs_rice.ccs.txt
                          
            # cp ${OUT_PATH}/${LOCATE_version}/${SAMPLE}/tmp_anno/result.txt ${OUT_PATH}/${LOCATE_version}/${SAMPLE}/${SAMPLE}.txt
            # cd ${OUT_PATH}/${LOCATE_version}/${SAMPLE} && tar -zcf tmp.tar.gz tmp* --remove-files
        fi
    done


done
        
# ${path2}/bin_inter/run_intersection_simu_dm6.sh

# #################
# #### div test ###
# #################

# ${path2}/bin/run_intersection_simu_hg38_for_div.sh -v ${LOCATE_version} -o ${OUT_PATH_dir}



# for SAMPLE in simulation_germ_clr simulation_germ_ccs simulation_germ_ont  # simulation_germ_clr ## simulation_soma_ont simulation_soma_clr simulation_soma_ccs #   # simulation_soma # HG001 # HG002 # SRR9685183.sra
# do
#     genome="hg38"
#     echo -e "\n## ${genome} | ${SAMPLE}" >> ${OUT_PATH_dir}/div.new.v3.txt
#     echo -e "## TLDR-pass & LOCATE-pass overlap\n" >> ${OUT_PATH_dir}/div.new.v3.txt
#     echo -e "Depth\tTLDR\t${LOCATE_version}" >> ${OUT_PATH_dir}/div.new.v3.txt


#     for depth in 1 2 3 4 5 10 20 30 40 50 # 20 30 40 50 # 1 2 3 4 5 10 20 30 40 50
#     do
#         data_path="${path2}/result/simulation/${genome}/intersection/${depth}/intersection/dis"
#         insertion_path="${path2}/result/simulation/${genome}/${depth}/${LOCATE_version}/${SAMPLE}"
#         # TLDR-pass-LOCATE.V2.no_polyA-pass
#         # cp ${data_path}/TLDR.raw.bed ${data_path}/TLDR.1.bed
#         awk -v pre="TLDR-pass-${LOCATE_version}-pass" '{if($2==pre){print $0}}' ${data_path}/${SAMPLE}.bp_dis.txt |  awk 'BEGIN{OFS="\t"}{print $7,$8,$9,$15,$21,$6}'> ${data_path}/TLDR.1.bed
#         awk -v pre="${LOCATE_version}-pass-TLDR-pass" '{if($2==pre){print $0}}' ${data_path}/${SAMPLE}.bp_dis.txt |  awk 'BEGIN{OFS="\t"}{print $7,$8,$9,$15,$21,$6}'> ${data_path}/${LOCATE_version}-pass.bed

#         bedtools intersect -a ${data_path}/TLDR.1.bed -b ${data_path}/${LOCATE_version}-pass.bed -wa -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6,$12,$4,toupper($5),$11}' > ${data_path}/TLDR_${LOCATE_version}-pass.bed
#         # bedtools intersect -a ${data_path}/LOCATE2.bed -b ${data_path}/LOCATE3.bed -wa -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6,$12,$4,$5,$11}' > ${data_path}/LOCATE2_LOCATE3.bed
#         bedtools intersect -a ${data_path}/TLDR_${LOCATE_version}-pass.bed -b ${insertion_path}/${SAMPLE}.txt -wa -wb > ${insertion_path}/div.bed # | awk '{if(and($27,2)){print $0}}'

#         div=`cut -f 1-5 ${insertion_path}/div.bed | grep -v "\-\_\-" |awk 'BEGIN{len1=0;div1=0; len2=0;div2=0}{split($4,a,"_");split($5,b,"_"); len1=len1+a[4]; div1=div1+a[5]*a[4];  len2=len2+b[4];  div2=div2+b[5]*b[4] }END{print div1/len1,div2/len2}'`
#         echo -e ${depth}"\t"${div}"\n" >> ${OUT_PATH_dir}/div.new.v3.txt
#     done

#     #################
#     #### div test ###
#     #################
# done
