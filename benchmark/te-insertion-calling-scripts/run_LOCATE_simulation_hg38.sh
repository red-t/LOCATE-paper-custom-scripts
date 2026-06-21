#!/bin/bash

path1="${path1:-path1}"
path2="${path2:-path2}"


# RAW_DATA_PATH="${path2}/rawdata/giab/"

DATA_PATH="${path2}/rawdata/simulation/"
OUT_PATH_dir="${path2}/result/simulation/new_ont_long" # "${path2}/simulation/result/"

# RAW_DATA_PATH="${path2}/rawdata/simulation/"
# DATA_PATH="${path2}/rawdata/simulation"
# OUT_PATH_dir="${path2}/result/simulation/hyplotype_LOCATE" # "${path2}/simulation/result/"



ANNO_PATH="${path1}" #"${path1}"
GENOME="hg38"
CPU="10"





for SAMPLE in simulation_germ_ont # simulation_germ_ccs simulation_germ_clr simulation_germ_ont #  HG003 HG004 # simulation_germ_clr # simulation_germ_clr simulation_germ_ont # simulation_germ_ccs simulation_germ_clr simulation_germ_ont  # simulation_germ_clr ## simulation_soma_ont simulation_soma_clr simulation_soma_ccs #   # simulation_soma # HG001 # HG002 # SRR9685183.sra
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
        
        # bam_file=${DATA_PATH}/${SAMPLE}.${depth}X.bam # ${DATA_PATH}/${GENOME}/tgs/${SAMPLE}_2_tgs.${depth}X.bam
        
        bam_file=${DATA_PATH}/${GENOME}/new_ont_long/${SAMPLE}.${depth}X.bam

        
        LOCATE_version="LOCATE" # "LOCATE.refinePolishing1"
        # LOCATE_version="v2.8_2"
        OUT_PATH=${OUT_PATH_dir}/${depth}
        
        [ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}
        [ ! -d ${OUT_PATH}/${LOCATE_version} ] && mkdir ${OUT_PATH}/${LOCATE_version}
        # [ -d ${OUT_PATH}/${LOCATE_version}/${SAMPLE} ] && rm -r ${OUT_PATH}/${LOCATE_version}/${SAMPLE}
        [ ! -d ${OUT_PATH}/${LOCATE_version}/${SAMPLE} ] && mkdir ${OUT_PATH}/${LOCATE_version}/${SAMPLE}

        # [ -d ${OUT_PATH}/TEMP3/${SAMPLE} ] && rm -r ${OUT_PATH}/TEMP3/${SAMPLE}

        # [ ! -d ${OUT_PATH}/TEMP3/${SAMPLE} ] && mkdir ${OUT_PATH}/TEMP3/${SAMPLE}

        if [ ! -f ${OUT_PATH}/${LOCATE_version}/${SAMPLE}/result.tsv ];then

            # start_TEMP3=$(date +%s)

            # rm -r ${path2}/for_SMS/tmp_*

            # rm -r ${OUT_PATH}/${LOCATE_version}/${SAMPLE}/*
            # cp -r ${OUT_PATH}/LOCATE.v2.7.old/${SAMPLE}/tmp_assm ${OUT_PATH}/${LOCATE_version}/${SAMPLE}/
            # rm ${OUT_PATH}/${LOCATE_version}/${SAMPLE}/tmp_assm/*assembled*

            # TEMP3_PATH="${path2}/TEMP3/refinePolishing1"
            # model_path="${path2}/bin/model"

            BAM=${bam_file} # "${path2}/rawdata/GIAB/tgs/HG002_2/bam_hs37d5/HG002_2.sort.bam" # "${path2}/temp/giab/bench/hs37d5/test.md.bam" # "${OUT_PATH}/${SAMPLE}_tgs.${depth}X.region.sorted.bam"       # BAM       Alignment 文件路径
            REPEAT="${ANNO_PATH}/hg38/ALSE.bed" # ${ANNO_PATH}/from_temp3/rmsk_200.bed "${ANNO_PATH}/hg38/ALSE.bed" 
            GAP="${ANNO_PATH}/from_temp3/gap.bed"
            BlackList="${ANNO_PATH}/from_temp3/BlackList.bed"
            TE_TEMP3="${ANNO_PATH}/from_temp3/GRCh38.transposon.fa" # ${ANNO_PATH}/${GENOME}/ALSE.fa # 
            GERM_MODEL="${ANNO_PATH}/from_temp3/GRCh38_G_1V1"
            SOMA_MODEL="${ANNO_PATH}/from_temp3/GRCh38_S_1V30"
            GENOME_FA="${ANNO_PATH}/from_temp3/GRCh38_no_alt.fa" # "${ANNO_PATH}/from_temp3/GRCh38_no_alt.fa" "${path1}/hg38/hg38.fa"
            REPEAT_CLASS="${path1}/${GENOME}/hg38.repeat.class"
            NPROCESS=7 # NPROCESS  最大进程数
            NTHREADS=30 # NTHREADS  每个进程最大线程数
            MINL=100    # MINL      Minimum fragment length


            source ${path2}/miniconda3/envs/O_O/etc/profile.d/conda.sh

            conda activate locate_pub && locate -b ${BAM} \
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
                                            # -t ${NTHREADS} \
# -R GRCh38_no_alt.fa -H GRCh38_HighFreq -L GRCh38_LowFreq -o output_path

            # cp ${OUT_PATH}/${LOCATE_version}/${SAMPLE}/tmp_anno/result.txt ${OUT_PATH}/${LOCATE_version}/${SAMPLE}/${SAMPLE}.txt
            cd ${OUT_PATH}/${LOCATE_version}/${SAMPLE} && tar -zcf tmp.tar.gz tmp* --remove-files
        fi
    done


done
        


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
