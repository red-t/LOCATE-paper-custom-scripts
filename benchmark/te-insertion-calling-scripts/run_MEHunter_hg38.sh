#!/bin/bash

path1="${path1:-path1}"
path2="${path2:-path2}"


OUT_PATH_dir="${path2}/result/simulation/new_ont_long"

for SAMPLE in simulation_germ_ont # simulation_germ_clr simulation_germ_ont simulation_soma_ccs simulation_soma_clr simulation_soma_ont # simulation_soma_ont simulation_soma_clr simulation_soma_ccs # simulation_germ_ont_0.02 simulation_soma_ont_0.02 #  # simulation_soma # HG001 # HG002 # SRR9685183.sra
do
    for depth in 1 2 3 4 5 10 20 30 40 50
    do

        #################
        # SHARED PARAMS #
        #################
        [ ! -d ${OUT_PATH_dir}/${depth}/MEHunter ] && mkdir -p ${OUT_PATH_dir}/${depth}/MEHunter
        # mv ${OUT_PATH_dir}/${depth}/MEHunter ${OUT_PATH_dir}/${depth}/MEHunter
        [ ! -d ${OUT_PATH_dir}/${depth}/MEHunter ] && mkdir ${OUT_PATH_dir}/${depth}/MEHunter
        [ ! -d ${OUT_PATH_dir}/${depth}/MEHunter/${SAMPLE} ] && mkdir ${OUT_PATH_dir}/${depth}/MEHunter/${SAMPLE}
        [ ! -d ${OUT_PATH_dir}/${depth}/MEHunter/${SAMPLE} ] && mkdir ${OUT_PATH_dir}/${depth}/MEHunter/${SAMPLE}
        [ ! -d ${OUT_PATH_dir}/${depth}/MEHunter/${SAMPLE}/cuteWork ] && mkdir ${OUT_PATH_dir}/${depth}/MEHunter/${SAMPLE}/cuteWork

        OUT_PATH=${OUT_PATH_dir}/${depth}/MEHunter/${SAMPLE}
        
        BAMPATH=${path2}/rawdata/simulation/hg38/new_ont_long/${SAMPLE}.${depth}X.bam # alignment file, .bam
        REFPATH=${path1}/from_temp3/GRCh38_no_alt.fa  # like hs37d5 or hg38, .fa
        CUTESV_WORKDIR=${OUT_PATH}/cuteWork # working directory for cuteSV intermediate files
        [ ! -d ${CUTESV_WORKDIR} ] && mkdir ${CUTESV_WORKDIR}
        CUTESV_OUTPUT=${OUT_PATH}/cuteWork/cuteSV.vcf # cuteSV output VCF
        PROCESS_NUM=10 # number of processes

        
        ##########################
        # CUTESV'S UNIQUE PARAMS #
        ##########################
        MIN_SUPPORT=1 # minimal number of supporting reads, for detailed info, see original cuteSV doc.
        cuteSVenvName=cuteSVenv # conda environment name for cuteSV

        [ -f ${OUT_PATH}/cuteWork/cuteSV_HiFi.vcf ] && rm -r ${OUT_PATH}/cuteWork/*

        if [[ ${SAMPLE}  =~ "clr" ]];then 
            max_cluster_bias_INS=100
		    diff_ratio_merging_INS=0.3
		    max_cluster_bias_DEL=200
		    diff_ratio_merging_DEL=0.5
        fi
        if [[ ${SAMPLE}  =~ "ont" ]];then 
            max_cluster_bias_INS=1000
		    diff_ratio_merging_INS=0.9
		    max_cluster_bias_DEL=1000
		    diff_ratio_merging_DEL=0.5
        fi
        if [[ ${SAMPLE}  =~ "ccs" ]];then 
            max_cluster_bias_INS=100
		    diff_ratio_merging_INS=0.3
		    max_cluster_bias_DEL=100
		    diff_ratio_merging_DEL=0.3
        fi

        start_MEHunter=$(date +%s)
        if [ ! -f ${CUTESV_OUTPUT} ];then
            source ${path2}/miniconda3/envs/O_O/etc/profile.d/conda.sh
            conda activate $cuteSVenvName
            cuteSV $BAMPATH $REFPATH $CUTESV_OUTPUT $CUTESV_WORKDIR --genotype \
                    -s $MIN_SUPPORT -t $PROCESS_NUM --report_readid --retain_work_dir --diff_ratio_filtering_TRA 1.1 \
                    --max_cluster_bias_INS ${max_cluster_bias_INS} \
                    --diff_ratio_merging_INS ${diff_ratio_merging_INS} \
                    --max_cluster_bias_DEL ${max_cluster_bias_DEL} \
                    --diff_ratio_merging_DEL ${diff_ratio_merging_DEL}
            conda deactivate
        fi


        echo "cuteSV Done"


        ############################
        # MEHunter'S UNIQUE PARAMS #
        ############################
        MEHUNTER_WORKDIR=${OUT_PATH} # working directory for MEHunter
        MEHUNTER_OUTPUT=${OUT_PATH}/MEHunter.vcf # path of MEHunter's final output, in .vcf format
        KNOWN_ME_PATH=${path2}/tools/MEHunter/ME_data/my.ME.fq # known mobile element sequences in FASTQ format
        DL_PATH=${path2}/tools/MEHunter_DL/ # See step2, path of the downloaded model.
        DEEPLEARNING_BATCHSIZE=32 # batch size of input, it highly affects memory usage, you may have to use another smaller number.
        MEHunterEnvName=MEHunterEnv # conda environment name for MEHunter



        source ${path2}/miniconda3/envs/O_O/etc/profile.d/conda.sh
        conda activate $MEHunterEnvName
        MEHunter $CUTESV_OUTPUT $BAMPATH \
                $CUTESV_WORKDIR $REFPATH $KNOWN_ME_PATH $MEHUNTER_WORKDIR $MEHUNTER_OUTPUT \
                --DL_module $DL_PATH --retain_work_dir -t $PROCESS_NUM --batch_size 32 --MAX_seqs 10
        conda deactivate

        echo "MEHunter Done"

        end_MEHunter=$(date +%s)
        time_MEHunter=$(( $end_MEHunter - $start_MEHunter ))


        echo -e "\n-----------\n" >> ${OUT_PATH_dir}/time_simulation_tgs_MEHunter.txt
        echo -e $(date) >> ${OUT_PATH_dir}/time_simulation_tgs_MEHunter.txt
        echo -e "sample\tdepth\tMEHunter" >> ${OUT_PATH_dir}/time_simulation_tgs_MEHunter.txt
        echo -e "${SAMPLE}\t${depth}\t${time_MEHunter}" >> ${OUT_PATH_dir}/time_simulation_tgs_MEHunter.txt

    done
done
