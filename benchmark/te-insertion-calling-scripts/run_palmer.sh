#!/bin/bash

path1="${path1:-path1}"
path2="${path2:-path2}"


PALMER_PATH="${path2}/tools/PALMER"
OUT_PATH_dir="${path2}/result/simulation/hg38/" #  "${path2}/temp/simulation/result/"
GENOME_FA="${path1}/from_temp3/GRCh38_no_alt.fa"

DATA_PATH="${path2}/rawdata/simulation/hg38/tgs/"


for SAMPLE in simulation_germ_ccs # simulation_soma_ont # simulation_soma_ont simulation_soma_clr simulation_soma_ccs # simulation_germ_ont_0.02 simulation_soma_ont_0.02 # simulation_germ_clr simulation_germ_ccs simulation_germ_ont simulation_soma_ont simulation_soma_clr simulation_soma_ccs  #   # simulation_soma # HG001 # HG002 # SRR9685183.sra
do
    
    for dep in 40 50 # 10 20 30 # 2 3 4 5 # 10 20 30 40 50 # 1 2 3 4 5 10 20 30 40 50 # 4 5 10 20 30 40 50
    do
        
        # break
        OUT_PATH=${OUT_PATH_dir}/${dep}
        [ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}
        [ ! -d ${OUT_PATH}/PALMER.time ] && mkdir ${OUT_PATH}/PALMER.time
        [ ! -d ${OUT_PATH}/PALMER.time/${SAMPLE} ] && mkdir ${OUT_PATH}/PALMER.time/${SAMPLE}
        
        for TE in ALU LINE SVA;do
            [ ! -d ${OUT_PATH}/PALMER.time/${SAMPLE}/${TE}/ ] && mkdir ${OUT_PATH}/PALMER.time/${SAMPLE}/${TE}/
            # sbatch < ${OUT_PATH}/PALMER/${TE}/
            cp ${path2}/bin/submit_slurm_job.sh ${OUT_PATH}/PALMER.time/${SAMPLE}/${TE}/${SAMPLE}.${dep}X.${TE}.sh
            sed -i '$a start_palmer=$(date +%s)' ${OUT_PATH}/PALMER.time/${SAMPLE}/${TE}/${SAMPLE}.${dep}X.${TE}.sh
            echo "PREFIX=${SAMPLE}_${dep}_${TE}" >> ${OUT_PATH}/PALMER.time/${SAMPLE}/${TE}/${SAMPLE}.${dep}X.${TE}.sh
            echo  " ${PALMER_PATH}/PALMER --input ${DATA_PATH}/${SAMPLE}.${dep}X.bam  --workdir ${OUT_PATH}/PALMER.time/${SAMPLE}/${TE}/ --ref_ver GRCh38 --output ${SAMPLE}.${dep}X.${TE} --type ${TE} --mode raw --ref_fa ${GENOME_FA} " >> ${OUT_PATH}/PALMER.time/${SAMPLE}/${TE}/${SAMPLE}.${dep}X.${TE}.sh
            sed -i '$a end_palmer=$(date +%s)\ntime_palmer=$(( $end_palmer - $start_palmer ))\necho ${PREFIX}\t${time_palmer} >> ${path2}/result/simulation/hg38/time.simulation.palmer.txt' ${OUT_PATH}/PALMER.time/${SAMPLE}/${TE}/${SAMPLE}.${dep}X.${TE}.sh

            sbatch < ${OUT_PATH}/PALMER.time/${SAMPLE}/${TE}/${SAMPLE}.${dep}X.${TE}.sh
        done
    done
done
