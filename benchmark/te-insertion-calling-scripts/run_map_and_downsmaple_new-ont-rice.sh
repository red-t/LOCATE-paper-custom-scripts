#!/bin/bash

path1="${path1:-path1}"
path2="${path2:-path2}"



OUT_PATH="${path2}/rawdata/simulation/IRGSPv1/new_ont_v2"   # "${path2}/simulation/result/"
[ ! -d ${OUT_PATH} ] && mkdir -p ${OUT_PATH}

ANNO_PATH="${path1}" #"${path1}"
GENOME_FA="${ANNO_PATH}/IRGSPv1/IRGSPv1.fa"


# TGS mapping


## ccs 

## ont
cp ${path2}/TestingDataset/IRGSPv1/1_ont/TGS.fastq.gz ${OUT_PATH}/simulation_germ_ont.fastq.gz
seqtk seq -A ${OUT_PATH}/simulation_germ_ont.fastq.gz > ${OUT_PATH}/simulation_germ_ont.50X.fa
samtools faidx ${OUT_PATH}/simulation_germ_ont.50X.fa


SEED=1
for sample in simulation_germ_ont 
do

    if [ ! -f ${OUT_PATH}/${simulation_germ_ccs}.50X.bam ];then
       minimap2 -aYx map-ont --MD -t 30 ${GENOME_FA} ${OUT_PATH}/${sample}.50X.fa > ${OUT_PATH}/${sample}.50X.sam
       samtools sort -@ 30 -O bam -o ${OUT_PATH}/${sample}.50X.bam ${OUT_PATH}/${sample}.50X.sam
       samtools index -@ 30 ${OUT_PATH}/${sample}.50X.bam
    fi

    for fraction in 0.8 0.6 0.4 0.2 0.1 0.08 0.06 0.04 0.02
    do
        
        depth=` echo "scale=0; 50*${fraction}/1" | bc -l `

        INPUT_FILE=${OUT_PATH}/${sample}.50X.bam
        OUT_FILE=${OUT_PATH}/${sample}.${depth}X.bam

        samtools view -@ 30 -s ${SEED}.${fraction##*.} ${INPUT_FILE} -b -o ${OUT_FILE}
        samtools index ${OUT_FILE}


        samtools fasta ${OUT_FILE} > ${OUT_PATH}/${sample}.${depth}X.fa
        
        sed -i "s/:/_/g" ${OUT_PATH}/${sample}.${depth}X.fa
        sed -i "s/;/_/g" ${OUT_PATH}/${sample}.${depth}X.fa
        sed -i "s/\//_/g" ${OUT_PATH}/${sample}.${depth}X.fa

        samtools faidx ${OUT_PATH}/${sample}.${depth}X.fa


        SEED=$(($SEED + 1))
    done
    
done 


