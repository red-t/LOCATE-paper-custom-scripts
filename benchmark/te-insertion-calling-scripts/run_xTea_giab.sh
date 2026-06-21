#!/bin/bash

path1="${path1:-path1}"
path2="${path2:-path2}"


WK_PATH="${path2}/temp/HG002_hg38"

# xtea_long -i sample_id.txt -b long_read_bam_list.txt -p ${WK_PATH}/output_hs37d5/ -o submit_jobs_hs37d5.sh --rmsk ${WK_PATH}/rep_lib_annotation/LINE/hg19/hg19_L1_larger_500_with_all_L1HS.out -r ${path1}/hs37d5/hs37d5.fa.gz --cns ${WK_PATH}/rep_lib_annotation/consensus/LINE1.fa --rep ${WK_PATH}/rep_lib_annotation/ --xtea ${path2}/bin/xTea/xtea_long/ -f 31 -y 15 -n 8 -m 32 --slurm -q 5days -t 2-00:00:00


# xtea_long -i sample_id_HG002_3.txt -b long_read_bam_list_HG002_3.txt -p ${WK_PATH}/output_hs37d5 -o submit_jobs_hs37d5_HG002_3.sh --rmsk ${WK_PATH}/rep_lib_annotation/LINE/hg19/hg19_L1_larger_500_with_all_L1HS.out -r ${path1}/hs37d5/hs37d5.fa.gz --cns ${WK_PATH}/rep_lib_annotation/consensus/LINE1.fa --rep ${WK_PATH}/rep_lib_annotation/ --xtea ${path2}/bin/xTea/xtea_long/ -f 31 -y 15 -n 8 -m 32 --slurm -q 5days -t 2-00:00:00

# xtea_long -i sample_id_HG002_1.txt -b long_read_bam_list_HG002_1.txt -p ${WK_PATH}/output_hs37d5 -o submit_jobs_hs37d5_HG002_1.sh --rmsk ${WK_PATH}/rep_lib_annotation/LINE/hg19/hg19_L1_larger_500_with_all_L1HS.out -r ${path1}/hs37d5/hs37d5.fa.gz --cns ${WK_PATH}/rep_lib_annotation/consensus/LINE1.fa --rep ${WK_PATH}/rep_lib_annotation/ --xtea ${path2}/bin/xTea/xtea_long/ -f 31 -y 15 -n 8 -m 32 --slurm -q 5days -t 2-00:00:00



# xtea_long -i sample_id_HG002_2_hg38.txt -b long_read_bam_list_HG002_2_hg38.txt -p ${WK_PATH}/output_hg38 -o submit_jobs_hg38_HG002_2.sh --rmsk ${WK_PATH}/rep_lib_annotation/LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out -r ${path1}/hg38/hg38.fa.gz --cns ${WK_PATH}/rep_lib_annotation/consensus/LINE1.fa --rep ${WK_PATH}/rep_lib_annotation/ --xtea ${path2}/bin/xTea/xtea_long/ -f 31 -y 15 -n 8 -m 32 --slurm -q 5days -t 2-00:00:00

# xtea_long -i sample_id_HG002_2_hg38.txt -b long_read_bam_list_HG002_2_hg38.txt -p ${WK_PATH}/output_hg38 -o submit_jobs_hg38_HG002_2.sh --rmsk ${WK_PATH}/rep_lib_annotation/LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out -r ${path1}/hg38/hg38.fa.gz --cns ${WK_PATH}/rep_lib_annotation/consensus/LINE1.fa --rep ${WK_PATH}/rep_lib_annotation/ --xtea ${path2}/bin/xTea/xtea_long/  -f 31 -y 15 -n 8 -m 32 --slurm -q 5days -t 2-00:00:00

# xtea_long -i sample_id_HG002_3.06.txt -b long_read_bam_list_HG002_3.06.txt -p ${WK_PATH}/output_hs37d5 -o submit_jobs_hs37d5_HG002_3.06.sh --rmsk ${WK_PATH}/rep_lib_annotation/LINE/hg19/hg19_L1_larger_500_with_all_L1HS.out -r ${path1}/hs37d5/hs37d5.fa.gz --cns ${WK_PATH}/rep_lib_annotation/consensus/LINE1.fa --rep ${WK_PATH}/rep_lib_annotation/ --xtea ${path2}/bin/xTea/xtea_long/ -f 31 -y 15 -n 8 -m 32 --slurm -q 5days -t 2-00:00:00

xtea_long -i sample_id_HG002_2_hg38.txt -b long_read_bam_list_HG002_2_hg38.txt -p ${WK_PATH}/output_hg38 -o submit_jobs_hg38_HG002_2.sh --rmsk ${WK_PATH}/rep_lib_annotation/LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out -r ${path1}/hg38/hg38.fa.gz --cns ${WK_PATH}/rep_lib_annotation/consensus/LINE1.fa --rep ${WK_PATH}/rep_lib_annotation/ --xtea ${path2}/bin/xTea/xtea_long/  -f 31 -y 15 -n 8 -m 32 --slurm -q 5days -t 2-00:00:00


