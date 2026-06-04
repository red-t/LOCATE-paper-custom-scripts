#!/bin/bash

user_path="${user_path:-${USER_PATH:-/path/to/user}}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
while getopts ":s:r:o:g:" OPTION;
do
        case $OPTION in
                s)  insert_seq=${OPTARG};;
                r)  ref_seq=${OPTARG};;
                o)  OUT_PATH=${OPTARG};;
                g)  GENOME=${OPTARG};;
        esac
done

# OUT_PATH="${user_path}/2022_long_reads/result/simulation"

echo -e ">insert_seq\n${insert_seq}" > ${OUT_PATH}/${GENOME}.test3.fa
echo -e ">ref_seq\n${ref_seq}" > ${OUT_PATH}/${GENOME}.ref3.fa


minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t 1 -aY -x map-ont ${OUT_PATH}/${GENOME}.ref3.fa ${OUT_PATH}/${GENOME}.test3.fa > ${OUT_PATH}/${GENOME}.test3.sam 2>${OUT_PATH}/${GENOME}.test3.sam.log
