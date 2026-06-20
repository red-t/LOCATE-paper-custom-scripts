#!/bin/bash

user_path="${user_path:-${USER_PATH:-/path/to/user}}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
help() {
    echo "usage: ./01_intersect_gold.sh -i input1,input2,input3 -g gold.insertion.bed -o out_path"
    echo -e "\t-i three bed format insertion file"
    echo -e "\t\t bed format: \n\t\t\t chrom start end TE meth:type:supporting_reads:frequency strand"
    echo -e "\t-g gold insertion transposon in BED format"
    echo -e "\t-l intersect length; default=5"
    echo -e "\t-s cutoff of supporting reads of germline insertion; default=0"
    echo -e "\t-f cutoff frequency of germline insertion; default=0.2"
    echo -e "\t-d sequencing depth label"
    echo -e "\t-c statistic file output path"
    echo -e "\t-r repeat file"
    echo -e "\t-b black list"
    echo -e "\t-p gap file"
    echo -e "\t-o output path"
    echo -e "\t-t data type"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    help
    exit 1
fi


while getopts ":i:l:s:g:f:d:c:r:b:p:o:t:h:" OPTION;
do
	case $OPTION in
        i)  INPUT_FILE=${OPTARG};;
        l)  intersect_len=${OPTARG};;
        s)  cutoff_supp=${OPTARG};;
        g)  gold_bed=${OPTARG};;
        f)  cutoff_freq=${OPTARG};;
        d)  depth=${OPTARG};;
        c)  static_file_path=${OPTARG};;
        r)  repeat_file=${OPTARG};;
        b)  black_list=${OPTARG};;
        p)  gap_file=${OPTARG};;
        o)  OUT_PATH=${OPTARG};;
        t)  data_type=${OPTARG};;
        h)  help && exit 1;;
	esac
done
echo ${static_file_path}


R_script_path="${user_path}/lrft2/scripts"

# rm ${OUT_PATH}/*

[ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}
[ -z ${intersect_len} ] && intersect_len=5
[ -z ${cutoff_supp} ] && cutoff_supp=0
[ -z ${cutoff_freq} ] && cutoff_freq=0
[ -z ${merge_dis} ] && merge_dis=5


# TEMP2
# sed 1d line_28.insertion.bed | awk -v l=50 -v p=TEMP2 'BEGIN{FS=OFS="\t"}{split($4,a,":");print $1,$2-l,$3+l,a[1],p":"$7":"$8":"$5,$6}'  | grep '1p1' > ${user_path}/lrft/result/simulation/intersection/TEMP2.bed
# simu
# awk -v l=50 -v p=simu 'BEGIN{FS=OFS="\t"}{split($4,a,"~");if($7==0.0001){type="soma"}else{type="germ"};print $1,$2-l,$3+l,a[2],p":"type":"5":"1,$6}' line_28.ins.summary > ${user_path}/lrft/result/simulation/intersection/simu.bed
# lrft
# sed 1d line_28_pacbio.lrft.table.txt | awk  -v l=50 -v p=lrft 'BEGIN{FS=OFS="\t"}{split($5,a,"|");split(a[3],b,"_");split($9,c,"|");split($6,d,"|");print $1,$2-l,$3+l,$4,p":"c[1]":"d[1]":"$8,b[1]}' | grep 'NI' > ${user_path}/lrft/result/simulation/intersection/lrft.bed
# tldr
# awk -v l=50 -v p=tldr 'BEGIN{FS=OFS="\t"}{print $2,$3-l,$4+l,$6,p":"$6":"$15":"$12,$5}' line_28_pacbio.table.txt | grep -v 'NA' | sed 1d > tldr.bed


# Required fields.
# chr   insertion_start insertion_end   strand  supporting_reads    frequency   


INPUT_FILES=${INPUT_FILE//,/ }
BED_NUM=0
[ -f ${OUT_PATH}/${data_type}_validation.ttl.all.bed ] && rm ${OUT_PATH}/${data_type}_validation.ttl.all.bed
for INPUT_FILE_TERM in ${INPUT_FILES[*]} # ${gold_bed}
do
    bed_file[${BED_NUM}]=${INPUT_FILE_TERM}
    bed_file_temp=${bed_file[${BED_NUM}]##*/}
    PREFIX[${BED_NUM}]=${bed_file_temp%%.bed*}
    cp ${INPUT_FILE_TERM} ${OUT_PATH}/
    cat ${INPUT_FILE_TERM} >> ${OUT_PATH}/${data_type}_validation.ttl.all.bed
    BED_NUM=$(($BED_NUM + 1))
done


LC_COLLATE=C sort -k1,1 -k2,2n ${OUT_PATH}/${data_type}_validation.ttl.all.bed  | sed "s/Alu/ALU/g" | bedtools merge -i - -c 4,5,6 -d ${merge_dis} -o collapse -delim ";"  > ${OUT_PATH}/${data_type}_validation.all.bed

# remove repeat masker
# bedtools intersect -v -a ${OUT_PATH}/${data_type}_validation.all.bed -b ${user_path}/annotation/dm3/dm3.rmsk.bed > ${OUT_PATH}/${data_type}_validation.all.ex_rmk.bed
# awk -v s=${cutoff_supp} -v f=${cutoff_freq} 'BEGIN{FS=OFS="\t"}{tag=0;split($5,a,";");for(i=1;i<=length(a);i++){split(a[i],b,":");if(b[3]>=s && b[4]>=f){tag=1}};if(tag==1){print $0}}' ${OUT_PATH}/${data_type}_validation.all.ex_rmk.bed > ${OUT_PATH}/${data_type}_validation.germline.bed

## filter by max supporting reads and max frequency
awk -v s=${cutoff_supp} -v f=${cutoff_freq} 'BEGIN{FS=OFS="\t"}{tag=0;split($5,a,";");for(i=1;i<=length(a);i++){split(a[i],b,":");if(b[3]>=s && b[4]>=f){tag=1}};if(tag==1){print $0}}' ${OUT_PATH}/${data_type}_validation.all.bed > ${OUT_PATH}/${data_type}_validation.germline.candi.bed

##### raw
cat ${OUT_PATH}/${data_type}_validation.germline.candi.bed ${gold_bed}| LC_COLLATE=C sort -k1,1 -k2,2n | bedtools merge -i - -c 4,5,6 -d ${merge_dis} -o collapse -delim ";" > ${OUT_PATH}/${data_type}_validation.germline.tmp.bed

##### test
# cat ${OUT_PATH}/${data_type}_validation.germline.candi.bed ${gold_bed}| LC_COLLATE=C sort -k1,1 -k2,2n | bedtools merge -i - -c 4,5,6 -d ${merge_dis} -o collapse -delim ";" > ${OUT_PATH}/${data_type}_validation.germline.1.bed
# bedtools intersect -a ${OUT_PATH}/${data_type}_validation.germline.1.bed -b ${user_path}/2022_long_reads/result/simulation/hg38/gold/ins.germ.internal_trunc.summary -v > ${OUT_PATH}/${data_type}_validation.germline.bed


#####   TE Class
awk 'BEGIN{a["NA"]="NA"}
     NR==FNR{a[$1]=$2}
     NR>FNR{split($4,b,";");
            for(i in b){count[b[i]]+=1};
            max_count=0;
            for (TE in count) {
                if (count[TE] > max_count) {
                    max_count = count[TE];
                    most_frequent_TE = TE;
                }
            };
            delete count;
            if(!a[most_frequent_TE]){a[most_frequent_TE]="."};print $0,a[most_frequent_TE] }' ${user_path}/annotation/repeat.class ${OUT_PATH}/${data_type}_validation.germline.tmp.bed > ${OUT_PATH}/${data_type}_validation.germline.class.bed
cp ${OUT_PATH}/${data_type}_validation.germline.class.bed ${OUT_PATH}/${data_type}_validation.germline.bed



meths=''
meths2=''
[ -f ${OUT_PATH}/${data_type}_validation.intersect.txt  ] && rm ${OUT_PATH}/${data_type}_validation.intersect.txt 
grep -no GOLD  ${OUT_PATH}/${data_type}_validation.germline.bed | awk -v me=GOLD 'BEGIN{FS=OFS="\t";b=0;i=1}{split($1,a,":");if(a[1]!=b){print me,a[1];b=a[1];i=1}else{print me,a[1]"."i;i=i+1;b=a[1]}}' >> ${OUT_PATH}/${data_type}_validation.intersect.txt 


echo ">>>${data_type}\t${depth}<<<<"

for meth in ${PREFIX[*]} #2} ${PREFIX3} ${PREFIX1}
do
    echo ${meth}
    # Use row names to identify overlaps.
    meth_name=${meth##*_}
    grep -no ${meth_name}  ${OUT_PATH}/${data_type}_validation.germline.bed | awk -v me=${meth_name} 'BEGIN{FS=OFS="\t";b=0;i=1}{split($1,a,":");if(a[1]!=b){print me,a[1];b=a[1];i=1}else{print me,a[1]"."i;i=i+1;b=a[1]}}' >> ${OUT_PATH}/${data_type}_validation.intersect.txt 
    if [ -z "${meths}" ]; then
        meths=${meth_name}
        meths2=${meth_name}
    else
        meths=${meths}" "${meth_name}
        meths2=${meths2}","${meth_name}
    fi

    if [ $((`cat ${OUT_PATH}/${meth}.bed | wc -l`)) -gt '0' ];then
        bedtools intersect -a ${OUT_PATH}/${meth}.bed -b ${gold_bed} | uniq > ${OUT_PATH}/${meth}.gold.bed
        ### caller false insertion 
        bedtools intersect -a ${OUT_PATH}/${meth}.bed -b ${gold_bed} -v > ${OUT_PATH}/${meth}.false.bed
        # bedtools intersect -a ${OUT_PATH}/${meth}.false.bed -b ${repeat_file} -c > ${OUT_PATH}/${meth}.false.n_repeat.bed
        bedtools intersect -a ${OUT_PATH}/${meth}.false.bed -b ${repeat_file} -wa -wb | awk -v meth=${meth_name} -v depth=${depth} -v data_type=${data_type} 'BEGIN{OFS="\t"}{split($10,a,"_");if($4~a[1]){ins=1}else{ins=2};split($5,b,":");print $1,$2,$3,$4,b[4],$6,data_type,meth,depth,a[1],ins}' > ${OUT_PATH}/${meth}.false.n_repeat.bed
        bedtools intersect -a ${OUT_PATH}/${meth}.false.bed -b ${repeat_file} -v | awk -v meth=${meth_name} -v depth=${depth} -v data_type=${data_type} 'BEGIN{OFS="\t"}{ins=0; split($5,b,":");print $1,$2,$3,$4,b[4],$6,data_type,meth,depth,"genome",ins}' >> ${OUT_PATH}/${meth}.false.n_repeat.bed

        bedtools intersect -a ${OUT_PATH}/${meth}.false.n_repeat.bed -b ${black_list} -c >  ${OUT_PATH}/${meth}.false.n_repeat.n_black.bed
        bedtools intersect -a ${OUT_PATH}/${meth}.false.n_repeat.n_black.bed -b ${gap_file} -c >  ${OUT_PATH}/${meth}.false.n_repeat.n_black.n_gap.bed
        cat ${OUT_PATH}/${meth}.false.n_repeat.n_black.n_gap.bed >> ${static_file_path}/false.bed
        # rm ${OUT_PATH}/${meth}.false.bed ${OUT_PATH}/${meth}.false.n_repeat.bed ${OUT_PATH}/${meth}.false.n_repeat.n_black.bed

        ### caller miss insertion  
        bedtools intersect -a ${gold_bed} -b ${OUT_PATH}/${meth}.bed -v > ${OUT_PATH}/${meth}.miss.bed
        # bedtools intersect -a ${OUT_PATH}/${meth}.miss.bed -b ${repeat_file} -c > ${OUT_PATH}/${meth}.miss.n_repeat.bed
        bedtools intersect -a ${OUT_PATH}/${meth}.miss.bed -b ${repeat_file} -wa -wb| awk -v meth=${meth_name} -v depth=${depth} -v data_type=${data_type} 'BEGIN{OFS="\t"}{split($10,a,"_");if($4~a[1]){ins=1}else{ins=2};split($5,b,":");print $1,$2,$3,$4,b[4],$6,data_type,meth,depth,a[1],ins}' > ${OUT_PATH}/${meth}.miss.n_repeat.bed
        bedtools intersect -a ${OUT_PATH}/${meth}.miss.bed -b ${repeat_file} -v | awk -v meth=${meth_name} -v depth=${depth} -v data_type=${data_type} 'BEGIN{OFS="\t"}{ins=0; split($5,b,":");print $1,$2,$3,$4,b[4],$6,data_type,meth,depth,"genome",ins}' >> ${OUT_PATH}/${meth}.miss.n_repeat.bed

        bedtools intersect -a ${OUT_PATH}/${meth}.miss.n_repeat.bed -b ${black_list} -c  > ${OUT_PATH}/${meth}.miss.n_repeat.n_black.bed
        bedtools intersect -a ${OUT_PATH}/${meth}.miss.n_repeat.n_black.bed -b ${gap_file} -c > ${OUT_PATH}/${meth}.miss.n_repeat.n_black.n_gap.bed
        cat ${OUT_PATH}/${meth}.miss.n_repeat.n_black.n_gap.bed >> ${static_file_path}/miss.bed
        # rm ${OUT_PATH}/${meth}.miss.bed ${OUT_PATH}/${meth}.miss.n_repeat.bed ${OUT_PATH}/${meth}.miss.n_repeat.n_black.bed 
    fi

    n_detected=`grep -o ${meth_name} ${OUT_PATH}/${data_type}_validation.germline.bed | wc -l`
    n_FP=`grep -o ${meth_name} ${OUT_PATH}/${data_type}_validation.germline.bed  | wc -l`
    n_intersect=`grep ${meth_name} ${OUT_PATH}/${data_type}_validation.germline.bed | grep "GOLD" | wc -l`
    n_groundTruth=`grep GOLD ${OUT_PATH}/${data_type}_validation.germline.bed | wc -l `

    TP=${n_intersect}
    FP=$(($n_detected-$n_intersect))
    FN=$(($n_groundTruth-$n_intersect))
    
    
    echo ">>>"${meth}"_"${TP}"_"${FP}"_"${FN}
    
    if [ "${n_detected}" == "0" ];then
        precision=0
    else
        precision=` echo "scale=8;$TP/$(($TP+$FP))" | bc `
    fi

    sensitivity=` echo "scale=8;$TP/$(($TP+$FN))" | bc `
    recall=` echo "scale=8;$TP/$(($TP+$FN))" | bc `
    pre_t_rec=` echo "scale=8;$precision*$recall" | bc`
    pre_a_rec=` echo "scale=8;$precision+$recall" | bc `
    
    tempt=`echo "$pre_t_rec > 0" | bc`
    if [ "${tempt}" != "1" ];then
        F1_score=0
    else
        F1_score=` echo "scale=8;2*$pre_t_rec/$pre_a_rec" | bc`
    fi
    
    # echo -e "${depth}\t${sensitivity}\t${precision}\t${recall}\t${F1_score}" >> ${static_file_path}/${meth}.statistic.txt
    # echo -e "${meth_name}\t${sensitivity}\t${precision}\t${F1_score}" >> ${static_file_path}/${data_type}.statistic.txt

    echo -e "${data_type}\t${depth}\t${meth_name}\t${sensitivity}\t${precision}\t${F1_score}" >> ${static_file_path}/performance.txt

    # echo -e "${sensitivity}\n${precision}\n${F1_score}" >> ${static_file_path}/${meth}.statistic.txt
    echo -e "${meth_name}\t${sensitivity}\t${precision}\t${F1_score}"

    # echo -e "${meth_name}\t${TP}\t${FP}" >> ${static_file_path}/${data_type}.TP_FP.count.txt
done



### TE Class

# OUT_PATH_all=${OUT_PATH}

# for TE_CLASS in LINE SINE RETROPOSON LTR DNA RC 
# do
#     echo -e "\n>>>>>>"${TE_CLASS}"\n"
#     [ ! -d ${OUT_PATH_all}/${TE_CLASS} ] && mkdir ${OUT_PATH_all}/${TE_CLASS}
#     OUT_PATH=${OUT_PATH_all}/${TE_CLASS}


#     meths=''
#     meths2=''
#     [ -f ${OUT_PATH}/${data_type}_validation.intersect.${TE_CLASS}.txt  ] && rm ${OUT_PATH}/${data_type}_validation.intersect.${TE_CLASS}.txt 
#     grep ${TE_CLASS} ${OUT_PATH_all}/${data_type}_validation.germline.bed > ${OUT_PATH}/${data_type}_validation.germline.${TE_CLASS}.bed
#     grep -no GOLD ${OUT_PATH}/${data_type}_validation.germline.${TE_CLASS}.bed | awk -v me=GOLD 'BEGIN{FS=OFS="\t";b=0;i=1}{split($1,a,":");if(a[1]!=b){print me,a[1];b=a[1];i=1}else{print me,a[1]"."i;i=i+1;b=a[1]}}' >> ${OUT_PATH}/${data_type}_validation.intersect.${TE_CLASS}.txt 


#     echo ">>>${data_type}\t${depth}<<<<"

#     for meth in ${PREFIX[*]} #2} ${PREFIX3} ${PREFIX1}
#     do
#         echo ${meth}
#         # Use row names to identify overlaps.
#         meth_name=${meth##*_}
#         grep -no ${meth_name} ${OUT_PATH}/${data_type}_validation.germline.${TE_CLASS}.bed | awk -v me=${meth_name} 'BEGIN{FS=OFS="\t";b=0;i=1}{split($1,a,":");if(a[1]!=b){print me,a[1];b=a[1];i=1}else{print me,a[1]"."i;i=i+1;b=a[1]}}' >> ${OUT_PATH}/${data_type}_validation.intersect.${TE_CLASS}.txt 
#         if [ -z "${meths}" ]; then
#             meths=${meth_name}
#             meths2=${meth_name}
#         else
#             meths=${meths}" "${meth_name}
#             meths2=${meths2}","${meth_name}
#         fi

#         n_detected=`grep -o ${meth_name} ${OUT_PATH}/${data_type}_validation.germline.${TE_CLASS}.bed | wc -l`
#         n_FP=`grep -o ${meth_name} ${OUT_PATH}/${data_type}_validation.germline.${TE_CLASS}.bed  | wc -l`
#         n_intersect=`grep ${meth_name} ${OUT_PATH}/${data_type}_validation.germline.${TE_CLASS}.bed | grep "GOLD" | wc -l`
#         n_groundTruth=`grep GOLD ${OUT_PATH}/${data_type}_validation.germline.${TE_CLASS}.bed | wc -l `

#         TP=${n_intersect}
#         FP=$(($n_detected-$n_intersect))
#         FN=$(($n_groundTruth-$n_intersect))
        
        
#         echo ">>>"${meth}"_"${TP}"_"${FP}"_"${FN}
        
#         if [ "${n_detected}" == "0" ];then
#             precision=0
#         else
#             precision=` echo "scale=8;$TP/$(($TP+$FP))" | bc `
#         fi

#         sensitivity=` echo "scale=8;$TP/$(($TP+$FN))" | bc `
#         recall=` echo "scale=8;$TP/$(($TP+$FN))" | bc `
#         pre_t_rec=` echo "scale=8;$precision*$recall" | bc`
#         pre_a_rec=` echo "scale=8;$precision+$recall" | bc `
        
#         tempt=`echo "$pre_t_rec > 0" | bc`
#         if [ "${tempt}" != "1" ];then
#             F1_score=0
#         else
#             F1_score=` echo "scale=8;2*$pre_t_rec/$pre_a_rec" | bc`
#         fi
        
#         # echo -e "${depth}\t${sensitivity}\t${precision}\t${recall}\t${F1_score}" >> ${static_file_path}/${meth}.statistic.txt
#         # echo -e "${meth_name}\t${sensitivity}\t${precision}\t${F1_score}" >> ${static_file_path}/${data_type}.statistic.txt

#         echo -e "${data_type}\t${depth}\t${meth_name}\t${sensitivity}\t${precision}\t${F1_score}" >> ${static_file_path}/performance.${TE_CLASS}.txt

#         # echo -e "${sensitivity}\n${precision}\n${F1_score}" >> ${static_file_path}/${meth}.statistic.txt
#         echo -e "${meth_name}\t${sensitivity}\t${precision}\t${F1_score}"

#         # echo -e "${meth_name}\t${TP}\t${FP}" >> ${static_file_path}/${data_type}.TP_FP.count.txt
#     done


# done



# echo ${meths2}
# awk -v me=${meths2}",GOLD" 'BEGIN{FS=OFS="\t"}{s="";split(me,a,",");for(i=1;i<=length(a);i++){b[i]=0;c=a[i];if($5~c){b[i]=1};s=s" "b[i]}print s}' ${OUT_PATH}/${data_type}_validation.germline.bed > ${OUT_PATH}/matrix.${depth}.txt

# awk 'BEGIN{OFS="\t"}{if($5~/GOLD/ && $5~/;/ && $5!~/TEMP2/ && $5!~/MELT/ ){print $0}}' ${OUT_PATH}/${data_type}_validation.germline.bed >  ${OUT_PATH}/${data_type}.tgs_spe.bed
# awk 'BEGIN{OFS="\t"}{if($5~/GOLD/ && $5~/TEMP3/ && $5!~/TEMP2/ && $5!~/MELT/ ){print $0}}' ${OUT_PATH}/${data_type}_validation.germline.bed >  ${OUT_PATH}/${data_type}.TEMP3_v_ngs.bed

# awk 'BEGIN{OFS="\t"}{if($5~/GOLD/ && $5~/;/ && $5!~/TEMP3/ && $5!~/xTea/ && $5!/TLDR/){print $0}}' ${OUT_PATH}/${data_type}_validation.germline.bed >  ${OUT_PATH}/${data_type}.ngs_spe.bed

# bedtools intersect -a ${OUT_PATH}/${data_type}.ngs_spe.bed -b ${repeat_file} -c > ${OUT_PATH}/${data_type}.ngs_spe.n_repeat.bed
# bedtools intersect -a ${OUT_PATH}/${data_type}.ngs_spe.n_repeat.bed -b ${black_list} -c  > ${OUT_PATH}/${data_type}.ngs_spe.n_repeat.n_black.bed
# bedtools intersect -a ${OUT_PATH}/${data_type}.ngs_spe.n_repeat.n_black.bed -b ${gap_file} -c > ${OUT_PATH}/${data_type}.ngs_spe.n_repeat.n_black.n_gap.bed
# rm ${OUT_PATH}/${data_type}.ngs_spe.bed ${OUT_PATH}/${data_type}.ngs_spe.n_repeat.bed ${OUT_PATH}/${data_type}.ngs_spe.n_repeat.n_black.bed 


# bedtools intersect -a ${OUT_PATH}/${data_type}.tgs_spe.bed -b ${repeat_file} -c > ${OUT_PATH}/${data_type}.tgs_spe.n_repeat.bed
# bedtools intersect -a ${OUT_PATH}/${data_type}.tgs_spe.n_repeat.bed -b ${black_list} -c  > ${OUT_PATH}/${data_type}.tgs_spe.n_repeat.n_black.bed
# bedtools intersect -a ${OUT_PATH}/${data_type}.tgs_spe.n_repeat.n_black.bed -b ${gap_file} -c > ${OUT_PATH}/${data_type}.tgs_spe.n_repeat.n_black.n_gap.bed
# rm ${OUT_PATH}/${data_type}.tgs_spe.bed ${OUT_PATH}/${data_type}.tgs_spe.n_repeat.bed ${OUT_PATH}/${data_type}.tgs_spe.n_repeat.n_black.bed 



# echo ""
