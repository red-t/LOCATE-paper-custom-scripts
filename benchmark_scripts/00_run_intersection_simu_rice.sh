#!/bin/bash
start=$(date +%s)
user_path="${user_path:-${USER_PATH:-/path/to/user}}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# while getopts ":v:" OPTION;
# do
#     case $OPTION in
#         v)  LOCATE_version=${OPTARG};;
#     esac
# done

LOCATE_version="LOCATE" # .refinePolishing1.3"

INSERTION_PATH_dir="${user_path}/2022_long_reads/result/simulation/IRGSPv1"
DATA_PATH_dir="${user_path}/2022_long_reads/result/simulation/IRGSPv1/insertion_revision"
OUT_PATH_dir="${user_path}/2022_long_reads/result/simulation/IRGSPv1/intersection_revision/" # "${user_path}/lrft2/result" # "${user_path}/lrft/result/simulation" 
lrft_BIN_PATH="${user_path}/2022_long_reads/scripts/"
ANNO_PATH="${user_path}/annotation"
GENOME="IRGSPv1"
SAMPLE="simulation"
CPU="15"

R_script_path="${user_path}/2022_long_reads/scripts/"

# intersection
# [ -d ${OUT_PATH_dir} ] && rm -r ${OUT_PATH_dir} 
[ ! -d ${INSERTION_PATH_dir}/figure ] && mkdir ${INSERTION_PATH_dir}/figure
[ ! -d ${DATA_PATH_dir} ] && mkdir -p ${DATA_PATH_dir}
[ ! -d ${OUT_PATH_dir} ] && mkdir -p ${OUT_PATH_dir}
[ ! -d ${OUT_PATH_dir}/pp_intersection ] && mkdir ${OUT_PATH_dir}/pp_intersection
[ ! -d ${OUT_PATH_dir}/pp_intersection/data ] && mkdir ${OUT_PATH_dir}/pp_intersection/data
[ ! -d ${OUT_PATH_dir}/pp_intersection/figure ] && mkdir ${OUT_PATH_dir}/pp_intersection/figure
# rm ${OUT_PATH_dir}/pp_intersection/data/*
[ -f ${OUT_PATH_dir}/pp_intersection/data/performance.txt ] && rm ${OUT_PATH_dir}/pp_intersection/data/performance*.txt
[ -f ${OUT_PATH_dir}/pp_intersection/data/miss.bed ] && rm ${OUT_PATH_dir}/pp_intersection/data/miss.bed
[ -f ${OUT_PATH_dir}/pp_intersection/data/false.bed ] && rm ${OUT_PATH_dir}/pp_intersection/data/false.bed

TRANSPOSON_SIZE="${user_path}/annotation/IRGSPv1/IRGSPv1.transposon.fa.fai"


simulation_sample="simulation_germ_clr"

ex_ch="_alt|random|"

# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;b[$1]=$2}NR>FNR{print $1,$2,$3,$4,$5,$6,$7,$8,$9,b[$4],a[$4]}' ins.soma.insSeq ins.soma.summary > ins.soma.gold.summary
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$3;b[$1]=$2}NR>FNR{print $1,$2,$3,$4,$5,$6,$7,$8,$9,b[$4],a[$4]}' ins.germ.insSeq ins.germ.summary > ins.germ.gold.summary

[ -f ${OUT_PATH}/intersection/dis/bp_dis.txt ] && rm ${OUT_PATH}/intersection/dis/bp_dis.txt

for simulation_sample in simulation_germ_ccs simulation_germ_clr simulation_germ_ont #  simulation_soma_ccs simulation_soma_ont simulation_soma_clr 
do
    for depth in 1 2 3 4 5 10 20 30 40 50
    do 

        OUT_PATH=${OUT_PATH_dir}/${depth}
        [ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}
        [ ! -d ${OUT_PATH}/intersection ] && mkdir ${OUT_PATH}/intersection
        [ ! -d ${OUT_PATH}/intersection/pdf ] && mkdir ${OUT_PATH}/intersection/pdf
        [ ! -d ${OUT_PATH}/intersection/dis ] && mkdir ${OUT_PATH}/intersection/dis



        DATA_PATH=${DATA_PATH_dir}/${depth}
        [ ! -d ${DATA_PATH} ] && mkdir ${DATA_PATH}
        [ ! -d ${DATA_PATH}/MELT ] && mkdir ${DATA_PATH}/MELT
        [ ! -d ${DATA_PATH}/TEMP2 ] && mkdir ${DATA_PATH}/TEMP2
        [ ! -d ${DATA_PATH}/LOCATE ] && mkdir ${DATA_PATH}/LOCATE
        [ ! -d ${DATA_PATH}/TEMP3_test_model ] && mkdir ${DATA_PATH}/TEMP3_test_model
        [ ! -d ${DATA_PATH}/TELR ] && mkdir ${DATA_PATH}/TELR
        [ ! -d ${DATA_PATH}/TrEMOLO ] && mkdir ${DATA_PATH}/TrEMOLO
        [ ! -d ${DATA_PATH}/tldr ] && mkdir ${DATA_PATH}/tldr

        [ ! -d ${DATA_PATH}/GraffiTE ] && mkdir ${DATA_PATH}/GraffiTE


        INSERTION_PATH=${INSERTION_PATH_dir}/${depth}


        if [[ $simulation_sample =~ "soma" ]];then
            ns_num=0
            insertion_type="soma"
        else
            
            if [ ${depth} -lt 10 ];then
                ns_num=1
            elif [ ${depth} -lt 20 ];then
                ns_num=2
            else
                ns_num=3
            fi
            ns_num=1
            insertion_type="germ"
        fi

        ### NGS copy
        if [ ! -d ${INSERTION_PATH}/melt/${simulation_sample} ];then
            ln -s ${INSERTION_PATH}/melt/simulation_${insertion_type} ${INSERTION_PATH}/melt/${simulation_sample}
            
        fi

        if [ ! -d ${INSERTION_PATH}/TEMP2/${simulation_sample} ];then
            cp -r ${INSERTION_PATH}/TEMP2/simulation_${insertion_type} ${INSERTION_PATH}/TEMP2/${simulation_sample}
            cp ${INSERTION_PATH}/TEMP2/${simulation_sample}/simulation_${insertion_type}_${depth}X.insertion.bed ${INSERTION_PATH}/TEMP2/${simulation_sample}/${simulation_sample}.insertion.bed
        fi 


        ## gold
        # simulation_bnech
        GOLD_INS_SUMM=${DATA_PATH}/../gold/ins.${insertion_type}.gold.summary 
        awk -v l=10 -v p=GOLD 'BEGIN{FS=OFS="\t"}{split($4,a,"~");print $1, $8-l, $9+l, a[2], p":gold:5:"$7":"$11, $6}' ${GOLD_INS_SUMM} > ${OUT_PATH}/intersection/${simulation_sample}_GOLD.bed
        
        gold_ins=${OUT_PATH}/intersection/${simulation_sample}_GOLD.bed

        ## MELT
        RESULT_PATH="${INSERTION_PATH}/melt/${simulation_sample}"
        [ -f ${RESULT_PATH}/melt.vcf.bed ] && rm ${RESULT_PATH}/melt.vcf.bed
        while read info; do
            te_info=(${info})
            melt_te=${te_info[0]%%_*}
            if [ -f ${RESULT_PATH}/${melt_te}.final_comp.vcf ];then
                # echo ${melt_te}
                # else
                if [ $((`cat ${RESULT_PATH}/${melt_te}.final_comp.vcf | wc -l`)) -gt '0' ];then
                    # bcftools query -f "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/LP\t%INFO/RP\t%INFO/SVLEN\t%INFO/TSD\t%INFO/MEINFO\t%INFO/SVLEN\n" ${RESULT_PATH}/${melt_te}.final_comp.vcf | awk -v l=10 'BEGIN{FS=OFS="\t"}{split($8,a,","); print $1,$2-l,$2+l,$3,"melt",a[4],$4+$5,$9,a[2],a[3],"-",$10}' | grep -v "chr_FBgn0003055" >> ${RESULT_PATH}/melt.vcf.bed
                    bcftools query -f "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/LP\t%INFO/RP\t%INFO/SVLEN\t%INFO/TSD\t%INFO/MEINFO\t%INFO/SVLEN\t[%GT]\n" ${RESULT_PATH}/${melt_te}.final_comp.vcf | awk -v l=10 'BEGIN{FS=OFS="\t"}{split($8,a,","); if($10=="1/1"){genotype="1"}else{genotype="0.5"}; print $1,$2-l,$2+l,$3,"melt",a[4],$4+$5,$9,a[2],a[3],"-",genotype}' | grep -v "chr_FBgn0003055" >> ${RESULT_PATH}/melt.vcf.bed

                fi
            fi
        done <<< "$( cat ${TRANSPOSON_SIZE} | awk '{split($1,a,"_");print a[1]}' )"
        
        if [ ! -f ${RESULT_PATH}/melt.vcf.bed ];then
            cat /dev/null > ${RESULT_PATH}/melt.vcf.bed
            cp ${RESULT_PATH}/melt.vcf.bed ${DATA_PATH}/MELT/${simulation_sample}.tmp.bed
        else
            cat ${RESULT_PATH}/melt.vcf.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools merge -i - -o first -c 4,5,6,7,8,9,10,11,12 > ${DATA_PATH}/MELT/${simulation_sample}.tmp.bed
        fi

        awk 'BEGIN{OFS="\t"}NR==FNR{split($1,a,"_");b[a[1]]=$1}NR>FNR{print $1,$2,$3,b[$4],$5,$6,$7,$12}' ${user_path}/annotation/IRGSPv1/IRGSPv1.transposon.2.class ${DATA_PATH}/MELT/${simulation_sample}.tmp.bed > ${DATA_PATH}/MELT/${simulation_sample}.bed
        
        cat ${DATA_PATH}/MELT/${simulation_sample}.bed | awk -v l=10 -v p=MELT 'BEGIN{FS=OFS="\t"}{print $1,$2-l,$3+l,$4,p":germ:"$7":"$8":-",$6}' > ${OUT_PATH}/intersection/${simulation_sample}_MELT.bed
        melt_ins=${OUT_PATH}/intersection/${simulation_sample}_MELT.bed


        ## TEMP2
        cp ${INSERTION_PATH}/TEMP2/${simulation_sample}/${simulation_sample}.insertion.bed ${DATA_PATH}/TEMP2/${simulation_sample}.insertion.bed
        if [[ $simulation_sample =~ "soma" ]];then
            sed 1d ${DATA_PATH}/TEMP2/${simulation_sample}.insertion.bed  | awk -v l=10 -v p=TEMP2 'BEGIN{FS=OFS="\t"}{split($4,a,":");if($8==1){print $1,$2-l,$3+l,a[1],p":"$7":"$8":"$5":-",$6}}' > ${OUT_PATH}/intersection/${simulation_sample}_TEMP2.bed
        else
            sed 1d ${DATA_PATH}/TEMP2/${simulation_sample}.insertion.bed | awk -v l=10 -v p=TEMP2 -v ns=${ns_num} 'BEGIN{FS=OFS="\t"}{split($4,a,":");if($8>ns){print $1,$2-l,$3+l,a[1],p":"$7":"$8":"$5":-",$6}}' > ${OUT_PATH}/intersection/${simulation_sample}_TEMP2.bed
        fi
        # sed 1d ${DATA_PATH}/TEMP2/${simulation_sample}_hg38X.insertion.bed  | awk -v l=10 -v p=TEMP2 -v ns=${ns_num} 'BEGIN{FS=OFS="\t"}{split($4,a,":");if($8>ns){print $1,$2-l,$3+l,a[1],p":"$7":"$8":"$5,$6}}' > ${OUT_PATH}/intersection/${simulation_sample}_TEMP2.bed
        temp2_ins=${OUT_PATH}/intersection/${simulation_sample}_TEMP2.bed
        echo 'cl1'


        #### tldr
        awk '{if($3>0){print $0}}' ${INSERTION_PATH}/tldr/${simulation_sample}.table.txt > ${DATA_PATH}/tldr/${simulation_sample}.table.txt
        if [[ $simulation_sample =~ "soma" ]];then
            sed 1d ${DATA_PATH}/tldr/${simulation_sample}.table.txt  | grep PASS | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"  | awk -v l=10 -v p=TLDR -v ns=${ns_num} 'BEGIN{FS=OFS="\t"}{if($15==1){if($1>0 && $6!= "A"){if($3>l){st=$3-l}else{st=0}; print $2,st,$4+l,$6,p":"$6":"$15":"$12":"$22,$5}}}' | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_TLDR.bed
            sed 1d ${DATA_PATH}/tldr/${simulation_sample}.table.txt  | grep PASS | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"  | awk -v l=10 -v p=TLDR-pass -v ns=${ns_num} 'BEGIN{FS=OFS="\t"}{if($15==1){if($1>0 && $6!= "A"){if($3>l){st=$3-l}else{st=0}; print $2,st,$4+l,$6,p":"$6":"$15":"$12":"$22,$5}}}' | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_TLDR-pass.bed
            sed 1d ${DATA_PATH}/tldr/${simulation_sample}.table.txt  | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"  | awk -v l=10 -v p=TLDR-all -v ns=${ns_num} 'BEGIN{FS=OFS="\t"}{if($15==1){if($1>0 && $6!= "A"){if($3>l){st=$3-l}else{st=0}; print $2,st,$4+l,$6,p":"$6":"$15":"$12":"$22,$5}}}' | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_TLDR-all.bed
        else
            sed 1d ${DATA_PATH}/tldr/${simulation_sample}.table.txt  | grep PASS | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"  | awk -v l=10 -v p=TLDR -v ns=${ns_num} 'BEGIN{FS=OFS="\t"}{if($15>ns){if($1>0 && $6!= "A"){if($3>l){st=$3-l}else{st=0}; print $2,st,$4+l,$6,p":"$6":"$15":"$12":"$22,$5}}}'  | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_TLDR.bed
            sed 1d ${DATA_PATH}/tldr/${simulation_sample}.table.txt  | grep PASS | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"  | awk -v l=10 -v p=TLDR-pass -v ns=${ns_num} 'BEGIN{FS=OFS="\t"}{if($15>ns){if($1>0 && $6!= "A"){if($3>l){st=$3-l}else{st=0}; print $2,st,$4+l,$6,p":"$6":"$15":"$12":"$22,$5}}}' | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_TLDR-pass.bed
            sed 1d ${DATA_PATH}/tldr/${simulation_sample}.table.txt  | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"  | awk -v l=10 -v p=TLDR-all -v ns=${ns_num} 'BEGIN{FS=OFS="\t"}{if($15>ns){if($1>0 && $6!= "A"){if($3>l){st=$3-l}else{st=0}; print $2,st,$4+l,$6,p":"$6":"$15":"$12":"$22,$5}}}' | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_TLDR-all.bed

        fi
        # sed 1d ${DATA_PATH}/tldr/simulation_hg38.table.txt | grep PASS | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]/|/g" | awk -v l=5 -v p=TLDR 'BEGIN{FS=OFS="\t"}{seq="";split($22,s,"|");for(i=2;i<=length(s)-1;i++){seq=seq""s[i]};print $2,$3-l,$4+l,$6,p":"$6":"$15":"$12,$5}'  > ${OUT_PATH}/intersection/TLDR.bed
        tldr_ins=${OUT_PATH}/intersection/${simulation_sample}_TLDR.bed
        tldr_pass_ins=${OUT_PATH}/intersection/${simulation_sample}_TLDR-pass.bed
        tldr_all_ins=${OUT_PATH}/intersection/${simulation_sample}_TLDR-all.bed
        echo 'cl2'


        # TELR
        vcf_result=` cat ${INSERTION_PATH}/TELR/${simulation_sample}/${simulation_sample}.${depth}X.telr.vcf | wc -l `
        if [ $vcf_result -gt 0 ];then
            sed 's/ID=[0-9]\+;//g' ${INSERTION_PATH}/TELR/${simulation_sample}/${simulation_sample}.${depth}X.telr.vcf | bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/FAMILY\t%INFO/AF\t%INFO/STRANDS\t%INFO/RE\t%INFO/TSD_SEQ\t%ALT\n" - | awk 'BEGIN{OFS="\t"}{split($4,a,"|");if($2>$3){st=$2;en=$2}else{st=$2;en=$3};print $1,st,en,a[1],$5,$6,$7,$8,$9}'> ${DATA_PATH}/TELR/${simulation_sample}.telr.bed
            awk -v l=10 -v p=TELR -v ns=${ns_num} 'BEGIN{OFS="\t"}{split($4,a,"|");if($2>$3){st=$3;en=$2}else{st=$2;en=$3};print $1,st,en+1,a[1],p":INS:"$7":"$5":"$9,$6}' ${DATA_PATH}/TELR/${simulation_sample}.telr.bed > ${OUT_PATH}/intersection/${simulation_sample}_TELR.bed

        else
            cat /dev/null > ${INSERTION_PATH}/TELR/${simulation_sample}/${simulation_sample}.${depth}X.telr.vcf
            cat /dev/null > ${DATA_PATH}/TELR/${simulation_sample}.telr.bed
            cat /dev/null > ${OUT_PATH}/intersection/${simulation_sample}_TELR.bed
        fi
        telr_ins=${OUT_PATH}/intersection/${simulation_sample}_TELR.bed
        echo 'cl3'


        # TrEMOLO
        
        bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/RE\t%ALT\t%INFO/AF\n"  ${INSERTION_PATH}/TrEMOLO/${simulation_sample}/OUTSIDER/VARIANT_CALLING/SV.vcf | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"sniffles."$5"."$4,$6,$7,$8}' > ${INSERTION_PATH}/TrEMOLO/${simulation_sample}/OUTSIDER/VARIANT_CALLING/SV.supp.bed
        sed 1d ${INSERTION_PATH}/TrEMOLO/${simulation_sample}/OUTSIDER/TrEMOLO_SV_TE/INS/INS_TREMOLO.csv | awk 'BEGIN{OFS="\t"}{split($2,a,":");print a[1],a[3],a[4],a[5],a[6],a[9],1}' >> ${INSERTION_PATH}/TrEMOLO/${simulation_sample}/OUTSIDER/VARIANT_CALLING/SV.supp.bed
        cat ${INSERTION_PATH}/TrEMOLO/${simulation_sample}/OUTSIDER/ET_FIND_FA/*fasta  | awk '{if($1~/>/){split($1,a,":");s=a[5]}else{print s, $1;s=""}}' > ${INSERTION_PATH}/TrEMOLO/${simulation_sample}/OUTSIDER/ET_FIND_FA/TE_INS.bed
        awk 'BEGIN{OFS="\t"}NR==FNR{a[$4]=$5"_"$7}NR>FNR{split($4,b,"|");split(a[b[2]], c, "_"); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,c[1],c[2] }' ${INSERTION_PATH}/TrEMOLO/${simulation_sample}/OUTSIDER/VARIANT_CALLING/SV.supp.bed ${INSERTION_PATH}/TrEMOLO/${simulation_sample}/TE_INFOS.bed | grep -v "DEL" > ${INSERTION_PATH}/TrEMOLO/${simulation_sample}/TE_INFOS.supp.bed

        grep chr ${INSERTION_PATH}/TrEMOLO/${simulation_sample}/TE_INFOS.supp.bed | sed 1d | awk 'BEGIN{OFS="\t"}{if($7>0){split($4,a,"|");if($2>$10){st=$10;en=$3}else{st=$2;en=$10};  print $1,st,en,a[1],$11,$5,$8,$9,1,$10,a[2],$12}}' | awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2}NR>FNR{if($11 in a){seq=a[$11]}else{seq="-"}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,seq,$12}' ${INSERTION_PATH}/TrEMOLO/${simulation_sample}/OUTSIDER/ET_FIND_FA/TE_INS.bed - > ${DATA_PATH}/TrEMOLO/${simulation_sample}.TrEMOLO.bed
        awk -v l=10 -v p=TrEMOLO -v ns=${ns_num} 'BEGIN{OFS="\t"}{split($4,a,"|"); if($2>$10){st=$10;en=$3}else{st=$2;en=$10};if($5>ns){print $1,st,en,a[1],p":INS:"$5":"$12":"$11,$6}}' ${DATA_PATH}/TrEMOLO/${simulation_sample}.TrEMOLO.bed > ${OUT_PATH}/intersection/${simulation_sample}_TrEMOLO.bed
        TrEMOLO_ins=${OUT_PATH}/intersection/${simulation_sample}_TrEMOLO.bed



        ### GraffiTE
        [ ! -d ${INSERTION_PATH}/GraffiTE ] && mkdir ${INSERTION_PATH}/GraffiTE
        # # cp -r ${user_path}/2022_long_reads/result/simulation/IRGSPv1/${depth}/GraffiTE/${simulation_sample}/ ${INSERTION_PATH}/GraffiTE/${simulation_sample}
        # bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/repeat_ids\t%INFO/SUPPORT\t%INFO/RM_hit_strands\t%ID\t%ALT\t[%GT]\n" ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/genotypes_repmasked_filtered.vcf > ${INSERTION_PATH}/GraffiTE/${simulation_sample}/GraffiTE.tmp.bed
        # # awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$9"_"$11"_"length($12)"_"$12}NR>FNR{if($7 in a){split(a[$7],b,"_"); if(b[1]=="0"){st=$2-b[3];en=$2}else{st=$2;en=$2+b[3]} }else{st=$2;en=$2+1}; split($4,c,",");split($6,d,",");if(d[1]=="C"){strand="-"}else{strand="+"}; print $1,st,en,c[1],$5,strand,$7,$8,a[$7] }' ${INSERTION_PATH}/GraffiTE/${simulation_sample}/3_TSD_search/TSD_summary.txt  ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/GraffiTE.tmp.bed > ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/${simulation_sample}_GraffiTE.bed
        # awk 'BEGIN{OFS="\t"}NR==FNR{if($14=="PASS"){a[$1]=$9"_"$11"_"length($12)"_"$12}}
        #                     NR>FNR{if($7 in a){split(a[$7],b,"_"); st=$2;en=$2+b[3];TSD=b[4];split($8,c,TSD);if(c[1]!=$8){ins_st=length(c[1])+b[3]+1;ins_len=length($8)-ins_st+1; ins_seq=substr($8,ins_st,ins_len)}else{ins_seq=substr($8,1,length(c[1]))} }
        #                            else{a[$7]=".";st=$2;en=$2+1;ins_seq=$8}; split($4,c,",");split($6,d,","); 
        #                                 if(d[1]=="C"){strand="-"}else{strand="+"}; 
        #                            if($9=="1|0"){AF=0.5}else{AF=1};
        #                            print $1,st,en,c[1],$5,strand,$7,ins_seq,a[$7],AF }' ${INSERTION_PATH}/GraffiTE/${simulation_sample}/3_TSD_search/TSD_summary.txt  ${INSERTION_PATH}/GraffiTE/${simulation_sample}/GraffiTE.tmp.bed | grep -v "DEL" > ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/${simulation_sample}_GraffiTE.bed
        
        # cp ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/${simulation_sample}_GraffiTE.bed ${DATA_PATH}/GraffiTE/${simulation_sample}.GraffiTE.bed
        # awk -v l=10 -v p=GraffiTE -v ns=${ns_num} 'BEGIN{OFS="\t"}{if($5>ns){print $1,$2,$3,$4,p":INS:"$5":"$10":"$8,$6}}' ${DATA_PATH}/GraffiTE/${simulation_sample}.GraffiTE.bed > ${OUT_PATH}/intersection/${simulation_sample}_GraffiTE.bed
        # GraffiTE_ins=${OUT_PATH}/intersection/${simulation_sample}_GraffiTE.bed


        # INSERTION_PATH_GT="${user_path}/2022_long_reads/result/simulation/IRGSPv1/test_graffiTE/${depth}"
        # genotype = true
        bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/repeat_ids\t%INFO/SUPPORT\t%INFO/RM_hit_strands\t%ID\t%ALT\t[%GT]\n" ${INSERTION_PATH}/GraffiTE/${simulation_sample}/4_Genotyping/GraffiTE.merged.genotypes.vcf.gz > ${INSERTION_PATH}/GraffiTE/${simulation_sample}/GraffiTE.tmp.bed

        # awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$9"_"$11"_"length($12)"_"$12}NR>FNR{if($7 in a){split(a[$7],b,"_"); if(b[1]=="0"){st=$2-b[3];en=$2}else{st=$2;en=$2+b[3]} }else{st=$2;en=$2+1}; split($4,c,",");split($6,d,",");if(d[1]=="C"){strand="-"}else{strand="+"}; print $1,st,en,c[1],$5,strand,$7,$8,a[$7] }' ${INSERTION_PATH}/GraffiTE/${simulation_sample}/3_TSD_search/TSD_summary.txt  ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/GraffiTE.tmp.bed > ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/${simulation_sample}_GraffiTE.bed
        awk 'BEGIN{OFS="\t"}NR==FNR{if($14=="PASS"){a[$1]=$9"_"$11"_"length($12)"_"$12}}
                            NR>FNR{if($7 in a){split(a[$7],b,"_"); st=$2;en=$2+b[3];TSD=b[4];split($8,c,TSD);if(c[1]!=$8){ins_st=length(c[1])+b[3]+1;ins_len=length($8)-ins_st+1; ins_seq=substr($8,ins_st,ins_len)}else{ins_seq=substr($8,1,length(c[1]))} }
                                   else{a[$7]=".";st=$2;en=$2+1;ins_seq=$8}; split($4,c,",");split($6,d,","); 
                                        if(d[1]=="C"){strand="-"}else{strand="+"}; 
                                   if($9=="0/1"){AF=0.5}else{ if($9=="1/1"){AF=1}else{AF=0} };
                                   print $1,st,en,c[1],$5,strand,$7,ins_seq,a[$7],AF }' ${INSERTION_PATH}/GraffiTE/${simulation_sample}/3_TSD_search/TSD_summary.txt  ${INSERTION_PATH}/GraffiTE/${simulation_sample}/GraffiTE.tmp.bed | grep -v "DEL" > ${INSERTION_PATH}/GraffiTE/${simulation_sample}/4_Genotyping/${simulation_sample}_GraffiTE.bed
        
        cp ${INSERTION_PATH}/GraffiTE/${simulation_sample}/4_Genotyping/${simulation_sample}_GraffiTE.bed ${DATA_PATH}/GraffiTE/${simulation_sample}.GraffiTE.bed
        awk -v l=10 -v p=GraffiTE -v ns=${ns_num} 'BEGIN{OFS="\t"}{if($5>ns){print $1,$2,$3,$4,p":INS:"$5":"$10":"$8,$6}}' ${DATA_PATH}/GraffiTE/${simulation_sample}.GraffiTE.bed > ${OUT_PATH}/intersection/${simulation_sample}_GraffiTE.bed
        GraffiTE_ins=${OUT_PATH}/intersection/${simulation_sample}_GraffiTE.bed


        # LOCATE
        # cat ${INSERTION_PATH}/${LOCATE_version}/${simulation_sample}/${simulation_sample}.txt > ${DATA_PATH}/LOCATE/${simulation_sample}.txt
        
        sed 1d ${INSERTION_PATH}/${LOCATE_version}/${simulation_sample}/result.tsv | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14,$15,$16,$17}' >  ${DATA_PATH}/LOCATE/${simulation_sample}.txt

        if [[ $simulation_sample =~ "soma" ]];then
            awk '{if(and($19,1)){print $0}}' ${DATA_PATH}/LOCATE/${simulation_sample}.txt | awk -v ns=${ns_num}  'BEGIN{FS=OFS="\t"}{if($10==1){gsub(/N*/,"",$13); print $1,$2,$3,$4,"LOCATE-pass:germ:"$10":"$5":"$16, $6}}'  > ${OUT_PATH}/intersection/${simulation_sample}_LOCATE-pass.bed
            cat ${DATA_PATH}/LOCATE/${simulation_sample}.txt | awk -v ns=${ns_num}  'BEGIN{OFS="\t"}{if($10==1){gsub(/N*/,"",$13); print $1,$2,$3,$4,"LOCATE-all:germ:"$10":"$5":"$16, $6}}'  > ${OUT_PATH}/intersection/${simulation_sample}_LOCATE-all.bed

        else
            awk '{if( $8=="True" ){print $0}}' ${DATA_PATH}/LOCATE/${simulation_sample}.txt | awk -v ns=${ns_num}  'BEGIN{FS=OFS="\t"}{if($11>ns){gsub(/N*/,"",$13); if($7=="1/1"){geno=1}else{if($7=="0/1"){geno=0.5}else{geno=0}}; print $1,$2+1,$3+1,$4,"LOCATE-pass:germ:"$11":"geno":"$13, $6}}'  | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_LOCATE-pass.bed
            cat ${DATA_PATH}/LOCATE/${simulation_sample}.txt | awk -v ns=${ns_num}  'BEGIN{OFS="\t"}{if($11>ns){ gsub(/N*/,"",$13);if($7=="1/1"){geno=1}else{if($7=="0/1"){geno=0.5}else{geno=0}};   print $1,$2+1,$3+1,$4,"LOCATE-all:germ:"$11":"geno":"$13, $6}}'  | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_LOCATE-all.bed

        fi
        # awk '{if(and($19,1)){print $0}}' ${DATA_PATH}/LOCATE/${simulation_sample}.txt | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,"TEMP3:germ:"$10":"$5":"$16, $6}'  > ${OUT_PATH}/intersection/${simulation_sample}_LOCATE.bed
        temp3_all_ins=${OUT_PATH}/intersection/${simulation_sample}_LOCATE-all.bed
        temp3_pass_ins=${OUT_PATH}/intersection/${simulation_sample}_LOCATE-pass.bed

        echo ">>>>intersection"
        # intersectte
        "${script_dir}/01_intersect_gold.sh" -i ${temp2_ins},${melt_ins},${tldr_all_ins},${tldr_pass_ins},${telr_ins},${temp3_all_ins},${temp3_pass_ins},${TrEMOLO_ins},${GraffiTE_ins} -g ${gold_ins} -d ${depth} -t ${simulation_sample} -c ${OUT_PATH_dir}/pp_intersection/data -o ${OUT_PATH}/intersection/pdf/ -r ${user_path}/annotation/IRGSPv1/IRGSPv1.rmsk.bed -b ${user_path}/annotation/IRGSPv1/IRGSPv1.blacklist.bed -p ${user_path}/annotation/IRGSPv1/IRGSPv1.gap.bed

        echo ">>>>distance"
        # distance type frequence
        # melt
        if [ $((`cat ${DATA_PATH}/MELT/${simulation_sample}.bed | wc -l`)) -gt '0' ];then
            cat ${DATA_PATH}/MELT/${simulation_sample}.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4":"$9":"$10,$12}' > ${DATA_PATH}/MELT/${simulation_sample}.melt.bed
            grep MELT ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{OFS="\t"}{if($9>$15){g_st=$15}else{g_st=$9};if($10>$16){g_en=$10}else{g_en=$16};  print $1,$2,$3,"GOLD",g_st,g_en,$11":"$12,$14,$18}' | bedtools intersect -a - -b ${DATA_PATH}/MELT/${simulation_sample}.melt.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,"MELT",$11-10,$11-10,$13,$14,"-"}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_MELT.dis.bed
        else
            cat /dev/null > ${OUT_PATH}/intersection/dis/${simulation_sample}_MELT.dis.bed
        fi 
        echo "ck1"
        
        # TEMP2
        grep TEMP2 ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{OFS="\t"}{if($9>$15){g_st=$15}else{g_st=$9};if($10>$16){g_en=$10}else{g_en=$16};  print $1,$2,$3,"GOLD",g_st,g_en,$11":"$12,$14,$18}' | bedtools intersect -a - -b ${DATA_PATH}/TEMP2/${simulation_sample}.insertion.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{split($13,a,":");print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TEMP2",$11+1,$12,a[1]":"a[2]-1":"a[3],$14,"-"}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_TEMP2.dis.bed


        echo "ck2"
        # tldr
        cat ${DATA_PATH}/tldr/${simulation_sample}.table.txt | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$7,$8,$9,$12,$22}' | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g" | sed 1d > ${DATA_PATH}/tldr/${simulation_sample}.table.bed
        grep TLDR-all ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{OFS="\t"}{if($9>$15){g_st=$15}else{g_st=$9};if($10>$16){g_en=$10}else{g_en=$16}; print $1,$2,$3,"GOLD",g_st,g_en,$11":"$12,$14,$18}'| bedtools intersect -a - -b ${DATA_PATH}/tldr/${simulation_sample}.table.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{str1="";split($17,a,""); for(i in a){if(a[i] ~ /[[:lower:]]/){str=a[i];str1=str1""str }}; freq=$16;if($16=="NA"){freq="-"}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TLDR-all",$11+1,$12+1,$13":"$14":"$15,$16,str1}' | grep -v NA > ${OUT_PATH}/intersection/dis/${simulation_sample}_TLDR-all.dis.bed
        grep TLDR-pass ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{OFS="\t"}{if($9>$15){g_st=$15}else{g_st=$9};if($10>$16){g_en=$10}else{g_en=$16}; print $1,$2,$3,"GOLD",g_st,g_en,$11":"$12,$14,$18}'| bedtools intersect -a - -b ${DATA_PATH}/tldr/${simulation_sample}.table.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{str1="";split($17,a,""); for(i in a){if(a[i] ~ /[[:lower:]]/){str=a[i];str1=str1""str }}; freq=$16;if($16=="NA"){freq="-"}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TLDR-pass",$11+1,$12+1,$13":"$14":"$15,$16,str1}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_TLDR-pass.dis.bed

        echo "ck3"

        # TELR
        grep TELR ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{OFS="\t"}{if($9>$15){g_st=$15}else{g_st=$9};if($10>$16){g_en=$10}else{g_en=$16}; print $1,$2,$3,"GOLD",g_st,g_en,$11":"$12,$14,$18}' | bedtools intersect -a - -b ${DATA_PATH}/TELR/${simulation_sample}.telr.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{if($2>0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TELR",$11,$12,$13":0:0",$14,$18}}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_TELR.dis.bed
        echo "ck4"
        
        # TrEMOLO
        grep TrEMOLO ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{OFS="\t"}{if($9>$15){g_st=$15}else{g_st=$9};if($10>$16){g_en=$10}else{g_en=$16};  print $1,$2,$3,"GOLD",g_st,g_en,$11":"$12,$14,$18}' | bedtools intersect -a - -b ${DATA_PATH}/TrEMOLO/${simulation_sample}.TrEMOLO.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{if($2>0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TrEMOLO",$11,$12,$13":0:0",$21,$20}}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_TrEMOLO.dis.bed
        echo "ck5"

        ### GraffiTE
        grep GraffiTE ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{OFS="\t"}{if($9>$15){g_st=$15}else{g_st=$9};if($10>$16){g_en=$10}else{g_en=$16};  print $1,$2,$3,"GOLD",g_st,g_en,$11":"$12,$14,$18}'  | bedtools intersect -a - -b ${DATA_PATH}/GraffiTE/${simulation_sample}.GraffiTE.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{seq=$17;if($17==""){seq="-"};  print $1,$2,$3,$4,$5,$6,$7,$8,$9,"GraffiTE",$11,$12,$13":"0":"0,$19,seq}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_GraffiTE.dis.bed


        # TEMP3
        grep LOCATE-all ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{OFS="\t"}{if($9>$15){g_st=$15}else{g_st=$9};if($10>$16){g_en=$10}else{g_en=$16};  print $1,$2,$3,"GOLD",g_st,g_en,$11":"$12,$14,$18}'  | bedtools intersect -a - -b ${DATA_PATH}/LOCATE/${simulation_sample}.txt -wa -wb | awk -v locate_version=${LOCATE_version} 'BEGIN{FS=OFS="\t"}{split($19,a,",");min_st = 10000;max_end = 0;for(i=1;i<=length(a);i++ ){if(a[i]!~/Poly/){split(a[i],b,":");split(b[2],c,"-");if(c[1]<min_st){min_st=c[1]};if(c[2]>max_end){max_end=c[2]} }};gsub(/N*/,"",$22);print $1,$2,$3,$4,$5,$6,$7,$8,$9,locate_version"-all",$11+1,$12+1,b[1]":"min_st":"max_end,$14,$22}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_${LOCATE_version}-all.dis.bed
        grep LOCATE-pass ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{OFS="\t"}{if($9>$15){g_st=$15}else{g_st=$9};if($10>$16){g_en=$10}else{g_en=$16};  print $1,$2,$3,"GOLD",g_st,g_en,$11":"$12,$14,$18}'  | bedtools intersect -a - -b ${DATA_PATH}/LOCATE/${simulation_sample}.txt -wa -wb | awk -v locate_version=${LOCATE_version} 'BEGIN{FS=OFS="\t"}{split($19,a,",");min_st = 10000;max_end = 0;for(i=1;i<=length(a);i++ ){if(a[i]!~/Poly/){split(a[i],b,":");split(b[2],c,"-");if(c[1]<min_st){min_st=c[1]};if(c[2]>max_end){max_end=c[2]} }};gsub(/N*/,"",$22);print $1,$2,$3,$4,$5,$6,$7,$8,$9,locate_version"-pass",$11+1,$12+1,b[1]":"min_st":"max_end,$14,$22}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_${LOCATE_version}-pass.dis.bed



        prefix=${simulation_sample}
        dis_bed_list=""
        # for meth in ${LOCATE_version}-pass-TLDR-pass TLDR-pass-${LOCATE_version}-pass # TLDR-pass-spe ${LOCATE_version}-pass-TLDR-pass ${LOCATE_version}-pass-spe ${LOCATE_version}-pass TLDR-pass ${LOCATE_version}-all MELT TEMP2 TLDR-all xTea-all xTea-pass
        
        # rm ${OUT_PATH}/intersection/dis/${simulation_sample}.${LOCATE_version}-all.bp_dis.bed
        # rm ${OUT_PATH}/intersection/dis/${simulation_sample}.${LOCATE_version}-pass.bp_dis.bed
        # rm ${OUT_PATH}/intersection/dis/${simulation_sample}.GraffiTE.bp_dis.bed
        rm ${OUT_PATH}/intersection/dis/${simulation_sample}.*.bp_dis.bed
        # rm ${OUT_PATH}/intersection/dis/${simulation_sample}.GraffiTE.bp_dis.bed

        for meth in ${LOCATE_version}-all ${LOCATE_version}-pass TLDR-pass TLDR-all MELT TEMP2 TELR TrEMOLO GraffiTE
        do
            echo $meth 
            # [ -f ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed ] && rm ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed
            if [ "${dis_bed_list}" == "" ];then
                dis_bed_list=${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed
            else
                dis_bed_list=${dis_bed_list}" "${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed
            fi

            # [ -f ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed ] && rm ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed
            
            if [ ! -f ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed ];then
                python "${script_dir}/02_cal_dis_div3.py" ${OUT_PATH}/intersection/dis/${simulation_sample}_${meth}.dis.bed ${user_path}/annotation/IRGSPv1/IRGSPv1.transposon.size  ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed ${simulation_sample} ${meth} ${depth} ${OUT_PATH_dir} IRGSPv1 "${script_dir}/03_cal_div3.sh"
            fi
        done 
        cat ${dis_bed_list} | LC_COLLATE=C sort -k1,1 -k2,2n - | bedtools merge -i - -c 4,5,6 -o first,first,collapse -delim ";" > ${OUT_PATH}/intersection/dis/${simulation_sample}.merge.bp_dis.bed


    done
done

# Rscript ${R_script_path}/intersection_precision_all.R TEMP2,MELT,TLDR,xTea,TEMP3 ${OUT_PATH_dir}/pp_intersection human

# Rscript ${R_script_path}/intersection_TP_FP_count_all.R TEMP2,MELT,TLDR,xTea,TEMP3 ${OUT_PATH_dir}/pp_intersection human
# Rscript ${R_script_path}/performance_dis.R ${OUT_PATH}/intersection/dis/bp_dis.txt ${OUT_PATH_dir}/pp_intersection human


sed -i "s/-/_/g" ${OUT_PATH_dir}/pp_intersection/data/performance*.txt



