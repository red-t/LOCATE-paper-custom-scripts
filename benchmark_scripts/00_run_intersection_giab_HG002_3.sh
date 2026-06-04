#!/bin/bash
start=$(date +%s)
user_path="${user_path:-${USER_PATH:-/path/to/user}}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INSERTION_PATH_dir="${user_path}/2022_long_reads/result/HG002_3"
DATA_PATH="${user_path}/2022_long_reads/result/HG002_3/insertion_revision"
OUT_PATH_dir="${user_path}/2022_long_reads/result/HG002_3/intersection_revision" # intersection_my_bench.r0" # "${user_path}/2022_long_reads/result/giab/intersection" # "${user_path}/2022_long_reads/result/GIAB/benchmark/intersection" # "${user_path}/lrft2/result" # "${user_path}/lrft/result/simulation" 
ANNO_PATH="${user_path}/annotation"
GENOME="hg38"
SAMPLE="HG002"
CPU="15"

R_script_path="${user_path}/2022_long_reads/scripts/"

# intersection
[ -d ${OUT_PATH_dir} ] && rm -r ${OUT_PATH_dir} # ${user_path}/2022_long_reads/result/GIAB/benchmark/tmp/
[ ! -d ${OUT_PATH_dir} ] && mkdir ${OUT_PATH_dir}
[ ! -d ${OUT_PATH_dir}/pp_intersection ] && mkdir ${OUT_PATH_dir}/pp_intersection
[ ! -d ${OUT_PATH_dir}/pp_intersection/data ] && mkdir ${OUT_PATH_dir}/pp_intersection/data
[ ! -d ${OUT_PATH_dir}/pp_intersection/figure ] && mkdir ${OUT_PATH_dir}/pp_intersection/figure

LOCATE_version="LOCATE"


sample="HG002_3"
for depth in 5 10 31 
do

    
    INSERTION_PATH=${INSERTION_PATH_dir}/${depth}

    OUT_PATH=${OUT_PATH_dir}/${depth}
    [ -d ${OUT_PATH}/intersection ] && rm -r ${OUT_PATH}/intersection/*
    [ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}
    [ ! -d ${OUT_PATH}/intersection ] && mkdir ${OUT_PATH}/intersection
    [ ! -d ${OUT_PATH}/intersection/pdf ] && mkdir ${OUT_PATH}/intersection/pdf
    [ ! -d ${OUT_PATH}/intersection/dis ] && mkdir ${OUT_PATH}/intersection/dis

    ##### gold
    # paper_bench from xTea
    # awk -v l=10 -v p=GOLD 'BEGIN{FS=OFS="\t"}{print $1, $2-l, $3+l, $4,p":gold:5:1",$5}' ${DATA_PATH}/paper_bench_set/bench.bed| grep -v ERV  > ${OUT_PATH}/intersection/GOLD.bed
    # awk -v l=10 -v p=GOLD 'BEGIN{FS=OFS="\t"}{print $1, $2-l, $3+l, $4,p":gold:5:1",$5}' ${DATA_PATH}/xTea_paper/bench.bed| grep -v ERV  > ${OUT_PATH}/intersection/GOLD.bed
    # awk -v l=10 -v p=GOLD 'BEGIN{FS=OFS="\t"}{print $1, $2-l, $3+l, $4,p":gold:5:1",$5}' ${DATA_PATH}/xTea_paper/bench_tprt.bed| grep -v ERV  > ${OUT_PATH}/intersection/GOLD.bed

    # my_bench
    # awk -v l=10 -v p=GOLD 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,p":gold:5:1",$16}' ${DATA_PATH}/gold/bench_set_check.bed | grep -v ERV  > ${OUT_PATH}/intersection/GOLD.bed


    awk -v l=10 -v p=GOLD 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,p":gold:5:1",$6}' ${user_path}/2022_long_reads/result/GIAB/benchmark/insertion/gold/bench_set_checked.v3.bed > ${OUT_PATH}/intersection/GOLD.tmp.bed
    bedtools intersect -a ${OUT_PATH}/intersection/GOLD.tmp.bed -b ${user_path}/2022_long_reads/result/GIAB/benchmark/insertion_revision/gold/HG002.insertions.hg38.min100bp.tsv -wa -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$15}' > ${OUT_PATH}/intersection/GOLD.bed
    bedtools intersect -a ${OUT_PATH}/intersection/GOLD.tmp.bed -b ${user_path}/2022_long_reads/result/GIAB/benchmark/insertion_revision/gold/HG002.insertions.hg38.min100bp.tsv -v | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,"-"}' >> ${OUT_PATH}/intersection/GOLD.bed


    # awk -v l=10 -v p=GOLD 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,p":gold:5:1",$5}'  ${user_path}/2022_long_reads/result/GIAB/hg38/gold/paper_bench.1.bed | grep -v ERV  > ${OUT_PATH}/intersection/GOLD.bed
    # awk -v l=10 -v p=GOLD 'BEGIN{FS=OFS="\t"}{split($4,a,";");split($6,b,";");print $1,$2,$3,a[1],p":gold:5:1",b[1]}' ${DATA_PATH}/matrix_v_use/one_sample/verify/gold.ex_del.bed | grep -v ERV   > ${OUT_PATH}/intersection/GOLD.bed
    # cp ${OUT_PATH_dir}/gold/HG002_gold.bed ${OUT_PATH}/intersection/gold.bed

    gold_ins=${OUT_PATH}/intersection/GOLD.bed



    # tldr
    sed 1d ${INSERTION_PATH}/tldr/${sample}_hg38_${depth}.table.txt | grep PASS | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]/|/g" | awk -v l=10 -v p=TLDR 'BEGIN{FS=OFS="\t"}{seq="";split($22,s,"|");for(i=2;i<=length(s)-1;i++){seq=seq""s[i]};print $2,$3-l,$4+l,$6,p":"$6":"$15":"$12,$5}' | grep -v random > ${OUT_PATH}/intersection/TLDR.bed
    tldr_ins=${OUT_PATH}/intersection/TLDR.bed

    # telr
    # bad sensitivity


    # xtea

    cat ${INSERTION_PATH}/xTea/classified_results.txt* | sort | uniq > ${INSERTION_PATH}/xTea/${sample}.xtea.bed
    cat ${INSERTION_PATH}/xTea/${sample}.xtea.bed | awk 'BEGIN{OFS="\t"}{split($5,a,":");strand=a[3]; print $1,$2-5,$2+5,$3, "xtea", a[3], 5, a[2]-a[1],a[1],a[2],$10 }' > ${INSERTION_PATH}/xTea/${sample}.xtea.insertion.bed
    cat ${INSERTION_PATH}/xTea/${sample}.xtea.bed | awk -v l=10 'BEGIN{FS=OFS="\t"}{split($5,a,":");split($2,b," ");if(a[2]>=0 && a[1]>=0 ){print $1,b[1]-l,b[1]+1+l,$3,"xtea_2",a[3],"5",a[2]-a[1],a[1],a[2],$10}}' \
        | LC_COLLATE=C sort -k1,1 -k2,2n - |  bedtools merge -i - -c 4,5,6,7,8,9,10,11 -o  distinct_only |  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,5,$8,$9,$10,$11}' \
        | awk -v l=10 -v p=xTea 'BEGIN{FS=OFS="\t"}{print $1,$2-l,$3+l,$4,p":germ:5:1",$6}' > ${OUT_PATH}/intersection/xTea.bed

    # cat ${DATA_PATH}/xtea/xtea_ccs.bed | awk -v l=5 -v p=xTea 'BEGIN{FS=OFS="\t"}{print $1,$2-l,$3+l,$4,p":germ:5:1",$6}' > ${OUT_PATH}/intersection/xTea.bed
    xtea_ins=${OUT_PATH}/intersection/xTea.bed




    #### LOCATE

    sed 1d ${INSERTION_PATH}/${LOCATE_version}/${sample}/result.tsv | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' >  ${INSERTION_PATH}/LOCATE/${sample}.txt
    # sed 1d ${INSERTION_PATH}/${LOCATE_version}/${sample}/result.tsv | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14,$15,$16,$17}' >  ${DATA_PATH}/LOCATE/${simulation_sample}.txt

    # awk '{if( $8=="True" ){print $0}}' ${INSERTION_PATH}/LOCATE/${sample}.txt | awk 'BEGIN{FS=OFS="\t"}{ gsub(/N*/,"",$13); print $1,$2+1,$3+1,$4,"LOCATE-pass:germ:"$11":"$5":"$13,$6}'  > ${OUT_PATH}/intersection/LOCATE.bed
    awk '{if( $8=="True" ){print $0}}' ${INSERTION_PATH}/LOCATE/${sample}.txt | awk -v ns=${ns_num}  'BEGIN{FS=OFS="\t"}{if($11>ns){gsub(/N*/,"",$13); if($7=="1/1"){geno=1}else{if($7=="0/1"){geno=0.5}else{geno=0}}; print $1,$2+1,$3+1,$4,"LOCATE:germ:"$11":"geno":"$13, $6}}'  > ${OUT_PATH}/intersection/LOCATE.bed

    LOCATE_ins=${OUT_PATH}/intersection/LOCATE.bed

            

    echo "ck"
    # TrEMOLO
    bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/RE\t%ALT\n"  ${INSERTION_PATH}/TrEMOLO/OUTSIDER/VARIANT_CALLING/SV.vcf | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"sniffles."$5"."$4,$6,$7}' > ${INSERTION_PATH}/TrEMOLO/OUTSIDER/VARIANT_CALLING/SV.supp.bed
    sed 1d ${INSERTION_PATH}/TrEMOLO/OUTSIDER/TrEMOLO_SV_TE/INS/INS_TREMOLO.csv | awk 'BEGIN{OFS="\t"}{split($2,a,":");print a[1],a[3],a[4],a[5],a[6],a[9]}' >> ${INSERTION_PATH}/TrEMOLO/OUTSIDER/VARIANT_CALLING/SV.supp.bed
    cat ${INSERTION_PATH}/TrEMOLO/OUTSIDER/ET_FIND_FA/*fasta  | awk '{if($1~/>/){split($1,a,":");s=a[5]}else{print s, $1;s=""}}' > ${INSERTION_PATH}/TrEMOLO/OUTSIDER/ET_FIND_FA/TE_INS.bed
    awk 'BEGIN{OFS="\t"}NR==FNR{a[$4]=$5}NR>FNR{split($4,b,"|");print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,a[b[2]] }' ${INSERTION_PATH}/TrEMOLO/OUTSIDER/VARIANT_CALLING/SV.supp.bed ${INSERTION_PATH}/TrEMOLO/TE_INFOS.bed | grep -v "DEL" > ${INSERTION_PATH}/TrEMOLO/TE_INFOS.supp.bed

    grep chr ${INSERTION_PATH}/TrEMOLO/TE_INFOS.supp.bed | sed 1d | awk 'BEGIN{OFS="\t"}{if($7>0){split($4,a,"|");if($2>$10){st=$10;en=$3}else{st=$2;en=$10};  print $1,st,en,a[1],$11,$5,$8,$9,1,$10,a[2]}}' | awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2}NR>FNR{if($11 in a){seq=a[$11]}else{seq="-"}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,seq}' ${INSERTION_PATH}/TrEMOLO/OUTSIDER/ET_FIND_FA/TE_INS.bed - > ${INSERTION_PATH}/TrEMOLO/TrEMOLO.bed
    awk -v l=10 -v p=TrEMOLO -v ns=3 'BEGIN{OFS="\t"}{split($4,a,"|"); if($2>$10){st=$10;en=$3}else{st=$2;en=$10};if($5>ns){print $1,st-l,en+l,a[1],p":INS:"$5":"1,$6}}' ${INSERTION_PATH}/TrEMOLO/TrEMOLO.bed > ${OUT_PATH}/intersection/TrEMOLO.bed
    TrEMOLO_ins=${OUT_PATH}/intersection/TrEMOLO.bed



    ### GraffiTE

    bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/repeat_ids\t%INFO/SUPPORT\t%INFO/RM_hit_strands\t%ID\t%ALT\n" ${INSERTION_PATH}/GraffiTE/2_Repeat_Filtering/genotypes_repmasked_filtered.vcf > ${INSERTION_PATH}/GraffiTE/GraffiTE.tmp.bed
    # awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$9"_"$11"_"length($12)"_"$12}NR>FNR{if($7 in a){split(a[$7],b,"_"); if(b[1]=="0"){st=$2-b[3];en=$2}else{st=$2;en=$2+b[3]} }else{st=$2;en=$2+1}; split($4,c,",");split($6,d,",");if(d[1]=="C"){strand="-"}else{strand="+"}; print $1,st,en,c[1],$5,strand,$7,$8,a[$7] }' ${INSERTION_PATH}/GraffiTE/${simulation_sample}/3_TSD_search/TSD_summary.txt  ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/GraffiTE.tmp.bed > ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/${simulation_sample}_GraffiTE.bed
    awk 'BEGIN{OFS="\t"}NR==FNR{if($14=="PASS"){a[$1]=$9"_"$11"_"length($12)"_"$12}}
                        NR>FNR{if($7 in a){split(a[$7],b,"_"); st=$2;en=$2+b[3];TSD=b[4];split($8,c,TSD);if(c[1]!=$8){ins_st=length(c[1])+b[3]+1;ins_len=length($8)-ins_st+1; ins_seq=substr($8,ins_st,ins_len)}else{ins_seq=substr($8,1,length(c[1]))} }
                                else{a[$7]=".";st=$2;en=$2+1;ins_seq=$8}; split($4,c,",");split($6,d,","); 
                                    if(d[1]=="C"){strand="-"}else{strand="+"}; 
                                print $1,st,en,c[1],$5,strand,$7,ins_seq,a[$7] }' ${INSERTION_PATH}/GraffiTE/3_TSD_search/TSD_summary.txt  ${INSERTION_PATH}/GraffiTE/GraffiTE.tmp.bed | grep -v "DEL" > ${INSERTION_PATH}/GraffiTE/2_Repeat_Filtering/GraffiTE.bed

    cp ${INSERTION_PATH}/GraffiTE/2_Repeat_Filtering/GraffiTE.bed ${INSERTION_PATH}/GraffiTE/GraffiTE.bed
    awk -v l=10 -v p=GraffiTE -v ns=0 'BEGIN{OFS="\t"}{if($5>ns){print $1,$2-l,$3+l,$4,p":INS:"$5":"1":"$8,$6}}' ${INSERTION_PATH}/GraffiTE/GraffiTE.bed > ${OUT_PATH}/intersection/GraffiTE.bed
    GraffiTE_ins=${OUT_PATH}/intersection/GraffiTE.bed


    echo ">>>>intersection"

    #### intersectte
    "${script_dir}/01_intersect_gold.sh" -i ${tldr_ins},${xtea_ins},${LOCATE_ins},${TrEMOLO_ins},${GraffiTE_ins} -g ${gold_ins} -d ${depth} -t giab -c ${OUT_PATH_dir}/pp_intersection/data -o ${OUT_PATH}/intersection/pdf/ -r ${user_path}/annotation/hg38/ALSE.bed -b /zata/zippy/boxu/annotation/from_temp3/BlackList.bed -p /zata/zippy/boxu/annotation/from_temp3/gap.bed



    # exit

    echo ">>>>distance"
    ## distance type frequence
    ## melt


    # tldr

    cat ${INSERTION_PATH}/tldr/${sample}_hg38_${depth}.table.txt | awk 'BEGIN{FS=OFS="\t"}{print $2,$3+1,$4+1,$7,$8,$9,$12,$22,$5}' | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g" > ${INSERTION_PATH}/tldr/${sample}_tgs_tldr.table.bed # ${simulation_sample}.table.bed
    # grep TLDR ${OUT_PATH}/intersection/pdf/giab_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${OUT_PATH}/intersection/GOLD.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"GOLD",$9,$10,$11":"$12,1,"-"}' | bedtools intersect -a - -b ${DATA_PATH}/tldr/${sample}_tgs_tldr.table.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{str1="";split($20,a,""); for(i in a){if(a[i] ~ /[[:lower:]]/){str=a[i];str1=str1""str }}; if(str1==""){str1="-"}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TLDR-all",$11,$12,$13":"$14":"$15,$16,str1}' > ${OUT_PATH}/intersection/dis/${sample}_tgs_TLDR-all.dis.bed
    grep TLDR ${OUT_PATH}/intersection/pdf/giab_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${OUT_PATH}/intersection/GOLD.bed -wa -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"GOLD",$9,$10,0":"0"~"$11,1,$13}' | bedtools intersect -a - -b ${INSERTION_PATH}/tldr/${sample}_tgs_tldr.table.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{split($17,a,""); str1=""; for(i in a){if(a[i] ~ /[[:lower:]]/){str=a[i];str1=str1""str }}; if (str1=="" ){str1="-"}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TLDR",$11,$12,$13":"$14":"$15,$16,str1}' | grep -v "NA:NA:NA" > ${OUT_PATH}/intersection/dis/${sample}_tgs_TLDR.dis.bed
    echo "ck3"

    # xTea
    # cat ${DATA_PATH}/xtea/${sample}_tgs_xtea_d5.bed | sort | uniq | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$2,$3":"$5,$10}'  > ${DATA_PATH}/xtea/${sample}_tgs_xtea_d5.table.bed
    # grep xTea-all ${OUT_PATH}/intersection/pdf/giab_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${OUT_PATH}/intersection/GOLD.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"GOLD",$9,$10,$11":"$12,1,"-"}' | bedtools intersect -a - -b ${DATA_PATH}/xtea/${sample}_tgs_xtea_d5.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"xTea-all",$11,$12,$13,0,$20}' > ${OUT_PATH}/intersection/dis/${sample}_tgs_xTea-all.dis.bed
    grep xTea ${OUT_PATH}/intersection/pdf/giab_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${OUT_PATH}/intersection/GOLD.bed -wa -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"GOLD",$9,$10,0":"0"~"$11,1,$13}' | bedtools intersect -a - -b ${INSERTION_PATH}/xTea/${sample}.xtea.insertion.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"xTea",$11,$12,$13":"$16":"$17,0,$20}' > ${OUT_PATH}/intersection/dis/${sample}_tgs_xTea.dis.bed

    echo "check4"


    # LOCATE
    # grep LOCATE-all ${OUT_PATH}/intersection/pdf/giab_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${OUT_PATH}/intersection/GOLD.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"GOLD",$9,$10,$11":"$12,1,"-"}' | bedtools intersect -a - -b ${DATA_PATH}/LOCATE/${sample}.txt -wa -wb | awk -v LOCATE_version=${LOCATE_version} 'BEGIN{FS=OFS="\t"}{split($18,a,",");min_st = 10000;max_end = 0;for(i=1;i<=length(a);i++ ){if(a[i]!~/Poly/){split(a[i],b,":");split(b[2],c,"-");if(c[1]<min_st){min_st=c[1]};if(c[2]>max_end){max_end=c[2]} }};gsub(/N*/,"",$25);print $1,$2,$3,$4,$5,$6,$7,$8,$9,LOCATE_version"-all",$11+1,$12+1,b[1]":"min_st":"max_end,0,$25}' > ${OUT_PATH}/intersection/dis/${sample}_tgs_${LOCATE_version}-all.dis.bed
    grep LOCATE ${OUT_PATH}/intersection/pdf/giab_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${OUT_PATH}/intersection/GOLD.bed -wa -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"GOLD",$9,$10,0":"0"~"$11,1,$13}' | bedtools intersect -a - -b ${INSERTION_PATH}/LOCATE/${sample}.txt -wa -wb | awk -v LOCATE_version=${LOCATE_version} 'BEGIN{FS=OFS="\t"}{split($19,a,",");min_st = 10000;max_end = 0;for(i=1;i<=length(a);i++ ){if(a[i]!~/Poly/){split(a[i],b,":");split(b[2],c,"-");if(c[1]<min_st){min_st=c[1]};if(c[2]>max_end){max_end=c[2]} }};gsub(/N*/,"",$22);print $1,$2,$3,$4,$5,$6,$7,$8,$9,LOCATE_version,$11+1,$12+1,b[1]":"min_st":"max_end,$14,$22}' > ${OUT_PATH}/intersection/dis/${sample}_tgs_${LOCATE_version}.dis.bed   



    prefix=${simulation_sample}
    dis_bed_list=""
    # for meth in ${LOCATE_version}-pass-TLDR-pass TLDR-pass-${LOCATE_version}-pass # TLDR-pass-spe ${LOCATE_version}-pass-TLDR-pass ${LOCATE_version}-pass-spe ${LOCATE_version}-pass TLDR-pass ${LOCATE_version}-all MELT TEMP2 TLDR-all xTea-all xTea-pass

    # rm ${OUT_PATH}/intersection/dis/${simulation_sample}.${LOCATE_version}-all.bp_dis.bed
    # rm ${OUT_PATH}/intersection/dis/${simulation_sample}.${LOCATE_version}-pass.bp_dis.bed
    # rm ${OUT_PATH}/intersection/dis/${simulation_sample}.GraffiTE.bp_dis.bed
    # rm ${OUT_PATH}/intersection/dis/${simulation_sample}.TrEMOLO.bp_dis.bed
    # rm ${OUT_PATH}/intersection/dis/${simulation_sample}.GraffiTE.bp_dis.bed


    for meth in ${LOCATE_version} TLDR xTea  # GraffiTE PALMER MEHunter TrEMOLO MELT TEMP2  # MELT TEMP2 xTea-all xTea-pass PALMER TELR TrEMOLO MEHunter
    do
        echo $meth 
        # [ -f ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed ] && rm ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed
        if [ "${dis_bed_list}" == "" ];then
            dis_bed_list=${OUT_PATH}/intersection/dis/${sample}_tgs.${meth}.bp_dis.bed
        else
            dis_bed_list=${dis_bed_list}" "${OUT_PATH}/intersection/dis/${sample}_tgs.${meth}.bp_dis.bed
        fi

        # [ -f ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed ] && rm ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed

        if [ ! -f ${OUT_PATH}/intersection/dis/${sample}_tgs.${meth}.bp_dis.bed ];then
            python "${script_dir}/02_cal_dis_div3_giab.py" ${OUT_PATH}/intersection/dis/${sample}_tgs_${meth}.dis.bed /zata/zippy/boxu/annotation/hg38/ALSE.size ${OUT_PATH}/intersection/dis/${sample}_tgs.${meth}.bp_dis.bed ${sample}_tgs ${meth} ${depth} ${OUT_PATH_dir} hg38 "${script_dir}/03_cal_div3_giab.sh"
        fi
    done 
    cat ${dis_bed_list} | LC_COLLATE=C sort -k1,1 -k2,2n - | bedtools merge -i - -c 4,5,6 -o first,first,collapse -delim ";" > ${OUT_PATH}/intersection/dis/giab.merge.bp_dis.bed

done


# "${script_dir}/01_intersect_gold.sh" -i ${temp2_ins},${melt_ins},${tldr_ins},${xtea_ins},${temp3_ins} -g ${gold_ins} -d ${depth} -t ${simulation_sample} -c ${OUT_PATH_dir}/pp_intersection/data -o ${OUT_PATH}/intersection/pdf/ -r ${user_path}/annotation/hg38/ALSE.bed -b /zata/zippy/boxu/annotation/from_temp3/BlackList.bed -p /zata/zippy/boxu/annotation/from_temp3/gap.bed

# prefix=bench_set
# # Rscript ${R_script_path}/intersection_precision_giab.R MELT,TEMP2,TLDR,xTea,TEMP3 ${OUT_PATH_dir}/pp_intersection ${prefix}
# Rscript ${R_script_path}/intersection_precision_all.R TEMP2,MELT,TLDR,xTea,TEMP3,MEHunter ${OUT_PATH_dir}/pp_intersection ${prefix}

# R --slave --no-restore --file=${R_script_path}/intersection_precision_one_depth.R --args MELT,TEMP2,SMS,TLDR,xTea,TEMP3 ${OUT_PATH_dir}/pp_intersection ${prefix}
# Rscript ${R_script_path}/intersection_upset.R ${OUT_PATH}/intersection/pdf/matrix.${depth}.txt ${OUT_PATH_dir}/pp_intersection/figure/${prefix}_upset.pdf MELT,TEMP2,TLDR,xTea,TEMP3,GOLD ${OUT_PATH_dir}/pp_intersection ${prefix}
# Rscript ${R_script_path}/intersection_venn.R ${OUT_PATH}/intersection/pdf/validation.intersect.txt ${OUT_PATH_dir}/pp_intersection/figure/venn_giab_hg38.pdf MELT,TEMP2,SMS,TLDR,xTea,TEMP3,GOLD # GOLD,SMS,xTea #  GOLD,TLDR,PALMER,SMS,xTea
