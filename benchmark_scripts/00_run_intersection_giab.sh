#!/bin/bash
start=$(date +%s)
user_path="${user_path:-${USER_PATH:-/path/to/user}}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DATA_PATH="${user_path}/2022_long_reads/result/GIAB/benchmark/insertion_revision"
OUT_PATH_dir="${user_path}/2022_long_reads/result/giab/intersection/intersection_revision" # intersection_my_bench.r0" # "${user_path}/2022_long_reads/result/giab/intersection" # "${user_path}/2022_long_reads/result/GIAB/benchmark/intersection" # "${user_path}/lrft2/result" # "${user_path}/lrft/result/simulation" 
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

dep=giab
depth=giab
OUT_PATH=${OUT_PATH_dir}
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
bedtools intersect -a ${OUT_PATH}/intersection/GOLD.tmp.bed -b ${user_path}/2022_long_reads/result/GIAB/benchmark/insertion_revision/gold/HG002.insertions.hg38.min100bp.tsv -v >> ${OUT_PATH}/intersection/GOLD.bed


# awk -v l=10 -v p=GOLD 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,p":gold:5:1",$5}'  ${user_path}/2022_long_reads/result/GIAB/hg38/gold/paper_bench.1.bed | grep -v ERV  > ${OUT_PATH}/intersection/GOLD.bed
# awk -v l=10 -v p=GOLD 'BEGIN{FS=OFS="\t"}{split($4,a,";");split($6,b,";");print $1,$2,$3,a[1],p":gold:5:1",b[1]}' ${DATA_PATH}/matrix_v_use/one_sample/verify/gold.ex_del.bed | grep -v ERV   > ${OUT_PATH}/intersection/GOLD.bed
# cp ${OUT_PATH_dir}/gold/HG002_gold.bed ${OUT_PATH}/intersection/gold.bed

gold_ins=${OUT_PATH}/intersection/GOLD.bed

### MELT
RESULT_PATH="${user_path}/2022_long_reads/result/GIAB/hg38/melt/HG002_2"
[ -f ${RESULT_PATH}/melt.vcf.bed ] && rm ${RESULT_PATH}/melt.vcf.bed

for melt_te in ALU LINE1 SVA HERVK;do
    if [ $((`cat ${RESULT_PATH}/${melt_te}.final_comp.vcf | wc -l`)) -gt '0' ];then
        bcftools query -f "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/LP\t%INFO/RP\t%INFO/SVLEN\t%INFO/TSD\t%INFO/MEINFO\t%INFO/SVLEN\n" ${RESULT_PATH}/${melt_te}.final_comp.vcf | awk -v l=10 'BEGIN{FS=OFS="\t"}{split($8,a,",");print $1,$2-l,$2+l,$3,"melt_2",a[4],$4+$5,$9,a[2],a[3],"-",$10}' >> ${RESULT_PATH}/melt.vcf.bed
    fi
done
cat ${RESULT_PATH}/*master.bed > ${RESULT_PATH}/melt.master1.bed

if [ ! -f ${RESULT_PATH}/melt.vcf.bed ];then
    cat /dev/null > ${RESULT_PATH}/melt.vcf.bed
    cp ${RESULT_PATH}/melt.vcf.bed ${DATA_PATH}/melt/HG002_2_ngs_melt.bed
else
    cat ${RESULT_PATH}/melt.vcf.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools merge -i - -o first -c 4,5,6,7,8,9,10,11 > ${DATA_PATH}/melt/HG002_2_ngs_melt.bed
fi

cat ${DATA_PATH}/melt/HG002_2_ngs_melt.bed  | grep -v ERV  | awk -v l=10 -v p=MELT 'BEGIN{FS=OFS="\t"}{print $1,$2-l,$3+l,$4,p":germ:5:1",$6}' > ${OUT_PATH}/intersection/MELT.bed
melt_ins=${OUT_PATH}/intersection/MELT.bed


# TEMP2
sed 1d ${DATA_PATH}/TEMP2/HG002_2_ngs_TEMP2_hg38X.insertion.bed | grep -v ERV | awk -v l=10 -v p=TEMP2 -v ns=2 'BEGIN{FS=OFS="\t"}{split($4,a,":");if($8>ns){print $1,$2-l,$3+l,a[1],p":"$7":"$8":"$5,$6}}' > ${OUT_PATH}/intersection/TEMP2.bed
temp2_ins=${OUT_PATH}/intersection/TEMP2.bed



# tldr
sed 1d ${DATA_PATH}/tldr/HG002_2_tgs_tldr_hg38X.insertion.bed | grep PASS | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]/|/g" | awk -v l=10 -v p=TLDR 'BEGIN{FS=OFS="\t"}{seq="";split($22,s,"|");for(i=2;i<=length(s)-1;i++){seq=seq""s[i]};print $2,$3-l,$4+l,$6,p":"$6":"$15":"$12,$5}'  > ${OUT_PATH}/intersection/TLDR.bed
tldr_ins=${OUT_PATH}/intersection/TLDR.bed

# telr
# bad sensitivity


# PALMER
#### PALMER
# cp ${user_path}/2022_long_reads/result/GIAB/hg38/xTea_paper/run_tools/Palmer/* ${DATA_PATH}/PALMER/
[ -f ${DATA_PATH}/PALMER/HG002_2_tgs_palmer_hg38X.tmp.bed ] && rm ${DATA_PATH}/PALMER/HG002_2_tgs_palmer_hg38X.tmp.bed 
for palmerTE in ALU LINE1 SVA
do
    # cat ${INSERTION_PATH}/PALMER/${simulation_sample}/${palmerTE}/${simulation_sample}.50X.${palmerTE}_calls.txt >> ${INSERTION_PATH}/PALMER/${simulation_sample}/${simulation_sample}.50X.palmer.bed
    if [ -f ${DATA_PATH}/PALMER/HG002_ccs_Palmer_${palmerTE}_calls.txt ];then
        cat ${DATA_PATH}/PALMER/HG002_ccs_Palmer_${palmerTE}_calls.txt | awk -v te=${palmerTE} -v dep=${depth} 'BEGIN{OFS="\t"}{if($10>=1){print $1,$2,$3,te,$10,$13,$6,$8,$14,$15,$16}}'   >> ${DATA_PATH}/PALMER/HG002_2_tgs_palmer_hg38X.tmp.bed
    fi
done
# awk -v te=${palmerTE} 'BEGIN{OFS="\t"}{print $2,$3,$4,te,$11,$14,$7,$9,$15,$16,$17}' ${INSERTION_PATH}/PALMER/${simulation_sample}/${simulation_sample}.50X.palmer.bed | sed 1d >> ${DATA_PATH}/PALMER/${simulation_sample}.palmer.bed # PALMER/${simulation_sample}/${simulation_sample}.50X.palmer.bed
# cp ${INSERTION_PATH}/PALMER/${simulation_sample}/${simulation_sample}.50X.palmer.bed ${DATA_PATH}/PALMER/${simulation_sample}.palmer.bed
LC_COLLATE=C sort -k1,1 -k2,2n ${DATA_PATH}/PALMER/HG002_2_tgs_palmer_hg38X.tmp.bed | bedtools merge -i - -c 4,5,6,7,8,9,10,11 -d 5 -o first  > ${DATA_PATH}/PALMER/HG002_2_tgs_palmer_hg38X.bed

cat ${DATA_PATH}/PALMER/HG002_2_tgs_palmer_hg38X.bed | awk -v l=10 -v p=PALMER -v ns=1 'BEGIN{FS=OFS="\t"}{if($5>ns){print $1,$2-l,$3+l,$4,p":germ:"$5":1",$6}}' > ${OUT_PATH}/intersection/PALMER.bed
# cat ${user_path}/2022_long_reads/result/GIAB/hg38/PALMER/HG002/HG002_2_tgs_palmer_hg38X.bed | awk -v l=10 -v p=PALMER -v ns=1 'BEGIN{FS=OFS="\t"}{if($5>ns){print $1,$2-l,$3+l,$4,p":germ:"$5":1",$6}}' > ${OUT_PATH}/intersection/PALMER.bed

palmer_ins=${OUT_PATH}/intersection/PALMER.bed



# xtea
if [ "0" == "0" ];then
    cat ${DATA_PATH}/xtea/HG002_2_tgs_xtea.bed  | grep -v ERV  | awk -v l=10 -v p=xTea 'BEGIN{FS=OFS="\t"}{print $1,$2-l,$3+l,$4,p":germ:5:1",$6}' > ${DATA_PATH}/xtea/HG002_2_tgs_xtea.1.bed
    bedtools intersect -a ${DATA_PATH}/xtea/HG002_2_tgs_xtea.1.bed -b ${DATA_PATH}/xtea/candidate_list_n_supp.txt -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,"xTea:"$4":"$11":1",$6}' > ${OUT_PATH}/intersection/xTea.bed
    bedtools intersect -a ${DATA_PATH}/xtea/HG002_2_tgs_xtea.1.bed -b ${DATA_PATH}/xtea/candidate_list_n_supp.txt -v | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,"xTea:"$4":"5":1",$6}' >> ${OUT_PATH}/intersection/xTea.bed
else
    cat ${DATA_PATH}/xtea/HG002_2_tgs_xtea_d5.bed | grep -v ERV  | awk -v l=10 -v p=xTea 'BEGIN{FS=OFS="\t"}{print $1,$2-l,$3+l,$4,p":germ:5:1",$6}' > ${OUT_PATH}/intersection/xTea.bed
fi

cat ${DATA_PATH}/xtea/HG002_2_tgs_xtea_hg38X.insertion.bed | sort | uniq | awk -v l=10 'BEGIN{FS=OFS="\t"}{split($5,a,":");split($2,b," ");if(a[2]>=0 && a[1]>=0 ){print $1,b[1]-l,b[1]+1+l,$3,"xtea_2",a[3],"5",a[2]-a[1],a[1],a[2],$10}}' \
     | LC_COLLATE=C sort -k1,1 -k2,2n - |  bedtools merge -i - -c 4,5,6,7,8,9,10,11 -o  distinct_only  | bedtools intersect -a - -b ${DATA_PATH}/xtea/output_hg38/HG002_2_hg38/candidate_list_n_supp.txt -wa -wb |  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$15,$8,$9,$10,$11}' \
     | awk -v l=10 -v p=xTea 'BEGIN{FS=OFS="\t"}{print $1,$2-l,$3+l,$4,p":germ:5:1",$6}' > ${OUT_PATH}/intersection/xTea.bed

# cat ${DATA_PATH}/xtea/xtea_ccs.bed | awk -v l=5 -v p=xTea 'BEGIN{FS=OFS="\t"}{print $1,$2-l,$3+l,$4,p":germ:5:1",$6}' > ${OUT_PATH}/intersection/xTea.bed
xtea_ins=${OUT_PATH}/intersection/xTea.bed


#### LOCATE
awk '{if(and($19,1)){print $0}}' ${DATA_PATH}/LOCATE/HG002_2.txt | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,"LOCATE:germ:"$10":1", $6}' > ${OUT_PATH}/intersection/LOCATE.bed
LOCATE_ins=${OUT_PATH}/intersection/LOCATE.bed

echo "ck"
# TrEMOLO
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/RE\t%ALT\n"  ${DATA_PATH}/TrEMOLO/OUTSIDER/VARIANT_CALLING/SV.vcf | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"sniffles."$5"."$4,$6,$7}' > ${DATA_PATH}/TrEMOLO/OUTSIDER/VARIANT_CALLING/SV.supp.bed
sed 1d ${DATA_PATH}/TrEMOLO/OUTSIDER/TrEMOLO_SV_TE/INS/INS_TREMOLO.csv | awk 'BEGIN{OFS="\t"}{split($2,a,":");print a[1],a[3],a[4],a[5],a[6],a[9]}' >> ${DATA_PATH}/TrEMOLO/OUTSIDER/VARIANT_CALLING/SV.supp.bed
cat ${DATA_PATH}/TrEMOLO/OUTSIDER/ET_FIND_FA/*fasta  | awk '{if($1~/>/){split($1,a,":");s=a[5]}else{print s, $1;s=""}}' > ${DATA_PATH}/TrEMOLO/OUTSIDER/ET_FIND_FA/TE_INS.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$4]=$5}NR>FNR{split($4,b,"|");print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,a[b[2]] }' ${DATA_PATH}/TrEMOLO/OUTSIDER/VARIANT_CALLING/SV.supp.bed ${DATA_PATH}/TrEMOLO/TE_INFOS.bed | grep -v "DEL" > ${DATA_PATH}/TrEMOLO/TE_INFOS.supp.bed

grep chr ${DATA_PATH}/TrEMOLO/TE_INFOS.supp.bed | sed 1d | awk 'BEGIN{OFS="\t"}{if($7>0){split($4,a,"|");if($2>$10){st=$10;en=$3}else{st=$2;en=$10};  print $1,st,en,a[1],$11,$5,$8,$9,1,$10,a[2]}}' | awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2}NR>FNR{if($11 in a){seq=a[$11]}else{seq="-"}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,seq}' ${DATA_PATH}/TrEMOLO/OUTSIDER/ET_FIND_FA/TE_INS.bed - > ${DATA_PATH}/TrEMOLO/TrEMOLO.bed
awk -v l=10 -v p=TrEMOLO -v ns=3 'BEGIN{OFS="\t"}{split($4,a,"|"); if($2>$10){st=$10;en=$3}else{st=$2;en=$10};if($5>ns){print $1,st-l,en+l,a[1],p":INS:"$5":"1,$6}}' ${DATA_PATH}/TrEMOLO/TrEMOLO.bed > ${OUT_PATH}/intersection/TrEMOLO.bed
TrEMOLO_ins=${OUT_PATH}/intersection/TrEMOLO.bed

### MEHunter
sed "s/'/ /g" ${DATA_PATH}/MEHunter/High_Quality_MEIs.txt | awk 'BEGIN{OFS="\t"}{print $4,$6,$6,$12,$10,".",$8}' > ${DATA_PATH}/MEHunter/MEHunter.bed
awk -v l=10 -v p=MEHunter  'BEGIN{OFS="\t"}{print $1, $2-l, $3+l, $4, p":INS:"$5":1", $6}' ${DATA_PATH}/MEHunter/MEHunter.bed > ${OUT_PATH}/intersection/MEHunter.bed
MEHunter_ins=${OUT_PATH}/intersection/MEHunter.bed


### GraffiTE

bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/repeat_ids\t%INFO/SUPPORT\t%INFO/RM_hit_strands\t%ID\t%ALT\n" ${DATA_PATH}/GraffiTE/2_Repeat_Filtering/genotypes_repmasked_filtered.vcf > ${DATA_PATH}/GraffiTE/GraffiTE.tmp.bed
# awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$9"_"$11"_"length($12)"_"$12}NR>FNR{if($7 in a){split(a[$7],b,"_"); if(b[1]=="0"){st=$2-b[3];en=$2}else{st=$2;en=$2+b[3]} }else{st=$2;en=$2+1}; split($4,c,",");split($6,d,",");if(d[1]=="C"){strand="-"}else{strand="+"}; print $1,st,en,c[1],$5,strand,$7,$8,a[$7] }' ${INSERTION_PATH}/GraffiTE/${simulation_sample}/3_TSD_search/TSD_summary.txt  ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/GraffiTE.tmp.bed > ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/${simulation_sample}_GraffiTE.bed
awk 'BEGIN{OFS="\t"}NR==FNR{if($14=="PASS"){a[$1]=$9"_"$11"_"length($12)"_"$12}}
                    NR>FNR{if($7 in a){split(a[$7],b,"_"); st=$2;en=$2+b[3];TSD=b[4];split($8,c,TSD);if(c[1]!=$8){ins_st=length(c[1])+b[3]+1;ins_len=length($8)-ins_st+1; ins_seq=substr($8,ins_st,ins_len)}else{ins_seq=substr($8,1,length(c[1]))} }
                            else{a[$7]=".";st=$2;en=$2+1;ins_seq=$8}; split($4,c,",");split($6,d,","); 
                                if(d[1]=="C"){strand="-"}else{strand="+"}; 
                            print $1,st,en,c[1],$5,strand,$7,ins_seq,a[$7] }' ${DATA_PATH}/GraffiTE/3_TSD_search/TSD_summary.txt  ${DATA_PATH}/GraffiTE/GraffiTE.tmp.bed | grep -v "DEL" > ${DATA_PATH}/GraffiTE/2_Repeat_Filtering/GraffiTE.bed

cp ${DATA_PATH}/GraffiTE/2_Repeat_Filtering/GraffiTE.bed ${DATA_PATH}/GraffiTE/GraffiTE.bed
awk -v l=10 -v p=GraffiTE -v ns=0 'BEGIN{OFS="\t"}{if($5>ns){print $1,$2-l,$3+l,$4,p":INS:"$5":"1":"$8,$6}}' ${DATA_PATH}/GraffiTE/GraffiTE.bed > ${OUT_PATH}/intersection/GraffiTE.bed
GraffiTE_ins=${OUT_PATH}/intersection/GraffiTE.bed


echo ">>>>intersection"

#### intersectte
"${script_dir}/01_intersect_gold.sh" -i ${melt_ins},${temp2_ins},${tldr_ins},${xtea_ins},${LOCATE_ins},${MEHunter_ins},${TrEMOLO_ins},${GraffiTE_ins},${palmer_ins} -g ${gold_ins} -d ${dep} -t giab -c ${OUT_PATH_dir}/pp_intersection/data -o ${OUT_PATH}/intersection/pdf/ -r ${user_path}/annotation/hg38/ALSE.bed -b /zata/zippy/boxu/annotation/from_temp3/BlackList.bed -p /zata/zippy/boxu/annotation/from_temp3/gap.bed


echo ">>>>distance"
## distance type frequence
## melt


# tldr
cat ${DATA_PATH}/tldr/HG002_2_tgs_tldr_hg38X.insertion.bed | awk 'BEGIN{FS=OFS="\t"}{print $2,$3+1,$4+1,$7,$8,$9,$12,$22,$5}' | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g" > ${DATA_PATH}/tldr/HG002_2_tgs_tldr.table.bed # ${simulation_sample}.table.bed
# grep TLDR ${OUT_PATH}/intersection/pdf/giab_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${OUT_PATH}/intersection/GOLD.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"GOLD",$9,$10,$11":"$12,1,"-"}' | bedtools intersect -a - -b ${DATA_PATH}/tldr/HG002_2_tgs_tldr.table.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{str1="";split($20,a,""); for(i in a){if(a[i] ~ /[[:lower:]]/){str=a[i];str1=str1""str }}; if(str1==""){str1="-"}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TLDR-all",$11,$12,$13":"$14":"$15,$16,str1}' > ${OUT_PATH}/intersection/dis/HG002_2_tgs_TLDR-all.dis.bed
grep TLDR ${OUT_PATH}/intersection/pdf/giab_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${OUT_PATH}/intersection/GOLD.bed -wa -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"GOLD",$9,$10,0":"0"~"$11,1,$13}' | bedtools intersect -a - -b ${DATA_PATH}/tldr/HG002_2_tgs_tldr.table.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{str1="-";split($17,a,""); for(i in a){if(a[i] ~ /[[:lower:]]/){str=a[i];str1=str1""str }}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TLDR",$11,$12,$13":"$14":"$15,$16,str1}' | grep -v "NA:NA:NA" > ${OUT_PATH}/intersection/dis/HG002_2_tgs_TLDR.dis.bed
echo "ck3"

# xTea
# cat ${DATA_PATH}/xtea/HG002_2_tgs_xtea_d5.bed | sort | uniq | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$2,$3":"$5,$10}'  > ${DATA_PATH}/xtea/HG002_2_tgs_xtea_d5.table.bed
# grep xTea-all ${OUT_PATH}/intersection/pdf/giab_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${OUT_PATH}/intersection/GOLD.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"GOLD",$9,$10,$11":"$12,1,"-"}' | bedtools intersect -a - -b ${DATA_PATH}/xtea/HG002_2_tgs_xtea_d5.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"xTea-all",$11,$12,$13,0,$20}' > ${OUT_PATH}/intersection/dis/HG002_2_tgs_xTea-all.dis.bed
grep xTea ${OUT_PATH}/intersection/pdf/giab_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${OUT_PATH}/intersection/GOLD.bed -wa -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"GOLD",$9,$10,0":"0"~"$11,1,$13}' | bedtools intersect -a - -b ${DATA_PATH}/xtea/HG002_2_tgs_xtea_d5.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"xTea",$11,$12,$13":"$16":"$17,0,$20}' > ${OUT_PATH}/intersection/dis/HG002_2_tgs_xTea.dis.bed
echo "check4"


# LOCATE
# grep LOCATE-all ${OUT_PATH}/intersection/pdf/giab_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${OUT_PATH}/intersection/GOLD.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"GOLD",$9,$10,$11":"$12,1,"-"}' | bedtools intersect -a - -b ${DATA_PATH}/LOCATE/HG002_2.txt -wa -wb | awk -v LOCATE_version=${LOCATE_version} 'BEGIN{FS=OFS="\t"}{split($18,a,",");min_st = 10000;max_end = 0;for(i=1;i<=length(a);i++ ){if(a[i]!~/Poly/){split(a[i],b,":");split(b[2],c,"-");if(c[1]<min_st){min_st=c[1]};if(c[2]>max_end){max_end=c[2]} }};gsub(/N*/,"",$25);print $1,$2,$3,$4,$5,$6,$7,$8,$9,LOCATE_version"-all",$11+1,$12+1,b[1]":"min_st":"max_end,0,$25}' > ${OUT_PATH}/intersection/dis/HG002_2_tgs_${LOCATE_version}-all.dis.bed
grep LOCATE ${OUT_PATH}/intersection/pdf/giab_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${OUT_PATH}/intersection/GOLD.bed -wa -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"GOLD",$9,$10,0":"0"~"$11,1,$13}'| bedtools intersect -a - -b ${DATA_PATH}/LOCATE/HG002_2.txt -wa -wb | awk -v LOCATE_version=${LOCATE_version} 'BEGIN{FS=OFS="\t"}{split($18,a,",");min_st = 10000;max_end = 0;for(i=1;i<=length(a);i++ ){if(a[i]!~/Poly/){split(a[i],b,":");split(b[2],c,"-");if(c[1]<min_st){min_st=c[1]};if(c[2]>max_end){max_end=c[2]} }};gsub(/N*/,"",$25);print $1,$2,$3,$4,$5,$6,$7,$8,$9,LOCATE_version,$11+1,$12+1,b[1]":"min_st":"max_end,0,$25}' > ${OUT_PATH}/intersection/dis/HG002_2_tgs_${LOCATE_version}.dis.bed



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
        dis_bed_list=${OUT_PATH}/intersection/dis/HG002_2_tgs.${meth}.bp_dis.bed
    else
        dis_bed_list=${dis_bed_list}" "${OUT_PATH}/intersection/dis/HG002_2_tgs.${meth}.bp_dis.bed
    fi

    # [ -f ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed ] && rm ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed

    if [ ! -f ${OUT_PATH}/intersection/dis/HG002_2_tgs.${meth}.bp_dis.bed ];then
        python "${script_dir}/02_cal_dis_div3_giab.py" ${OUT_PATH}/intersection/dis/HG002_2_tgs_${meth}.dis.bed /zata/zippy/boxu/annotation/hg38/ALSE.size ${OUT_PATH}/intersection/dis/HG002_2_tgs.${meth}.bp_dis.bed HG002_2_tgs ${meth} ${depth} ${OUT_PATH_dir} hg38 "${script_dir}/03_cal_div3_giab.sh"
    fi
done 
cat ${dis_bed_list} | LC_COLLATE=C sort -k1,1 -k2,2n - | bedtools merge -i - -c 4,5,6 -o first,first,collapse -delim ";" > ${OUT_PATH}/intersection/dis/giab.merge.bp_dis.bed



# "${script_dir}/01_intersect_gold.sh" -i ${temp2_ins},${melt_ins},${tldr_ins},${xtea_ins},${temp3_ins} -g ${gold_ins} -d ${dep} -t ${simulation_sample} -c ${OUT_PATH_dir}/pp_intersection/data -o ${OUT_PATH}/intersection/pdf/ -r ${user_path}/annotation/hg38/ALSE.bed -b /zata/zippy/boxu/annotation/from_temp3/BlackList.bed -p /zata/zippy/boxu/annotation/from_temp3/gap.bed

# prefix=bench_set
# # Rscript ${R_script_path}/intersection_precision_giab.R MELT,TEMP2,TLDR,xTea,TEMP3 ${OUT_PATH_dir}/pp_intersection ${prefix}
# Rscript ${R_script_path}/intersection_precision_all.R TEMP2,MELT,TLDR,xTea,TEMP3,MEHunter ${OUT_PATH_dir}/pp_intersection ${prefix}

# R --slave --no-restore --file=${R_script_path}/intersection_precision_one_depth.R --args MELT,TEMP2,SMS,TLDR,xTea,TEMP3 ${OUT_PATH_dir}/pp_intersection ${prefix}
# Rscript ${R_script_path}/intersection_upset.R ${OUT_PATH}/intersection/pdf/matrix.${dep}.txt ${OUT_PATH_dir}/pp_intersection/figure/${prefix}_upset.pdf MELT,TEMP2,TLDR,xTea,TEMP3,GOLD ${OUT_PATH_dir}/pp_intersection ${prefix}
# Rscript ${R_script_path}/intersection_venn.R ${OUT_PATH}/intersection/pdf/validation.intersect.txt ${OUT_PATH_dir}/pp_intersection/figure/venn_giab_hg38.pdf MELT,TEMP2,SMS,TLDR,xTea,TEMP3,GOLD # GOLD,SMS,xTea #  GOLD,TLDR,PALMER,SMS,xTea
