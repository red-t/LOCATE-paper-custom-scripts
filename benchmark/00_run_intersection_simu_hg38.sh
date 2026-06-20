
user_path="${user_path:-${USER_PATH:-/path/to/user}}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# while getopts ":v:" OPTION;
# do
#     case $OPTION in
#         v)  LOCATE_version=${OPTARG};;
#     esac
# done

LOCATE_version="LOCATE" #"LOCATE.ra.new" # "LOCATE.new.v3" # "LOCATE.rp"

INSERTION_PATH_dir="${user_path}/2022_long_reads/result/simulation/hg38/"
DATA_PATH_dir="${user_path}/2022_long_reads/result/simulation/hg38/insertion_revision"
OUT_PATH_dir="${user_path}/2022_long_reads/result/simulation/hg38/intersection_revision/" # "${user_path}/lrft2/result" # "${user_path}/lrft/result/simulation" 
ANNO_PATH="${user_path}/annotation"
GENOME="hg38"
SAMPLE="simulation"
CPU="15"



# intersection
# [ -d ${OUT_PATH_dir} ] && rm -r ${OUT_PATH_dir} 
[ ! -d ${OUT_PATH_dir} ] && mkdir ${OUT_PATH_dir}
[ ! -d ${OUT_PATH_dir}/pp_intersection ] && mkdir ${OUT_PATH_dir}/pp_intersection
[ ! -d ${OUT_PATH_dir}/pp_intersection/data ] && mkdir ${OUT_PATH_dir}/pp_intersection/data
[ ! -d ${OUT_PATH_dir}/pp_intersection/figure ] && mkdir ${OUT_PATH_dir}/pp_intersection/figure

[ -f ${OUT_PATH_dir}/pp_intersection/data/performance.txt ] && rm ${OUT_PATH_dir}/pp_intersection/data/performance*.txt
[ -f ${OUT_PATH_dir}/pp_intersection/data/false.bed ] && rm ${OUT_PATH_dir}/pp_intersection/data/false.bed
[ -f ${OUT_PATH_dir}/pp_intersection/data/miss.bed ] && rm ${OUT_PATH_dir}/pp_intersection/data/miss.bed


ex_ch="_alt|random|"


for simulation_sample in simulation_germ_ccs simulation_germ_clr simulation_germ_ont # simulation_germ_ont # simulation_germ_ccs simulation_germ_clr simulation_germ_ont # simulation_soma_ccs simulation_soma_ont simulation_soma_clr   
do
    for depth in 1 2 3 4 5 10 20 30 40 50
    do 

        OUT_PATH=${OUT_PATH_dir}/${depth}
        [ ! -d ${OUT_PATH} ] && mkdir ${OUT_PATH}
        [ ! -d ${OUT_PATH}/intersection ] && mkdir ${OUT_PATH}/intersection
        # [ -d ${OUT_PATH}/intersection/pdf ] && rm -r ${OUT_PATH}/intersection/pdf
        [ ! -d ${OUT_PATH}/intersection/pdf ] && mkdir ${OUT_PATH}/intersection/pdf
        [ ! -d ${OUT_PATH}/intersection/dis ] && mkdir ${OUT_PATH}/intersection/dis



        DATA_PATH=${DATA_PATH_dir}/${depth}
        [ ! -d ${DATA_PATH} ] && mkdir ${DATA_PATH}
        [ ! -d ${DATA_PATH}/MELT ] && mkdir ${DATA_PATH}/MELT
        [ ! -d ${DATA_PATH}/TEMP2 ] && mkdir ${DATA_PATH}/TEMP2
        [ ! -d ${DATA_PATH}/LOCATE ] && mkdir ${DATA_PATH}/LOCATE
        [ ! -d ${DATA_PATH}/xTea ] && mkdir ${DATA_PATH}/xTea
        [ ! -d ${DATA_PATH}/tldr ] && mkdir ${DATA_PATH}/tldr
        [ ! -d ${DATA_PATH}/PALMER ] && mkdir ${DATA_PATH}/PALMER
        [ ! -d ${DATA_PATH}/MEHunter ] && mkdir ${DATA_PATH}/MEHunter

        [ ! -d ${DATA_PATH}/TELR ] && mkdir ${DATA_PATH}/TELR
        [ ! -d ${DATA_PATH}/TrEMOLO ] && mkdir ${DATA_PATH}/TrEMOLO

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
            # ns_num=0
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
        GOLD_INS_SUMM="${user_path}/2022_long_reads/result/simulation/hg38/gold/ins.${insertion_type}.gold.summary"

        awk -v l=10 -v p=GOLD 'BEGIN{FS=OFS="\t"}{split($4,t,"~");if($5~/\[/){split($5,a,"[");split(a[2],b,".")}else{b[1]=1};if(b[1]>=0){print $1, $8-l, $9+l, t[2], p":gold:5:"$7":"$11, $6}}' ${GOLD_INS_SUMM}  | grep -v ${ex_ch} | bedtools intersect -a - -b /zata/zippy/boxu/annotation/from_temp3/GRCh38.gap.bed -v > ${OUT_PATH}/intersection/${simulation_sample}_GOLD.bed
        # awk -v l=10 -v p=GOLD 'BEGIN{FS=OFS="\t"}{split($4,a,"~");if($5~/\[/){split($5,c,"[");split(c[2],b,".")}else{b[1]=1};if(b[1]<=1){print $1, $8-l, $9+l, a[2], p":gold:5:"$7":"$11, $6}}' ${GOLD_INS_SUMM}  | grep -v ${ex_ch} | bedtools intersect -a - -b /zata/zippy/boxu/annotation/from_temp3/GRCh38.gap.bed -v > ${OUT_PATH}/intersection/${simulation_sample}_GOLD.bed
        
        gold_ins=${OUT_PATH}/intersection/${simulation_sample}_GOLD.bed

        ## MELT
        RESULT_PATH="${INSERTION_PATH}/melt/${simulation_sample}"
        [ -f ${RESULT_PATH}/melt.vcf.bed ] && rm ${RESULT_PATH}/melt.vcf.bed
        
        for melt_te in ALU LINE1 SVA HERVK;do
            if [ $((`cat ${RESULT_PATH}/${melt_te}.final_comp.vcf | wc -l`)) -gt '0' ];then
                # add frequency / genotype
                bcftools query -f "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/LP\t%INFO/RP\t%INFO/SVLEN\t%INFO/TSD\t%INFO/MEINFO\t%INFO/SVLEN\t[%GT]\n" ${RESULT_PATH}/${melt_te}.final_comp.vcf | awk -v l=10 'BEGIN{FS=OFS="\t"}{split($8,a,","); if($10=="1/1"){genotype="1"}else{genotype="0.5"}; print $1,$2-l,$2+l,$3,"melt_2",a[4],$4+$5,$9,a[2],a[3],"-",genotype}' >> ${RESULT_PATH}/melt.vcf.bed
            fi
        done
        cat ${RESULT_PATH}/*master.bed > ${RESULT_PATH}/melt.master1.bed

        if [ ! -f ${RESULT_PATH}/melt.vcf.bed ];then
            cat /dev/null > ${RESULT_PATH}/melt.vcf.bed
            cp ${RESULT_PATH}/melt.vcf.bed ${DATA_PATH}/MELT/${simulation_sample}.bed
        else
            cat ${RESULT_PATH}/melt.vcf.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools merge -i - -o first -c 4,5,6,7,8,9,10,11,12 > ${DATA_PATH}/MELT/${simulation_sample}.bed
        fi

        cat ${DATA_PATH}/MELT/${simulation_sample}.bed | awk -v l=10 -v p=MELT 'BEGIN{FS=OFS="\t"}{print $1,$2-l,$3+l,$4,p":germ:"$7":"$12":-",$6}' | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_MELT.bed
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

        # tldr
        sed 's/ERVK/HERVK/g' ${INSERTION_PATH}/tldr/${simulation_sample}_${depth}.table.txt > ${DATA_PATH}/tldr/${simulation_sample}.table.txt
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


        #### xtea
        cat ${INSERTION_PATH}/xTea/${simulation_sample}_${depth}X/classified_results.txt* | uniq | sed 's/HERV-K/HERVK/g' > ${DATA_PATH}/xTea/${simulation_sample}.txt
        cat ${DATA_PATH}/xTea/${simulation_sample}.txt | sort | uniq | awk 'BEGIN{FS=OFS="\t"}{split($5,a,":");split($2,b," ");print $1,b[1],b[1]+1,$3,"xTea:"$3":5:1",a[3],a[2]-a[1], $10}' | LC_COLLATE=C sort -k1,1 -k2,2n  | awk -v l=5 -v p=xTea 'BEGIN{FS=OFS="\t"}{print $1,$2-l,$3+l,$4,p":germ:5:1:"$8,$6}' | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_xTea.bed
        xtea_ins=${OUT_PATH}/intersection/${simulation_sample}_xTea.bed

        # xtea-pass
        [ -f ${DATA_PATH}/xTea/${simulation_sample}.pass.txt ] && rm ${DATA_PATH}/xTea/${simulation_sample}.pass.txt
        for TE in LINE1 SVA ALU HERV;do 
            cat ${INSERTION_PATH}/xTea/${simulation_sample}_${depth}X/classified_results.txt.${TE}.txt >> ${DATA_PATH}/xTea/${simulation_sample}.pass.txt
        done
        # cat ${INSERTION_PATH}/xTea/${simulation_sample}_${depth}X/classified_results.txt.ALU.txt ${INSERTION_PATH}/xTea/${simulation_sample}_${depth}X/classified_results.txt.LINE1.txt    > ${DATA_PATH}/xTea/${simulation_sample}.txt
        cat ${DATA_PATH}/xTea/${simulation_sample}.pass.txt | sort | uniq | sed 's/HERV-K/HERVK/g' | awk 'BEGIN{FS=OFS="\t"}{split($5,a,":");split($2,b," ");print $1,b[1],b[1]+1,$3,"xTea-pass:"$3":5:1",a[3],a[2]-a[1], $10}' | LC_COLLATE=C sort -k1,1 -k2,2n  | awk -v l=5 -v p=xTea-pass 'BEGIN{FS=OFS="\t"}{print $1,$2-l,$3+l,$4,p":germ:5:1:"$8,$6}'| grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_xTea-pass.bed
        xtea_pass_ins=${OUT_PATH}/intersection/${simulation_sample}_xTea-pass.bed

        [ -f ${DATA_PATH}/xTea/${simulation_sample}.all.txt ] && rm ${DATA_PATH}/xTea/${simulation_sample}.all.txt
        cat ${INSERTION_PATH}/xTea/${simulation_sample}_${depth}X/classified_results.txt* | sort | uniq | sed 's/HERV-K/HERVK/g' > ${DATA_PATH}/xTea/${simulation_sample}.all.txt
        cat ${DATA_PATH}/xTea/${simulation_sample}.all.txt | awk 'BEGIN{FS=OFS="\t"}{split($5,a,":");split($2,b," ");print $1,b[1],b[1]+1,$3,"xTea-all:"$3":5:1",a[3],a[2]-a[1], $10}' | LC_COLLATE=C sort -k1,1 -k2,2n | awk -v l=5 -v p=xTea-all 'BEGIN{FS=OFS="\t"}{print $1,$2-l,$3+l,$4,p":germ:5:1:"$8,$6}' | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_xTea-all.bed
        xtea_all_ins=${OUT_PATH}/intersection/${simulation_sample}_xTea-all.bed



        #### TELR
        
        vcf_result=` cat ${INSERTION_PATH}/TELR/${simulation_sample}/${simulation_sample}.${depth}X.telr.vcf | wc -l `
        if [ $vcf_result -gt 0 ];then
            sed 's/ID=[0-9]\+;//g' ${INSERTION_PATH}/TELR/${simulation_sample}/${simulation_sample}.${depth}X.telr.vcf | bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/FAMILY\t%INFO/AF\t%INFO/STRANDS\t%INFO/RE\t%INFO/TSD_SEQ\t%ALT\n" - | awk 'BEGIN{OFS="\t"}{split($4,a,"|");if($2>$3){st=$2;en=$2}else{st=$2;en=$3};print $1,st,en,a[1],$5,$6,$7,$8,$9}'> ${DATA_PATH}/TELR/${simulation_sample}.telr.bed
            awk -v l=10 -v p=TELR -v ns=${ns_num} 'BEGIN{OFS="\t"}{split($4,a,"|");if($2>$3){st=$3;en=$2}else{st=$2;en=$3}; if($7>ns){print $1,st,en+1,a[1],p":INS:"$7":"$5":"$9,$6}}' ${DATA_PATH}/TELR/${simulation_sample}.telr.bed > ${OUT_PATH}/intersection/${simulation_sample}_TELR.bed
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


        ### MEHunter
        # Extract AF from cuteSV results.
        bcftools query -f "%CHROM\t%POS\t%INFO/AF\t[%GT]\t%ID\n" ${INSERTION_PATH}/MEHunter_v2/${simulation_sample}/MEHunter.vcf > ${INSERTION_PATH}/MEHunter_v2/${simulation_sample}/MEHunter.fre.bed
        sed "s/'/ /g" ${INSERTION_PATH}/MEHunter_v2/${simulation_sample}/High_Quality_MEIs.txt | awk 'BEGIN{OFS="\t"}NR==FNR{ split($5,ins_id,"_"); a[ins_id[1]]=$3}NR>FNR{print $4,$6,$6,$12,$10,".",$8,a[$2]}' ${INSERTION_PATH}/MEHunter_v2/${simulation_sample}/MEHunter.fre.bed - > ${DATA_PATH}/MEHunter/${simulation_sample}.MEHunter.bed
        awk -v l=10 -v p=MEHunter  'BEGIN{OFS="\t"}{print $1, $2-l, $3+l, $4, p":INS:"$5":"$8":"$7, $6}' ${DATA_PATH}/MEHunter/${simulation_sample}.MEHunter.bed > ${OUT_PATH}/intersection/${simulation_sample}_MEHunter.bed
        MEHunter_ins=${OUT_PATH}/intersection/${simulation_sample}_MEHunter.bed


        #### PALMER
        # AF is not available.
        [ ! -d ${INSERTION_PATH}/PALMER/${simulation_sample} ] && mkdir -p ${INSERTION_PATH}/PALMER/${simulation_sample}
        cp ${user_path}/2022_long_reads/result/simulation/hg38/${depth}/PALMER/${simulation_sample}/* ${INSERTION_PATH}/PALMER/${simulation_sample}
        [ -f ${INSERTION_PATH}/PALMER/${simulation_sample}/${simulation_sample}.${depth}X.palmer.bed ] && rm ${INSERTION_PATH}/PALMER/${simulation_sample}/${simulation_sample}.${depth}X.palmer.bed
        for palmerTE in ALU LINE SVA HERVK
        do
            # cat ${INSERTION_PATH}/PALMER/${simulation_sample}/${palmerTE}/${simulation_sample}.50X.${palmerTE}_calls.txt >> ${INSERTION_PATH}/PALMER/${simulation_sample}/${simulation_sample}.50X.palmer.bed
            if [ -f ${INSERTION_PATH}/PALMER/${simulation_sample}/${simulation_sample}.${depth}X.${palmerTE}_calls.txt ];then
                sed 1d ${INSERTION_PATH}/PALMER/${simulation_sample}/${simulation_sample}.${depth}X.${palmerTE}_calls.txt | awk -v te=${palmerTE} -v dep=${depth} 'BEGIN{OFS="\t"}{if($11>=1 && $12 >= 0.1*dep){print $2,$3,$4,te,$11,$14,$7,$9,$15,$16,$17}}'   >> ${INSERTION_PATH}/PALMER/${simulation_sample}/${simulation_sample}.${depth}X.palmer.bed
            fi
        done
        # awk -v te=${palmerTE} 'BEGIN{OFS="\t"}{print $2,$3,$4,te,$11,$14,$7,$9,$15,$16,$17}' ${INSERTION_PATH}/PALMER/${simulation_sample}/${simulation_sample}.50X.palmer.bed | sed 1d >> ${DATA_PATH}/PALMER/${simulation_sample}.palmer.bed # PALMER/${simulation_sample}/${simulation_sample}.50X.palmer.bed
        # cp ${INSERTION_PATH}/PALMER/${simulation_sample}/${simulation_sample}.50X.palmer.bed ${DATA_PATH}/PALMER/${simulation_sample}.palmer.bed
        LC_COLLATE=C sort -k1,1 -k2,2n ${INSERTION_PATH}/PALMER/${simulation_sample}/${simulation_sample}.${depth}X.palmer.bed | bedtools merge -i - -c 4,5,6,7,8,9,10,11 -d 5 -o first  > ${DATA_PATH}/PALMER/${simulation_sample}.palmer.bed
        
        cat ${DATA_PATH}/PALMER/${simulation_sample}.palmer.bed | awk -v l=10 -v p=PALMER -v ns=${ns_num} 'BEGIN{FS=OFS="\t"}{if($5>ns){print $1,$2-l,$3+l,$4,p":germ:"$5":1:-",$6}}' > ${OUT_PATH}/intersection/${simulation_sample}_PALMER.bed
        palmer_ins=${OUT_PATH}/intersection/${simulation_sample}_PALMER.bed


        ### GraffiTE
        
        # cp -r ${user_path}/2022_long_reads/result/simulation/hg38/${depth}/GraffiTE/${simulation_sample}/ ${INSERTION_PATH}/GraffiTE/${simulation_sample}
        # genotype = false 
        # bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/repeat_ids\t%INFO/SUPPORT\t%INFO/RM_hit_strands\t%ID\t%ALT\t[%GT]\n" ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/genotypes_repmasked_filtered.vcf > ${INSERTION_PATH}/GraffiTE/${simulation_sample}/GraffiTE.tmp.bed
        
        INSERTION_PATH_GT="${user_path}/2022_long_reads/result/simulation/hg38/test_graffiTE/${depth}"
        # genotype = true
        bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/repeat_ids\t%INFO/SUPPORT\t%INFO/RM_hit_strands\t%ID\t%ALT\t[%GT]\n" ${INSERTION_PATH_GT}/GraffiTE/${simulation_sample}/4_Genotyping/GraffiTE.merged.genotypes.vcf.gz > ${INSERTION_PATH_GT}/GraffiTE/${simulation_sample}/GraffiTE.tmp.bed

        # awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$9"_"$11"_"length($12)"_"$12}NR>FNR{if($7 in a){split(a[$7],b,"_"); if(b[1]=="0"){st=$2-b[3];en=$2}else{st=$2;en=$2+b[3]} }else{st=$2;en=$2+1}; split($4,c,",");split($6,d,",");if(d[1]=="C"){strand="-"}else{strand="+"}; print $1,st,en,c[1],$5,strand,$7,$8,a[$7] }' ${INSERTION_PATH}/GraffiTE/${simulation_sample}/3_TSD_search/TSD_summary.txt  ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/GraffiTE.tmp.bed > ${INSERTION_PATH}/GraffiTE/${simulation_sample}/2_Repeat_Filtering/${simulation_sample}_GraffiTE.bed
        awk 'BEGIN{OFS="\t"}NR==FNR{if($14=="PASS"){a[$1]=$9"_"$11"_"length($12)"_"$12}}
                            NR>FNR{if($7 in a){split(a[$7],b,"_"); st=$2;en=$2+b[3];TSD=b[4];split($8,c,TSD);if(c[1]!=$8){ins_st=length(c[1])+b[3]+1;ins_len=length($8)-ins_st+1; ins_seq=substr($8,ins_st,ins_len)}else{ins_seq=substr($8,1,length(c[1]))} }
                                   else{a[$7]=".";st=$2;en=$2+1;ins_seq=$8}; split($4,c,",");split($6,d,","); 
                                        if(d[1]=="C"){strand="-"}else{strand="+"}; 
                                   if($9=="0/1"){AF=0.5}else{ if($9=="1/1"){AF=1}else{AF=0} };
                                   print $1,st,en,c[1],$5,strand,$7,ins_seq,a[$7],AF }' ${INSERTION_PATH_GT}/GraffiTE/${simulation_sample}/3_TSD_search/TSD_summary.txt  ${INSERTION_PATH_GT}/GraffiTE/${simulation_sample}/GraffiTE.tmp.bed | grep -v "DEL" > ${INSERTION_PATH_GT}/GraffiTE/${simulation_sample}/4_Genotyping/${simulation_sample}_GraffiTE.bed
        
        cp ${INSERTION_PATH_GT}/GraffiTE/${simulation_sample}/4_Genotyping/${simulation_sample}_GraffiTE.bed ${DATA_PATH}/GraffiTE/${simulation_sample}.GraffiTE.bed
        awk -v l=10 -v p=GraffiTE -v ns=${ns_num} 'BEGIN{OFS="\t"}{if($5>ns){print $1,$2,$3,$4,p":INS:"$5":"$10":"$8,$6}}' ${DATA_PATH}/GraffiTE/${simulation_sample}.GraffiTE.bed > ${OUT_PATH}/intersection/${simulation_sample}_GraffiTE.bed
        GraffiTE_ins=${OUT_PATH}/intersection/${simulation_sample}_GraffiTE.bed


        ##### LOCATE
        # cat ${INSERTION_PATH}/${LOCATE_version}/${simulation_sample}/${simulation_sample}.txt | sed 's/ERVK/HERVK/g' > ${DATA_PATH}/LOCATE/${simulation_sample}.txt

        # New version with genotype support.
        # cp ${user_path}/2022_long_reads/result/simulation/hg38/${depth}/${LOCATE_version}/${simulation_sample}/${simulation_sample}.txt ${DATA_PATH}/LOCATE/${simulation_sample}.txt
        sed 1d ${user_path}/2022_long_reads/result/simulation/hyplotype_LOCATE/${depth}/${LOCATE_version}/${simulation_sample}/result.tsv | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14,$15,$16,$17}' >  ${DATA_PATH}/LOCATE/${simulation_sample}.txt


        if [[ $simulation_sample =~ "soma" ]];then
            # awk '{if(and($19,1)){print $0}}' ${DATA_PATH}/LOCATE/${simulation_sample}.txt | awk -v ns=${ns_num}  'BEGIN{FS=OFS="\t"}{if($10==1){print $1,$2,$3,$4,"LOCATE:germ:"$10":"$5":"$16, $6}}'  | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_LOCATE.bed

            awk '{if(and($19,1)){print $0}}' ${DATA_PATH}/LOCATE/${simulation_sample}.txt | awk -v ns=${ns_num}  'BEGIN{FS=OFS="\t"}{if($12==1){gsub(/N*/,"",$13); print $1,$2,$3,$4,"LOCATE-pass:germ:"$10":"$5":"$16, $6}}'  | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_LOCATE-pass.bed
            cat ${DATA_PATH}/LOCATE/${simulation_sample}.txt | awk -v ns=${ns_num}  'BEGIN{FS=OFS="\t"}{if($10==1){gsub(/N*/,"",$13); print $1,$2,$3,$4,"LOCATE-all:germ:"$10":"$5":"$16, $6}}'  | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_LOCATE-all.bed

        else
            awk '{if( $8=="True" ){print $0}}' ${DATA_PATH}/LOCATE/${simulation_sample}.txt | awk -v ns=${ns_num}  'BEGIN{FS=OFS="\t"}{if($11>ns){gsub(/N*/,"",$13); if($7=="1/1"){geno=1}else{if($7=="0/1"){geno=0.5}else{geno=0}}; print $1,$2+1,$3+1,$4,"LOCATE-pass:germ:"$11":"geno":"$13, $6}}'  | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_LOCATE-pass.bed
            cat ${DATA_PATH}/LOCATE/${simulation_sample}.txt | awk -v ns=${ns_num}  'BEGIN{FS=OFS="\t"}{if($11>ns){ gsub(/N*/,"",$13);if($7=="1/1"){geno=1}else{if($7=="0/1"){geno=0.5}else{geno=0}};   print $1,$2+1,$3+1,$4,"LOCATE-all:germ:"$11":"geno":"$13, $6}}'  | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_LOCATE-all.bed
        fi
        
        # awk '{if(and($19,1)){print $0}}' ${DATA_PATH}/TEMP3/${simulation_sample}.txt | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,"TEMP3:germ:"$10":"$5":"$16, $6}'  | grep -v ${ex_ch} > ${OUT_PATH}/intersection/${simulation_sample}_TEMP3.bed
        temp3_all_ins=${OUT_PATH}/intersection/${simulation_sample}_LOCATE-all.bed
        temp3_pass_ins=${OUT_PATH}/intersection/${simulation_sample}_LOCATE-pass.bed





        echo ">>>>intersection"
        # ${user_path}/2022_long_reads/bin_inter/intersect_gold.sh -i ${temp2_ins},${melt_ins},${tldr_pass_ins},${tldr_all_ins},${xtea_pass_ins},${xtea_all_ins},${telr_ins},${TrEMOLO_ins},${palmer_ins},${MEHunter_ins},${temp3_all_ins},${temp3_pass_ins} -g ${gold_ins} -d ${depth} -t ${simulation_sample} -c ${OUT_PATH_dir}/pp_intersection/data -o ${OUT_PATH}/intersection/pdf/ -r ${user_path}/annotation/hg38/ALSE.bed -b /zata/zippy/boxu/annotation/from_temp3/BlackList.bed -p /zata/zippy/boxu/annotation/from_temp3/gap.bed
        "${script_dir}/01_intersect_gold.sh" -i ${temp2_ins},${melt_ins},${tldr_pass_ins},${tldr_all_ins},${xtea_pass_ins},${xtea_all_ins},${telr_ins},${TrEMOLO_ins},${palmer_ins},${MEHunter_ins},${temp3_all_ins},${temp3_pass_ins},${GraffiTE_ins} -g ${gold_ins} -d ${depth} -t ${simulation_sample} -c ${OUT_PATH_dir}/pp_intersection/data -o ${OUT_PATH}/intersection/pdf/ -r ${user_path}/annotation/hg38/RepeatMasker/new.cns.with_polyA/repeat.cons.with_polyA.bed -b /zata/zippy/boxu/annotation/from_temp3/BlackList.bed -p /zata/zippy/boxu/annotation/from_temp3/gap.bed

        echo ">>>>distance"
        ## distance type frequence
        ## melt
        if [ $((`cat ${DATA_PATH}/MELT/${simulation_sample}.bed | wc -l`)) -gt '0' ];then

            cat ${DATA_PATH}/MELT/${simulation_sample}.bed | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4":"$9":"$10,$12}' > ${DATA_PATH}/MELT/${simulation_sample}.melt.bed
            grep MELT ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"GOLD",$14,$15,$10":"$11,$13,$17}' | bedtools intersect -a - -b ${DATA_PATH}/MELT/${simulation_sample}.melt.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,"MELT",$11-10,$11-10,$13,$14,"-"}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_MELT.dis.bed
        else
            cat /dev/null > ${OUT_PATH}/intersection/dis/${simulation_sample}_MELT.dis.bed
        fi 
        echo "ck1"
        
        # TEMP2
        grep TEMP2 ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"GOLD",$14,$15,$10":"$11,$13,$17}' | bedtools intersect -a - -b ${DATA_PATH}/TEMP2/${simulation_sample}.insertion.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{split($13,a,":");print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TEMP2",$11+1,$12,a[1]":"a[2]-1":"a[3],$14,"-"}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_TEMP2.dis.bed
        echo "ck2"


        # tldr
        cat ${DATA_PATH}/tldr/${simulation_sample}.table.txt | awk 'BEGIN{FS=OFS="\t"}{print $2,$3+1,$4+1,$7,$8,$9,$12,$22}' | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g" > ${DATA_PATH}/tldr/${simulation_sample}.table.bed
        grep TLDR-all ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"GOLD",$14,$15,$10":"$11,$13,$17}' | bedtools intersect -a - -b ${DATA_PATH}/tldr/${simulation_sample}.table.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{str1="";split($17,a,""); for(i in a){if(a[i] ~ /[[:lower:]]/){str=a[i];str1=str1""str }}; if(str1==""){str1="-"}; freq=$16;if($16=="NA"){freq="-"}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TLDR-all",$11,$12,$13":"$14":"$15,freq,str1}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_TLDR-all.dis.bed
        grep TLDR-pass ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"GOLD",$14,$15,$10":"$11,$13,$17}' | bedtools intersect -a - -b ${DATA_PATH}/tldr/${simulation_sample}.table.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{str1="";split($17,a,""); for(i in a){if(a[i] ~ /[[:lower:]]/){str=a[i];str1=str1""str }}; freq=$16;if($16=="NA"){freq="-"}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TLDR-pass",$11,$12,$13":"$14":"$15,$16,str1}' | grep -v "NA:NA:NA" > ${OUT_PATH}/intersection/dis/${simulation_sample}_TLDR-pass.dis.bed
        echo "ck3"

        # xTea
        cat ${DATA_PATH}/xTea/${simulation_sample}.txt | sort | uniq | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$2,$3":"$5,$10}'  > ${DATA_PATH}/xTea/${simulation_sample}.bed
        grep xTea-all ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"GOLD",$14,$15,$10":"$11,$13,$17}' | bedtools intersect -a - -b ${DATA_PATH}/xTea/${simulation_sample}.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"xTea-all",$11,$12,$13,"-",$14}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_xTea-all.dis.bed
        grep xTea-pass ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"GOLD",$14,$15,$10":"$11,$13,$17}' | bedtools intersect -a - -b ${DATA_PATH}/xTea/${simulation_sample}.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"xTea-pass",$11,$12,$13,"-",$14}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_xTea-pass.dis.bed

        echo "check4"



        # TELR
        grep TELR ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{FS=OFS="\t"}{if($8>$14){g_st=$14}else{g_st=$8};if($9>$15){g_en=$9}else{g_en=$15}; print $1,$2,$3,"GOLD",g_st,g_en,$10":"$11,$13,$17}' | bedtools intersect -a - -b ${DATA_PATH}/TELR/${simulation_sample}.telr.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{if($2>0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TELR",$11,$12,$13":0:0",$14,$18}}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_TELR.dis.bed
        echo "ck4"
        
        # TrEMOLO
        grep TrEMOLO ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{FS=OFS="\t"}{if($8>$14){g_st=$14}else{g_st=$8};if($9>$15){g_en=$9}else{g_en=$15}; print $1,$2,$3,"GOLD",g_st,g_en,$10":"$11,$13,$17}' | bedtools intersect -a - -b ${DATA_PATH}/TrEMOLO/${simulation_sample}.TrEMOLO.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{if($2>0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TrEMOLO",$11,$12,$13":0:0",$21,$20}}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_TrEMOLO.dis.bed
        echo "ck5"

        ## MEHunter
        grep MEHunter ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{FS=OFS="\t"}{if($8>$14){g_st=$14}else{g_st=$8};if($9>$15){g_en=$9}else{g_en=$15}; print $1,$2,$3,"GOLD",g_st,g_en,$10":"$11,$13,$17}' | bedtools intersect -a - -b ${DATA_PATH}/MEHunter/${simulation_sample}.MEHunter.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{if($2>0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,"MEHunter",$11,$12,$13":0:0",$17,$16}}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_MEHunter.dis.bed
        echo "ck6"

        ### PALMER
        grep PALMER ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb |  awk 'BEGIN{FS=OFS="\t"}{if($8>$14){g_st=$14}else{g_st=$8};if($9>$15){g_en=$9}else{g_en=$15}; print $1,$2,$3,"GOLD",g_st,g_en,$10":"$11,$13,$17}' | bedtools intersect -a - -b ${DATA_PATH}/PALMER/${simulation_sample}.palmer.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"PALMER",$11,$12,$13":"$16":"$17,"-","-"}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_PALMER.dis.bed
        echo "ck7"
        
        ### GraffiTE
        grep GraffiTE ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{FS=OFS="\t"}{if($8>$14){g_st=$14}else{g_st=$8};if($9>$15){g_en=$9}else{g_en=$15}; print $1,$2,$3,"GOLD",g_st,g_en,$10":"$11,$13,$17}' | bedtools intersect -a - -b ${DATA_PATH}/GraffiTE/${simulation_sample}.GraffiTE.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"}{ seq=$17;if($17==""){seq="-"}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,"GraffiTE",$11,$12,$13":"0":"0,$19,seq}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_GraffiTE.dis.bed
        echo "ck8"
        
        # LOCATE
        grep LOCATE-all ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{FS=OFS="\t"}{if($8>$14){g_st=$14}else{g_st=$8};if($9>$15){g_en=$9}else{g_en=$15}; print $1,$2,$3,"GOLD",g_st,g_en,$10":"$11,$13,$17}' | bedtools intersect -a - -b ${DATA_PATH}/LOCATE/${simulation_sample}.txt -wa -wb | awk -v LOCATE_version=${LOCATE_version} 'BEGIN{FS=OFS="\t"}{split($19,a,",");min_st = 10000;max_end = 0;for(i=1;i<=length(a);i++ ){if(a[i]!~/Poly/){split(a[i],b,":");split(b[2],c,"-");if(c[1]<min_st){min_st=c[1]};if(c[2]>max_end){max_end=c[2]} }};gsub(/N*/,"",$22); if($16=="1/1"){geno=1}else{if($16=="0/1"){geno=0.5}else{geno=0}};  print $1,$2,$3,$4,$5,$6,$7,$8,$9,LOCATE_version"-all",$11+1,$12+1,b[1]":"min_st":"max_end,geno,$22}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_${LOCATE_version}-all.dis.bed
        grep LOCATE-pass ${OUT_PATH}/intersection/pdf/${simulation_sample}_validation.germline.bed | grep GOLD | bedtools intersect -a - -b ${GOLD_INS_SUMM} -wa -wb | awk 'BEGIN{FS=OFS="\t"}{if($8>$14){g_st=$14}else{g_st=$8};if($9>$15){g_en=$9}else{g_en=$15}; print $1,$2,$3,"GOLD",g_st,g_en,$10":"$11,$13,$17}'| bedtools intersect -a - -b ${DATA_PATH}/LOCATE/${simulation_sample}.txt -wa -wb | awk -v LOCATE_version=${LOCATE_version} 'BEGIN{FS=OFS="\t"}{split($19,a,",");min_st = 10000;max_end = 0;for(i=1;i<=length(a);i++ ){if(a[i]!~/Poly/){split(a[i],b,":");split(b[2],c,"-");if(c[1]<min_st){min_st=c[1]};if(c[2]>max_end){max_end=c[2]} }};gsub(/N*/,"",$22); if($16=="1/1"){geno=1}else{if($16=="0/1"){geno=0.5}else{geno=0}};  print $1,$2,$3,$4,$5,$6,$7,$8,$9,LOCATE_version"-pass",$11+1,$12+1,b[1]":"min_st":"max_end,geno,$22}' > ${OUT_PATH}/intersection/dis/${simulation_sample}_${LOCATE_version}-pass.dis.bed
        echo "ck9"


        prefix=${simulation_sample}
        dis_bed_list=""
        # for meth in ${LOCATE_version}-pass-TLDR-pass TLDR-pass-${LOCATE_version}-pass # TLDR-pass-spe ${LOCATE_version}-pass-TLDR-pass ${LOCATE_version}-pass-spe ${LOCATE_version}-pass TLDR-pass ${LOCATE_version}-all MELT TEMP2 TLDR-all xTea-all xTea-pass
        
        # rm ${OUT_PATH}/intersection/dis/${simulation_sample}.${LOCATE_version}-all.bp_dis.bed
        # rm ${OUT_PATH}/intersection/dis/${simulation_sample}.${LOCATE_version}-pass.bp_dis.bed
        # rm ${OUT_PATH}/intersection/dis/${simulation_sample}.GraffiTE.bp_dis.bed
        rm ${OUT_PATH}/intersection/dis/${simulation_sample}.*.bp_dis.bed


        for meth in ${LOCATE_version}-all ${LOCATE_version}-pass TLDR-pass TLDR-all xTea-all xTea-pass GraffiTE PALMER MEHunter TrEMOLO TELR MELT TEMP2  # MELT TEMP2 xTea-all xTea-pass PALMER TELR TrEMOLO MEHunter
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
                python "${script_dir}/02_cal_dis_div3.py" ${OUT_PATH}/intersection/dis/${simulation_sample}_${meth}.dis.bed /zata/zippy/boxu/annotation/hg38/ALSE.size ${OUT_PATH}/intersection/dis/${simulation_sample}.${meth}.bp_dis.bed ${simulation_sample} ${meth} ${depth} ${OUT_PATH_dir} hg38 "${script_dir}/03_cal_div3.sh"
            fi
        done 
        
        cat ${dis_bed_list} | LC_COLLATE=C sort -k1,1 -k2,2n - | bedtools merge -i - -c 4,5,6 -o first,first,collapse -delim ";" > ${OUT_PATH}/intersection/dis/${simulation_sample}.merge.bp_dis.bed


    done
done

# Rscript ${R_script_path}/intersection_precision_all.R TEMP2,MELT,TLDR,xTea,TEMP3 ${OUT_PATH_dir}/pp_intersection human

# Rscript ${R_script_path}/intersection_TP_FP_count_all.R TEMP2,MELT,TLDR,xTea,TEMP3 ${OUT_PATH_dir}/pp_intersection human
# Rscript ${R_script_path}/performance_dis.R ${OUT_PATH}/intersection/dis/bp_dis.txt ${OUT_PATH_dir}/pp_intersection human

# python ${user_path}/2022_long_reads/bin/cal_dis_div2.py ${OUT_PATH}/intersection/dis/${simulation_sample}_${meth}.dis.bed /zata/zippy/boxu/annotation/hg38/ALSE.size ${OUT_PATH}/intersection/dis/${simulation_sample}.bp_dis.txt ${simulation_sample} ${meth} ${depth}

# python ${user_path}/2022_long_reads/bin/cal_dis_div2.py ./simulation_germ_clr_TLDR-pass.dis.bed  /zata/zippy/boxu/annotation/hg38/ALSE.size t.txt simulation_germ_clr TLDR-pass 5

sed -i "s/-/_/g" ${OUT_PATH_dir}/pp_intersection/data/performance*.txt

