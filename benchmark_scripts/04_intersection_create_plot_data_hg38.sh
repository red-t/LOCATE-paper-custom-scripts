#!/bin/bash

user_path="${user_path:-${USER_PATH:-/path/to/user}}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOME="hg38"
intersection_data_path="${user_path}/2022_long_reads/result/simulation/${GENOME}/intersection_revision"
echo -e "sample\tdepth\tmeth\tcaller_count\tgold_count" > ${intersection_data_path}/all.merge.count.txt

for SAMPLE in simulation_germ_ccs simulation_germ_clr simulation_germ_ont # simulation_germ_clr simulation_germ_ccs # simulation_soma_ont simulation_soma_clr simulation_soma_ccs #   # simulation_soma # HG001 # HG002 # SRR9685183.sra
do
    for depth in 1 2 3 4 5 10 20 30 40 50
    do 
        gold_count=` cat ${intersection_data_path}/${depth}/intersection/${SAMPLE}_GOLD.bed | wc -l `

        for meth in LOCATE-all LOCATE-pass TLDR-pass TLDR-all xTea-all xTea-pass GraffiTE PALMER MEHunter TrEMOLO TELR MELT TEMP2  # MELT TEMP2 xTea-all xTea-pass PALMER TELR TrEMOLO MEHunter
        do
            caller_count=` cat ${intersection_data_path}/${depth}/intersection/${SAMPLE}_${meth}.bed | wc -l `
            
            echo -e "${SAMPLE}\t${depth}\t${meth}\t${caller_count}\t${gold_count}" >> ${intersection_data_path}/all.merge.count.txt

        done

    done
done



GENOME="hg38"

intersection_data_path="${user_path}/2022_long_reads/result/simulation/${GENOME}/intersection_revision"

# MEHunter.bp_dis\tMEHunter.te_dis\tMEHunter.te_identity\tMEHunter.te_len\tMEHunter.div
echo -e "insID\tsample\tdepth\tTE\ttype\tfrequency\tLOCATE_all.bp_dis\tLOCATE_all.te_dis\tLOCATE_all.te_identity\tLOCATE_all.te_len\tLOCATE_all.div\tLOCATE_all.genotype\tLOCATE_pass.bp_dis\tLOCATE_pass.te_dis\tLOCATE_pass.te_identity\tLOCATE_pass.te_len\tLOCATE_pass.div\tLOCATE_pass.genotype\tMELT.bp_dis\tMELT.te_dis\tMELT.te_identity\tMELT.te_len\tMELT.div\tMELT.genotype\tTEMP2.bp_dis\tTEMP2.te_dis\tTEMP2.te_identity\tTEMP2.te_len\tTEMP2.div\tTEMP2.genotype\tTLDR_all.bp_dis\tTLDR_all.te_dis\tTLDR_all.te_identity\tTLDR_all.te_len\tTLDR_all.div\tTLDR_all.genotype\tTLDR_pass.bp_dis\tTLDR_pass.te_dis\tTLDR_pass.te_identity\tTLDR_pass.te_len\tTLDR_pass.div\tTLDR_pass.genotype\txTea_all.bp_dis\txTea_all.te_dis\txTea_all.te_identity\txTea_all.te_len\txTea_all.div\txTea_all.genotype\txTea_pass.bp_dis\txTea_pass.te_dis\txTea_pass.te_identity\txTea_pass.te_len\txTea_pass.div\txTea_pass.genotype\tPALMER.bp_dis\tPALMER.te_dis\tPALMER.te_identity\tPALMER.te_len\tPALMER.div\tPALMER.genotype\tTELR.bp_dis\tTELR.te_dis\tTELR.te_identity\tTELR.te_len\tTELR.div\tTELR.genotype\tTrEMOLO.bp_dis\tTrEMOLO.te_dis\tTrEMOLO.te_identity\tTrEMOLO.te_len\tTrEMOLO.div\tTrEMOLO.genotype\tMEHunter.bp_dis\tMEHunter.te_dis\tMEHunter.te_identity\tMEHunter.te_len\tMEHunter.div\tMEHunter.genotype\tGraffiTE.bp_dis\tGraffiTE.te_dis\tGraffiTE.te_identity\tGraffiTE.te_len\tGraffiTE.div\tGraffiTE.genotype\tchrom\tstart\tend" > ${intersection_data_path}/all.merge.dis.txt


for SAMPLE in simulation_germ_ccs simulation_germ_clr simulation_germ_ont # simulation_germ_clr simulation_germ_ccs # simulation_soma_ont simulation_soma_clr simulation_soma_ccs #   # simulation_soma # HG001 # HG002 # SRR9685183.sra
do
    for depth in 1 2 3 4 5 10 20 30 40 50
    do 

        DATA_PATH="${intersection_data_path}/${depth}/intersection/dis"

        awk '{split($6,a,";");LOCATE_ALL=".";LOCATE_PASS=".";MELT=".";PALMER=".";TEMP2=".";TLDR_ALL=".";TLDR_PASS=".";xTea_ALL=".";xTea_PASS=".";PALMER=".";TELR=".";TrEMOLO=".";MEHunter=".";GraffiTE=".";
            for(i in a){if(a[i]~/LOCATE-all/){LOCATE_ALL=a[i]};if(a[i]~/LOCATE-pass/){ LOCATE_PASS=a[i] };if(a[i]~/MELT/){MELT=a[i]};if(a[i]~/PALMER/){PALMER=a[i]};if(a[i]~/TEMP2/){TEMP2=a[i]};if(a[i]~/TLDR-all/){TLDR_ALL=a[i]};if(a[i]~/TLDR-pass/){TLDR_PASS=a[i]};
                        if(a[i]~/xTea-all/){xTea_ALL=a[i]};if(a[i]~/xTea-pass/){xTea_PASS=a[i]};if(a[i]~/TELR/){TELR=a[i]};if(a[i]~/TrEMOLO/){TrEMOLO=a[i]};if(a[i]~/MEHunter/){MEHunter=a[i]};if(a[i]~/GraffiTE/){GraffiTE=a[i]} };
            split($5,b,"_"); insID="ID"NR; print $1,$2,$3,$4,b[4],b[5],LOCATE_ALL,LOCATE_PASS,MELT,TEMP2,TLDR_ALL,TLDR_PASS,xTea_ALL,xTea_PASS,PALMER,TELR,TrEMOLO,MEHunter,GraffiTE,insID}' ${DATA_PATH}/${SAMPLE}.merge.bp_dis.bed \
            | awk -v sample=${SAMPLE} 'BEGIN{OFS="\t"}{split($5,a,":");split(a[1],b,"~"); TE=b[2]; type=b[3]; depth=$4; frequency=$6;insID="ID"NR;
                                    if($7=="."){LOCATE_all_t="-_-_-_-_-_-_-_-_-_-"}else{LOCATE_all_t=$7};split(LOCATE_all_t,LOCATE_ALL,"_");
                                    if($8=="."){LOCATE_pass_t="-_-_-_-_-_-_-_-_-_-"}else{LOCATE_pass_t=$8};split(LOCATE_pass_t,LOCATE_pass,"_");
                                    if($9=="."){MELT_t="-_-_-_-_-_-_-_-_-_-"}else{MELT_t=$9};split(MELT_t,MELT,"_");
                                    if($10=="."){TEMP2_t="-_-_-_-_-_-_-_-_-_-"}else{TEMP2_t=$10};split(TEMP2_t,TEMP2,"_");
                                    if($11=="."){TLDR_all_t="-_-_-_-_-_-_-_-_-_-"}else{TLDR_all_t=$11};split(TLDR_all_t,TLDR_all,"_");
                                    if($12=="."){TLDR_pass_t="-_-_-_-_-_-_-_-_-_-"}else{TLDR_pass_t=$12};split(TLDR_pass_t,TLDR_pass,"_");
                                    if($13=="."){xTea_all_t="-_-_-_-_-_-_-_-_-_-"}else{xTea_all_t=$13};split(xTea_all_t,xTea_all,"_");
                                    if($14=="."){xTea_pass_t="-_-_-_-_-_-_-_-_-_-"}else{xTea_pass_t=$14};split(xTea_pass_t,xTea_pass,"_");
                                    if($15=="."){PALMER_t="-_-_-_-_-_-_-_-_-_-"}else{PALMER_t=$15};split(PALMER_t,PALMER,"_");
                                    if($16=="."){TELR_t="-_-_-_-_-_-_-_-_-_-"}else{TELR_t=$16};split(TELR_t,TELR,"_");
                                    if($17=="."){TrEMOLO_t="-_-_-_-_-_-_-_-_-_-"}else{TrEMOLO_t=$17};split(TrEMOLO_t,TrEMOLO,"_");
                                    if($18=="."){MEHunter_t="-_-_-_-_-_-_-_-_-_-"}else{MEHunter_t=$18};split(MEHunter_t,MEHunter,"_");
                                    if($19=="."){GraffiTE_t="-_-_-_-_-_-_-_-_-_-"}else{GraffiTE_t=$19};split(GraffiTE_t,GraffiTE,"_");


                                    print insID,sample,depth,TE,type,frequency,
                                        LOCATE_ALL[2],LOCATE_ALL[3],LOCATE_ALL[4],LOCATE_ALL[8],LOCATE_ALL[9],LOCATE_ALL[10],
                                        LOCATE_pass[2],LOCATE_pass[3],LOCATE_pass[4],LOCATE_pass[8],LOCATE_pass[9],LOCATE_pass[10],
                                        MELT[2],MELT[3],MELT[4],MELT[8],MELT[9],MELT[10],
                                        TEMP2[2],TEMP2[3],TEMP2[4],TEMP2[8], TEMP2[9], TEMP2[10], 
                                        TLDR_all[2],TLDR_all[3],TLDR_all[4],TLDR_all[8],TLDR_all[9],TLDR_all[10],
                                        TLDR_pass[2],TLDR_pass[3],TLDR_pass[4],TLDR_pass[8],TLDR_pass[9],TLDR_pass[10],
                                        xTea_all[2],xTea_all[3],xTea_all[4],xTea_all[8],xTea_all[9],xTea_all[10],
                                        xTea_pass[2],xTea_pass[3],xTea_pass[4],xTea_pass[8],xTea_pass[9],xTea_pass[10],
                                        PALMER[2],PALMER[3],PALMER[4],PALMER[8],PALMER[9],PALMER[10],
                                        TELR[2],TELR[3],TELR[4],TELR[8],TELR[9],TELR[10],
                                        TrEMOLO[2],TrEMOLO[3],TrEMOLO[4],TrEMOLO[8],TrEMOLO[9],TrEMOLO[10],
                                        MEHunter[2],MEHunter[3],MEHunter[4],MEHunter[8],MEHunter[9],MEHunter[10],
                                        GraffiTE[2],GraffiTE[3],GraffiTE[4],GraffiTE[8],GraffiTE[9],GraffiTE[10],
                                        $1,$2,$3   }' >> ${intersection_data_path}/all.merge.dis.txt
    done
done

head ${intersection_data_path}/all.merge.dis.txt


