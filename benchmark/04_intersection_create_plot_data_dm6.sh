#!/bin/bash

user_path="${user_path:-${USER_PATH:-/path/to/user}}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOME="dm6"

# LOCATE  melt  TELR  TEMP2  tldr  TrEMOLO

echo -e "insID\tsample\tdepth\tTE\ttype\tfrequency\tLOCATE_all.bp_dis\tLOCATE_all.te_dis\tLOCATE_all.te_identity\tLOCATE_all.te_len\tLOCATE_all.div\tLOCATE_all.genotype\tLOCATE_pass.bp_dis\tLOCATE_pass.te_dis\tLOCATE_pass.te_identity\tLOCATE_pass.te_len\tLOCATE_pass.div\tLOCATE_pass.genotype\tMELT.bp_dis\tMELT.te_dis\tMELT.te_identity\tMELT.te_len\tMELT.div\tMELT.genotype\tTEMP2.bp_dis\tTEMP2.te_dis\tTEMP2.te_identity\tTEMP2.te_len\tTEMP2.div\tTEMP2.genotype\tTLDR_all.bp_dis\tTLDR_all.te_dis\tTLDR_all.te_identity\tTLDR_all.te_len\tTLDR_all.div\tTLDR_all.genotype\tTLDR_pass.bp_dis\tTLDR_pass.te_dis\tTLDR_pass.te_identity\tTLDR_pass.te_len\tTLDR_pass.div\tTLDR_pass.genotype\tTELR.bp_dis\tTELR.te_dis\tTELR.te_identity\tTELR.te_len\tTELR.div\tTELR.genotype\tTrEMOLO.bp_dis\tTrEMOLO.te_dis\tTrEMOLO.te_identity\tTrEMOLO.te_len\tTrEMOLO.div\tTrEMOLO.genotype\tGraffiTE.bp_dis\tGraffiTE.te_dis\tGraffiTE.te_identity\tGraffiTE.te_len\tGraffiTE.div\tGraffiTE.genotype\tchrom\tstart\tend" > ${user_path}/2022_long_reads/result/simulation/${GENOME}/intersection_revision/all.merge.dis.txt

for SAMPLE in simulation_germ_ont simulation_germ_clr simulation_germ_ccs # simulation_soma_ont simulation_soma_clr simulation_soma_ccs #   # simulation_soma # HG001 # HG002 # SRR9685183.sra
do
    for depth in 1 2 3 4 5 10 20 30 40 50
    do 
        DATA_PATH="${user_path}/2022_long_reads/result/simulation/${GENOME}/intersection_revision/${depth}/intersection/dis"

        awk '{split($6,a,";");LOCATE_all=".";LOCATE_PASS=".";MELT=".";TEMP2=".";TLDR_ALL=".";TLDR_PASS=".";TELR=".";TrEMOLO="."; GraffiTE=".";
            for(i in a){if(a[i]~/LOCATE-all/){LOCATE_all=a[i]};if(a[i]~/LOCATE-pass/){ LOCATE_PASS=a[i] };if(a[i]~/MELT/){MELT=a[i]};if(a[i]~/TEMP2/){TEMP2=a[i]};if(a[i]~/TLDR-all/){TLDR_ALL=a[i]};if(a[i]~/TLDR-pass/){TLDR_PASS=a[i]};if(a[i]~/TELR/){TELR=a[i]};if(a[i]~/TrEMOLO/){TrEMOLO=a[i]};if(a[i]~/GraffiTE/){GraffiTE=a[i]}  };
            split($5,b,"~"); split(b[3],c,"_");   print $1,$2,$3,$4,$5,c[2],LOCATE_all,LOCATE_PASS,MELT,TEMP2,TLDR_ALL,TLDR_PASS,TELR,TrEMOLO,GraffiTE}' ${DATA_PATH}/${SAMPLE}.merge.bp_dis.bed \
            | awk -v sample=${SAMPLE} 'BEGIN{OFS="\t"}{split($5,a,":");split(a[1],b,"~"); TE=b[2]; type=b[3]; depth=$4; frequency=$6;insID="ID"NR;
                                    if($7=="."){LOCATE_all_t="-_-_-_-_-_-_-_-_-_-"}else{LOCATE_all_t=$7};split(LOCATE_all_t,LOCATE_all,"_");
                                    if($8=="."){LOCATE_pass_t="-_-_-_-_-_-_-_-_-_-"}else{LOCATE_pass_t=$8};split(LOCATE_pass_t,LOCATE_pass,"_");
                                    if($9=="."){MELT_t="-_-_-_-_-_-_-_-_-_-"}else{MELT_t=$9};split(MELT_t,MELT,"_");
                                    if($10=="."){TEMP2_t="-_-_-_-_-_-_-_-_-_-"}else{TEMP2_t=$10};split(TEMP2_t,TEMP2,"_");
                                    if($11=="."){TLDR_all_t="-_-_-_-_-_-_-_-_-_-"}else{TLDR_all_t=$11};split(TLDR_all_t,TLDR_all,"_");
                                    if($12=="."){TLDR_pass_t="-_-_-_-_-_-_-_-_-_-"}else{TLDR_pass_t=$12};split(TLDR_pass_t,TLDR_pass,"_");
                                    if($13=="."){TELR_t="-_-_-_-_-_-_-_-_-_-"}else{TELR_t=$13};split(TELR_t,TELR,"_");
                                    if($14=="."){TrEMOLO_t="-_-_-_-_-_-_-_-_-_-"}else{TrEMOLO_t=$14};split(TrEMOLO_t,TrEMOLO,"_");
                                    if($15=="."){GraffiTE_t="-_-_-_-_-_-_-_-_-_-"}else{GraffiTE_t=$15};split(GraffiTE_t,GraffiTE,"_");


                                    print insID,sample,depth,TE,type,frequency,
                                        LOCATE_all[2],LOCATE_all[3],LOCATE_all[4],LOCATE_all[8],LOCATE_all[9],LOCATE_all[10],
                                        LOCATE_pass[2],LOCATE_pass[3],LOCATE_pass[4],LOCATE_pass[8],LOCATE_pass[9],LOCATE_pass[10],
                                        MELT[2],MELT[3],MELT[4],MELT[8],MELT[9],MELT[10],
                                        TEMP2[2],TEMP2[3],TEMP2[4],TEMP2[8], TEMP2[9], TEMP2[10], 
                                        TLDR_all[2],TLDR_all[3],TLDR_all[4],TLDR_all[8],TLDR_all[9],TLDR_all[10],
                                        TLDR_pass[2],TLDR_pass[3],TLDR_pass[4],TLDR_pass[8],TLDR_pass[9],TLDR_pass[10],
                                        TELR[2],TELR[3],TELR[4],TELR[8],TELR[9],TELR[10],
                                        TrEMOLO[2],TrEMOLO[3],TrEMOLO[4],TrEMOLO[8],TrEMOLO[9],TrEMOLO[10],
                                        GraffiTE[2],GraffiTE[3],GraffiTE[4],GraffiTE[8],GraffiTE[9],GraffiTE[10],
                                        $1,$2,$3  }' >> ${user_path}/2022_long_reads/result/simulation/${GENOME}/intersection_revision/all.merge.dis.txt
    done
done

head ${user_path}/2022_long_reads/result/simulation/${GENOME}/intersection_revision/all.merge.dis.txt


