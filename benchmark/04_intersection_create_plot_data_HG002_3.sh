#!/bin/bash


user_path="${user_path:-${USER_PATH:-/path/to/user}}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOME="hg38"

# MEHunter.bp_dis\tMEHunter.te_dis\tMEHunter.te_identity\tMEHunter.te_len\tMEHunter.div
echo -e "insID\tsample\tdepth\tTE\ttype\tfrequency\tLOCATE.bp_dis\tLOCATE.te_dis\tLOCATE.te_identity\tLOCATE.te_len\tLOCATE.div\tTLDR.bp_dis\tTLDR.te_dis\tTLDR.te_identity\tTLDR.te_len\tTLDR.div\txTea.bp_dis\txTea.te_dis\txTea.te_identity\txTea.te_len\txTea.div\tchrom\tstart\tend" > ${user_path}/2022_long_reads/result/giab/intersection/intersection_revision/all.merge.dis.txt


for SAMPLE in giab # simulation_germ_clr simulation_germ_ccs # simulation_soma_ont simulation_soma_clr simulation_soma_ccs #   # simulation_soma # HG001 # HG002 # SRR9685183.sra
do
    for depth in hg38
    do 

        DATA_PATH="${user_path}/2022_long_reads/result/giab/intersection/intersection_revision/intersection/dis"

        awk '{split($6,a,";");LOCATE=".";TLDR=".";xTea=".";
            for(i in a){if(a[i]~/LOCATE/){LOCATE=a[i]};if(a[i]~/TLDR/){TLDR=a[i]};if(a[i]~/xTea/){xTea=a[i]};
            split($5,b,"_"); insID="ID"NR; print $1,$2,$3,$4,b[4],b[5],LOCATE,TLDR,xTea,insID}}' ${DATA_PATH}/${SAMPLE}.merge.bp_dis.bed \
            | awk -v sample=${SAMPLE} -v depth=${depth} 'BEGIN{OFS="\t"}{split($5,a,":");split(a[1],b,"~"); split(a[2],c,"~"); split(c[2],d,"_"); TE=d[1]; type="1"; frequency=$6;insID="ID"NR;
                                    if($7=="."){LOCATE_t="-_-_-_-_-_-_-_-_-"}else{LOCATE_t=$7};split(LOCATE_t,LOCATE,"_");
                                    if($8=="."){TLDR_t="-_-_-_-_-_-_-_-_-"}else{TLDR_t=$8};split(TLDR_t,TLDR,"_");
                                    if($9=="."){xTea_t="-_-_-_-_-_-_-_-_-"}else{xTea_t=$9};split(xTea_t,xTea,"_");


                                    print insID,sample,depth,TE,type,frequency,
                                        LOCATE[2],LOCATE[3],LOCATE[4],LOCATE[8],LOCATE[9],
                                        TLDR[2],TLDR[3],TLDR[4],TLDR[8],TLDR[9],
                                        xTea[2],xTea[3],xTea[4],xTea[8],xTea[9],
                                        $1,$2,$3   }' >> ${user_path}/2022_long_reads/result/giab/intersection/intersection_revision/all.merge.dis.txt
    done
done

head ${user_path}/2022_long_reads/result/giab/intersection/intersection_revision/all.merge.dis.txt


