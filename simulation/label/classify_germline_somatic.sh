#!/bin/bash
# 分类 germline 和 somatic insertion

[ -z $FILELIST ] && echo "ERROR: FILELIST not set" && exit 1
[ -z $SOMATIC_FREQ ] && SOMATIC_FREQ=0.01

script_dir=`dirname $0`
classify_py="${script_dir}/classify_germline_somatic.py"

for i in $(seq 1 $(wc -l < $FILELIST))
do
  str=`sed "${i}q;d" $FILELIST`
  toks=($str)
  BAM=${toks[0]}
  WORKDIR=`dirname $BAM`
  OUTPUTDIR=${WORKDIR}/result_no_secondary
  WORKDIR=`dirname $WORKDIR`

  cd $OUTPUTDIR
  cat $WORKDIR/*/*summary > tmp.summary

  python $classify_py \
    --summary tmp.summary \
    --input TP.bed \
    --output tmp_TP.bed \
    --freq ${SOMATIC_FREQ}

  awk 'NR==FNR {freq[NR]=$5; next} {$5=freq[FNR]; print}' tmp_TP.bed TP_clt.txt > tmp_TP_clt.txt
  awk -v freq=${SOMATIC_FREQ} '$5>freq' tmp_TP_clt.txt > TP_clt_G.txt
  awk -v freq=${SOMATIC_FREQ} '$5<=freq' tmp_TP_clt.txt > TP_clt_S.txt

  rm tmp.summary tmp_TP.bed tmp_TP_clt.txt
  echo "Done: $OUTPUTDIR"
done
