#! /bin/bash

[ -z $FILELIST ] && FILELIST=your_path_to/simulation/simulate_training_data/dm6/filelist_examples/label_filelist_dm6
REPEAT='your_path_to/Dm6.rmsk.bed'
GAP='your_path_to/Dm6.gap.bed'
TE='your_path_to/GRCh38.transposon.fa'
TEMP3='your_path_to/TEMP3'
BLACKLIST=''
germModel=''
somaModel=''


for i in $(seq 1 900)
do
  str=`sed "${i}q;d" $FILELIST`
  toks=($str)
  BAM=${toks[0]}
  WORKDIR=`dirname $BAM`
  OUTPUTDIR=${WORKDIR}/result_no_secondary
  WORKDIR=`dirname $WORKDIR`
  [ ! -d $OUTPUTDIR ] && mkdir $OUTPUTDIR
  cd $OUTPUTDIR && echo "Changing wd: " $(pwd) && echo "Labeling for ${BAM}"

  your_path_to/simulation/label/label_protocol.sh \
    -b $BAM \
    -d $WORKDIR \
    -r $REPEAT \
    -g $GAP \
    -P 10 \
    -T 5 \
    --refTe $TE \
    --blacklist $BLACKLIST \
    --minLen 100 \
    --TEMP3 $TEMP3 \
    --subSize 5 \
    --germModel $germModel \
    --somaModel $somaModel

  echo "Done." && echo ""
done
