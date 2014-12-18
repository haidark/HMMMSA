#!/bin/bash
RV=RV12
HMMSUMPATH=$HOME/HMMSUM
BALIPATH=$HMMSUMPATH/bb3_release/bali_score_src
msfFiles=( $HMMSUMPATH/bb3_release/$RV/*.msf )
msaDIR=$HMMSUMPATH/HMM_$RV
outFile=$msaDIR/scores.txt
rm $outFile
touch $outFile
for msfFile in "${msfFiles[@]}"
do
  echo "processing $msfFile"
  #parse the filename from the path
  fileName="${msfFile##*/}"
  name="${fileName%.*}"
  msaFile="$msaDIR/$name.msa"
  #run bali_score
  $BALIPATH/bali_score $msfFile $msaFile | tail -1 >> $outFile
done
