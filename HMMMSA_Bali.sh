#!/bin/bash

RV=RV11
HMMSUMPATH=$HOME/HMMSUM
tfaFiles=( $HMMSUMPATH/bb3_release/$RV/*.tfa )
outDIR=$HMMSUMPATH/HMM_$RV
#make the directory if it does not exist
mkdir $outDIR
for tfaFile in "${tfaFiles[@]}"
do
  echo "processing $tfaFile"
  #parse the filename from the path
  fileName="${tfaFile##*/}"
  name="${fileName%.*}"
  outFile="$outDIR/$name.msa"
  #echo $outFile
  #run HMMMSA
  $HMMSUMPATH/xHMMMSA2 $tfaFile $outFile
done
