#!/bin/bash

if [ $# -eq 3 ]
then
  pid=$1 # particle name (Phi, KStar, Rho)
  flag=$2 # SE, ME
  location=$3
  filename=${pid}_${flag}_Forest.list
  rm ${filename}

  find /gpfs01/star/${location}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/${pid}/Forest/${pid}${flag} -name *${pid}*${flag}*.root > ${filename}
  #ls -d -1 /gpfs01/star/scratch/gwilks3/${pid}/Forest/*_${pid}_${flag}*.root > ${filename}
fi
