#!/bin/bash

if [ $# -eq 2 ]
then
  pid=$1 # particle name (Phi, KStar, Rho)
  flag=$2 # SE, ME

  filename=${pid}_${flag}_Forest.list
  rm ${filename}

  ls -d -1 /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/Forest/${pid}${flag}/*_${pid}_${flag}_*.root > ${filename}
fi
