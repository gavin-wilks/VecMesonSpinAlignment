#!/bin/bash

if [ $# -eq 3 ]
then
  pid=$1 # particle name (Phi, KStar, Rho)
  flag=$2 # SE, ME
  jobid=$3
  filename=${pid}_${flag}_Forest.list
  rm ${filename}

  ls -d -1 /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/${pid}/Forest/*_${pid}_${flag}_${jobid}*.root > ${filename}
  #ls -d -1 /gpfs01/star/scratch/gwilks3/${pid}/Forest/*_${pid}_${flag}*.root > ${filename}
fi
