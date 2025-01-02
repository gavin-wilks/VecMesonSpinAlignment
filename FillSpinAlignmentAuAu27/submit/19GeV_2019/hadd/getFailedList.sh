#!/bin/bash

if [ $# -eq 4 ]
then
  jobid=$1
  energy=$2 #19p6
  year=$3 #2019
  njobs=$4
  rm Error.list
  
  for((j=0; j<${njobs}; j++ ))
  do
    #${j} >> Error.list
    echo "${j}"
    grep "hadd exiting due to error in" /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu${energy}GeV_${year}/hadd/Log/RunQA/${jobid}_${j}.log >> Error.list
  done
fi
