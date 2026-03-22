#!/bin/bash

if [ $# -eq 3 ]
part=$1
loc=$2
jobid=$3
then
  rm hadd.list

  ls -d -1 /gpfs01/star/${loc}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEff/${part}/*${jobid}*.root > hadd.list
fi
