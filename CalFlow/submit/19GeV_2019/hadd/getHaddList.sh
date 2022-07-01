#!/bin/bash

if [ $# -eq 2 ]
then

SM=$1
job=$2
  rm hadd.list

  ls -d -1 /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/FlowYields/Phi${SM}/*${job}*.root > hadd.list
fi
