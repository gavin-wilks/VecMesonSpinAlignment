#!/bin/bash

if [ $# -eq 4 ]
then
  location=$1
  pid=$2
  SM=$3
  jobid=$4

  rm EventPlane_hadd.list

  find /gpfs01/star/${location}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/DirectFillEventMixing/${pid}/Forest/${SM}/*${jobid}* -name "*.root" -size +1k > EventPlane_hadd.list
fi
