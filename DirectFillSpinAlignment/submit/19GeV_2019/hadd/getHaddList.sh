#!/bin/bash

if [ $# -eq 1 ]
then
  location=scratch
  pid=KStar
  SM=$1
  rm EventPlane_hadd.list

  find /gpfs01/star/${location}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/DirectFill/${pid}/Forest/${SM}/ -name "*.root" -size +1k > EventPlane_hadd.list
fi
