#!/bin/bash

if [ $# -eq 3 ]
part=$1
SM=$2
loc=$3
then
  rm hadd.list

  ls -d -1 /gpfs01/star/${loc}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignmentYields/${part}${SM}/*.root > hadd.list
fi
