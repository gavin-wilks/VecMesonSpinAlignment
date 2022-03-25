#!/bin/bash

if [ $# -eq 2 ]
part=$1
SM=$2
then
  rm hadd.list

  ls -d -1 /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignmentYields/${part}${SM}/*.root > hadd.list
fi
