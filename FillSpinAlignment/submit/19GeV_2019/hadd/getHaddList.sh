#!/bin/bash

if [ $# -eq 0 ]
then
  rm hadd.list

  ls -d -1 /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignmentYields/PhiSE/*.root > hadd.list
fi
