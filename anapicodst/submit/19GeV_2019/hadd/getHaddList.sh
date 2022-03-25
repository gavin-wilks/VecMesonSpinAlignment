#!/bin/bash

if [ $# -eq 1 ]
part=$1
then
  rm hadd.list

  ls -d -1 /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/Embedding/Phi/${part}/*.root > hadd.list
fi
