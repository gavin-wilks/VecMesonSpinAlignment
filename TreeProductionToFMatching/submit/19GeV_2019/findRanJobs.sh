#!/bin/bash
date

if [ $# -ne 2 ]
 then
  echo -e "\033[31m Please input your trigger, and try a again ! bye. \033[0m"
  exit 1
fi

jobid=$1
numjobs=$2

rm ranjobs.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`ls -d -1 /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/Resolution/file_19GeV_Resolution_${jobid}_${j}.root`" == "/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/Resolution/file_19GeV_Resolution_${jobid}_${j}.root" ]
#if [ "`ls /gpfs01/star/scratch/gwilks3/limitedetaphi_timing/outtree${style}_${reco}/ | grep _${j}.root`" != "${Style}Tree_${j}.root" ]

then
  echo "${j}" | tee -a  ranjobs.log
fi
done

