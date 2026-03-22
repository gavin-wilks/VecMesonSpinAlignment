#!/bin/bash
date

#if [ $# -ne 3 ]
# then
#  echo -e "\033[31m Please input your trigger, and try a again ! bye. \033[0m"
#  exit 1
#fi

jobid=$1
numjobs=`ls JOBS\/list\/*${jobid}* | wc -l`
echo "${numjobs}"
#step=$2
step=ReCenterParameter_TPCOnly
loc=pwg
mode=ReCenterPar

if [ ${step} == ReCenterParameter ] 
  then
  mode=ReCenterPar
fi

if [ ${step} == ShiftParameter ] 
  then
  mode=ShiftPar
fi

if [ ${step} == Resolution ] 
  then 
  mode=Resolution
fi


rm missingjobs.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`ls -d -1 /gpfs01/star/${loc}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/${step}/file_19GeV_${mode}_${jobid}_${j}.root`" !=  "/gpfs01/star/${loc}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/${step}/file_19GeV_${mode}_${jobid}_${j}.root" ]

then
  echo "${j}" | tee -a  missingjobs.log
fi
done

