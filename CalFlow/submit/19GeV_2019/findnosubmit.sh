#!/bin/bash
date

if [ $# -ne 3 ]
 then
  echo -e "\033[31m Please input your trigger, and try a again ! bye. \033[0m"
  exit 1
fi

jobid=$1
part=$2
SM=$3
numjobs=`ls JOBS\/list\/*${jobid}* | wc -l`
echo "${numjobs}"



rm missingjobs.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`ls -d -1 /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/FlowYields/${part}${SM}_eta1p5_NoEToF/Yields_${part}_${SM}_19GeV_${jobid}_${j}.root`" !=  "/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/FlowYields/${part}${SM}_eta1p5_NoEToF/Yields_${part}_${SM}_19GeV_${jobid}_${j}.root" ]

then
  echo "${j}" | tee -a  missingjobs.log
fi
done

