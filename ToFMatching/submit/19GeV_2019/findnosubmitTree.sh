#!/bin/bash
date

if [ $# -ne 3 ]
 then
  echo -e "\033[31m Please input your trigger, and try a again ! bye. \033[0m"
  exit 1
fi

jobid=$1
part=$2 # Phi, KStar
#SM=$3 # SE, ME
where=$3
numjobs=`ls JOBS\/list\/*${jobid}* | wc -l`
echo "${numjobs}"



rm missingjobs.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`find /gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/ToFMatching/${part}/ -name "*_${jobid}_${j}.root"`" !=  "/gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/ToFMatching/${part}/file_19GeV_ToFMatch_${jobid}_${j}.root" ]

then
  echo "${j}" | tee -a  missingjobs.log
fi
done

