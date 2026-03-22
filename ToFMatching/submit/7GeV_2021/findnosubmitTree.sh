#!/bin/bash
date

if [ $# -ne 2 ]
 then
  echo -e "\033[31m Please input your trigger, and try a again ! bye. \033[0m"
  exit 1
fi

jobid=$1
part=$2 # Phi, KStar
#SM=$3 # SE, ME
where=pwg
folder=PhiEff_Fixed
numjobs=`ls JOBS\/list\/*${jobid}* | wc -l`
echo "${numjobs}"



rm missingjobs.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`find /gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu7GeV_2021/OutPut/ToFMatching/${part}_${folder}/ -name "*_${jobid}_${j}.root"`" !=  "/gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu7GeV_2021/OutPut/ToFMatching/${part}_${folder}/file_7GeV_ToFMatch_${jobid}_${j}.root" ]

then
  echo "${j}" | tee -a  missingjobs.log
fi
done

