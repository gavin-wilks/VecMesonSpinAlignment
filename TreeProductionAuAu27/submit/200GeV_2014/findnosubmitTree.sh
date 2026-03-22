#!/bin/bash
date

jobid=$1
part=Phi # Phi, KStar
SM=$2 # SE, ME
where=pwg
library=$3
prod=$4
lum=$5
folder=Phi${SM}_eta1_20240823_resubmit_${library}_${prod}_${lum}
numjobs=`ls \/gpfs01\/star\/pwg\/gwilks3\/JOBS\/list\/*${jobid}* | wc -l`
echo "${numjobs}"



rm missingjobs_${SM}_${library}_${prod}_${lum}.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`ls -d -1 /gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu200GeV_2014/OutPut/SpinAlignment/${part}/Forest/${folder}/file_200GeV_${part}_${SM}_${jobid}_${j}.root`" !=  "/gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu200GeV_2014/OutPut/SpinAlignment/${part}/Forest/${folder}/file_200GeV_${part}_${SM}_${jobid}_${j}.root" ]

then
  echo "${j}" | tee -a  missingjobs_${SM}_${library}_${prod}_${lum}.log
fi
done

