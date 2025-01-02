#!/bin/bash
date

jobid=$1
part=Phi
SM=SE
loc=pwg
folder=EtaMode0_20240124
numjobs=`ls JOBS\/list\/*${jobid}* | wc -l`
echo "${numjobs}"



rm missingjobs.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`ls -d -1 /gpfs01/star/${loc}/gwilks3/VectorMesonSpinAlignment/AuAu14GeV_2019/OutPut/SpinAlignmentYields/${part}${SM}_${folder}/EP_${part}_${SM}_14GeV_${jobid}_${j}.txt`" !=  "/gpfs01/star/${loc}/gwilks3/VectorMesonSpinAlignment/AuAu14GeV_2019/OutPut/SpinAlignmentYields/${part}${SM}_${folder}/EP_${part}_${SM}_14GeV_${jobid}_${j}.txt" ]

then
  echo "${j}" | tee -a  missingjobs.log
fi
done

