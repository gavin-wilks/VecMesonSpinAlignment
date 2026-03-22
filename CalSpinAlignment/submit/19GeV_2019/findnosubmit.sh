#!/bin/bash
date

jobid=$1
part=Phi
SM=$2
loc=pwg
folder=EtaMode0_20240911_ProcessingKaonTTrees_PtYPhi
numjobs=`ls JOBS\/list\/*${jobid}* | wc -l`
echo "${numjobs}"



rm missingjobs.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`ls -d -1 /gpfs01/star/${loc}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/ProcessKaonTTrees/${part}${SM}_${folder}/output_${SM}_${jobid}_${j}.root`" !=  "/gpfs01/star/${loc}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/ProcessKaonTTrees/${part}${SM}_${folder}/output_${SM}_${jobid}_${j}.root" ]

then
  echo "${j}" | tee -a  missingjobs.log
fi
done

