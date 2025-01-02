#!/bin/bash
date

jobid=$1
part=Phi
SM=$2
loc=pwg
folder=EtaMode0_20240910_cent2060_kaonspectra_TTrees_TPCOnly_Retry
numjobs=`ls JOBS\/list\/*${jobid}* | wc -l`
echo "${numjobs}"



rm missingjobs.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`ls -d -1 /gpfs01/star/${loc}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignmentYields/${part}${SM}_${folder}/Yields_${part}_${SM}_19GeV_${jobid}_${j}.root`" !=  "/gpfs01/star/${loc}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignmentYields/${part}${SM}_${folder}/Yields_${part}_${SM}_19GeV_${jobid}_${j}.root" ]

then
  echo "${j}" | tee -a  missingjobs.log
fi
done

