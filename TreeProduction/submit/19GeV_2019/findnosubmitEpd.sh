#!/bin/bash
date

jobid=$1
part=Phi # Phi, KStar
SM=$2 # SE, ME
where=pwg
#folder=$3
numjobs=`ls JOBS\/list\/*${jobid}* | wc -l`
echo "${numjobs}"



rm missingjobs.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`ls -d -1 /gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/${part}/Forest/FirstOrderPhi${SM}/file_19GeV_EpdFlow_4_${jobid}_${j}.root`" !=  "/gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/${part}/Forest/FirstOrderPhi${SM}/file_19GeV_EpdFlow_4_${jobid}_${j}.root" ]

then
  echo "${j}" | tee -a  missingjobs.log
fi
done

