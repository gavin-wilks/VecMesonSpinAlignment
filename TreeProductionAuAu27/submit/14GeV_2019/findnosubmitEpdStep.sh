#!/bin/bash
date

jobid=$1
part=Phi # Phi, KStar
#SM=$2 # SE, ME
where=scratch
#folder=$3
numjobs=`ls JOBS\/list\/*${jobid}* | wc -l`
step=5
echo "${numjobs}"



rm missingjobs.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`ls -d -1 /gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu14GeV_2019/OutPut/SpinAlignment/${part}/EpdStep${step}/file_14GeV_EpdCorrections_${step}_${jobid}_${j}.root`" !=  "/gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu14GeV_2019/OutPut/SpinAlignment/${part}/EpdStep${step}/file_14GeV_EpdCorrections_${step}_${jobid}_${j}.root" ]

then
  echo "${j}" | tee -a  missingjobs.log
fi
done

