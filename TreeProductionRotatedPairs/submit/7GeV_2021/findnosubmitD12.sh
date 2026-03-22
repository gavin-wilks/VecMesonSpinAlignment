#!/bin/bash
date

jobid=$1
part=Phi # Phi, KStar
#SM=$2 # SE, ME
where=pwg
#folder=$3
numjobs=`ls JOBS\/list\/*${jobid}* | wc -l`
epdstep=5
epdfolder=${epdstep}_RetryAgain
echo "${numjobs}"



rm missingjobs.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`ls -d -1 /gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu7GeV_2021/OutPut/SpinAlignment/Epd_Step${epdfolder}/file_7GeV_EpdCorrections_${epdstep}_${jobid}_${j}.root`" !=  "/gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu7GeV_2021/OutPut/SpinAlignment/Epd_Step${epdfolder}/file_7GeV_EpdCorrections_${epdstep}_${jobid}_${j}.root" ]

then
  echo "${j}" | tee -a  missingjobs.log
fi
done

