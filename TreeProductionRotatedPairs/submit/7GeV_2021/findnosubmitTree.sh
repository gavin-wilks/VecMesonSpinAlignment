#!/bin/bash
date

jobid=$1
part=Phi # Phi, KStar
SM=$2 # SE, ME
where=pwg
folder=Phi${SM}
numjobs=`ls JOBS\/list\/*${jobid}* | wc -l`
echo "${numjobs}"



rm missingjobs.log
for((j=0; j<${numjobs}; j++ ))
do
if [ "`ls -d -1 /gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu7GeV_2019/OutPut/SpinAlignment/${part}/Forest/${folder}/file_7GeV_${part}_${SM}_${jobid}_${j}.root`" !=  "/gpfs01/star/${where}/gwilks3/VectorMesonSpinAlignment/AuAu7GeV_2019/OutPut/SpinAlignment/${part}/Forest/${folder}/file_7GeV_${part}_${SM}_${jobid}_${j}.root" ]

then
  echo "${j}" | tee -a  missingjobs.log
fi
done

