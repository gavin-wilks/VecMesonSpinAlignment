#!/bin/bash

if [ $# -eq 4 ]
then
  jobid=$1
  energy=$2 #19
  year=$3 #2019
  step=$4 #recenter, shift, resolution, chargedflow
  rm EventPlane_hadd.list

  ls -d -1 /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu${energy}GeV_${year}/OutPut/SpinAlignment/${step}/*${jobid}*.root > EventPlane_hadd.list
fi
