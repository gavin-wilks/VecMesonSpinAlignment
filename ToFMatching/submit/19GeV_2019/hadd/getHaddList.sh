#!/bin/bash

if [ $# -eq 2 ]
then
  jobid=$1
  step=$2 #recenter, shift, resolution, chargedflow
  rm EventPlane_hadd.list

  ls -d -1 /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu${energy}GeV_${year}/OutPut/SpinAlignment/${step}/*${jobid}*.root > EventPlane_hadd.list
fi
