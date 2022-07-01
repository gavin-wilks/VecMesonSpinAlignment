#!/bin/bash

if [ $# -eq 2 ]
then
  jobid=$1
  energy=19 #19
  year=2019 #2019
  step=$2 #recenter, shift, resolution, chargedflow
  rm EventPlane_hadd.list

  find /gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu${energy}GeV_${year}/OutPut/SpinAlignment/${step}/ -name "*${jobid}*.root" > EventPlane_hadd.list
fi
