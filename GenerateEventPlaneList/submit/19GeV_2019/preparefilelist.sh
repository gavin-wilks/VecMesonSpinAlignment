#!/bin/bash
date

if [ $# -eq 1 ]

jobid=$1

then
  echo ${jobid}
  rm resubmit.list
  
  for FILE in `cat missingjobs.log`
  do
       echo "This is working"
       echo ${FILE}
       cat /star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/GenerateEventPlaneList/submit/19GeV_2019/JOBS/list/sched${jobid}_${FILE}.list >> resubmit.list
  done
fi
