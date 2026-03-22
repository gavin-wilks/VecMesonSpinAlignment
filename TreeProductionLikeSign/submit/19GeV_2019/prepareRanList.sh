#!/bin/bash
date

if [ $# -eq 1 ]

jobid=$1

then
  echo ${jobid}
  rm Run.list
  
  for FILE in `cat ranjobs.log`
  do
       echo "This is working"
       echo ${FILE}
       cat /star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/TreeProduction/submit/19GeV_2019/JOBS/list/sched${jobid}_${FILE}.list >> Run.list
  done
fi
