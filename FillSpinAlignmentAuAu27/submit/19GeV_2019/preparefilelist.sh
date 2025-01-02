#!/bin/bash
date

if [ $# -eq 5 ]

jobid=$1
library=$2
prod=$3
lum=$4
SM=$5

then
  echo ${jobid}
  rm resubmit_${SM}_${library}_${prod}_${lum}.list
  
  for FILE in `cat missingjobs_${SM}_${library}_${prod}_${lum}.log`
  do
       echo "This is working"
       echo ${FILE}
       cat /gpfs01/star/pwg/gwilks3/JOBS/list/sched${jobid}_${FILE}.list >> resubmit_${SM}_${library}_${prod}_${lum}.list
  done
fi
