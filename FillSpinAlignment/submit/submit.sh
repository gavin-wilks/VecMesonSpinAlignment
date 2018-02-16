#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-9 --time=01:00:00 FillSpinAlignment.slr # test mode
  sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-32 --time=08:00:00 FillSpinAlignment.slr # production mode

  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-9 --time=02:00:00 FillSpinAlignmentZdcSmd.slr # test mode for ZDCSMD
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-32 --time=08:00:00 FillSpinAlignmentZdcSmd.slr # production mode for ZDCSMD
else
  echo "Wrong number of parameters!!"
fi

