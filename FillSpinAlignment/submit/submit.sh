#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  ####submit macro####
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-2 --time=08:00:00 FillSpinAlignment.slr # 11 GeV production mode
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-6 --time=08:00:00 FillSpinAlignment.slr # 19 GeV production mode
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-10 --time=08:00:00 FillSpinAlignment.slr # 27 GeV production mode
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-24 --time=08:00:00 FillSpinAlignment.slr # 39 GeV production mode
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-24 --time=08:00:00 FillSpinAlignment.slr # 62 GeV production mode
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-80 --time=08:00:00 FillSpinAlignment.slr # 200 GeV production mode
  ####submit macro####

  ####test macro####
  sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-3 --time=08:00:00 FillSpinAlignment.slr # test mode
  ####test macro####

  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-9 --time=02:00:00 FillSpinAlignmentZdcSmd.slr # test mode for ZDCSMD
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-32 --time=08:00:00 FillSpinAlignmentZdcSmd.slr # production mode for ZDCSMD
else
  echo "Wrong number of parameters!!"
fi

