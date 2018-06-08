#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  ####submit macro#####
  sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-99 --time=4:00:00 toyMcPhiDecay.slr
  ####submit macro#####

  ####resubmit macro#####
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=2000-2000 --time=2:00:00 toyMcPhiDecay.slr # resubmit mode of phi TTree
  ####resubmit macro#####

  ####test macro#####
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-9 --time=2:00:00 toyMcPhiDecay.slr # test mode of phi TTree
  ####test macro#####
else
  echo "Wrong number of parameters!!"
fi

