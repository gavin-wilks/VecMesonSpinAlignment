#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  ####submit macro#####
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-58 --time=2:00:00 VecMesonTree.slr # phi TTree mode => 11 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-147 --time=2:00:00 VecMesonTree.slr # phi TTree mode => 19 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-216 --time=2:00:00 VecMesonTree.slr # phi TTree mode => 27 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-496 --time=2:00:00 VecMesonTree.slr # phi TTree mode => 39 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-509 --time=2:00:00 VecMesonTree.slr # phi TTree mode => 62 GeV
  sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-1613 --time=4:00:00 VecMesonTree.slr # phi TTree mode => 200 GeV
  ####submit macro#####

  ####resubmit macro#####
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=2000-2000 --time=2:00:00 VecMesonTree.slr # resubmit mode of phi TTree
  ####resubmit macro#####

  ####test macro#####
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-19 --time=2:00:00 VecMesonTree.slr # test mode of phi TTree
  ####test macro#####

  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-1613 --time=8:00:00 ZdcSmdTree.slr # ZDCSMD phi TTree mode
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=10000-10001 --time=8:00:00 ZdcSmdTree.slr # resubmit mode of ZDCSMD phi TTree 
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-23 --time=1:00:00 ZdcSmdTree.slr # test mode of ZDCSMD phi TTree
else
  echo "Wrong number of parameters!!"
fi

