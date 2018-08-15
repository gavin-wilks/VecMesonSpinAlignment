#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  ####submit macro#####
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-58 --time=2:00:00 ToFMatching.slr # phi TTree mode => 11 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-147 --time=2:00:00 ToFMatching.slr # phi TTree mode => 19 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-216 --time=2:00:00 ToFMatching.slr # phi TTree mode => 27 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-496 --time=2:00:00 ToFMatching.slr # phi TTree mode => 39 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-509 --time=2:00:00 ToFMatching.slr # phi TTree mode => 62 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-1613 --time=4:00:00 ToFMatching.slr # phi TTree mode => 200 GeV
  ####submit macro#####

  ####resubmit macro#####
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=2000-2000 --time=2:00:00 ToFMatching.slr # resubmit mode of phi TTree
  ####resubmit macro#####

  ####test macro#####
  sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-19 --time=2:00:00 ToFMatching.slr # test mode of phi TTree
  ####test macro#####
else
  echo "Wrong number of parameters!!"
fi

