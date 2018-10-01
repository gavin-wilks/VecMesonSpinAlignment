#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  ####submit macro#####
  sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-58 --time=2:00:00 ToFMatching.slr # Kaon ToFMatching mode => 11 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-147 --time=2:00:00 ToFMatching.slr # Kaon ToFMatching mode => 19 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-216 --time=2:00:00 ToFMatching.slr # Kaon ToFMatching mode => 27 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-496 --time=2:00:00 ToFMatching.slr # Kaon ToFMatching mode => 39 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-509 --time=2:00:00 ToFMatching.slr # Kaon ToFMatching mode => 62 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-1613 --time=4:00:00 ToFMatching.slr # Kaon ToFMatching mode => 200 GeV
  ####submit macro#####

  ####resubmit macro#####
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=2000-2000 --time=2:00:00 ToFMatching.slr # resubmit mode of Kaon ToFMatching
  ####resubmit macro#####

  ####test macro#####
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-19 --time=2:00:00 ToFMatching.slr # test mode of Kaon ToFMatching
  ####test macro#####
else
  echo "Wrong number of parameters!!"
fi

