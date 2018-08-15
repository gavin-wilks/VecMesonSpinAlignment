#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  ####submit macro#####
  sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-2377 --time=1:00:00 KaonEmbedding.slr # Kaon embedding => 11 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-147 --time=2:00:00 KaonEmbedding.slr # Kaon embedding => 19 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-216 --time=2:00:00 KaonEmbedding.slr # Kaon embedding => 27 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-496 --time=2:00:00 KaonEmbedding.slr # Kaon embedding => 39 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-509 --time=2:00:00 KaonEmbedding.slr # Kaon embedding => 62 GeV
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-1613 --time=4:00:00 KaonEmbedding.slr # Kaon embedding => 200 GeV
  ####submit macro#####

  ####resubmit macro#####
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=2000-2000 --time=1:00:00 KaonEmbedding.slr # resubmit mode of Kaon embedding
  ####resubmit macro#####

  ####test macro#####
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-19 --time=1:00:00 KaonEmbedding.slr # test mode of Kaon embedding
  ####test macro#####
else
  echo "Wrong number of parameters!!"
fi

