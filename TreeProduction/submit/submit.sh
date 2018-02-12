#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-19 --time=8:00:00 VecMesonTree.slr # phi TTree mode
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=10000-10035 --time=8:00:00 VecMesonTree.slr # phi TTree resubmit mode

  sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-1613 --time=8:00:00 ZdcSmdTree.slr # ZDCSMD phi TTree mode
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=10000-10001 --time=8:00:00 ZdcSmdTree.slr # resubmit mode of ZDCSMD phi TTree 
  # sbatch -p shared-chos --account rhstar --ntasks=1 --array=0-23 --time=1:00:00 ZdcSmdTree.slr # test mode of ZDCSMD phi TTree
else
  echo "Wrong number of parameters!!"
fi

