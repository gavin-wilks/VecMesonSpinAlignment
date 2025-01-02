#!/bin/sh

pid=0
StartEvent=0
StopEvent=10000000000
# StopEvent=10000 # test
codePath=/star/u/sunxuhit/STAR/VecMesonSpinAlignment

##########Energy Selection##########
# energy=0  # 7.7GeV
# Energy=7GeV
 
# energy=1  # 11.5GeV
# outPath=/star/data01/pwg/sunxuhit/AuAu11GeV/
 
# energy=2  # 19.6GeV
# outPath=/star/data01/pwg/sunxuhit/AuAu19GeV/
 
# energy=3  # 27GeV
# outPath=/star/data01/pwg/sunxuhit/AuAu27GeV/
 
# energy=4  # 39GeV
# outPath=/star/data01/pwg/sunxuhit/AuAu39GeV/
 
# energy=5  # 62.4GeV
# outPath=/star/data01/pwg/sunxuhit/AuAu62GeV/
 
energy=6  # 200GeV
outPath=/star/data01/pwg/sunxuhit/AuAu200GeV/
##########Energy Selection##########

flag_ME=0 # 0 for SE | 1 for ME
SM=SE
# flag_ME=1 # 0
# SM=ME

cd JOBS

# for jobID in `seq 0 2` # 11.5 GeV
# for jobID in `seq 0 6` # 19.6 GeV
# for jobID in `seq 0 10` # 27 GeV
# for jobID in `seq 0 24` # 39 GeV
# for jobID in `seq 0 24` # 62.4 GeV
for jobID in `seq 0 80` # 200 GeV
do
  # echo $jobID
  # echo $energy
  # echo $outPath
  star-submit-template -template ../multiTreeReaderJobTempl.xml -entities energy=$energy,flag_ME=$flag_ME,jobID=$jobID,StartEvent=$StartEvent,StopEvent=$StopEvent,pid=$pid,SM=$SM,codePath=$codePath,outPath=$outPath
done
