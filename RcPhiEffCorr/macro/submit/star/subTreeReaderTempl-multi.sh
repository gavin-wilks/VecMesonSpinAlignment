#!/bin/sh

pid=0
NMAX=1000000
# NMAX=10000 # test
codePath=/star/u/sunxuhit/STAR/VecMesonSpinAlignment

##########Energy Selection##########
# energy=0  # 7.7GeV
# Energy=7GeV
 
# energy=1  # 11.5GeV
# year=1
# cut=0
# outPath=/star/data01/pwg/sunxuhit/AuAu11GeV/
 
# energy=2  # 19.6GeV
# year=0
# cut=0
# outPath=/star/data01/pwg/sunxuhit/AuAu19GeV/
 
# energy=3  # 27GeV
# year=0
# cut=0
# outPath=/star/data01/pwg/sunxuhit/AuAu27GeV/
 
# energy=4  # 39GeV
# year=1
# cut=0
# outPath=/star/data01/pwg/sunxuhit/AuAu39GeV/
 
# energy=5  # 62.4GeV
# year=1
# cut=0
# outPath=/star/data01/pwg/sunxuhit/AuAu62GeV/
 
energy=6  # 200GeV
year=0
cut=0
outPath=/star/data01/pwg/sunxuhit/AuAu200GeV/
##########Energy Selection##########

cd JOBS

counter=0
for jobID in `seq 1 100`
do
  # echo $jobID
  # echo $energy
  # echo $year
  # echo $cut
  # echo $outPath
  star-submit-template -template ../multiTreeReaderJobTempl.xml -entities energy=$energy,year=$year,cut=$cut,pid=$pid,NMAX=$NMAX,jobID=$jobID,codePath=$codePath,outPath=$outPath
done
