#!/bin/sh

codePath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/ToFMatching

##########Energy Selection##########
# energy=0  # 200GeV
# library=SL18h
# listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/200GeV_2014
# outPath=/star/data01/pwg/sunxuhit/AuAu200GeV_2014
 
# energy=1  # 54.0GeV
# library=SL18c
# listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/54GeV_2017
# outPath=/star/data01/pwg/sunxuhit/AuAu54GeV_2017

energy=0  # 19GeV
library=SL22b
listPath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/PidFlow/FileList/7p7GeV_2021
outPath=/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu7GeV_2021
##########Energy Selection##########

##########Mode Selection##########
#mode=0
#outDir=ReCenterParameter

#mode=1
#outDir=ShiftParameter

#mode=2
#outDir=Resolution

#mode=3
#outDir=Phi/Forest
##########Mode Selection##########

##########Mixed Event Selection##########
#flag_ME=0 # 0 for SE | 1 for ME
#SM=SE
#flag_ME=1 # 0
#SM=ME
##########Mixed Event Selection##########

#flag_PID=2 # 0 for phi | 1 for rho | 2 for kstar
pid=0
outDir=ToFMatching/Phi_PhiEff_Fixed
#chunk=10

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/${outDir}
mkdir -p ${outPath}/OutPut/${outDir}

##########Test Production##########
#star-submit-template -template testProductionTemp.xml -entities pid=$pid,energy=$energy,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
star-submit-template -template resubmitProductionTemp.xml -entities pid=$pid,energy=$energy,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir

##########Re-Submit##########
