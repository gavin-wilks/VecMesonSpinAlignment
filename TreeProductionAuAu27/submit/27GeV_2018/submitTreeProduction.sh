#!/bin/sh

codePath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/TreeProductionAuAu27

##########Energy Selection##########
# energy=0  # 200GeV
# library=SL18h
# listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/200GeV_2014
# outPath=/star/data01/pwg/sunxuhit/AuAu200GeV_2014
 
# energy=1  # 54.0GeV
# library=SL18c
# listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/54GeV_2017
# outPath=/star/data01/pwg/sunxuhit/AuAu54GeV_2017

energy=5  # 19GeV
library=SL19b
listPath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/PidFlow/FileList/27GeV_2018
outPath=/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu27GeV_2018
##########Energy Selection##########

##########Mode Selection##########
#mode=0
#outDir=ReCenterParameter_Retry_NoVzVzVpdCut

#mode=1
#outDir=ShiftParameter_Retry_NoVzVzVpdCut

#mode=2
#outDir=Resolution_Retry_NoVzVzVpdCut

mode=3
#outDir=Phi/Forest
##########Mode Selection##########

##########Mixed Event Selection##########
#flag_ME=0 # 0 for SE | 1 for ME
#SM=SE
flag_ME=1 # 0
SM=ME
##########Mixed Event Selection##########

flag_PID=0 # 0 for phi | 1 for rho | 2 for kstar


outDir=Phi/Forest/Phi${SM}_eta1_20240712_WithVzVzVpdCut_NoPsi2Bin

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/SpinAlignment/${outDir}
mkdir -p ${outPath}/OutPut/SpinAlignment/${outDir}

##########Test Production##########
star-submit-template -template ${codePath}/submit/27GeV_2018/testProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,flag_PID=$flag_PID,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
#star-submit-template -template resubmitProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,flag_PID=$flag_PID,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir

##########Re-Submit##########