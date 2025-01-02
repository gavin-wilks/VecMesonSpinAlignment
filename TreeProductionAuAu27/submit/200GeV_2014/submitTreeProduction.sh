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

energy=6  # 200GeV
library=$1
prod=$2
lum=$3
listPath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/PidFlow/FileList/200GeV_2014
outPath=/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu200GeV_2014
##########Energy Selection##########

##########Mode Selection##########
#mode=0
#outDir=ReCenterParameter_20240813

#mode=1
#outDir=ShiftParameter_20240814

#mode=2
#outDir=Resolution_20240815

mode=4
#outDir=Phi/Forest
##########Mode Selection##########

##########Mixed Event Selection##########
flag_ME=$4 # 0 for SE | 1 for ME
SM=$5
#flag_ME=1 # 0
#SM=ME
##########Mixed Event Selection##########

flag_PID=0 # 0 for phi | 1 for rho | 2 for kstar


outDir=Phi/Forest/Phi${SM}_eta1_20241107_EP

#mkdir -p JOBS/report
#mkdir -p JOBS/csh
#mkdir -p JOBS/list
mkdir -p /gpfs01/star/pwg/gwilks3/JOBS/files
mkdir -p /gpfs01/star/pwg/gwilks3/JOBS/report
mkdir -p /gpfs01/star/pwg/gwilks3/JOBS/csh
mkdir -p /gpfs01/star/pwg/gwilks3/JOBS/list

mkdir -p ${outPath}/Log/SpinAlignment/${outDir}_${library}_${prod}_${lum}
mkdir -p ${outPath}/OutPut/SpinAlignment/${outDir}_${library}_${prod}_${lum}

##########Test Production##########
star-submit-template -template ${codePath}/submit/200GeV_2014/testProductionTemp.xml -entities lum=$lum,prod=$prod,mode=$mode,energy=$energy,flag_ME=$flag_ME,flag_PID=$flag_PID,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
#star-submit-template -template resubmitProductionTemp.xml -entities lum=$lum,prod=$prod,mode=$mode,energy=$energy,flag_ME=$flag_ME,flag_PID=$flag_PID,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir

##########Re-Submit##########

rm schedTemplateExp.xml 
