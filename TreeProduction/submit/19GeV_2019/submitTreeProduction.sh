#!/bin/sh

codePath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/TreeProduction

##########Energy Selection##########
# energy=0  # 200GeV
# library=SL18h
# listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/200GeV_2014
# outPath=/star/data01/pwg/sunxuhit/AuAu200GeV_2014
 
# energy=1  # 54.0GeV
# library=SL18c
# listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/54GeV_2017
# outPath=/star/data01/pwg/sunxuhit/AuAu54GeV_2017

energy=4  # 19GeV
library=SL23d
listPath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/PidFlow/FileList/19p6GeV_2019
outPath=/gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019
##########Energy Selection##########

##########Mode Selection##########
#mode=0
#outDir=ReCenterParameter

#mode=1
#outDir=ShiftParameter

#mode=2
#outDir=Resolution

#mode=3
mode=4
#outDir=Phi/Forest
##########Mode Selection##########

##########Mixed Event Selection##########
flag_ME=$1 # 0 for SE | 1 for ME
SM=$2
#flag_ME=1 # 0
#SM=ME
##########Mixed Event Selection##########

flag_PID=0 # 0 for phi | 1 for rho | 2 for kstar


#outDir=Phi/Forest/Phi${SM}_eta1_20240724_pt0p9m2only
#outDir=Phi/Forest/Phi${SM}_eta1_20240904_TPConly
#outDir=Phi/Forest/Phi${SM}_eta1_20241021_PIDEff_DCA2_20Pt_24Eta_24Phi_AllCent
outDir=Phi/Forest/Phi${SM}_eta1p5_20250102_PIDEff_DCA2_20Pt_36Eta_24Phi_AllCent
#outDir=Phi/Forest/Phi_KaonPID_20240724

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/SpinAlignment/${outDir}
mkdir -p ${outPath}/OutPut/SpinAlignment/${outDir}

##########Test Production##########
star-submit-template -template ${codePath}/submit/19GeV_2019/testProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,flag_PID=$flag_PID,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
#star-submit-template -template resubmitProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,flag_PID=$flag_PID,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir

##########Re-Submit##########
