#!/bin/sh

codePath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/AcceptanceCorrection

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
library=SL21c
listPath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/FileLists/19GeV_2019
outPath=/gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019
##########Energy Selection##########

##########Mode Selection##########
pid=0 #0 = phi
ptbin=$1
nJobs=500
nEvents=20000
rho00=$2
rapidity=$3

#mode=2
#pid=KStar

# mode=2
# outDir=Resolution

# mode=3
# outDir=Phi/Forest
##########Mode Selection##########

outDir=AcceptanceCorrection/VaryInput_rho00

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/${outDir}
mkdir -p ${outPath}/OutPut/${outDir}

##########Test Production##########
star-submit-template -template testProductionTemp.xml -entities rapidity=$rapidity,rho00=$rho00,ptbin=$ptbin,pid=$pid,energy=$energy,nJobs=$nJobs,nEvents=$nEvents,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
#star-submit-template -template resubmitProductionTemp.xml -entities pid=$pid,mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########
