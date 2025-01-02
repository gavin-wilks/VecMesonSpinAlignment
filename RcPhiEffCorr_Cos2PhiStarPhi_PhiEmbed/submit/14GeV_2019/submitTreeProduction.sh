#!/bin/sh

codePath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/RcPhiEffCorr

##########Energy Selection##########
# energy=0  # 200GeV
# library=SL18h
# listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/200GeV_2014
# outPath=/star/data01/pwg/sunxuhit/AuAu200GeV_2014
 
# energy=1  # 54.0GeV
# library=SL18c
# listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/54GeV_2017
# outPath=/star/data01/pwg/sunxuhit/AuAu54GeV_2017

energy=3  # 14GeV
library=SL21c
listPath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/FileLists/14GeV_2019
outPath=/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu14GeV_2019
##########Energy Selection##########

##########Mode Selection##########
pid=0 #0 = phi
inputpt=$1
startpt=$2
stoppt=$3
mode=$4
etamode=$5
order=$6

#mode=2
#pid=KStar

# mode=2
# outDir=Resolution

# mode=3
# outDir=Phi/Forest
##########Mode Selection##########
#outDir=CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order${order}
outDir=CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order${order}_NoRapiditySpectra_FixedFirstEP

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/${outDir}
mkdir -p ${outPath}/OutPut/${outDir}

##########Test Production##########
star-submit-template -template testProductionTemp.xml -entities order=$order,pid=$pid,inputpt=$inputpt,startpt=$startpt,stoppt=$stoppt,mode=$mode,etamode=$etamode,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir,energy=$energy
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
#star-submit-template -template resubmitProductionTemp.xml -entities pid=$pid,mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########
