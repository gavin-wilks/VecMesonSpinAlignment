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

energy=4  # 19GeV
library=SL21c
listPath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/FileLists/19GeV_2019
outPath=/gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019
##########Energy Selection##########

##########Mode Selection##########
pid=0 #0 = phi
particle=Phi
cut=0 #idk
year=0 #idk
pt=$1
cent=$2
mode=$3
fitmode=$4
nJobs=$5
nEvents=12500
#mode=2
#pid=KStar

# mode=2
# outDir=Resolution

# mode=3
# outDir=Phi/Forest
##########Mode Selection##########

outDir=${particle}_yabs1_EPRes_yspectraflat_fixedEPRes_20240724_nov2_flatpT_flatRP_phiEff_PhiEmbed_MaxFunctionCalls50000

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/CosEff/${outDir}
mkdir -p ${outPath}/OutPut/CosEff/${outDir}

##########Test Production##########
star-submit-template -template testProductionTemp.xml -entities fitmode=$fitmode,pt=$pt,cent=$cent,mode=$mode,particle=$particle,pid=$pid,energy=$energy,nJobs=$nJobs,nEvents=$nEvents,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
#star-submit-template -template testKStarProductionTemp.xml -entities fitmode=$fitmode,pt=$pt,cent=$cent,mode=$mode,particle=$particle,pid=$pid,energy=$energy,nJobs=$nJobs,nEvents=$nEvents,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
#star-submit-template -template resubmitProductionTemp.xml -entities pid=$pid,mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########
