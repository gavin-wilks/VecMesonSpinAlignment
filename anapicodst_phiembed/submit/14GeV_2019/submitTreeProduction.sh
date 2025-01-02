#!/bin/sh

codePath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/anapicodst

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

pid=Phi
mode=0 # 0 for phi, 1 for rho, 2 for KStar
particle=Kplus
gid=11 # pi+ = 8 | pi- = 9 | K+ = 11 | K- = 12   
##########Mode Selection##########

outDir=${pid}/${particle}_eta1

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/Embedding/${outDir}
mkdir -p ${outPath}/OutPut/Embedding/${outDir}

##########Test Production##########
#star-submit-template -u ie -template testProductionTemp.xml -entities mode=$mode,gid=$gid,particle=$particle,energy=$energy,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir 
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
star-submit-template -template resubmitProductionTemp.xml -entities mode=$mode,gid=$gid,particle=$particle,energy=$energy,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########
