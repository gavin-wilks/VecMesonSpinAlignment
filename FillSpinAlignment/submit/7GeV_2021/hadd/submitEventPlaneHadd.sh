#!/bin/sh

#codePath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/PidFlow/RunQA

##########Energy Selection##########
#energy=0  # 14p5GeV
#library=SL21c
#listPath=/star/u/gwilks3/Workspace/global_spin_alignment/PidFlow/FileList/14p5GeV_2019
#outPath=/gpfs01/star/scratch/gwilks3/global_spin_alignment/AuAu14p5GeV_2019
 
#energy=1  # 19p6GeV
#library=SL21c
listPath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/FillSpinAlignment/submit/7GeV_2021/hadd
outPath=/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu7GeV_2021
##########Energy Selection##########

##########Mode Selection##########
#mode=0
#outDir=ReCenterParameter
etamode=$1
part=Phi
SM=$2
order=$3

outDir=${part}${SM}_etamode${etamode}_order${order}_20240511
outFile=Yields_${part}_${SM}_7GeV
##########Mode Selection##########

mkdir -p ${outPath}/hadd/SpinAlignmentYields/${outDir}
mkdir -p ${outPath}/hadd/SpinAlignmentYieldsLog/${outDir}

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

##########Full Production##########
star-submit-template -template myhadd.xml -entities order=$order,etamode=${etamode},SM=$SM,outFile=$outFile,outPath=$outPath,outDir=$outDir,listPath=$listPath
##########Full Production##########

#star-submit-template -template myhadd_resubmit.xml -entities outFile=$outFile,outPath=$outPath,outDir=$outDir,listPath=$listPath
