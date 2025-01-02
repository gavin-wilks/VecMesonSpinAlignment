#!/bin/sh

#codePath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/PidFlow/RunQA

##########Energy Selection##########
#energy=0  # 14p5GeV
#library=SL21c
#listPath=/star/u/gwilks3/Workspace/global_spin_alignment/PidFlow/FileList/14p5GeV_2019
#outPath=/gpfs01/star/scratch/gwilks3/global_spin_alignment/AuAu14p5GeV_2019
 
#energy=1  # 19p6GeV
#library=SL21c
listPath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/TreeProduction/submit/14GeV_2019/hadd
outPath=/gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu14GeV_2019
##########Energy Selection##########

##########Mode Selection##########
#mode=0
#outDir=ReCenterPar
#outDir=ShiftPar
outDir=Resolution
#outDir=EpdFlow
outFile=file_14GeV_${outDir}
##########Mode Selection##########

mkdir -p ${outPath}/hadd/TreeProductionLog/${outDir}
mkdir -p ${outPath}/hadd/TreeProduction/${outDir}

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

##########Full Production##########
star-submit-template -u ie -template myhadd.xml -entities outFile=$outFile,outPath=$outPath,outDir=$outDir,listPath=$listPath
##########Full Production##########

#star-submit-template -template myhadd_resubmit.xml -entities outFile=$outFile,outPath=$outPath,outDir=$outDir,listPath=$listPath
