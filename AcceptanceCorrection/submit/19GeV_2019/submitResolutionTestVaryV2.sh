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
rho00=$1
rapidity=1.5
cent=$2
mode=$3
etamode=$4
res=$5
fullEP=0
v2=$6
v2val=$7
ptbin=$8 
nJobs=100
#EP_res=$7
nEvents=25000

#mode=2
#pid=KStar

# mode=2
# outDir=Resolution

# mode=3
# outDir=Phi/Forest
##########Mode Selection##########

#outDir=AcceptanceTupleProduction_fullEP${fullEP}_noDelta
#outDir=RerunPt_${res}_EtaCut${etamode}_0p2y0p4_v2${v2}_withBeta_MANYPLOTS_EDGE_v2val${v2val}
#outDir=AcceptanceQA_${res}_EtaCut${etamode}_v2study_0p8yabs1_v2val${v2val}
outDir=AcceptanceTupleProduction_subEP_withBeta_FirstOrder
mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/${outDir}
mkdir -p ${outPath}/OutPut/${outDir}

##########Test Production##########
star-submit-template -template testResolutionVaryV2.xml -entities v2val=$v2val,v2=$v2,fullEP=$fullEP,res=$res,etamode=$etamode,cent=$cent,mode=$mode,rapidity=$rapidity,rho00=$rho00,ptbin=$ptbin,pid=$pid,energy=$energy,nJobs=$nJobs,nEvents=$nEvents,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
#star-submit-template -template resubmitProductionTemp.xml -entities pid=$pid,mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########
