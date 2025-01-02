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

energy=0  # 19GeV
library=SL22b
listPath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/FileLists/7GeV_2021
outPath=/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu7GeV_2021
##########Energy Selection##########

##########Mode Selection##########
pid=0 #0 = phi
inputpt=$1
startpt=$2
stoppt=$3
mode=$4
#mode=0
#etamode=$5
etamode=0
#order=$6
order=$5
rerho1n1=0.0
ysigma=1000
rho00=0.3333
real=0.0
imag=0.0
imrho1n1=0.0

#mode=2
#pid=KStar

# mode=2
# outDir=Resolution

# mode=3
# outDir=Phi/Forest
##########Mode Selection##########
#outDir=CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order${order}_NoRapiditySpectra_InputRho0p4_Fixed
outDir=CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order${order}_NoRapiditySpectra_FixedFirstEP_RCBins_officialv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240514_fixed2Dfunc_rerho1n1${rerho1n1}_ysigma${ysigma}_rho00${rho00}_r${real}_i${imag}_imrho1n1${imrho1n1}

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/${outDir}
mkdir -p ${outPath}/OutPut/${outDir}

##########Test Production##########
star-submit-template -template testProductionTemp.xml -entities real=$real,imag=$imag,imrho1n1=$imrho1n1,rho00=$rho00,ysigma=$ysigma,rerho1n1=$rerho1n1,order=$order,pid=$pid,inputpt=$inputpt,startpt=$startpt,stoppt=$stoppt,mode=$mode,etamode=$etamode,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir,energy=$energy
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
#star-submit-template -template resubmitProductionTemp.xml -entities pid=$pid,mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########