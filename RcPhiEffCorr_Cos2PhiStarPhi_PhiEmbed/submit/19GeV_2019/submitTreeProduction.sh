#!/bin/sh

codePath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/RcPhiEffCorr_Cos2PhiStarPhi_PhiEmbed

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
outPath=/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019
##########Energy Selection##########

##########Mode Selection##########
pid=0 #0 = phi
inputpt=$1
startpt=$2
stoppt=$3
#mode=$4
mode=0
#etamode=$5
etamode=0
#order=$6
order=2
ysigma=1000
rho00=0.3333
rerho1n1=0.0
imrho1n1=0.0
real=0.0
imag=0.0
hrho00=$4
hrerho1n1=$5
himrho1n1=$6
hreal=$7
himag=$8
v2=$9

#mode=2
#pid=KStar

# mode=2
# outDir=Resolution

# mode=3
# outDir=Phi/Forest
##########Mode Selection##########
#outDir=CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order${order}_NoRapiditySpectra_InputRho0p4_Fixed
outDir=CEPT_O${order}_NoYSpec_v2${v2}_phieff_phiTPC_20240815he_ysig${ysigma}_rho${rho00}_rerho${rerho1n1}_imrho${imrho1n1}_r${real}_i${imag}_hrho${hrho00}_hrerho${hrerho1n1}_himrho${himrho1n1}_hr${hreal}_hi${himag}

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/${outDir}
mkdir -p ${outPath}/OutPut/${outDir}

##########Test Production##########
star-submit-template -template testProductionTemp.xml -entities v2=$v2,hrerho1n1=$hrerho1n1,himrho1n1=$himrho1n1,hreal=$hreal,himag=$himag,hrho00=$hrho00,real=$real,imag=$imag,imrho1n1=$imrho1n1,rho00=$rho00,ysigma=$ysigma,rerho1n1=$rerho1n1,order=$order,pid=$pid,inputpt=$inputpt,startpt=$startpt,stoppt=$stoppt,mode=$mode,etamode=$etamode,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir,energy=$energy
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
#star-submit-template -template resubmitProductionTemp.xml -entities pid=$pid,mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########
