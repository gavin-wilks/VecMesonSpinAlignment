#!/bin/sh

codePath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/RcPhiEffCorr_Cos2PhiStarPhi

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
ptset=$1
yset=$2
bincos=$3
binphi=$4
pid=0 #0 = phi
inputpt=$5
setnum=$6
startpt=$7
stoppt=5
#mode=$4
mode=0
#etamode=$5
etamode=0
#order=$6
order=2
rerho1n1=0.0
ysigma=1000
rho00=0.3333
real=0.0
imag=0.0
imrho1n1=0.0
helicity=0.3333


#mode=2
#pid=KStar

# mode=2
# outDir=Resolution

# mode=3
# outDir=Phi/Forest
##########Mode Selection##########
#outDir=CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order${order}_NoRapiditySpectra_InputRho0p4_Fixed
outDir=ptset_${ptset}_yset_${yset}_bincos${bincos}_binphi${binphi}_ptbin${inputpt}_set${setnum}_20250204_multiplejobs_v2onoff_fixedyweight_tpceffcentdep_explicitptselection_RC

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/${outDir}
mkdir -p ${outPath}/OutPut/${outDir}

##########Test Production##########
star-submit-template -template testProductionTemp.xml -entities setnum=$setnum,ptset=${ptset},yset=${yset},bincos=${bincos},binphi=${binphi},helicity=$helicity,real=$real,imag=$imag,imrho1n1=$imrho1n1,rho00=$rho00,ysigma=$ysigma,rerho1n1=$rerho1n1,order=$order,pid=$pid,inputpt=$inputpt,startpt=$startpt,stoppt=$stoppt,mode=$mode,etamode=$etamode,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir,energy=$energy
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
#star-submit-template -template resubmitProductionTemp.xml -entities pid=$pid,mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########
