#!/bin/sh

codePath=/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/anapicodst_phiembed

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

pid=Phi
mode=0 # 0 for phi, 1 for rho, 2 for KStar
particle=Phi
gid=12 # pi+ = 8 | pi- = 9 | K+ = 11 | K- = 12   
##########Mode Selection##########

#rho=$1
#real=$2
#imag=$3
#re=$4
#im=$5
#hrho=$6
#hreal=$7
#himag=$8
#hre=$9
#him=$10
rho=0.3333
real=0.0
imag=0.0
re=0.0
im=0.0
hrho=0.3333
hreal=0.0
himag=0.0
hre=0.0
him=0.0
phiswitch=$1

#outDir=${pid}/${particle}_eta1_WithWeights_PrelimV2_20240820_RapidityIncluded_RC0NonNullMc_FixedHelicityAngle
#outDir=${pid}/${particle}_eta1_trees_notofmatch_EP_finerPID_20240925_48phi_phiswitch${phiswitch}_rho${rho}_real${real}_imag${imag}_re${re}_im${im}_hrho${hrho}_hreal${hreal}_himag${himag}_hre${hre}_him${him}
outDir=${pid}/${particle}_20241212_phiswitch${phiswitch}_noep_fixtrackmatching_addnsig
mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

mkdir -p ${outPath}/Log/Embedding/${outDir}
mkdir -p ${outPath}/OutPut/Embedding/${outDir}

##########Test Production##########
star-submit-template -template testProductionTemp.xml -entities phiswitch=$phiswitch,rho=$rho,real=$real,imag=$imag,re=$re,im=$im,hrho=$hrho,hreal=$hreal,himag=$himag,hre=$hre,him=$him,mode=$mode,gid=$gid,particle=$particle,energy=$energy,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
#star-submit-template -template resubmitProductionTemp.xml -entities mode=$mode,gid=$gid,particle=$particle,energy=$energy,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########
