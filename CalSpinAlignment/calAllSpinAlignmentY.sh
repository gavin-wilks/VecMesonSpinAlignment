#!/bin/bash

#for i in 0 1 2 3 4 
#do 
#  root -l -b -q oldsys_calSpinAlignmentSysPhiPt_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,1,\"eta1_eta1\",0,${i}\)
#  root -l -b -q oldsys_calSpinAlignmentSysPhiPt_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,2,\"eta1_eta1\",0,${i}\)
#  #root -l -b -q calSpinAlignmentSysPhiCent_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,1,\"eta1_eta1\",0,${i}\)
#  #root -l -b -q calSpinAlignmentSysPhiCent_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,2,\"eta1_eta1\",0,${i}\)
#  #root -l -b -q calSpinAlignmentSysPhiRapidity_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,1,\"eta1_eta1\",0,${i}\)
#  #root -l -b -q calSpinAlignmentSysPhiRapidity_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,2,\"eta1_eta1\",0,${i}\)
#done
#for i in {0..9}
flag=_m105_2Gamma

for i in {0..0}
do 
  #root -l -b -q calSpinAlignmentSysPhiInt_Global_2D_OffDiag_Sys_Exp_OnlyBWInt_Voigt.C\(4,0,0,0,2,\"eta1_eta1\",0,${i},\"${flag}\"\)
  root -l -b -q calSpinAlignmentSysPhiInt_Global_2D_OffDiag_Sys_Exp_OnlyBWInt.C\(4,0,0,0,2,\"eta1_eta1\",0,${i},\"${flag}\"\)
  root -l -b -q calSpinAlignmentSysPhiInt_Global_2D_OffDiag_Sys_Exp_OnlyBWInt.C\(4,0,0,0,1,\"eta1_eta1\",0,${i},\"${flag}\"\)
  #root -l -b -q calSpinAlignmentSysPhiInt_Global_2D_OffDiag_Sys_Exp_OnlyBWInt_Voigt.C\(4,0,0,0,1,\"eta1_eta1\",0,${i},\"${flag}\"\)
  #root -l -b -q calSpinAlignmentSysPhiInt_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,2,\"eta1_eta1\",0,${i},\"${flag}\"\)
  #root -l -b -q calSpinAlignmentSysPhiInt_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,1,\"eta1_eta1\",0,${i},\"${flag}\"\)
  #root -l -b -q calSpinAlignmentSysPhiPt_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,1,\"eta1_eta1\",0,${i},\"${flag}\"\)
  #root -l -b -q calSpinAlignmentSysPhiPt_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,2,\"eta1_eta1\",0,${i},\"${flag}\"\)
  #root -l -b -q calSpinAlignmentSysPhiPt_Global_2D_OffDiag_Sys_Exp_Rotated.C\(4,0,0,0,1,\"eta1_eta1\",0,${i}\)
  #root -l -b -q calSpinAlignmentSysPhiPt_Global_2D_OffDiag_Sys_Exp_Rotated.C\(4,0,0,0,2,\"eta1_eta1\",0,${i}\)
  #root -l -b -q calSpinAlignmentSysPhiPt_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,1,\"eta1_eta1\",0,${i}\)
  #root -l -b -q calSpinAlignmentSysPhiPt_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,2,\"eta1_eta1\",0,${i}\)
  #root -l -b -q calSpinAlignmentSysPhiPt_Global_2D_OffDiag_Sys_Exp_Voigt.C\(4,0,0,0,2,\"eta1_eta1\",0,${i}\)
  #root -l -b -q calSpinAlignmentSysPhiPt_Global_2D_OffDiag_Sys_Exp_Voigt.C\(4,0,0,0,1,\"eta1_eta1\",0,${i}\)
  #root -l -b -q calSpinAlignmentSysPhiPt_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,2,\"eta1_eta1\",0,${i}\)
  #root -l -b -q calSpinAlignmentSysPhiCent_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,1,\"eta1_eta1\",0,${i}\)
  #root -l -b -q calSpinAlignmentSysPhiCent_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,2,\"eta1_eta1\",0,${i}\)
  #root -l -b -q calSpinAlignmentSysPhiRapidity_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,1,\"eta1_eta1\",0,${i}\)
  #root -l -b -q calSpinAlignmentSysPhiRapidity_Global_2D_OffDiag_Sys_Exp.C\(4,0,0,0,2,\"eta1_eta1\",0,${i}\)
done

#root -l -b -q calSpinAlignmentSysPhiEtaIntegratedPoly3.C\(3,0,0,0,6,1,2,0,0,0,0,1,\"eta0p4\"\)
#root -l -b -q calSpinAlignmentSysPhiEtaIntegratedPoly3.C\(3,0,0,0,6,1,2,0,0,0,0,2,\"eta0p6\"\)
#root -l -b -q calSpinAlignmentSysPhiEtaIntegratedPoly3.C\(3,0,0,0,6,1,2,0,0,0,0,1,\"eta0p8\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys.C\(3,0,0,0,6,1,2,0,0,0,0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys.C\(3,0,0,0,6,1,2,0,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys.C\(4,0,0,0,6,1,2,0,0,0,0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys.C\(4,0,0,0,6,1,2,0,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiEtaIntegratedPoly3.C\(3,0,0,0,6,1,2,0,0,0,0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiEtaIntegratedPoly3.C\(3,0,0,0,6,1,2,0,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiEtaIntegratedPoly3.C\(4,0,0,0,6,1,2,0,0,0,0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiEtaIntegratedPoly3.C\(4,0,0,0,6,1,2,0,0,0,0,2,\"eta1_eta1\"\)


#root -l -b -q calSpinAlignmentSysPhiY_PolySys.C\(0,0,0,0,6,1,2,0,0,0,0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys.C\(0,0,0,0,6,1,2,0,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys.C\(4,0,0,0,6,1,2,0,0,0,0,1,\"eta1_eta1\"\)


#void calSpinAlignmentSysPhiY_PolySys_TPCTOF(int energy = 4, int pid = 0, int year = 0, bool random3D = false, int etaQA = 9, int ptQA = 1, int centQA = 2, int dcaQA = 0, int nsigQA = 0, int normQA = 0, int sigQA = 0, int order = 2, std::string etamode = "eta1_eta1")
##root -l -b -q calSpinAlignmentSysPhiY_PolySys_TPCTOF_LargeEdgeBin.C\(4,0,0,0,6,0,0,0,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys_TPCTOF_LargeEdgeBin_DipAngleComp.C\(4,0,0,0,6,0,0,0,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys_TPCTOF_LargeEdgeBin_Voigt.C\(4,0,0,0,6,0,0,0,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys_TPCTOF.C\(4,0,0,0,6,0,0,0,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys_TPCOnly.C\(4,0,0,0,6,0,0,0,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys_TPCOnly.C\(4,0,0,0,6,0,0,0,0,0,0,2,\"eta1p5_eta1p5\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys_TPCTOF.C\(4,0,0,0,6,0,0,0,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiY_PolySys_TPCTOF.C\(4,0,0,0,6,0,0,0,0,0,0,2,\"eta1p5_eta1p5\"\)

#root -l -b -q calSpinAlignmentSysPhiRapidity_Global_2D_OffDiag.C\(4,0,0,0,2,\"eta1_eta1\",0\)
#root -l -b -q calSpinAlignmentSysPhiRapidity_Global_2D_OffDiag.C\(4,0,0,0,1,\"eta1_eta1\",0\)

#root -l -b -q calSpinAlignmentSysPhiPt_Global_2D_OffDiag.C\(4,0,0,0,2,\"eta1_eta1\",0\)

#root -l -b -q calSpinAlignmentSysPhiPt_Global_1D_OffDiag.C\(4,0,0,0,2,\"eta1_eta1\",0\)
#root -l -b -q calSpinAlignmentSysPhiPt_Global_1D_OffDiag.C\(4,0,0,1,1,\"eta1_eta1\",0\)

#root -l -b -q calSpinAlignmentSysPhiPt_Global_2D_OffDiag_NonFixedWidth.C\(4,0,0,0,2,\"eta1_eta1\",0\)
#root -l -b -q calSpinAlignmentSysPhiPt_Global_2D_OffDiag_NonFixedWidth.C\(4,0,0,0,1,\"eta1_eta1\",0\)

#root -l -b -q calSpinAlignmentSysPhiCent_Global_2D_OffDiag.C\(4,0,0,0,2,\"eta1_eta1\",0\)
#root -l -b -q calSpinAlignmentSysPhiCent_Global_2D_OffDiag.C\(4,0,0,0,1,\"eta1_eta1\",0\)
#root -l -b -q calSpinAlignmentSysPhiCent_Global_1D_OffDiag.C\(4,0,0,0,2,\"eta1_eta1\",0\)
#root -l -b -q calSpinAlignmentSysPhiCent_Global_1D_OffDiag.C\(4,0,0,1,1,\"eta1_eta1\",0\)

#root -l -b -q calSpinAlignmentSysPhiRapidity_Global_2D_OffDiag.C\(4,0,0,0,2,\"eta1_eta1\",0\)
#root -l -b -q calSpinAlignmentSysPhiRapidity_Global_2D_OffDiag.C\(4,0,0,0,1,\"eta1_eta1\",0\)
#root -l -b -q calSpinAlignmentSysPhiRapidity_Global_1D_OffDiag.C\(4,0,0,0,2,\"eta1_eta1\",0\)
#root -l -b -q calSpinAlignmentSysPhiRapidity_Global_1D_OffDiag.C\(4,0,0,1,1,\"eta1_eta1\",0\)


#root -l -b -q calSpinAlignmentSysPhiPt_Global_1D_OffDiag.C\(4,0,0,1,1,\"eta1_eta1\",0\)
#root -l -b -q calSpinAlignmentSysPhiRapidity_Global_1D_OffDiag.C\(4,0,0,1,1,\"eta1_eta1\",0\)
#root -l -b -q calSpinAlignmentSysPhiCent_Global_1D_OffDiag.C\(4,0,0,1,1,\"eta1_eta1\",0\)
