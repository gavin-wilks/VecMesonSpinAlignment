#!/bin/bash


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
#root -l -b -q calSpinAlignmentSysPhiY_PolySys_TPCTOF.C\(4,0,0,0,6,0,0,0,0,0,0,2,\"eta1_eta1\"\)
root -l -b -q calSpinAlignmentSysPhiY_PolySys_TPCOnly.C\(4,0,0,0,6,0,0,0,0,0,0,2,\"eta1_eta1\"\)
