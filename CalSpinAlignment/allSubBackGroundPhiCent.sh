#!/bin/bash

# inputs ==> energy, eta, pid, year, date, random3D
#root -l -b -q subBackGroundPhiCent.C\(4,0,0,\"20230720\",0,1,\"eta0p4\"\)
#root -l -b -q calSpinAlignmentSysPhiCentPoly3Integration.C\(4,0,0,0,4,0,0,0,1,\"eta0p4\"\)
#root -l -b -q subBackGroundPhiCent.C\(4,0,0,\"20230720\",0,1,\"eta0p6\"\)
#root -l -b -q calSpinAlignmentSysPhiCentPoly3Integration.C\(4,0,0,0,4,0,0,0,1,\"eta0p6\"\)
#root -l -b -q subBackGroundPhiCent.C\(4,0,0,\"20230720\",0,1,\"eta0p8\"\)
#root -l -b -q calSpinAlignmentSysPhiCentPoly3Integration.C\(4,0,0,0,4,0,0,0,1,\"eta0p8\"\)
#root -l -b -q subBackGroundPhiCent.C\(4,0,0,\"20230720\",0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiCentPoly3Integration.C\(4,0,0,0,4,0,0,0,1,\"eta1_eta1\"\)

#root -l -b -q subBackGroundPhiCent.C\(3,0,0,\"20230803\",0,1,\"eta0p4\"\)
#root -l -b -q calSpinAlignmentSysPhiCentPoly3Integration.C\(3,0,0,0,4,0,0,0,1,\"eta0p4\"\)
#root -l -b -q subBackGroundPhiCent.C\(3,0,0,\"20230803\",0,1,\"eta0p6\"\)
#root -l -b -q calSpinAlignmentSysPhiCentPoly3Integration.C\(3,0,0,0,4,0,0,0,1,\"eta0p6\"\)
#root -l -b -q subBackGroundPhiCent.C\(3,0,0,\"20230803\",0,1,\"eta0p8\"\)
#root -l -b -q calSpinAlignmentSysPhiCentPoly3Integration.C\(3,0,0,0,4,0,0,0,1,\"eta0p8\"\)
#root -l -b -q subBackGroundPhiCent.C\(3,0,0,\"20230803\",0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiCentPoly3Integration.C\(3,0,0,0,4,0,0,0,1,\"eta1_eta1\"\)

mkdir -p figures/Phi/7GeV/centralitystudy
mkdir -p figures/Phi/14GeV/centralitystudy
mkdir -p figures/Phi/19GeV/centralitystudy

#root -l -b -q subBackGroundPhiCent.C\(3,0,0,\"20230803\",0,1,\"eta1_eta1\"\)
#root -l -b -q subBackGroundPhiCent.C\(3,0,0,\"20230625\",0,2,\"eta1_eta1\"\)
#root -l -b -q subBackGroundPhiCent.C\(4,0,0,\"20230720\",0,1,\"eta1_eta1\"\)
#root -l -b -q subBackGroundPhiCent.C\(4,0,0,\"20230531\",0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiCentPoly3Integration.C\(3,0,0,0,4,0,0,0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiCentPoly3Integration.C\(3,0,0,0,4,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiCentPoly3Integration.C\(4,0,0,0,4,0,0,0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiCentPoly3Integration.C\(4,0,0,0,4,0,0,0,2,\"eta1_eta1\"\)

root -l -b -q subBackGroundPhiCent.C\(0,0,0,\"20240511\",0,1,\"eta1_eta1\"\)
root -l -b -q subBackGroundPhiCent.C\(0,0,0,\"20240511\",0,2,\"eta1_eta1\"\)
root -l -b -q calSpinAlignmentSysPhiCent_PolySys.C\(0,0,0,0,4,0,0,0,1,\"eta1_eta1\"\)
root -l -b -q calSpinAlignmentSysPhiCent_PolySys.C\(0,0,0,0,4,0,0,0,2,\"eta1_eta1\"\)

#root -l -b -q calSpinAlignmentSysPhiCent_PolySys.C\(3,0,0,0,4,0,0,0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiCent_PolySys.C\(3,0,0,0,4,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiCent_PolySys.C\(4,0,0,0,4,0,0,0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhiCent_PolySys.C\(4,0,0,0,4,0,0,0,2,\"eta1_eta1\"\)
