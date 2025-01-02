#!/bin/bash

#root -l -b -q subBackGround.C\(4,0,0,\"20230531\",0,2,\"eta0p4\"\)
#root -l -b -q subBackGround.C\(4,0,0,\"20230531\",0,2,\"eta0p6\"\)
#root -l -b -q subBackGround.C\(4,0,0,\"20230531\",0,2,\"eta0p8\"\)
#root -l -b -q subBackGround.C\(4,0,0,\"20230531\",0,2,\"eta1_eta1\"\)
#root -l -b -q subBackGround.C\(3,0,0,\"20230625\",0,2,\"eta1_eta1\"\)
#root -l -b -q subBackGround.C\(4,0,0,\"20230720\",0,1,\"eta0p4\"\)
#root -l -b -q subBackGround.C\(4,0,0,\"20230720\",0,1,\"eta0p6\"\)
#root -l -b -q subBackGround.C\(4,0,0,\"20230720\",0,1,\"eta0p8\"\)
#root -l -b -q subBackGround.C\(4,0,0,\"20230720\",0,1,\"eta1_eta1\"\)


#root -l -b -q calSpinAlignmentSysPhi.C\(4,0,0,0,2,\"eta0p4\"\)
#root -l -b -q calSpinAlignmentSysPhi.C\(4,0,0,0,2,\"eta0p6\"\)
#root -l -b -q calSpinAlignmentSysPhi.C\(4,0,0,0,2,\"eta0p8\"\)
#root -l -b -q calSpinAlignmentSysPhi.C\(4,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhi.C\(3,0,0,0,2,\"eta1_eta1\"\)


#root -l -b -q subBackGround.C\(3,0,0,\"20230803\",0,1,\"eta0p4\"\)
#root -l -b -q subBackGround.C\(3,0,0,\"20230803\",0,1,\"eta0p6\"\)
#root -l -b -q subBackGround.C\(3,0,0,\"20230803\",0,1,\"eta0p8\"\)
#root -l -b -q subBackGround.C\(3,0,0,\"20230803\",0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhi.C\(3,0,0,0,1,\"eta0p4\"\)
#root -l -b -q calSpinAlignmentSysPhi.C\(3,0,0,0,1,\"eta0p6\"\)
#root -l -b -q calSpinAlignmentSysPhi.C\(3,0,0,0,1,\"eta0p8\"\)
#root -l -b -q calSpinAlignmentSysPhi.C\(3,0,0,0,1,\"eta1_eta1\"\)

#root -l -b -q subBackGround.C\(3,0,0,\"20230803\",0,2,\"eta0p4\"\)
#root -l -b -q subBackGround.C\(3,0,0,\"20230803\",0,2,\"eta0p6\"\)
#root -l -b -q subBackGround.C\(3,0,0,\"20230803\",0,2,\"eta0p8\"\)

mkdir -p figures/Phi/7GeV/pTstudy
mkdir -p figures/Phi/14GeV/pTstudy
mkdir -p figures/Phi/19GeV/pTstudy

#root -l -b -q subBackGround.C\(3,0,0,\"20230625\",0,2,\"eta1_eta1\"\)
#root -l -b -q subBackGround.C\(4,0,0,\"20230531\",0,2,\"eta1_eta1\"\)
#root -l -b -q subBackGround.C\(3,0,0,\"20230803\",0,1,\"eta1_eta1\"\)
#root -l -b -q subBackGround.C\(4,0,0,\"20230720\",0,1,\"eta1_eta1\"\)
#root -l -b -q subBackGround.C\(5,0,0,\"20240424\",0,2,\"eta1_eta1\"\)

#root -l -b -q calSpinAlignmentSysPhi_PolySys.C\(3,0,0,0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhi_PolySys.C\(3,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhi_PolySys.C\(4,0,0,0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhi_PolySys.C\(4,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhi_PolySys.C\(5,0,0,0,2,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhi_PolySys.C\(5,0,0,0,2,\"eta1_eta1\"\)

#root -l -b -q subBackGround.C\(0,0,0,\"20240511\",0,1,\"eta1_eta1\"\)
#root -l -b -q subBackGround.C\(0,0,0,\"20240511\",0,2,\"eta1_eta1\"\)

root -l -b -q calSpinAlignmentSysPhi_PolySys.C\(0,0,0,0,1,\"eta1_eta1\"\)
#root -l -b -q calSpinAlignmentSysPhi_PolySys.C\(0,0,0,0,2,\"eta1_eta1\"\)
