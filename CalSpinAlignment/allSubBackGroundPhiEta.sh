#!/bin/bash

mkdir -p figures/Phi/14GeV/rapiditystudy
mkdir -p figures/Phi/19GeV/rapiditystudy
mkdir -p figures/Phi/27GeV/rapiditystudy
mkdir -p ../output/AuAu27GeV/Phi
mkdir -p ../output/AuAu27GeV/Phi/Poly
mkdir -p figures/Phi/7GeV/rapiditystudy
mkdir -p ../output/AuAu7GeV/Phi
mkdir -p ../output/AuAu7GeV/Phi/Poly

# inputs ==> energy, eta, pid, year, date, random3D
for i in {0..13}
do
  #for j in {0..6}
  #do
    #root -l -b -q subBackGroundPhiEta.C\(5,${i},0,0,\"20240424\",0,2,\"eta1_eta1\"\)
    #root -l -b -q subBackGroundPhiEta.C\(0,${i},0,0,\"20240511\",0,1,\"eta1_eta1\"\)
    #root -l -b -q subBackGroundPhiEta.C\(0,${i},0,0,\"20240511\",0,2,\"eta1_eta1\"\)
    #root -l -b -q subBackGroundPhiEta_PhiPsi.C\(5,${i},${j},0,0,\"20240425\",0,2,\"eta1_eta1\"\)
    #root -l -b -q subBackGroundPhiEta_TPCTOF.C\(4,${i},0,0,\"20230531\",0,2,\"eta1_eta1\"\)
    root -l -b -q subBackGroundPhiEta_TPCOnly.C\(4,${i},0,0,\"20241003\",0,2,\"eta1p5_eta1p5\"\)
    #root -l -b -q subBackGroundPhiEta.C\(4,${i},0,0,\"20230720\",0,1,\"eta1_eta1\"\)
  #done
done
