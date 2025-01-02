#!/bin/sh

# void calculateFPtY_BinbyBin(const int energy = 4, const int pid = 0, bool doall = true, bool isBesI = false, bool random3D = false, int mode = 2, int etamode = 0, int iy = 3) {

root -l -b -q calculateFPtY_BinbyBin_EPRes.C\(4,0,1,0,0,2,0,1,\"sub\"\)
root -l -b -q calculateFPtY_BinbyBin_EPRes.C\(4,0,1,0,0,2,0,2,\"sub\"\)
root -l -b -q calculateFPtY_BinbyBin_EPRes.C\(4,0,1,0,0,2,0,3,\"sub\"\)
root -l -b -q calculateFPtY_BinbyBin_EPRes.C\(4,0,1,0,0,2,0,4,\"sub\"\)
root -l -b -q calculateFPtY_BinbyBin_EPRes.C\(4,0,1,0,0,2,0,5,\"sub\"\)

