#!/bin/sh

# void calculateFCent(const int energy = 4, const int pid = 0, int ipt = 0, bool doall = true, bool isBesI = false, bool random3D = false, int mode = 0, int etamode = 4) {
order=$1
#root -l -b -q calculateFCent_EPRes_EffAcc.C\(4,0,0,1,0,0,0,4,\"EffAcc_EPSmear\"\)
#root -l -b -q calculateFCent_EPRes_EffAcc.C\(4,0,1,1,0,0,0,4,\"EffAcc_EPSmear\"\)
#root -l -b -q calculateFCent_EPRes_EffAcc.C\(4,0,2,1,0,0,0,4,\"EffAcc_EPSmear\"\)
root -l -b -q calculateFCent_EPRes_EffAcc.C\(4,0,0,1,0,0,0,4,\"EffAcc_EPSmear\",${order}\)
root -l -b -q calculateFCent_EPRes_EffAcc.C\(4,0,1,1,0,0,0,4,\"EffAcc_EPSmear\",${order}\)
root -l -b -q calculateFCent_EPRes_EffAcc.C\(4,0,2,1,0,0,0,4,\"EffAcc_EPSmear\",${order}\)

#root -l -b -q calculateFCent_EPRes.C\(4,0,0,1,0,0,0,4,\"sub\"\)
#root -l -b -q calculateFCent_EPRes.C\(4,0,1,1,0,0,0,4,\"sub\"\)
#root -l -b -q calculateFCent_EPRes.C\(4,0,2,1,0,0,0,4,\"sub\"\)

