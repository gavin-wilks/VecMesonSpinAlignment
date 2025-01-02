#!/bin/sh

# void calculateFY_BinbyBin(const int energy = 4, const int pid = 0, bool doall = true, bool isBesI = false, bool random3D = false, int mode = 1, int etamode = 0, int ypadding = 2, int ipt = 0, int icent = 0, std::string etastring = "eta1_eta1") {
order=$1
root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,5,3,0,0,\"EffAcc_EPSmear_RapidityV2\",${order}\)
root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,5,3,0,1,\"EffAcc_EPSmear_RapidityV2\",${order}\)
root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,5,3,0,2,\"EffAcc_EPSmear_RapidityV2\",${order}\)
root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,5,3,1,0,\"EffAcc_EPSmear_RapidityV2\",${order}\)
root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,5,3,1,1,\"EffAcc_EPSmear_RapidityV2\",${order}\)
root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,5,3,1,2,\"EffAcc_EPSmear_RapidityV2\",${order}\)
#root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,5,3,0,0,\"EffAcc_EPSmear_RapidityV2\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,5,3,0,1,\"EffAcc_EPSmear_RapidityV2\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,5,3,0,2,\"EffAcc_EPSmear_RapidityV2\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,5,3,1,0,\"EffAcc_EPSmear_RapidityV2\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,5,3,1,1,\"EffAcc_EPSmear_RapidityV2\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,5,3,1,2,\"EffAcc_EPSmear_RapidityV2\"\)

#root -l -b -q calculateFY_BinbyBin_EPRes.C\(4,0,1,0,0,1,5,3,0,0,\"sub\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes.C\(4,0,1,0,0,1,5,3,0,1,\"sub\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes.C\(4,0,1,0,0,1,5,3,0,2,\"sub\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes.C\(4,0,1,0,0,1,5,3,1,0,\"sub\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes.C\(4,0,1,0,0,1,5,3,1,1,\"sub\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes.C\(4,0,1,0,0,1,5,3,1,2,\"sub\"\)
