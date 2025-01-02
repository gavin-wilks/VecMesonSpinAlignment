#!/bin/sh

# void calculateFY_BinbyBin(const int energy = 4, const int pid = 0, bool doall = true, bool isBesI = false, bool random3D = false, int mode = 1, int etamode = 0, int ypadding = 2, int ipt = 0, int icent = 0, std::string etastring = "eta1_eta1") {
#void calculateFY_BinbyBin_EPRes(const int energy = 4, const int pid = 0, bool doall = false, bool isBesI = false, bool random3D = false, int mode = 1, int etamode = 3, int ypadding = 5, int ipt = 1, int icent = 1, std::string ep = "noEPsmear")

#order=$1

#root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,0,2,0,0,\"EffAcc_EPSmear_RapidityV2\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,0,2,0,1,\"EffAcc_EPSmear_RapidityV2\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,0,2,0,2,\"EffAcc_EPSmear_RapidityV2\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,0,2,1,0,\"EffAcc_EPSmear_RapidityV2\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,0,2,1,1,\"EffAcc_EPSmear_RapidityV2\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc.C\(4,0,1,0,0,1,0,2,1,2,\"EffAcc_EPSmear_RapidityV2\"\)

#for energy in 3 4 
#do
#  for order in 1 2 
#  do 
#    root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_19GeV.C\(${energy},0,1,0,0,1,0,2,0,0,\"EffAcc_EPSmear_RapidityV2\",${order}\)
#    root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_19GeV.C\(${energy},0,1,0,0,1,0,2,0,1,\"EffAcc_EPSmear_RapidityV2\",${order}\)
#    root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_19GeV.C\(${energy},0,1,0,0,1,0,2,0,2,\"EffAcc_EPSmear_RapidityV2\",${order}\)
#    root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_19GeV.C\(${energy},0,1,0,0,1,0,2,1,0,\"EffAcc_EPSmear_RapidityV2\",${order}\)
#    root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_19GeV.C\(${energy},0,1,0,0,1,0,2,1,1,\"EffAcc_EPSmear_RapidityV2\",${order}\)
#    root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_19GeV.C\(${energy},0,1,0,0,1,0,2,1,2,\"EffAcc_EPSmear_RapidityV2\",${order}\)
#  done
#done

#3d random
for energy in 0
do
  for order in 1 2
  do 
    for poly in 0 1 
    do
      #root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi.C\(${energy},0,1,0,0,1,0,2,0,0,${poly},\"EffAcc_EPSmear_RapidityV2\",${order},-2,0,1\)
      #root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi.C\(${energy},0,1,0,0,1,0,2,0,1,${poly},\"EffAcc_EPSmear_RapidityV2\",${order},-2,0,1\)
      #root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi.C\(${energy},0,1,0,0,1,0,2,0,2,${poly},\"EffAcc_EPSmear_RapidityV2\",${order},-2,0,1\)
      #root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi.C\(${energy},0,1,0,0,1,0,2,1,0,${poly},\"EffAcc_EPSmear_RapidityV2\",${order},-2,0,1\)
      #root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi.C\(${energy},0,1,0,0,1,0,2,1,1,${poly},\"EffAcc_EPSmear_RapidityV2\",${order},-2,0,1\)
      #root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi.C\(${energy},0,1,0,0,1,0,2,1,2,${poly},\"EffAcc_EPSmear_RapidityV2\",${order},-2,0,1\)

#void calculateFY_BinbyBin_EPRes_EffAcc_Poly(const int energy = 4, const int pid = 0, bool doall = false, bool isBesI = false, bool random3D = false, int mode = 1, int etamode = 0, int ypadding = 2, int ipt = 1, int icent = 1, int i_poly = 0, std::string ep = "sub", int order = 1, int yspectra = 0, int EP = 1, int v2 = 0) {
      root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi.C\(${energy},0,1,0,0,1,0,2,0,0,${poly},\"\",${order},-10,1,1\)
      root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi.C\(${energy},0,1,0,0,1,0,2,0,1,${poly},\"\",${order},-10,1,1\)
      root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi.C\(${energy},0,1,0,0,1,0,2,0,2,${poly},\"\",${order},-10,1,1\)
      root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi.C\(${energy},0,1,0,0,1,0,2,1,0,${poly},\"\",${order},-10,1,1\)
      root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi.C\(${energy},0,1,0,0,1,0,2,1,1,${poly},\"\",${order},-10,1,1\)
      root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi.C\(${energy},0,1,0,0,1,0,2,1,2,${poly},\"\",${order},-10,1,1\)
    done
  done
done

#root -l -b -q calculateFY_BinbyBin_EPRes.C\(4,0,1,0,0,1,0,2,0,0,\"sub\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes.C\(4,0,1,0,0,1,0,2,0,1,\"sub\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes.C\(4,0,1,0,0,1,0,2,0,2,\"sub\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes.C\(4,0,1,0,0,1,0,2,1,0,\"sub\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes.C\(4,0,1,0,0,1,0,2,1,1,\"sub\"\)
#root -l -b -q calculateFY_BinbyBin_EPRes.C\(4,0,1,0,0,1,0,2,1,2,\"sub\"\)

