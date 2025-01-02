#!/bin/bash

for order in 1 2 
do 

  #root -l -b -q calculateFPt_EPRes_EffAcc_14GeV.C\(3,0,1,0,0,0,${order}\)

  root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_14GeV.C\(3,0,1,0,0,1,0,2,0,0,\"\",${order}\)
  root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_14GeV.C\(3,0,1,0,0,1,0,2,0,1,\"\",${order}\)
  root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_14GeV.C\(3,0,1,0,0,1,0,2,0,2,\"\",${order}\)
  root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_14GeV.C\(3,0,1,0,0,1,0,2,1,0,\"\",${order}\)
  root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_14GeV.C\(3,0,1,0,0,1,0,2,1,1,\"\",${order}\)
  root -l -b -q calculateFY_BinbyBin_EPRes_EffAcc_14GeV.C\(3,0,1,0,0,1,0,2,1,2,\"\",${order}\)
  
  root -l -b -q calculateFCent_EPRes_EffAcc_14GeV.C\(3,0,0,1,0,0,0,0,\"\",${order}\)
  root -l -b -q calculateFCent_EPRes_EffAcc_14GeV.C\(3,0,1,1,0,0,0,0,\"\",${order}\)
  root -l -b -q calculateFCent_EPRes_EffAcc_14GeV.C\(3,0,2,1,0,0,0,0,\"\",${order}\)


done
