#!/bin/sh

# void calculateFCent(const int energy = 4, const int pid = 0, int ipt = 0, bool doall = true, bool isBesI = false, bool random3D = false, int mode = 0, int etamode = 4) {
#root -l -b -q calculateFCent_EPRes_EffAcc.C\(4,0,0,1,0,0,0,0,\"EffAcc_EPSmear\"\)
#root -l -b -q calculateFCent_EPRes_EffAcc.C\(4,0,1,1,0,0,0,0,\"EffAcc_EPSmear\"\)
#root -l -b -q calculateFCent_EPRes_EffAcc.C\(4,0,2,1,0,0,0,0,\"EffAcc_EPSmear\"\)

#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(3,0,0,1,0,0,0,0,\"EffAcc_EPSmear\",${order},0\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(3,0,1,1,0,0,0,0,\"EffAcc_EPSmear\",${order},0\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(3,0,2,1,0,0,0,0,\"EffAcc_EPSmear\",${order},0\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(3,0,0,1,0,0,0,0,\"EffAcc_EPSmear\",1,-1,0,1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(3,0,1,1,0,0,0,0,\"EffAcc_EPSmear\",1,-1,0,1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(3,0,2,1,0,0,0,0,\"EffAcc_EPSmear\",1,-1,0,1\)
#
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(4,0,0,1,0,0,0,0,\"EffAcc_EPSmear\",1,-1,0,1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(4,0,1,1,0,0,0,0,\"EffAcc_EPSmear\",1,-1,0,1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(4,0,2,1,0,0,0,0,\"EffAcc_EPSmear\",1,-1,0,1\)
#
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(3,0,0,1,0,0,0,0,\"EffAcc_EPSmear\",2,-1,0,1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(3,0,1,1,0,0,0,0,\"EffAcc_EPSmear\",2,-1,0,1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(3,0,2,1,0,0,0,0,\"EffAcc_EPSmear\",2,-1,0,1\)

#root -l -b -q calculateFCent_EPRes_EffAcc_Poly_PhiPsi.C\(4,0,0,1,0,0,0,0,\"EffAcc_EPSmear\",2,-2,0,1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly_PhiPsi.C\(4,0,1,1,0,0,0,0,\"EffAcc_EPSmear\",2,-2,0,1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly_PhiPsi.C\(4,0,2,1,0,0,0,0,\"EffAcc_EPSmear\",2,-2,0,1\)

#root -l -b -q calculateFCent_EPRes_EffAcc_Poly_PhiPsi.C\(4,0,0,1,0,0,0,0,\"EffAcc_EPSmear\",2,-3,1,1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly_PhiPsi.C\(4,0,1,1,0,0,0,0,\"EffAcc_EPSmear\",2,-3,1,1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly_PhiPsi.C\(4,0,2,1,0,0,0,0,\"EffAcc_EPSmear\",2,-3,1,1\)
root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(0,0,0,1,0,0,0,0,\"\",1,-10,0,1\)
root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(0,0,1,1,0,0,0,0,\"\",1,-10,0,1\)
root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(0,0,2,1,0,0,0,0,\"\",1,-10,0,1\)
root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(0,0,0,1,0,0,0,0,\"\",2,-10,0,1\)
root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(0,0,1,1,0,0,0,0,\"\",2,-10,0,1\)
root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(0,0,2,1,0,0,0,0,\"\",2,-10,0,1\)


#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(3,0,0,1,0,0,0,0,\"EffAcc_EPSmear\",${order},1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(3,0,1,1,0,0,0,0,\"EffAcc_EPSmear\",${order},1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(3,0,2,1,0,0,0,0,\"EffAcc_EPSmear\",${order},1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(4,0,0,1,0,0,0,0,\"EffAcc_EPSmear\",${order},1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(4,0,1,1,0,0,0,0,\"EffAcc_EPSmear\",${order},1\)
#root -l -b -q calculateFCent_EPRes_EffAcc_Poly.C\(4,0,2,1,0,0,0,0,\"EffAcc_EPSmear\",${order},1\)

#root -l -b -q calculateFCent_EPRes.C\(4,0,0,1,0,0,0,0,\"sub\"\)
#root -l -b -q calculateFCent_EPRes.C\(4,0,1,1,0,0,0,0,\"sub\"\)
#root -l -b -q calculateFCent_EPRes.C\(4,0,2,1,0,0,0,0,\"sub\"\)


