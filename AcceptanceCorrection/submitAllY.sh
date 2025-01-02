#!/bin/sh
order=$1
# void calculateFY_BinbyBin(const int energy = 4, const int pid = 0, bool doall = true, bool isBesI = false, bool random3D = false, int mode = 1, int etamode = 0, int ypadding = 2, int ipt = 0, int icent = 0, std::string etastring = "eta1_eta1") {
./submitAllY_EtaMode0.sh ${order} 
./submitAllY_EtaMode3.sh ${order}
./submitAllY_EtaMode4.sh ${order}
./submitAllY_EtaMode5.sh ${order}
