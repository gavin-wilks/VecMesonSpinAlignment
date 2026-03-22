#!/bin/bash

#void EffMcPhi(const int Energy = 4, const long StartEvent = 0, const long StopEvent = 1000000000, const int PID = 0, const int Year = 0, const int Mode = 0, const int inputpt = 4, const int startpt = 2, const int stoppt = 5, const char* Setting = "noToF", const int etamode = 0, const int order = 2, const int iter = 0, const int study = 0, const float rerho1n1 = 0.0, const char *sigmay = "1000", const float rho00 = 1./3., const float reterms = 0.0, const float imterms = 0.0, const float imrho1n1 = 0.0, const float rhohelicity = 1./3., const float ptfixed = 1, const float yfixed = 0, int bincos = 10, int binphi = 10, int tofflag = 0, int singlekaonflag = 0) 


#void EffMcPhi(const int Energy = 4, const long StartEvent = 0, const long StopEvent = 1000000000, const int PID = 0, const int Year = 0, const int Mode = 0, const int inputpt = 4, const int startpt = 2, const int stoppt = 5, const char* Setting = "noToF", const int etamode = 0, const int order = 2, const int iter = 0, const int study = 0, const int method = 2, const float rerho1n1 = 0.0, const char *sigmay = "1000", const float rho00 = 1./3., const float reterms = 0.0, const float imterms = 0.0, const float imrho1n1 = 0.0, const float rhohelicity = 1./3., const float ptfixed = 1, const float yfixed = 0, int bincos = 10, int binphi = 10, int tofflag = 0, int singlekaonflag = 0) 

root -l -b -q EffMcPhi.C\(4,0,10000000,0,1,0,5,10,200,\"noToF\",0,1,0,0,2,0.01,\"1000\",0.01,0.0,0.0,0.0,0.3333333,100.0,100.0,50,10,0,1\)
#root4star -l -b -q EffMcPhi.C\(4,0,100000,0,1,0,5,2,0,\"noToF\",0,1,0,2,2,0.0,\"1000\",0.3333333,0.0,0.0,0.0,0.3333333,100.0,100.0,10,10,0,1\)
#root4star -l -b -q EffMcPhi.C\(4,0,100000,0,1,0,5,2,0,\"noToF\",0,1,0,1,1,0.0,\"1000\",0.3333333,0.0,0.0,0.0,0.3333333,100.0,100.0,10,10,0,1\)
#root4star -l -b -q EffMcPhi.C\(4,0,100000,0,1,0,5,2,0,\"noToF\",0,1,0,2,1,0.0,\"1000\",0.3333333,0.0,0.0,0.0,0.3333333,100.0,100.0,10,10,0,1\)

#root4star -l -b -q EffMcPhi.C\(4,0,1000000000,0,1,0,5,3,1000,\"noToF\",0,2,0.0,\"1000\",0.3333333,0.0,0.0,0.0,0.3333333,100.0,100.0,10,10,0,1\)
#root4star -l -b -q EffMcPhi.C\(4,0,1000000000,0,1,0,5,4,1000,\"noToF\",0,2,0.0,\"1000\",0.3333333,0.0,0.0,0.0,0.3333333,100.0,100.0,10,10,0,1\)
#root4star -l -b -q EffMcPhi.C\(4,0,1000000000,0,1,0,5,5,1000,\"noToF\",0,2,0.0,\"1000\",0.3333333,0.0,0.0,0.0,0.3333333,100.0,100.0,10,10,0,1\)
