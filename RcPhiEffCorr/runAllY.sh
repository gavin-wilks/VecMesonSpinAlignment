#!/bin/bash


#void EffMcPhi(const int Energy = 4, const long StartEvent = 0, const long StopEvent = 1000000000, const int PID = 0, const int Year = 0, const int Mode = 0, const int inputpt = 0, const int startpt = 2, const int stoppt = 3, const char* Setting = "noToF", const int etamode = 0, const int order = 2) // Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV, 6: 200 GeV | PID = 0: phi, 2 Kstar | Mode = 0: cent = 20-60% | Mode = 1: centrality dependence

#starver dev

root4star -l -b -q EffMcPhi.C\(0,0,10000000000,0,0,2,4,0,1,\"noToF\",0,1,0.0,\"1000\",0.3333333,0.0,0.0,0.0\)
mv Eff_7GeV_SingleParticle_noToF_Mode2_EtaMode0_PtBins0_1.root order1/.
root4star -l -b -q EffMcPhi.C\(0,0,10000000000,0,0,2,4,0,1,\"noToF\",0,2,0.0,\"1000\",0.3333333,0.0,0.0,0.0\)
mv Eff_7GeV_SingleParticle_noToF_Mode2_EtaMode0_PtBins0_1.root order2/.
#root4star -l -b -q EffMcPhi.C\(0,0,10000000000,0,0,1,4,0,2,\"noToF\",0,2,0.0,\"1000\",0.3333333,0.0,0.0,0.0\)
#root4star -l -b -q EffMcPhi.C\(0,0,10000000000,0,0,2,4,0,1,\"noToF\",0,2,0.0,\"1000\",0.3333333,0.0,0.0,0.0\)
#root4star -l -b -q EffMcPhi.C\(4,0,1000000000,0,0,0,1,4,4,\"noToF\",0,2\)
#root4star -l -b -q EffMcPhi.C\(4,0,1000000000,0,0,0,2,5,5,\"noToF\",0,2\)
