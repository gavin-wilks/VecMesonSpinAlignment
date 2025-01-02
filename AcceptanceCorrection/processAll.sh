#!/bin/bash

#void EffMcPhi(const char* list ="/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt1_subEP_withBeta.list" , const int Energy = 4, const long StartEvent = 0, const long StopEvent = 100000000, const int PID = 0, const int Mode = 1, const int etamode = 3, const int inputpt = 1, const int startpt = 0, const int stoppt = 0) // Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV, 6: 200 GeV | PID = 0: phi, 2 Kstar | Mode = 0: cent = 20-60% | Mode = 1: centrality dependence

root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt1_subEP_withBeta.list\",4,0,100000000,0,1,0,1,0,0\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt4_subEP_withBeta.list\",4,0,100000000,0,1,0,4,1,1\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt1_subEP_withBeta.list\",4,0,100000000,0,1,3,1,0,0\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt4_subEP_withBeta.list\",4,0,100000000,0,1,3,4,1,1\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt1_subEP_withBeta.list\",4,0,100000000,0,1,4,1,0,0\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt4_subEP_withBeta.list\",4,0,100000000,0,1,4,4,1,1\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt1_subEP_withBeta.list\",4,0,100000000,0,1,5,1,0,0\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt4_subEP_withBeta.list\",4,0,100000000,0,1,5,4,1,1\)

root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt1_subEP_withBeta.list\",4,0,100000000,0,0,0,1,0,1\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt4_subEP_withBeta.list\",4,0,100000000,0,0,0,4,2,2\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt1_subEP_withBeta.list\",4,0,100000000,0,0,3,1,0,1\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt4_subEP_withBeta.list\",4,0,100000000,0,0,3,4,2,2\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt1_subEP_withBeta.list\",4,0,100000000,0,0,4,1,0,1\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt4_subEP_withBeta.list\",4,0,100000000,0,0,4,4,2,2\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt1_subEP_withBeta.list\",4,0,100000000,0,0,5,1,0,1\)
root4star -l -b -q EffMcPhi.C\(\"/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt4_subEP_withBeta.list\",4,0,100000000,0,0,5,4,2,2\)
