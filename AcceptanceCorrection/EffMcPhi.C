#include <TSystem>

// root4star -b -q EffMcPhi.C\(2,0,100000024,0,0,0\)
void EffMcPhi(const char* list ="/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Acceptance/Lists/Phi_19GeV_Tuples_pt1_subEP_withBeta.list" , const int Energy = 4, const long StartEvent = 0, const long StopEvent = 100000000, const int PID = 0, const int Mode = 1, const int etamode = 3, const int inputpt = 1, const int startpt = 0, const int stoppt = 0) // Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV, 6: 200 GeV | PID = 0: phi, 2 Kstar | Mode = 0: cent = 20-60% | Mode = 1: centrality dependence
{
  gSystem->Load("libStEffMcPhi.so");

  cout << "All libraries are loaded!!!!" << endl;
  cout << "Start to Calculate Efficiency!!" << endl;

  StEffMcPhi *mEffMcPhi = new StEffMcPhi(list,Energy,StartEvent,StopEvent,PID,Mode,etamode,inputpt,startpt,stoppt);
  mEffMcPhi->Init();
  mEffMcPhi->Make();
  mEffMcPhi->Finish();
  delete mEffMcPhi;

  cout << "End of the Calculation!!" << endl;
}
