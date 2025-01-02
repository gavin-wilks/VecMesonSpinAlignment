#include <TSystem>

// root4star -b -q EffMcPhi.C\(2,0,100000024,0,0,0\)
void EffMcPhi(const int Energy = 0, const long StartEvent = 0, const long StopEvent = 100000, const int PID = 0, const int Year = 0, const int Mode = 0, const int inputpt = 4, const int startpt = 2, const int stoppt = 3, const char* Setting = "noToF", const int etamode = 0, const int order = 2, const float rerho1n1 = 0.0, const char *sigmay = "100", const float rho00 = 1./3., const float reterms = 0.0, const float imterms = 0.0, const float imrho1n1 = 0.2) // Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV, 6: 200 GeV | PID = 0: phi, 2 Kstar | Mode = 0: cent = 20-60% | Mode = 1: centrality dependence
{
  gSystem->Load("libStEffMcPhi.so");

  cout << "All libraries are loaded!!!!" << endl;
  cout << "Start to Calculate Efficiency!!" << endl;

  StEffMcPhi *mEffMcPhi = new StEffMcPhi(Energy,StartEvent,StopEvent,PID,Year,Mode,inputpt,startpt,stoppt,Setting,etamode,order,rerho1n1,rho00,reterms,imterms,imrho1n1);
  mEffMcPhi->setSigmay(sigmay);
  //mEffMcPhi->setRho00(rho00);
  mEffMcPhi->Init();
  mEffMcPhi->Make();
  mEffMcPhi->Finish();
  delete mEffMcPhi;

  cout << "End of the Calculation!!" << endl;
}
