#include <TSystem>

// root4star -b -q EffMcPhi.C\(2,0,100000024,0,0,0\)
void EffMcPhi(const int Energy = 0, const long StartEvent = 0, const long StopEvent = 100000, const int PID = 0, const int Year = 0, const int Mode = 0, const int inputpt = 4, const int startpt = 2, const int stoppt = 5, const char* Setting = "noToF", const int etamode = 0, const int order = 2, const char *sigmay = "1000", const float v2 = 0.1, const float rho00 = 1./3., const float rerho1n1 = 0.0, const float imrho1n1 = 0.0, const float reterms = 0.0, const float imterms = 0.0, const float hrho00 = 1./3., const float hrerho1n1 = 0.0, const float himrho1n1 = 0.0, const float hreterms = 0.0, const float himterms = 0.0) // Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV, 6: 200 GeV | PID = 0: phi, 2 Kstar | Mode = 0: cent = 20-60% | Mode = 1: centrality dependence
{
  gSystem->Load("libStEffMcPhi.so");

  cout << "All libraries are loaded!!!!" << endl;
  cout << "Start to Calculate Efficiency!!" << endl;


  StEffMcPhi *mEffMcPhi = new StEffMcPhi(Energy,StartEvent,StopEvent,PID,Year,Mode,inputpt,startpt,stoppt,Setting,etamode,order);

  mEffMcPhi->setSigmay(sigmay);
  mEffMcPhi->setV2(v2);

  mEffMcPhi->setRho00(rho00);
  mEffMcPhi->setReRho1n1(rerho1n1);
  mEffMcPhi->setImRho1n1(imrho1n1);
  mEffMcPhi->setReal(reterms);
  mEffMcPhi->setImag(imterms);

  mEffMcPhi->setHRho00(hrho00);
  mEffMcPhi->setHReRho1n1(hrerho1n1);
  mEffMcPhi->setHImRho1n1(himrho1n1);
  mEffMcPhi->setHReal(hreterms);
  mEffMcPhi->setHImag(himterms);

  mEffMcPhi->Init();
  mEffMcPhi->Make();
  mEffMcPhi->Finish();
  delete mEffMcPhi;

  cout << "End of the Calculation!!" << endl;
}
