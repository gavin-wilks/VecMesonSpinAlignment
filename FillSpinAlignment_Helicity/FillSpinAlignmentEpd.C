#include <TSystem>
#include "TStopwatch.h"

void FillSpinAlignmentEpd(const char* list = "../FileLists/19GeV_2019/Phi_SE_FirstOrder_eta1_test.list", const char *jobId = "testing", const Int_t energy = 3, const Int_t X_flag = 0, const Int_t mode = 0, const Int_t etamode = 0)
{
  // mBeamEnergy[NumBeamEnergy] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  // X_flag: 0 for Same Event, 1 for Mixed Event
  // List: different number for different TTree list
  // mode: 0 for phi meson, 1 for K*, 2 for K0S

  TStopwatch *stopWatch = new TStopwatch();
  stopWatch->Start();

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StRefMultCorr");
  //gSystem->Load("StPicoEvent");
  //gSystem->Load("StPicoDstMaker");
  gSystem->Load("StMesonEvent");
  gSystem->Load("StVecMesonAna");
  //gSystem->Load("StRunIdEventsDb");

  cout << "All libraries are loaded!!!!" << endl;

  cout << "Start to Read Trees!" << endl;

  StVecMesonAnaEpd *mVecMesonAna = new StVecMesonAnaEpd(list,jobId,energy,X_flag,mode,etamode);
  mVecMesonAna->Init();
  mVecMesonAna->Make();
  mVecMesonAna->Finish();

  stopWatch->Stop();
  stopWatch->Print();

  cout << "End of the Calculation!!" << endl;
}
