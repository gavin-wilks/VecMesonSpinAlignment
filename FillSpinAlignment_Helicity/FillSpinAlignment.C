#include <TSystem>
#include "TStopwatch.h"

void FillSpinAlignment(const char* list = "submit/19GeV_2019/JOBS/list/schedF07E1BA3957BC4F0978B6B95317ABA0E_8.list", const char *jobId = "resubmit", const Int_t energy = 4, const Int_t X_flag = 1, const Int_t mode = 0, const Int_t etamode = 0, const Int_t study = 0)
{
  // mBeamEnergy[NumBeamEnergy] = {"7GeV","9GeV","11GeV","14GeV","19GeV"};
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

  StVecMesonAna *mVecMesonAna = new StVecMesonAna(list,jobId,energy,X_flag,mode,etamode,study);
  mVecMesonAna->Init();
  mVecMesonAna->Make();
  mVecMesonAna->Finish();

  stopWatch->Stop();
  stopWatch->Print();

  cout << "End of the Calculation!!" << endl;
}
