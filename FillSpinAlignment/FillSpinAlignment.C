#include <TSystem>
#include "TStopwatch.h"

void FillSpinAlignment(const char* list = "../FileLists/19GeV_2019/Phi_SE_Forest.list", const char *jobId = "Final", const Int_t energy = 4, const Int_t X_flag = 0, const Int_t mode = 0)
{
  // mBeamEnergy[NumBeamEnergy] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  // X_flag: 0 for Same Event, 1 for Mixed Event
  // List: different number for different TTree list
  // mode: 0 for phi meson, 1 for K*, 2 for K0S

  TStopwatch *stopWatch = new TStopwatch();
  stopWatch->Start();

  gSystem->Load("StRefMultCorr");
  // gSystem->Load("StPicoDstMaker");
  gSystem->Load("StMesonEvent");
  gSystem->Load("StVecMesonAna");
  gSystem->Load("StRunIdEventsDb");

  cout << "All libraries are loaded!!!!" << endl;

  cout << "Start to Read Trees!" << endl;

  StVecMesonAna *mVecMesonAna = new StVecMesonAna(list,jobId,energy,X_flag,mode);
  mVecMesonAna->Init();
  mVecMesonAna->Make();
  mVecMesonAna->Finish();

  stopWatch->Stop();
  stopWatch->Print();

  cout << "End of the Calculation!!" << endl;
}
