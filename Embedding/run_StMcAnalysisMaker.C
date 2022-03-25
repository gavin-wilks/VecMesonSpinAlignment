#include <iostream>
#include <string>

void Load();

using namespace std;

void run_StMcAnalysisMaker(const char* file = "/star/embed/embedding/production_19GeV_2019/Kplus_200_20214203/P21ic.SL21c/2019/057/20057003/st_physics_adc_20057003_raw_4000002.geant.root", std::string outFile = "test", int energy = 4)
{
  //Check STAR Library. Please set SL_version to the original star library used
  // in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl

  string SL_version = "SL21c"; // 19.6 GeV
  string env_SL = getenv("STAR");

  if (env_SL.find(SL_version) == string::npos)
  {
    cout << SL_version << " " << env_SL << endl;
    cout << "Environment Star Library does not match the requested library in run_st_etree.C. Exiting..." << endl;
    return;
  }

  // Load shared libraries
  Load();

  // Create chain
  StChain* chain = new StChain;

  // I/O maker
  StIOMaker* ioMaker = new StIOMaker;
  ioMaker->SetFile(file);
  ioMaker->SetIOMode("r");
  ioMaker->SetBranch("*", 0, "0");
  //ioMaker->SetBranch("McEventBranch",0,"r");
  ioMaker->SetBranch("geantBranch", 0, "r");
  ioMaker->SetBranch("eventBranch", 0, "r");

  TString picodstfile = file;

  if(picodstfile.First("$") != -1)
  {
    picodstfile.ReplaceAll("$","");
    picodstfile = getenv(picodstfile.Data());
  }

  picodstfile.ReplaceAll(".event.root", ".picoDst.root");
  picodstfile.ReplaceAll(".geant.root", ".picoDst.root");
  cout << "Reading PicoDst file " << picodstfile << endl;
  //StPicoDstMaker* picoDstMaker = new StPicoDstMaker(0, 0, "", picodstfile.Data(), "", 100000, "PicoDst");

  StPicoDstMaker *picoDstMaker = new StPicoDstMaker(2,picodstfile.Data(),"picoDst");
  StMcEventMaker *mcEventMaker = new StMcEventMaker();
  mcEventMaker->doPrintEventInfo = false;
  mcEventMaker->doPrintMemoryInfo = false;

  StAssociationMaker* assoc = new StAssociationMaker;
  assoc->useInTracker();
  assoc->SetDebug();

  //.. see example in CVS: StRoot/macros/mudst/exampleEmc.C
  // Need St_db_Maker for Emc calibration
  // St_db_Maker* dbMk = new St_db_Maker("StarDb", "MySQL:StarDb");
  //dbMk->SetMaxEntryTime(20100301,0);
  //dbMk->SetDateTime(20080101,000001);
  //int refMultCorrLoad = gSystem->Load("StRefMultCorr");
  //StRefMultCorr* grefmultCorrUtil = NULL;

  //if (refMultCorrLoad == -1)
  //{
  //  cout << "StRefMultCorr library is not available" << endl;
  //}
  //else
  //{
  //  grefmultCorrUtil = CentralityMaker::instance()->getRefMultCorr();
    // grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr();
    // grefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
    // grefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14.txt");
    //
    // for (int i = 0; i < 6; ++i)
    // {
    //   cout << i << " " << grefmultCorrUtil->get(i, 0) << endl;
    // }
  //}

  // Monte Carlo event maker
  StMcAnalysisMaker* analysis = new StMcAnalysisMaker;
  analysis->setOutFileName(outFile);
  analysis->setEnergy(energy);
  //analysis->setRefMultCorr(grefmultCorrUtil);

  // Initialize chain
  chain->Init();
  chain->EventLoop(1e7);
  // chain->EventLoop(100);
  chain->Finish();

  //delete chain;
}

void Load()
{
   
   gROOT->Macro("LoadLogger.C");
   gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
   loadSharedLibraries();
   gSystem->Load("StChain");
   gSystem->Load("StPicoEvent");
   gSystem->Load("StPicoDstMaker");
   gSystem->Load("StIOMaker");
   gSystem->Load("StTreeMaker");
   gSystem->Load("StMcEvent");
   gSystem->Load("StMcEventMaker");
   gSystem->Load("StAssociationMaker");
   gSystem->Load("StDbLib.so");
   gSystem->Load("StDbBroker.so");
   gSystem->Load("libglobal_Tables.so");
   gSystem->Load("St_db_Maker.so");
   gSystem->Load("StDetectorDbMaker");
   gSystem->Load("StTpcDb");
   gSystem->Load("StDbUtilities");
   gSystem->Load("StMcEvent");
   gSystem->Load("StMcEventMaker");
   gSystem->Load("StAssociationMaker");
   gSystem->Load("StEmcRawMaker");
   gSystem->Load("StEmcADCtoEMaker");
   gSystem->Load("StPreEclMaker");
   gSystem->Load("StEpcMaker");
   gSystem->Load("StEmcSimulatorMaker");
   gSystem->Load("StEmcUtil");
   gSystem->Load("StEEmcUtil");
   gSystem->Load("StEEmcDbMaker");
   gSystem->Load("StEmcTriggerMaker");
   gSystem->Load("StDaqLib");
   gSystem->Load("StMcAnalysisMaker");
}
