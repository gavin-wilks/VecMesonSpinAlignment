//wrapper to makePicoDstQA.C
//for ACLiC mode

void runmakePicoDstQA(TString InputFileList,Int_t nevents,TString OutputFile, TString jobId, Int_t mEnergy, Int_t mGid)
{
   gROOT->Reset();
   gROOT->Macro("loadMuDst.C");
   gSystem->Load("StRefMultCorr");
   gSystem->Load("StPicoEvent");
   gSystem->Load("StPicoDstMaker");
   gSystem->Load("StVecMesonMaker");
   gSystem->Load("StEffKaon");
   gROOT->LoadMacro("makePicoDstQA.C+");
   makePicoDstQA(InputFileList,nevents,OutputFile,jobId,mEnergy,mGid); 
}
