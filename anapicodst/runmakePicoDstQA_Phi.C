//wrapper to makePicoDstQA.C
//for ACLiC mode

void runmakePicoDstQA_Phi(TString InputFileList,Int_t nevents,TString OutputFile, TString jobId, Int_t mEnergy, Int_t mGid, Int_t mPid)
{
   gROOT->Reset();
   gROOT->Macro("loadMuDst.C");
   gSystem->Load("StRefMultCorr");
   gSystem->Load("StPicoEvent");
   gSystem->Load("StPicoDstMaker");
   gSystem->Load("StVecMesonMaker");
   gSystem->Load("StEffKaon");
   gROOT->LoadMacro("makePicoDstQA_Phi.C+");
   makePicoDstQA_Phi(InputFileList,nevents,OutputFile,jobId,mEnergy,mGid,mPid); 
}
