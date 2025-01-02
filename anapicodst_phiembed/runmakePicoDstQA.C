//wrapper to makePicoDstQA.C
//for ACLiC mode

void runmakePicoDstQA(TString InputFileList,Int_t nevents,TString OutputFile, TString jobId, Int_t mEnergy, Int_t mGid, Int_t mPid)
{
   gROOT->Reset();
   gROOT->Macro("loadMuDst.C");
   gSystem->Load("StRefMultCorr");
   gSystem->Load("StPicoEvent");
   gSystem->Load("StPicoDstMaker");
   gSystem->Load("StVecMesonMaker");
   gSystem->Load("StEffKaon");
   gROOT->LoadMacro("makePicoDstQA.C+");
   makePicoDstQA(InputFileList,nevents,OutputFile,jobId,mEnergy,mGid,mPid,1./3.,0.0,0.0,0.0,0.0,1./3.,0.0,0.0,0.0,0.0); 
}
