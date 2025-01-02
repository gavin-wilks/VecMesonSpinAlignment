//wrapper to makePicoDstQA.C
//for ACLiC mode

//#include "makePicoDstQA_Phi.C"

void runmakePicoDstQA_Phi(TString InputFileList,Int_t nevents,TString OutputFile, TString jobId, Int_t mEnergy, Int_t mGid, Int_t mPid, bool phiswitch, float rho = 0.3333, float real = 0.0, float imag = 0.0, float re = 0.0, float im = 0.0, float hrho = 0.0, float hreal = 0.0, float himag = 0.0, float hre = 0.0, float him = 0.0)
{
   gROOT->Reset();
   gROOT->Macro("loadMuDst.C");
   gSystem->Load("StRefMultCorr");
   gSystem->Load("StPicoEvent");
   gSystem->Load("StPicoDstMaker");
   gSystem->Load("StVecMesonMaker");
   gSystem->Load("StEffKaon");
   gROOT->LoadMacro("makePicoDstQA_Phi.C+");
   makePicoDstQA_Phi(InputFileList,nevents,OutputFile,jobId,mEnergy,mGid,mPid,phiswitch,rho,real,imag,re,im,hrho,hreal,himag,hre,him); 
}
