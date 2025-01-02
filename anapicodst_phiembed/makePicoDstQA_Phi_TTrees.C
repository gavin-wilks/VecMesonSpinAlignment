//Xianglei Zhu 
//Skeleton embedding PicoDst analysis macro with StPicoDstMaker 
//Run it with the wrapper in ACLIC mode, CINT mode for debug ONLY

#ifndef __CINT__
#include "TRandom.h" 
#include "TROOT.h"
#include "TMath.h"
#include "TSystem.h"
#include <iostream>
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH3.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeHelper.h"
#include "TNtuple.h"
#include "TDatime.h"
#include "StarRoot/TUnixTime.h"
#include "StChain.h"
#include "StMessMgr.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoMcTrack.h"
#include "StPicoEvent/StPicoMcVertex.h"
#include "StPicoEvent/StPicoArrays.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StBTofHeader.h"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "StRoot/StVecMesonMaker/StVecMesonCut.h"
#include "StRoot/StVecMesonMaker/StUtility.h"
#include "StRoot/StEffKaon/StEffHistManger.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StRoot/Utility/functions.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include <map>
#include "TGraphAsymmErrors.h"
//#include <pair>
#include <algorithm>
typedef std::map<std::string,TH3F*> TH3FMap;
typedef std::map<std::string,TH1D*> TH1DMap;
typedef std::map<std::string,TF1*> TF1Map;
TH1DMap h_TofKplus;
TH1DMap h_TofKminus;
TH1D *h_FrameEta_ToF[2];
TH1D *h_FramePhi_ToF[2];

TF1Map f_TofKplus;
TF1Map f_TofKminus;

TF1Map f_nsig_PID;
TF1Map f_m2_PID;
TH3FMap h_m2_PID;
TH3FMap h_nsig_PID;

TF1 *f_mRhoPt[6];
TF1 *f_mRhoPt_Helicity[6];
TF2 *f_mRhoPt_2D[6];
TF2 *f_mRhoPt_Helicity_2D[6];
TF2* f_mRhoY_Helicity_2D[10];

float mMaxData[6];
float mMaxHelicity[6];
float mMax2D[6];
float mMaxHelicity2D[6];
float mMaxHelicity2DY[10];

//void readTofEff(int energy);
//void readTofEffFit(int energy);
//void findHist_ToF(TLorentzVector const& lKaon, int iParticleIndex, float Psi2, int& EtaBin, int& PhiBin);
//bool tofReconstructed(int iParticleIndex, int cent, float Psi2, TLorentzVector const& lKaon);

#endif

#include "TLorentzVector.h"
#include "TF1.h"
#include "TF2.h"

bool Sampling(TF1 *f_rhoPhy, float CosThetaStar, float wMax)
{
  return !(gRandom->Rndm() > f_rhoPhy->Eval(CosThetaStar)/wMax);
}

bool Sampling2D(TF2 *f_rhoPhy, float ThetaStar, float PhiPrime, float wMax)
{
  return !(gRandom->Rndm() > f_rhoPhy->Eval(ThetaStar,PhiPrime)/wMax);
}

void readm2PID_hist()
{
  string InPut = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/PID/m2_PID_hist_19GeV.root");
  TFile *File = TFile::Open(InPut.c_str());
  std::cout << "m2 PID file: " << InPut << endl; 

  int icent = 9;
  for(int ic = 0; ic < 2; ic++)
  {
    string histname = Form("h_m2_PID_cent%d_charge%d",icent,ic);
    h_m2_PID[histname] = (TH3F*)File->Get(histname.c_str());
  } 
}

void findm2HistPID(TLorentzVector const& lKaon, int icharge, int &GlobalBin)
{
  double pt  = lKaon.Pt();
  double eta = lKaon.PseudoRapidity();
  double phi = lKaon.Phi();
  while(phi < -TMath::Pi()) phi += 2.0*TMath::Pi();
  while(phi >  TMath::Pi()) phi -= 2.0*TMath::Pi();

  int icent = 9;
  string histname = Form("h_m2_PID_cent%d_charge%d",icent,icharge);
  GlobalBin = h_m2_PID[histname]->FindBin(pt,eta,phi);
}

bool passhist_m2_PID(int icharge, int icent, TLorentzVector const& lKaon, int GlobalBin)
{
  TH3F *hist = NULL;

  string KEY; 
  KEY = Form("h_m2_PID_cent%d_charge%d",icent,icharge);
  hist = h_m2_PID[KEY]; // only 20-60%

  bool pass_m2 = false;
  double prob_m2 = 0.0;
  prob_m2 = hist->GetBinContent(GlobalBin);
  pass_m2 = gRandom->Rndm() < prob_m2;
   
  return pass_m2;
}

void readnsigPID_hist()
{
  string InPut = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/PID/nsig_PID_hist_19GeV.root");
  TFile *File = TFile::Open(InPut.c_str());
  std::cout << "m2 PID file: " << InPut << endl; 

  int icent = 9;
  for(int ic = 0; ic < 2; ic++)
  {
    string histname = Form("h_nsig_PID_cent%d_charge%d",icent,ic);
    h_nsig_PID[histname] = (TH3F*)File->Get(histname.c_str());
  } 
}

void findnsigHistPID(TLorentzVector const& lKaon, int icharge, int &GlobalBin)
{
  double pt  = lKaon.Pt();
  double eta = lKaon.PseudoRapidity();
  double phi = lKaon.Phi();
  while(phi < -TMath::Pi()) phi += 2.0*TMath::Pi();
  while(phi >  TMath::Pi()) phi -= 2.0*TMath::Pi();

  int icent = 9;
  string histname = Form("h_nsig_PID_cent%d_charge%d",icent,icharge);
  GlobalBin = h_nsig_PID[histname]->FindBin(pt,eta,phi);
}

bool passhist_nsig_PID(int icharge, int icent, TLorentzVector const& lKaon, int GlobalBin)
{
  TH3F *hist = NULL;

  string KEY; 
  KEY = Form("h_nsig_PID_cent%d_charge%d",icent,icharge);
  hist = h_nsig_PID[KEY]; // only 20-60%

  bool pass_nsig = false;
  double prob_nsig = 0.0;
  prob_nsig = hist->GetBinContent(GlobalBin);
  pass_nsig = gRandom->Rndm() < prob_nsig;
   
  return pass_nsig;
}


void readm2PID()
{
  string InPut = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/PID/m2_PID_funcparameters_19GeV.root");
  TFile *File = TFile::Open(InPut.c_str());
  std::cout << "m2 PID file: " << InPut << endl; 

  int icent = 9;
  for(int ic = 0; ic < 2; ic++)
  {
    for(int ieta = 0; ieta < 10; ieta++)
    {
      for(int iphi = 0; iphi < 12; iphi++)
      {
        TH1F *parameterHist; 
        string funcname = Form("m2_PID_cent%d_charge%d_eta%d_phi%d",icent,ic,ieta,iphi);
        f_m2_PID[funcname] = new TF1(funcname.c_str(),m2_PID_tanh,0.1,5.0,4); 

        string histname = Form("m2_parameters_cent%d_charge%d_eta%d_phi%d",icent,ic,ieta,iphi);
        parameterHist = (TH1F*)File->Get(histname.c_str());

        for(int i_par = 0; i_par < 4; i_par++)
        {
          f_m2_PID[funcname]->SetParameter(i_par,parameterHist->GetBinContent(1+i_par));
        }
      }
    }
  } 
}

void readnsigPID()
{
  string InPut = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/PID/nsig_PID_funcparameters_19GeV.root");
  TFile *File = TFile::Open(InPut.c_str());
  std::cout << "nsig PID file: " << InPut << endl; 

  int icent = 9;
  for(int ic = 0; ic < 2; ic++)
  {
    for(int ieta = 0; ieta < 10; ieta++)
    {
      for(int iphi = 0; iphi < 12; iphi++)
      {
        TH1F *parameterHist; 
        string funcname = Form("nsig_PID_cent%d_charge%d_eta%d_phi%d",icent,ic,ieta,iphi);
        f_nsig_PID[funcname] = new TF1(funcname.c_str(),nsig_PID_tanh,0.1,5.0,3); 

        string histname = Form("nsig_parameters_cent%d_charge%d_eta%d_phi%d",icent,ic,ieta,iphi);
        parameterHist = (TH1F*)File->Get(histname.c_str());

        for(int i_par = 0; i_par < 3; i_par++)
        {
          f_nsig_PID[funcname]->SetParameter(i_par,parameterHist->GetBinContent(1+i_par));
        }
      }
    }
  } 
}

void findFuncPID(TLorentzVector const& lKaon, int icharge, int &EtaBin, int &PhiBin)
{
  double eta = lKaon.PseudoRapidity();
  double phi = lKaon.Phi();
  while(phi < -TMath::Pi()) phi += 2.0*TMath::Pi();
  while(phi >  TMath::Pi()) phi -= 2.0*TMath::Pi();

  //cout << "Eta = " << eta << endl;

  for(int ieta = 0; ieta < 10; ieta++)
  {
    if(eta >= double(ieta-5)/5.0 && eta < double(ieta-4)/5.0) 
    {
      EtaBin = ieta;      
      break;
    }
  }

  for(int iphi = 0; iphi < 12; iphi++)
  {
    if(phi >= TMath::Pi()*(double(iphi)-6.)/6. && phi < TMath::Pi()*(double(iphi)-5.)/6.) 
    {
      PhiBin = iphi;      
      break;
    }
  }
}

bool pass_m2_PID(int icharge, int icent, TLorentzVector const& lKaon, int EtaBin, int PhiBin)
{
  TF1 *f_m2 = NULL;

  string KEY; 
  KEY = Form("m2_PID_cent%d_charge%d_eta%d_phi%d",icent,icharge,EtaBin,PhiBin);
  //cout << KEY << endl;
  f_m2 = f_m2_PID[KEY]; // only 20-60%
  //f_m2->Print();

  double pt = lKaon.Perp();

  bool pass_m2 = false;
  double prob_m2 = 0.0;
  prob_m2 = f_m2->Eval(pt);
  pass_m2 = gRandom->Rndm() < prob_m2;
   
  return pass_m2;
}

bool pass_nsig_PID(int icharge, int icent, TLorentzVector const& lKaon, int EtaBin, int PhiBin)
{
  TF1 *f_nsig = NULL;

  string KEY; 
  KEY = Form("nsig_PID_cent%d_charge%d_eta%d_phi%d",icent,icharge,EtaBin,PhiBin);
  //cout << KEY << endl;
  f_nsig = f_nsig_PID[KEY]; // only 20-60%
  //f_nsig->Print();

  double pt = lKaon.Perp();

  bool pass_nsig = false;
  double prob_nsig = 0.0;
  prob_nsig = f_nsig->Eval(pt);
  pass_nsig = gRandom->Rndm() < prob_nsig;
   
  return pass_nsig;
}

//TF1* readm2PID(int icent, int icharge, float ieta, float iphi)
//{
//  TH1F *parameterHist; 
//  string funcname = Form("m2_PID_cent%d_charge%d_eta%d_phi%d",icent,icharge,ieta,iphi);
//  TF1 *m2_PID_pT = new TF1(funcname.c_str(),m2_PID_tanh,0.1,5.0,4); 
//
//  string InPut = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/PID/m2_PID_funcparameters_19GeV.root");
//  TFile *File = TFile::Open(InPut.c_str());
//  std::cout << "m2 PID file: " << InPut << endl; 
//
//  string histname = Form("m2_parameters_cent%d_charge%d_eta%d_phi%d",icent,icharge,ieta,iphi);
//  parameterHist = (TH1F*)File->Get(histname.c_str());
//
//  for(int i_par = 0; i_par < 5; i_par++)
//  {
//    m2_PID_pT->SetParameter(parameterHist->GetBinContent(1+i_par));
//  }
// 
//  return m2_PID_pT;
//}

//TF1* readnsigPID(int icent, int icharge, float ieta, float iphi)
//{
//  string InPut = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/PID/nsig_PID_funcparameters_19GeV.root");
//  TFile *File = TFile::Open(InPut.c_str());
//  std::cout << "nsig PID file: " << InPut << endl; 
//
//  int icent = 9;
//  for(int 
//  TH1F *parameterHist; 
//  string funcname = Form("nsig_PID_cent%d_charge%d_eta%d_phi%d",icent,icharge,ieta,iphi);
//  f_nsig_PID[funcname] = new TF1(funcname.c_str(),nsig_PID_tanh,0.1,5.0,3); 
//
//  string histname = Form("nsig_parameters_cent%d_charge%d_eta%d_phi%d",icent,icharge,ieta,iphi);
//  parameterHist = (TH1F*)File->Get(histname.c_str());
//
//  for(int i_par = 0; i_par < 3; i_par++)
//  {
//    nsig_PID_pT->SetParameter(parameterHist->GetBinContent(1+i_par));
//  }
// 
//  return nsig_PID_pT;
//}


void readTofEff(int energy)
{
  // string inputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  //string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/ToF/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str());
  //string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/ToF/Eff_%s_ToFMatch_DCA2_nsigma0p5.root",vmsa::mBeamEnergy[energy].c_str());
  string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/ToF/Eff_%s_ToFMatch_FinerEta_240bins_rebinto48phibins.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  cout << "OPEN Efficiency File for K+ and K-: " << inputfile.c_str() << endl;

  h_FrameEta_ToF[0] = (TH1D*)File_InPut->Get("h_FrameEta_ToF")->Clone();
  h_FramePhi_ToF[0] = (TH1D*)File_InPut->Get("h_FramePhi_ToF")->Clone();
  std::cout << "NUMBER OF ETABINS = " << h_FrameEta_ToF[0]->GetNbinsX() << endl;;
  std::cout << "NUMBER OF PHIBINS = " << h_FramePhi_ToF[0]->GetNbinsX() << endl;;
  int neta = h_FrameEta_ToF[0]->GetNbinsX();
  int nphi = h_FramePhi_ToF[0]->GetNbinsX();
  h_FrameEta_ToF[1] = (TH1D*)File_InPut->Get("h_FrameEta_ToF")->Clone();
  h_FramePhi_ToF[1] = (TH1D*)File_InPut->Get("h_FramePhi_ToF")->Clone();

  for(int i_cent = 9; i_cent < 10; ++i_cent)
  {
    for(int i_eta = 0; i_eta < neta; ++i_eta)
    {
      string KEY;
      KEY = Form("h_mEfficiency_Kplus_Cent_%d_Eta_%d",i_cent,i_eta);
      h_TofKplus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());
      cout << "Kplus KEY: " << KEY.c_str() << endl;

      KEY = Form("h_mEfficiency_Kminus_Cent_%d_Eta_%d",i_cent,i_eta);
      h_TofKminus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());
      cout << "Kminus KEY: " << KEY.c_str() << endl;
      for(int i_phi = 0; i_phi < nphi; ++i_phi)
      {
	// string KEY;
	KEY = Form("h_mEfficiency_Kplus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_TofKplus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());
	cout << "Kplus KEY: " << KEY.c_str() << endl;

	KEY = Form("h_mEfficiency_Kminus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_TofKminus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());
	cout << "Kminus KEY: " << KEY.c_str() << endl;
      }	
    }
  }

}

void readTofEffFit(int energy)
{
  // string inputKplus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/FitPar_AuAu%s_Kplus_first.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());

  // string inputKminus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/FitPar_AuAu%s_Kminus_first.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());

  TFile *File_Kplus[10];
  TFile *File_Kminus[10];

  int neta = h_FrameEta_ToF[0]->GetNbinsX();
  int nphi = h_FramePhi_ToF[0]->GetNbinsX();
  for(int i_cent = 9; i_cent < 10; ++i_cent)
  {

    //string inputKplus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/ToF/FitPar_AuAu%s_Kplus_cent%d_DCA2_nsigma0p5.root",vmsa::mBeamEnergy[energy].c_str(),i_cent);
    string inputKplus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/ToF/FitPar_AuAu%s_Kplus_cent%d_FinerEta_240bins_rebinto48phibins.root",vmsa::mBeamEnergy[energy].c_str(),i_cent);
    //string inputKplus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/ToF/FitPar_AuAu%s_Kplus_cent%d.root",vmsa::mBeamEnergy[energy].c_str(),i_cent);
    File_Kplus[i_cent] = TFile::Open(inputKplus.c_str());
    cout << "OPEN ToF Matching Efficiency Fit File for K+: " << inputKplus.c_str() << endl;

    //string inputKminus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/ToF/FitPar_AuAu%s_Kminus_cent%d_DCA2_nsigma0p5.root",vmsa::mBeamEnergy[energy].c_str(),i_cent);
    string inputKminus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/ToF/FitPar_AuAu%s_Kminus_cent%d_FinerEta_240bins_rebinto48phibins.root",vmsa::mBeamEnergy[energy].c_str(),i_cent);
    //string inputKminus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/ToF/FitPar_AuAu%s_Kminus_cent%d.root",vmsa::mBeamEnergy[energy].c_str(),i_cent);
    File_Kminus[i_cent] = TFile::Open(inputKminus.c_str());
    cout << "OPEN ToF Matching Efficiency Fit File for K-: " << inputKminus.c_str() << endl;
      
    

    for(int i_eta = 0; i_eta < neta; ++i_eta) // eta
    {
      TH1D *h_Kplus = NULL;
      string KEY;
      KEY = Form("h_mFitParameters_Kplus_Cent_%d_Eta_%d",i_cent,i_eta);
      h_Kplus = (TH1D*)File_Kplus[i_cent]->Get(KEY.c_str());
      //cout << "Kplus KEY: " << KEY.c_str() << endl;
      KEY = Form("f_mToFMatch_Kplus_Cent_%d_Eta_%d",i_cent,i_eta);
      f_TofKplus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
      for(int i_par = 0; i_par < 7; ++i_par)
      {
        f_TofKplus[KEY]->FixParameter(i_par,h_Kplus->GetBinContent(i_par+1));
        // if(i_par < 2) cout << "Kplus: i_par = " << i_par << ", par = " << f_TofKplus[KEY]->GetParameter(i_par) << endl;
      }

      TH1D *h_Kminus = NULL;
      KEY = Form("h_mFitParameters_Kminus_Cent_%d_Eta_%d",i_cent,i_eta);
      h_Kminus = (TH1D*)File_Kminus[i_cent]->Get(KEY.c_str());
      //cout << "Kminus KEY: " << KEY.c_str() << endl;
      KEY = Form("f_mToFMatch_Kminus_Cent_%d_Eta_%d",i_cent,i_eta);
      f_TofKminus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
      for(int i_par = 0; i_par < 7; ++i_par)
      {
        f_TofKminus[KEY]->FixParameter(i_par,h_Kminus->GetBinContent(i_par+1));
        // if(i_par < 2) cout << "Kminus: i_par = " << i_par << ", par = " << f_TofKminus[KEY]->GetParameter(i_par) << endl;
      }
    }

    for(int i_eta = 0; i_eta < neta; ++i_eta) // eta & phi
    {
      for(int i_phi = 0; i_phi < nphi; ++i_phi)
      {
        TH1D *h_Kplus = NULL;
        string KEY;
        KEY = Form("h_mFitParameters_Kplus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
        h_Kplus = (TH1D*)File_Kplus[i_cent]->Get(KEY.c_str());
        //cout << "Kplus KEY: " << KEY.c_str() << endl;
        KEY = Form("f_mToFMatch_Kplus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
        f_TofKplus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
        for(int i_par = 0; i_par < 7; ++i_par)
        {
          f_TofKplus[KEY]->FixParameter(i_par,h_Kplus->GetBinContent(i_par+1));
          // if(i_par < 2) cout << "Kplus: i_par = " << i_par << ", par = " << f_TofKplus[KEY]->GetParameter(i_par) << endl;
        }

        TH1D *h_Kminus = NULL;
        KEY = Form("h_mFitParameters_Kminus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
        h_Kminus = (TH1D*)File_Kminus[i_cent]->Get(KEY.c_str());
        //cout << "Kminus KEY: " << KEY.c_str() << endl;
        KEY = Form("f_mToFMatch_Kminus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
        f_TofKminus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
        for(int i_par = 0; i_par < 7; ++i_par)
        {
          f_TofKminus[KEY]->FixParameter(i_par,h_Kminus->GetBinContent(i_par+1));
          // if(i_par < 2) cout << "Kminus: i_par = " << i_par << ", par = " << f_TofKminus[KEY]->GetParameter(i_par) << endl;
        }
      }	
    }
  }
}

void findHist_ToF(TLorentzVector const& lKaon, int iParticleIndex, float Psi2, int& EtaBin, int& PhiBin)
{
  float eta = lKaon.Eta();
  EtaBin = h_FrameEta_ToF[iParticleIndex]->FindBin(eta)-1;
  // cout << "eta = " << eta << ", EtaBin = " << EtaBin << endl;

  //float phi = lKaon.Phi()-Psi2;
  //float phi_shift = AngleShift(phi);
  //PhiBin = h_FramePhi_ToF[iParticleIndex]->FindBin(phi_shift)-1;

  float phi = lKaon.Phi();
  //cout << "kaon phi = " << phi << endl;
  PhiBin = h_FramePhi_ToF[iParticleIndex]->FindBin(phi)-1;
  // cout << "phi = " << phi << ", PhiBin = " << PhiBin << endl;
}

bool tofReconstructed(int iParticleIndex, int cent, float Psi2, TLorentzVector const& lKaon, bool phiswitch)
{
   if(fabs(lKaon.Eta()) >= vmsa::mEtaMax) return false;

   TH1D *h_ToF = NULL;
   TF1 *f_ToF = NULL;
   int EtaBin_ToF = -1;
   int PhiBin_ToF = -1;
  
   if(fabs(lKaon.Eta()) <= 1.0) findHist_ToF(lKaon,iParticleIndex,Psi2,EtaBin_ToF,PhiBin_ToF);

   if(phiswitch == 1)
   {
     if (iParticleIndex == 0)
     {
       if(fabs(lKaon.Eta()) <= 1.0)
       {
         string KEY_ToF = Form("h_mEfficiency_Kplus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
         //cout << KEY_ToF << endl;
         h_ToF = h_TofKplus[KEY_ToF];

         string KEY_ToFFit; 
         KEY_ToFFit = Form("f_mToFMatch_Kplus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
         //cout << KEY_ToFFit << endl;
         f_ToF = f_TofKplus[KEY_ToFFit]; // only 20-60%
       }
     }
     else
     {
       if(fabs(lKaon.Eta()) <= 1.0)
       {
         string KEY_ToF = Form("h_mEfficiency_Kminus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
         //cout << KEY_ToF << endl;
         h_ToF = h_TofKminus[KEY_ToF];

         string KEY_ToFFit;
         KEY_ToFFit = Form("f_mToFMatch_Kminus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
         //cout << KEY_ToFFit << endl;
         f_ToF = f_TofKminus[KEY_ToFFit]; // only 20-60%
       }
     }
   }
   if(phiswitch == 0)
   {
     if (iParticleIndex == 0)
     {
       if(fabs(lKaon.Eta()) <= 1.0)
       {
         string KEY_ToF = Form("h_mEfficiency_Kplus_Cent_%d_Eta_%d",cent,EtaBin_ToF); // get ToF eff @ 20-60%
         h_ToF = h_TofKplus[KEY_ToF];

         string KEY_ToFFit; 
         KEY_ToFFit = Form("f_mToFMatch_Kplus_Cent_%d_Eta_%d",cent,EtaBin_ToF); // get ToF eff @ 20-60%
         f_ToF = f_TofKplus[KEY_ToFFit]; // only 20-60%
       }
     }
     else
     {
       if(fabs(lKaon.Eta()) <= 1.0)
       {
         string KEY_ToF = Form("h_mEfficiency_Kminus_Cent_%d_Eta_%d",cent,EtaBin_ToF); // get ToF eff @ 20-60%
         h_ToF = h_TofKminus[KEY_ToF];

         string KEY_ToFFit;
         KEY_ToFFit = Form("f_mToFMatch_Kminus_Cent_%d_Eta_%d",cent,EtaBin_ToF); // get ToF eff @ 20-60%
         f_ToF = f_TofKminus[KEY_ToFFit]; // only 20-60%
       }
     }
   }
   double pt = lKaon.Perp();

   bool is_ToF = false;
   if(fabs(lKaon.Eta()) <= 1.0)
   { 
     int const bin_ToF = h_ToF->FindBin(pt); // tof fit and hist combined
     double prob_tof = 0.0;
     if(pt > 0.3) prob_tof = f_ToF->Eval(pt); // donot use fit extrapolation
     //if(pt > 0.2) prob_tof = f_ToF->Eval(pt); // donot use fit extrapolation
     else prob_tof = h_ToF->GetBinContent(bin_ToF);
     is_ToF = gRandom->Rndm() < prob_tof;
   }

   if(fabs(lKaon.Eta()) < 1.0) return is_ToF;
}



void makePicoDstQA_Phi_TTrees(TString InputFileList, Int_t nEvents = 0, TString OutputFile = "test.histo", TString jobId = "1", Int_t mEnergy = 4, Int_t mGid = 11, Int_t mPid = 0, bool phiswitch = 1, Float_t rho00 = 1./3., Float_t real = 0.0, Float_t imag = 0.0, Float_t rerho1n1 = 0.0, Float_t imrho1n1 = 0.0, Float_t hrho00 = 1./3., Float_t hreal = 0.0, Float_t himag = 0.0, Float_t hrerho1n1 = 0.0, Float_t himrho1n1 = 0.0);

void makePicoDstQA_Phi_TTrees(TString InputFileList, Int_t nEvents, TString OutputFile, TString jobId, Int_t mEnergy, Int_t mGid, Int_t mPid, bool phiswitch, Float_t rho00, Float_t real, Float_t imag, Float_t rerho1n1, Float_t imrho1n1, Float_t hrho00, Float_t hreal, Float_t himag, Float_t hrerho1n1, Float_t himrho1n1) 
{
  int mCent;
  float mEpFull;
  float mRefMult;
  int mNToFMatch;
  float mVx;
  float mVy;
  float mVz;
  float mWeight;

  TLorentzVector mLMeson;
  TLorentzVector mLKp;
  TLorentzVector mLKm; 
  Bool_t mRcExists;
  TLorentzVector mLRcKp;
  TLorentzVector mLRcKm;

  int mNHitsFitP;
  int mNHitsMaxP;
  float mDEdxP;
  float mDcaP;
  int mNHitsFitM;
  int mNHitsMaxM;
  float mDEdxM;
  float mDcaM;

  float mMcPhiStar;
  float mMcCosThetaStar;
  float mMcPhiPrime;
  float mMcCosTheta;
  float mMcHelicityAngle;
  
  float mRcPhiStar;
  float mRcCosThetaStar;
  float mRcPhiPrime;
  float mRcCosTheta;
  float mRcHelicityAngle;
  
  Bool_t mPassTpcKp;
  Bool_t mPassTofKp;
  Bool_t mPassM2Kp;
  Bool_t mPassNsigKp;

  Bool_t mPassTpcKm;
  Bool_t mPassTofKm;
  Bool_t mPassM2Km;
  Bool_t mPassNsigKm;
 
  float mNSigmaKaonP;
  float mNSigmaKaonM;
 
  TTree *mKaonTree = new TTree("kaontree","Tree for kaon info");

  mKaonTree->Branch("cent", &mCent, "cent/I");
  mKaonTree->Branch("epfull", &mEpFull, "epfull/F");
  mKaonTree->Branch("refmult", &mRefMult, "refmult/F");
  mKaonTree->Branch("ntofmatch", &mNToFMatch, "ntofmatch/I");
  mKaonTree->Branch("vx", &mVx, "vz/F");
  mKaonTree->Branch("vy", &mVy, "vz/F");
  mKaonTree->Branch("vz", &mVz, "vz/F");
  mKaonTree->Branch("weight", &mWeight, "weight/F");

  mKaonTree->Branch("lmeson", &mLMeson);
  mKaonTree->Branch("lkp", &mLKp);
  mKaonTree->Branch("lkm", &mLKm);
  mKaonTree->Branch("rcexists", &mRcExists, "rcexists/O");
  mKaonTree->Branch("rlkp", &mLRcKp);
  mKaonTree->Branch("rlkm", &mLRcKm);

  mKaonTree->Branch("nhitsfitp", &mNHitsFitP, "nhitsfitp/I");
  mKaonTree->Branch("nhitsmaxp", &mNHitsMaxP, "nhitsmaxp/I");
  mKaonTree->Branch("dedxp", &mDEdxP, "dedxp/F");
  mKaonTree->Branch("dcap", &mDcaP, "dcap/F");
  mKaonTree->Branch("nhitsfitm", &mNHitsFitM, "nhitsfitm/I");
  mKaonTree->Branch("nhitsmaxm", &mNHitsMaxM, "nhitsmaxm/I");
  mKaonTree->Branch("dedxm", &mDEdxM, "dedxm/F");
  mKaonTree->Branch("dcam", &mDcaM, "dcam/F");

  mKaonTree->Branch("mcphistar", &mMcPhiStar, "mcphistar/F");
  mKaonTree->Branch("mccosthetastar", &mMcCosThetaStar, "mccosthetastar/F");
  mKaonTree->Branch("mcphiprime", &mMcPhiPrime, "mcphiprime/F");
  mKaonTree->Branch("mccostheta", &mMcCosTheta, "mccostheta/F");
  mKaonTree->Branch("mchelicityangle", &mMcHelicityAngle, "mchelicityangle/F");
  
  mKaonTree->Branch("rcphistar", &mRcPhiStar, "rcphistar/F");
  mKaonTree->Branch("rccosthetastar", &mRcCosThetaStar, "rccosthetastar/F");
  mKaonTree->Branch("rcphiprime", &mRcPhiPrime, "rcphiprime/F");
  mKaonTree->Branch("rccostheta", &mRcCosTheta, "rccostheta/F");
  mKaonTree->Branch("rchelicityangle", &mRcHelicityAngle, "rchelicityangle/F");
  
  mKaonTree->Branch("passtpckp", &mPassTpcKp, "passtpckp/O");
  mKaonTree->Branch("passtofkp", &mPassTofKp, "passtofkp/O");
  mKaonTree->Branch("passm2kp", &mPassM2Kp, "passm2kp/O");
  mKaonTree->Branch("passnskp", &mPassNsigKp, "passnskp/O");

  mKaonTree->Branch("passtpckm", &mPassTpcKm, "passtpckm/O");
  mKaonTree->Branch("passtofkm", &mPassTofKm, "passtofkm/O");
  mKaonTree->Branch("passm2km", &mPassM2Km, "passm2km/O");
  mKaonTree->Branch("passnskm", &mPassNsigKm, "passnskm/O");

  mKaonTree->Branch("nsigmakaonp", &mNSigmaKaonP, "nsigmakaonp/F");
  mKaonTree->Branch("nsigmakaonm", &mNSigmaKaonM, "nsigmakaonm/F");
  
  // Load libraries for CINT mode
  #ifdef __CINT__
    gROOT  -> Macro("loadMuDst.C");
    gSystem->Load("StRefMultCorr");
    gSystem->Load("StPicoEvent");
    gSystem->Load("StPicoDstMaker");
    gSystem->Load("StVecMesonMaker");
    gSystem->Load("StEffKaon");
    //gSystem->Load("StUtility");
  #endif

  std::string inputfilepthelicity = Form("StRoot/Utility/Rho/%s/AccPhiPtSys_eta1_eta1_PolySys_Helicity_2D_OffDiag_Embed.root",vmsa::mBeamEnergy[mEnergy].c_str());
  TFile *filepthelicity = TFile::Open(inputfilepthelicity.c_str());
  TGraphAsymmErrors *g_rho00_helicity = (TGraphAsymmErrors*) filepthelicity->Get("rhoRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_real_helicity = (TGraphAsymmErrors*) filepthelicity->Get("realRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_imag_helicity = (TGraphAsymmErrors*) filepthelicity->Get("imagRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_rerho1n1_helicity = (TGraphAsymmErrors*) filepthelicity->Get("rerho1n1Raw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_imrho1n1_helicity = (TGraphAsymmErrors*) filepthelicity->Get("imrho1n1Raw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");

  g_rho00_helicity->Print();
  g_real_helicity->Print(); 
  g_imag_helicity->Print(); 
  g_rerho1n1_helicity->Print(); 
  g_imrho1n1_helicity->Print(); 



  std::string yinputfilepthelicity = Form("StRoot/Utility/Rho/%s/AccPhiPtSys_eta1_eta1_PolySys_Helicity_2D_OffDiag_Embed_Rapidity.root",vmsa::mBeamEnergy[mEnergy].c_str());
  TFile *yfilepthelicity = TFile::Open(yinputfilepthelicity.c_str());
  TGraphAsymmErrors *g_rho00_helicity_y = (TGraphAsymmErrors*) yfilepthelicity->Get("rhoRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_real_helicity_y = (TGraphAsymmErrors*) yfilepthelicity->Get("realRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_imag_helicity_y = (TGraphAsymmErrors*) yfilepthelicity->Get("imagRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_rerho1n1_helicity_y = (TGraphAsymmErrors*) yfilepthelicity->Get("rerho1n1Raw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_imrho1n1_helicity_y = (TGraphAsymmErrors*) yfilepthelicity->Get("imrho1n1Raw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");


  double helicityrho00[6] = {0.333333,0.333333,0.307032,0.295729,0.288274,0.246577};
  for(int ipt = 2; ipt < 6; ipt++)
  { 
    f_mRhoPt[ipt] = new TF1(Form("f_mRho_pt%d",ipt),SpinDensity,-1.0,1.0,2);
    f_mRhoPt_Helicity[ipt] = new TF1(Form("f_mRhoHelicity_pt%d",ipt),SpinDensity,-1.0,1.0,2);
    f_mRhoPt_2D[ipt] = new TF2(Form("f_mRho2D_pt%d",ipt),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),5);
    f_mRhoPt_Helicity_2D[ipt] = new TF2(Form("f_mRhoHelicity2D_pt%d",ipt),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),5);
    TF2 *temp = new TF2("temp",SpinDensity2Dcosneg,-1.0,1.0,0.0,2.0*TMath::Pi(),5);
    TF2 *htemp = new TF2("temp_helicity",SpinDensity2Dcosneg,-1.0,1.0,0.0,2.0*TMath::Pi(),5);
    double pt, rho; 
    //g_rho_pt->GetPoint(ipt,pt,rho);
    //cout << "ipt = " << ipt << ", pt = " << pt << ", rho = " << rho << endl;
    f_mRhoPt[ipt]->FixParameter(0,rho00);// set by user
    f_mRhoPt[ipt]->FixParameter(1,0.75);

    mMaxData[ipt] = f_mRhoPt[ipt]->GetMaximum(-1.0,1.0);
    cout << "Maximum of data rh00 = " << mMaxData[ipt] << endl; 

    //f_mRhoPt_Helicity[ipt]->FixParameter(0,helicityrho00[ipt]); // set by data
    f_mRhoPt_Helicity[ipt]->FixParameter(0,hrho00);       // set by user
    f_mRhoPt_Helicity[ipt]->FixParameter(1,0.75);
  
    mMaxHelicity[ipt] = f_mRhoPt_Helicity[ipt]->GetMaximum(-1.0,1.0);
    cout << "Maximum of helicity rh00 = " << mMaxHelicity[ipt] << endl; 

    double hpt;
    double hrho00_data, hreal_data, himag_data, hrerho1n1_data, himrho1n1_data;
    double hrho00_data_err, hreal_data_err, himag_data_err, hrerho1n1_data_err, himrho1n1_data_err;
    g_rho00_helicity->GetPoint(ipt-2,hpt,hrho00_data); 
    g_real_helicity->GetPoint(ipt-2,hpt,hreal_data);
    g_imag_helicity->GetPoint(ipt-2,hpt,himag_data);
    g_rerho1n1_helicity->GetPoint(ipt-2,hpt,hrerho1n1_data);
    g_imrho1n1_helicity->GetPoint(ipt-2,hpt,himrho1n1_data);

    //hrho00_data_err = g_rho00_helicity->GetErrorYhigh(ipt-2); 
    //hreal_data_err = g_real_helicity->GetErrorYhigh(ipt-2);
    //himag_data_err = g_imag_helicity->GetErrorYhigh(ipt-2);
    //hrerho1n1_data_err = g_rerho1n1_helicity->GetErrorYhigh(ipt-2);
    //himrho1n1_data_err = g_imrho1n1_helicity->GetErrorYhigh(ipt-2);

    f_mRhoPt_Helicity_2D[ipt]->FixParameter(0,hrho00_data);
    f_mRhoPt_Helicity_2D[ipt]->FixParameter(1,hreal_data);
    f_mRhoPt_Helicity_2D[ipt]->FixParameter(2,himag_data);
    f_mRhoPt_Helicity_2D[ipt]->FixParameter(3,hrerho1n1_data);
    f_mRhoPt_Helicity_2D[ipt]->FixParameter(4,himrho1n1_data);

    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(0,hrho00_data);//+hrho00_data_err);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(1,hreal_data);//-hreal_data_err);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(2,himag_data);//-himag_data_err);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(3,hrerho1n1_data);//-hrerho1n1_data_err);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(4,himrho1n1_data);//-himrho1n1_data_err);

    cout << "ipt = " << ipt << endl;
    cout << "helicity rho00 = " << hrho00_data << endl;
    cout << "helicity real = " << hreal_data << endl;
    cout << "helicity imag = " << himag_data << endl;
    cout << "helicity rerho1n1 = " << hrerho1n1_data << endl;
    cout << "helicity imrho1n1 = " << himrho1n1_data << endl;

    htemp->FixParameter(0,hrho00_data);
    htemp->FixParameter(1,hreal_data);
    htemp->FixParameter(2,himag_data);
    htemp->FixParameter(3,hrerho1n1_data);
    htemp->FixParameter(4,himrho1n1_data);


    //f_mRhoPt_2D[ipt]->FixParameter(0,rho);
    f_mRhoPt_2D[ipt]->FixParameter(0,rho00);
    f_mRhoPt_2D[ipt]->FixParameter(1,real);
    f_mRhoPt_2D[ipt]->FixParameter(2,imag);
    f_mRhoPt_2D[ipt]->FixParameter(3,rerho1n1);
    f_mRhoPt_2D[ipt]->FixParameter(4,imrho1n1);

    temp->FixParameter(0,rho00);
    temp->FixParameter(1,real);
    temp->FixParameter(2,imag);
    temp->FixParameter(3,rerho1n1);
    temp->FixParameter(4,imrho1n1);

    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(0,hrho00);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(1,hreal);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(2,himag);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(3,hrerho1n1);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(4,himrho1n1);

    //htemp->FixParameter(0,hrho00);
    //htemp->FixParameter(1,hreal);
    //htemp->FixParameter(2,himag);
    //htemp->FixParameter(3,hrerho1n1);
    //htemp->FixParameter(4,himrho1n1);

    //mMax = f_mRhoPt_2D[ipt]->GetMaximum();
    //double x[2] = {TMath::Pi()/2.,TMath::Pi()/2.};
    double x,y;
    temp->GetMinimumXY(x,y);
    cout << "Maximum x,y = " << x << "," << y << endl;
    mMax2D[ipt] = TMath::Abs(temp->Eval(x,y));
    cout << "Maximum of the function is " << mMax2D[ipt] << endl;

    double hx,hy;
    htemp->GetMinimumXY(hx,hy);
    cout << "Maximum hx,hy = " << hx << "," << hy << endl;
    mMaxHelicity2D[ipt] = TMath::Abs(htemp->Eval(hx,hy));
    cout << "Maximum of the function is " << mMaxHelicity2D[ipt] << endl;

    //cout << "x,y both = TMath::Pi()/2, eval func = " << temp->Eval(1.57113,1.5708) << endl;
  }  
  for(int iy = 0; iy < 10; iy++)
  { 
    f_mRhoY_Helicity_2D[iy] = new TF2(Form("f_mRhoHelicity2D_y%d",iy),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),5);
    TF2 *htemp = new TF2("temp_helicity",SpinDensity2Dcosneg,-1,1,0.0,2.0*TMath::Pi(),5);

    double hpt;
    double hrho00_data, hreal_data, himag_data, hrerho1n1_data, himrho1n1_data;
    double hrho00_data_err, hreal_data_err, himag_data_err, hrerho1n1_data_err, himrho1n1_data_err;
    g_rho00_helicity_y->GetPoint(iy,hpt,hrho00_data); 
    g_real_helicity_y->GetPoint(iy,hpt,hreal_data);
    g_imag_helicity_y->GetPoint(iy,hpt,himag_data);
    g_rerho1n1_helicity_y->GetPoint(iy,hpt,hrerho1n1_data);
    g_imrho1n1_helicity_y->GetPoint(iy,hpt,himrho1n1_data);

    hrho00_data_err = g_rho00_helicity_y->GetErrorYhigh(iy); 
    hreal_data_err = g_real_helicity_y->GetErrorYhigh(iy);
    himag_data_err = g_imag_helicity_y->GetErrorYhigh(iy);
    hrerho1n1_data_err = g_rerho1n1_helicity_y->GetErrorYhigh(iy);
    himrho1n1_data_err = g_imrho1n1_helicity_y->GetErrorYhigh(iy);

    f_mRhoY_Helicity_2D[iy]->FixParameter(0,hrho00_data);
    f_mRhoY_Helicity_2D[iy]->FixParameter(1,hreal_data);
    f_mRhoY_Helicity_2D[iy]->FixParameter(2,himag_data);
    f_mRhoY_Helicity_2D[iy]->FixParameter(3,hrerho1n1_data);
    f_mRhoY_Helicity_2D[iy]->FixParameter(4,himrho1n1_data);

    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(0,hrho00_data+hrho00_data_err);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(1,hreal_data-hreal_data_err);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(2,himag_data-himag_data_err);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(3,hrerho1n1_data-hrerho1n1_data_err);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(4,himrho1n1_data-himrho1n1_data_err);

    cout << "iy = " << iy << endl;
    cout << "helicity rho00 = " << hrho00_data << endl;
    cout << "helicity real = " << hreal_data << endl;
    cout << "helicity imag = " << himag_data << endl;
    cout << "helicity rerho1n1 = " << hrerho1n1_data << endl;
    cout << "helicity imrho1n1 = " << himrho1n1_data << endl;

    htemp->FixParameter(0,hrho00_data);
    htemp->FixParameter(1,hreal_data);
    htemp->FixParameter(2,himag_data);
    htemp->FixParameter(3,hrerho1n1_data);
    htemp->FixParameter(4,himrho1n1_data);

    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(0,hrho00);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(1,hreal);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(2,himag);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(3,hrerho1n1);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(4,himrho1n1);

    //htemp->FixParameter(0,hrho00);
    //htemp->FixParameter(1,hreal);
    //htemp->FixParameter(2,himag);
    //htemp->FixParameter(3,hrerho1n1);
    //htemp->FixParameter(4,himrho1n1);

    double hx,hy;
    htemp->GetMinimumXY(hx,hy);
    cout << "Maximum hx,hy = " << hx << "," << hy << endl;
    mMaxHelicity2DY[iy] = TMath::Abs(htemp->Eval(hx,hy));
    cout << "Maximum of the function is " << mMaxHelicity2DY[iy] << endl;

    //cout << "x,y both = TMath::Pi()/2, eval func = " << temp->Eval(1.57113,1.5708) << endl;
  }  

  readTofEff(mEnergy);
  readTofEffFit(mEnergy);

  readm2PID();
  readnsigPID();

  readm2PID_hist();
  readnsigPID_hist();

  // List of member links in the chain
  StChain*                    chain  =  new StChain ;

  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(StPicoDstMaker::IoRead,InputFileList,"picoDst");

  StUtility* mUtility = new StUtility(mEnergy);
  mUtility->initRunIndex(); // initialize std::map for run index 
  //mUtility->initEventPlane(); // initialize std::map for event planes 

  StRefMultCorr *mRefMultCorr = new StRefMultCorr("refmult");

  StVecMesonCut* mVecMesonCut = new StVecMesonCut(mEnergy);

  StEffHistManger* mEffHistManger = new StEffHistManger();

  if ( nEvents == 0 )  nEvents = 1000000000 ;       // Take all events in nFiles if nEvents = 0

  // ---------------- modify here according to your QA purpose --------------------------
  //TFile *tags_output = new TFile( "/gpfs01/star/scratch/gwilks3/"+OutputFile+"_"+jobId+".root" , "recreate" ) ;
  TFile *tags_output = new TFile(OutputFile+"_"+jobId+".root" , "recreate" ) ;
  tags_output->cd();

  // ---------------- end of histogram and tree booking --------------------------------

  // Loop over the links in the chain
  Int_t iInit = chain -> Init() ;
  if (iInit) chain->FatalErr(iInit,"on init");

  //mEffHistManger->InitHist();
  //mEffHistManger->InitHistQA();
  int nMC[5] = {0};
  int nRC[5] = {0};
  
  // chain -> EventLoop(1,nEvents) ;  //will output lots of useless debugging info.
  Int_t istat = 0, i = 1;
  while (i <= nEvents && istat != 2) {
     if(i%1000==0)cout << endl << "== Event " << i << " start ==" << endl;
     chain->Clear();
     //cout << "Event = " << i << endl;
     istat = chain->Make(i);

     if (istat == 2)
	  cout << "Last  event processed. Status = " << istat << endl;
     if (istat == 3)
	  cout << "Error event processed. Status = " << istat << endl;
     i++;

     if(istat != kStOK)continue; //skip those suspectible events
     
  // ---------------- modify here according to your QA purpose --------------------------
     //cout<<"In event #. "<<i-1<<" Maker status "<<istat<<endl;

     StPicoDst* mPicoDst = picoDstMaker->picoDst();
     if(!mPicoDst) {
	  LOG_WARN << " No PicoDst " << endm; continue;
     }

     StPicoEvent* mPicoEvent = mPicoDst->event();
     if(!mPicoEvent) {
	  LOG_WARN << " No PicoEvent " << endm; continue;
     }

     //triggerdd
     //if ( ! mPicoEvent->isTrigger(640001) && ! mPicoEvent->isTrigger(640011) && ! mPicoEvent->isTrigger(640021) && ! mPicoEvent->isTrigger(640031) && ! mPicoEvent->isTrigger(640041) && ! mPicoEvent->isTrigger(640051) ) continue ; // BES-II 19.6 GeV
    
     Int_t runId = mPicoEvent->runId();
     Int_t eventId = mPicoEvent->eventId();
     Int_t refMult = mPicoEvent->refMult();
     Float_t vx = mPicoEvent->primaryVertex().x();
     Float_t vy = mPicoEvent->primaryVertex().y();
     Float_t vz = mPicoEvent->primaryVertex().z();
     Float_t zdcX = mPicoEvent->ZDCx();
     const unsigned short nBTofMatched = mPicoEvent->nBTOFMatch();

     mRefMultCorr->init(runId);

     if(mRefMultCorr->isBadRun( runId ))
     {
       LOG_ERROR << "Bad Run! Skip!" << endm; continue;
     }
   
     bool isPileUpEvent = false;
     // IMPORTANT: vertex position is needed for Au+Au 19.6 GeV 2019
     if (mRefMultCorr->isPileUpEvent( refMult, nBTofMatched, vz ) ) isPileUpEvent = true;
     mRefMultCorr->initEvent(refMult,vz,zdcX);

     const Int_t cent9 = mRefMultCorr->getCentralityBin9();       // 0: 70-80%, 1: 60-70%, ..., 6: 10-20%, 7: 5-10%, 8: 0-5%
     const Double_t reweight = mRefMultCorr->getWeight();           // Retrieve weight
 
     if(cent9 < -0.5) continue;

     //TEMPORARY CENT CUT FOR 20-60% CENTRALITY
     //if(cent9 < 2 || cent9 > 5) continue;

     if(isPileUpEvent) continue;
     if(!mVecMesonCut->passEventCut(mPicoDst)) continue;

     //const int runIndex = mUtility->findRunIndex(runId);
     // 
     //if(runIndex < 0)
     //{
     //  LOG_ERROR << "Could not find this run Index from StUtility! Skip!" << endm; continue;
     //}

     const float ep_west = 0.0;
     const float ep_east = 0.0;
     //const float ep_full = 0.0;
     //const float ep_west = mUtility->findEPwest(eventId);
     //const float ep_east = mUtility->findEPeast(eventId);
     const float ep_full = 0.0;//mUtility->findEPfull(eventId);
     if( ep_west < -900.0 || ep_east < -900.0 || ep_full < -900.0 )
     {
       LOG_ERROR << "Could not find this event plane from StUtility! Skip!" << endm; continue;
     }

     //mEffHistManger->FillEventHistQA(reweight, cent9, vx, vy, vz, nBTofMatched, refMult);

     //Vz
     //if ( fabs(mPicoEvent->primaryVertex().Z()) > 70.0 ) continue ;
     //Vr
     //if ( mPicoEvent->primaryVertex().Perp() > 2.0 ) continue ;
     
     //fill MC histograms
     //The MC arrays in PicoDst
     Int_t NoMuMcVertices = mPicoDst->numberOfMcVertices();
     Int_t NoMuMcTracks = mPicoDst->numberOfMcTracks();
     //mPicoDst->printMcVertices();
     //mPicoDst->printMcTracks();
     //LOG_INFO <<"# of MC tracks = "<< NoMuMcTracks <<" # of MC vertices = "<< NoMuMcVertices << endm;
     if (! NoMuMcVertices || ! NoMuMcTracks) {
	  //LOG_WARN << "Ev. " << i  << " has no MC information ==> skip it" << endm;
	  continue;
     }
     Int_t nMc = 0;

     std::map< int, std::pair<StPicoMcTrack*, StPicoMcTrack*> > phiMesonMcDaughters;
     std::map< int, StPicoMcTrack*> phiMeson;
     std::map< int, int> nMCDaughters;
     std::map< int, int> nRCDaughters;
     // Loop for MC tracks
     for(Int_t itrk=0; itrk<NoMuMcTracks; itrk++){
	  StPicoMcTrack *mcTrack = (StPicoMcTrack *) mPicoDst->mcTrack(itrk);
	  if (! mcTrack) continue;

	  // Select only Triggered Mc Vertex, i.e. the MC track should originate from PV (IdVx=1)
	  Int_t idMcVx = mcTrack->idVtxStart();
	  //if (IdVx != 1) continue;
	  //if(mcTrack->geantId()==10151 && idMcVx == 1)// all phi-mesons (no channel selected)
          //{
          //  TLorentzVector lphi = mcTrack->fourMomentum();
          //  //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,lphi.Pt(),lphi.Phi(),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0);
          //  
          //}

          //Only consider K+ and K- for selected channel
	  if(mcTrack->geantId() != 11 && mcTrack->geantId() != 12) continue; //geantId cut for Kplus = 11 and Kminus = 12

          Int_t idMcTrack = -1;
          Int_t parentGeantId = -1;
	  while (idMcVx != 1) {
	     StPicoMcVertex *mcVertex = (StPicoMcVertex *) mPicoDst->mcVertex(idMcVx-1);
	     idMcTrack = mcVertex->idOfParentTrack();
	     if (! idMcTrack) break;
	     StPicoMcTrack *mcTrackP = (StPicoMcTrack *) mPicoDst->mcTrack(idMcTrack-1);
             parentGeantId = mcTrackP->geantId();
             //cout << "Parent geantid = " << parentGeantId << endl;
             //cout << "idMcTrack = " << idMcTrack << endl;
             if(parentGeantId == 10151)
             {
               phiMeson[idMcTrack] = mcTrackP;
             }
 	     idMcVx = mcTrackP->idVtxStart();
	     if (! idMcVx) break;
	  }
	  if (idMcVx != 1) continue; //this MC track is not eventually originated from PV


          if(parentGeantId != 10151) continue;
          if(idMcTrack == -1) continue;

          if(mcTrack->geantId() == 11) 
          {
            //cout << "phi meson idMcTrack = " << idMcTrack << endl; 
            //cout << "kaon+ idVtxStart = " << mcTrack->idVtxStart() << endl;
            //cout << "phi   idVtxStop  = " << phiMeson[idMcTrack]->idVtxStop() << endl << endl;

            if( mcTrack->idVtxStart() != phiMeson[idMcTrack]->idVtxStop() ) continue;
            phiMesonMcDaughters[idMcTrack].first = mcTrack;
            nMCDaughters[idMcTrack]++;
            //phiMesonRcDaughters[idMcTrack].first = ptrack;
            //cout << "index = " << idMcTrack << " kplus mc" << endl;
          }
          if(mcTrack->geantId() == 12) 
          {
            //cout << "phi meson idMcTrack = " << idMcTrack << endl; 
            //cout << "kaon- idVtxStart = " << mcTrack->idVtxStart() << endl;
            //cout << "phi   idVtxStop  = " << phiMeson[idMcTrack]->idVtxStop() << endl << endl;
            if( mcTrack->idVtxStart() != phiMeson[idMcTrack]->idVtxStop() ) continue;
            phiMesonMcDaughters[idMcTrack].second = mcTrack;
            nMCDaughters[idMcTrack]++;
            //phiMesonRcDaughters[idMcTrack].second = ptrack;
            //cout << "index = " << idMcTrack << " kminus mc" << endl;
          }

     }

     //fill Event QA histograms
     Int_t nTracks = mPicoDst->numberOfTracks();

     //std::map< int, std::pair<StPicoMcTrack*, StPicoMcTrack*> > phiMesonMcDaughters;
     std::map< int, std::pair<StPicoTrack*, StPicoTrack*> > phiMesonRcDaughters;
 
     for(Int_t i=0; i<nTracks; i++)                // Main loop for Iterating over tracks
     {
          StPicoTrack* ptrack ; 
	  ptrack = mPicoDst->track(i);  // Pointer to a track
          //cout << "ptrack->isPrimary() = " << ptrack->isPrimary() << endl;
	  if(!ptrack->isPrimary()) continue;

	  if (ptrack->idTruth() <= 0 || ptrack->idTruth() > NoMuMcTracks) {
	     //cout << "Illegal idTruth " << ptrack->idTruth() << " The track is ignored" << endl;
	     continue;
	  }
	  StPicoMcTrack *mcTrack = (StPicoMcTrack *) mPicoDst->mcTrack(ptrack->idTruth()-1);
	  if (!mcTrack) {
	     LOG_WARN << "Inconsistency in mcArray(1), ignored" << endm;
	     continue;
	  }
          //mcTrack->Print();
	  if (mcTrack->id() != ptrack->idTruth()) {
	     LOG_WARN << "Mismatched idTruth " << ptrack->idTruth() << " and mcTrack Id " <<  mcTrack->id() 
		  << " this track is ignored" <<  endm;
	  }
	  Int_t idMcVx = mcTrack->idVtxStart();
	  Int_t idMcVxStop = mcTrack->idVtxStop();
	  Int_t idMcVxItrmd = mcTrack->idVtxItrmd();
          //cout << "idMcVx = " << idMcVx << endl;
          //cout << "idMcVxStop = " << idMcVxStop << endl;
          //cout << "idMcVxItrmd = " << idMcVxItrmd << endl;
          Int_t parentGeantId = -1;
          Int_t pdgId = -1;
          Int_t idMcTrack = -1;
          Int_t idMcTrackStop = -1;
	  //if(idMcVxStop > 0) 
          //{
          //   StPicoMcVertex *mcVertexStop = (StPicoMcVertex *) mPicoDst->mcVertex(idMcVxStop-1);
	  //   idMcTrackStop = mcVertexStop->idOfParentTrack();
	  //   StPicoMcTrack *mcTrackPstop = (StPicoMcTrack *) mPicoDst->mcTrack(idMcTrackStop-1);
          //   cout << "parentGeantId STOP= " << mcTrackPstop->geantId() << endl;
          //   cout << "idMcTrackStop = " << idMcTrackStop << endl;
          //}
	  while (idMcVx != 1) {
	     StPicoMcVertex *mcVertex = (StPicoMcVertex *) mPicoDst->mcVertex(idMcVx-1);
	     //StPicoMcVertex *mcVertexStop = (StPicoMcVertex *) mPicoDst->mcVertex(idMcVxStop-1);
	     //StPicoMcVertex *mcVertexStop = (StPicoMcVertex *) mPicoDst->mcVertex(idMcVxStop-1);
	     //StPicoMcVertex *mcVertex = (StPicoMcVertex *) mPicoDst->mcVertex(idMcVx);
	     //StPicoMcVertex *mcVertexStop = (StPicoMcVertex *) mPicoDst->mcVertex(idMcVxStop);
	     idMcTrack = mcVertex->idOfParentTrack();
	     //idMcTrackStop = mcVertexStop->idOfParentTrack();
	     //StPicoMcTrack *mcTrackPstop = (StPicoMcTrack *) mPicoDst->mcTrack(idMcTrackStop-1);
             //cout << "parentGeantId STOP= " << mcTrackPstop->geantId() << endl;
             //cout << "idMcTrackStop = " << idMcTrackStop << endl;
	     if (! idMcTrack) break;
	     StPicoMcTrack *mcTrackP = (StPicoMcTrack *) mPicoDst->mcTrack(idMcTrack-1);
             parentGeantId = mcTrackP->geantId();
             pdgId = mcTrackP->pdgId();
	     idMcVx = mcTrackP->idVtxStart();
	     if (! idMcVx) break;
	  }
	  if (idMcVx != 1) continue; //this MC track is not eventually originated from PV

          //cout << "parentGeantId = " << parentGeantId << endl;      
          //cout << "GeantId = " << mcTrack->geantId() << endl;      
	  if(mcTrack->geantId() != 11 && mcTrack->geantId() != 12) continue; //geantId cut for Kplus = 11 and Kminus = 12
 
          if (parentGeantId != 10151) continue;
          if (idMcTrack == -1) continue;

	  //if(mcTrack->geantId() != 11) continue; //geantId cut for Kplus = 11 and Kminus = 12
	  //if(mcTrack->idVtxStart() != 1) {
	  //   LOG_WARN<<"mc track may not directly originate from PV!"<<endm;
	  //}

          //double P   = mcTrack->p().Mag();
          //double Pt  = mcTrack->p().Perp();
          //double Phi = mcTrack->p().Phi();
          //double Eta = mcTrack->p().PseudoRapidity();
          //double Y = mcTrack->p().Rapidity();
   
          //double gRcPt = ptrack->gMom().Perp();
          //double pRcPt = ptrack->pMom().Perp();
 
          if(ptrack->qaTruth()<50.) continue;

          //cout << "nsigma_kaon = " << ptrack->nSigmaKaon() << endl;    

          //-------------------------McKaon-----------------------------------------------------
          //if( fabs(Eta) > vmsa::mEtaMax ) continue; // eta cut
          //if(!(Pt > vmsa::mGlobPtMin && P < vmsa::mPrimMomMax) ) continue; // pt cut

          if(mcTrack->geantId() == 11) 
          {
            //phiMesonMcDaughters[idMcTrack].first = mcTrack;
            if( mcTrack->idVtxStart() != phiMeson[idMcTrack]->idVtxStop() ) continue;
            phiMesonRcDaughters[idMcTrack].first = ptrack;
            nRCDaughters[idMcTrack]++;
            //cout << "index = " << idMcTrack << " kplus mc" << endl;
          }
          if(mcTrack->geantId() == 12) 
          {
            //phiMesonMcDaughters[idMcTrack].second = mcTrack;
            if( mcTrack->idVtxStart() != phiMeson[idMcTrack]->idVtxStop() ) continue;
            phiMesonRcDaughters[idMcTrack].second = ptrack;
            nRCDaughters[idMcTrack]++;
            //mVecMesonCut->getPrimaryMass2(ptrack, mPicoDst);
            //cout << "index = " << idMcTrack << " kminus mc" << endl;
          }

          //mEffHistManger->FillHistMc(cent9,Pt,Eta,Phi,ep_west,ep_east);
          //-------------------------McKaon----------------------------------------------------- 

          //if( !mVecMesonCut->passTrackMeson(ptrack, mPicoEvent, mPid) ) continue;
          //mEffHistManger->FillHistRc(cent9,Pt,Eta,Phi,ep_west,ep_east);
          //mEffHistManger->FillHistPt(cent9,Pt,gRcPt,pRcPt);   
       

	  //if(ptrack->qaTruth()<50.) continue; // Quality of MC track (probably out of 100)

	  //if(ptrack->nHits()<=15)continue; // nHits cut
	  //if(ptrack->flag()<=0)continue;
	  //if(abs(ptrack->charge())!=1) continue; //if the charge of the track is not equal to 1

	  //TVector3 p = ptrack->pMom();
	  //hPhi->Fill(p.Phi());
	  //hPt->Fill(p.Perp());
	  //hEta->Fill(p.PseudoRapidity());
	  //if(fabs(p.PseudoRapidity())<1.5)hSelPt->Fill(p.Perp());
	  //end of the filling  
     }
   
     //std::sort( phiIndex.begin(), phiIndex.end() );
     //phiIndex.erase( std::unique( phiIndex.begin(), phiIndex.end() ), phiIndex.end() );

 
     for(std::map< int, std::pair<StPicoMcTrack*,StPicoMcTrack*> >::iterator phi = phiMesonMcDaughters.begin(); phi != phiMesonMcDaughters.end(); phi++) 
     {
       int index = (phi->first);
       StPicoMcTrack *mcTrackKP = phiMesonMcDaughters[index].first;
       StPicoMcTrack *mcTrackKM = phiMesonMcDaughters[index].second;
       //cout << "We have " << nMCDaughters[index] << " MC Daughters " << endl;
       //cout << "We have " << nRCDaughters[index] << " RC Daughters " << endl;
       nMC[nMCDaughters[index]]++;
       nRC[nRCDaughters[index]]++;
       if(mcTrackKP == nullptr || mcTrackKM == nullptr) continue;
       if(nMCDaughters[index] != 2) continue;
       //cout << "We have MC tracks" << endl;

       double KP_P   = mcTrackKP->p().Mag();
       double KP_Pt  = mcTrackKP->p().Perp();
       double KP_Phi = mcTrackKP->p().Phi();
       double KP_Eta = mcTrackKP->p().PseudoRapidity();
       double KM_P   = mcTrackKM->p().Mag();
       double KM_Pt  = mcTrackKM->p().Perp();
       double KM_Phi = mcTrackKM->p().Phi();
       double KM_Eta = mcTrackKM->p().PseudoRapidity();


       TLorentzVector l_kp = mcTrackKP->fourMomentum();
       TLorentzVector l_km = mcTrackKM->fourMomentum();
      
       //cout << "K+ Mass = " << l_kp.M() << endl;
       //cout << "K- Mass = " << l_km.M() << endl << endl;

       //l_kp.SetPtEtaPhiM(static_cast<Double_t>(KP_Pt),static_cast<Double_t>(KP_Eta),static_cast<Double_t>(KP_Phi),static_cast<Double_t>(vmsa::mMassKaon));
       //l_km.SetPtEtaPhiM(static_cast<Double_t>(KM_Pt),static_cast<Double_t>(KM_Eta),static_cast<Double_t>(KM_Phi),static_cast<Double_t>(vmsa::mMassKaon));
       TLorentzVector lphi = l_kp + l_km;        
       double KP_Y   = l_kp.Rapidity();
       double KM_Y   = l_km.Rapidity();
       double KP_M2  = l_kp.M2();
       double KM_M2  = l_km.M2();

       TVector3 vMcPhiBeta = -1.0*lphi.BoostVector();
       TLorentzVector lMcKP;
       //lMcKP.SetPtEtaPhiM(KP_Pt,KP_Eta,KP_Phi,vmsa::mMassKaon);
       lMcKP = l_kp;//.SetPtEtaPhiM(KP_Pt,KP_Eta,KP_Phi,vmsa::mMassKaon);
       lMcKP.Boost(vMcPhiBeta);
       TVector3 vMcKP = lMcKP.Vect().Unit(); // direction of K+ momentum in phi-meson rest frame
       double phistar = vMcKP.Phi();

       TVector3 phiMomentumLabUnit = lphi.Vect().Unit();
        
       TLorentzVector lBeamPos;
       lBeamPos.SetPxPyPzE(0.0,0.0,9.75,9.796);
       lBeamPos.Boost(vMcPhiBeta);
       TLorentzVector lBeamNeg;
       lBeamNeg.SetPxPyPzE(0.0,0.0,-9.75,9.796);
       lBeamNeg.Boost(vMcPhiBeta);
       
       TLorentzVector lBeamTot = lBeamPos + lBeamNeg;
  
       TVector3 vBeamPos = lBeamPos.Vect();
       TVector3 vBeamNeg = lBeamNeg.Vect();

       TVector3 yaxisfrombeam = vBeamPos.Cross(vBeamNeg).Unit();
   
       TVector3 vBeamTotalUnit = lBeamTot.Vect().Unit();

       TVector3 zaxisfrombeam = vBeamTotalUnit;
       TVector3 xaxisfrombeam = yaxisfrombeam.Cross(zaxisfrombeam).Unit();

       double costheta = vMcKP.Dot(zaxisfrombeam); 

       double xprojection = vMcKP.Dot(xaxisfrombeam);
       double yprojection = vMcKP.Dot(yaxisfrombeam);
  
       double helicityangle = TMath::ATan2(yprojection,xprojection);
       if(helicityangle < 0.0) helicityangle += 2.0*TMath::Pi();
  
       //TVector3 QVector(TMath::Sin(Psi),-1.0*TMath::Cos(Psi),0.0);
       TVector3 QVectorMc(TMath::Sin(ep_full),-1.0*TMath::Cos(ep_full),0.0);
       TVector3 mcxprime(TMath::Cos(ep_full),TMath::Sin(ep_full),0.0);
       TVector3 mczprime(0.0,0.0,1.0); // beam direction
       TVector3 mcycalc = mczprime.Cross(phiMomentumLabUnit);
       TVector3 mcyhelicity = mcycalc.Unit();
       TVector3 mcxcalc = mcyhelicity.Cross(mczprime);
       TVector3 mcxhelicity = mcxcalc.Unit();
       

       Double_t mcproj_xprime = vMcKP.Dot(mcxprime);
       Double_t mcproj_zprime = vMcKP.Dot(mczprime);
       Double_t mczxangle = TMath::ATan2(TMath::Abs(mcproj_xprime),TMath::Abs(mcproj_zprime));
       Float_t mcphiprime = 0.0;
       if(mcproj_zprime > 0.0)
       {
         if(mcproj_xprime > 0.0)  mcphiprime = 2.0*TMath::Pi()-mczxangle; 
         if(mcproj_xprime < 0.0)  mcphiprime = mczxangle;
         if(mcproj_xprime == 0.0) mcphiprime = 0.0; 
       }               
       if(mcproj_zprime < 0.0)
       {
         if(mcproj_xprime > 0.0)  mcphiprime = TMath::Pi()+mczxangle;
         if(mcproj_xprime < 0.0)  mcphiprime = TMath::Pi()-mczxangle;
         if(mcproj_xprime == 0.0) mcphiprime = TMath::Pi();
       }
       if(mcproj_zprime == 0.0)
       {
         if(mcproj_xprime > 0.0)  mcphiprime = 3.0*TMath::Pi()/2.0;
         if(mcproj_xprime < 0.0)  mcphiprime = TMath::Pi()/2.0;
         if(mcproj_xprime == 0.0) mcphiprime = 0.0;
       }

       Double_t mcproj_xhelicity = vMcKP.Dot(mcxhelicity);
       Double_t mcproj_yhelicity = vMcKP.Dot(mcyhelicity);
       Double_t mcxyangle_helicity = TMath::ATan2(mcproj_yhelicity,mcproj_xhelicity);
       Double_t mcphihelicity = 0.0;
       if(mcproj_yhelicity > 0) mcphihelicity = mcxyangle_helicity - TMath::Pi();
       if(mcproj_yhelicity < 0) mcphihelicity = mcxyangle_helicity + TMath::Pi();
       if(mcphihelicity < 0) mcphihelicity += TMath::Pi()*2.0;

       TVector3 nQMc = QVectorMc.Unit(); // direction of QVector
       float McCosThetaStar = vMcKP.Dot(nQMc);
       //cout << "index = " << index << " mc phi candidate" << endl;
    
       //if(TMath::Abs(lphi.Rapidity()) > 1.0) continue;
       //if(lphi.Pt() < 1.2 || lphi.Pt() >= 4.2) continue;

       //int ptbin = -1; 
       //for(int ipt = 2; ipt < 6; ipt++)
       //{
       //  if(lphi.Pt() >= vmsa::pt_low[mEnergy][ipt] && lphi.Pt() < vmsa::pt_up[mEnergy][ipt]) 
       //  {
       //    ptbin = ipt; 
       //    break;
       //  }
       //}
      
       //if(!Sampling(f_mRhoPt[ptbin],TMath::Abs(McCosThetaStar),mMaxData)) continue;
       //if(!Sampling(f_mRhoPt[ptbin],McCosThetaStar,mMaxData[ptbin])) continue;
       //if(!Sampling2D(f_mRhoPt_2D[ptbin],McCosThetaStar,mcphiprime,mMax)) continue;
       //if(!Sampling2D(f_mRhoPt_Helicity_2D[ptbin],costheta,helicityangle,mMaxHelicity2D[ptbin])) continue;
       //cout << "ptbin = " << ptbin << ", BeforeSampling" << endl;
       //if(!Sampling2D(f_mRhoPt_Helicity_2D[ptbin],costheta,helicityangle,mMaxHelicity2D[ptbin])) continue;
       //cout << "After Sampling" << endl;
       //if(!Sampling2D(f_mRhoY_Helicity_2D[ybin] ,costheta,helicityangle,mMaxHelicity2DY[ybin])) continue;
       //if(!SamplingHelicity(f_mRhoPt_Helicity[ptbin],costheta,mMaxHelicity[ptbin])) continue;

       //cout << "index = " << index << " mc phi candidate passed rapidity cut" << endl;

       //cout << "Before filling MC hist " << endl;
       //mEffHistManger->FillHistMc(cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,phistar,McCosThetaStar,mcphiprime,costheta,mcphihelicity);

       StPicoMcTrack *mcMesonTrack = phiMeson[index];
       TLorentzVector l_meson = mcMesonTrack->fourMomentum();
       //l_meson->Print();
       //mEffHistManger->FillHistRc(reweight,cent9,l_meson.Pt(),l_meson.Eta(),l_meson.Rapidity(),l_meson.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,KP_M2,KM_M2,0.0,0.0,0);

       //if(!(KP_Pt > vmsa::mGlobPtMin && KP_P < vmsa::mPrimMomMax)) continue; 
       //if(!(KM_Pt > vmsa::mGlobPtMin && KM_P < vmsa::mPrimMomMax)) continue; 
       //if(fabs(KP_Eta) > vmsa::mEtaMax) continue;
       //if(fabs(KM_Eta) > vmsa::mEtaMax) continue;
       //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,KP_M2,KM_M2,0.0,0.0,1);
       //cout << "After filling MC hist " << endl;

       StPicoTrack *ptrackKP = phiMesonRcDaughters[index].first;
       StPicoTrack *ptrackKM = phiMesonRcDaughters[index].second;
       if(ptrackKP == nullptr || ptrackKM == nullptr || nRCDaughters[index] != 2)
       {
         mCent = cent9;
         mEpFull = ep_full;
         mRefMult = refMult;
         mNToFMatch = nBTofMatched;
         mVx = vx;
         mVy = vy;
         mVz = vz;
         mWeight = reweight;

         mLMeson = l_meson;
         mLKp = l_kp;
         mLKm = l_km;
         mRcExists = false;
         mLRcKp = TLorentzVector(-999,-999,-999,-999);
         mLRcKm = TLorentzVector(-999,-999,-999,-999);

         mNHitsFitP = -999;
         mNHitsMaxP = -999;
         mDEdxP = -999;
         mDcaP  = -999;
         mNHitsFitM = -999;
         mNHitsMaxM = -999;
         mDEdxM = -999;
         mDcaM  = -999;

         mMcPhiStar = phistar;
         mMcCosThetaStar = McCosThetaStar;
         mMcPhiPrime = mcphiprime;
         mMcCosTheta = costheta;
         mMcHelicityAngle = helicityangle;
         
         mRcPhiStar = -999;
         mRcCosThetaStar = -999;
         mRcPhiPrime = -999;
         mRcCosTheta = -999;
         mRcHelicityAngle = -999;
         
         mPassTpcKp = false;
         mPassTofKp = false;
         mPassM2Kp = false;
         mPassNsigKp = false; 

         mPassTpcKm = false;
         mPassTofKm = false;
         mPassM2Km = false;
         mPassNsigKm = false; 

         mNSigmaKaonP = -999.;
         mNSigmaKaonM = -999.;
  
  
         //cout << "There is a nullptr and we are going to fill the TTree" << endl;       
         mKaonTree->Fill();
         //cout << "There is a nullptr and we FILLED the TTree" << endl;       

         continue;
       }

       //cout << "We have RC tracks" << endl;

       double rKP_P   = ptrackKP->pMom().Mag();
       double rKP_Pt  = ptrackKP->pMom().Perp();
       double rKP_Phi = ptrackKP->pMom().Phi();
       double rKP_Eta = ptrackKP->pMom().PseudoRapidity();
       int    rKP_nhits = ptrackKP->nHitsFit();
       int    rKP_nhitsmax = ptrackKP->nHitsMax();
       float  rKP_nhitsratio = (float)ptrackKP->nHitsFit()/(float)ptrackKP->nHitsMax();
       float  rKP_dca = ptrackKP->gDCA(vx,vy,vz);
       float  rKP_dedx = ptrackKP->dEdx();
       double rKM_P   = ptrackKM->pMom().Mag();
       double rKM_Pt  = ptrackKM->pMom().Perp();
       double rKM_Phi = ptrackKM->pMom().Phi();
       double rKM_Eta = ptrackKM->pMom().PseudoRapidity();
       int    rKM_nhits = ptrackKM->nHitsFit();
       int    rKM_nhitsmax = ptrackKM->nHitsMax();
       float  rKM_nhitsratio = (float)ptrackKM->nHitsFit()/(float)ptrackKM->nHitsMax();
       float  rKM_dca = ptrackKM->gDCA(vx,vy,vz);
       float  rKM_dedx = ptrackKM->dEdx();

       //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,KP_M2,KM_M2,0.0,0.0,2);

       //cout << "index = " << index << " rc phi candidate" << endl;


       //cout << "index = " << index << " rc phi pass" << endl;

       TLorentzVector rl_kp;
       TLorentzVector rl_km;
       rl_kp.SetPtEtaPhiM(static_cast<Double_t>(rKP_Pt),static_cast<Double_t>(rKP_Eta),static_cast<Double_t>(rKP_Phi),static_cast<Double_t>(vmsa::mMassKaon));
       rl_km.SetPtEtaPhiM(static_cast<Double_t>(rKM_Pt),static_cast<Double_t>(rKM_Eta),static_cast<Double_t>(rKM_Phi),static_cast<Double_t>(vmsa::mMassKaon));
       TLorentzVector rlphi = rl_kp + rl_km;
       double rKP_Y   = rl_kp.Rapidity();
       double rKM_Y   = rl_km.Rapidity();

       TVector3 vRcPhiBeta = -1.0*rlphi.BoostVector();
       TLorentzVector lRcKP;
       //lRcKP.SetPtEtaPhiM(rKP_Pt,rKP_Eta,rKP_Phi,vmsa::mMassKaon);
       lRcKP = rl_kp;//.SetPtEtaPhiM(rKP_Pt,rKP_Eta,rKP_Phi,vmsa::mMassKaon);
       lRcKP.Boost(vRcPhiBeta);
       TVector3 vRcKP = lRcKP.Vect().Unit(); // direction of K+ momentum in phi-meson rest frame
       double phistarRc = vRcKP.Phi();
       TVector3 phiMomentumLabUnitRc = rlphi.Vect().Unit();

       TLorentzVector lBeamPosRc;
       lBeamPosRc.SetPxPyPzE(0.0,0.0,9.75,9.796);
       lBeamPosRc.Boost(vRcPhiBeta);
       TLorentzVector lBeamNegRc;
       lBeamNegRc.SetPxPyPzE(0.0,0.0,-9.75,9.796);
       lBeamNegRc.Boost(vRcPhiBeta);
       
       TLorentzVector lBeamTotRc = lBeamPosRc + lBeamNegRc;
  
       TVector3 vBeamPosRc = lBeamPosRc.Vect();
       TVector3 vBeamNegRc = lBeamNegRc.Vect();

       TVector3 yaxisfrombeamRc = vBeamPosRc.Cross(vBeamNegRc).Unit();
   
       TVector3 vBeamTotalUnitRc = lBeamTotRc.Vect().Unit();
       
       TVector3 zaxisfrombeamRc = vBeamTotalUnitRc;
       TVector3 xaxisfrombeamRc = yaxisfrombeamRc.Cross(zaxisfrombeamRc).Unit();

       double costhetaRc = vRcKP.Dot(zaxisfrombeamRc); 

       double xprojectionRc = vRcKP.Dot(xaxisfrombeamRc);
       double yprojectionRc = vRcKP.Dot(yaxisfrombeamRc);
  
       double helicityangleRc = TMath::ATan2(yprojectionRc,xprojectionRc);
       if(helicityangleRc < 0.0) helicityangleRc += 2.0*TMath::Pi();

       TVector3 QVectorRc(TMath::Sin(ep_full),-1.0*TMath::Cos(ep_full),0.0);
       TVector3 xprime(TMath::Cos(ep_full),TMath::Sin(ep_full),0.0);
       TVector3 zprime(0.0,0.0,1.0);
       TVector3 ycalc = zprime.Cross(phiMomentumLabUnitRc);
       TVector3 yhelicity = ycalc.Unit();
       TVector3 xcalc = yhelicity.Cross(zprime);
       TVector3 xhelicity = xcalc.Unit();

       Double_t proj_xprime = vRcKP.Dot(xprime);
       Double_t proj_zprime = vRcKP.Dot(zprime);
       Double_t zxangle = TMath::ATan2(TMath::Abs(proj_zprime),TMath::Abs(proj_xprime));
       Float_t phiprimeRc = 0.0;
       if(proj_zprime > 0.0)
       {
         if(proj_xprime > 0.0)  phiprimeRc = 2.0*TMath::Pi()-zxangle; 
         if(proj_xprime < 0.0)  phiprimeRc = zxangle;
         if(proj_xprime == 0.0) phiprimeRc = 0.0; 
       }               
       if(proj_zprime < 0.0)
       {
         if(proj_xprime > 0.0)  phiprimeRc = TMath::Pi()+zxangle;
         if(proj_xprime < 0.0)  phiprimeRc = TMath::Pi()-zxangle;
         if(proj_xprime == 0.0) phiprimeRc = TMath::Pi();
       }
       if(proj_zprime == 0.0)
       {
         if(proj_xprime > 0.0)  phiprimeRc = 3.0*TMath::Pi()/2.0;
         if(proj_xprime < 0.0)  phiprimeRc = TMath::Pi()/2.0;
         if(proj_xprime == 0.0) phiprimeRc = 0.0;
       }

       Double_t proj_xhelicity = vRcKP.Dot(xhelicity);
       Double_t proj_yhelicity = vRcKP.Dot(yhelicity);
       Double_t xyangle_helicity = TMath::ATan2(proj_yhelicity,proj_xhelicity);
       Double_t phihelicity = 0.0;
       if(proj_yhelicity > 0) phihelicity = xyangle_helicity - TMath::Pi();
       if(proj_yhelicity < 0) phihelicity = xyangle_helicity + TMath::Pi();
       if(phihelicity < 0) phihelicity += TMath::Pi()*2.0;

       TVector3 nQRc = QVectorRc.Unit(); // direction of QVector
       float RcCosThetaStar = vRcKP.Dot(nQRc);

       float rkp_ns = ptrackKP->nSigmaKaon();
       float rkm_ns = ptrackKM->nSigmaKaon();


       //mEffHistManger->FillHistRc(reweight,cent9,rlphi.Pt(),rlphi.Eta(),rlphi.Rapidity(),rlphi.Phi(),ep_full,rKP_Pt,rKP_Eta,rKP_Y,rKP_Phi,rKM_Pt,rKM_Eta,rKM_Y,rKM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistarRc,RcCosThetaStar,phiprimeRc,costhetaRc,helicityangleRc,KP_M2,KM_M2,rkp_ns,rkm_ns,3);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 0, rKP_dca, rKP_nhits, rKP_nhitsratio, rKP_P, rKP_dedx, 3);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 1, rKM_dca, rKM_nhits, rKM_nhitsratio, rKM_P, rKM_dedx, 3);
       //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,KP_M2,KM_M2,0.0,0.0,2);
       //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,3);

       //mEffHistManger->FillHistRc(reweight,cent9,rlphi.Pt(),rlphi.Eta(),rlphi.Rapidity(),rlphi.Phi(),ep_full,rKP_Pt,rKP_Eta,rKP_Y,rKP_Phi,rKM_Pt,rKM_Eta,rKM_Y,rKM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistarRc,RcCosThetaStar,phiprimeRc,costhetaRc,helicityangleRc,KP_M2,KM_M2,rkp_ns,rkm_ns,4);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 0, rKP_dca, rKP_nhits, rKP_nhitsratio, rKP_P, rKP_dedx, 4);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 1, rKM_dca, rKM_nhits, rKM_nhitsratio, rKM_P, rKM_dedx, 4);
//       mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,KP_M2,KM_M2,0.0,0.0,3);
       //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,4);
       //mEffHistManger->FillHistRc(reweight,cent9,rlphi.Pt(),rlphi.Eta(),rlphi.Rapidity(),rlphi.Phi(),ep_full,rKP_Pt,rKP_Eta,rKP_Y,rKP_Phi,rKM_Pt,rKM_Eta,rKM_Y,rKM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistarRc,RcCosThetaStar,phiprimeRc,costhetaRc,helicityangleRc,KP_M2,KM_M2,rkp_ns,rkm_ns,5);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 0, rKP_dca, rKP_nhits, rKP_nhitsratio, rKP_P, rKP_dedx, 5);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 1, rKM_dca, rKM_nhits, rKM_nhitsratio, rKM_P, rKM_dedx, 5);
//       mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,KP_M2,KM_M2,0.0,0.0,4);

       //if( TMath::Abs(rKP_Eta) > 1.0 || TMath::Abs(rKM_Eta) > 1.0) continue; 
       //mEffHistManger->FillHistRc(reweight,cent9,rlphi.Pt(),rlphi.Eta(),rlphi.Rapidity(),rlphi.Phi(),ep_full,rKP_Pt,rKP_Eta,rKP_Y,rKP_Phi,rKM_Pt,rKM_Eta,rKM_Y,rKM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistarRc,RcCosThetaStar,phiprimeRc,costhetaRc,helicityangleRc,KP_M2,KM_M2,rkp_ns,rkm_ns,4);

       //cout << "Before TPC" << endl;
       Bool_t passTpcKp = true;
       Bool_t passTpcKm = true;
       //if(rlphi.Pt() < 1.2 || rlphi.Pt() >= 4.2) continue;
       //if(TMath::Abs(rlphi.Rapidity()) > 1.0) continue;
       //if(fabs(rKP_Eta) > vmsa::mEtaMax) continue;
       //if(fabs(rKM_Eta) > vmsa::mEtaMax) continue;
       //if(!(rKP_Pt > vmsa::mGlobPtMin && rKP_P < vmsa::mPrimMomMax)) continue; 
       //if(!(rKM_Pt > vmsa::mGlobPtMin && rKM_P < vmsa::mPrimMomMax)) continue; 
       //if( !mVecMesonCut->passTrackMeson(ptrackKP, mPicoEvent, mPid) ) continue;
       //if( !mVecMesonCut->passTrackMeson(ptrackKM, mPicoEvent, mPid) ) continue;

       if( !mVecMesonCut->passTrackMeson(ptrackKP, mPicoEvent, mPid) ) passTpcKp = false;
       if( !mVecMesonCut->passTrackMeson(ptrackKM, mPicoEvent, mPid) ) passTpcKm = false;
       //if( (rlphi.Pt() < 1.2 || rlphi.Pt() >= 4.2)
       //|| TMath::Abs(rlphi.Rapidity()) > 1.0
       //|| fabs(rKP_Eta) > vmsa::mEtaMax
       //|| fabs(rKM_Eta) > vmsa::mEtaMax
       //|| !(rKP_Pt > vmsa::mGlobPtMin && rKP_P < vmsa::mPrimMomMax)
       //|| !(rKM_Pt > vmsa::mGlobPtMin && rKM_P < vmsa::mPrimMomMax)
       //||  !mVecMesonCut->passTrackMeson(ptrackKP, mPicoEvent, mPid) 
       //||  !mVecMesonCut->passTrackMeson(ptrackKM, mPicoEvent, mPid) ) 
       //{
       //  //mVecMesonCut->getPrimaryMass2(ptrackKP, mPicoDst);
       //  //mVecMesonCut->getPrimaryMass2(ptrackKM, mPicoDst);
       //  passTpc = false;
       //}
       //cout << "After TPC" << endl;

       //if(fabs(rkm_ns) > 2.5 || fabs(rkp_ns) > 2.5) continue;
       //mEffHistManger->FillHistRc(reweight,cent9,rlphi.Pt(),rlphi.Eta(),rlphi.Rapidity(),rlphi.Phi(),ep_full,rKP_Pt,rKP_Eta,rKP_Y,rKP_Phi,rKM_Pt,rKM_Eta,rKM_Y,rKM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistarRc,RcCosThetaStar,phiprimeRc,costhetaRc,helicityangleRc,KP_M2,KM_M2,rkp_ns,rkm_ns,6);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 0, rKP_dca, rKP_nhits, rKP_nhitsratio, rKP_P, rKP_dedx, 6);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 1, rKM_dca, rKM_nhits, rKM_nhitsratio, rKM_P, rKM_dedx, 6);
       //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,KP_M2,KM_M2,0.0,0.0,5);
       //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,5);
 
       //cout << "Before TOF Matching" << endl;
       //if( !tofReconstructed(0,cent9,ep_full,rl_kp,phiswitch) ) continue;
       //if( !tofReconstructed(1,cent9,ep_full,rl_km,phiswitch) ) continue;
       //if( !tofReconstructed(0,9,ep_full,rl_kp,phiswitch) ) continue;
       //if( !tofReconstructed(1,9,ep_full,rl_km,phiswitch) ) continue;
       
       //cout << "Before TOF" << endl;
       Bool_t passTofKp = true;
       Bool_t passTofKm = true;
       if( !tofReconstructed(0,9,ep_full,rl_kp,phiswitch) ) passTofKp = false;
       if( !tofReconstructed(1,9,ep_full,rl_km,phiswitch) ) passTofKm = false;
       //cout << "After TOF" << endl;
       //cout << "After TOF Matching" << endl;
       //mEffHistManger->FillHistRc(reweight,cent9,rlphi.Pt(),rlphi.Eta(),rlphi.Rapidity(),rlphi.Phi(),ep_full,rKP_Pt,rKP_Eta,rKP_Y,rKP_Phi,rKM_Pt,rKM_Eta,rKM_Y,rKM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistarRc,RcCosThetaStar,phiprimeRc,costhetaRc,helicityangleRc,KP_M2,KM_M2,rkp_ns,rkm_ns,7);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 0, rKP_dca, rKP_nhits, rKP_nhitsratio, rKP_P, rKP_dedx, 7);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 1, rKM_dca, rKM_nhits, rKM_nhitsratio, rKM_P, rKM_dedx, 7);
       //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,KP_M2,KM_M2,0.0,0.0,6);
       //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,6);

       //cout << "index = " << index << " mc phi candidate passed rapidity cut" << endl;

       //cout << "Before PID cuts " << endl;

       //int EtaBinM = -1;
       //int PhiBinM = -1;
       //findFuncPID(rl_km,0,EtaBinM,PhiBinM);
       //int EtaBinP = -1;
       //int PhiBinP = -1;
       //findFuncPID(rl_kp,1,EtaBinP,PhiBinP);

       int GlobalBinM = -1;
       int GlobalBinP = -1;
       findm2HistPID(rl_km,0,GlobalBinM);
       findm2HistPID(rl_kp,1,GlobalBinP);
  
       //cout << "EtaBinM = " << EtaBinM << ", PhiBinM = " << PhiBinM << endl;
       //cout << "EtaBinP = " << EtaBinP << ", PhiBinP = " << PhiBinP << endl;
       //if( !pass_m2_PID(0,9,rl_km,EtaBinM,PhiBinM) ) continue;
       //if( !pass_m2_PID(1,9,rl_kp,EtaBinP,PhiBinP) ) continue;
       //cout << "Before M2" << endl;
       Bool_t passm2kp = true;
       Bool_t passm2km = true;
     
       //if(EtaBinM < 0 || EtaBinP < 0 || PhiBinM < 0 || PhiBinM < 0) passm2 = false;
       //else
       //{
       //  if( !pass_m2_PID(0,9,rl_km,EtaBinM,PhiBinM)
       //  ||  !pass_m2_PID(1,9,rl_kp,EtaBinP,PhiBinP) ) passm2 = false; 
       //}
       if(GlobalBinM < 0 || GlobalBinP < 0) 
       {
         passm2kp = false;
         passm2km = false;
       }
       else
       {
         if( !passhist_m2_PID(0,9,rl_km,GlobalBinM) ) passm2km = false;
         if( !passhist_m2_PID(1,9,rl_kp,GlobalBinP) ) passm2kp = false; 
       }
       //cout << "After M2" << endl;
       //cout << "K- passed m2 PID cut" << endl;
       //cout << "K- passed nsig PID cut" << endl;
 
       //mEffHistManger->FillHistRc(reweight,cent9,rlphi.Pt(),rlphi.Eta(),rlphi.Rapidity(),rlphi.Phi(),ep_full,rKP_Pt,rKP_Eta,rKP_Y,rKP_Phi,rKM_Pt,rKM_Eta,rKM_Y,rKM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistarRc,RcCosThetaStar,phiprimeRc,costhetaRc,helicityangleRc,KP_M2,KM_M2,rkp_ns,rkm_ns,8);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 0, rKP_dca, rKP_nhits, rKP_nhitsratio, rKP_P, rKP_dedx, 8);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 1, rKM_dca, rKM_nhits, rKM_nhitsratio, rKM_P, rKM_dedx, 8);
       //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,KP_M2,KM_M2,0.0,0.0,7);
       //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,7);

       //cout << "EtaBin = " << EtaBin << ", PhiBin = " << PhiBin << endl;
       //cout << "K+ passed m2 PID cut" << endl;
       //if( !pass_nsig_PID(0,9,rl_km,EtaBinM,PhiBinM) ) continue;
       //if( !pass_nsig_PID(1,9,rl_kp,EtaBinP,PhiBinP) ) continue;
       //cout << "Before Nsig" << endl;
       Bool_t passnskp = true;
       Bool_t passnskm = true;
       //if(EtaBinM < 0 || EtaBinP < 0 || PhiBinM < 0 || PhiBinM < 0) passm2 = false;
       //else
       //{
       //  if( !pass_nsig_PID(0,9,rl_km,EtaBinM,PhiBinM)
       //  ||  !pass_nsig_PID(1,9,rl_kp,EtaBinP,PhiBinP) ) passns = false;
       //}
       GlobalBinM = -1;
       GlobalBinP = -1;
       findnsigHistPID(rl_km,0,GlobalBinM);
       findnsigHistPID(rl_kp,1,GlobalBinP);
       if(GlobalBinM < 0 || GlobalBinP < 0) 
       {
         passnskp = false;
         passnskm = false;
       }
       else
       {
         if( !passhist_nsig_PID(0,9,rl_km,GlobalBinM) ) passnskm = false;
         if( !passhist_nsig_PID(1,9,rl_kp,GlobalBinP) ) passnskp = false;
       }
       //cout << "K+ passed nsig PID cut" << endl;
       //cout << "After PID cuts " << endl;

       //cout << "Before filling RC hist " << endl;
       //mEffHistManger->FillHistRc(reweight,cent9,rlphi.Pt(),rlphi.Eta(),rlphi.Rapidity(),rlphi.Phi(),ep_full,rKP_Pt,rKP_Eta,rKP_Y,rKP_Phi,rKM_Pt,rKM_Eta,rKM_Y,rKM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistarRc,RcCosThetaStar,phiprimeRc,costhetaRc,helicityangleRc,KP_M2,KM_M2,rkp_ns,rkm_ns,9);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 0, rKP_dca, rKP_nhits, rKP_nhitsratio, rKP_P, rKP_dedx, 9);
       //mEffHistManger->FillTrackHistQA(reweight, cent9, 1, rKM_dca, rKM_nhits, rKM_nhitsratio, rKM_P, rKM_dedx, 9);
//       mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,helicityangle,KP_M2,KM_M2,0.0,0.0,8);
       //mEffHistManger->FillHistRc(reweight,cent9,lphi.Pt(),lphi.Eta(),lphi.Rapidity(),lphi.Phi(),ep_full,KP_Pt,KP_Eta,KP_Y,KP_Phi,KM_Pt,KM_Eta,KM_Y,KM_Phi,l_meson.Pt(),l_meson.Phi(),l_meson.Rapidity(),phistar,McCosThetaStar,mcphiprime,costheta,mcphihelicity,8);
       //cout << "After filling RC hist " << endl;

       //cout << "azimuthal angle = " << lphi.Phi() << endl;
       
       mCent = cent9;
       mEpFull = ep_full;
       mRefMult = refMult;
       mNToFMatch = nBTofMatched;
       mVx = vx;
       mVy = vy;
       mVz = vz;
       mWeight = reweight;

       mLMeson = l_meson;
       mLKp = l_kp;
       mLKm = l_km;
       mRcExists = true;
       mLRcKp = rl_kp;
       mLRcKm = rl_km;

       mNHitsFitP = rKP_nhits;
       mNHitsMaxP = rKP_nhitsmax;
       mDEdxP = rKP_dedx;
       mDcaP  = rKP_dca;
       mNHitsFitM = rKM_nhits;
       mNHitsMaxM = rKM_nhitsmax;
       mDEdxM = rKM_dedx;
       mDcaM  = rKM_dca;

       mMcPhiStar = phistar;
       mMcCosThetaStar = McCosThetaStar;
       mMcPhiPrime = mcphiprime;
       mMcCosTheta = costheta;
       mMcHelicityAngle = helicityangle;
       
       mRcPhiStar = phistarRc;
       mRcCosThetaStar = RcCosThetaStar;
       mRcPhiPrime = phiprimeRc;
       mRcCosTheta = costhetaRc;
       mRcHelicityAngle = helicityangleRc;
       
       mPassTpcKp = passTpcKp;
       mPassTofKp = passTofKp;
       mPassM2Kp = passm2kp;
       mPassNsigKp = passnskp;

       mPassTpcKm = passTpcKm;
       mPassTofKm = passTofKm;
       mPassM2Km = passm2km;
       mPassNsigKm = passnskm;
  
       mNSigmaKaonP = ptrackKP->nSigmaKaon();
       mNSigmaKaonM = ptrackKM->nSigmaKaon();

 
       //cout << "There is RC information and we are going to fill the TTree" << endl;       
       mKaonTree->Fill();
       //cout << "There is RC information we FILLED the TTree" << endl;       

     }


     phiMesonMcDaughters.clear();
     phiMesonRcDaughters.clear();
     phiMeson.clear();
     nMCDaughters.clear();
     nRCDaughters.clear();
  }
  for(int i = 0; i < 5; i++)
  {
    cout << nMC[i] << " mesons with " << i << " MC daughters " << endl;
  } 

  cout << endl;

  for(int i = 0; i < 5; i++)
  {
    cout << nRC[i] << " mesons with " << i << " RC daughters " << endl;
  }
  //mEffHistManger->CalEfficiency();
  //mEffHistManger->CalEffPtEtaPhi();

  tags_output->cd();
  //mEffHistManger->WriteHist();
  //mEffHistManger->WriteHistQA();
  mKaonTree->Write();

  //hEffPt->Divide(hSelPt,hSelPtMc,1,1,"B");

  if (nEvents > 1) chain -> Finish() ;

  if(tags_output!=NULL) tags_output -> Write() ;
  if(tags_output!=NULL) tags_output -> Close() ;
  //flush(tags_output);
  delete tags_output;

  // Cleanup
  delete chain ;
}

