#include <iostream>
#include <string> 
#include <map>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TPythia6.h"
#include "TPythia6Decayer.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TNtuple.h"
#include "TF1.h"
// #include "/global/homes/x/xusun/STAR/VecMesonSpinAlignment/Utility/StSpinAlignmentCons.h"
// #include "/global/homes/x/xusun/STAR/VecMesonSpinAlignment/Utility/functions.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/functions.h"

using namespace std;

typedef std::map<std::string,TH1D*> TH1DMap;
typedef std::map<std::string,TF1*> TF1Map;

void readEfficiency(int energy, int year, int cut, int jobID);
void readTofEff(int energy);
void readTofEffFit(int energy);
void getKinematics(TLorentzVector& lPhi, double const mass);
void setDecayChannels(int const pid);
void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters);
void fill(int const kf, TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus);
bool tpcReconstructed(int iParticleIndex, int cent, TLorentzVector const& lKaon);
void findHist(TLorentzVector const& lKaon, int iParticleIndex, int& EtaBin, int& PhiBin); // iParticleIndex = 0 => K+
void findHist_ToF(TLorentzVector const& lKaon, int iParticleIndex, int& EtaBin, int& PhiBin); // iParticleIndex = 0 => K+
void write();

TPythia6Decayer* pydecay;
TNtuple* McPhiMeson;

TH1DMap h_EffKplus;
TH1DMap h_EffKminus;
TH1D *h_FrameEta[2];
TH1D *h_FramePhi[2];

TH1DMap h_TofKplus;
TH1DMap h_TofKminus;
TH1D *h_FrameEta_ToF[2];
TH1D *h_FramePhi_ToF[2];

TF1Map f_TofKplus;
TF1Map f_TofKminus;

TFile *File_OutPut;

void toyMcPhiDecay(const int energy = 4, const int pid = 0, const int year = 1, const int cut = 0, const int NMax = 50, const int jobID = 9)
{
  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();
  gRandom->SetSeed();
  readEfficiency(energy,year,cut,jobID);
  readTofEff(energy);
  readTofEffFit(energy);

  pydecay = TPythia6Decayer::Instance();
  pydecay->Init();
  setDecayChannels(pid); // phi--> K+K-

  TClonesArray ptl("TParticle", 10);
  TLorentzVector *lPhi = new TLorentzVector();
  for(int i_ran = 0; i_ran < NMax; ++i_ran)
  {
    if (floor(10.0*i_ran/ static_cast<float>(NMax)) > floor(10.0*(i_ran-1)/ static_cast<float>(NMax)))
    cout << "=> processing data: " << 100.0*i_ran/ static_cast<float>(NMax) << "%" << endl;

    getKinematics(*lPhi,vmsa::InvMass[pid]);
    // if (fabs(lPhi->Phi()) > TMath::Pi()) continue;
    decayAndFill(vmsa::decayMother[pid],lPhi,ptl);

    if (i_ran % 1000 == 1) McPhiMeson->AutoSave("SaveSelf");
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  write();

  stopWatch->Stop();   
  stopWatch->Print();
}

void getKinematics(TLorentzVector& lPhi, double const mass)
{
   double const pt = gRandom->Uniform(vmsa::ptMin, vmsa::ptEffMax);
   // double const pt = gRandom->Uniform(0.8, 2.4);
   double const y = gRandom->Uniform(-vmsa::acceptanceRapidity, vmsa::acceptanceRapidity);
   double const phi = TMath::TwoPi() * gRandom->Rndm();

   double const mT = sqrt(mass * mass + pt * pt);
   double const pz = mT * sinh(y);
   double const E = mT * cosh(y);

   lPhi.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

void setDecayChannels(int const pid)
{
  int const mdme = vmsa::decayChannels[pid];
  cout << "mdme = " << mdme << endl;
  for (int idc = vmsa::decayChannelsFirst[pid]; idc < vmsa::decayChannelsSecond[pid] + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0); // close all decay channel
  TPythia6::Instance()->SetMDME(mdme, 1, 1); // open the one we need
  int *PYSeed = new int;
  TPythia6::Instance()->SetMRPY(1,(int)PYSeed); // Random seed
}

void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters)
{
  pydecay->Decay(kf, lPhi);
  pydecay->ImportParticles(&daughters);

  TLorentzVector lKplus;
  TLorentzVector lKminus;

  int nTrk = daughters.GetEntriesFast();
  for (int iTrk = 0; iTrk < nTrk; ++iTrk)
  {
    TParticle* ptl0 = (TParticle*)daughters.At(iTrk);

    switch (ptl0->GetPdgCode())
    {
      case 321:
	ptl0->Momentum(lKplus);
	break;
      case -321:
	ptl0->Momentum(lKminus);
	break;
      default:
	break;
    }
  }
  daughters.Clear("C");

  fill(kf,lPhi,lKplus,lKminus);
}

void fill(int const kf, TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus)
{
   int const centrality = floor(vmsa::NCentMax * gRandom->Rndm());
   // int const centrality = floor(2 * gRandom->Rndm());
   TLorentzVector lRcPhi = lKplus + lKminus; // phi meson reconstruction
   // cout << "lPhi.pt = " << lPhi->Pt() << ", lPhi.eta = " << lPhi->Eta() << ", lPhi.phi = " << lPhi->Phi() << ", lPhi.m = " << lPhi->M() << endl;
   // cout << "lRcPhi.pt = " << lRcPhi.Pt() << ", lRcPhi.eta = " << lRcPhi.Eta() << ", lRcPhi.phi = " << lRcPhi.Phi() << ", lRcPhi.m = " << lRcPhi.M() << endl;
   // cout << endl;

   float arr[110];
   int iArr = 0;
   arr[iArr++] = centrality; // McPhi
   arr[iArr++] = lPhi->Pt();
   arr[iArr++] = lPhi->P();
   arr[iArr++] = lPhi->PseudoRapidity();
   arr[iArr++] = lPhi->Rapidity();
   arr[iArr++] = lPhi->Phi();
   arr[iArr++] = lPhi->M();
   arr[iArr++] = kf;

   arr[iArr++] = lKplus.Pt();
   arr[iArr++] = lKplus.PseudoRapidity();
   arr[iArr++] = lKplus.Rapidity();
   arr[iArr++] = lKplus.Phi();
   arr[iArr++] = lKplus.M();
   arr[iArr++] = 321;

   arr[iArr++] = lKplus.Pt();
   arr[iArr++] = lKplus.PseudoRapidity();
   arr[iArr++] = lKplus.Rapidity();
   arr[iArr++] = lKplus.Phi();
   arr[iArr++] = lKplus.M();
   arr[iArr++] = tpcReconstructed(0,centrality,lKplus);

   arr[iArr++] = lKminus.Pt();
   arr[iArr++] = lKminus.PseudoRapidity();
   arr[iArr++] = lKminus.Rapidity();
   arr[iArr++] = lKminus.Phi();
   arr[iArr++] = lKminus.M();
   arr[iArr++] = -321;

   arr[iArr++] = lKminus.Pt();
   arr[iArr++] = lKminus.PseudoRapidity();
   arr[iArr++] = lKminus.Rapidity();
   arr[iArr++] = lKminus.Phi();
   arr[iArr++] = lKminus.M();
   arr[iArr++] = tpcReconstructed(1,centrality,lKminus);

   arr[iArr++] = lRcPhi.Pt();
   arr[iArr++] = lRcPhi.P();
   arr[iArr++] = lRcPhi.PseudoRapidity();
   arr[iArr++] = lRcPhi.Rapidity();
   arr[iArr++] = lRcPhi.Phi();
   arr[iArr++] = lRcPhi.M();

   McPhiMeson->Fill(arr);
   // if(lRcPhi.Pt() < 10e-4) cout << "lPhi->Pt = " << lPhi->Pt() << ", lRcPhi.Pt = " << lRcPhi.Pt() << endl;
}

bool tpcReconstructed(int iParticleIndex, int cent, TLorentzVector const& lKaon)
{
   if(fabs(lKaon.Eta()) >= 1.0) return false;

   TH1D *h_TPC = NULL;
   int EtaBin_TPC = -1;
   int PhiBin_TPC = -1;
   findHist(lKaon,iParticleIndex,EtaBin_TPC,PhiBin_TPC);

   TH1D *h_ToF = NULL;
   TF1 *f_ToF = NULL;
   int EtaBin_ToF = -1;
   int PhiBin_ToF = -1;
   findHist_ToF(lKaon,iParticleIndex,EtaBin_ToF,PhiBin_ToF);

   // cout << KEY.c_str() << endl;
   if (iParticleIndex == 0)
   {
     // string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
     string KEY_TPC = Form("h_mEff_Cent_9_Eta_%d_Phi_%d",EtaBin_TPC,PhiBin_TPC); // get TPC eff @ 20-60%
     h_TPC = h_EffKplus[KEY_TPC];

     // string KEY_ToF = Form("h_mEfficiency_Kplus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff
     string KEY_ToF = Form("h_mEfficiency_Kplus_Cent_9_Eta_%d_Phi_%d",EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
     h_ToF = h_TofKplus[KEY_ToF];

     // string KEY_ToFFit = Form("f_mToFMatch_Kplus_Cent_9_Eta_%d_Phi_%d",EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
     string KEY_ToFFit = Form("f_mToFMatch_Kplus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
     f_ToF = f_TofKplus[KEY_ToFFit]; // only 20-60%
   }
   else
   {
     // string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
     string KEY_TPC = Form("h_mEff_Cent_9_Eta_%d_Phi_%d",EtaBin_TPC,PhiBin_TPC); // get TPC eff @ 20-60%
     h_TPC = h_EffKminus[KEY_TPC];

     // string KEY_ToF = Form("h_mEfficiency_Kminus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff
     string KEY_ToF = Form("h_mEfficiency_Kminus_Cent_9_Eta_%d_Phi_%d",EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
     h_ToF = h_TofKminus[KEY_ToF];

     // string KEY_ToFFit = Form("f_mToFMatch_Kminus_Cent_9_Eta_%d_Phi_%d",EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
     string KEY_ToFFit = Form("f_mToFMatch_Kminus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
     f_ToF = f_TofKminus[KEY_ToFFit]; // only 20-60%
   }

   double pt = lKaon.Perp();
   int const bin_TPC = h_TPC->FindBin(pt);
   bool is_TPC = gRandom->Rndm() < h_TPC->GetBinContent(bin_TPC);

   int const bin_ToF = h_ToF->FindBin(pt); // tof fit and hist combined
   double prob_tof = 0.0;
   if(pt > 0.3) prob_tof = f_ToF->Eval(pt); // donot use fit extrapolation
   else prob_tof = h_ToF->GetBinContent(bin_ToF);
   bool is_ToF = gRandom->Rndm() < prob_tof;

   // int const bin_ToF = h_ToF->FindBin(pt); // tof hist only
   // bool is_ToF = gRandom->Rndm() < h_ToF->GetBinContent(bin_ToF);

   // bool is_ToF = gRandom->Rndm() < f_ToF->Eval(pt); // tof fit only

   // cout << "is_TPC: " << is_TPC << ", is_ToF: " << is_ToF << ", is_TPC && is_ToF: " << (is_TPC && is_ToF) << endl;
   // cout << "pt = " << pt << ", h_ToF = " << h_ToF->GetBinContent(bin_ToF) << ", f_ToF = " << f_ToF->Eval(pt) << ", prob_tof = " << prob_tof << endl; // comparison between eff from hist and func

   return is_TPC && is_ToF;
}

void findHist(TLorentzVector const& lKaon, int iParticleIndex, int& EtaBin, int& PhiBin)
{
  float eta = lKaon.Eta();
  EtaBin = h_FrameEta[iParticleIndex]->FindBin(eta)-1;
  // cout << "eta = " << eta << ", EtaBin = " << EtaBin << endl;
  float phi = lKaon.Phi();
  PhiBin = h_FramePhi[iParticleIndex]->FindBin(phi)-1;
}

void findHist_ToF(TLorentzVector const& lKaon, int iParticleIndex, int& EtaBin, int& PhiBin)
{
  float eta = lKaon.Eta();
  EtaBin = h_FrameEta_ToF[iParticleIndex]->FindBin(eta)-1;
  // cout << "eta = " << eta << ", EtaBin = " << EtaBin << endl;
  float phi = lKaon.Phi();
  PhiBin = h_FramePhi_ToF[iParticleIndex]->FindBin(phi)-1;
  // cout << "phi = " << phi << ", PhiBin = " << PhiBin << endl;
}

void readEfficiency(int energy, int year, int cut, int jobID)
{
  // string inputKplus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[0].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  string inputKplus = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_%s_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[0].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  TFile *File_Kplus = TFile::Open(inputKplus.c_str());
  cout << "OPEN Efficiency File for K+: " << inputKplus.c_str() << endl;

  // string inputKminus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[1].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  string inputKminus = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_%s_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[1].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  TFile *File_Kminus = TFile::Open(inputKminus.c_str());
  cout << "OPEN Efficiency File for K-: " << inputKminus.c_str() << endl;

  h_FrameEta[0] = (TH1D*)File_Kplus->Get("h_FrameEta");
  h_FramePhi[0] = (TH1D*)File_Kplus->Get("h_FramePhi");
  h_FrameEta[1] = (TH1D*)File_Kminus->Get("h_FrameEta");
  h_FramePhi[1] = (TH1D*)File_Kminus->Get("h_FramePhi");

  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
    {
      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
	string KEY = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	// cout << "KEY = " << KEY.c_str() << endl;
	h_EffKplus[KEY] = (TH1D*)File_Kplus->Get(KEY.c_str());
	h_EffKminus[KEY] = (TH1D*)File_Kminus->Get(KEY.c_str());
      }	
    }
  }

  // string outputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/Efficiency/Eff_%s_SingleKaon_%s_%s_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str(),jobID);
  string outputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/Phi/Efficiency/Eff_%s_SingleKaon_%s_%s_%d_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str(),jobID);
  cout << "OutPut File set to: " << outputfile.c_str() << endl;
  File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();

  int BufSize = (int)pow(2., 16.);
  // int Split = 1;

  const char* varlist = "Centrality:McPt:McP:McEta:McY:McPhi:McInvMass:McPid:" // MC phi 
                        "KpMcPt:KpMcEta:KpMcY:KpMcPhi:KpMcM:KpMcPid:" // MC K+ 
                        "KpRcPt:KpRcEta:KpRcY:KpRcPhi:KpRcM:KpRcTpc:" // RC K+
                        "KmMcPt:KmMcEta:KmMcY:KmMcPhi:KmMcM:KmMcPid:" // MC K-
                        "KmRcPt:KmRcEta:KmRcY:KmRcPhi:KmRcM:KmRcTpc:" // RC K-
                        "RcPt:RcP:RcEta:RcY:RcPhi:RcInvMass"; // reconstructed phi

  McPhiMeson = new TNtuple("McPhiMeson", "McPhiMeson", varlist, BufSize);
}

void readTofEff(int energy)
{
  // string inputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/ToFMatch/Eff_%s_ToFMatch_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  cout << "OPEN Efficiency File for K+ and K-: " << inputfile.c_str() << endl;

  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_eta = 0; i_eta < vmsa::BinEta+2; ++i_eta)
    {
      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
	string KEY;
	KEY = Form("h_mEfficiency_Kplus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_TofKplus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());
	// cout << "Kplus KEY: " << KEY.c_str() << endl;

	KEY = Form("h_mEfficiency_Kminus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_TofKminus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());
	// cout << "Kminus KEY: " << KEY.c_str() << endl;
      }	
    }
  }

  h_FrameEta_ToF[0] = (TH1D*)File_InPut->Get("h_FrameEta_ToF")->Clone();
  h_FramePhi_ToF[0] = (TH1D*)File_InPut->Get("h_FramePhi_ToF")->Clone();
  h_FrameEta_ToF[1] = (TH1D*)File_InPut->Get("h_FrameEta_ToF")->Clone();
  h_FramePhi_ToF[1] = (TH1D*)File_InPut->Get("h_FramePhi_ToF")->Clone();
}

void readTofEffFit(int energy)
{
  // string inputKplus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/FitPar_AuAu%s_Kplus_first.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputKplus = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/ToFMatch/FitPar_AuAu%s_Kplus_second_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Kplus = TFile::Open(inputKplus.c_str());
  cout << "OPEN ToF Matching Efficiency Fit File for K+: " << inputKplus.c_str() << endl;

  // string inputKminus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/FitPar_AuAu%s_Kminus_first.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputKminus = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/ToFMatch/FitPar_AuAu%s_Kminus_second_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Kminus = TFile::Open(inputKminus.c_str());
  cout << "OPEN ToF Matching Efficiency Fit File for K-: " << inputKminus.c_str() << endl;

  for(int i_eta = 0; i_eta < vmsa::BinEta+2; ++i_eta) // eta
  {
    TH1D *h_Kplus = NULL;
    string KEY;
    KEY = Form("h_mFitParameters_Kplus_Cent_9_Eta_%d",i_eta);
    h_Kplus = (TH1D*)File_Kplus->Get(KEY.c_str());
    // cout << "Kplus KEY: " << KEY.c_str() << endl;
    KEY = Form("f_mToFMatch_Kplus_Cent_9_Eta_%d",i_eta);
    f_TofKplus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.2,10,7);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_TofKplus[KEY]->FixParameter(i_par,h_Kplus->GetBinContent(i_par+1));
      // if(i_par < 2) cout << "Kplus: i_par = " << i_par << ", par = " << f_TofKplus[KEY]->GetParameter(i_par) << endl;
    }

    TH1D *h_Kminus = NULL;
    KEY = Form("h_mFitParameters_Kminus_Cent_9_Eta_%d",i_eta);
    h_Kminus = (TH1D*)File_Kminus->Get(KEY.c_str());
    // cout << "Kminus KEY: " << KEY.c_str() << endl;
    KEY = Form("f_mToFMatch_Kminus_Cent_9_Eta_%d",i_eta);
    f_TofKminus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.2,10,7);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_TofKminus[KEY]->FixParameter(i_par,h_Kminus->GetBinContent(i_par+1));
      // if(i_par < 2) cout << "Kminus: i_par = " << i_par << ", par = " << f_TofKminus[KEY]->GetParameter(i_par) << endl;
    }
  }

  for(int i_eta = 0; i_eta < vmsa::BinEta+2; ++i_eta) // eta & phi
  {
    for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
    {
      TH1D *h_Kplus = NULL;
      string KEY;
      KEY = Form("h_mFitParameters_Kplus_Cent_9_Eta_%d_Phi_%d",i_eta,i_phi);
      h_Kplus = (TH1D*)File_Kplus->Get(KEY.c_str());
      // cout << "Kplus KEY: " << KEY.c_str() << endl;
      KEY = Form("f_mToFMatch_Kplus_Cent_9_Eta_%d_Phi_%d",i_eta,i_phi);
      f_TofKplus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.2,10,7);
      for(int i_par = 0; i_par < 7; ++i_par)
      {
	f_TofKplus[KEY]->FixParameter(i_par,h_Kplus->GetBinContent(i_par+1));
	// if(i_par < 2) cout << "Kplus: i_par = " << i_par << ", par = " << f_TofKplus[KEY]->GetParameter(i_par) << endl;
      }

      TH1D *h_Kminus = NULL;
      KEY = Form("h_mFitParameters_Kminus_Cent_9_Eta_%d_Phi_%d",i_eta,i_phi);
      h_Kminus = (TH1D*)File_Kminus->Get(KEY.c_str());
      // cout << "Kminus KEY: " << KEY.c_str() << endl;
      KEY = Form("f_mToFMatch_Kminus_Cent_9_Eta_%d_Phi_%d",i_eta,i_phi);
      f_TofKminus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.2,10,7);
      for(int i_par = 0; i_par < 7; ++i_par)
      {
	f_TofKminus[KEY]->FixParameter(i_par,h_Kminus->GetBinContent(i_par+1));
	// if(i_par < 2) cout << "Kminus: i_par = " << i_par << ", par = " << f_TofKminus[KEY]->GetParameter(i_par) << endl;
      }
    }	
  }
}

void write()
{
  File_OutPut->cd();
  McPhiMeson->Write();
  File_OutPut->Close();
}
