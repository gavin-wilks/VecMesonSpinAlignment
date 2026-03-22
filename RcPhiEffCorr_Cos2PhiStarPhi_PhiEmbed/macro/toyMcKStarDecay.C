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
#include "TH3F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StRoot/Utility/functions.h"

using namespace std;

typedef std::map<std::string,TH1D*> TH1DMap;
typedef std::map<std::string,TF1*> TF1Map;

void readEfficiency(int energy, int year, int cut, const char *jobID);
void readTofEff(int energy);
void readTofEffFit(int energy);
void getKinematics(TLorentzVector& lKStar, double const mass);
void setDecayChannels(int const pid);
void decayAndFill(int const kf, int const cent, TLorentzVector* lKStar, TClonesArray& daughters);
void fill(int const kf, int const centrality, TLorentzVector* lKStar, TLorentzVector const& lK, TLorentzVector const& lPi);
bool tpcReconstructed(int iParticleIndex, int cent, float Psi2, TLorentzVector const& lKaon);
void findHist(TLorentzVector const& lKaon, int iParticleIndex, float Psi2, int& EtaBin, int& PhiBin); // iParticleIndex = 0 => K+
float AngleShift(float phi);
void findHist_ToF(TLorentzVector const& lKaon, int iParticleIndex, int& EtaBin, int& PhiBin); // iParticleIndex = 0 => K+
TF1* readv2(int energy, int pid, int centrality);
TF1* readspec(int energy, int pid, int centrality);
void write();

TPythia6Decayer* pydecay;
TNtuple* McKStarMeson;

TH1DMap h_EffKplus;
TH1DMap h_EffKminus;
TH1DMap h_EffPiplus;
TH1DMap h_EffPiminus;
TH1D *h_FrameEta[4];
TH1D *h_FramePhi[4];

TH1DMap h_TofKplus;
TH1DMap h_TofKminus;
TH1DMap h_TofPiplus;
TH1DMap h_TofPiminus;
TH1D *h_FrameEta_ToF[4];
TH1D *h_FramePhi_ToF[4];

TF1Map f_TofKplus;
TF1Map f_TofKminus;
TF1Map f_TofPiplus;
TF1Map f_TofPiminus;

// sampling functions and histograms
TF1 *f_v2, *f_spec, *f_flow;
TH3F *h_Tracks;
TH2F *h_phiRP;

TFile *File_OutPut;

Double_t pt_set[3] = {1.0, 1.5, 5.0};
int pt_bin;

void toyMcKStarDecay(const int energy = 4, const int pt = 0, const int centrality = 9, const int pid = 2, const int year = 0, const int cut = 0, const int NMax = 1000, const char* jobID = "testing")
{
  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();
  gRandom->SetSeed();
  readEfficiency(energy,year,cut,jobID);
  //readTofEff(energy);
  //readTofEffFit(energy);

  pt_bin = pt;
  // v2 & spectra implementation
  f_v2   = readv2(energy,pid,centrality);
  f_spec = readspec(energy,pid,centrality);
  f_flow = new TF1("f_flow",flowSample,-TMath::Pi(),TMath::Pi(),1);

  h_Tracks = new TH3F("h_Tracks","h_Tracks",20,vmsa::ptMin,vmsa::ptMax,vmsa::BinY,-1.0,1.0,36,-TMath::Pi(),TMath::Pi());
  h_phiRP = new TH2F("h_phiRP","h_phiRP",20,vmsa::ptMin,vmsa::ptMax,36,-TMath::Pi(),TMath::Pi()); // QA histogram for v2 sample

  pydecay = TPythia6Decayer::Instance();
  pydecay->Init();
  setDecayChannels(pid); // K*0 --> Kpi

  TClonesArray ptl("TParticle", 10);
  TLorentzVector *lKStar = new TLorentzVector();
  for(int i_ran = 0; i_ran < NMax; ++i_ran)
  {
    if (floor(10.0*i_ran/ static_cast<float>(NMax)) > floor(10.0*(i_ran-1)/ static_cast<float>(NMax)))
    cout << "=> processing data: " << 100.0*i_ran/ static_cast<float>(NMax) << "%" << endl;
  //  setDecayChannels(pid); // K*0 --> Kpi

    getKinematics(*lKStar,vmsa::InvMass[pid]);
    // if (fabs(lKStar->Phi()) > TMath::Pi()) continue;
    decayAndFill(vmsa::decayMother[pid],centrality,lKStar,ptl);

    if (i_ran % 1000 == 1) McKStarMeson->AutoSave("SaveSelf");
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  write();

  stopWatch->Stop();   
  stopWatch->Print();
}

void getKinematics(TLorentzVector& lKStar, double const mass)
{
  //double const pt = gRandom->Uniform(vmsa::ptMin, vmsa::ptEffMax); // sample with flat distribution
  double const pt = f_spec->GetRandom(pt_set[pt_bin], pt_set[pt_bin+1]); // sample with measured spectra
  double const y = gRandom->Uniform(-vmsa::acceptanceRapidity, vmsa::acceptanceRapidity);
  // double const phi = TMath::TwoPi() * gRandom->Rndm(); // sample flat distribution
  //double const phi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
  f_flow->ReleaseParameter(0);
  f_flow->SetParameter(0,f_v2->Eval(pt));
  double const phi = f_flow->GetRandom(); // sample with measured v2

  double const mT = sqrt(mass * mass + pt * pt);
  double const pz = mT * sinh(y);
  double const E = mT * cosh(y);

  lKStar.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

void setDecayChannels(int const pid)
{
  //float const randAnti = gRandom->Uniform(0.0,1.0);
  //bool anti = false;
  //if(randAnti < 0.5) anti = true;
  int const mdme = vmsa::decayChannels[pid];

  //if(!anti)
 // {
    for (int idc = vmsa::decayChannelsFirst[pid]; idc < vmsa::decayChannelsSecond[pid] + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0); // close all decay channel
    TPythia6::Instance()->SetMDME(mdme, 1, 1); // open the one we need
    int *PYSeed = new int;
    TPythia6::Instance()->SetMRPY(1,(int)PYSeed); // Random seed
    //cout << "mdme = " << mdme << endl;
  //}
  //if(anti)
 // {
  //  for (int idc = vmsa::decayChannelsFirst[pid]; idc < vmsa::decayChannelsSecond[pid] + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0); // close all decay channel
    //TPythia6::Instance()->SetMDME(mdme, 1, 0); // open the one we need
 //   TPythia6::Instance()->SetMDME(mdme, 1, 3); // open the one we need
 //   int *PYSeed = new int;
  //  TPythia6::Instance()->SetMRPY(1,(int)PYSeed); // Random seed
  //  cout << "mdme = " << -mdme << endl;
 // }
}

void decayAndFill(int const kf, int const cent, TLorentzVector* lKStar, TClonesArray& daughters)
{
  float const randAnti = gRandom->Uniform(0.0,1.0);
  bool anti = false;
  if(randAnti < 0.5) anti = true;
  
  if(!anti) pydecay->Decay(kf, lKStar);
  if(anti)  pydecay->Decay(-kf, lKStar);
  
  pydecay->ImportParticles(&daughters);

  TLorentzVector lK;
  TLorentzVector lPi;

  int nTrk = daughters.GetEntriesFast();
  for (int iTrk = 0; iTrk < nTrk; ++iTrk)
  {
    TParticle* ptl0 = (TParticle*)daughters.At(iTrk);

    switch (ptl0->GetPdgCode())
    {
      case 321:
	ptl0->Momentum(lK);
	break;
      case -321:
	ptl0->Momentum(lK);
	break;
      case 211:
        ptl0->Momentum(lPi); 
        break;
      case -211:
        ptl0->Momentum(lPi);
        break;
      default:
	break;
    }
  }
  daughters.Clear("C");

  fill(kf,cent,lKStar,lK,lPi);
}

void fill(int const kf, int const centrality, TLorentzVector* lKStar, TLorentzVector const& lK, TLorentzVector const& lPi)
{
  //int const centrality = floor(vmsa::NCentMax * gRandom->Rndm());
  //float const Psi2 = gRandom->Uniform(-0.5*TMath::Pi(),0.5*TMath::Pi()); // random event plane angle
  float const Psi2 = 0.0; // fixed event plane angle
  // cout << "centrality = " << centrality << ", Psi2 = " << Psi2 << endl;
  // int const centrality = floor(2 * gRandom->Rndm());
  TLorentzVector lRcKStar = lK + lPi; // phi meson reconstruction
  // cout << "lKStar.pt = " << lKStar->Pt() << ", lKStar.eta = " << lKStar->Eta() << ", lKStar.phi = " << lKStar->Phi() << ", lKStar.m = " << lKStar->M() << endl;
  // cout << "lRcKStar.pt = " << lRcPhi.Pt() << ", lRcPhi.eta = " << lRcPhi.Eta() << ", lRcPhi.phi = " << lRcPhi.Phi() << ", lRcPhi.m = " << lRcPhi.M() << endl;
  // cout << endl;

  float arr[110];
  int iArr = 0;
  arr[iArr++] = centrality; // McPhi
  arr[iArr++] = Psi2; 
  arr[iArr++] = lKStar->Pt();
  arr[iArr++] = lKStar->P();
  arr[iArr++] = lKStar->PseudoRapidity();
  arr[iArr++] = lKStar->Rapidity();
  arr[iArr++] = lKStar->Phi();
  arr[iArr++] = lKStar->M();
  arr[iArr++] = kf;

  arr[iArr++] = lK.Pt();
  arr[iArr++] = lK.PseudoRapidity();
  arr[iArr++] = lK.Rapidity();
  arr[iArr++] = lK.Phi();
  arr[iArr++] = lK.M();
  arr[iArr++] = (kf > 0) ? 321 : -321 ;

  arr[iArr++] = lK.Pt();
  arr[iArr++] = lK.PseudoRapidity();
  arr[iArr++] = lK.Rapidity();
  arr[iArr++] = lK.Phi();
  arr[iArr++] = lK.M();
  arr[iArr++] = (kf > 0) ? tpcReconstructed(0,centrality,Psi2,lK) : tpcReconstructed(1,centrality,Psi2,lK);

  arr[iArr++] = lPi.Pt();
  arr[iArr++] = lPi.PseudoRapidity();
  arr[iArr++] = lPi.Rapidity();
  arr[iArr++] = lPi.Phi();
  arr[iArr++] = lPi.M();
  arr[iArr++] = (kf > 0) ? -211 : 211 ;

  arr[iArr++] = lPi.Pt();
  arr[iArr++] = lPi.PseudoRapidity();
  arr[iArr++] = lPi.Rapidity();
  arr[iArr++] = lPi.Phi();
  arr[iArr++] = lPi.M();
  arr[iArr++] = (kf > 0) ? tpcReconstructed(3,centrality,Psi2,lPi) : tpcReconstructed(2,centrality,Psi2,lPi);

  arr[iArr++] = lRcKStar.Pt();
  arr[iArr++] = lRcKStar.P();
  arr[iArr++] = lRcKStar.PseudoRapidity();
  arr[iArr++] = lRcKStar.Rapidity();
  arr[iArr++] = lRcKStar.Phi();
  arr[iArr++] = lRcKStar.M();

  McKStarMeson->Fill(arr);
  // if(lRcPhi.Pt() < 10e-4) cout << "lKStar->Pt = " << lKStar->Pt() << ", lRcPhi.Pt = " << lRcPhi.Pt() << endl;

  h_Tracks->Fill(lKStar->Pt(),lKStar->Rapidity(),lKStar->Phi()); // fill QA histograms
  h_phiRP->Fill(lKStar->Pt(),lKStar->Phi());
}

bool tpcReconstructed(int iParticleIndex, int cent, float Psi2, TLorentzVector const& lKaon)
{
   if(fabs(lKaon.Eta()) >= vmsa::mEtaMax) return false;

   TH1D *h_TPC = NULL;
   int EtaBin_TPC = -1;
   int PhiBin_TPC = -1;
   findHist(lKaon,iParticleIndex,Psi2,EtaBin_TPC,PhiBin_TPC);

   //TH1D *h_ToF = NULL;
   //TF1 *f_ToF = NULL;
   //int EtaBin_ToF = -1;
   //int PhiBin_ToF = -1;
   //findHist_ToF(lKaon,iParticleIndex,EtaBin_ToF,PhiBin_ToF);

   // cout << KEY.c_str() << endl;
   if (iParticleIndex == 0)
   {
     // string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
     string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff @ 20-60%
     h_TPC = h_EffKplus[KEY_TPC];

     //string KEY_ToF = Form("h_mEfficiency_Kplus_Cent_9_Eta_%d_Phi_%d",EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
     ////string KEY_ToF = Form("h_mEfficiency_Kplus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
     //h_ToF = h_TofKplus[KEY_ToF];

     //string KEY_ToFFit = Form("f_mToFMatch_Kplus_Cent_9_Eta_%d_Phi_%d",EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
     ////string KEY_ToFFit = Form("f_mToFMatch_Kplus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
     //f_ToF = f_TofKplus[KEY_ToFFit]; // only 20-60%
   }
   if (iParticleIndex == 2)
   {
     // string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
     string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff @ 20-60%
     h_TPC = h_EffPiplus[KEY_TPC];

     //string KEY_ToF = Form("h_mEfficiency_Piplus_Cent_9_Eta_%d_Phi_%d",EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
     ////string KEY_ToF = Form("h_mEfficiency_Piplus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
     //h_ToF = h_TofPiplus[KEY_ToF];

     //string KEY_ToFFit = Form("f_mToFMatch_Piplus_Cent_9_Eta_%d_Phi_%d",EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
     ////string KEY_ToFFit = Form("f_mToFMatch_Piplus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
     //f_ToF = f_TofPiplus[KEY_ToFFit]; // only 20-60%
   }
   if (iParticleIndex == 1)
   {
     // string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
     string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff @ 20-60%
     h_TPC = h_EffKminus[KEY_TPC];

     //string KEY_ToF = Form("h_mEfficiency_Kminus_Cent_9_Eta_%d_Phi_%d",EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
     ////string KEY_ToF = Form("h_mEfficiency_Kminus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
     //h_ToF = h_TofKminus[KEY_ToF];

     //string KEY_ToFFit = Form("f_mToFMatch_Kminus_Cent_9_Eta_%d_Phi_%d",EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
     ////string KEY_ToFFit = Form("f_mToFMatch_Kminus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
     //f_ToF = f_TofKminus[KEY_ToFFit]; // only 20-60%
   }
   if (iParticleIndex == 3)
   {
     // string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
     string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff @ 20-60%
     h_TPC = h_EffPiminus[KEY_TPC];

     //string KEY_ToF = Form("h_mEfficiency_Piminus_Cent_9_Eta_%d_Phi_%d",EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
     ////string KEY_ToF = Form("h_mEfficiency_Piminus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
     //h_ToF = h_TofPiminus[KEY_ToF];

     //string KEY_ToFFit = Form("f_mToFMatch_Piminus_Cent_9_Eta_%d_Phi_%d",EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
     ////string KEY_ToFFit = Form("f_mToFMatch_Piminus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
     //f_ToF = f_TofPiminus[KEY_ToFFit]; // only 20-60%
   }

   double pt = lKaon.Perp();
   int const bin_TPC = h_TPC->FindBin(pt);
   bool is_TPC = gRandom->Rndm() < h_TPC->GetBinContent(bin_TPC);

   //int const bin_ToF = h_ToF->FindBin(pt); // tof fit and hist combined
   //double prob_tof = 0.0;
   //if(pt > 0.3) prob_tof = f_ToF->Eval(pt); // donot use fit extrapolation
   //else prob_tof = h_ToF->GetBinContent(bin_ToF);
   //bool is_ToF = gRandom->Rndm() < prob_tof;

   // int const bin_ToF = h_ToF->FindBin(pt); // tof hist only
   // bool is_ToF = gRandom->Rndm() < h_ToF->GetBinContent(bin_ToF);

   // bool is_ToF = gRandom->Rndm() < f_ToF->Eval(pt); // tof fit only

   // cout << "is_TPC: " << is_TPC << ", is_ToF: " << is_ToF << ", is_TPC && is_ToF: " << (is_TPC && is_ToF) << endl;
   // cout << "pt = " << pt << ", h_ToF = " << h_ToF->GetBinContent(bin_ToF) << ", f_ToF = " << f_ToF->Eval(pt) << ", prob_tof = " << prob_tof << endl; // comparison between eff from hist and func

   return is_TPC; //&& is_ToF;
   //return true;
}

void findHist(TLorentzVector const& lKaon, int iParticleIndex, float Psi2, int& EtaBin, int& PhiBin)
{
  float eta = lKaon.Eta();
  EtaBin = h_FrameEta[iParticleIndex]->FindBin(eta)-1;
  // cout << "eta = " << eta << ", EtaBin = " << EtaBin << endl;
  //float phi = lKaon.Phi();
  float phi = lKaon.Phi()-Psi2;
  float phi_shift = AngleShift(phi);
  PhiBin = h_FramePhi[iParticleIndex]->FindBin(phi_shift)-1;
  //PhiBin = h_FramePhi[iParticleIndex]->FindBin(phi)-1;
}

float AngleShift(float phi)
{
  double const Psi2_low[3] = {-3.0*TMath::Pi()/2.0,-1.0*TMath::Pi()/2.0,1.0*TMath::Pi()/2.0};
  double const Psi2_up[3]  = {-1.0*TMath::Pi()/2.0, 1.0*TMath::Pi()/2.0,3.0*TMath::Pi()/2.0};

  float phi_shift = -999.0;
  for(int psi_bin = 0; psi_bin < 3; ++psi_bin)
  {
    if(phi >= Psi2_low[psi_bin] && phi < Psi2_up[psi_bin])
    {
      phi_shift = phi - (psi_bin-1)*2.0*TMath::Pi()/2.0;
    }
  }

  return phi_shift;
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

void readEfficiency(int energy, int year, int cut, const char *jobID)
{
  // string inputKplus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[0].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  string inputKplus = Form("/star/data01/pwg/gwilks3/VectorMesonSpinAlignment/Data/KStar/Efficiency/TPC/Eff_%s_%s.root",vmsa::mParType[0].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Kplus = TFile::Open(inputKplus.c_str());
  cout << "OPEN Efficiency File for K+: " << inputKplus.c_str() << endl;

  // string inputKminus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[1].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  string inputKminus = Form("/star/data01/pwg/gwilks3/VectorMesonSpinAlignment/Data/KStar/Efficiency/TPC/Eff_%s_%s.root",vmsa::mParType[1].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Kminus = TFile::Open(inputKminus.c_str());
  cout << "OPEN Efficiency File for K-: " << inputKminus.c_str() << endl;

  // string inputKplus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[0].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  string inputPiplus = Form("/star/data01/pwg/gwilks3/VectorMesonSpinAlignment/Data/KStar/Efficiency/TPC/Eff_%s_%s.root",vmsa::mParType[2].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Piplus = TFile::Open(inputPiplus.c_str());
  cout << "OPEN Efficiency File for pi+: " << inputPiplus.c_str() << endl;

  // string inputPiminus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[1].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  string inputPiminus = Form("/star/data01/pwg/gwilks3/VectorMesonSpinAlignment/Data/KStar/Efficiency/TPC/Eff_%s_%s.root",vmsa::mParType[3].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Piminus = TFile::Open(inputPiminus.c_str());
  cout << "OPEN Efficiency File for pi-: " << inputPiminus.c_str() << endl;

  h_FrameEta[0] = (TH1D*)File_Kplus->Get("h_FrameEta");
  h_FramePhi[0] = (TH1D*)File_Kplus->Get("h_FramePhi");
  h_FrameEta[1] = (TH1D*)File_Kminus->Get("h_FrameEta");
  h_FramePhi[1] = (TH1D*)File_Kminus->Get("h_FramePhi");
  h_FrameEta[2] = (TH1D*)File_Piplus->Get("h_FrameEta");
  h_FramePhi[2] = (TH1D*)File_Piplus->Get("h_FramePhi");
  h_FrameEta[3] = (TH1D*)File_Piminus->Get("h_FrameEta");
  h_FramePhi[3] = (TH1D*)File_Piminus->Get("h_FramePhi");

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
	h_EffPiplus[KEY] = (TH1D*)File_Piplus->Get(KEY.c_str());
	h_EffPiminus[KEY] = (TH1D*)File_Piminus->Get(KEY.c_str());
      }	
    }
  }
 
  // string outputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/Efficiency/Eff_%s_SingleKaon_%s_%s_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str(),jobID);
  string outputfile = Form("./Eff_%s_SingleParticle_2060_%s_KStar.root",vmsa::mBeamEnergy[energy].c_str(),jobID);
  cout << "OutPut File set to: " << outputfile.c_str() << endl;
  File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();

  int BufSize = (int)pow(2., 16.);
  // int Split = 1;

  const char* varlist = "Centrality:Psi2:McPt:McP:McEta:McY:McPhi:McInvMass:McPid:" // MC phi 
                        "KpMcPt:KpMcEta:KpMcY:KpMcPhi:KpMcM:KpMcPid:" // MC K+ 
                        "KpRcPt:KpRcEta:KpRcY:KpRcPhi:KpRcM:KpRcTpc:" // RC K+
                        "KmMcPt:KmMcEta:KmMcY:KmMcPhi:KmMcM:KmMcPid:" // MC K-
                        "KmRcPt:KmRcEta:KmRcY:KmRcPhi:KmRcM:KmRcTpc:" // RC K-
                        "RcPt:RcP:RcEta:RcY:RcPhi:RcInvMass"; // reconstructed phi

  McKStarMeson = new TNtuple("McKStarMeson", "McKStarMeson", varlist, BufSize);
}

void readTofEff(int energy)
{
  // string inputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Efficiency/ToF/Eff_%s_ToFMatch_2060.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  cout << "OPEN Efficiency File for K+, K-, pi+, pi-: " << inputfile.c_str() << endl;

  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
    {
      string KEY;
      KEY = Form("h_mEfficiency_Kplus_Cent_%d_Eta_%d",i_cent,i_eta);
      h_TofKplus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());
      // cout << "Kplus KEY: " << KEY.c_str() << endl;

      KEY = Form("h_mEfficiency_Kminus_Cent_%d_Eta_%d",i_cent,i_eta);
      h_TofKminus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());
      // cout << "Kminus KEY: " << KEY.c_str() << endl;

      KEY = Form("h_mEfficiency_Piplus_Cent_%d_Eta_%d",i_cent,i_eta);
      h_TofPiplus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());

      KEY = Form("h_mEfficiency_Piminus_Cent_%d_Eta_%d",i_cent,i_eta);
      h_TofPiminus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());

      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
	// string KEY;
	KEY = Form("h_mEfficiency_Kplus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_TofKplus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());
	// cout << "Kplus KEY: " << KEY.c_str() << endl;

	KEY = Form("h_mEfficiency_Kminus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_TofKminus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());
	// cout << "Kminus KEY: " << KEY.c_str() << endl;

	KEY = Form("h_mEfficiency_Piplus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_TofPiplus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());
	// cout << "Kplus KEY: " << KEY.c_str() << endl;

	KEY = Form("h_mEfficiency_Piminus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_TofPiminus[KEY] = (TH1D*)File_InPut->Get(KEY.c_str());
	// cout << "Kminus KEY: " << KEY.c_str() << endl;
      }	
    }
  }

  h_FrameEta_ToF[0] = (TH1D*)File_InPut->Get("h_FrameEta_ToF")->Clone();
  h_FramePhi_ToF[0] = (TH1D*)File_InPut->Get("h_FramePhi_ToF")->Clone();
  h_FrameEta_ToF[1] = (TH1D*)File_InPut->Get("h_FrameEta_ToF")->Clone();
  h_FramePhi_ToF[1] = (TH1D*)File_InPut->Get("h_FramePhi_ToF")->Clone();
  h_FrameEta_ToF[2] = (TH1D*)File_InPut->Get("h_FrameEta_ToF")->Clone();
  h_FramePhi_ToF[2] = (TH1D*)File_InPut->Get("h_FramePhi_ToF")->Clone();
  h_FrameEta_ToF[3] = (TH1D*)File_InPut->Get("h_FrameEta_ToF")->Clone();
  h_FramePhi_ToF[3] = (TH1D*)File_InPut->Get("h_FramePhi_ToF")->Clone();
}

void readTofEffFit(int energy)
{
  string inputKplus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/KStar/Efficiency/ToF/FitPar_AuAu%s_Kplus_2060.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Kplus = TFile::Open(inputKplus.c_str());
  cout << "OPEN ToF Matching Efficiency Fit File for K+: " << inputKplus.c_str() << endl;

  string inputKminus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/KStar/Efficiency/ToF/FitPar_AuAu%s_Kminus_2060.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Kminus = TFile::Open(inputKminus.c_str());
  cout << "OPEN ToF Matching Efficiency Fit File for K-: " << inputKminus.c_str() << endl;

  string inputPiplus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/KStar/Efficiency/ToF/FitPar_AuAu%s_Piplus_2060.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Piplus = TFile::Open(inputPiplus.c_str());
  cout << "OPEN ToF Matching Efficiency Fit File for pi+: " << inputPiplus.c_str() << endl;

  string inputPiminus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/KStar/Efficiency/ToF/FitPar_AuAu%s_Piminus_2060.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Piminus = TFile::Open(inputPiminus.c_str());
  cout << "OPEN ToF Matching Efficiency Fit File for pi-: " << inputPiminus.c_str() << endl;

  for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta) // eta
  {
    TH1D *h_Kplus = NULL;
    string KEY;
    KEY = Form("h_mFitParameters_Kplus_Cent_9_Eta_%d",i_eta);
    h_Kplus = (TH1D*)File_Kplus->Get(KEY.c_str());
    // cout << "Kplus KEY: " << KEY.c_str() << endl;
    KEY = Form("f_mToFMatch_Kplus_Cent_9_Eta_%d",i_eta);
    f_TofKplus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
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
    f_TofKminus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_TofKminus[KEY]->FixParameter(i_par,h_Kminus->GetBinContent(i_par+1));
      // if(i_par < 2) cout << "Kminus: i_par = " << i_par << ", par = " << f_TofKminus[KEY]->GetParameter(i_par) << endl;
    }

    TH1D *h_Piplus = NULL;
    KEY = Form("h_mFitParameters_Piplus_Cent_9_Eta_%d",i_eta);
    h_Piplus = (TH1D*)File_Piplus->Get(KEY.c_str());
    // cout << "Kplus KEY: " << KEY.c_str() << endl;
    KEY = Form("f_mToFMatch_Piplus_Cent_9_Eta_%d",i_eta);
    f_TofPiplus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_TofPiplus[KEY]->FixParameter(i_par,h_Piplus->GetBinContent(i_par+1));
      // if(i_par < 2) cout << "Kplus: i_par = " << i_par << ", par = " << f_TofKplus[KEY]->GetParameter(i_par) << endl;
    }

    TH1D *h_Piminus = NULL;
    KEY = Form("h_mFitParameters_Piminus_Cent_9_Eta_%d",i_eta);
    h_Piminus = (TH1D*)File_Piminus->Get(KEY.c_str());
    // cout << "Kminus KEY: " << KEY.c_str() << endl;
    KEY = Form("f_mToFMatch_Piminus_Cent_9_Eta_%d",i_eta);
    f_TofPiminus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_TofPiminus[KEY]->FixParameter(i_par,h_Piminus->GetBinContent(i_par+1));
      // if(i_par < 2) cout << "Kminus: i_par = " << i_par << ", par = " << f_TofKminus[KEY]->GetParameter(i_par) << endl;
    }
  }

  for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta) // eta & phi
  {
    for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
    {
      TH1D *h_Kplus = NULL;
      string KEY;
      KEY = Form("h_mFitParameters_Kplus_Cent_9_Eta_%d_Phi_%d",i_eta,i_phi);
      h_Kplus = (TH1D*)File_Kplus->Get(KEY.c_str());
      // cout << "Kplus KEY: " << KEY.c_str() << endl;
      KEY = Form("f_mToFMatch_Kplus_Cent_9_Eta_%d_Phi_%d",i_eta,i_phi);
      f_TofKplus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
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
      f_TofKminus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
      for(int i_par = 0; i_par < 7; ++i_par)
      {
	f_TofKminus[KEY]->FixParameter(i_par,h_Kminus->GetBinContent(i_par+1));
	// if(i_par < 2) cout << "Kminus: i_par = " << i_par << ", par = " << f_TofKminus[KEY]->GetParameter(i_par) << endl;
      }
      
      TH1D *h_Piplus = NULL;
      string KEY;
      KEY = Form("h_mFitParameters_Piplus_Cent_9_Eta_%d_Phi_%d",i_eta,i_phi);
      h_Piplus = (TH1D*)File_Piplus->Get(KEY.c_str());
      // cout << "Kplus KEY: " << KEY.c_str() << endl;
      KEY = Form("f_mToFMatch_Piplus_Cent_9_Eta_%d_Phi_%d",i_eta,i_phi);
      f_TofPiplus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
      for(int i_par = 0; i_par < 7; ++i_par)
      {
	f_TofPiplus[KEY]->FixParameter(i_par,h_Piplus->GetBinContent(i_par+1));
	// if(i_par < 2) cout << "Kplus: i_par = " << i_par << ", par = " << f_TofKplus[KEY]->GetParameter(i_par) << endl;
      }

      TH1D *h_Piminus = NULL;
      KEY = Form("h_mFitParameters_Piminus_Cent_9_Eta_%d_Phi_%d",i_eta,i_phi);
      h_Piminus = (TH1D*)File_Piminus->Get(KEY.c_str());
      // cout << "Kminus KEY: " << KEY.c_str() << endl;
      KEY = Form("f_mToFMatch_Piminus_Cent_9_Eta_%d_Phi_%d",i_eta,i_phi);
      f_TofPiminus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
      for(int i_par = 0; i_par < 7; ++i_par)
      {
	f_TofPiminus[KEY]->FixParameter(i_par,h_Piminus->GetBinContent(i_par+1));
	// if(i_par < 2) cout << "Kminus: i_par = " << i_par << ", par = " << f_TofKminus[KEY]->GetParameter(i_par) << endl;
      }
    }	
  }
}

//TF1* readv2(int energy, int pid)
//{
//  string InPutV2 = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/KStar/MonteCarlo/Data/KStar_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
//  TFile *File_v2 = TFile::Open(InPutV2.c_str());
//  TGraphAsymmErrors *g_v2 = (TGraphAsymmErrors*)File_v2->Get("g_v2");
//  TF1 *f_v2 = new TF1("f_v2",v2_pT_FitFunc,vmsa::ptMin,vmsa::ptMax,5);
//  f_v2->FixParameter(0,2);
//  f_v2->SetParameter(1,0.1);
//  f_v2->SetParameter(2,0.1);
//  f_v2->SetParameter(3,0.1);
//  f_v2->SetParameter(4,0.1);
//  f_v2->SetLineColor(2);
//  f_v2->SetLineWidth(2);
//  f_v2->SetLineStyle(2);
//  g_v2->Fit(f_v2,"N");
//
//  /*
//  TCanvas *c_v2 = new TCanvas("c_v2","c_v2",10,10,800,800);
//  c_v2->cd()->SetLeftMargin(0.15);
//  c_v2->cd()->SetBottomMargin(0.15);
//  c_v2->cd()->SetTicks(1,1);
//  c_v2->cd()->SetGrid(0,0);
//  TH1F *h_v2 = new TH1F("h_v2","h_v2",100,0.0,10.0);
//  for(int i_bin = 1; i_bin < 101; ++i_bin)
//  {
//    h_v2->SetBinContent(i_bin,-10.0);
//    h_v2->SetBinError(i_bin,1.0);
//  }
//  h_v2->SetTitle("");
//  h_v2->SetStats(0);
//  h_v2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
//  h_v2->GetXaxis()->CenterTitle();
//  h_v2->GetYaxis()->SetTitle("v_{2}");
//  h_v2->GetYaxis()->CenterTitle();
//  h_v2->GetYaxis()->SetRangeUser(0.0,0.2);
//  h_v2->Draw("pE");
//  g_v2->Draw("pE same");
//  f_v2->Draw("l same");
//  */
//
//  return f_v2;
//}
//
//TF1* readspec(int energy, int pid)
//{
//  string InPutSpec = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/KStar/MonteCarlo/Data/KStar_Spec.root",vmsa::mBeamEnergy[energy].c_str());
//  TFile *File_Spec = TFile::Open(InPutSpec.c_str());
//  TGraphAsymmErrors *g_spec = (TGraphAsymmErrors*)File_Spec->Get("g_spec");
//  TF1 *f_Levy = new TF1("f_Levy",Levy,vmsa::ptMin,vmsa::ptMax,3);
//  f_Levy->SetParameter(0,1);
//  f_Levy->SetParameter(1,10);
//  f_Levy->SetParameter(2,0.1);
//  f_Levy->SetLineStyle(2);
//  f_Levy->SetLineColor(4);
//  f_Levy->SetLineWidth(2);
//  g_spec->Fit(f_Levy,"N");
//
//  TF1 *f_spec = new TF1("f_spec",pTLevy,vmsa::ptMin,vmsa::ptMax,3);
//  f_spec->SetParameter(0,f_Levy->GetParameter(0));
//  f_spec->SetParameter(1,f_Levy->GetParameter(1));
//  f_spec->SetParameter(2,f_Levy->GetParameter(2));
//  f_spec->SetLineStyle(2);
//  f_spec->SetLineColor(2);
//  f_spec->SetLineWidth(2);
//
//  /*
//  TCanvas *c_spec = new TCanvas("c_spec","c_spec",10,10,800,800);
//  c_spec->cd()->SetLeftMargin(0.15);
//  c_spec->cd()->SetBottomMargin(0.15);
//  c_spec->cd()->SetTicks(1,1);
//  c_spec->cd()->SetGrid(0,0);
//  c_spec->SetLogy();
//  TH1F *h_spec = new TH1F("h_spec","h_spec",100,0.0,10.0);
//  for(int i_bin = 1; i_bin < 101; ++i_bin)
//  {
//    h_spec->SetBinContent(i_bin,-10.0);
//    h_spec->SetBinError(i_bin,1.0);
//  }
//  h_spec->SetTitle("");
//  h_spec->SetStats(0);
//  h_spec->GetXaxis()->SetTitle("p_{T} (GeV/c)");
//  h_spec->GetXaxis()->CenterTitle();
//  h_spec->GetYaxis()->SetTitle("dN/p_{T}dp_{T}");
//  h_spec->GetYaxis()->CenterTitle();
//  h_spec->GetYaxis()->SetRangeUser(1E-6,10);
//  h_spec->Draw("pE");
//  g_spec->Draw("pE same");
//  f_Levy->Draw("l same");
//  f_spec->Draw("l same");
//  */
//
//  return f_spec;
//}

TF1* readv2(int energy, int pid, int centrality)
{
  string InPutV2 = Form("/gpfs01/star/pwg/subhash/Public/Kstar_spec_data/kstarv2_200_Run11.root");
  cout << InPutV2 << endl;
  TFile *File_v2 = TFile::Open(InPutV2.c_str());
  TH1F *g_v2 = (TH1F*)File_v2->Get("hkstarv2_200");
  //TF1 *f_v2 = (TF1*)File_v2->Get("v2Fit");
 // TF1 *f_v2 = new TF1("f_v2",v2_pT_FitFunc_Poly3,vmsa::ptMin,vmsa::ptMax,9);
 // f_v2->FixParameter(0,2);
 // f_v2->SetParameter(1,0.1);
 // f_v2->SetParameter(2,0.1);
 // f_v2->SetParameter(3,0.1);
 // f_v2->SetParameter(4,0.1);
 // f_v2->SetLineColor(2);
 // f_v2->SetLineWidth(2);
 // f_v2->SetLineStyle(2);
 // g_v2->Fit(f_v2,"N");

  TF1 *f_v2 = f_Kstarv2_200 = new TF1("f_Kstarv2_200","(2*7.74602e-02/(1+exp(-( (x/2) - 5.31219e-01)/1.40637e-01) )) - 2*2.43233e-03", 0, 10.0); 
 
  TCanvas *c_v2 = new TCanvas("c_v2","c_v2",10,10,800,800);
  c_v2->cd()->SetLeftMargin(0.15);
  c_v2->cd()->SetBottomMargin(0.15);
  c_v2->cd()->SetTicks(1,1);
  c_v2->cd()->SetGrid(0,0);
  TH1F *h_v2 = new TH1F("h_v2","h_v2",100,0.0,10.0);
  for(int i_bin = 1; i_bin < 101; ++i_bin)
  {
    h_v2->SetBinContent(i_bin,-10.0);
    h_v2->SetBinError(i_bin,1.0);
  }
  h_v2->SetTitle("");
  h_v2->SetStats(0);
  h_v2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_v2->GetXaxis()->CenterTitle();
  h_v2->GetYaxis()->SetTitle("v_{2}");
  h_v2->GetYaxis()->CenterTitle();
  h_v2->GetYaxis()->SetRangeUser(0.0,0.2);
  h_v2->Draw("pE");
  g_v2->Draw("pE same");
  f_v2->Draw("l same");
  c_v2->SaveAs("v2.pdf");



  return f_v2;
}

TF1* readspec(int energy, int pid, int centrality)
{
  TCanvas *c1 = new TCanvas();
  c1->SetFillColor(0);
  c1->SetGrid(0,0);
  c1->SetTitle(0);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.15);

  string centString[9] = {"6080","6080","4060","4060","3040","2030","2010","010","010"};
  string InPutSpec;
  TH1F *g_spec;
  if(centrality < 9)
  {
    string InPutSpec = Form("/gpfs01/star/pwg/subhash/Public/Kstar_spec_data/Both_kstarmeson_PhPsi_Bg1_f0.77_1.04_H_c%s_Levy_%s_def_FW.root",centString[centrality].c_str(),vmsa::mBeamEnergy[energy].c_str());
    cout << InPutSpec << endl;
    TFile *File_Spec = TFile::Open(InPutSpec.c_str());
    g_spec = (TH1F*)File_Spec->Get("hcorr_yield");
  }
  if(centrality == 9)
  {
    //TGraphAsymmErrors *g_temp[3];
    for(int i = 3; i < 6; i++)
    {
      string InPutSpec = Form("/gpfs01/star/pwg/subhash/Public/Kstar_spec_data/Both_kstarmeson_PhPsi_Bg1_f0.77_1.04_H_c%s_Levy_%s_def_FW.root",centString[i].c_str(),vmsa::mBeamEnergy[energy].c_str());
      cout << InPutSpec << endl;
      TFile *File_Spec = TFile::Open(InPutSpec.c_str());
      if(i == 3) g_spec = (TH1F*)File_Spec->Get("hcorr_yield")->Clone();
      else g_spec->Add((TH1F*)File_Spec->Get("hcorr_yield"));
    }
    //const int nPoints = g_temp[i]->GetN();
    //double xval[nPoints]  = {0.0};
    //double sum[nPoints]   = {0.0};
    //double errxl[nPoints] = {0.0};
    //double errxh[nPoints] = {0.0};
    //double erryl[nPoints] = {0.0};
    //double erryh[nPoints] = {0.0};
    //for(int i = 3; i < 6; i++)
    //{
    //  for(int ip = 0; ip < nPoints; ip++)
    //  {  
    //    double x,y;
    //    g_temp[i]->GetPoint(j,x,y);
    //    cout << "i = " << i << ", j = " << j << ", x = " << x << "< y = " << y << endl;
    //    xval[j]  += x*y; 
    //    sum[j]   += y; 
    //    errxl[j] += (g_temp[i]->GetErrorXlow(j)*g_temp[i]->GetErrorXlow(j));
    //    errxh[j] += (g_temp[i]->GetErrorXhigh(j)*g_temp[i]->GetErrorXhigh(j));
    //    erryl[j] += (g_temp[i]->GetErrorYlow(j)*g_temp[i]->GetErrorYlow(j));
    //    erryh[j] += (g_temp[i]->GetErrorYhigh(j)*g_temp[i]->GetErrorYhigh(j));
    //  }
    //}
    //for(int ip = 0; ip < nPoints; ip++)
    //{
    //  xval[j] /= sum[j];
    //  errxl[j] = TMath::Sqrt(errxl[j]); 
    //  errxh[j] = TMath::Sqrt(errxh[j]);
    //  erryl[j] = TMath::Sqrt(erryl[j]);
    //  erryh[j] = TMath::Sqrt(erryh[j]); 
    //  g_spec->SetPoint(ip,xval[j],sum[j]);
    //  g_spec->SetPointError(ip,errxl[j],errxh[j],erryl[j],erryh[j]);
    //}  
  }

  TF1 *f_Levy = new TF1("f_Levy",LevyKStar,vmsa::ptMin,vmsa::ptMax,3);
  f_Levy->SetParameter(0,10);
  f_Levy->SetParameter(1,1000);
  f_Levy->SetParameter(2,0.05);
  f_Levy->SetLineStyle(2);
  f_Levy->SetLineColor(4);
  f_Levy->SetLineWidth(2);
  g_spec->Fit(f_Levy,"N");


//  TF1 *f_spec = new TF1("f_spec",pTLevy,vmsa::ptMin,vmsa::ptMax,3);
  TF1 *f_spec = new TF1("f_spec",pTLevyKStar,pt_set[pt_bin], pt_set[pt_bin+1], 3);
  f_spec->SetParameter(0,f_Levy->GetParameter(0));
  f_spec->SetParameter(1,f_Levy->GetParameter(1));
  f_spec->SetParameter(2,f_Levy->GetParameter(2));
  f_spec->SetLineStyle(2);
  f_spec->SetLineColor(2);
  f_spec->SetLineWidth(2);



/*
  TF1 *f_Expo = new TF1("f_Expo",Expo,vmsa::ptMin,vmsa::ptMax,2);
//  TF1 *f_Expo = new TF1("f_Expo","gaus(0)",vmsa::ptMin,vmsa::ptMax);
  f_Expo->SetParameter(0,10);
  f_Expo->SetParameter(1,1);
  g_spec->Fit(f_Expo,"N");

  TF1 *f_spec = new TF1("f_spec",Expo,0.01,5.1,2);
//  TF1 *f_spec = new TF1("f_spec","gaus(0)",0.01,5.1);
  f_spec->SetParameter(0,f_Expo->GetParameter(0));
  f_spec->SetParameter(1,f_Expo->GetParameter(1));
  f_spec->SetParameter(2,f_Expo->GetParameter(2));
*/


  g_spec->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  g_spec->GetXaxis()->SetLabelSize(0.05);
  g_spec->GetXaxis()->SetTitleSize(0.05);
  g_spec->GetXaxis()->SetTitleOffset(1.0);
  g_spec->GetYaxis()->SetTitle("d^{2}N/2#pip_{T}dp_{T}dy");
  g_spec->GetYaxis()->SetLabelSize(0.05);
  g_spec->GetYaxis()->SetTitleSize(0.05);
  g_spec->GetYaxis()->SetTitleOffset(1.1);
  g_spec->Draw("ap");
  f_Levy->Draw("same");
  c1->SetLogy();
  c1->SaveAs("pt.pdf");


  
  TCanvas *c_spec = new TCanvas("c_spec","c_spec",10,10,800,800);
  c_spec->cd()->SetBottomMargin(0.15);
  c_spec->cd()->SetTicks(1,1);
  c_spec->cd()->SetGrid(0,0);
  c_spec->SetLogy();
  TH1F *h_spec = new TH1F("h_spec","h_spec",100,0.0,10.0);
  for(int i_bin = 1; i_bin < 101; ++i_bin)
  {
    h_spec->SetBinContent(i_bin,-10.0);
    h_spec->SetBinError(i_bin,1.0);
  }
  h_spec->SetTitle("");
  h_spec->SetStats(0);
  h_spec->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_spec->GetXaxis()->CenterTitle();
  h_spec->GetYaxis()->SetTitle("dN/p_{T}dp_{T}");
  h_spec->GetYaxis()->CenterTitle();
  h_spec->GetYaxis()->SetRangeUser(1E-6,10);
  h_spec->Draw("pE");
  g_spec->Draw("pE same");
  f_Levy->Draw("l same");
  f_spec->Draw("l same");
  
  TCanvas *c10 = new TCanvas("c_ptdist","c_ptdist",10,10,800,800);
  c10->cd();
  c10->SetFillColor(0);
  c10->SetGrid(0,0);
  c10->SetTitle(0);
  c10->SetBottomMargin(0.15);
  c10->SetLeftMargin(0.15);
  f_spec->Draw();
  f_spec->SaveAs("indivitualptspec.pdf");

  return f_spec;
}

void write()
{
  File_OutPut->cd();
  McKStarMeson->Write();
  h_Tracks->Write();
  h_phiRP->Write();
  File_OutPut->Close();
}
