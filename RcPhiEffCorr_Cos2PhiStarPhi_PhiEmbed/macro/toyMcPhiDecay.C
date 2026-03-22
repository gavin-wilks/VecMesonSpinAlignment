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
#include "StRoot/Utility/v2/phiv2_14GeV.h";

using namespace std;

typedef std::map<std::string,TH1D*> TH1DMap;
typedef std::map<std::string,TF1*> TF1Map;

void readEfficiency(int energy, int year, int cut, const char *jobID);
void readTofEff(int energy);
void readTofEffFit(int energy, int centrality);
void getKinematics(TLorentzVector& lPhi, double const mass, int energy, int centrality);
void setDecayChannels(int const pid);
void decayAndFill(int cent, int const kf, TLorentzVector* lPhi, TClonesArray& daughters);
void fill(int cent, int const kf, TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus);
bool tpcReconstructed(int iParticleIndex, int cent, float Psi2, TLorentzVector const& lKaon);
bool tofReconstructed(int iParticleIndex, int cent, float Psi2, TLorentzVector const& lKaon);
void findHist(TLorentzVector const& lKaon, int iParticleIndex, float Psi2, int& EtaBin, int& PhiBin); // iParticleIndex = 0 => K+
float AngleShift(float phi);
void findHist_ToF(TLorentzVector const& lKaon, int iParticleIndex, float Psi2, int& EtaBin, int& PhiBin); // iParticleIndex = 0 => K+
TF1* readv2(int energy, int pid, int cent);
TF1* readspec(int energy, int pid, int cent);
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

// sampling functions and histograms
TF1 *f_v2, *f_spec, *f_flow, *f_y, *f_mass;
TH3F *h_Tracks;
TH2F *h_phiRP;

TFile *File_OutPut;

const Double_t v2SysErr010[10]     = {0.018410215, 0.002158695, 0.001272407, 0.005236625, 0.009963413, 0.010982662, 0.01098977, 0.020200927, 0.048425845};
const Double_t v2SysErr1040[10]    = {0.006775937, 0.001888813, 0.000728774, 0.000981321, 0.002251885, 0.002495498, 0.010080815, 0.002546398, 0.010536196};
const Double_t v2SysErr4080[10]    = {0.009521688, 0.008394166, 0.003205859, 0.010239825, 0.009919353, 0.019296126, 0.017187247, 0.01696486, 0.092181837};

int mMode;
int tableNum[5][9] = {{6,6,5,5,4,3,2,1,1},
                      {0,0,0,0,0,0,0,0,0},
                      {12,12,11,11,10,9,8,7,7},
                      {0,0,0,0,0,0,0,0,0},
                      {18,18,17,17,16,15,14,13,13}};

int tableNumV2[5][9] = {{0,0,0,0,0,0,0,0,0},
                        {0,0,0,0,0,0,0,0,0},
                        {0,0,0,0,0,0,0,0,0},
                        {0,0,0,0,0,0,0,0,0},
                        {239,239,239,239,141,141,141,43,43}};

double v2_7GeV[9] = {0.125818,0.125818,0.125818,0.125818,0.039951,0.039951,0.039951,0.013462,0.013462};

int mPtBin = 0;
double pt_set[5] = {1.2,1.8,2.4,3.2,4.2};

double pt_low[5]  = {1.0,2.0,3.0,2.0,1.0};
double pt_high[5] = {3.0,4.0,5.0,5.0,5.0};

int mFitMode = 0;
int mCentrality;

double mPsiRP = 0.0;
double mPsi2 = 0.0;
//double mPsi2RP = 0.0;
double mRes = 0.0;
double mChi = 0.0;
double mDelta = 0.0;

double mPsi1 = 0.0;
//double mPsi1RP = 0.0;
double mRes1 = 0.0;
double mChi1 = 0.0;
double mDelta1 = 0.0;

TF1 *f_pDel;
TF1 *f_pDel1;
TF1 *f_res;
TF1 *f_res1;

TFile *mInPutFile_Res;
TFile *mInPutFile_Res1;

int mNMax = 0;
double cos2delta = 0;
double cosdelta1 = 0;
double cos2delta1 = 0;

void toyMcPhiDecay(const int energy = 0, int ptbin = 4, const int centrality = 5, const int pid = 0, const int NMax = 10000, int mode = 1, int fitmode = 0, const char* jobID = "20221209") 
{
  mNMax = NMax;
  // mode = 0 --> pT depedence (20-60%), mode = 1 --> centrality depedence (1.0 < pT < 5.0) GeV/c
  // fitmode = 0 --> fit to plateau, fitmode = 1 --> fit to eta1
  mPtBin = ptbin;
  mMode = mode;
  mCentrality = centrality;

  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();
  gRandom->SetSeed();
  readEfficiency(energy,0,0,jobID);
  readTofEff(energy);
  readTofEffFit(energy,centrality);


  TString InPutFile_Res = Form("Utility/Resolution/file_%s_Resolution.root",vmsa::mBeamEnergy[energy].c_str());
  mInPutFile_Res = TFile::Open(InPutFile_Res.Data());

  TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get("p_mRes2_Sub");
  Float_t Res_raw = p_res2->GetBinContent(p_res2->FindBin(mCentrality));
  mRes = TMath::Sqrt(Res_raw);
  cout << "Resolution2 = " << mRes << endl;
  f_res = new TF1("resolution",EventPlaneResolution,0,80,0);
  mChi = f_res->GetX(mRes);

  cout << "mChi = " << mChi << endl;
  f_pDel = new TF1("deltaPsi",EventPlaneDist,-TMath::Pi()/2.0,TMath::Pi()/2.0,2);
  f_pDel->FixParameter(0,mChi);
  f_pDel->FixParameter(1,1.0/(2.0*TMath::Pi()));


  TString InPutFile_Res1 = Form("Utility/EpdResolution/Resolution_file_%s_EpdCorrections_4.root",vmsa::mBeamEnergy[energy].c_str());
  if(energy == 3 || energy == 0) InPutFile_Res1 = Form("Utility/EpdResolution/Resolution_file_%s_EpdCorrections_6.root",vmsa::mBeamEnergy[energy].c_str());
  mInPutFile_Res1 = TFile::Open(InPutFile_Res1.Data());

  TProfile *p_res1 = (TProfile*)mInPutFile_Res1->Get("AveCosDeltaPsi1");
  Float_t Res_raw1 = p_res1->GetBinContent(p_res1->FindBin(mCentrality));
  mRes1 = TMath::Sqrt(Res_raw1);
  cout << "Resolution1 = " << mRes1 << endl;
  f_res1 = new TF1("resolution1",EventPlaneResolution,0,80,0);
  mChi1 = f_res->GetX(mRes1); // This is for sub  event plane resolution
  mChi1 *= TMath::Sqrt(2.0);  // This is for full event plane resolution
  cout << "mChi1 = " << mChi1 << endl;
  f_pDel1 = new TF1("deltaPsi1",EventPlaneDist1st,-TMath::Pi(),TMath::Pi(),2);
  f_pDel1->FixParameter(0,mChi1);
  f_pDel1->FixParameter(1,1.0/(2.0*TMath::Pi()));
  


  //TCanvas *c_del = new TCanvas("c_del","c_del",10,10,800,800);
  //c_del->cd()->SetLeftMargin(0.15);
  //c_del->cd()->SetBottomMargin(0.15);
  //c_del->cd()->SetTicks(1,1);
  //c_del->cd()->SetGrid(0,0);
  //TH1F *h_del = new TH1F("h_del","h_del",100,-10.0,10.0);
  //for(int i_bin = 1; i_bin < 101; ++i_bin)
  //{
  //  h_del->SetBinContent(i_bin,-10.0);
  //  h_del->SetBinError(i_bin,1.0);
  //}
  //h_del->GetXaxis()->SetRangeUser(-TMath::Pi()/2.0,TMath::Pi()/2.0);
  //h_del->SetTitle("");
  //h_del->SetStats(0);
  //h_del->GetXaxis()->SetTitle("#Delta");
  //h_del->GetXaxis()->CenterTitle();
  //h_del->GetYaxis()->SetTitle("arb. units");
  //h_del->GetYaxis()->CenterTitle();
  //h_del->GetYaxis()->SetRangeUser(0.0,1.0);
  //h_del->Draw("pE");
  ////g_del->Draw("pE same");
  //f_pDel->Draw("l same");
  //c_del->SaveAs("del.pdf");


  // v2 & spectra implementation
  if(energy != 0) f_v2   = readv2(energy,pid,centrality);
  f_spec = readspec(energy,pid,centrality);
  f_flow = new TF1("f_flow",flowSample,-TMath::Pi(),TMath::Pi(),1);
  //f_y = new TF1("f_y","exp(-0.2375*x*x)",-5.0,5.0);
  f_y = new TF1("f_y","exp(-0.95*x*x)",-5.0,5.0);
  //f_mass = new TF1("f_mass",BreitWigner,vmsa::InvMass[pid]-0.00426*2,vmsa::InvMass[pid]+0.00426*2,3);
  //f_mass->SetParameter(0,vmsa::InvMass[pid]);
  //f_mass->SetParameter(1,0.00426);
  //f_mass->SetParameter(2,1.0);

  h_Tracks = new TH3F("h_Tracks","h_Tracks",20,vmsa::ptMin,vmsa::ptMax,vmsa::BinY,-1.0,1.0,36,-TMath::Pi(),TMath::Pi());
  h_phiRP = new TH2F("h_phiRP","h_phiRP",20,vmsa::ptMin,vmsa::ptMax,36,-TMath::Pi(),TMath::Pi()); // QA histogram for v2 sample

  pydecay = TPythia6Decayer::Instance();
  pydecay->Init();
  setDecayChannels(pid); // phi--> K+K-

  TClonesArray ptl("TParticle", 10);
  TLorentzVector *lPhi = new TLorentzVector();
  for(int i_ran = 0; i_ran < NMax; ++i_ran)
  {
    if (floor(10.0*i_ran/ static_cast<float>(NMax)) > floor(10.0*(i_ran-1)/ static_cast<float>(NMax)))
    cout << "=> processing data: " << 100.0*i_ran/ static_cast<float>(NMax) << "%" << endl;

    //double invariantmass = f_mass->GetRandom(vmsa::InvMass[pid]-0.012,vmsa::InvMass[pid]+0.012);
    //getKinematics(*lPhi,invariantmass);
    getKinematics(*lPhi,vmsa::InvMass[pid],energy,centrality);
    // if (fabs(lPhi->Phi()) > TMath::Pi()) continue;
    decayAndFill(centrality,vmsa::decayMother[pid],lPhi,ptl);

    if (i_ran % 1000 == 1) McPhiMeson->AutoSave("SaveSelf");
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  write();

  stopWatch->Stop();   
  stopWatch->Print();
}

void getKinematics(TLorentzVector& lPhi, double const mass, int energy, int centrality)
{
  //double const pt = gRandom->Uniform(vmsa::ptMin, vmsa::ptEffMax); // sample with flat distribution
  if (mMode == 0 ) double const pt = f_spec->GetRandom(pt_set[mPtBin], pt_set[mPtBin+1]); // sample with measured spectra
  if (mMode == 1 ) double const pt = f_spec->GetRandom(pt_low[mPtBin], pt_high[mPtBin]); // sample with measured spectra
  //double const y = f_y->GetRandom(-1.0,1.0);//Uniform(-vmsa::acceptanceRapidity, vmsa::acceptanceRapidity);
  double const y = gRandom->Uniform(-1.0,1.0);//Uniform(-vmsa::acceptanceRapidity, vmsa::acceptanceRapidity);
  // double const phi = TMath::TwoPi() * gRandom->Rndm(); // sample flat distribution
  //double const phi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
 
  //double Psi2 = gRandom->Uniform(-TMath::Pi()/2.0,TMath::Pi()/2.0); 
  double Psi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
  mPsiRP = Psi;

  double delta = f_pDel->GetRandom();
  mDelta = delta;

  mPsi2 = Psi + delta; // EP, the RP is smeared by delta
  //mPsi2RP = Psi2;       // RP

  //Angle wrapping
  while(mPsi2 > 0.5*TMath::Pi())
  {
    mPsi2 -= TMath::Pi();
  }
  while(mPsi2 < -0.5*TMath::Pi())
  {
    mPsi2 += TMath::Pi();
  }

  double delta1 = f_pDel1->GetRandom();
  mDelta1 = delta1;

  mPsi1 = Psi + delta1; // EP, the RP is smeared by delta
  //mPsi1RP = Psi1;       // RP

  //Angle wrapping
  while(mPsi1 > TMath::Pi())
  {
    mPsi1 -= 2.0*TMath::Pi();
  }
  while(mPsi1 < -TMath::Pi())
  {
    mPsi1 += 2.0*TMath::Pi();
  }

  //mPsi2 = gRandom->Uniform(-0.5*TMath::Pi(),0.5*TMath::Pi()); // random event plane angle
  //f_flow->ReleaseParameter(0);
  //if(energy != 0) f_flow->SetParameter(0,f_v2->Eval(pt));
  //if(energy == 0) f_flow->SetParameter(0,v2_7GeV[centrality]);
  //double const phi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
  //double const phi = f_flow->GetRandom() + mPsiRP; // sample with measured v2
  //if(phi > TMath::Pi())  phi -= 2.0*TMath::Pi();
  //if(phi < -TMath::Pi()) phi += 2.0*TMath::Pi();
  double const phi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
  double const mT = sqrt(mass * mass + pt * pt);
  double const pz = mT * sinh(y);
  double const E = mT * cosh(y);

  lPhi.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);

  //cos2delta += TMath::Cos(;
  //cosdelta1 = 0;
  //cos2delta1 = 0;
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

void decayAndFill(int centrality, int const kf, TLorentzVector* lPhi, TClonesArray& daughters)
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

  fill(centrality,kf,lPhi,lKplus,lKminus);
}

void fill(int centrality, int const kf, TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus)
{
  //int const centrality = floor(vmsa::NCentMax * gRandom->Rndm());
  //float const Psi2 = gRandom->Uniform(-0.5*TMath::Pi(),0.5*TMath::Pi()); // random event plane angle
  //float const Psi2 = 0.0; // fixed event plane angle
  // cout << "centrality = " << centrality << ", Psi2 = " << Psi2 << endl;
  // int const centrality = floor(2 * gRandom->Rndm());
  TLorentzVector lRcPhi = lKplus + lKminus; // phi meson reconstruction
  // cout << "lPhi.pt = " << lPhi->Pt() << ", lPhi.eta = " << lPhi->Eta() << ", lPhi.phi = " << lPhi->Phi() << ", lPhi.m = " << lPhi->M() << endl;
  // cout << "lRcPhi.pt = " << lRcPhi.Pt() << ", lRcPhi.eta = " << lRcPhi.Eta() << ", lRcPhi.phi = " << lRcPhi.Phi() << ", lRcPhi.m = " << lRcPhi.M() << endl;
  // cout << endl;

  float arr[110];
  int iArr = 0;
  //cout << "Centrality = " << centrality << endl;
  arr[iArr++] = centrality; // McPhi
  arr[iArr++] = mPsiRP; 
  arr[iArr++] = mPsi1; 
  arr[iArr++] = mPsi2; 
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
  arr[iArr++] = tpcReconstructed(0,centrality,mPsi2,lKplus);
  arr[iArr++] = tofReconstructed(0,centrality,mPsi2,lKplus);

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
  arr[iArr++] = tpcReconstructed(1,centrality,mPsi2,lKminus);
  arr[iArr++] = tofReconstructed(1,centrality,mPsi2,lKminus);

  arr[iArr++] = lRcPhi.Pt();
  arr[iArr++] = lRcPhi.P();
  arr[iArr++] = lRcPhi.PseudoRapidity();
  arr[iArr++] = lRcPhi.Rapidity();
  arr[iArr++] = lRcPhi.Phi();
  arr[iArr++] = lRcPhi.M();

  McPhiMeson->Fill(arr);
  // if(lRcPhi.Pt() < 10e-4) cout << "lPhi->Pt = " << lPhi->Pt() << ", lRcPhi.Pt = " << lRcPhi.Pt() << endl;

  h_Tracks->Fill(lPhi->Pt(),lPhi->Rapidity(),lPhi->Phi()); // fill QA histograms
  h_phiRP->Fill(lPhi->Pt(),lPhi->Phi());
}

//bool tpcReconstructed(int iParticleIndex, int cent, float Psi2, TLorentzVector const& lKaon)
//{
//   if(fabs(lKaon.Eta()) >= vmsa::mEtaMax) return false;
//
//   TH1D *h_TPC = NULL;
//   int EtaBin_TPC = -1;
//   int PhiBin_TPC = -1;
//   findHist(lKaon,iParticleIndex,Psi2,EtaBin_TPC,PhiBin_TPC);
//
//   TH1D *h_ToF = NULL;
//   TF1 *f_ToF = NULL;
//   int EtaBin_ToF = -1;
//   int PhiBin_ToF = -1;
//  
//   if(fabs(lKaon.Eta()) <= 1.0) findHist_ToF(lKaon,iParticleIndex,Psi2,EtaBin_ToF,PhiBin_ToF);
//   //findHist_ToF(lKaon,iParticleIndex,EtaBin_ToF,PhiBin_ToF);
//
//   // cout << KEY.c_str() << endl;
//   if (iParticleIndex == 0)
//   {
//     string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
//     //string KEY_TPC = Form("h_mEff_Cent_9_Eta_%d_Phi_%d",EtaBin_TPC,PhiBin_TPC); // get TPC eff @ 20-60%
//     h_TPC = h_EffKplus[KEY_TPC];
//
//     if(fabs(lKaon.Eta()) <= 1.0)
//     {
//       string KEY_ToF = Form("h_mEfficiency_Kplus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
//       //string KEY_ToF = Form("h_mEfficiency_Kplus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
//       h_ToF = h_TofKplus[KEY_ToF];
//
//       string KEY_ToFFit; 
//       if(mFitMode == 0) KEY_ToFFit = Form("f_mToFMatch_Kplus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
//       if(mFitMode == 1) KEY_ToFFit = Form("f_mToFMatch_Kplus_Cent_%d_Eta_%d",cent,EtaBin_ToF); // get ToF eff with eta only @ 20-60%
//       f_ToF = f_TofKplus[KEY_ToFFit]; // only 20-60%
//     }
//   }
//   else
//   {
//     string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
//     //string KEY_TPC = Form("h_mEff_Cent_9_Eta_%d_Phi_%d",EtaBin_TPC,PhiBin_TPC); // get TPC eff @ 20-60%
//     h_TPC = h_EffKminus[KEY_TPC];
//
//     if(fabs(lKaon.Eta()) <= 1.0)
//     {
//       string KEY_ToF = Form("h_mEfficiency_Kminus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
//       //string KEY_ToF = Form("h_mEfficiency_Kminus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
//       h_ToF = h_TofKminus[KEY_ToF];
//
//       string KEY_ToFFit;
//       if(mFitMode == 0) KEY_ToFFit = Form("f_mToFMatch_Kminus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
//       if(mFitMode == 1) KEY_ToFFit = Form("f_mToFMatch_Kminus_Cent_%d_Eta_%d",cent,EtaBin_ToF); // get ToF eff with eta only @ 20-60%
//       f_ToF = f_TofKminus[KEY_ToFFit]; // only 20-60%
//     }
//   }
//
//   double pt = lKaon.Perp();
//   int const bin_TPC = h_TPC->FindBin(pt);
//   bool is_TPC = gRandom->Rndm() < h_TPC->GetBinContent(bin_TPC);
//
//   bool is_ToF = false;
//   if(fabs(lKaon.Eta()) <= 1.0)
//   { 
//     int const bin_ToF = h_ToF->FindBin(pt); // tof fit and hist combined
//     double prob_tof = 0.0;
//     if(pt > 0.3) prob_tof = f_ToF->Eval(pt); // donot use fit extrapolation
//     else prob_tof = h_ToF->GetBinContent(bin_ToF);
//     is_ToF = gRandom->Rndm() < prob_tof;
//   }
//
//   // int const bin_ToF = h_ToF->FindBin(pt); // tof hist only
//   // bool is_ToF = gRandom->Rndm() < h_ToF->GetBinContent(bin_ToF);
//
//   // bool is_ToF = gRandom->Rndm() < f_ToF->Eval(pt); // tof fit only
//
//   // cout << "is_TPC: " << is_TPC << ", is_ToF: " << is_ToF << ", is_TPC && is_ToF: " << (is_TPC && is_ToF) << endl;
//   // cout << "pt = " << pt << ", h_ToF = " << h_ToF->GetBinContent(bin_ToF) << ", f_ToF = " << f_ToF->Eval(pt) << ", prob_tof = " << prob_tof << endl; // comparison between eff from hist and func
//
//   if(fabs(lKaon.Eta()) <= 1.0) return is_TPC && is_ToF;
//   if(fabs(lKaon.Eta()) > 1.0) return is_TPC;
//   
//}

bool tpcReconstructed(int iParticleIndex, int cent, float Psi2, TLorentzVector const& lKaon)
{
   if(fabs(lKaon.Eta()) >= vmsa::mEtaMax) return false;

   TH1D *h_TPC = NULL;
   int EtaBin_TPC = -1;
   int PhiBin_TPC = -1;
   findHist(lKaon,iParticleIndex,Psi2,EtaBin_TPC,PhiBin_TPC);

   TH1D *h_ToF = NULL;
   TF1 *f_ToF = NULL;
   int EtaBin_ToF = -1;
   int PhiBin_ToF = -1;
  
   if (iParticleIndex == 0)
   {
     string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
     //string KEY_TPC = Form("h_mEff_Cent_9_Eta_%d_Phi_%d",EtaBin_TPC,PhiBin_TPC); // get TPC eff @ 20-60%
     h_TPC = h_EffKplus[KEY_TPC];
   }
   else
   {
     string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
     //string KEY_TPC = Form("h_mEff_Cent_9_Eta_%d_Phi_%d",EtaBin_TPC,PhiBin_TPC); // get TPC eff @ 20-60%
     h_TPC = h_EffKminus[KEY_TPC];
   }

   double pt = lKaon.Perp();
   int const bin_TPC = h_TPC->FindBin(pt);
   bool is_TPC = gRandom->Rndm() < h_TPC->GetBinContent(bin_TPC);

   return is_TPC;
   
}
bool tofReconstructed(int iParticleIndex, int cent, float Psi2, TLorentzVector const& lKaon)
{
   if(fabs(lKaon.Eta()) >= vmsa::mEtaMax) return false;

   TH1D *h_ToF = NULL;
   TF1 *f_ToF = NULL;
   int EtaBin_ToF = -1;
   int PhiBin_ToF = -1;
  
   if(fabs(lKaon.Eta()) <= 1.0) findHist_ToF(lKaon,iParticleIndex,Psi2,EtaBin_ToF,PhiBin_ToF);

   if (iParticleIndex == 0)
   {
     if(fabs(lKaon.Eta()) <= 1.0)
     {
       string KEY_ToF = Form("h_mEfficiency_Kplus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
       //string KEY_ToF = Form("h_mEfficiency_Kplus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
       h_ToF = h_TofKplus[KEY_ToF];

       string KEY_ToFFit; 
       if(mFitMode == 0) KEY_ToFFit = Form("f_mToFMatch_Kplus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
       if(mFitMode == 1) KEY_ToFFit = Form("f_mToFMatch_Kplus_Cent_%d_Eta_%d",cent,EtaBin_ToF); // get ToF eff with eta only @ 20-60%
       f_ToF = f_TofKplus[KEY_ToFFit]; // only 20-60%
     }
   }
   else
   {
     if(fabs(lKaon.Eta()) <= 1.0)
     {
       string KEY_ToF = Form("h_mEfficiency_Kminus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
       //string KEY_ToF = Form("h_mEfficiency_Kminus_Cent_9_Eta_%d",EtaBin_ToF); // get ToF eff with eta only @ 20-60%
       h_ToF = h_TofKminus[KEY_ToF];

       string KEY_ToFFit;
       if(mFitMode == 0) KEY_ToFFit = Form("f_mToFMatch_Kminus_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_ToF,PhiBin_ToF); // get ToF eff @ 20-60%
       if(mFitMode == 1) KEY_ToFFit = Form("f_mToFMatch_Kminus_Cent_%d_Eta_%d",cent,EtaBin_ToF); // get ToF eff with eta only @ 20-60%
       f_ToF = f_TofKminus[KEY_ToFFit]; // only 20-60%
     }
   }

   double pt = lKaon.Perp();

   bool is_ToF = false;
   if(fabs(lKaon.Eta()) <= 1.0)
   { 
     int const bin_ToF = h_ToF->FindBin(pt); // tof fit and hist combined
     double prob_tof = 0.0;
     if(pt > 0.3) prob_tof = f_ToF->Eval(pt); // donot use fit extrapolation
     else prob_tof = h_ToF->GetBinContent(bin_ToF);
     is_ToF = gRandom->Rndm() < prob_tof;
   }

   if(fabs(lKaon.Eta()) < 1.0) return is_ToF;
}

void findHist(TLorentzVector const& lKaon, int iParticleIndex, float Psi2, int& EtaBin, int& PhiBin)
{
  float eta = lKaon.Eta();
  EtaBin = h_FrameEta[iParticleIndex]->FindBin(eta)-1;

  //float phi = lKaon.Phi()-Psi2;
  //float phi_shift = AngleShift(phi);
  //PhiBin = h_FramePhi[iParticleIndex]->FindBin(phi_shift)-1;

  float phi = lKaon.Phi();
  PhiBin = h_FramePhi[iParticleIndex]->FindBin(phi)-1;
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

void findHist_ToF(TLorentzVector const& lKaon, int iParticleIndex, float Psi2, int& EtaBin, int& PhiBin)
{
  float eta = lKaon.Eta();
  EtaBin = h_FrameEta_ToF[iParticleIndex]->FindBin(eta)-1;
  // cout << "eta = " << eta << ", EtaBin = " << EtaBin << endl;

  //float phi = lKaon.Phi()-Psi2;
  //float phi_shift = AngleShift(phi);
  //PhiBin = h_FramePhi_ToF[iParticleIndex]->FindBin(phi_shift)-1;

  float phi = lKaon.Phi();
  PhiBin = h_FramePhi_ToF[iParticleIndex]->FindBin(phi)-1;
  // cout << "phi = " << phi << ", PhiBin = " << PhiBin << endl;
}

void readEfficiency(int energy, int year, int cut, const char *jobID)
{
  // string inputKplus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[0].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  string inputKplus = Form("/star/data01/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/TPC/Eff_%s_%s.root",vmsa::mParType[0].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Kplus = TFile::Open(inputKplus.c_str());
  cout << "OPEN Efficiency File for K+: " << inputKplus.c_str() << endl;

  // string inputKminus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[1].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  string inputKminus = Form("/star/data01/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/TPC/Eff_%s_%s.root",vmsa::mParType[1].c_str(),vmsa::mBeamEnergy[energy].c_str());
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
  string outputfile = Form("./Eff_%s_SingleKaon_cent%d_pt%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),mCentrality,mPtBin,jobID);
  cout << "OutPut File set to: " << outputfile.c_str() << endl;
  File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();

  int BufSize = (int)pow(2., 16.);
  // int Split = 1;

  const char* varlist = "Centrality:PsiRP:Psi1:Psi2:McPt:McP:McEta:McY:McPhi:McInvMass:McPid:" // MC phi 
                        "KpMcPt:KpMcEta:KpMcY:KpMcPhi:KpMcM:KpMcPid:" // MC K+ 
                        "KpRcPt:KpRcEta:KpRcY:KpRcPhi:KpRcM:KpRcTpc:KpRcTof:" // RC K+
                        "KmMcPt:KmMcEta:KmMcY:KmMcPhi:KmMcM:KmMcPid:" // MC K-
                        "KmRcPt:KmRcEta:KmRcY:KmRcPhi:KmRcM:KmRcTpc:KmRcTof:" // RC K-
                        "RcPt:RcP:RcEta:RcY:RcPhi:RcInvMass"; // reconstructed phi

  McPhiMeson = new TNtuple("McPhiMeson", "McPhiMeson", varlist, BufSize);
}

void readTofEff(int energy)
{
  // string inputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/ToF/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  cout << "OPEN Efficiency File for K+ and K-: " << inputfile.c_str() << endl;

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
      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
	// string KEY;
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

void readTofEffFit(int energy, int i_cent)
{
  // string inputKplus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/FitPar_AuAu%s_Kplus_first.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputKplus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/ToF/FitPar_AuAu%s_Kplus_cent%d.root",vmsa::mBeamEnergy[energy].c_str(),i_cent);
  TFile *File_Kplus = TFile::Open(inputKplus.c_str());
  cout << "OPEN ToF Matching Efficiency Fit File for K+: " << inputKplus.c_str() << endl;

  // string inputKminus = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/FitPar_AuAu%s_Kminus_first.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputKminus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/ToF/FitPar_AuAu%s_Kminus_cent%d.root",vmsa::mBeamEnergy[energy].c_str(),i_cent);
  TFile *File_Kminus = TFile::Open(inputKminus.c_str());
  cout << "OPEN ToF Matching Efficiency Fit File for K-: " << inputKminus.c_str() << endl;

  //for(int i_cent = 0; i_cent < 10; ++i_cent)
  //{
  for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta) // eta
  {
    TH1D *h_Kplus = NULL;
    string KEY;
    KEY = Form("h_mFitParameters_Kplus_Cent_%d_Eta_%d",i_cent,i_eta);
    h_Kplus = (TH1D*)File_Kplus->Get(KEY.c_str());
    // cout << "Kplus KEY: " << KEY.c_str() << endl;
    KEY = Form("f_mToFMatch_Kplus_Cent_%d_Eta_%d",i_cent,i_eta);
    f_TofKplus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_TofKplus[KEY]->FixParameter(i_par,h_Kplus->GetBinContent(i_par+1));
      // if(i_par < 2) cout << "Kplus: i_par = " << i_par << ", par = " << f_TofKplus[KEY]->GetParameter(i_par) << endl;
    }

    TH1D *h_Kminus = NULL;
    KEY = Form("h_mFitParameters_Kminus_Cent_%d_Eta_%d",i_cent,i_eta);
    h_Kminus = (TH1D*)File_Kminus->Get(KEY.c_str());
    // cout << "Kminus KEY: " << KEY.c_str() << endl;
    KEY = Form("f_mToFMatch_Kminus_Cent_%d_Eta_%d",i_cent,i_eta);
    f_TofKminus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_TofKminus[KEY]->FixParameter(i_par,h_Kminus->GetBinContent(i_par+1));
      // if(i_par < 2) cout << "Kminus: i_par = " << i_par << ", par = " << f_TofKminus[KEY]->GetParameter(i_par) << endl;
    }
  }

  for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta) // eta & phi
  {
    for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
    {
      TH1D *h_Kplus = NULL;
      string KEY;
      KEY = Form("h_mFitParameters_Kplus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
      h_Kplus = (TH1D*)File_Kplus->Get(KEY.c_str());
      // cout << "Kplus KEY: " << KEY.c_str() << endl;
      KEY = Form("f_mToFMatch_Kplus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
      f_TofKplus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
      for(int i_par = 0; i_par < 7; ++i_par)
      {
        f_TofKplus[KEY]->FixParameter(i_par,h_Kplus->GetBinContent(i_par+1));
        // if(i_par < 2) cout << "Kplus: i_par = " << i_par << ", par = " << f_TofKplus[KEY]->GetParameter(i_par) << endl;
      }

      TH1D *h_Kminus = NULL;
      KEY = Form("h_mFitParameters_Kminus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
      h_Kminus = (TH1D*)File_Kminus->Get(KEY.c_str());
      // cout << "Kminus KEY: " << KEY.c_str() << endl;
      KEY = Form("f_mToFMatch_Kminus_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
      f_TofKminus[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.1,10,7);
      for(int i_par = 0; i_par < 7; ++i_par)
      {
        f_TofKminus[KEY]->FixParameter(i_par,h_Kminus->GetBinContent(i_par+1));
        // if(i_par < 2) cout << "Kminus: i_par = " << i_par << ", par = " << f_TofKminus[KEY]->GetParameter(i_par) << endl;
      }
    }	
  }
  //}
}

//TF1* readv2(int energy, int pid, int centrality){
//  //string centFile[9] = {"4080","4080","4080","4080","1040","1040","1040","0010","0010"};
//  
//  string centlabel = "4080";
//  if(centrality >= 4 && centrality <= 6) centlabel = "1040";
//  if(centrality >= 7 && centrality <= 8) centlabel = "0010";
//
//  if(energy != 3)
//  {
//    string InPutV2 = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
//    if((energy == 4 || energy == 2 || energy == 0) && mMode == 1) InPutV2 = "/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/v2/HEPData-ins1395151-v2-root.root";
//    TFile *File_v2 = TFile::Open(InPutV2.c_str());
//    std::cout << "v2 file: " << InPutV2 << endl;
//
//    TGraphAsymmErrors *g_v2;
//
//    if(mMode == 0) g_v2 = (TGraphAsymmErrors*)File_v2->Get("g_v2");
//
//    if((energy == 4 || energy == 2 || energy == 0) && mMode == 1) 
//    {
//      TDirectory *dir = (TDirectory*) File_v2->Get(Form("Table %d",tableNumV2[energy][centrality]));
//      dir->cd(); 
//      g_v2 = (TGraphAsymmErrors*)dir->Get("Graph1D_y1");
//    }
//  }
//  if( energy == 3 && mMode == 1 )
//  {
//    int centidx = 0;
//    if(centrality >= 0 && centrality <= 3) centidx = 3;
//    if(centrality >= 4 && centrality <= 6) centidx = 2;
//    if(centrality >= 7 && centrality <= 8) centidx = 1;
//    g_v2 = new TGraphAsymmErrors();
//    for(int ipt = 0; ipt < phiv2_14::ptbins; ipt++)
//    {
//      g_v2->SetPoint(ipt, phiv2_14::pt[centidx][ipt], phiv2_14::v2[centidx][ipt]);
//      double stat = phiv2_14::stat[centidx][ipt];
//      double sys  = phiv2_14::sys[centidx][ipt];
//      double totalerr = TMath::Sqrt( stat*stat + sys*sys );
//      g_v2->SetPointError(ipt, 0.0, 0.0, totalerr, totalerr);  
//    }
//
//  }
//
//
//  TF1 *f_v2 = new TF1("f_v2",v2_pT_FitFunc,vmsa::ptMin,vmsa::ptMax,5);
//  f_v2->FixParameter(0,2);
//  f_v2->SetParameter(1,0.1);
//  f_v2->SetParameter(2,0.1);
//  f_v2->SetParameter(3,0.1);
//  f_v2->SetParameter(4,0.1);
//  f_v2->SetLineColor(2);
//  f_v2->SetLineWidth(2);
//  f_v2->SetLineStyle(2);
//  cout << "Fitting v2" << endl;
//  g_v2->Fit(f_v2,"N");
//
//  
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
//  c_v2->SaveAs(Form("v2_%s_cent%s.pdf",vmsa::mBeamEnergy[energy].c_str(),centlabel.c_str()));
//
//
//
//  return f_v2;
//}

TF1* readv2(int energy, int pid, int centrality){
  //string centFile[9] = {"4080","4080","4080","4080","1040","1040","1040","0010","0010"};
  
  string centlabel = "4080";
  if(centrality >= 4 && centrality <= 6) centlabel = "1040";
  if(centrality >= 7 && centrality <= 8) centlabel = "0010";

  if(energy != 3)
  {
    string InPutV2 = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
    if((energy == 2 || energy == 0) && mMode == 1) InPutV2 = "/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/v2/HEPData-ins1395151-v2-root.root";
    if(energy == 4 && mMode == 1) InPutV2 = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/v2/OutPhi_v2_Cent%s.root",centlabel.c_str());
    TFile *File_v2 = TFile::Open(InPutV2.c_str());
    std::cout << "v2 file: " << InPutV2 << endl;

    TGraphAsymmErrors *g_v2;

    if(mMode == 0) g_v2 = (TGraphAsymmErrors*)File_v2->Get("g_v2");

    if((energy == 2 || energy == 0) && mMode == 1) 
    {
      TDirectory *dir = (TDirectory*) File_v2->Get(Form("Table %d",tableNumV2[energy][centrality]));
      dir->cd(); 
      g_v2 = (TGraphAsymmErrors*)dir->Get("Graph1D_y1");
    }
    if(energy == 4 && mMode == 1) 
    {
      g_v2 = (TGraphAsymmErrors*) File_v2->Get("Graph");
      g_v2->Print();
      //for(int ipoint = 0; ipoint < g_v2->GetN(); ipoint++)
      //{
      //  double v2errStat = g_v2->GetErrorYhigh(ipoint);
      //  cout << v2errStat << endl;
      //  double v2errSys;
      //  if(centlabel == "0010") v2errSys = v2SysErr010[ipoint];  
      //  if(centlabel == "1040") v2errSys = v2SysErr1040[ipoint]; 
      //  if(centlabel == "4080") v2errSys = v2SysErr4080[ipoint]; 
      //  cout << v2errSys << endl;        
 
      //  double v2err = TMath::Sqrt(v2errStat*v2errStat+v2errSys*v2errSys);
      //  cout << v2err << endl;
      //  //g_v2->SetPointError(ipoint,0.0,0.0,v2err,v2err);
      //  g_v2->SetPointEYlow(ipoint,v2err);
      //  g_v2->SetPointEYhigh(ipoint,v2err);
      //  cout << ipoint << "    Set the point error" << endl;
      //  cout << "point error = " << g_v2->GetErrorYhigh(ipoint) << endl;
      //}  
    }
  }
  
  if( energy == 3 && mMode == 1 )
  {
    int centidx = 0;
    if(centrality >= 0 && centrality <= 3) centidx = 3;
    if(centrality >= 4 && centrality <= 6) centidx = 2;
    if(centrality >= 7 && centrality <= 8) centidx = 1;
    g_v2 = new TGraphAsymmErrors();
    for(int ipt = 0; ipt < phiv2_14::ptbins; ipt++)
    {
      g_v2->SetPoint(ipt, phiv2_14::pt[centidx][ipt], phiv2_14::v2[centidx][ipt]);
      double stat = phiv2_14::stat[centidx][ipt];
      double sys  = phiv2_14::sys[centidx][ipt];
      double totalerr = TMath::Sqrt( stat*stat + sys*sys );
      g_v2->SetPointError(ipt, 0.0, 0.0, totalerr, totalerr);  
    }

  }


  TF1 *f_v2 = new TF1("f_v2",v2_pT_FitFunc,vmsa::ptMin,vmsa::ptMax,5);
  f_v2->FixParameter(0,2);
  f_v2->SetParameter(1,0.1);
  f_v2->SetParameter(2,0.1);
  f_v2->SetParameter(3,0.1);
  f_v2->SetParameter(4,0.1);
  f_v2->SetLineColor(2);
  f_v2->SetLineWidth(2);
  f_v2->SetLineStyle(2);
  cout << "Fitting v2" << endl;
  g_v2->Fit(f_v2,"N");

  
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
  h_v2->GetYaxis()->SetRangeUser(0.0,0.35);
  h_v2->Draw("pE");
  g_v2->Draw("pE same");
  f_v2->Draw("l same");
  c_v2->SaveAs(Form("v2_%s_cent%s.pdf",vmsa::mBeamEnergy[energy].c_str(),centlabel.c_str()));



  return f_v2;
}


//TF1* readv2(int energy, int pid)
//{
//  string InPutV2 = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
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

TF1* readspec(int energy, int pid, int centrality)
{
  TCanvas *c1 = new TCanvas("c_pt","c_pt",10,10,800,800);
  c1->SetFillColor(0);
  c1->SetGrid(0,0);
  c1->SetTitle(0);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.15);
 
  double ptlow[5][11] = { {0.4,0.5,0.6,0.7,0.9,1.0,1.3,0.,0.,0.,0.}, 
                          {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                          {0.4,0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,0.},
                          {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                          {0.4,0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.0} };
 
  double pthigh[5][11] = { {0.5,0.6,0.7,0.9,1.0,1.3,1.7,0.,0.,0.,0.}, 
                           {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                           {0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.5,0.},
                           {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                           {0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.0,4.0} };

  double ptvalues[5][9][11] = { { {0.4541,0.5510,0.6500,0.7973,0.9490,1.1399,1.4818,0.,0.,0.,0.},
                                  {0.4541,0.5510,0.6500,0.7973,0.9490,1.1399,1.4818,0.,0.,0.,0.},
                                  {0.4498,0.5562,0.6518,0.8012,0.9497,1.1433,1.4818,0.,0.,0.,0.},
                                  {0.4498,0.5562,0.6518,0.8012,0.9497,1.1433,1.4818,0.,0.,0.,0.},
                                  {0.4437,0.5493,0.6523,0.8020,0.9498,1.1446,1.4839,0.,0.,0.,0.},
                                  {0.4430,0.5493,0.6521,0.8018,0.9498,1.1442,1.4834,0.,0.,0.,0.},                 
                                  {0.4449,0.5494,0.6531,0.8027,0.9500,1.1456,1.4857,0.,0.,0.,0.},
                                  {0.4461,0.5496,0.6543,0.8037,0.9502,1.1467,1.4876,0.,0.,0.,0.},  
                                  {0.4461,0.5496,0.6543,0.8037,0.9502,1.1467,1.4876,0.,0.,0.,0.} },
                                { {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},                 
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},  
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.} },
                                { {0.4500,0.5556,0.6518,0.8013,0.9497,1.1438,1.4835,1.8392,2.2188,2.8795,0.},
                                  {0.4500,0.5556,0.6518,0.8013,0.9497,1.1438,1.4835,1.8392,2.2188,2.8795,0.},
                                  {0.4500,0.5551,0.6516,0.8008,0.9496,1.1427,1.4813,1.8377,2.2137,2.8589,0.},
                                  {0.4500,0.5551,0.6516,0.8008,0.9496,1.1427,1.4813,1.8377,2.2137,2.8589,0.},
                                  {0.4437,0.5494,0.6523,0.8021,0.9498,1.1446,1.4840,1.8390,2.2165,2.8633,0.},
                                  {0.4441,0.5494,0.6525,0.8022,0.9499,1.1449,1.4845,1.8393,2.2173,2.8660,0.},                 
                                  {0.4459,0.5494,0.6539,0.8035,0.9501,1.1464,1.4871,1.8408,2.2214,2.8798,0.},
                                  {0.4457,0.5494,0.6535,0.8033,0.9501,1.1463,1.4868,1.8406,2.2210,2.8785,0.},  
                                  {0.4457,0.5494,0.6535,0.8033,0.9501,1.1463,1.4868,1.8406,2.2210,2.8785,0.} },
                                { {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},                 
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},  
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.} },
                                { {0.4500,0.5547,0.6515,0.8006,0.9496,1.1427,1.4814,1.8380,2.2150,2.7142,3.3666},
                                  {0.4500,0.5547,0.6515,0.8006,0.9496,1.1427,1.4814,1.8380,2.2150,2.7142,3.3666},
                                  {0.4434,0.5495,0.6523,0.8020,0.9498,1.1446,1.4845,1.8396,2.2191,2.7178,3.3767},
                                  {0.4434,0.5495,0.6523,0.8020,0.9498,1.1446,1.4845,1.8396,2.2191,2.7178,3.3767},
                                  {0.4455,0.5494,0.6535,0.8031,0.9500,1.1460,1.4864,1.8404,2.2203,2.7178,3.3718},
                                  {0.4453,0.5494,0.6531,0.8029,0.9500,1.1458,1.4860,1.8402,2.2198,2.7173,3.3699},                 
                                  {0.4463,0.5494,0.6543,0.8039,0.9502,1.1469,1.4879,1.8412,2.2227,2.7202,3.3798},
                                  {0.4457,0.5494,0.6535,0.8033,0.9501,1.1463,1.4869,1.8407,2.2211,2.7186,3.3745},  
                                  {0.4457,0.5494,0.6535,0.8033,0.9501,1.1463,1.4869,1.8407,2.2211,2.7186,3.3745} } };

  string centlabel = "6080";
  if(centrality >= 2 && centrality <= 3) centlabel = "4060";
  if(centrality >= 4 && centrality <= 4) centlabel = "3040";
  if(centrality >= 5 && centrality <= 5) centlabel = "2030";
  if(centrality >= 6 && centrality <= 6) centlabel = "1020";
  if(centrality >= 7 && centrality <= 8) centlabel = "0010";

  if(energy != 3)
  {
    string InPutSpec = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_Spec.root",vmsa::mBeamEnergy[energy].c_str());
    if((energy == 4 || energy == 2 || energy == 0) && mMode == 1) InPutSpec = "/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/pTspectra/HEPData-ins1378002-v1-root.root";
    TFile *File_Spec = TFile::Open(InPutSpec.c_str());
    cout << "Input spectra" << InPutSpec << endl;
 
    TGraphAsymmErrors *g_spec;
    if(mMode == 0) g_spec = (TGraphAsymmErrors*)File_Spec->Get("g_spec");
    if((energy == 4 || energy == 2 || energy == 0) && mMode == 1) 
    {
      TDirectory *dir = (TDirectory*) File_Spec->Get(Form("Table %d",tableNum[energy][centrality]));
      dir->cd(); 
      g_spec = (TGraphAsymmErrors*)dir->Get("Graph1D_y1");
      for(int i = 0; i < g_spec->GetN(); i++)
      {
        double x,y;
        g_spec->GetPoint(i,x,y);
        g_spec->SetPoint(i,ptvalues[energy][centrality][i],y);
        //g_spec->SetPointError(i,fabs(ptvalues[energy][centrality][i]-ptlow[energy][i]),fabs(ptvalues[energy][centrality][i]-pthigh[energy][i]),g_spec->GetErrorYlow(i),g_spec->GetErrorYhigh(i));
        g_spec->SetPointError(i,0.0,0.0,g_spec->GetErrorYlow(i),g_spec->GetErrorYhigh(i));
      }
    }

    TF1 *f_Levy = new TF1("f_Levy",Levy,vmsa::ptMin,vmsa::ptMax,3);
    f_Levy->SetParameter(0,1);
    if(energy == 4 && centrality >= 0 && centrality <= 1) f_Levy->SetParameter(0,0.01);
    if(energy == 0 ) f_Levy->SetParameter(0,0.01);
    f_Levy->SetParameter(1,10);
    f_Levy->SetParameter(2,0.1);
    f_Levy->SetLineStyle(2);
    f_Levy->SetLineColor(4);
    f_Levy->SetLineWidth(2);
    cout << "Fitting full pT distribution" << endl;
    g_spec->Fit(f_Levy,"N");


    //TF1 *f_spec = new TF1("f_spec",pTLevy,vmsa::ptMin,vmsa::ptMax,3);
    TF1 *f_spec = new TF1("f_spec",pTLevy,pt_low[mPtBin], pt_high[mPtBin], 3);
    f_spec->SetParameter(0,f_Levy->GetParameter(0));
    f_spec->SetParameter(1,f_Levy->GetParameter(1));
    f_spec->SetParameter(2,f_Levy->GetParameter(2));
    f_spec->SetLineStyle(2);
    f_spec->SetLineColor(2);
    f_spec->SetLineWidth(2);


  

/*    TF1 *f_spec = new TF1("f_spec",Expo,0.01,5.1,2);
//    TF1 *f_spec = new TF1("f_spec","gaus(0)",0.01,5.1);
    f_spec->SetParameter(0,f_Expo->GetParameter(0));
    f_spec->SetParameter(1,f_Expo->GetParameter(1));
    f_spec->SetParameter(2,f_Expo->GetParameter(2));
*/  



  }

  double dNdy14[9] = {0.009126,0.009126,0.042299,0.042299,0.095731,0.152881,0.223975,0.329114,0.329114};    
  double Texp14[9] = {0.211483,0.211483,0.228253,0.228253,0.241851,0.245252,0.266930,0.268719,0.268719};

  TCanvas *c10 = new TCanvas("c_ptdist","c_ptdist",10,10,800,800);
  c10->cd();
  c10->SetFillColor(0);
  c10->SetGrid(0,0);
  c10->SetTitle(0);
  c10->SetBottomMargin(0.15);
  c10->SetLeftMargin(0.15);
  TF1 *f_Expo = new TF1("f_Expo",specExp,vmsa::ptMin,vmsa::ptMax,2);
  if(energy == 3)
  {
    f_spec = new TF1("f_specExp",pTspecExp, pt_low[mPtBin], pt_high[mPtBin], 2);
    f_spec->SetParameter(0,dNdy14[centrality]);
    f_spec->SetParameter(1,Texp14[centrality]);
//    TF1 *f_Expo = new TF1("f_Expo","gaus(0)",vmsa::ptMin,vmsa::ptMax);
    f_Expo->SetParameter(0,dNdy14[centrality]);
    f_Expo->SetParameter(1,Texp14[centrality]);
  }
  f_spec->Draw();
  c10->SaveAs(Form("indivitualptspec_%s_cent%s.pdf",vmsa::mBeamEnergy[energy].c_str(),centlabel.c_str()));
 
  if(energy != 3)
  {
    c1->cd(); 
    c1->SetLogy();
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
    c1->SaveAs(Form("ptpec_%s_cent%s.pdf",vmsa::mBeamEnergy[energy].c_str(),centlabel.c_str()));
  }
  if(energy == 3)
  {
    c1->cd(); 
    c1->SetLogy();
    f_Expo->Draw();
    c1->SaveAs(Form("ptpec_%s_cent%s.pdf",vmsa::mBeamEnergy[energy].c_str(),centlabel.c_str()));
  }
  /*
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
  */


  return f_spec;
}

//`TF1* readspec(int energy, int pid)
//`{
//`  string InPutSpec = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_Spec.root",vmsa::mBeamEnergy[energy].c_str());
//`  TFile *File_Spec = TFile::Open(InPutSpec.c_str());
//`  TGraphAsymmErrors *g_spec = (TGraphAsymmErrors*)File_Spec->Get("g_spec");
//`  TF1 *f_Levy = new TF1("f_Levy",Levy,vmsa::ptMin,vmsa::ptMax,3);
//`  f_Levy->SetParameter(0,1);
//`  f_Levy->SetParameter(1,10);
//`  f_Levy->SetParameter(2,0.1);
//`  f_Levy->SetLineStyle(2);
//`  f_Levy->SetLineColor(4);
//`  f_Levy->SetLineWidth(2);
//`  g_spec->Fit(f_Levy,"N");
//`
//`  TF1 *f_spec = new TF1("f_spec",pTLevy,vmsa::ptMin,vmsa::ptMax,3);
//`  f_spec->SetParameter(0,f_Levy->GetParameter(0));
//`  f_spec->SetParameter(1,f_Levy->GetParameter(1));
//`  f_spec->SetParameter(2,f_Levy->GetParameter(2));
//`  f_spec->SetLineStyle(2);
//`  f_spec->SetLineColor(2);
//`  f_spec->SetLineWidth(2);
//`
//`  /*
//`  TCanvas *c_spec = new TCanvas("c_spec","c_spec",10,10,800,800);
//`  c_spec->cd()->SetLeftMargin(0.15);
//`  c_spec->cd()->SetBottomMargin(0.15);
//`  c_spec->cd()->SetTicks(1,1);
//`  c_spec->cd()->SetGrid(0,0);
//`  c_spec->SetLogy();
//`  TH1F *h_spec = new TH1F("h_spec","h_spec",100,0.0,10.0);
//`  for(int i_bin = 1; i_bin < 101; ++i_bin)
//`  {
//`    h_spec->SetBinContent(i_bin,-10.0);
//`    h_spec->SetBinError(i_bin,1.0);
//`  }
//`  h_spec->SetTitle("");
//`  h_spec->SetStats(0);
//`  h_spec->GetXaxis()->SetTitle("p_{T} (GeV/c)");
//`  h_spec->GetXaxis()->CenterTitle();
//`  h_spec->GetYaxis()->SetTitle("dN/p_{T}dp_{T}");
//`  h_spec->GetYaxis()->CenterTitle();
//`  h_spec->GetYaxis()->SetRangeUser(1E-6,10);
//`  h_spec->Draw("pE");
//`  g_spec->Draw("pE same");
//`  f_Levy->Draw("l same");
//`  f_spec->Draw("l same");
//`  */
//`
//`  return f_spec;
//`}

void write()
{
  File_OutPut->cd();
  McPhiMeson->Write();
  h_Tracks->Write();
  h_phiRP->Write();
  File_OutPut->Close();
}
