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
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "../Utility/functions.h"
#include "../Utility/StSpinAlignmentCons.h"

using namespace std;

TF1* readv2(int energy, int pid, int centrality);
TF1* readspec(int energy, int pid, int centrality);
void getKinematics(TLorentzVector& lPhi, double const mass);
void setDecayChannels(int const pid);
void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters);
void fill(TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus);
void write(int energy,int Nrho);
TVector3 CalBoostedVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec);
bool Sampling(TF1 *f_rhoPhy,float CosThetaStar);
float GausSmearing(TF1 *f_gaus);
float getChi(float Resolution);
float EventPlaneSmearing(TF1 *f_gaus);
double FuncAD(double *x_val, double *par);

// histograms
TH2F *h_cos, *h_eff, *h_ptcut;
TH2F *h_cos_narrow, *h_eff_narrow, *h_ptcut_narrow;
TProfile *resolution;
TH1F *h_theta, *h_theta_star;
TProfile *cos2phi;
TF1 *cut_pt = new TF1("cut_pt","1",0,1);

Double_t pt_set[7] = {0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};
int pt_bin;
double eta_gap;


TFile *eff_file;
TH1D *eff_hist;

// sampling functions
TF1 *f_v2, *f_spec, *f_flow, *f_rhoPhy, *f_y, *f_gaus;

TPythia6Decayer* pydecay;

void McEtaD(int Nrho = 33,int SetPt=3, int energy = 4, int pid = 0, int cent = 9, double eta_cut=1.0,  int const NMax = 1000000)
{
  int   const BinPt    = vmsa::BinPt;
  int   const BinY     = vmsa::BinY;
  int   const BinPhi   = vmsa::BinPhi;
  float const rhoDelta = vmsa::rhoDelta; //rhoDelta = 0.01

  pt_bin = SetPt;
  eta_gap = eta_cut;

//  string Info = Form("sampling rhophy = %.2f with %d tracks!!!!",rhoDelta*Nrho,NMax);
  string Info = Form("sampling rhophy = 1/3 with %d tracks!!!!",NMax);
  cout << Info.c_str() << endl;

  string HistName;

  HistName = Form("h_cos_%d",Nrho);
  h_cos = new TH2F(HistName.c_str(), HistName.c_str(), 6, 0.5, 6.5, 7, 0, 1.0);
//  resolution = new TProfile("resolution","resolution", 2, 0.5,2.5);
  cos2phi = new TProfile("cos2phi","cos2phi",2, 0.5,2.5);

  h_theta = new TH1F("h_theta","h_theta", 7, 0, 1);
  h_theta_star = new TH1F("h_theta_star","h_theta_star", 7, 0, 1);

  f_v2   = readv2(energy,pid,cent);
  f_spec = readspec(energy,pid,cent);
  f_y = new TF1("f_y","1 + 0*exp(-x*x)", -2,2);
  f_gaus = new TF1("f_gaus", "exp(-x*x/(2*[0]*[0]))/(sqrt(2.*TMath::Pi()*[0]*[0]))", -TMath::Pi(), TMath::Pi());


  f_flow = new TF1("f_flow",flowSample,-TMath::Pi(),TMath::Pi(),1);

//  float rhoPhy = Nrho*rhoDelta;
  float rhoPhy = 1./3.;
  f_rhoPhy = new TF1("f_rhoPhy",SpinDensity,-1.0,1.0,2);
  f_rhoPhy->FixParameter(0,rhoPhy);
  f_rhoPhy->FixParameter(1,0.75);

  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();
  if(gRandom) delete gRandom;
  gRandom = new TRandom3();
  gRandom->SetSeed();

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
    decayAndFill(vmsa::decayMother[pid],lPhi,ptl);
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  write(energy,Nrho);

  stopWatch->Stop();   
  stopWatch->Print();
}


TF1* readv2(int energy, int pid, int centrality)
{
  string InPutV2 = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_v2 = TFile::Open(InPutV2.c_str());
  TGraphAsymmErrors *g_v2 = (TGraphAsymmErrors*)File_v2->Get("g_v2");
  TF1 *f_v2 = new TF1("f_v2",v2_pT_FitFunc,vmsa::ptMin,vmsa::ptMax,5);
  f_v2->FixParameter(0,2);
  f_v2->SetParameter(1,0.1);
  f_v2->SetParameter(2,0.1);
  f_v2->SetParameter(3,0.1);
  f_v2->SetParameter(4,0.1);
  f_v2->SetLineColor(2);
  f_v2->SetLineWidth(2);
  f_v2->SetLineStyle(2);
  g_v2->Fit(f_v2,"N");

  /*
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
  */

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

  string InPutSpec = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_Spec.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Spec = TFile::Open(InPutSpec.c_str());
  TGraphAsymmErrors *g_spec = (TGraphAsymmErrors*)File_Spec->Get("g_spec");

  TF1 *f_Levy = new TF1("f_Levy",Levy,vmsa::ptMin,vmsa::ptMax,3);
  f_Levy->SetParameter(0,1);
  f_Levy->SetParameter(1,10);
  f_Levy->SetParameter(2,0.1);
  f_Levy->SetLineStyle(2);
  f_Levy->SetLineColor(4);
  f_Levy->SetLineWidth(2);
  g_spec->Fit(f_Levy,"N");


//  TF1 *f_spec = new TF1("f_spec",pTLevy,vmsa::ptMin,vmsa::ptMax,3);
  TF1 *f_spec = new TF1("f_spec",pTLevy,pt_set[pt_bin-1], pt_set[pt_bin], 3);
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
  c1->SaveAs("pt.eps");


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

void getKinematics(TLorentzVector& lPhi, double const mass)
{
  f_flow->ReleaseParameter(0);
//  double const pt = f_spec->GetRandom(vmsa::ptMin, vmsa::ptMax);
  double const pt = f_spec->GetRandom(pt_set[pt_bin-1], pt_set[pt_bin]);
//  double const y = gRandom->Uniform(-vmsa::acceptanceRapidity, vmsa::acceptanceRapidity);
  double const y = gRandom->Uniform(-2., 2);
//  double const y = f_y->GetRandom(-1,1);
  f_flow->SetParameter(0,f_v2->Eval(pt));
  double const phi = f_flow->GetRandom();

  // double const pt = gRandom->Uniform(vmsa::ptMin, vmsa::ptMax);
  // double const y = gRandom->Uniform(-vmsa::acceptanceRapidity, vmsa::acceptanceRapidity);
  // double const phi = TMath::TwoPi() * gRandom->Rndm();

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

  fill(lPhi,lKplus,lKminus);
}

void fill(TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus)
{
  TVector3 vMcKpBoosted = CalBoostedVector(lKplus,lPhi); // boost Kplus back to phi-meson rest frame

  double PhiPt = lPhi->Pt();
  double KplusEta = lKplus.Eta();
  double KminusEta = lKminus.Eta();

  TVector3 nQ(0.0,-1.0,0.0);
  TVector3 zQ(0.0,0.0,1.0);

  float CosThetaStarEP = vMcKpBoosted.Dot(nQ);
  float CosThetaStarZP = vMcKpBoosted.Dot(zQ);

  if(TMath::Abs(KplusEta)<=eta_gap && TMath::Abs(KminusEta)<eta_gap){
    h_theta->Fill(TMath::Abs(CosThetaStarZP));
    h_theta_star->Fill(TMath::Abs(CosThetaStarEP));
  }

  return;

}

float GausSmearing(TF1 *f_gaus)
{
  float Psi2 = f_gaus->GetRandom(-TMath::PiOver2(),TMath::PiOver2());
  return Psi2;
}

float getChi(float Resolution)
{
  TF1 *f_res = new TF1("f_res",EventPlaneResolution,0,10,0);
  double chi = f_res->GetX(Resolution);

  return chi;
}

float EventPlaneSmearing(TF1 *f_EP)
{
  float Psi2 = f_EP->GetRandom(-TMath::PiOver2(),TMath::PiOver2());
  return Psi2;
}

TVector3 CalBoostedVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec)
{
  TVector3 vMcBeta = -1.0*lMcVec->BoostVector(); // boost vector

  TLorentzVector lKaon = lMcDau;
  lKaon.Boost(vMcBeta); // boost Kplus back to phi-meson rest frame
  TVector3 vMcDauStar = lKaon.Vect().Unit(); // momentum direction of Kplus in phi-meson rest frame

  return vMcDauStar;
}

bool Sampling(TF1 *f_rhoPhy, float CosThetaStar)
{
  float wMax;
  if(f_rhoPhy->GetParameter(0) <= 1.0/3.0) wMax = f_rhoPhy->Eval(0.0);
  else wMax = f_rhoPhy->Eval(1.0);
  return !(gRandom->Rndm() > f_rhoPhy->Eval(CosThetaStar)/wMax);
}

void write(int energy, int Nrho)
{
  TF1 *Func_rho = new TF1("Func_rho","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
  TF1 *Func_A = new TF1("Func_A","[0]*(1.+[1]*(x*x))",0,1);
  TF1 *Func_AD = new TF1("Func_AD",FuncAD,0,1,4);

  double D_theta, D_theta_error;
  double D, D_error;
  double rho_before, error_before;
  double rho_after, error_after;

  Func_A->SetParameter(0,h_theta->GetBinContent(1));
  Func_A->SetParameter(1,0);
  h_theta->Fit(Func_A,"ER");
  D_theta = Func_A->GetParameter(1);
  D_theta_error = Func_A->GetParError(1);

  Func_rho->SetParameter(0,h_theta_star->GetBinContent(1));
  Func_rho->SetParameter(1,1./3.);
  h_theta_star->Fit(Func_rho,"ER");
  rho_before = Func_rho->GetParameter(1);
  error_before = Func_rho->GetParError(1);

  Func_AD->SetParameter(0,h_theta_star->GetBinContent(1));
  Func_AD->SetParameter(1,1./3.);
  Func_AD->FixParameter(2,-D_theta/(2.+D_theta));
  Func_AD->FixParameter(3,1);
  h_theta_star->Fit(Func_AD,"ER");
  D = Func_AD->GetParameter(2);
  D_error = Func_AD->GetParError(2);
  rho_after = Func_AD->GetParameter(1);
  error_after = Func_AD->GetParError(1);

  cout<<rho_before<<" +/- "<<error_before<<endl;
  cout<<rho_after<<" +/- "<<error_after<<endl;

  cout<<"pT: "<<pt_set[pt_bin-1]<<" ~ "<<pt_set[pt_bin]<<endl;
  cout<<"F':"<<D_theta<<" +/- "<<D_theta_error<<endl;
  cout<<"F: "<<D<<endl;

}

double FuncAD(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double D = par[2];
  double R = par[3];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double result = (1+Bs*D/2) + (As+D)*CosTheta*CosTheta + (As*D-Bs*D/2)*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}


