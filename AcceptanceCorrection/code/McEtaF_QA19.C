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
#include "./Utility/functions.h"
#include "./Utility/StSpinAlignmentCons.h"

using namespace std;

TF1* readv2(int energy, int pid, int centrality);
TF1* readspec(int energy, int pid, int centrality);
void getKinematics(TLorentzVector& lPhi, double const mass);
void setDecayChannels(int const pid);
void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters);
void fill(TLorentzVector* lPhi, TLorentzVector lKplus, TLorentzVector lKminus);
void write(int energy,int Nrho);
TVector3 CalBoostedVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec);
bool Sampling(TF1 *f_rhoPhy,float CosThetaStar);
float GausSmearing(TF1 *f_gaus);
float getChi(float Resolution);
float EventPlaneSmearing(TF1 *f_gaus);
double FuncAD(double *x_val, double *par);

// histograms
TProfile *resolution;
TH1F *h_eta, *h_res, *h_both;
TH1F *h_theta, *h_theta_star;
TProfile *cos2phi;
TF1 *cut_pt = new TF1("cut_pt","1",0,1);
TProfile* CosBeta = new TProfile("CosBeta","",1,0,1);

Double_t pt_set[7] = {0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};
int pt_bin;

TFile *eff_file;
TH1D *eff_hist;

// sampling functions
TF1 *f_v2, *f_spec, *f_flow, *f_rhoPhy, *f_y, *f_gaus, *f_EP;

TPythia6Decayer* pydecay;

float const Resolution = 0.1;

void McEtaF_QA19(int Nrho = 33,int SetPt=2, int energy = 2, int pid = 0, int cent = 4, int const NMax = 1000000)
{
  int   const BinPt    = vmsa::BinPt;
  int   const BinY     = vmsa::BinY;
  int   const BinPhi   = vmsa::BinPhi;
  float const rhoDelta = vmsa::rhoDelta; //rhoDelta = 0.01

  pt_bin = SetPt;

  string Info = Form("sampling rhophy = %.2f with %d tracks!!!!",rhoDelta*Nrho,NMax);
  cout << Info.c_str() << endl;

  string HistName;

  resolution = new TProfile("resolution","resolution", 2, 0.5,2.5);
  cos2phi = new TProfile("cos2phi","cos2phi",2, 0.5,2.5);

  h_theta = new TH1F("h_theta","h_theta", 7, 0, 1);
  h_theta_star = new TH1F("h_theta_star","h_theta_star", 7, 0, 1);
  h_eta = new TH1F("h_eta","h_eta", 7, 0, 1);

  f_v2   = readv2(energy,pid,cent);
  f_spec = readspec(energy,pid,cent);
  f_y = new TF1("f_y","1 - 0*exp(-x*x)", -2,2);
  f_gaus = new TF1("f_gaus", "exp(-x*x/(2*[0]*[0]))/(sqrt(2.*TMath::Pi()*[0]*[0]))", -TMath::Pi(), TMath::Pi());


  f_flow = new TF1("f_flow",flowSample,-TMath::Pi(),TMath::Pi(),1);

  float rhoPhy = Nrho*rhoDelta;
  f_rhoPhy = new TF1("f_rhoPhy",SpinDensity,-1.0,1.0,2);
  f_rhoPhy->FixParameter(0,rhoPhy);
  f_rhoPhy->FixParameter(1,0.75);

  // float const Resolution = 0.81; // xin's resolution
  cout << "InPut Resolution = " << Resolution << endl;

  float const chi = getChi(Resolution);
  f_EP = new TF1("f_EP",EventPlaneDist,-TMath::PiOver2(),TMath::PiOver2(),2);
  f_EP->FixParameter(0,chi);
  f_EP->FixParameter(1,0.5/TMath::Pi());

  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();
  if(gRandom) delete gRandom;
  gRandom = new TRandom3();
  gRandom->SetSeed(0);

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
  string InPutV2 = Form("../res_phi/Data/Phi_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
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

  string InPutSpec = Form("../res_phi/Data/Phi_Spec.root",vmsa::mBeamEnergy[energy].c_str());
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
  double const y = gRandom->Uniform(-1., 1.);
//  double const y = f_y->GetRandom(-1,1);

//  f_flow->SetParameter(0,f_v2->Eval(pt));
//  double const phi = f_flow->GetRandom();
  double const phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());

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

void fill(TLorentzVector* lPhi, TLorentzVector lKplus, TLorentzVector lKminus)
{
  double PhiPt = lPhi->Pt();
  double PhiEta = lPhi->Eta();

  TVector3 nQ(0.0,-1.0,0.0); // direction of angular momentum with un-smeared EP
  TVector3 zQ(0.0,0.0,1.0);

  f_gaus->SetParameter(0, 0.5);
  float Psi2Gaus = GausSmearing(f_gaus);
  //resolution->Fill(1,cos(Psi2Gaus));
  //resolution->Fill(2,cos(2.*Psi2Gaus));

  TVector3 sQ(TMath::Sin(Psi2Gaus), -1.*TMath::Cos(Psi2Gaus), 0);

  float PsiEP = EventPlaneSmearing(f_EP);
  TVector3 nQSmear(sin(PsiEP),-cos(PsiEP),0.0); // direction of angular momentum with EP Smearing
  resolution->Fill(1,cos(PsiEP));
  resolution->Fill(2,cos(2.*PsiEP));

  double CosThetaStar = gRandom->Uniform(-1,1);
  double thetaStar = TMath::ACos(CosThetaStar);
  TVector3 boost = -(lKplus+lKminus).BoostVector();
  lKplus.Boost(boost); //boost to rest frame
  lKminus.Boost(boost); //boost to rest frame
  double py_cms = lKplus.P()*cos(thetaStar);
  double pxz_cms = lKplus.P()*sin(thetaStar);
  double phi = gRandom->Uniform(0,2.*TMath::Pi());
  double px_cms = pxz_cms*cos(phi);
  double pz_cms = pxz_cms*sin(phi);

  //now rotate to J direction, rotated angle around z is jrot = azimuthOfJ - pi/2
  double PI = TMath::Pi();
  double psi = nQ.Phi() + PI/2.;  
  double cosjrot = sin(psi-PI/2.);
  double sinjrot = -cos(psi-PI/2.);
  double px_cms_new = px_cms*cosjrot - py_cms*sinjrot;
  double py_cms_new = px_cms*sinjrot + py_cms*cosjrot;

  lKplus.SetPx(px_cms_new);
  lKplus.SetPy(py_cms_new);
  lKplus.SetPz(pz_cms);
  
  lKminus.SetPx(-px_cms_new);
  lKminus.SetPy(-py_cms_new);
  lKminus.SetPz(-pz_cms);

  TVector3 vMcKpBoosted = lKplus.Vect().Unit();
  
  //the beta prime in x-z, x refer to eventplane x
  psi = nQSmear.Phi() + PI/2.; 
  cosjrot = sin(psi-PI/2.);
  sinjrot = -cos(psi-PI/2.);
  px_cms_new = px_cms*cosjrot - py_cms*sinjrot;
  double beta = atan2(pz_cms,px_cms_new);

  lKplus.Boost(-boost);
  lKminus.Boost(-boost);

  double KplusEta = lKplus.Eta();
  double KminusEta = lKminus.Eta();

  double eta_gap = 1.0;

  float CosThetaStarRP = vMcKpBoosted.Dot(nQ);
  float CosThetaStarSP = vMcKpBoosted.Dot(nQSmear);
  float CosThetaStarZP = vMcKpBoosted.Dot(zQ);

  //if(TMath::Abs(KplusEta)<=eta_gap && TMath::Abs(KminusEta)<eta_gap){
  if(TMath::Abs(KplusEta)<=eta_gap && TMath::Abs(KminusEta)<=eta_gap){
    h_theta->Fill(TMath::Abs(CosThetaStarZP));
    //h_theta_star->Fill(TMath::Abs(CosThetaStarSP));
    h_theta_star->Fill(TMath::Abs(CosThetaStarRP));
    CosBeta->Fill(0.5,cos(2.*beta));
  }

  if(!Sampling(f_rhoPhy,CosThetaStarRP)) return;

  if(TMath::Abs(KplusEta)<=eta_gap && TMath::Abs(KminusEta)<=eta_gap){
    h_eta->Fill(TMath::Abs(CosThetaStarSP));
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
TVector3 CalBoostedVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec) {
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

  double res = resolution->GetBinContent(2);

  Func_A->SetParameter(0,h_theta->GetBinContent(1));
  Func_A->SetParameter(1,0);
  h_theta->Fit(Func_A,"ER");
  double F_theta = Func_A->GetParameter(1);
  double F_theta_error = Func_A->GetParError(1);
  
  Func_A->SetParameter(0,h_theta_star->GetBinContent(1));
  Func_A->SetParameter(1,0);
  h_theta_star->Fit(Func_A,"ER");
  double F = Func_A->GetParameter(1);
  double F_error = Func_A->GetParError(1);

  Func_rho->SetParameter(0,h_theta->GetBinContent(1));
  Func_rho->SetParameter(1, 1./3.);
  h_eta->Fit(Func_rho,"ER");
  double rho_obs = Func_rho->GetParameter(1);
  double rho_obs_error = Func_rho->GetParError(1);

  Func_AD->SetParameter(0,h_theta->GetBinContent(1));
  Func_AD->SetParameter(1,1./3.);
  Func_AD->FixParameter(2, F);
  Func_AD->FixParameter(3, Resolution);
  h_eta->Fit(Func_AD,"ER");
  double rho_rec = Func_AD->GetParameter(1);
  double rho_rec_error = Func_AD->GetParError(1);
  cout<<"printing parameters ::"<<endl;
  cout<<F<<" +/- "<<F_error<<endl;
  cout<<rho_obs<<" +/- "<<rho_obs_error<<endl;
  cout<<rho_rec<<" +/- "<<rho_rec_error<<endl;
  cout<<"cos2beta: "<<CosBeta->GetBinContent(1)<<" +/- "<<CosBeta->GetBinError(1)<<endl;
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


int main(int argc, char*argv[]){
int Nrho = atoi(argv[1]);
//McEtaF_QA19(Nrho,int SetPt=2, int energy = 2, int pid = 0, int cent = 4, int const NMax = 200000);
McEtaF_QA19(Nrho, 2, 2, 0, 4, 200000);
TFile* f = new TFile(Form("Hist_%s_%s.root",argv[1],argv[2]),"RECREATE");
f->cd();
h_theta->Write();
h_theta_star->Write();
h_eta->Write();
resolution->Write();
CosBeta->Write();
f->Write();
f->Close();
}
