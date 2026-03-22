#include <iostream>
#include <string> 
#include <map>
#include "TFile.h"
#include "TLorentzVector.h"
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
#include "./functions.h"
#include "./StSpinAlignmentCons.h"
#include "TStyle.h"

using namespace std;
float get_Mass();
TF1* read_f_rho(double par0, double par1, double par2, double par3);
TF1* readv2(int energy, int pid, int centrality);
TF1* readspec(int energy, int pid, int centrality);
void getKinematics(TLorentzVector& lPhi, double const mass);
void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters, double const Mass_rho,int const SetPt);
void fill(TLorentzVector* lPhi, TLorentzVector const& lPiplus, TLorentzVector const& lPiminus, int const SetPt);
void write(int energy,int Nrho, double Resolution,int NMax);
TVector3 CalBoostedVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec);
TLorentzVector Cal_labVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec);
bool Sampling(TF1 *f_rhoPhy,float CosThetaStar);
float GausSmearing(TF1 *f_gaus);
float getChi(float Resolution);
float EventPlaneSmearing(TF1 *f_gaus);
double FuncAD(double *x_val, double *par);
double Func4A(double *x_val, double *par);
bool x_Sam(double S_plus);
//void gauss(double *g1, double *g2);
// histograms
TH2F *h_cos, *h_eff, *h_ptcut;
TH2F *h_cos_narrow, *h_eff_narrow, *h_ptcut_narrow;
TProfile *resolution;
TH1F *h_theta, *h_theta_star, *h_out1, *h_out2, *h_phi, *h_theta_star_before;
TH1D *h_pi_plus, *h_pi_minus,  *h_rho_phi, *h_pi_py, *h_pi_pz, *h_pi_px_unit, *h_pi_py_unit, *h_pi_pz_unit , *h_pi_pz1, *h_pi_pz2, *h_pi_pz3, *h_pi_pz4, *h_pi_pz5;
TProfile *cos2phi;
TF1 *cut_pt = new TF1("cut_pt","1",0,1);

//Double_t pt_set[29] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0};
double pt_set[5] = {1.2, 1.8, 2.4, 3.0, 4.0};
  double const m_pi = 0.49368;
TFile *eff_file;
TH1F *eff_hist;

double pt;
int pt_bin;

// sampling functions
TF1 *f_v2, *f_spec, *f_flow, *f_rhoPhy, *f_y, *f_gaus, *f_gaus_Mass, *f_rho_mass, *f_EP;


void P_phi_eta06_v2(int Nrho = 33,int SetPt=2, int energy = 2, int pid = 0, int cent = 4, int const NMax = 100000)
{
gStyle->SetOptStat(0000);
  int   const BinPt    = vmsa::BinPt;
  int   const BinY     = vmsa::BinY;
  int   const BinPhi   = vmsa::BinPhi;
  float const rhoDelta = vmsa::rhoDelta; //rhoDelta = 0.01

  pt_bin = SetPt;
  pt = (pt_set[pt_bin-1] + pt_set[pt_bin])/2.;
  string Info = Form("sampling rhophy = %.2f with %d tracks!!!!",rhoDelta*Nrho,NMax);
  cout << Info.c_str() << endl;

  string HistName;

  HistName = Form("h_cos_%d",Nrho);
  h_cos = new TH2F(HistName.c_str(), HistName.c_str(), 6, 0.5, 6.5, 7, 0, 1.0);
  resolution = new TProfile("resolution","resolution", 2, 0.5,2.5);
  cos2phi = new TProfile("cos2phi","cos2phi",2, 0.5,2.5);

  h_theta = new TH1F("h_theta","h_theta", 7, 0, 1);
  h_theta_star = new TH1F("h_theta_star","h_theta_star", 7, 0, 1);
  h_theta_star_before = new TH1F("h_theta_star_before","h_theta_star_before", 7, 0, 1);

  h_out1 = new TH1F("h_out1","h_out1", 7, 0, 1);
  h_out2 = new TH1F("h_out2","h_out2", 7, 0, 1);
  h_phi = new TH1F("h_phi","h_phi", 100, -5, 5);

  h_pi_plus = new TH1D("h_pi_plus","h_pi_plus", 15, 0, 3);
  h_pi_minus = new TH1D("h_pi_minus","h_pi_minus", 1000, 0, 3);

  h_rho_phi = new TH1D("h_rho_phi","h_rho_phi", 320, -TMath::Pi(),TMath::Pi());

  TString *Title;
  Title = new TString("phi_");
  *Title +="pt_";
  *Title +=pt_bin;
  *Title += "_hhh_pi.root";
  TFile *outputFile = new TFile(Title->Data(),"RECREATE");
  delete Title;


  TCanvas *c1 = new TCanvas();
  c1->SetFillColor(0);
  c1->SetGrid(0,0);
  c1->SetTitle(0);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.15);

  f_v2   = readv2(energy,pid,cent);
  f_v2->Draw();
  c1->SaveAs("v2.eps");
//  f_v2   = readv2(6,pid,cent);
  f_spec = readspec(energy,pid,cent);
f_rho_mass = read_f_rho(10000.,0.76850, 0.15100, pt);
  f_y = new TF1("f_y","0 + exp(-10.*x*x)", -1,1);

  TCanvas *c2 = new TCanvas();
  c2->SetFillColor(0);
  c2->SetGrid(0,0);
  c2->SetTitle(0);
  c2->SetBottomMargin(0.15);
  c2->SetLeftMargin(0.15);

  f_y->Draw();
  c2->SaveAs("f_y.eps");

  f_gaus = new TF1("f_gaus", "exp(-x*x/(2*[0]*[0]))/(sqrt(2.*TMath::Pi()*[0]*[0]))", -TMath::Pi(), TMath::Pi());
 // f_gaus_Mass = new TF1("f_gaus_Mass", "exp(-(x-0.775)*(x-0.775)/(2*0.16*0.16))/(sqrt(2.*TMath::Pi()*0.16*0.16))", 0.455, 1.095);
  f_gaus_Mass = new TF1("f_gaus_Mass", "exp(-(x-1.01946)*(x-1.01946)/(2.*0.00426*0.00426))/(sqrt(2.*TMath::Pi()*0.00426*0.00426))", 0.98, 1.05);
 // f_rho_mass = new TF1("f_rho_mass", "exp(-sqrt(x*x+0.15*0.15)/0.12)*x/sqrt(x*x+0.15*0.15)*0.775*x*(pow(((x*x - 4.*0.13957*0.13957)/(0.775*0.775 - 4.*0.13957*0.13957)),1./2.)*0.775*0.16/x)/(pow((x*x-0.775*0.775),2.) + 0.775*0.775*(pow(((x*x - 4.*0.13957*0.13957)/(0.775*0.775 - 4.*0.13957*0.13957)),1./2.)*0.775*0.16/x)*(pow(((x*x - 4.*0.13957*0.13957)/(0.775*0.775 - 4.*0.13957*0.13957)),1./2.)*0.775*0.16/x))", 0.455, 1.095);

  f_flow = new TF1("f_flow",flowSample,-TMath::Pi(),TMath::Pi(),1);

  float rhoPhy = Nrho*rhoDelta;
cout<<"rhoPhy="<<rhoPhy<<endl;
  f_rhoPhy = new TF1("f_rhoPhy",SpinDensity,-1.0,1.0,2);
  f_rhoPhy->FixParameter(0,rhoPhy);
  f_rhoPhy->FixParameter(1,0.75);

//res
  float const Resolution = 0.1;//readRes(energy,pid,cent);
  cout << "InPut Resolution = " << Resolution << endl;

  float const chi = getChi(Resolution);
  f_EP = new TF1("f_EP",EventPlaneDist,-TMath::PiOver2(),TMath::PiOver2(),2);
  f_EP->FixParameter(0,chi);
  f_EP->FixParameter(1,0.5/TMath::Pi());


  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();
  if(gRandom) delete gRandom;
  gRandom = new TRandom3();
  gRandom->SetSeed();


  TClonesArray ptl("TParticle", 10);
  TLorentzVector *lPhi = new TLorentzVector();
  for(int i_ran = 0; i_ran < NMax; ++i_ran)
  {
    if (floor(10.0*i_ran/ static_cast<float>(NMax)) > floor(10.0*(i_ran-1)/ static_cast<float>(NMax)))
    cout << "=> processing data: " << 100.0*i_ran/ static_cast<float>(NMax) << "%" << endl;
float aa_Mass = get_Mass();

    getKinematics(*lPhi,aa_Mass);
    decayAndFill(vmsa::decayMother[pid],lPhi,ptl, aa_Mass, SetPt);
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

//write output root file
  outputFile->cd();

  h_pi_plus->Write();
  h_theta_star->Write();
  h_rho_phi->Write();

  outputFile->Close();


  write(energy,Nrho,Resolution,NMax);
  stopWatch->Stop();   
  stopWatch->Print();
}

float get_Mass()
{
 // float Mass = f_gaus_Mass->GetRandom(0.455, 1.095);
 // float Mass = f_rho_mass->GetRandom(0.455, 1.095);
  float Mass = f_gaus_Mass->GetRandom(1.01094, 1.02798);
//float Mass = 0.775;
  return Mass;
/*  TF1 *f_mass = new TF1("f_mass", "[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",0.455, 1.095);
  f_mass->SetParameter(0,1);
  f_mass->SetParameter(1,0.775);
  f_mass->SetParameter(2,0.16);
  return f_mass;
*/
}

TF1* read_f_rho(double par0, double par1, double par2, double par3)
{
TF1 *f_rho_mass = new TF1("d_frho_in",d_frho_in,0.4, 1.6,4);
f_rho_mass->FixParameter(0,par0);
f_rho_mass->FixParameter(1,par1);
f_rho_mass->FixParameter(2,par2);
f_rho_mass->FixParameter(3,par3);

  return f_rho_mass;
}


TF1* readv2(int energy, int pid, int centrality)
{
  //string InPutV2 = Form("./Phi_v2_4080.root");
 // TFile *File_v2 = TFile::Open("./Rho_v2_4080.root");//InPutV2);
  TFile *File_v2 = TFile::Open("./Data/Phi_v2_1040.root");//InPutV2);
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

  //string InPutSpec = Form("./Phi_spec_4080.root");
 // TFile *File_Spec = TFile::Open("./Rho_spec_4080.root");//InPutSpec);
  TFile *File_Spec = TFile::Open("./Data/Phi_Spec.root");//InPutSpec);
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

  return f_spec;
}

void getKinematics(TLorentzVector& lPhi, double const mass)
{
  f_flow->ReleaseParameter(0);
  double const pt = f_spec->GetRandom(pt_set[pt_bin-1], pt_set[pt_bin]);
//  double const pt = gRandom->Uniform(pt_set[pt_bin-1], pt_set[pt_bin]);
  double const y = gRandom->Uniform(-1.0, 1.0);
//  double const y = f_y->GetRandom(-0.8, 0.8);
  f_flow->SetParameter(0, f_v2->Eval(pt));
//  double const phi = f_flow->GetRandom();
  double const phi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());

//  h_rho_phi->Fill(phi);
  double const mT = sqrt(mass * mass + pt * pt);
  double const pz = mT * sinh(y);
  double const E = mT * cosh(y);
//cout<<"pz= "<<pz<<"E= "<<E<<endl;
  lPhi.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}


void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters, double const Mass_rho,int const SetPt)
{
  TLorentzVector lPiplus;
  TLorentzVector lPiminus;
  double M_rho = Mass_rho;
  double pp1 = sqrt( pow((M_rho*M_rho - 2.*m_pi*m_pi),2) - pow((2.*m_pi*m_pi),2)) / (2.*M_rho);

  double const pi_phi = TMath::TwoPi() * gRandom->Rndm();
  double const pi_theta = acos(1.-2.*gRandom->Rndm());
  double pi_E = sqrt(pp1*pp1 + m_pi*m_pi);
  lPiplus.SetPxPyPzE(pp1*sin(pi_theta)*cos(pi_phi), pp1*sin(pi_theta)*sin(pi_phi), pp1*cos(pi_theta), pi_E);
  lPiminus.SetPxPyPzE(-1.*pp1*sin(pi_theta)*cos(pi_phi), -1.*pp1*sin(pi_theta)*sin(pi_phi), -1.*pp1*cos(pi_theta), pi_E);

  fill(lPhi,lPiplus,lPiminus, SetPt);
}

void fill(TLorentzVector* lPhi, TLorentzVector const& lPiplus, TLorentzVector const& lPiminus, int const SetPt)
{
//  TVector3 vMcKpBoosted = CalBoostedVector(lPiplus,lPhi); // boost Piplus back to phi-meson rest frame
  TLorentzVector lPi = lPiplus;
  TVector3 unit_Piplus = lPi.Vect().Unit();

  TLorentzVector hh_Piplus = Cal_labVector(lPiplus, lPhi);
  TLorentzVector hh_Piminus = Cal_labVector(lPiminus, lPhi);

 // double PhiPt = lPhi->Pt();
  double PiplusEta = hh_Piplus.Eta();
  double PiminusEta = hh_Piminus.Eta();
  double PiplusPt = hh_Piplus.Pt();
  double PiminusPt = hh_Piminus.Pt();

  double PiplusPy = hh_Piplus.Py();
  double PiplusPz = hh_Piplus.Pz();
//  h_pi_plus->Fill(PiplusPt);

  TVector3 nQ(0.0,-1.0,0.0); // direction of angular momentum with un-smeared EP
  TVector3 zQ(0.0,0.0,1.0);

//direction of angular momentum with EP Smearing
  float PsiEP = EventPlaneSmearing(f_EP);
  TVector3 nQSmear(sin(PsiEP),-cos(PsiEP),0.0); // direction of angular momentum with EP Smearing
  float CosThetaStarEP = unit_Piplus.Dot(nQSmear);
  float phiSmear = lPhi->Phi()-PsiEP;
  if(phiSmear > TMath::Pi())  phiSmear -= TMath::TwoPi();
  if(phiSmear < -TMath::Pi()) phiSmear += TMath::TwoPi();

  float CosThetaStarRP = unit_Piplus.Dot(nQ);
 // float CosThetaStarZP = unit_Piplus.Dot(zQ);
  if(!Sampling(f_rhoPhy,CosThetaStarRP)) return;
  h_theta_star_before->Fill(TMath::Abs(CosThetaStarEP));

  double eta_gap = 0.6;
  double pt_gap = 0.2;
  if(TMath::Abs(PiplusEta)<=eta_gap && TMath::Abs(PiminusEta)<=eta_gap) h_theta_star->Fill(TMath::Abs(CosThetaStarEP));
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

  TLorentzVector lPion = lMcDau;
  lPion.Boost(vMcBeta); // boost Piplus back to phi-meson rest frame
  TVector3 vMcDauStar = lPion.Vect().Unit(); // momentum direction of Piplus in phi-meson rest frame

  return vMcDauStar;
}

TLorentzVector Cal_labVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec)
{
  TVector3 vMcBeta = lMcVec->BoostVector(); // boost vector

  TLorentzVector lPion = lMcDau;
  lPion.Boost(vMcBeta); // boost Piplus back to lab frame


  return lPion;
}

bool Sampling(TF1 *f_rhoPhy, float CosThetaStar)
{
  float wMax;
  if(f_rhoPhy->GetParameter(0) <= 1.0/3.0) wMax = f_rhoPhy->Eval(0.0);
  else wMax = f_rhoPhy->Eval(1.0);
  return !(gRandom->Rndm() > f_rhoPhy->Eval(CosThetaStar)/wMax);
}
bool x_Sam(double S_plus)
{
  return !(gRandom->Rndm() > S_plus);
}

void write(int energy, int Nrho, double Resolution, int NMax)
{
  TF1 *Func_rho = new TF1("Func_rho","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
  TF1 *Func_A = new TF1("Func_A","[0]*(1.+[1]*(x*x))",0,1);
  TF1 *Func_AD = new TF1("Func_AD",FuncAD,0,1,4);//2-order
  double D, D_error, D_theta, D_theta_error;
//  cout<<"pT: "<<pt_set[pt_bin-1]<<" ~ "<<pt_set[pt_bin]<<endl;
  Func_rho->SetParameter(0,h_theta_star->GetBinContent(1));
  Func_rho->SetParameter(1,1./3.);
  h_theta_star->Fit(Func_rho,"ERQ");
  cout<<"picture_37, before:    "<<Func_rho->GetParameter(1)<<" ,  +/-  "<<Func_rho->GetParError(1)<<",      pT: "<<pt_set[pt_bin-1]<<" ~ "<<pt_set[pt_bin]<<endl;



  TString *Title;
      Title = new TString("./rho_pt");
        *Title += "_";
        *Title += pt_bin;
//        *Title += "_Nrho_";
//        *Title += Nrho;
        *Title += ".pdf";

  TCanvas *c12 = new TCanvas();
  c12->SetFillColor(0);
  c12->SetGrid(0,0);
  c12->SetTitle(0);
  c12->SetBottomMargin(0.15);
  c12->SetLeftMargin(0.15);
  h_theta_star_before->Draw();
  h_theta_star_before->Sumw2();
  Func_rho->SetLineColor(4);
  Func_rho->Draw("same");
  c12->SaveAs(Title->Data());
      delete Title;

//Here we start to Divide:
  TH1D *h_theta_star_clone = (TH1D*)h_theta_star->Clone("h_theta_star_clone");
  h_theta_star->Sumw2();
  h_theta_star_before->Sumw2();
  h_theta_star->Divide(h_theta_star_before);

  Func_A->SetParameter(0,h_theta_star->GetBinContent(1));
  Func_A->SetParameter(1,0);
  h_theta_star->Fit(Func_A,"ER");
  D = Func_A->GetParameter(1);
  D_error = Func_A->GetParError(1);

  cout<<"pT: "<<pt_set[pt_bin-1]<<" ~ "<<pt_set[pt_bin]<<endl;
  cout<<"D: "<<D<<" +/- "<<D_error<<endl;

  
  Func_AD->SetParameter(0, h_theta_star_clone->GetBinContent(1));
  Func_AD->SetParameter(1,1./3.);
  Func_AD->FixParameter(2, D);
//  Func_AD->FixParameter(2, -D_theta/(2.+D_theta));
  Func_AD->FixParameter(3, Resolution);
  h_theta_star_clone->Fit(Func_AD,"ERQ");
h_theta_star_clone->Draw();
  cout<<"output:before correction,   "<<Func_rho->GetParameter(1)<<" +/- "<<Func_rho->GetParError(1)<<" counts:"<<(h_theta_star->GetEntries())/NMax<<endl;
  cout<<"output:after correction,    "<<Func_AD->GetParameter(1)<<" +/- "<<Func_AD->GetParError(1)<<endl;



}
//2-order accept correction
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
/*
//4-order accept correction
double FuncAD(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double G = par[3];
  double R = par[4];

  
  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));


  double result = (2.0 + F - Bs*F/2.0 + 3.0*G/4.0 - Bs*G/2.0) 
  + (2.0*As - F*(1.0 - As - Bs) - G*(3.0/2.0 - 3.0*As/4.0 - 3.0*Bs/2.0))*CosTheta*CosTheta 
  + (G*(3.0/4.0 - 3.0*As/2.0 - 3.0*Bs/2.0) - F*(As + Bs/2.0))*CosTheta*CosTheta*CosTheta*CosTheta 
  + (G*(3.0*As/4.0 + Bs/2.0))*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta;
  return N*result;

}
*/
double Func4A(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double F = par[1];
  double G = par[2];

  double result = 1. + (4.*F+3.*G)/8. - (2.*F+3.*G)/4.*CosTheta*CosTheta + 3.*G/8.*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}
