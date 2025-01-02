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
#include "Utility/functions.h"
#include "Utility/StSpinAlignmentCons.h"

using namespace std;

TF1* readv2(int energy, int pid, int centrality);
TF1* readspec(int energy, int pid, int centrality);
void getKinematics(TLorentzVector& lPhi, double const mass);
void setDecayChannels(int const pid);
void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters, int cent);
void fill(TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus, int centrality);
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
TH1F *h_theta[6], *h_theta_before[6], *h_theta_star[6], *h_theta_star_before[6], *h_theta_star_RP[6], *h_theta_star_before_RP[6], *h_out1, *h_out2;
TH2F *h2_theta_star[6], *h2_theta_star_before[6], *h2_theta_star_RP[6], *h2_theta_star_RP_before[6];
TProfile *cos2phi;
TF1 *cut_pt = new TF1("cut_pt","1",0,1);

Double_t pt_set[7] = {0.4,0.8,1.2,1.8,2.4,3.0,4.2};
Double_t pt_centset[7] = {1.0,1.6,2.4,3.4,5.0};
Double_t pt_tuple_low[4] = {1.0,2.0,3.0,2.0};
Double_t pt_tuple_high[4] = {3.0,4.0,5.0,5.0};

int pt_bin;
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

double CentralityValues[9] = {75.0,65.0,55.0,45.0,35.0,25.0,15.0,7.5,2.5};

double yinput;

double cos2betaBefore = 0.0;
double cos4betaBefore = 0.0;
double cos2betaAfter = 0.0;
double cos4betaAfter = 0.0;

double cos2betaBeforeP = 0.0;
double cos4betaBeforeP = 0.0;
double cos2betaAfterP = 0.0;
double cos4betaAfterP = 0.0;

TH1F* h_cosbeta[6];      
TH1F* h_cos2beta[6];     
TH1F* h_cos4beta[6];     
   
TH1F* h_cosbetaP[6];     
TH1F* h_cos2betaP[6];    
TH1F* h_cos4betaP[6];    

TProfile* p_cosbeta[6];  
TProfile* p_cos2beta[6]; 
TProfile* p_cos4beta[6]; 

TProfile* p_cosbetaP[6]; 
TProfile* p_cos2betaP[6];
TProfile* p_cos4betaP[6];
   
TH1F* h_beta[6];         
TH1F* h_tstar[6];        

TH1F* h_betaP[6];        
TH1F* h_tstarP[6];       

TH2F* h_betatstar[6];    
TH2F* h_betatstarP[6];   


double cos2r2 = 0.0;
double cos2del = 0.0;
double cosdel = 0.0;
double sin2del = 0.0;
double sindel = 0.0;

double mDelta = 0.0;

TFile *eff_file;
TH1D *eff_hist;

// sampling functions
TF1 *f_v2, *f_spec, *f_flow, *f_rhoPhy, *f_y, *f_gaus;
TF1 *f_pDel, *f_res;

TPythia6Decayer* pydecay;
TNtuple* McPhiMeson;

TFile *File_OutPut;

int mEtaMode = 0;

double mChi = 0.0;

double mPsi2 = 0.0;
double mRP = 0.0;
double mPsi2RP = 0.0;

double mNEvents = 0;

bool v2flag;

double v2value;

void McEtaF_Resolution_Beta(double Nrho=3333, int SetPt=1, int energy = 4, int pid = 0, int cent = 4, int const NMax = 500, double rapidity = 1.0, int mode = 1, int etamode = 3, double Res = 0.4, bool fullEP = 0, bool v2on = 0, const char* jobID = "1")
{
  v2flag = v2on;
  mNEvents = NMax;
  mEtaMode = etamode;
  string outputfile;
  if(mode == 0) outputfile = Form("McAcceptanceOutput_pt%d_energy%d_pid%d_cent%d_EtaMode_%d_%s.root",SetPt,energy,pid,cent,mEtaMode,jobID);
  if(mode == 1) outputfile = Form("McAcceptanceOutput_pt%d_energy%d_pid%d_cent%d_EtaMode_%d_nrho%d_%s.root",SetPt,energy,pid,cent,mEtaMode,int(Nrho),jobID);
  if(mode != 0 && mode != 1) outputfile = Form("McAcceptanceOutput_pt%d_energy%d_pid%d_cent%d_%s.root",SetPt,energy,pid,cent,jobID);
  File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  cout << "Output set to" << outputfile << endl;

  int   const BinPt    = vmsa::BinPt;
  int   const BinY     = vmsa::BinY;
  int   const BinPhi   = vmsa::BinPhi;
  float const rhoDelta = 0.0001; //rhoDelta = 0.0001


  int BufSize = (int)pow(2., 16.);
  // int Split = 1;

  const char* varlist = "Centrality:McPt:McP:McEta:McY:McPhi:McCos:McCosTheta:McKpEta:McKmEta:McCosRP"; // reconstructed phi

  McPhiMeson = new TNtuple("McPhiMeson", "McPhiMeson", varlist, BufSize);


  mMode = mode;

  pt_bin = SetPt; 
  //if(mode == 3) 
  //{
  //  pt_bin = 1;
  //  pt_set[1] = 5.0;
  //  pt_set[0] = 1.0;
  //}

  yinput = rapidity;

  //TString InPutFile_Res = Form("Utility/Resolution/file_%s_Resolution.root",vmsa::mBeamEnergy[energy].c_str());
  //mInPutFile_Res = TFile::Open(InPutFile_Res.Data());

  //TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get("p_mRes2_Sub");
  //Float_t Res_raw = p_res2->GetBinContent(p_res2->FindBin(cent));
  //Res = TMath::Sqrt(Res_raw);
  cout << "Resolution = " << Res << endl;
  f_res = new TF1("resolution",EventPlaneResolution,0,80,0);
  mChi = f_res->GetX(Res);
  //if(fullEP) mChi = mChi * TMath::Sqrt(2);
  //for( int i = 1; i <= 100; i++){
  //  cout << "res = " <<  double(i)/100. << ",   chi = " << f_res->GetX(double(i)/100.) << endl;
  //}
  cout << "mChi = " << mChi << endl;
  f_pDel = new TF1("deltaPsi",EventPlaneDist,-TMath::Pi()/2.0,TMath::Pi()/2.0,2);
  f_pDel->FixParameter(0,mChi);
  f_pDel->FixParameter(1,1.0/(2.0*TMath::Pi()));

  
  TCanvas *c_del = new TCanvas("c_del","c_del",10,10,800,800);
  c_del->cd()->SetLeftMargin(0.15);
  c_del->cd()->SetBottomMargin(0.15);
  c_del->cd()->SetTicks(1,1);
  c_del->cd()->SetGrid(0,0);
  TH1F *h_del = new TH1F("h_del","h_del",100,-10.0,10.0);
  for(int i_bin = 1; i_bin < 101; ++i_bin)
  {
    h_del->SetBinContent(i_bin,-10.0);
    h_del->SetBinError(i_bin,1.0);
  }
  h_del->GetXaxis()->SetRangeUser(-TMath::Pi()/2.0,TMath::Pi()/2.0);
  h_del->SetTitle("");
  h_del->SetStats(0);
  h_del->GetXaxis()->SetTitle("#Delta");
  h_del->GetXaxis()->CenterTitle();
  h_del->GetYaxis()->SetTitle("arb. units");
  h_del->GetYaxis()->CenterTitle();
  h_del->GetYaxis()->SetRangeUser(0.0,1.0);
  h_del->Draw("pE");
  //g_del->Draw("pE same");
  f_pDel->Draw("l same");
  c_del->SaveAs("del.pdf");

  cout << "ptbin = " <<  pt_bin << "     yinput = " << yinput << endl;
  if(mMode == 0) cout << "pt range = [" << pt_set[pt_bin-1] << "," << pt_set[pt_bin] << "] GeV/c" << endl;  
  if(mMode == 1) cout << "pt range = [" << 1.2 << "," << 1.8 << "] GeV/c" << endl;  
  if(mMode == 3) cout << "pt range = [" << pt_tuple_low[pt_bin-1] << "," << pt_tuple_high[pt_bin-1] << "] GeV/c" << endl;  

  string Info = Form("sampling rhophy = %.2f with %d tracks!!!!",rhoDelta*Nrho,NMax);
  cout << Info.c_str() << endl;

  string HistName;

  HistName = Form("h_cos_%d",Nrho);
  h_cos = new TH2F(HistName.c_str(), HistName.c_str(), 6, 0.5, 6.5, 7, 0, 1.0);
  resolution = new TProfile("resolution","resolution", 2, 0.5,2.5);
  cos2phi = new TProfile("cos2phi","cos2phi",2, 0.5,2.5);


  std::string histname;
  for(int ie; ie < 6; ie++) //index corresponds to eta cut, 0 = none, 1 = 1.0, 2 = 0.8, 3 = 0.6, 4 = 0.4, 5 = 0.2
  {
    histname = Form("h_theta_%d",ie);
    h_theta[ie]           = new TH1F(histname.c_str(),histname.c_str(), 10, 0, 1);
    histname = Form("h_theta_star_%d",ie);
    h_theta_star[ie]      = new TH1F(histname.c_str(),histname.c_str(), 7, 0, 1);
    histname = Form("h_theta_star_RP_%d",ie);
    h_theta_star_RP[ie]   = new TH1F(histname.c_str(),histname.c_str(), 7, 0, 1);

    histname = Form("h2_theta_star_%d",ie);
    h2_theta_star[ie]     = new TH2F(histname.c_str(),histname.c_str(), 7, 0, 1, 30, 0.0, 2*TMath::Pi());
    histname = Form("h2_theta_star_RP_%d",ie);
    h2_theta_star_RP[ie]  = new TH2F(histname.c_str(),histname.c_str(), 7, 0, 1, 30, 0.0, 2*TMath::Pi());


    histname = Form("h_cosbeta_%d",ie);
    h_cosbeta[ie]      = new TH1F(histname.c_str(),histname.c_str(), 50, -1, 1);
    histname = Form("h_cos2beta_%d",ie);
    h_cos2beta[ie]     = new TH1F(histname.c_str(),histname.c_str(), 50, -1, 1);
    histname = Form("h_cos4beta_%d",ie);
    h_cos4beta[ie]     = new TH1F(histname.c_str(),histname.c_str(), 50, -1, 1); 

    histname = Form("h_cosbetaP_%d",ie);
    h_cosbetaP[ie]     = new TH1F(histname.c_str(),histname.c_str(), 50, -1, 1);
    histname = Form("h_cos2betaP_%d",ie);
    h_cos2betaP[ie]    = new TH1F(histname.c_str(),histname.c_str(), 50, -1, 1);
    histname = Form("h_cos4betaP_%d",ie);
    h_cos4betaP[ie]    = new TH1F(histname.c_str(),histname.c_str(), 50, -1, 1);
    
    histname = Form("p_cosbeta_%d",ie);
    p_cosbeta[ie]      = new TProfile(histname.c_str(),histname.c_str(), 7, 0.0, 1.0);
    histname = Form("p_cos2beta_%d",ie);
    p_cos2beta[ie]     = new TProfile(histname.c_str(),histname.c_str(), 7, 0.0, 1.0);
    histname = Form("p_cos4beta_%d",ie);
    p_cos4beta[ie]     = new TProfile(histname.c_str(),histname.c_str(), 7, 0.0, 1.0);
    
    histname = Form("p_cosbetaP_%d",ie);
    p_cosbetaP[ie]     = new TProfile(histname.c_str(),histname.c_str(), 7, 0.0, 1.0);
    histname = Form("p_cos2betaP_%d",ie);
    p_cos2betaP[ie]    = new TProfile(histname.c_str(),histname.c_str(), 7, 0.0, 1.0);
    histname = Form("p_cos4betaP_%d",ie);
    p_cos4betaP[ie]    = new TProfile(histname.c_str(),histname.c_str(), 7, 0.0, 1.0);
       
    histname = Form("h_beta_%d",ie);
    h_beta[ie]         = new TH1F(histname.c_str(),histname.c_str(), 50, 0.0, 2.0*TMath::Pi());
    histname = Form("h_tstar_%d",ie);
    h_tstar[ie]        = new TH1F(histname.c_str(),histname.c_str(), 50, 0.0, TMath::Pi());
    
    histname = Form("h_betaP_%d",ie);
    h_betaP[ie]        = new TH1F(histname.c_str(),histname.c_str(), 50, 0.0, 2.0*TMath::Pi());
    histname = Form("h_tstarP_%d",ie);
    h_tstarP[ie]       = new TH1F(histname.c_str(),histname.c_str(), 50, 0.0, TMath::Pi());
    
    histname = Form("h_betatstar_%d",ie);
    h_betatstar[ie]    = new TH2F(histname.c_str(),histname.c_str(), 50, 0.0, TMath::Pi(), 50, 0.0, 2.0*TMath::Pi());
    histname = Form("h_betatstarP_%d",ie);
    h_betatstarP[ie]   = new TH2F(histname.c_str(),histname.c_str(), 50, 0.0, TMath::Pi(), 50, 0.0, 2.0*TMath::Pi());
  }


  f_v2   = readv2(energy,pid,cent);
  v2value = f_v2->Eval(1.5); 
  cout << "v2 value = " << v2value << endl;
 // TCanvas *c2 = new TCanvas("c2","c2");
 // c2->SetFillColor(0);
 // c2->SetGrid(0,0);
 // c2->SetTitle(0);
 // c2->SetBottomMargin(0.15);
 // c2->SetLeftMargin(0.15);
 // f_v2->Draw();
 // c1->SaveAs("v2.pdf");
  //f_v2   = readv2(6,pid,cent);
  f_spec = readspec(energy,pid,cent);
  //f_y = new TF1("f_y","0 + exp(-100*x*x)", -6.0,6.0);
  //f_gaus = new TF1("f_gaus", "exp(-x*x/(2*[0]*[0]))/(sqrt(2.*TMath::Pi()*[0]*[0]))", -TMath::Pi(), TMath::Pi());

  //f_gaus = new TF1("f_gaus", EventPlaneGaus, -TMath::Pi()/2., TMath::Pi(), 1);
  //f_gaus->FixParameter(0, EP_res);

  f_flow = new TF1("f_flow",flowSample,-TMath::Pi(),TMath::Pi(),1);

  float rhoPhy = Nrho*rhoDelta;
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
    decayAndFill(vmsa::decayMother[pid],lPhi,ptl,cent);

    if (mMode == 3 && i_ran % 1000 == 1) McPhiMeson->AutoSave("SaveSelf");
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  write(energy,Nrho);

  stopWatch->Stop();   
  stopWatch->Print();
}


TF1* readv2(int energy, int pid, int centrality)
{
  //string centFile[9] = {"4080","4080","4080","4080","1040","1040","1040","0010","0010"};
  string InPutV2 = Form("/gpfs01/star/pwg/sunxuhit/BackupData/SpinAlignment/data/AuAu%s/Phi/MonteCarlo/Data/Phi_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
  //if((energy == 4 || energy == 3 || energy == 2 || energy == 0) && (mMode == 1 || mMode == 3)) InPutV2 = "/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/v2/HEPData-ins1395151-v2-root.root";
  TFile *File_v2 = TFile::Open(InPutV2.c_str());
  std::cout << "v2 file: " << InPutV2 << endl;

  TGraphAsymmErrors *g_v2;

  if(mMode == 0) g_v2 = (TGraphAsymmErrors*)File_v2->Get("g_v2");
  if(mMode == 1) g_v2 = (TGraphAsymmErrors*)File_v2->Get("g_v2");

  //if((energy == 4 || energy == 3 || energy == 2 || energy == 0) && (mMode == 1 || mMode == 3)) 
  //{
  //  TDirectory *dir = (TDirectory*) File_v2->Get(Form("Table %d",tableNumV2[energy][centrality]));
  //  dir->cd(); 
  //  g_v2 = (TGraphAsymmErrors*)dir->Get("Graph1D_y1");
  //}

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
  h_v2->GetYaxis()->SetRangeUser(0.0,0.2);
  h_v2->Draw("pE");
  g_v2->Draw("pE same");
  f_v2->Draw("l same");
  c_v2->SaveAs("v2.pdf");



  return f_v2;
}

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

  string InPutSpec = Form("/gpfs01/star/pwg/sunxuhit/BackupData/SpinAlignment/data/AuAu%s/Phi/MonteCarlo/Data/Phi_Spec.root",vmsa::mBeamEnergy[energy].c_str());
  //if((energy == 4 || energy == 2 || energy == 0) && (mMode == 1 || mMode == 3)) InPutSpec = "/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/pTspectra/HEPData-ins1378002-v1-root.root";
  TFile *File_Spec = TFile::Open(InPutSpec.c_str());
  cout << "Input spectra" << InPutSpec << endl;
 
  TGraphAsymmErrors *g_spec;
  if(mMode == 0) g_spec = (TGraphAsymmErrors*)File_Spec->Get("g_spec");
  if(mMode == 1) g_spec = (TGraphAsymmErrors*)File_Spec->Get("g_spec");
  //if((energy == 4 || energy == 2 || energy == 0) && (mMode == 1 || mMode == 3)) 
  //{
  //  TDirectory *dir = (TDirectory*) File_Spec->Get(Form("Table %d",tableNum[energy][centrality]));
  //  dir->cd(); 
  //  g_spec = (TGraphAsymmErrors*)dir->Get("Graph1D_y1");
  //  for(int i = 0; i < g_spec->GetN(); i++)
  //  {
  //    double x,y;
  //    g_spec->GetPoint(i,x,y);
  //    g_spec->SetPoint(i,ptvalues[energy][centrality][i],y);
  //    g_spec->SetPointError(i,fabs(ptvalues[energy][centrality][i]-ptlow[energy][i]),fabs(ptvalues[energy][centrality][i]-pthigh[energy][i]),g_spec->GetErrorYlow(i),g_spec->GetErrorYhigh(i));
  //  }
  //}

  TF1 *f_Levy = new TF1("f_Levy",Levy,vmsa::ptMin,vmsa::ptMax,3);
  f_Levy->SetParameter(0,1);
  f_Levy->SetParameter(1,10);
  f_Levy->SetParameter(2,0.1);
  f_Levy->SetLineStyle(2);
  f_Levy->SetLineColor(4);
  f_Levy->SetLineWidth(2);
  cout << "Fitting full pT distribution" << endl;
  g_spec->Fit(f_Levy,"N");


//  TF1 *f_spec = new TF1("f_spec",pTLevy,vmsa::ptMin,vmsa::ptMax,3);
  TF1 *f_spec;
  if(mMode == 0) f_spec = new TF1("f_spec",pTLevy,pt_set[pt_bin-1], pt_set[pt_bin], 3);
  //if(mMode == 1) f_spec = new TF1("f_spec",pTLevy,pt_tuple_low[pt_bin-1], pt_tuple_high[pt_bin-1], 3);
  if(mMode == 1) f_spec = new TF1("f_spec",pTLevy,1.2, 1.8, 3);
  //if(mMode == 3) f_spec = new TF1("f_spec",pTLevy, 1.0, 5.0, 3);
  if(mMode == 3) f_spec = new TF1("f_spec",pTLevy,pt_tuple_low[pt_bin-1], pt_tuple_high[pt_bin-1], 3);
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
  c1->SaveAs("pt.pdf");

  TCanvas *c10 = new TCanvas("c_ptdist","c_ptdist",10,10,800,800);
  c10->cd();
  c10->SetFillColor(0);
  c10->SetGrid(0,0);
  c10->SetTitle(0);
  c10->SetBottomMargin(0.15);
  c10->SetLeftMargin(0.15);
  f_spec->Draw();
  f_spec->SaveAs("indivitualptspec.pdf");
  

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
  double pt; 
  if(mMode == 0) pt = f_spec->GetRandom(pt_set[pt_bin-1], pt_set[pt_bin]);
  //if(mMode == 1) pt = f_spec->GetRandom(1.2, 1.8);
  if(mMode == 1) pt = 1.5;
  //if(mMode == 3) pt = f_spec->GetRandom(1.0, 5.0);
  if(mMode == 3) pt = f_spec->GetRandom(pt_tuple_low[pt_bin-1],pt_tuple_high[pt_bin-1]);
  //double const pt = gRandom->Uniform(pt_set[pt_bin-1], pt_set[pt_bin]);
//  double const y = gRandom->Uniform(-vmsa::acceptanceRapidity, vmsa::acceptanceRapidity);
  //double const y = gRandom->Uniform(-1., 1.);
  double y, Psi2 = 0.0;
  mPsi2 = 0.0;
  mPsi2RP = 0.0;
  if(mMode == 0 || mMode == 3) y = gRandom->Uniform(-yinput, yinput);
  if(mMode == 1)
  {
    if(mEtaMode == 0) y = gRandom->Uniform(0.8,1.0);
    if(mEtaMode == 3) y = gRandom->Uniform(0.2,0.4);
    if(mEtaMode == 4) y = gRandom->Uniform(0.4,0.6);
    if(mEtaMode == 5) y = gRandom->Uniform(0.6,0.8);
  }    
  //if(mMode == 1)
  //{
  //  if(gRandom->Uniform() > 0.5) y = gRandom->Uniform(0.8,1.0);
  //  else y = gRandom->Uniform(-1.0,-0.8);
  //}
  //cout << " y = " << y << endl;

  //Let's assume for the time being that the RP angle is always 0;
  Psi2 = 0.0; 

  double delta = f_pDel->GetRandom();
  mDelta = delta;

  mPsi2 = Psi2 + delta; // EP, the RP is smeared by delta
  mPsi2RP = Psi2;       // RP

  //Angle wrapping
  while(mPsi2 > 0.5*TMath::Pi())
  {
    mPsi2 -= TMath::Pi();
  }
  while(mPsi2 < -0.5*TMath::Pi())
  {
    mPsi2 += TMath::Pi();
  }

  double phi = 0.0;  
  if(v2flag)
  {
    f_flow->SetParameter(0,v2value);//f_v2->Eval(pt));
    //f_flow->SetParameter(0,f_v2->Eval(pt));
    //f_flow->SetParameter(0,1.5*f_v2->Eval(pt)); //scale up by 50%
    phi = f_flow->GetRandom(); // This is actually phi - RP
    phi += mPsi2RP; //Therefore we must add RP
  }  
  else 
  {
    phi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
  }


  //phi = phi + Psi2;// + delta;
  //phi = phi + mPsi2;// + delta;

  cos2r2  += TMath::Cos(2.0*(mPsi2-mPsi2RP));
  cos2del += TMath::Cos(2.0*delta);
  cosdel  += TMath::Cos(delta);
  sin2del += TMath::Sin(2.0*delta);
  sindel  += TMath::Sin(delta);

  //cout << "mRP = " << mRP << ",     mPsi2 = " << mPsi2 << ",    delta = " << delta << endl;

  //mPsi2RP = Psi2;
  //cout << "Original Psi2 EP = " << mPsi2 << endl;


  //cout << "Wrapped Psi2 EP = " << mPsi2 << endl;

  //cout << "y = " << y << ", Psi2 = " << Psi2 << ", delta = " << delta << ", phi = " << defaultphi << ", final phi = " << phi << endl;

  //cout << "pt = " << pt << "     y = " << y << "     phi = " << phi << endl;

  //double const phi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
  // double const pt = gRandom->Uniform(vmsa::ptMin, vmsa::ptMax);
  // double const y = gRandom->Uniform(-vmsa::acceptanceRapidity, vmsa::acceptanceRapidity);
  //double const phi = TMath::TwoPi() * gRandom->Rndm();

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

void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters, int cent)
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

  fill(lPhi,lKplus,lKminus,cent);
}

void fill(TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus, int centrality)
{
  TVector3 vMcKpBoosted = CalBoostedVector(lKplus,lPhi); // boost Kplus back to phi-meson rest frame

  double PhiRapidity = lPhi->Rapidity();
  double PhiPt = lPhi->Pt();
  double KplusPhi = vMcKpBoosted.Phi();
  double KplusEta = lKplus.Eta();
  double KminusEta = lKminus.Eta();
  //double KplusPt = lKplus.Pt();
  //double KminusPt = lKminus.Pt();

  //cout << "K+ eta = " << KplusEta << "        K+ boosted eta = " << vMcKpBoosted.Eta() << endl;

  TVector3 nQ(TMath::Sin(mPsi2),-1.0*TMath::Cos(mPsi2),0.0); // direction of angular momentum with un-smeared EP
  TVector3 nQRP(TMath::Sin(mPsi2RP),-1.0*TMath::Cos(mPsi2RP),0.0); // direction of angular momentum with un-smeared EP
  //TVector3 nQRP(TMath::Sin(mRP),-1.0*TMath::Cos(mRP),0.0); // direction of angular momentum with un-smeared EP
  TVector3 zQ(0.0,0.0,1.0);

  float Theta = 2.*TMath::ATan(TMath::Exp(-KplusEta));
  float CosThetaInLab = TMath::Cos(Theta);

  TVector3 nPP = (lPhi->Vect().Unit()).Cross(zQ).Unit();

  TVector3 vMcBeta = -1.0*lPhi->BoostVector(); // boost vector
  
  TLorentzVector lPsi2RP, lPsi2EP, lZ;

  lPsi2RP.SetXYZT(nQRP.X(), nQRP.Y(), nQRP.Z(), 0.0);
  lPsi2EP.SetXYZT(nQ.X(), nQ.Y(), nQ.Z(), 0.0);
  lZ.SetXYZT(zQ.X(), zQ.Y(), zQ.Z(), 0.0);

  lPsi2RP.Boost(vMcBeta); // boost normal to RP to phi-meson rest frame
  lPsi2EP.Boost(vMcBeta); // boost normal to EP to phi-meson rest frame
  lZ.Boost(vMcBeta); // boost normal to EP to phi-meson rest frame

  TVector3 vNRP = lPsi2RP.Vect().Unit(); // momentum direction of Kplus in phi-meson rest frame
  TVector3 vNEP = lPsi2EP.Vect().Unit(); // momentum direction of Kplus in phi-meson rest frame
  TVector3 vZ = lZ.Vect().Unit(); // momentum direction of Kplus in phi-meson rest frame

  float CosThetaStarEP = vMcKpBoosted.Dot(nQ);
  float CosThetaStarRP = vMcKpBoosted.Dot(nQRP);
  float CosThetaStarZP = vMcKpBoosted.Dot(zQ);
  //float CosThetaStarEP = vMcKpBoosted.Dot(vNEP);
  //float CosThetaStarRP = vMcKpBoosted.Dot(vNRP);
  //float CosThetaStarZP = vMcKpBoosted.Dot(vZ);
  
  //cout << "cos(theta) in lab = " << CosThetaInLab << ",       cos(theta) in phi meson rest frame = " << CosThetaStarZP << endl;
  //cout << "Theta in lab = " << Theta << ",        Theta in phi meson rest frame = " << TMath::ACos(CosThetaStarZP) << endl;

  if(!Sampling(f_rhoPhy,TMath::Abs(CosThetaStarRP))) return;

  float SinThetaStarEP = TMath::Sin(TMath::ACos(CosThetaStarEP));
  float SinThetaStarRP = TMath::Sin(TMath::ACos(CosThetaStarRP));
  float SinThetaStarZP = TMath::Sin(TMath::ACos(CosThetaStarZP));

  // Need to calculate sin(beta)
  TVector3 xVecRP(TMath::Cos(mPsi2RP),TMath::Sin(mPsi2RP),0.0); // direction of x axis in RP frame
  //TVector3 xVecRP(TMath::Cos(mRP),TMath::Sin(mRP),0.0); // direction of x axis in RP frame
  TVector3 xVecEP(TMath::Cos(mPsi2),TMath::Sin(mPsi2),0.0); // direction of x axis in EP frame

  TLorentzVector xPsi2RP, xPsi2EP;

  xPsi2RP.SetXYZT(xVecRP.X(), xVecRP.Y(), xVecRP.Z(), 0.0);
  xPsi2EP.SetXYZT(xVecEP.X(), xVecEP.Y(), xVecEP.Z(), 0.0);

  xPsi2RP.Boost(vMcBeta); // boost normal to RP to phi-meson rest frame
  xPsi2EP.Boost(vMcBeta); // boost normal to EP to phi-meson rest frame

  TVector3 xNRP = xPsi2RP.Vect().Unit(); // momentum direction of Kplus in phi-meson rest frame
  TVector3 xNEP = xPsi2EP.Vect().Unit(); // momentum direction of Kplus in phi-meson rest frame

  float SinBeta  = -(xVecRP.Dot(nQRP.Cross(vMcKpBoosted)))/(nQRP.Cross(vMcKpBoosted)).Mag();
  float SinBetaP = -(xVecEP.Dot(nQ.Cross(vMcKpBoosted)))/(nQ.Cross(vMcKpBoosted)).Mag();
  float Beta  = TMath::ASin(SinBeta);
  float BetaP = TMath::ASin(SinBetaP);

  //float SinBeta  = -(xVecRP.Dot(nQRP.Cross(vMcKpBoosted)))/(nQRP.Cross(vMcKpBoosted)).Mag();
  //float SinBetaP = -(xVecEP.Dot(nQ.Cross(vMcKpBoosted)))/(nQ.Cross(vMcKpBoosted)).Mag();
  //float Beta  = TMath::ASin(SinBeta);
  //float BetaP = TMath::ASin(SinBetaP);

  //float SinBeta  = -(xNRP.Dot(vNRP.Cross(vMcKpBoosted)))/(vNRP.Cross(vMcKpBoosted)).Mag();
  //float SinBetaP = -(xNEP.Dot(vNEP.Cross(vMcKpBoosted)))/(vNEP.Cross(vMcKpBoosted)).Mag();
  //float Beta  = TMath::ASin(SinBeta);
  //float BetaP = TMath::ASin(SinBetaP);


  //cout << endl << endl << "Beta = " << Beta << "       Beta' = " << BetaP << endl;
  //cout << "K+ direction in phi meson frame" << endl;
  //vMcKpBoosted.Print();

  TVector3 vMcKpBoostedInRP = vMcKpBoosted;
           vMcKpBoostedInRP.RotateZ(-mPsi2RP);
           //vMcKpBoostedInRP.RotateZ(-mRP);
  TVector3 vMcKpBoostedInEP = vMcKpBoosted;
           vMcKpBoostedInEP.RotateZ(-mPsi2);

  //cout << "Psi2 RP = " << mPsi2RP << endl;
  //cout << "K+ direction in phi meson RP frame" << endl;
  //vMcKpBoostedInRP.Print();
  //cout << "Psi2 EP = " << mPsi2 << endl;
  //cout << "K+ direction in phi meson EP frame" << endl;
  //vMcKpBoostedInEP.Print();

  //if( vMcKpBoostedInRP.GetX() >= 0.0 && vMcKpBoostedInRP.GetZ() >= 0.0) 
  double c2betaBW = TMath::Cos(2.0*Beta);
  double c2betapBW = TMath::Cos(2.0*BetaP);

  if( vMcKpBoostedInRP.X() >= 0.0 && vMcKpBoostedInRP.Z() <  0.0) Beta = Beta + 2*TMath::Pi();
  if( vMcKpBoostedInRP.X() <  0.0 && vMcKpBoostedInRP.Z() >= 0.0) Beta = TMath::Pi() - Beta;
  if( vMcKpBoostedInRP.X() <  0.0 && vMcKpBoostedInRP.Z() <  0.0) Beta = TMath::Pi() - Beta;

  if( vMcKpBoostedInEP.X() >= 0.0 && vMcKpBoostedInEP.Z() <  0.0) BetaP = BetaP + 2*TMath::Pi();
  if( vMcKpBoostedInEP.X() <  0.0 && vMcKpBoostedInEP.Z() >= 0.0) BetaP = TMath::Pi() - BetaP;
  if( vMcKpBoostedInEP.X() <  0.0 && vMcKpBoostedInEP.Z() <  0.0) BetaP = TMath::Pi() - BetaP;

  //cout << "cos(2beta)  difference =  " << c2betaBW-TMath::Cos(2.0*Beta) << endl;
  //cout << "cos(2beta') difference =  " << c2betapBW-TMath::Cos(2.0*BetaP) << endl;

  h_cosbeta[0]->Fill(TMath::Cos(Beta));
  h_cosbetaP[0]->Fill(TMath::Cos(BetaP));
  h_cos2beta[0]->Fill(TMath::Cos(2.0*Beta));
  h_cos2betaP[0]->Fill(TMath::Cos(2.0*BetaP));
  h_cos4beta[0]->Fill(TMath::Cos(4.0*Beta));
  h_cos4betaP[0]->Fill(TMath::Cos(4.0*BetaP));

  p_cosbeta[0]->Fill(TMath::Abs(CosThetaStarRP),TMath::Cos(Beta));
  p_cosbetaP[0]->Fill(TMath::Abs(CosThetaStarEP),TMath::Cos(BetaP));
  p_cos2beta[0]->Fill(TMath::Abs(CosThetaStarRP),TMath::Cos(2.0*Beta));
  p_cos2betaP[0]->Fill(TMath::Abs(CosThetaStarEP),TMath::Cos(2.0*BetaP));
  p_cos4beta[0]->Fill(TMath::Abs(CosThetaStarRP),TMath::Cos(4.0*Beta));
  p_cos4betaP[0]->Fill(TMath::Abs(CosThetaStarEP),TMath::Cos(4.0*BetaP));

  h_beta[0]->Fill(Beta);    
  h_tstar[0]->Fill(TMath::ACos(CosThetaStarRP));   
                  
  h_betaP[0]->Fill(BetaP);   
  h_tstarP[0]->Fill(TMath::ACos(CosThetaStarEP));  
  
  h_betatstar[0]->Fill(TMath::ACos(CosThetaStarRP),Beta);                           
  h_betatstarP[0]->Fill(TMath::ACos(CosThetaStarEP),BetaP);   

  //cout << "True Beta = " << Beta << "       True Beta' = " << BetaP << endl;

  //TMatrixD RotIntoRPFrame(3,3);
  //RotIntoRPFrame(0,0) = TMath::Cos(mPsi2RP);          RotIntoRPFrame(0,1) = 0.0;              RotIntoRPFrame(0,2) = -TMath::Sin(mPsi2RP);
  //RotIntoRPFrame(1,0) = 0.0;                          RotIntoRPFrame(1,1) = 1.0;              RotIntoRPFrame(1,2) = 0.0;
  //RotIntoRPFrame(2,0) = TMath::Sin(mPsi2RP);          RotIntoRPFrame(2,1) = 0.0;              RotIntoRPFrame(2,2) = TMath::Cos(mPsi2RP);
  //  
  //TMatrixD RotIntoEPFrame(3,3);
  //RotIntoEPFrame(0,0) = TMath::Cos(mPsi2);          RotIntoEPFrame(0,1) = 0.0;              RotIntoEPFrame(0,2) = -TMath::Sin(mPsi2);
  //RotIntoEPFrame(1,0) = 0.0;                        RotIntoEPFrame(1,1) = 1.0;              RotIntoEPFrame(1,2) = 0.0;
  //RotIntoEPFrame(2,0) = TMath::Sin(mPsi2);          RotIntoEPFrame(2,1) = 0.0;              RotIntoEPFrame(2,2) = TMath::Cos(mPsi2);


  

  //float SinBeta = -(zQ.Dot(nQRP.Cross(vMcKpBoosted)))/(nQRP.Cross(vMcKpBoosted)).Mag();
  //float SinBetaP = -(zQ.Dot(nQ.Cross(vMcKpBoosted)))/(nQ.Cross(vMcKpBoosted)).Mag();

  //cout << "Mag of vMcKpBoosted = " << vMcKpBoosted.Mag() << "                 Mag of cross = " << (zQ.Cross(vMcKpBoosted)).Mag() << endl;
  //cout << "vMcKpBoosted " << endl;
  //vMcKpBoosted.Print();
  //cout << "Cross Product of zQ and vMcKpBoosted " << endl;
  //(zQ.Cross(vMcKpBoosted)).Print();

  //cout << "CosThetaStarRP (from dot product) = " << CosThetaStarRP << "   cos(theta*) from sin(theta)sin(Psi-phi) = " << SinThetaStarZP*TMath::Sin(mRP - KplusPhi) << "   cos(theta*) from sin(theta)sin(phi-Psi2) = " << SinThetaStarZP*TMath::Sin(KplusPhi-mRP) << endl;
  //cout << "CosThetaStarEP (from dot product) = " << CosThetaStarEP << "   cos(theta*') from sin(theta)sin(Psi2-phi) = " << SinThetaStarZP*TMath::Sin(mPsi2 - KplusPhi) << "   cos(theta*') from sin(theta)sin(phi-Psi2) = " << SinThetaStarZP*TMath::Sin(KplusPhi-mPsi2) << endl;
  //cout << "CosThetaStarRP (from dot product) = " << CosThetaStarRP << "   cos(theta*) from sin(theta)sin(Psi2' - delta -phi) = " << SinThetaStarZP*TMath::Sin(mPsi2 -mDelta - KplusPhi) << "   cos(theta*) from sin(theta)sin(phi-Psi2'+delta) = " << SinThetaStarZP*TMath::Sin(KplusPhi-mPsi2+mDelta) << endl;
  
   
  //cout << "CosThetaZP (from dot product) = " << CosThetaStarZP << "   cos(theta) from sin(theta*)sin(beta) = " << SinThetaStarRP*TMath::Sin(Beta) << endl;
  //cout << "CosThetaZP (from dot product) = " << CosThetaStarZP << "   cos(theta) from sin(theta*')sin(beta') = " << SinThetaStarEP*TMath::Sin(BetaP) << endl;
  



  double eta_gap[5] = {1.0,0.8,0.6,0.4,0.2};
  //double pt_gap = 0.2;
  //double pt_plus=lKplus.Pt(), pt_minus=lKminus.Pt();
  //double ratio = TMath::Abs(pt_plus-pt_minus);
  float arr[12];
  int iArr = 0;
  if(mMode == 3)
  {
    arr[iArr++] = centrality; // McPhi
    arr[iArr++] = lPhi->Pt();
    arr[iArr++] = lPhi->P();
    arr[iArr++] = lPhi->PseudoRapidity();
    arr[iArr++] = lPhi->Rapidity();
    arr[iArr++] = lPhi->Phi();
    arr[iArr++] = CosThetaStarZP;
    arr[iArr++] = CosThetaStarEP;
    arr[iArr++] = KplusEta;
    arr[iArr++] = KminusEta;
    arr[iArr++] = CosThetaStarRP;
    
    McPhiMeson->Fill(arr);
  }

  if(mMode != 3)
  {
    h_theta[0]->Fill(TMath::Abs(CosThetaStarZP));
    h_theta_star[0]->Fill(TMath::Abs(CosThetaStarEP));
    h_theta_star_RP[0]->Fill(TMath::Abs(CosThetaStarRP));

    h2_theta_star[0]->Fill(TMath::Abs(CosThetaStarEP),BetaP);
    h2_theta_star_RP[0]->Fill(TMath::Abs(CosThetaStarRP),Beta);
  }
  //cout << "Before:  K+ eta = " << KplusEta << "     K- eta = " << KminusEta << "     phi rapidity = " << PhiRapidity << endl;
//  if(TMath::Abs(KplusEta)<=eta_gap && TMath::Abs(KminusEta)<=eta_gap)
   // cout << "AFTER:  K+ eta = " << KplusEta << "     K- eta = " << KminusEta << "     phi rapidity = " << PhiRapidity << endl;
  if(mMode != 3)
  {
    for(int ie = 1; ie < 6; ie++)
    {
      if(TMath::Abs(KplusEta)<=eta_gap[ie-1] && TMath::Abs(KminusEta)<=eta_gap[ie-1])// && KplusPt > 0.1 && KminusPt > 0.1 && KminusPt < 10.0 && KplusPt < 10.0)
      {
        h_theta[ie]->Fill(TMath::Abs(CosThetaStarZP));
        h_theta_star[ie]->Fill(TMath::Abs(CosThetaStarEP));
        h_theta_star_RP[ie]->Fill(TMath::Abs(CosThetaStarRP));

        h2_theta_star[ie]->Fill(TMath::Abs(CosThetaStarEP),BetaP);
        h2_theta_star_RP[ie]->Fill(TMath::Abs(CosThetaStarRP),Beta);
 
        h_cosbeta[ie]->Fill(TMath::Cos(Beta));
        h_cosbetaP[ie]->Fill(TMath::Cos(BetaP));
        h_cos2beta[ie]->Fill(TMath::Cos(2.0*Beta));
        h_cos2betaP[ie]->Fill(TMath::Cos(2.0*BetaP));
        h_cos4beta[ie]->Fill(TMath::Cos(4.0*Beta));
        h_cos4betaP[ie]->Fill(TMath::Cos(4.0*BetaP));

        p_cosbeta[ie]->Fill(TMath::Abs(CosThetaStarRP),TMath::Cos(Beta));
        p_cosbetaP[ie]->Fill(TMath::Abs(CosThetaStarEP),TMath::Cos(BetaP));
        p_cos2beta[ie]->Fill(TMath::Abs(CosThetaStarRP),TMath::Cos(2.0*Beta));
        p_cos2betaP[ie]->Fill(TMath::Abs(CosThetaStarEP),TMath::Cos(2.0*BetaP));
        p_cos4beta[ie]->Fill(TMath::Abs(CosThetaStarRP),TMath::Cos(4.0*Beta));
        p_cos4betaP[ie]->Fill(TMath::Abs(CosThetaStarEP),TMath::Cos(4.0*BetaP));

        h_beta[ie]->Fill(Beta);    
        h_tstar[ie]->Fill(TMath::ACos(CosThetaStarRP));   
                        
        h_betaP[ie]->Fill(BetaP);   
        h_tstarP[ie]->Fill(TMath::ACos(CosThetaStarEP));  
                             
        h_betatstar[ie]->Fill(TMath::ACos(CosThetaStarRP),Beta);                           
        h_betatstarP[ie]->Fill(TMath::ACos(CosThetaStarEP),BetaP);   
      }
    }
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
  File_OutPut->cd();

  if(mMode != 3)
  {
    for(int ie = 0; ie < 6; ie++)
    {
      h_theta[ie]->Write();
      h_theta_star[ie]->Write();
      h_theta_star_RP[ie]->Write();

      h2_theta_star[ie]->Write();
      h2_theta_star_RP[ie]->Write();

      h_cosbeta[ie]->Write();  
      h_cos2beta[ie]->Write();  
      h_cos4beta[ie]->Write();  
      h_cosbetaP[ie]->Write(); 
      h_cos2betaP[ie]->Write(); 
      h_cos4betaP[ie]->Write(); 

      p_cosbeta[ie]->Write();  
      p_cos2beta[ie]->Write();  
      p_cos4beta[ie]->Write();  
      p_cosbetaP[ie]->Write(); 
      p_cos2betaP[ie]->Write(); 
      p_cos4betaP[ie]->Write(); 

      h_beta[ie]->Write();     
      h_tstar[ie]->Write();    
                      
      h_betaP[ie]->Write();    
      h_tstarP[ie]->Write();   
                         
      h_betatstar[ie]->Write();   
      h_betatstarP[ie]->Write();  
    }
  }

  if(mMode == 3) McPhiMeson->Write();
  //h_out1->Write();
  //h_out2->Write();
  File_OutPut->Close();

  cos2r2  /= mNEvents;
  cos2del /= mNEvents;
  cosdel  /= mNEvents;
  sin2del /= mNEvents;
  sindel  /= mNEvents;

  cout << "cos(2*(Psi2-RP)) = " << cos2r2 << endl;
  cout << "cos(2*delta) = " << cos2del << endl;
  cout << "cos(delta)   = " << cosdel << endl;
  cout << "sin(2*delta) = " << sin2del << endl;
  cout << "sin(delta)   = " << sindel << endl;
  
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


