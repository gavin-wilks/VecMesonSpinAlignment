#include "StEffMcPhiHelicityGlobal.h"
#include "StEffHistMangerHelicityGlobal.h"
#include "StEffCut.h"
#include <string>
#include "TFile.h"
#include "TNtuple.h"
#include "TBranch.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRotation.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StRoot/Utility/functions.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TRotation.h"
#include <algorithm>
#include <unordered_set>

typedef std::map<std::string,TH3F*> TH3FMap;
typedef std::map<std::string,TH1D*> TH1DMap;
typedef std::map<std::string,TF1*> TF1Map;
TH3FMap h_m2_PID;
TH3FMap h_nsig_PID;

TH1DMap h_EffKplus;
TH1DMap h_EffKminus;
TH1D *h_FrameEta[2];
TH1D *h_FramePhi[2];

void readEfficiency(int energy)
{
  string inputKplus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/TPC/Eff_%s_%s_SingleKaon.root",vmsa::mParType[0].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Kplus = TFile::Open(inputKplus.c_str());
  cout << "OPEN Efficiency File for K+: " << inputKplus.c_str() << endl;

  string inputKminus = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/TPC/Eff_%s_%s_SingleKaon.root",vmsa::mParType[1].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Kminus = TFile::Open(inputKminus.c_str());
  cout << "OPEN Efficiency File for K-: " << inputKminus.c_str() << endl;

  //h_FrameEta[0] = (TH1D*)File_Kplus->Get("h_FrameEta");
  h_FrameEta[0] = new TH1D("h_FrameEta0","h_FrameEta0",vmsa::BinEta,-vmsa::mEtaMax,vmsa::mEtaMax);
  h_FramePhi[0] = (TH1D*)File_Kplus->Get("h_FramePhi");
  h_FrameEta[1] = new TH1D("h_FrameEta1","h_FrameEta1",vmsa::BinEta,-vmsa::mEtaMax,vmsa::mEtaMax);
  //h_FrameEta[1] = (TH1D*)File_Kminus->Get("h_FrameEta");
  h_FramePhi[1] = (TH1D*)File_Kminus->Get("h_FramePhi");

  h_FrameEta[0]->SetDirectory(0); 
  h_FramePhi[0]->SetDirectory(0);
  h_FrameEta[1]->SetDirectory(0);
  h_FramePhi[1]->SetDirectory(0);

  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
    {
      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
	string KEY = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_EffKplus[KEY] = (TH1D*)File_Kplus->Get(KEY.c_str());
	h_EffKminus[KEY] = (TH1D*)File_Kminus->Get(KEY.c_str());
        //h_EffKplus[KEY]->Print(); 
        //h_EffKminus[KEY]->Print();
        h_EffKplus[KEY]->SetDirectory(0);; 
        h_EffKminus[KEY]->SetDirectory(0);;
      }	
    }
  }
}

void findHist(TLorentzVector const& lKaon, int iParticleIndex, int& EtaBin, int& PhiBin)
{
  float eta = lKaon.Eta();
  //cout << "eta = " << eta  << endl;
  EtaBin = h_FrameEta[iParticleIndex]->FindBin(eta)-1;

  //float phi = lKaon.Phi()-Psi2;
  //float phi_shift = AngleShift(phi);
  //PhiBin = h_FramePhi[iParticleIndex]->FindBin(phi_shift)-1;

  float phi = lKaon.Phi();
  while(phi < -TMath::Pi()) phi += 2.*TMath::Pi();
  while(phi >  TMath::Pi()) phi += 2.*TMath::Pi();
  PhiBin = h_FramePhi[iParticleIndex]->FindBin(phi)-1;
}

bool tpcReconstructed(int iParticleIndex, int cent, TLorentzVector const& lKaon)
{
   if(fabs(lKaon.Eta()) >= vmsa::mEtaMax) return false;
    
   //cout << "Before FindHist" << endl;
    
   TH1D *h_TPC = NULL;
   int EtaBin_TPC = -1;
   int PhiBin_TPC = -1;
   findHist(lKaon,iParticleIndex,EtaBin_TPC,PhiBin_TPC);

   //cout << "Eta bin = " << EtaBin_TPC << ", Phi bin = " << PhiBin_TPC << endl;
   //cout << "After FindHist" << endl;

   if (iParticleIndex == 0)
   {
     string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
     h_TPC = h_EffKplus[KEY_TPC];
   }
   else
   {
     string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
     h_TPC = h_EffKminus[KEY_TPC];
   }
   //cout << "Grab Hist" << endl;

   //h_TPC->Print();

   double pt = lKaon.Perp();
   if(pt < 0.1) return false;
   int const bin_TPC = h_TPC->FindBin(pt);
   bool is_TPC = gRandom->Rndm() < h_TPC->GetBinContent(bin_TPC);
   //cout << "Finished function" << endl;

   //delete h_TPC;
   return is_TPC;
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



ClassImp(StEffMcPhiHelicityGlobal)

int StEffMcPhiHelicityGlobal::mInput_flag = 1;

bool Sampling(TF1 *f_rhoPhy, float CosThetaStar, float max);
bool SamplingHelicity(TF1 *f_rhoPhy, float CosThetaStar, float max);
//bool Sampling2D(TF2 *f_rhoPhy, float ThetaStar, float beta, float max);
bool Sampling2D(TF2 *f_rhoPhy, float ThetaStar, float beta, float max);
bool Sampling2DWeight(TF2 *f_rhoPhy, float ThetaStar, float beta, float max);

TF1* readv2(int energy, int pid, int centrality){
  //string centFile[9] = {"4080","4080","4080","4080","1040","1040","1040","0010","0010"};
  int tableNumV2[5][9] = {{0,0,0,0,0,0,0,0,0},
                        {0,0,0,0,0,0,0,0,0},
                        {0,0,0,0,0,0,0,0,0},
                        {0,0,0,0,0,0,0,0,0},
                        {239,239,239,239,141,141,141,43,43}};
  
  string centlabel = "4080";
  if(centrality >= 4 && centrality <= 6) centlabel = "1040";
  if(centrality >= 7 && centrality <= 8) centlabel = "0010";
  TGraphAsymmErrors *g_v2;
  if(energy != 3)
  {
    string InPutV2 = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
    if((energy == 2 || energy == 0) ) InPutV2 = "/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/v2/HEPData-ins1395151-v2-root.root";
    if(energy == 4 ) InPutV2 = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/v2/OutPhi_v2_Cent%s.root",centlabel.c_str());
    TFile *File_v2 = TFile::Open(InPutV2.c_str());
    std::cout << "v2 file: " << InPutV2 << endl;

    g_v2 = (TGraphAsymmErrors*)File_v2->Get("g_v2");

    if((energy == 2 || energy == 0) ) 
    {
      TDirectory *dir = (TDirectory*) File_v2->Get(Form("Table %d",tableNumV2[energy][centrality]));
      dir->cd(); 
      g_v2 = (TGraphAsymmErrors*)dir->Get("Graph1D_y1");
    }
    if(energy == 4) 
    {
      g_v2 = (TGraphAsymmErrors*) File_v2->Get("Graph");
      g_v2->Print();
    }
  }
  
  //if( energy == 3 && mMode == 1 )
  //{
  //  int centidx = 0;
  //  if(centrality >= 0 && centrality <= 3) centidx = 3;
  //  if(centrality >= 4 && centrality <= 6) centidx = 2;
  //  if(centrality >= 7 && centrality <= 8) centidx = 1;
  //  g_v2 = new TGraphAsymmErrors();
  //  for(int ipt = 0; ipt < phiv2_14::ptbins; ipt++)
  //  {
  //    g_v2->SetPoint(ipt, phiv2_14::pt[centidx][ipt], phiv2_14::v2[centidx][ipt]);
  //    double stat = phiv2_14::stat[centidx][ipt];
  //    double sys  = phiv2_14::sys[centidx][ipt];
  //    double totalerr = TMath::Sqrt( stat*stat + sys*sys );
  //    g_v2->SetPointError(ipt, 0.0, 0.0, totalerr, totalerr);  
  //  }

  //}


  TF1 *f_v2 = new TF1("f_v2",v2_pT_FitFunc,vmsa::ptMin,vmsa::ptMax,5);
  f_v2->FixParameter(0,2);
  f_v2->SetParameter(1,0.1);
  f_v2->SetParameter(2,0.1);
  f_v2->SetParameter(3,0.1);
  f_v2->SetParameter(4,0.1);
  cout << "Fitting v2" << endl;
  g_v2->Fit(f_v2,"N");

  return f_v2;
}

StEffMcPhiHelicityGlobal::StEffMcPhiHelicityGlobal(int Energy, long StartEvent, long StopEvent, int PID, int year, int mode, int inputpt, int startpt, int stoppt, const char* setting, int etamode, int order = 2, float realrho1n1 = 0.0, float rho00 = 1./3., float reterms = 0.0, float imterms = 0.0, float imagrho1n1 = 0.0, float rhohelicity = 1./3., float ptfixed = 0.0, float yfixed = 0.0) 
{
  mOrder = order;
  energy = Energy;
  pid = PID;
  mInputPt = inputpt;
  mStartPt = startpt;
  mStopPt = stoppt;
  mEtaMode = etamode;
  mMode = mode;
  rerho1n1 = realrho1n1;
  imrho1n1 = imagrho1n1;
  real = reterms;
  imag = imterms;
  mrho00 = rho00;

  mrho00helicity = rhohelicity; 

  for(int i = 0; i < 9; i++)
  {
    f_mV2[i] = readv2(energy,pid,i);
  }
  
  std::string EP[2] = {"","2nd"};

  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_prelimv2_EPeff_randomRP/%s_%s_yabs1_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_PhiEffFinerBins_prelimv2_flatRP_yabs1p1/%s_%s_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEff/Phi_pt%d_y%d.root",ptfixed,yfixed);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/GlobalHelicityStudies/19GeV/Phi_pt%1.1f_y%1.1f.root",ptfixed,yfixed);
  string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/GlobalHelicityStudies/19GeV/Phi_pt%1.1f_y%1.1f_pt%d_set%d/cent%d_%d.root",ptfixed,yfixed,inputpt,year,startpt,stoppt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_nov2_flatRP_yabs1_PhiEmbed/%s_%s_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_PhiEffFinerBins_nov2_flatRP_yabs1p1_flatpt/%s_%s_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_nov2_flatRP_yabs1p1/%s_%s_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_PhiEff/%s_%s_yabs1_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_v2times3/%s_%s_yabs1_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_EP/%s_%s_eta1_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);

  SetInPutFile(InPutFile); // set input list

  SetStartEvent(StartEvent); // set start event
  SetStopEvent(StopEvent); // set stop event

  string OutPutFile = Form("Eff_%s_SingleParticle_%s_Mode%d_EtaMode%d_pt%1.1f_y%1.1f_cent%d_%d.root",vmsa::mBeamEnergy[energy].c_str(),setting,mode,etamode,ptfixed,yfixed,startpt,stoppt);
  SetOutPutFile(OutPutFile); // set output file

  mEffCut = new StEffCut();
  mEffHistManger = new StEffHistMangerHelicityGlobal(energy, pid, mode, startpt, stoppt, ptfixed, yfixed);

  cout << "Created histogram manager" << endl;

  ////////////// LOAD 1ST ORDER EP INFORMATION //////////////////////////////
  //TF1* f_res = new TF1("resolution",EventPlaneResolution,0,80,0);

  //TString InPutFile_Res1 = Form("StRoot/Utility/EpdResolution/Resolution_file_%s_EpdCorrections_4.root",vmsa::mBeamEnergy[energy].c_str());
  //mInPutFile_Res1 = TFile::Open(InPutFile_Res1.Data());

  //TProfile *p_res1 = (TProfile*)mInPutFile_Res1->Get("AveCosDeltaPsi1");
  //for(int icent = 0; icent < 9; icent++)
  //{
  //  float Res_raw1 = p_res1->GetBinContent(p_res1->FindBin(icent));
  //  float Res1 = TMath::Sqrt(Res_raw1);
  //  float Chi1 = f_res->GetX(Res1); // This is for sub  event plane resolution
  //  Chi1 *= TMath::Sqrt(2.0);  // This is for full event plane resolution
  //  cout << "Centrality = " << icent
  //       << ",   Resolution1 = " << Res1
  //       << ",   Chi1 = " << Chi1 << endl;

  //  //Now store chi values
  //  mChi[0][icent] = Chi1;
  //}
  //f_pDel1 = new TF1("deltaPsi1",EventPlaneDist1st,-TMath::Pi(),TMath::Pi(),2);
  //f_pDel1->FixParameter(1,1.0/(2.0*TMath::Pi()));
  ////////////// LOAD 1ST ORDER EP INFORMATION //////////////////////////////

  ////////////// LOAD 2ND ORDER EP INFORMATION //////////////////////////////
  //TString InPutFile_Res2 = Form("StRoot/Utility/Resolution/file_%s_Resolution.root",vmsa::mBeamEnergy[energy].c_str());
  //mInPutFile_Res2 = TFile::Open(InPutFile_Res2.Data());

  //TProfile *p_res2 = (TProfile*)mInPutFile_Res2->Get("p_mRes2_Sub");
  //for(int icent = 0; icent < 9; icent++)
  //{
  //  float Res_raw2 = p_res2->GetBinContent(p_res2->FindBin(icent));
  //  float Res2 = TMath::Sqrt(Res_raw2);
  //  float Chi2 = f_res->GetX(Res2);

  //  cout << "Centrality = " << icent
  //       << ",   Resolution2 = " << Res2
  //       << ",   Chi2 = " << Chi2 << endl;

  //  //Now store chi values
  //  mChi[1][icent] = Chi2;
  //}
  //f_pDel2 = new TF1("deltaPsi2",EventPlaneDist,-TMath::Pi()/2.0,TMath::Pi()/2.0,2);
  //f_pDel2->FixParameter(1,1.0/(2.0*TMath::Pi()));
  ////////////// LOAD 2ND ORDER EP INFORMATION //////////////////////////////

 // cout << "Trying to load rho" << endl;
 // std::string inputfilept     = Form("StRoot/Utility/Rho/%s/Rho_AccResSysErrors_F_0_eta1_eta1_PolySys_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
 // if(mOrder == 1) inputfilept = Form("StRoot/Utility/Rho/%s/Rho_AccResSysErrors_F_0_eta1_eta1_PolySys_FirstOrder_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
 // std::string inputfilecent;
 // std::string inputfiley;
 // if(mOrder == 1) inputfilecent = Form("StRoot/Utility/Rho/%s/RhoCent_AccResSysErrors_eta1_eta1_PolySys_FirstOrder_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
 // if(mOrder == 2) inputfilecent = Form("StRoot/Utility/Rho/%s/RhoCent_AccResSysErrors_eta1_eta1_PolySys_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
 // if(mOrder == 1) inputfiley = Form("StRoot/Utility/Rho/%s/RhoEta_AccResSysErrors_eta1_eta1_PolySys_FirstOrder_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
 // if(mOrder == 2) inputfiley = Form("StRoot/Utility/Rho/%s/RhoEta_AccResSysErrors_eta1_eta1_PolySys_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());

 // //std::string inputfilecent = "StRoot/Utility/Rho/";
 // //std::string inputfiley = "StRoot/Utility/Rho/";
 // 
 // cout << "Trying to open files" << endl;
 // TFile *filept   = TFile::Open(inputfilept.c_str());
 // TFile *filecent = TFile::Open(inputfilecent.c_str());
 // TFile *filey    = TFile::Open(inputfiley.c_str());
 // cout << "Opening files" << endl;
 

 // TGraphAsymmErrors *g_rho_pt = (TGraphAsymmErrors*) filept->Get(Form("g_rho00_order%d_%s_%s_StatError",mOrder,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()));
 // double helicityrho00[6] = {0.333333,0.333333,0.307032,0.295729,0.288274,0.246577};

  for(int j = 0; j < 5; j++)
  {
    f_mRhoPt_2D[0][j] = new TF2(Form("f_mRho2D_0_%d",j),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),5);
    double rhoinput = 1./3.+(j-2)*0.1;
    cout << "rho input = " << rhoinput << endl;
    f_mRhoPt_2D[0][j]->FixParameter(0,rhoinput);
    f_mRhoPt_2D[0][j]->FixParameter(1,0.0);
    f_mRhoPt_2D[0][j]->FixParameter(2,0.0);
    f_mRhoPt_2D[0][j]->FixParameter(3,0.0);
    f_mRhoPt_2D[0][j]->FixParameter(4,0.0);
  }
  for(int irho = 1; irho < 5; irho++)
  {
    for(int j = 0; j < 5; j++)
    {
      f_mRhoPt_2D[irho][j] = new TF2(Form("f_mRho2D_%d_%d",irho,j),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),5);
      double rhoinput = (j-2)*0.1;
      cout << "rho input = " << rhoinput << endl;
      f_mRhoPt_2D[irho][j]->FixParameter(0,1./3.);
      f_mRhoPt_2D[irho][j]->FixParameter(1,0.0);
      f_mRhoPt_2D[irho][j]->FixParameter(2,0.0);
      f_mRhoPt_2D[irho][j]->FixParameter(3,0.0);
      f_mRhoPt_2D[irho][j]->FixParameter(4,0.0);
    
      f_mRhoPt_2D[irho][j]->FixParameter(irho,rhoinput);
    }
  }

 // for(int ipt = 2; ipt < 6; ipt++)
 // { 
 //   f_mRhoPt[ipt] = new TF1(Form("f_mRho_pt%d",ipt),SpinDensity,-1.0,1.0,2);
 //   f_mRhoPt_Helicity[ipt] = new TF1(Form("f_mRhoHelicity_pt%d",ipt),SpinDensity,-1.0,1.0,2);
 //   f_mRhoPt_2D[ipt] = new TF2(Form("f_mRho2D_pt%d",ipt),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),5);
 //   TF2 *temp = new TF2("temp",SpinDensity2Dcosneg,0.0,TMath::Pi(),0.0,2.0*TMath::Pi(),5);
 //   double pt, rho; 
 //   g_rho_pt->GetPoint(ipt,pt,rho);
 //   cout << "ipt = " << ipt << ", pt = " << pt << ", rho = " << rho << endl;
 //   //cout << "ipt = " << ipt << ", pt = " << pt << ", rho = " << rho << endl;
 //   //f_mRhoPt[ipt]->FixParameter(0,0.5);
 //   //f_mRhoPt[ipt]->FixParameter(0,rho); // set by data
 //   f_mRhoPt[ipt]->FixParameter(0,mrho00);// set by user
 //   f_mRhoPt[ipt]->FixParameter(1,0.75);

 //   mMaxData[ipt] = f_mRhoPt[ipt]->GetMaximum(-1.0,1.0);
 //   cout << "Maximum of data rh00 = " << mMaxData[ipt] << endl; 

 //   //f_mRhoPt_Helicity[ipt]->FixParameter(0,helicityrho00[ipt]); // set by data
 //   f_mRhoPt_Helicity[ipt]->FixParameter(0,mrho00helicity);       // set by user
 //   f_mRhoPt_Helicity[ipt]->FixParameter(1,0.75);
 // 
 //   mMaxHelicity[ipt] = f_mRhoPt_Helicity[ipt]->GetMaximum(-1.0,1.0);
 //   cout << "Maximum of helicity rh00 = " << mMaxHelicity[ipt] << endl; 

 //   cout << "rho00 = " << mrho00 << ", rerho1n1 = " << rerho1n1 << endl;

 //   f_mRhoPt_2D[ipt]->FixParameter(0,mrho00);
 //   //f_mRhoPt_2D[ipt]->FixParameter(0,rho);
 //   f_mRhoPt_2D[ipt]->FixParameter(1,real);
 //   f_mRhoPt_2D[ipt]->FixParameter(2,imag);
 //   f_mRhoPt_2D[ipt]->FixParameter(3,rerho1n1);
 //   f_mRhoPt_2D[ipt]->FixParameter(4,imrho1n1);

 //   temp->FixParameter(0,mrho00);
 //   temp->FixParameter(1,real);
 //   temp->FixParameter(2,imag);
 //   temp->FixParameter(3,rerho1n1);
 //   temp->FixParameter(4,imrho1n1);

 //   //mMax = f_mRhoPt_2D[ipt]->GetMaximum();
 //   //double x[2] = {TMath::Pi()/2.,TMath::Pi()/2.};
 //   double x,y;
 //   temp->GetMinimumXY(x,y);
 //   cout << "Maximum x,y = " << x << "," << y << endl;
 //   mMax = TMath::Abs(temp->Eval(x,y));
 //   cout << "Maximum of the function is " << mMax << endl;

 //   //cout << "x,y both = TMath::Pi()/2, eval func = " << temp->Eval(1.57113,1.5708) << endl;
 // }  
 
  //`for(int ipt = 0; ipt < vmsa::pt_rebin_cent; ipt++)
  //`{
  //`  TGraphAsymmErrors *g_rho_cent = (TGraphAsymmErrors*) filecent->Get(Form("rhoRawStat_pt_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",ipt,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str()));
  //`  for(int icent = 0; icent < 9; icent++)
  //`  {
  //`    f_mRhoCent[ipt][icent] = new TF1(Form("f_mRho_pt%d_cent%d",ipt,icent),SpinDensity,-1.0,1.0,2);
  //`    double cent, rho; 
  //`    g_rho_cent->GetPoint(icent,cent,rho);
  //`    cout << "ipt = " << ipt << ", cent = " << cent << ", rho = " << rho << endl;
  //`    f_mRhoCent[ipt][icent]->FixParameter(0,rho);
  //`    f_mRhoCent[ipt][icent]->FixParameter(1,0.75);
  //`  }
  //`}    

  //`for(int ipt = 0; ipt < vmsa::pt_rebin_y; ipt++)
  //`{
  //`  for(int icent = 0; icent < vmsa::cent_rebin_total; icent++)
  //`  {
  //`    TGraphAsymmErrors *g_rho_y = (TGraphAsymmErrors*) filey->Get(Form("rhoRawStat_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,icent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1));
  //`    for(int iy = 0; iy < 10; iy++)
  //`    { 
  //`      f_mRhoY[ipt][icent][iy] = new TF1(Form("f_mRho_pt%d_cent%d_y%d",ipt,icent,iy),SpinDensity,-1.0,1.0,2);
  //`      double y, rho; 
  //`      g_rho_y->GetPoint(iy,y,rho);
  //`      cout << "ipt = " << ipt << ", icent = " << icent << ", y = " << y << ", rho = " << rho << endl;
  //`      f_mRhoY[ipt][icent][iy]->FixParameter(0,rho);
  //`      f_mRhoY[ipt][icent][iy]->FixParameter(1,0.75);
  //`    }
  //`  }
  //`}    
 
}

StEffMcPhiHelicityGlobal::~StEffMcPhiHelicityGlobal()
{
}

//------------------------------------------------------------
void StEffMcPhiHelicityGlobal::SetInPutFile(const string inputfile)
{
  mInPutFile = inputfile;
  cout << "Input file was set to: " << mInPutFile.c_str() << endl;
}

void StEffMcPhiHelicityGlobal::SetOutPutFile(const string outputfile)
{
  mOutPutFile = outputfile;
  cout << "Output file was set to: " << mOutPutFile.c_str() << endl;
}

void StEffMcPhiHelicityGlobal::SetStartEvent(const long StartEvent)
{
  mStartEvent = StartEvent;
  cout << "nStartEvent = " << mStartEvent << endl;
}

void StEffMcPhiHelicityGlobal::SetStopEvent(const long StopEvent)
{
  mStopEvent = StopEvent;
  cout << "nStopEvent = " << mStopEvent << endl;
}
//------------------------------------------------------------

void StEffMcPhiHelicityGlobal::Init()
{
  mEffHistManger->InitHist(mBinCos,mBinPhi);
  //mEffHistManger->InitKaonHist();
  //mEffHistManger->InitPhiHist();

  readm2PID_hist();
  cout << "Read in m2 PID " << endl;
  readnsigPID_hist();
  cout << "Read in nsig PID " << endl;
  readEfficiency(energy);
  cout << "Read in TPC Efficiency " << endl;

  ToFFile = TFile::Open(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/ToFMatching/ToFMatching_19GeV_dipangle.root"));
  ToFHist[0] = (TH3F*) ToFFile->Get("kplus_ratio");
  ToFHist[0]->Print();
  ToFHist[1] = (TH3F*) ToFFile->Get("kminus_ratio");
  ToFHist[1]->Print();
  cout << "Read in TOF Efficiency " << endl;


  // initialize the TNtuple
  cout << "Right before inputfile" << endl;
  cout << mInPutFile << endl;
  mFile_InPut = TFile::Open(mInPutFile.c_str());
  cout << "OPEN InPut File: " << mInPutFile.c_str() << endl;

  if(pid == 0) mNtuple = (TNtuple*)mFile_InPut->Get("McPhiMeson");

  // initialize Ntuple
  mNtuple->SetBranchAddress("Centrality",&mCentrality);
  mNtuple->SetBranchAddress("PsiRP",&mPsi);
  mNtuple->SetBranchAddress("Psi1",&mPsi1);
  mNtuple->SetBranchAddress("Psi2",&mPsi2);
  mNtuple->SetBranchAddress("McPt",&mMcPt);
  mNtuple->SetBranchAddress("McP",&mMcP);
  mNtuple->SetBranchAddress("McEta",&mMcEta);
  mNtuple->SetBranchAddress("McY",&mMcY);
  mNtuple->SetBranchAddress("McPhi",&mMcPhi);
  mNtuple->SetBranchAddress("McInvMass",&mMcInvMass);
  mNtuple->SetBranchAddress("McPid",&mMcPid);

  mNtuple->SetBranchAddress("KpMcPt",&mKpMcPt);
  mNtuple->SetBranchAddress("KpMcEta",&mKpMcEta);
  mNtuple->SetBranchAddress("KpMcY",&mKpMcY);
  mNtuple->SetBranchAddress("KpMcPhi",&mKpMcPhi);
  mNtuple->SetBranchAddress("KpMcM",&mKpMcM);
  mNtuple->SetBranchAddress("KpMcPid",&mKpMcPid);

  mNtuple->SetBranchAddress("KpRcPt",&mKpRcPt);
  mNtuple->SetBranchAddress("KpRcEta",&mKpRcEta);
  mNtuple->SetBranchAddress("KpRcY",&mKpRcY);
  mNtuple->SetBranchAddress("KpRcPhi",&mKpRcPhi);
  mNtuple->SetBranchAddress("KpRcM",&mKpRcM);
  mNtuple->SetBranchAddress("KpRcTpc",&mKpRcTpc);
  mNtuple->SetBranchAddress("KpRcTof",&mKpRcTof);
  
  mNtuple->SetBranchAddress("KmMcPt",&mKmMcPt);
  mNtuple->SetBranchAddress("KmMcEta",&mKmMcEta);
  mNtuple->SetBranchAddress("KmMcY",&mKmMcY);
  mNtuple->SetBranchAddress("KmMcPhi",&mKmMcPhi);
  mNtuple->SetBranchAddress("KmMcM",&mKmMcM);
  mNtuple->SetBranchAddress("KmMcPid",&mKmMcPid);

  mNtuple->SetBranchAddress("KmRcPt",&mKmRcPt);
  mNtuple->SetBranchAddress("KmRcEta",&mKmRcEta);
  mNtuple->SetBranchAddress("KmRcY",&mKmRcY);
  mNtuple->SetBranchAddress("KmRcPhi",&mKmRcPhi);
  mNtuple->SetBranchAddress("KmRcM",&mKmRcM);
  mNtuple->SetBranchAddress("KmRcTpc",&mKmRcTpc);
  mNtuple->SetBranchAddress("KmRcTof",&mKmRcTof);

  mNtuple->SetBranchAddress("RcPt",&mRcPt);
  mNtuple->SetBranchAddress("RcP",&mRcP);
  mNtuple->SetBranchAddress("RcEta",&mRcEta);
  mNtuple->SetBranchAddress("RcY",&mRcY);
  mNtuple->SetBranchAddress("RcPhi",&mRcPhi);
  mNtuple->SetBranchAddress("RcInvMass",&mRcInvMass);

  int num_tracks = mNtuple->GetEntriesFast();
  cout << "Number of tracks in McPhiMeson = " << num_tracks<< endl;

  if(mStartEvent > num_tracks) mStartEvent = num_tracks;
  if(mStopEvent  > num_tracks) mStopEvent  = num_tracks;
  cout << "New nStartEvent = " << mStartEvent << ", new nStopEvent = " << mStopEvent << endl;

  mFile_OutPut = new TFile(mOutPutFile.c_str(),"RECREATE");
  //f_y = new TF1("f_y",Form("exp(-x*x/2/%s)",mSigmay.c_str()),-1.0,1.0);

  const double  pythialow[19] = { 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.6, 4.2};;
  const double pythiahigh[19] = { 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.6, 4.2, 5.0};;

  double norm[19] = {1.14659, 1.15282, 1.16102, 1.16754, 1.17421, 1.20569, 1.18756, 1.20748, 1.20482, 1.24736, 1.23579, 1.20088, 1.21209, 1.26302, 1.22473, 1.29652, 1.29134, 1.67825, 1.67232};
  double mean[19] = {0.000950466, -0.00205036, 0.00635211, -0.00776218, -0.00329925, -0.00254986, 0.0106082, 0.00417584, 0.0110738, 0.0024981, 0.0287173, -0.0166698, -6.69877e-05, -0.0367374, -0.0154559, 0.0763402, -0.0882059, 0.113747, -0.142898};
  double sigma[19] = {0.959611, 0.94694, 0.929603, 0.912229, 0.899472, 0.838369, 0.875892, 0.843554, 0.839229, 0.771095, 0.791428, 0.844824, 0.825385, 0.743799, 0.787513, 0.698393, 0.707714, 0.425395, 0.742268};

  for(int i = 0; i < 19; i++)
  { 
    pythiaflat[i]= new TF1(Form("pythiaflat_%d",i),"[0]*exp(-(x-[1])*(x-[1])/2/[2]/[2])",-1.0,1.0);
    pythiaflat[i]->SetParameter(0,norm[i]);
    pythiaflat[i]->SetParameter(1,mean[i]);
    pythiaflat[i]->SetParameter(2,sigma[i]);
  }


}

void StEffMcPhiHelicityGlobal::Make()
{
  long totalmc = 0; 
  long totalprocessed = 0; 
  long totalbeforemc = 0; 

  long start_event_use = mStartEvent;
  long stop_event_use  = mStopEvent;
  
  double v2_7GeV[9] = {0.125818,0.125818,0.125818,0.125818,0.039951,0.039951,0.039951,0.013462,0.013462};

  int yptstart[4] = {4,10,13,16};
  int yptstop[4] = {9,12,15,17};

  float ptbinedges[5] = {1.2,1.8,2.4,3.0,4.2};

  gRandom->SetSeed();
  mNtuple->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry
  TF1* f_flow[11];

  TF1 *f_flowreal = new TF1("f_flowreal",flowSampleNorm,-TMath::Pi(),TMath::Pi(),1);
  for(int i = 0; i < 11; i++)
  {
    double v2input = 1./30.*i;
    cout << "v2 input = " << v2input << endl;
    f_flow[i] = new TF1("f_flow",flowSampleNorm,-TMath::Pi(),TMath::Pi(),1);
    f_flow[i]->SetParameter(0,v2input);
  }
 
  //string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/ptyspectra_datarcratio_%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  //TFile *File_Spectra = TFile::Open(inputfile.c_str());
  //TH2D *h_mSpectraRatio = (TH2D*) ((TH2D*) File_Spectra->Get("pty_datarcratio"))->Clone();
  //double spectramax = h_mSpectraRatio->GetMaximum();
  //h_mSpectraRatio->Scale(1./spectramax);
  //cout << "spectra max = " << spectramax << ", normalzied max " << h_mSpectraRatio->GetMaximum() << endl;

  //int eventtarget[4] = {1282940,301140,49281,7521};
  int eventtarget[4] = {1282940,301140,49281,100};
  int eventcounts[5][5] = {0};

  //std::vector<Long64_t> indices(stop_event_use);
  //for (Long64_t i = 0; i < stop_event_use; ++i) {
  //    indices[i] = i;
  //}
  std::unordered_set<Long64_t> usedIndices;

  //double maxRho[5][5] = {0.0};
  //for(int irho = 0; irho < 5; irho++)
  //{
  //  for(int j = 0; j < 5; j++)
  //  {
  //    maxRho[irho][j] = f_mRhoPt_2D[irho][j]->GetMaximum();
  //  }
  //}
  //double maxRapidity[19] = {0.0};
  //for(int i = 0; i < 19; i++) maxRapidity[i] = pythiaflat[i]->GetMaximum();
  //std::random_shuffle(indices.begin(), indices.end(), [&](int stop_event_use) { return gRandom->Integer(stop_event_use); });

  TRandom3 randomGen;

  for(long i_track = start_event_use; i_track < stop_event_use; ++i_track) 
  { 
    //Long64_t randomIndex = randomGen.Integer(stop_event_use);
    
    //if(usedIndices.find(randomIndex) == usedIndices.end()) 
    //{
    //  if(!mNtuple->GetEntry(randomIndex)) break;
    //  //cout << "About to grab a random index" << endl;
    //  usedIndices.insert(randomIndex);
    //  //cout << "Grabbed a random index" << endl;
    //}
    //else 
    //{
    //  continue;
    //}

    if (!mNtuple->GetEntry(i_track)) // take track information
      break;  // end of data chunk
    if (floor(10.0*i_track/ static_cast<float>(stop_event_use)) > floor(10.0*(i_track-1)/ static_cast<float>(stop_event_use)))
    {
      cout << "=> processing data: " << 100.0*i_track/ static_cast<float>(stop_event_use) << "%" << endl;
    }
    totalprocessed++;
    //if(i_track%1000 == 0)
    //{
    //  int sum = 0;
    //  for(int irho = 0; irho < 5; irho++)
    //  {
    //    for(int j = 0; j < 5; j++)
    //    {
    //      if(eventcounts[irho][j] > eventtarget[mInputPt]) sum++;
    //      cout << irho << " " << j << " " << eventcounts[irho][j] << endl; 
    //   }
    //  }
    //  if(sum == 25) break;
    //  cout << "Our check is okay " << sum << " < 25" << endl;
    //}

    //cout << "TRACK " << i_track << endl;
    McVecMeson McPhi; // initialize McPhi
    McPhi.Centrality = mCentrality;
    //if(mMode == 0 && (( McPhi.Centrality > 5 || McPhi.Centrality < 2 ))) continue;
    McPhi.Psi        = mPsi;
    McPhi.Psi1       = mPsi1;
    McPhi.Psi2       = mPsi2;
    //cout << "McPhi.Psi2 =" << McPhi.Psi2 << endl;;
    McPhi.McPt       = mMcPt;
    McPhi.McP        = mMcP;
    McPhi.McEta      = mMcEta;
    McPhi.McY        = mMcY;
    McPhi.McPhi      = mMcPhi;
    //cout << "phi = " << McPhi.McPhi << endl;
    McPhi.McInvMass  = mMcInvMass;
    McPhi.McPid      = mMcPid;

    McDecayDau McKP; // initialize McKplus
    McKP.McPt  = mKpMcPt;
    McKP.McEta = mKpMcEta;
    McKP.McY   = mKpMcY;
    McKP.McPhi = mKpMcPhi;
    McKP.McM   = mKpMcM;
    McKP.McPid = mKpMcPid;

    RcDecayDau RcKP; // initialize RcKplus
    RcKP.RcPt  = mKpRcPt;
    RcKP.RcEta = mKpRcEta;
    RcKP.RcY   = mKpRcY;
    RcKP.RcPhi = mKpRcPhi;
    RcKP.RcM   = mKpRcM;
    RcKP.RcTpc = mKpRcTpc;
    RcKP.RcTof = mKpRcTof;

    McDecayDau McKM; // initialize McKminus
    McKM.McPt  = mKmMcPt;
    McKM.McEta = mKmMcEta;
    McKM.McY   = mKmMcY;
    McKM.McPhi = mKmMcPhi;
    McKM.McM   = mKmMcM;
    McKM.McPid = mKmMcPid;

    RcDecayDau RcKM; // initialize RcKminus
    RcKM.RcPt  = mKmRcPt;
    RcKM.RcEta = mKmRcEta;
    RcKM.RcY   = mKmRcY;
    RcKM.RcPhi = mKmRcPhi;
    RcKM.RcM   = mKmRcM;
    RcKM.RcTpc = mKmRcTpc;
    RcKM.RcTof = mKmRcTof;

    RcVecMeson RcPhi; // initialize RcPhi
    RcPhi.RcPt      = mRcPt;
    RcPhi.RcP       = mRcP;
    RcPhi.RcEta     = mRcEta;
    RcPhi.RcY       = mRcY;
    RcPhi.RcPhi     = mRcPhi;
    RcPhi.RcInvMass = mRcInvMass;


    //f_flow->ReleaseParameter(0);
    ////f_flow->SetParameter(0,mv2);//f_mV2[McPhi.Centrality]->Eval(McPhi.McPt));
    //f_flow->SetParameter(0,f_mV2[int(McPhi.Centrality)]->Eval(McPhi.McPt));
    ////if(energy == 0) f_flow->SetParameter(0,v2_7GeV[int(McPhi.Centrality)]);//f_mV2[McPhi.Centrality]->Eval(McPhi.McPt));
    ////cout << "Set the parameter" << endl;
    //double maximum = f_flow->Eval(0.0);
    ////cout << "Maximum = " << maximum << endl;
    //double phiminusPsi = McPhi.McPhi-McPhi.Psi; 
    ////cout << "grabbed maximum and phiminusPsi = " << phiminusPsi<< endl;
    //while(phiminusPsi < -TMath::Pi()) phiminusPsi += 2.0*TMath::Pi();
    //while(phiminusPsi >  TMath::Pi()) phiminusPsi -= 2.0*TMath::Pi();
    ////cout << "Angle Wrapping = " << phiminusPsi << endl;
    //double valuev2 = f_flow->Eval(phiminusPsi);
    //double randomv2 = gRandom->Uniform(0,1);
    //if(randomv2 > valuev2/maximum) continue; 

    //cout << "maximum = " << maximum << endl;
    //cout << "valuev2/maximum = " << valuev2/maximum << endl;    

    //cout << "mMode = " << mMode << endl;
    //cout << "Centrality = " << McPhi.Centrality << endl;

    //cout << "McPhi.McPt = " << McPhi.McPt << ", RcPhi.RcPt = " << RcPhi.RcPt << endl;
    //cout << "McPhi.McY = " << McPhi.McY << ", RcPhi.RcY = " << RcPhi.RcY << endl;

/* Gavin comment out 11/09/2024
    if( mMode == 0 && (McPhi.McPt <= vmsa::pt_low[energy][mStartPt] || McPhi.McPt > vmsa::pt_up[energy][mStopPt]))    continue;
    //if( mMode == 0 && (RcPhi.RcPt <= vmsa::pt_low[energy][mStartPt] || RcPhi.RcPt > vmsa::pt_up[energy][mStopPt]))    continue; //Swap to RC level cut
    if( mMode == 3 && (McPhi.McPt <= vmsa::pt_low[energy][mStartPt] || McPhi.McPt > vmsa::pt_up[energy][mStopPt]))    continue;

    //double val = f_y->Eval(McPhi.McY);
    //double random = gRandom->Uniform(0.0,1.0);
    //if(random > val) continue;
 

    //if( mMode == 3 && (McPhi.Centrality < 2 || McPhi.Centrality > 5) ) continue;
    if( mMode == 3) cout << "McPhi.McPt = " << McPhi.McPt << "    is in range [" << vmsa::pt_low[energy][mStartPt] << "," << vmsa::pt_up[energy][mStopPt] << "]" << endl;
  
    if( mMode == 1 && (McPhi.McPt <= vmsa::pt_low_cent[energy][mStartPt] || McPhi.McPt > vmsa::pt_up_cent[energy][mStopPt]) ) continue;
    //if( mMode == 1) cout << "McPhi.McPt = " << McPhi.McPt << "    is in range [" << vmsa::pt_low_cent[energy][mStartPt] << "," << vmsa::pt_up_cent[energy][mStopPt] << "]" << endl;

    if( mMode == 2 && (McPhi.McPt <= vmsa::pt_low_y[energy][mStartPt]    || McPhi.McPt > vmsa::pt_up_y[energy][mStopPt]   ) ) continue;
    //if( mMode == 2) cout << "McPhi.McPt = " << McPhi.McPt << "    is in range [" << vmsa::pt_low_y[energy][mStartPt] << "," << vmsa::pt_up_y[energy][mStopPt] << "]" << endl;
     

    ///if( mEtaMode == 0 && TMath::Abs(McPhi.McY) > 1.0 ) continue;
    if( mEtaMode == 0 && TMath::Abs(McPhi.McY) > 1.0 ) continue; //Swap to RC level cut
    //if( mEtaMode == 3 && TMath::Abs(McPhi.McY) > 0.4 ) continue;
    //if( mEtaMode == 4 && TMath::Abs(McPhi.McY) > 0.6 ) continue;
    //if( mEtaMode == 5 && TMath::Abs(McPhi.McY) > 0.8 ) continue;
   
    //cout << "BEFORE CUTS" << endl;
    //cout << "McKM.Phi = " << McKM.McPhi << "    RcKM.Phi = " << RcKM.RcPhi << endl;
    //cout << "McKP.Phi = " << McKP.McPhi << "    RcKP.Phi = " << RcKP.RcPhi << endl;

    //if( !mEffCut->passTrackCutPhi(McPhi) ) continue; // eta cuts for McPhi 
    //if( !mEffCut->passTrackCut(McKP) ) continue; // eta cut for McKplus
    //if( !mEffCut->passTrackCut(McKM) ) continue; // eta cut for McKminus

    //cout << "AFTER CUTS" << endl;
    //cout << "McKM.Phi = " << McKM.McPhi << "    RcKM.Phi = " << RcKM.RcPhi << endl;
    //cout << "McKP.Phi = " << McKP.McPhi << "    RcKP.Phi = " << RcKP.RcPhi << endl;


    double Psi = 0.0;

    if(mOrder == 1) 
    {
      Psi = McPhi.Psi1;
      //f_pDel1->FixParameter(0,mChi[0][int(McPhi.Centrality)]);
      //double delta = f_pDel1->GetRandom();
      //Psi += delta; // smear the reconstructed event plane
 
      ////Angle wrapping
      //while(Psi > TMath::Pi())
      //{
      //  Psi -= 2.0*TMath::Pi();
      //}
      //while(Psi < -TMath::Pi())
      //{
      //  Psi += 2.0*TMath::Pi();
      //}
    }

    if(mOrder == 2) 
    {
      Psi = McPhi.Psi2;
      //f_pDel2->FixParameter(0,mChi[1][int(McPhi.Centrality)]);
      //double delta = f_pDel2->GetRandom();
      //Psi += delta; // smear the reconstructed event plane
 
      ////Angle wrapping
      //while(Psi > 0.5*TMath::Pi())
      //{
      //  Psi -= TMath::Pi();
      //}
      //while(Psi < -0.5*TMath::Pi())
      //{
      //  Psi += TMath::Pi();
      //}
    }
 Gavin comment out 11/09/2024*/
    double deviationP = 0.0;//gRandom->Gaus(0,0.05*RcKP.RcPt);
    double deviationM = 0.0;//gRandom->Gaus(0,0.05*RcKM.RcPt);

    double posPt = RcKP.RcPt + deviationP;
    double negPt = RcKM.RcPt + deviationM;
 
    // Unsmeared 
    TLorentzVector lMcKP;
    lMcKP.SetPtEtaPhiM(RcKP.RcPt,RcKP.RcEta,RcKP.RcPhi,RcKP.RcM);
    TLorentzVector rl_kp = lMcKP;
    TLorentzVector lMcKM;
    lMcKM.SetPtEtaPhiM(RcKM.RcPt,RcKM.RcEta,RcKM.RcPhi,RcKM.RcM);
    TLorentzVector rl_km = lMcKM;

    TLorentzVector lMcPhi;
    lMcPhi = lMcKP + lMcKM;
    //lMcPhi.SetPtEtaPhiM(RcPhi.RcPt,RcPhi.RcEta,RcPhi.RcPhi,RcPhi.RcInvMass);
    TVector3 vMcPhiBeta = -1.0*lMcPhi.BoostVector();
    TVector3 phiMomentumLabUnit = lMcPhi.Vect().Unit();
 
    lMcKP.Boost(vMcPhiBeta);
    TVector3 vMcKP = lMcKP.Vect().Unit(); // direction of K+ momentum in phi-meson rest frame
 
    double costheta = vMcKP.Dot(phiMomentumLabUnit);

    TLorentzVector lBeamPos;
    lBeamPos.SetPxPyPzE(0.0,0.0,9.75,9.796);
    lBeamPos.Boost(vMcPhiBeta);
    TLorentzVector lBeamNeg;
    lBeamNeg.SetPxPyPzE(0.0,0.0,-9.75,9.796);
    lBeamNeg.Boost(vMcPhiBeta);
    
    TLorentzVector lBeamTot = lBeamPos + lBeamNeg;
  
    TVector3 vBeamPos = lBeamPos.Vect();
    TVector3 vBeamNeg = lBeamNeg.Vect();

    TVector3 vBeamTotalUnit = lBeamTot.Vect().Unit();
    TVector3 yaxisfrombeam = vBeamPos.Cross(vBeamNeg).Unit();

    TVector3 zaxisfrombeam = vBeamTotalUnit;
    TVector3 xaxisfrombeam = yaxisfrombeam.Cross(zaxisfrombeam).Unit();
    
    Double_t CosThetaStarH = vMcKP.Dot(zaxisfrombeam);

    double xprojection = vMcKP.Dot(xaxisfrombeam);
    double yprojection = vMcKP.Dot(yaxisfrombeam);
  
    double helicityangle = TMath::ATan2(yprojection,xprojection);
    while(helicityangle < 0.0) helicityangle += 2.0*TMath::Pi();
    while(helicityangle > 2.0*TMath::Pi()) helicityangle -= 2.0*TMath::Pi();

    TVector3 QVectorMc(TMath::Sin(McPhi.Psi),-1.0*TMath::Cos(McPhi.Psi),0.0);

    double phistar = vMcKP.Phi();
    TVector3 mcyprime(TMath::Cos(TMath::Pi()+McPhi.Psi),TMath::Sin(TMath::Pi()+McPhi.Psi),0.0);
    TVector3 mcxprime(0.0,0.0,1.0);
    Double_t mcproj_yprime = vMcKP.Dot(mcyprime);
    Double_t mcproj_xprime = vMcKP.Dot(mcxprime);
    Float_t mcphiprime = TMath::ATan2(mcproj_yprime,mcproj_xprime);
    while(mcphiprime <  0.0) mcphiprime += 2.0*TMath::Pi();
    while(mcphiprime > 2.0*TMath::Pi()) mcphiprime -= 2.0*TMath::Pi();

    TVector3 nQMc = QVectorMc.Unit(); // direction of QVector
    float McCosThetaStar = vMcKP.Dot(nQMc);

    // Smeared
    TLorentzVector slMcKP;
    slMcKP.SetPtEtaPhiM(posPt,RcKP.RcEta,RcKP.RcPhi,RcKP.RcM);
    TLorentzVector srl_kp = slMcKP;
    TLorentzVector slMcKM;
    slMcKM.SetPtEtaPhiM(negPt,RcKM.RcEta,RcKM.RcPhi,RcKM.RcM);
    TLorentzVector srl_km = slMcKM;

    TLorentzVector slMcPhi;
    slMcPhi = slMcKP + slMcKM;
    //lMcPhi.SetPtEtaPhiM(RcPhi.RcPt,RcPhi.RcEta,RcPhi.RcPhi,RcPhi.RcInvMass);
    TVector3 svMcPhiBeta = -1.0*slMcPhi.BoostVector();
    TVector3 sphiMomentumLabUnit = slMcPhi.Vect().Unit();
 
    slMcKP.Boost(svMcPhiBeta);
    TVector3 svMcKP = slMcKP.Vect().Unit(); // direction of K+ momentum in phi-meson rest frame
 
    double scostheta = svMcKP.Dot(sphiMomentumLabUnit);

    TLorentzVector slBeamPos;
    slBeamPos.SetPxPyPzE(0.0,0.0,9.75,9.796);
    slBeamPos.Boost(svMcPhiBeta);
    TLorentzVector slBeamNeg;
    slBeamNeg.SetPxPyPzE(0.0,0.0,-9.75,9.796);
    slBeamNeg.Boost(svMcPhiBeta);
    
    TLorentzVector slBeamTot = slBeamPos + slBeamNeg;
  
    TVector3 svBeamPos = slBeamPos.Vect();
    TVector3 svBeamNeg = slBeamNeg.Vect();

    TVector3 svBeamTotalUnit = slBeamTot.Vect().Unit();
    TVector3 syaxisfrombeam = svBeamPos.Cross(svBeamNeg).Unit();

    TVector3 szaxisfrombeam = svBeamTotalUnit;
    TVector3 sxaxisfrombeam = syaxisfrombeam.Cross(szaxisfrombeam).Unit();
    
    Double_t sCosThetaStarH = svMcKP.Dot(szaxisfrombeam);

    double sxprojection = svMcKP.Dot(sxaxisfrombeam);
    double syprojection = svMcKP.Dot(syaxisfrombeam);
  
    double shelicityangle = TMath::ATan2(syprojection,sxprojection);
    while(shelicityangle < 0.0) shelicityangle += 2.0*TMath::Pi();
    while(shelicityangle > 2.0*TMath::Pi()) shelicityangle -= 2.0*TMath::Pi();

    double sphistar = svMcKP.Phi();
    Double_t smcproj_yprime = svMcKP.Dot(mcyprime);
    Double_t smcproj_xprime = svMcKP.Dot(mcxprime);
    Float_t smcphiprime = TMath::ATan2(smcproj_yprime,smcproj_xprime);
    while(smcphiprime <  0.0) smcphiprime += 2.0*TMath::Pi();
    while(smcphiprime > 2.0*TMath::Pi()) smcphiprime -= 2.0*TMath::Pi();

    float sMcCosThetaStar = svMcKP.Dot(nQMc);
/*
    TLorentzVector lMcPhi;
    lMcPhi.SetPtEtaPhiM(McPhi.McPt,McPhi.McEta,McPhi.McPhi,vmsa::InvMass[pid]);
    TVector3 vMcPhiBeta = -1.0*lMcPhi.BoostVector();
    TVector3 phiMomentumLabUnit = lMcPhi.Vect().Unit();
 
    TLorentzVector lMcKP;
    lMcKP.SetPtEtaPhiM(McKP.McPt,McKP.McEta,McKP.McPhi,vmsa::mMassKaon);
    lMcKP.Boost(vMcPhiBeta);
    TVector3 vMcKP = lMcKP.Vect().Unit(); // direction of K+ momentum in phi-meson rest frame
 
    double costheta = vMcKP.Dot(phiMomentumLabUnit);

    TLorentzVector lBeamPos;
    lBeamPos.SetPxPyPzE(0.0,0.0,9.75,9.796);
    lBeamPos.Boost(vMcPhiBeta);
    TLorentzVector lBeamNeg;
    lBeamNeg.SetPxPyPzE(0.0,0.0,-9.75,9.796);
    lBeamNeg.Boost(vMcPhiBeta);
    
    TLorentzVector lBeamTot = lBeamPos + lBeamNeg;
  
    TVector3 vBeamPos = lBeamPos.Vect();
    TVector3 vBeamNeg = lBeamNeg.Vect();

    TVector3 vBeamTotalUnit = lBeamTot.Vect().Unit();
    TVector3 yaxisfrombeam = vBeamPos.Cross(vBeamNeg).Unit();

    TVector3 zaxisfrombeam = vBeamTotalUnit;
    TVector3 xaxisfrombeam = yaxisfrombeam.Cross(zaxisfrombeam).Unit();
    
    Double_t CosThetaStarH = vMcKP.Dot(zaxisfrombeam);

    double xprojection = vMcKP.Dot(xaxisfrombeam);
    double yprojection = vMcKP.Dot(yaxisfrombeam);
  
    double helicityangle = TMath::ATan2(yprojection,xprojection);
    while(helicityangle < 0.0) helicityangle += 2.0*TMath::Pi();
    while(helicityangle > 2.0*TMath::Pi()) helicityangle -= 2.0*TMath::Pi();

    TVector3 QVectorMc(TMath::Sin(McPhi.Psi),-1.0*TMath::Cos(McPhi.Psi),0.0);

    double phistar = vMcKP.Phi();
    TVector3 mcyprime(TMath::Cos(TMath::Pi()+McPhi.Psi),TMath::Sin(TMath::Pi()+McPhi.Psi),0.0);
    TVector3 mcxprime(0.0,0.0,1.0);
    Double_t mcproj_yprime = vMcKP.Dot(mcyprime);
    Double_t mcproj_xprime = vMcKP.Dot(mcxprime);
    Float_t mcphiprime = TMath::ATan2(mcproj_yprime,mcproj_xprime);
    while(mcphiprime <  0.0) mcphiprime += 2.0*TMath::Pi();
    while(mcphiprime > 2.0*TMath::Pi()) mcphiprime -= 2.0*TMath::Pi();

    TVector3 nQMc = QVectorMc.Unit(); // direction of QVector
    float McCosThetaStar = vMcKP.Dot(nQMc);
*/

/*
    TLorentzVector lRcPhi;
    lRcPhi.SetPtEtaPhiM(RcPhi.RcPt,RcPhi.RcEta,RcPhi.RcPhi,RcPhi.RcInvMass);
    TVector3 vRcPhiBeta = -1.0*lRcPhi.BoostVector();
    //TVector3 phiMomentumLabRc = lRcPhi.Vect();
    TVector3 phiMomentumLabUnitRc = lRcPhi.Vect().Unit();
    //TVector3 zAxisRc(0, 0, 1);
    //TVector3 rotationAxisRc = phiMomentumLabRc.Cross(zAxisRc);
    //double rotationAxisRc = phiMomentumLabRc.Angle(zAxisRc);

    //TRotation rotationRc; 
    //rotationRc.Rotate(rotationAngleRc, rotationAxisRc);

    TLorentzVector lRcKP;
    lRcKP.SetPtEtaPhiM(RcKP.RcPt,RcKP.RcEta,RcKP.RcPhi,vmsa::mMassKaon);
    lRcKP.Boost(vRcPhiBeta);
    TVector3 vRcKP = lRcKP.Vect().Unit(); // direction of K+ momentum in phi-meson rest frame

    //TLorentzVector KplusRc = lRcKP;
    //TVector3 KplusVecRc = lRcKP.Vect(); // direction of K+ momentum in phi-meson rest frame
    //KplusVecRc.Transform(rotationRc);
    //KplusRc.SetVect(KplusVecRc);
 
    //double thetaRc = KplusVecRc.Theta();
    //double costhetaRc = TMath::Cos(thetaRc);
    double costhetaRc = vRcKP.Dot(phiMomentumLabUnitRc);

    TVector3 QVectorRc(TMath::Sin(Psi),-1.0*TMath::Cos(Psi),0.0);
    double phistarRc = vRcKP.Phi();
    TVector3 xprime(TMath::Cos(Psi),TMath::Sin(Psi),0.0);
    TVector3 zprime(0.0,0.0,1.0);
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
*/

    //cout << "phistar - phistarRc = " << phistar - phistarRc << endl;

            
    //cout << "McPhi.McPhi = " << McPhi.McPhi << endl;
    //cout << "atan2(y,x) = " << TMath::ATan2(vMcKP.Y(),vMcKP.X()) << endl;    
    
    //if(mOrder == 3) 
    //{
    //  double randtheta = gRandom->Uniform(0.,TMath::Pi());
    //  double randPsi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
    //  //cout << randtheta << endl;
    //  //cout << randPsi << endl;
    //  QVector.SetXYZ(TMath::Sin(randPsi)*TMath::Cos(randtheta),-1.0*TMath::Cos(randPsi)*TMath::Cos(randtheta),TMath::Sin(randtheta));     
    //}
    
//    TVector3 nQRc = QVectorRc.Unit(); // direction of QVector
    //float McThetaStar = TMath::ACos(vMcKP.Dot(nQMc));
    //if(McThetaStar < 0.0) McThetaStar+=TMath::Pi();
    //if(McThetaStar > TMath::Pi()) McThetaStar-=TMath::Pi();
//    float RcCosThetaStar = vRcKP.Dot(nQRc);

/*Gavin comment out 11/09/24
    double helicityinglobal = -1.0*McCosThetaStar*TMath::Sqrt(1.-McPhi.McPt*McPhi.McPt/McPhi.McP/McPhi.McP*TMath::Sin(McPhi.McPhi-McPhi.Psi)*TMath::Sin(McPhi.McPhi-McPhi.Psi)) + TMath::Sqrt(1.-McCosThetaStar*McCosThetaStar)*TMath::Cos(McPhi.McPhi-McPhi.Psi)*McPhi.McPt/McPhi.McP*TMath::Sin(McPhi.McPhi-McPhi.Psi);
 
    //cout << "Helicity cos(theta*) = " << CosThetaStarH << ",    global frame relation = " << helicityinglobal << endl;

    //double globalinhelicity =  -1.0*CosThetaStarH*TMath::Sin(Psi)*lMcPhi.Pz()/lMcPhi.P() + CosThetaStarH*TMath::Cos(Psi)*lMcPhi.Py()/lMcPhi.P() + (vMcKP.Cross(zaxisfrombeam)).Dot(QVectorMc);
    //double globalinhelicity =  CosThetaStarH*zaxisfrombeam.Dot(QVectorMc) + (vMcKP-CosThetaStarH*zaxisfrombeam).Dot(QVectorMc);

    double th = TMath::ATan2(-lMcPhi.Pz(),lMcPhi.Perp());
    //double CosTh = TMath::Cos(th);
    //double SinTh = TMath::Sin(th);
    double CosTh = lMcPhi.Perp()/lMcPhi.P();
    double SinTh = -lMcPhi.Pz()/lMcPhi.P();
    double CosPsi = TMath::Cos(McPhi.Psi);
    double SinPsi = TMath::Sin(McPhi.Psi);
    double CosPhi = TMath::Cos(lMcPhi.Phi());
    double SinPhi = TMath::Sin(lMcPhi.Phi());
    double CosPhiPsi = TMath::Cos(lMcPhi.Phi()-McPhi.Psi);
    double SinPhiPsi = TMath::Sin(lMcPhi.Phi()-McPhi.Psi);
    //double globalinhelicity = -SinTh*TMath::Sqrt(1.-CosThetaStarH*CosThetaStarH)*TMath::Cos(lMcPhi.Phi()-helicityangle)-CosTh*CosThetaStarH;
    //double globalinhelicity = -TMath::Sin(lMcPhi.Phi()-McPhi.Psi)*SinTh*TMath::Sqrt(1.-CosThetaStarH*CosThetaStarH)*TMath::Cos(helicityangle)
    //                          -TMath::Cos(lMcPhi.Phi()-McPhi.Psi)*TMath::Sqrt(1.-CosThetaStarH*CosThetaStarH)*TMath::Sin(helicityangle)
    //                          +TMath::Sin(lMcPhi.Phi()-McPhi.Psi)*CosTh*CosThetaStarH;
    //double globalinhelicity = -SinTh*TMath::Sqrt(1.-CosThetaStarH*CosThetaStarH)*TMath::Cos(helicityangle)
    //                          -CosTh*TMath::Cos(lMcPhi.Phi()-McPhi.Psi)*TMath::Sqrt(1.-CosThetaStarH*CosThetaStarH)*TMath::Sin(helicityangle)
    //                          +TMath::Sin(lMcPhi.Phi()-McPhi.Psi)*CosTh*CosThetaStarH;

    double globalinhelicity =  SinTh*SinPhiPsi*TMath::Sqrt(1.-CosThetaStarH*CosThetaStarH)*TMath::Cos(helicityangle)
                              +CosPhiPsi*TMath::Sqrt(1.-CosThetaStarH*CosThetaStarH)*TMath::Sin(helicityangle)
                              +CosTh*SinPhiPsi*CosThetaStarH;
    TMatrixD Rot(3,3);
    //Rot(0, 0) = -SinPsi*CosTh*CosPhi+CosPsi*SinPhi;
    //Rot(0, 1) = -SinPsi*CosTh*SinPhi-CosPsi*CosPhi;
    //Rot(0, 2) = SinPsi*SinTh;
    //Rot(1, 0) = CosPsi*CosTh*CosPhi+SinPsi*SinPhi;
    //Rot(1, 1) = CosPsi*CosTh*SinPhi-SinPsi*CosPhi;
    //Rot(1, 2) = -CosPsi*SinTh;
    //Rot(2, 0) = SinTh*CosPhi;
    //Rot(2, 1) = SinTh*SinPhi;
    //Rot(2, 2) = CosTh;

    //Rot(0, 0) = -SinPsi*CosTh*CosPhi-CosPsi*SinPhi;
    //Rot(0, 1) = -SinPsi*CosTh*SinPhi+CosPsi*CosPhi;
    //Rot(0, 2) = SinPsi*SinTh;
    //Rot(1, 0) = CosPsi*CosTh*CosPhi-SinPsi*SinPhi;
    //Rot(1, 1) = CosPsi*CosTh*SinPhi+SinPsi*CosPhi;
    //Rot(1, 2) = -CosPsi*SinTh;
    //Rot(2, 0) = -SinTh*CosPhi;
    //Rot(2, 1) = -SinTh*SinPhi;
    //Rot(2, 2) = -CosTh;

    Rot(0, 0) = -CosPsi*SinPhi-CosPhi*CosTh*SinPsi;
    Rot(0, 1) =  CosPhi*CosPsi-CosTh*SinPhi*SinPsi;
    Rot(0, 2) =  SinPsi*SinTh;
    Rot(1, 0) =  CosPhi*CosPsi*CosTh-SinPhi*SinPsi;
    Rot(1, 1) =  CosPsi*CosTh*SinPhi+CosPhi*SinPsi;
    Rot(1, 2) = -CosPsi*SinTh;
    Rot(2, 0) = -SinTh*CosPhi;
    Rot(2, 1) = -SinTh*SinPhi;
    Rot(2, 2) = -CosTh;
  
   // Rot.Invert();

    TMatrixD RotH(3,3);
   // RotH(0, 0) = -CosTh;
   // RotH(0, 1) = 0.0;
   // RotH(0, 2) = -SinTh;
   // RotH(1, 0) = CosPhiPsi*SinTh;
   // RotH(1, 1) = -SinPhiPsi;
   // RotH(1, 2) = -CosPhiPsi*CosTh;
   // RotH(2, 0) = -SinPhiPsi*SinTh;
   // RotH(2, 1) = -CosPhiPsi;
   // RotH(2, 2) = SinPhiPsi*CosTh;

    RotH(0, 0) = -CosTh;
    RotH(0, 1) = 0.0;
    RotH(0, 2) = SinTh;
    RotH(1, 0) = CosPhiPsi*SinTh;
    RotH(1, 1) = -SinPhiPsi;
    RotH(1, 2) = CosPhiPsi*CosTh;
    RotH(2, 0) = SinPhiPsi*SinTh;
    RotH(2, 1) = CosPhiPsi;
    RotH(2, 2) = SinPhiPsi*CosTh;

    TMatrixD RotZ1(3,3);
    RotZ1(0, 0) = -CosPhi; 
    RotZ1(0, 1) = -SinPhi;
    RotZ1(0, 2) = 0.0;
    RotZ1(1, 0) = SinPhi;
    RotZ1(1, 1) = -CosPhi;
    RotZ1(1, 2) = 0.0;
    RotZ1(2, 0) = 0.0;
    RotZ1(2, 1) = 0.0;
    RotZ1(2, 2) = 1.0;

    TMatrixD RotZ3(3,3);
    RotZ3(0, 0) = SinPsi; 
    RotZ3(0, 1) = CosPsi;
    RotZ3(0, 2) = 0.0;
    RotZ3(1, 0) = -CosPsi;
    RotZ3(1, 1) = SinPsi;
    RotZ3(1, 2) = 0.0;
    RotZ3(2, 0) = 0.0;
    RotZ3(2, 1) = 0.0;
    RotZ3(2, 2) = 1.0;

    TMatrixD RotY2(3,3);
    RotY2(0, 0) = CosTh; 
    RotY2(0, 1) = 0.0;
    RotY2(0, 2) = SinTh;
    RotY2(1, 0) = 0.0;
    RotY2(1, 1) = 1.0;
    RotY2(1, 2) = 0.0;
    RotY2(2, 0) = -SinTh;
    RotY2(2, 1) = 0.0;
    RotY2(2, 2) = CosTh;

    TMatrixD totalRot = RotZ3 * RotY2 * RotZ1;
     
    TMatrixD vh(1,3);
    vh(0,0) = TMath::Sqrt(1.-CosThetaStarH*CosThetaStarH)*TMath::Cos(helicityangle);
    vh(0,1) = TMath::Sqrt(1.-CosThetaStarH*CosThetaStarH)*TMath::Sin(helicityangle);
    vh(0,2) = CosThetaStarH;

    TMatrixD zh(1,3);
    zh(0,0) = zaxisfrombeam.X();
    zh(0,1) = zaxisfrombeam.Y();
    zh(0,2) = zaxisfrombeam.Z();

    TMatrixD xh(1,3);
    xh(0,0) = xaxisfrombeam.X();
    xh(0,1) = xaxisfrombeam.Y();
    xh(0,2) = xaxisfrombeam.Z();

    TMatrixD yh(1,3);
    yh(0,0) = yaxisfrombeam.X();
    yh(0,1) = yaxisfrombeam.Y();
    yh(0,2) = yaxisfrombeam.Z();

    TVector3 zh2(zh(0,0),zh(0,1),zh(0,2));
    TVector3 yh2(yh(0,0),yh(0,1),yh(0,2));
    TVector3 xh2(xh(0,0),xh(0,1),xh(0,2));
    TVector3 vh2(vh(0,0),vh(0,1),vh(0,2));
    TVector3 vg2(TMath::Sqrt(1.-McCosThetaStar*McCosThetaStar)*TMath::Cos(mcphiprime),TMath::Sqrt(1.-McCosThetaStar*McCosThetaStar)*TMath::Sin(mcphiprime),McCosThetaStar);
    zh2.RotateZ(TMath::Pi()-lMcPhi.Phi());
    //zh2.RotateZ(-zaxisfrombeam.Phi());
    zh2.RotateY(th);
    //zh2.Print();
    zh2.RotateX(TMath::Pi());
    zh2.RotateZ(McPhi.Psi-TMath::Pi()/2.0);
    xh2.RotateZ(TMath::Pi()-lMcPhi.Phi());
    //zh2.RotateZ(-zaxisfrombeam.Phi());
    xh2.RotateY(th);
    //zh2.Print();
    xh2.RotateX(TMath::Pi());
    xh2.RotateZ(McPhi.Psi-TMath::Pi()/2.0);
    yh2.RotateZ(TMath::Pi()-lMcPhi.Phi());
    //zh2.RotateZ(-zaxisfrombeam.Phi());
    yh2.RotateY(th);
    //zh2.Print();
    yh2.RotateX(TMath::Pi());
    yh2.RotateZ(McPhi.Psi-TMath::Pi()/2.0);
    vh2.RotateZ(TMath::Pi()-lMcPhi.Phi());
    //zh2.RotateZ(-zaxisfrombeam.Phi());
    vh2.RotateY(th);
    //zh2.Print();
    vh2.RotateX(TMath::Pi());
    vh2.RotateZ(McPhi.Psi-TMath::Pi()/2.0);
    //cout << "Psi = " << McPhi.Psi << endl;
    //cout << "Psi-TMath::Pi()/2.0 = " << McPhi.Psi-TMath::Pi()/2.0 << endl;

    TMatrixD kp(1,3);
    kp(0,0) = vMcKP.X();
    kp(0,1) = vMcKP.Y();
    kp(0,2) = vMcKP.Z();

    zh.T();
    xh.T();
    yh.T();
    vh.T();
    kp.T();
    TMatrixD product = kp * (Rot * zh);

    cout << "Global cos(theta*) = " << McCosThetaStar << ", Helicity cos(theta*) = " << CosThetaStarH << ",  Helicity frame relation to get global cos(theta*) = " << globalinhelicity << endl;
    //cout << "Global cos(theta*) = " << McCosThetaStar << ", Helicity cos(theta*) = " << CosThetaStarH << ",  Helicity frame relation to get global cos(theta*) = " << endl;
    //product.Print();

    TMatrixD productRot = Rot * zh;
    TMatrixD productRotX = Rot * xh;
    TMatrixD productRotY = Rot * yh;
    TMatrixD productRotV = Rot * vh;
    TMatrixD productRotkp = Rot * kp; 
    
    TVector3 kp2(productRotkp(0,0),productRotkp(1,0),productRotkp(2,0));
    TVector3 kpdef(kp(0,0),kp(1,0),kp(2,0));
  
    //cout << "zg.Dot(zh) = " << nQMc.Dot(zaxisfrombeam) << endl;

    productRot.Print();
    zh2.Print();  
    QVectorMc.Print();
    productRotX.Print();
    xh2.Print();  
    productRotY.Print();
    yh2.Print();
    //QVectorMc.RotateZ(-TMath::Pi()/2.);
    //QVectorMc.Print();

    cout << "1. Rotation of helicity frame vector for K+ to the global frame" << endl;
    productRotV.Print();
    cout << "2. Before the rotation" << endl;
    vh.Print();
    cout << "3. Same thing as 1, but using Tvector3 rotations" << endl;
    vh2.Print();  
    cout << "4. Global kaon vector" << endl;
    vg2.Print();

    
    RotH.Print();
    //RotH.Invert();
    RotH.Print();
    TMatrixD productInvRotV = RotH * vh;
    productInvRotV.Print();
    

    cout << "kp2.Dot(QVectorMc) = " << kp2.Dot(QVectorMc) << endl;
    cout << "kp2.Dot(zh2) = " << kp2.Dot(zh2) << endl;
    cout << "kpdef.Dot(zh2) = " << kpdef.Dot(zh2) << endl;
    
    //cout << "Phi difference = " << zh2.Phi() - QVectorMc.Phi() << endl;
    //cout << endl;

    //cout << "helicity cos(theta*) from phi momentum = " << costheta << ", from boosts and sum = " << CosThetaStarH << endl; 
 
    //cout << "mc costheta* - rc costheta* = " << McCosThetaStar-RcCosThetaStar << endl;
Gavin comment out 11/09/24*/
    
//    int ptbin = -1; 
//    int centbin = -1;
//    int ybin = -1;
//    if(mMode == 0)
//    {
//      for(int ipt = 2; ipt < 6; ipt++)
//      {
//        //if(RcPhi.RcPt > vmsa::pt_low[energy][ipt] && RcPhi.RcPt <= vmsa::pt_up[energy][ipt]) 
//
//        if(McPhi.McPt > vmsa::pt_low[energy][ipt] && McPhi.McPt <= vmsa::pt_up[energy][ipt]) 
//        {
//          ptbin = ipt; 
//          //cout << "ptbin = " << ptbin << ", pt = " << McPhi.McPt << endl;
//          break;
//        }
//      }
//    }
    //if(mMode == 1)
    //{
    //  for(int ipt = 0; ipt < vmsa::pt_rebin_cent; ipt++)
    //  {
    //    if(McPhi.McPt > vmsa::pt_low_cent[energy][ipt] && McPhi.McPt <= vmsa::pt_up_cent[energy][ipt]) 
    //    {
    //      ptbin = ipt; 
    //      //cout << "ptbin = " << ptbin << ", pt = " << McPhi.McPt << endl;
    //      break;
    //    }
    //  }
    //  for(int icent = 0; icent < 9; icent++)
    //  {
    //    if(int(McPhi.Centrality) == icent) 
    //    {
    //      centbin = icent;
    //      break;
    //    }
    //  }
    //}
    //if(mMode == 2)
    //{
    //  for(int ipt = 0; ipt < vmsa::pt_rebin_y; ipt++)
    //  {
    //    if(McPhi.McPt > vmsa::pt_low_y[energy][ipt] && McPhi.McPt <= vmsa::pt_up_y[energy][ipt]) 
    //    {
    //      ptbin = ipt; 
    //      //cout << "ptbin = " << ptbin << ", pt = " << McPhi.McPt << endl;
    //      break;
    //    }
    //  }
    //  for(int icent = 0; icent < vmsa::cent_rebin_total; icent++)
    //  {
    //    if(int(McPhi.Centrality) >= vmsa::cent_rebin[icent] && int(McPhi.Centrality) < vmsa::cent_rebin[icent+1]) 
    //    {
    //      centbin = icent;
    //     // cout << "cent = " << McPhi.Centrality << ",  centbin = " << centbin << endl;
    //      break;
    //    }
    //  }
    //  for(int iy = 0; iy < 10; iy++)
    //  {
    //    if(McPhi.McY >= float(iy-5)/5. && McPhi.McY < float(iy-4)/5.)
    //    {
    //      ybin = iy;
    //      //cout << "y = " << McPhi.McY << ", ybin = " << ybin << endl;
    //      break;
    //    }
    //    else if ( McPhi.McY ==  1.)
    //    {
    //      ybin = 9; 
    //    }
    //  }
    //}

    ////cout << "Before Sampling" << endl;
    ////cout << "ptbin = " << ptbin << ", centbin = " << centbin << ", ybin = " << ybin << endl;
    //cout << "Before Sampling" << endl;
    
    //if(mMode == 0)
   // {
      //if(!Sampling(f_mRhoPt[ptbin],TMath::Abs(McCosThetaStar),mMaxData)) continue;
      //if(!Sampling(f_mRhoPt[ptbin],McCosThetaStar,mMaxData[ptbin])) continue;
      //cout << "Are we sampling? " << endl;
      //cout << "McCosThetaStar = " << McCosThetaStar << " mcphiprime = " << mcphiprime << endl;
      //if(!Sampling2D(f_mRhoPt_2D[ptbin],McCosThetaStar,mcphiprime,mMax)) continue;
      //if(!SamplingHelicity(f_mRhoPt_Helicity[ptbin],costheta,mMaxHelicity[ptbin])) continue;
      //cout << "pass sampling" << endl;
   // }
    //if(mMode == 1)
    //{
    //  if(!Sampling(f_mRhoCent[ptbin][centbin],TMath::Abs(McCosThetaStar))) continue;
    //}
    //if(mMode == 2)
    //{
    //  if(!Sampling(f_mRhoY[ptbin][centbin][ybin],TMath::Abs(McCosThetaStar))) continue;
    //}
    //double spectrarandom = gRandom->Uniform(0.0,1.0);
    //double spectravalue  = h_mSpectraRatio->GetBinContent(h_mSpectraRatio->FindBin(McPhi.McPt,McPhi.McY));
    //if(spectrarandom > spectravalue) continue;  
    
    //cout << "After Sampling" << endl;

    //mEffHistManger->FillPhiHistMc(McPhi.Centrality,McPhi.McPt,McPhi.McPhi,McPhi.McY,phistar,McCosThetaStar,costheta,McKP.McPt,McKP.McY);
    //cout << "After MC Fill 1" << endl;

    //flow loop 

    double phiPsi = RcPhi.RcPhi - McPhi.Psi;
    //cout << "phiPsi = " << phiPsi/TMath::Pi() << endl;
    while(phiPsi < -TMath::Pi()) phiPsi += 2.0*TMath::Pi();
    while(phiPsi >  TMath::Pi()) phiPsi -= 2.0*TMath::Pi();

    double valRapidity = -1;
    const double  pythialow[19] = { 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.6, 4.2};;
    const double pythiahigh[19] = { 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.6, 4.2, 5.0};;
    
    //cout << "|y| = " << TMath::Abs(lMcPhi.Rapidity()) << endl;
    //cout << "pT  = " << lMcPhi.Pt() << endl;
    //cout << "ptbinedges[mInputPt] = " << ptbinedges[mInputPt] << endl;
    //cout << "ptbinedges[mInputPt+1] = " << ptbinedges[mInputPt+1] << endl;

    if(TMath::Abs(lMcPhi.Rapidity()) > 1.0) continue; 
    if(!(lMcPhi.Pt() >= ptbinedges[mInputPt] && lMcPhi.Pt() < ptbinedges[mInputPt+1])) continue; 


    double weighty = 1.0;
    //for(int i = 0; i < 19; i++)

    for(int i = yptstart[mInputPt]; i <= yptstop[mInputPt]; i++)
    {
      if(lMcPhi.Pt() >= pythialow[i] && lMcPhi.Pt() < pythiahigh[i]) 
      {
        //valRapidity = pythiaflat[i]->Eval(slMcPhi.Rapidity())/pythiaflat[i]->GetMaximum();
        weighty = pythiaflat[i]->Eval(lMcPhi.Rapidity());
        break;
      } 
    }  
    //if(valRapidity < 0) continue;

    //double randRapidity = randomGen.Uniform(0,1);
    //if(valRapidity > randRapidity) continue;

    //double valV2 = f_flowreal->Eval(phiPsi)/f_flowreal->GetMaximum();
    //double randV2 = randomGen.Uniform(0,1);
    //if(valV2 > randV2) continue;

    totalbeforemc++;

    float etacuts[10] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    for(int i = 0; i < 1; i++)
    {
      double v2real = f_mV2[int(McPhi.Centrality)]->Eval(lMcPhi.Pt());
      f_flowreal->SetParameter(0,v2real);
      //if(i == 0) f_flowreal->SetParameter(0,0.0);
      //if(i == 1) f_flowreal->SetParameter(0,v2real);
      double weightflow = weightflow = f_flowreal->Eval(phiPsi);
      //double weightflow = f_flowreal->Eval(phiPsi);
      //helicity spin alignment loop
      for(int irho = 0; irho < 5; irho++)
      {
        for(int j = 0; j < 5; j++)
        {
          //double weight = f_mRhoPt_2D[irho][j]->Eval(CosThetaStarH,helicityangle) * f_flow[i]->Eval(phiPsi);
          //double weight = f_mRhoPt_2D[irho][j]->Eval(McCosThetaStar,mcphiprime);// * f_flow[i]->Eval(phiPsi);
          //double weight = f_mRhoPt_2D[irho][j]->Eval(McCosThetaStar,mcphiprime) * f_flowreal->Eval(phiPsi);
          //double valRho = f_mRhoPt_2D[irho][j]->Eval(McCosThetaStar,mcphiprime)/f_mRhoPt_2D[irho][j]->GetMaximum();
  
          //double valRho = f_mRhoPt_2D[irho][j]->Eval(McCosThetaStar,mcphiprime)/maxRho[irho][j];
        
          double weightrho = f_mRhoPt_2D[irho][j]->Eval(McCosThetaStar,mcphiprime);

          double weight = weighty * weightflow * weightrho;
          //double randRho = randomGen.Uniform(0,1);

          //if(valRho > randRho) continue;

          //mEffHistManger->FillHistMc(McPhi.Centrality,McPhi.McPt,McPhi.McY,McPhi.McPhi,McCosThetaStar,McPhi.Psi,phistar,mcphiprime,McKP.McPt,McKP.McY,CosThetaStarH,helicityangle,weight,i,irho,j);
          for(int ipt = 0; ipt < 5; ipt+=4)
          {
            //double deviationP = gRandom->Gaus(0,0.05*McKP.McPt);
            //double deviationM = gRandom->Gaus(0,0.05*McKM.McPt);
            if(RcKP.RcPt < float(ipt)*0.025 || RcKM.RcPt < float(ipt)*0.025) continue;
            //if(posPt < float(ipt)*0.025 || negPt < float(ipt)*0.025) continue;
            //mEffHistManger->FillHistMc(McPhi.Centrality,slMcPhi.Pt(),slMcPhi.Rapidity(),slMcPhi.Phi(),sMcCosThetaStar,McPhi.Psi,sphistar,smcphiprime,McKP.McPt,McKP.McY,sCosThetaStarH,shelicityangle,weight,i,irho,j,0,ipt);
            mEffHistManger->FillHistMc(McPhi.Centrality,lMcPhi.Pt(),lMcPhi.Rapidity(),lMcPhi.Phi(),McCosThetaStar,McPhi.Psi,phistar,mcphiprime,McKP.McPt,McKP.McY,CosThetaStarH,helicityangle,weight,i,irho,j,0,ipt);
            //mEffHistManger->FillHistMc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,sMcCosThetaStar,McPhi.Psi,sphistar,smcphiprime,McKP.McPt,McKP.McY,sCosThetaStarH,shelicityangle,weight,i,irho,j,0,ipt);
          }

          if(TMath::Abs(slMcPhi.Rapidity()) > 1.0) continue; 
          if(!(slMcPhi.Pt() >= ptbinedges[mInputPt] && slMcPhi.Pt() < ptbinedges[mInputPt+1])) continue; 
              
          for(int ieta = 9; ieta < 10; ieta++)
          {
            if(fabs(RcKP.RcEta) > etacuts[ieta] || fabs(RcKM.RcEta) > etacuts[ieta]) continue;
            for(int ipt = 0; ipt < 5; ipt+=4)
            {
              //double deviationP = gRandom->Gaus(0,0.05*McKP.McPt);
              //double deviationM = gRandom->Gaus(0,0.05*McKM.McPt);
              //if(RcKP.RcPt < float(ipt)*0.025 || RcKM.RcPt < float(ipt)*0.025) continue;
              if(posPt < float(ipt)*0.025 || negPt < float(ipt)*0.025) continue;
              //cout << "Before TPC Reconstructed" << endl;

              if(!tpcReconstructed(0,McPhi.Centrality,rl_kp)) continue;
              if(!tpcReconstructed(1,McPhi.Centrality,rl_km)) continue;
              //cout << "After TPC Reconstructed" << endl;

              bool passToF = false;

              int GlobalBinP = -1;
              int GlobalBinM = -1;
              GlobalBinP = ToFHist[0]->FindBin(srl_kp.Pt(),srl_kp.Eta(),srl_kp.Phi());
              GlobalBinM = ToFHist[1]->FindBin(srl_km.Pt(),srl_km.Eta(),srl_km.Phi());
              if(GlobalBinP < 0 || GlobalBinM < 0) passToF = false;   
              float valToFP = ToFHist[0]->GetBinContent(GlobalBinP);
              float valToFM = ToFHist[1]->GetBinContent(GlobalBinM);
      
              float probP = randomGen.Uniform(0,1);
              float probM = randomGen.Uniform(0,1);

              if(probP < valToFP && probM < valToFM) passToF = true;


              if(!passToF) continue;             
              //cout << "After TOF" << endl;
          
              Bool_t passnskp = true;
              Bool_t passnskm = true;
              GlobalBinM = -1;
              GlobalBinP = -1;
              findnsigHistPID(srl_km,0,GlobalBinM);
              findnsigHistPID(srl_kp,1,GlobalBinP);
              if(GlobalBinM < 0 || GlobalBinP < 0) 
              {
                passnskp = false;
                passnskm = false;
              }
              else
              {
                if( !passhist_nsig_PID(0,9,srl_km,GlobalBinM) ) passnskm = false;
                if( !passhist_nsig_PID(1,9,srl_kp,GlobalBinP) ) passnskp = false;
              }

              if(!passnskm || !passnskp) continue;
              //cout << "After nsig" << endl;


              Bool_t passm2kp = true;
              Bool_t passm2km = true;
              GlobalBinM = -1;
              GlobalBinP = -1;
              findm2HistPID(srl_km,0,GlobalBinM);
              findm2HistPID(srl_kp,1,GlobalBinP);
              if(GlobalBinM < 0 || GlobalBinP < 0) 
              {
                passm2kp = false;
                passm2km = false;
              }
              else
              {
                if( !passhist_m2_PID(0,9,srl_km,GlobalBinM) ) passm2km = false;
                if( !passhist_m2_PID(1,9,srl_kp,GlobalBinP) ) passm2kp = false; 
              }

              if(!passm2km || !passm2kp) continue;
              //cout << "After nsig" << endl;


              //mEffHistManger->FillHistMc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,McCosThetaStar,McPhi.Psi,phistar,mcphiprime,McKP.McPt,McKP.McY,CosThetaStarH,helicityangle,weight,i,irho,j,ieta+1,ipt);
              //if(eventcounts[irho][j] < eventtarget[mInputPt-2]) mEffHistManger->FillHistMc(McPhi.Centrality,slMcPhi.Pt(),slMcPhi.Rapidity(),slMcPhi.Phi(),sMcCosThetaStar,McPhi.Psi,sphistar,smcphiprime,McKP.McPt,McKP.McY,sCosThetaStarH,shelicityangle,weight,i,irho,j,ieta+1,ipt);
              //if(eventcounts[irho][j] < eventtarget[mInputPt]) mEffHistManger->FillHistMc(McPhi.Centrality,slMcPhi.Pt(),slMcPhi.Rapidity(),slMcPhi.Phi(),sMcCosThetaStar,McPhi.Psi,sphistar,smcphiprime,McKP.McPt,McKP.McY,sCosThetaStarH,shelicityangle,1.0,i,irho,j,ieta+1,ipt);
              //mEffHistManger->FillHistMc(McPhi.Centrality,lMcPhi.Pt(),lMcPhi.Rapidity(),lMcPhi.Phi(),McCosThetaStar,McPhi.Psi,phistar,mcphiprime,McKP.McPt,McKP.McY,CosThetaStarH,helicityangle,weight,i,irho,j,ieta+1,ipt);
              mEffHistManger->FillHistMc(McPhi.Centrality,slMcPhi.Pt(),slMcPhi.Rapidity(),slMcPhi.Phi(),sMcCosThetaStar,McPhi.Psi,sphistar,smcphiprime,McKP.McPt,McKP.McY,sCosThetaStarH,shelicityangle,weight,i,irho,j,ieta+1,ipt);
               
              //eventcounts[irho][j]++;
            } 
          }
        }
      }
    }
    //cout << "After MC Fill 2" << endl;
    //mEffHistManger->FillHistMc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc);
    //mEffHistManger->FillKaonHistMc(McPhi.Centrality,McPhi.McPhi,McPhi.McPt,McPhi.McY,McKP.McPt,McKP.McY,McKP.McEta,McKP.McPhi,phistar,McCosThetaStar,costheta); // K+
    //cout << "After MC Fill 3" << endl;
    //mEffHistManger->FillKaonHistMc(McPhi.Centrality,McPhi.McPhi,McPhi.McPt,McPhi.McY,McKM.McPt,McKM.McY,McKM.McEta,McKM.McPhi,phistar,McCosThetaStar,costheta); // K-
    //cout << "After MC Fill 4" << endl;
    //mEffHistManger->FillKaonDeltaHistMc(McPhi.Centrality,McPhi.McPhi,McPhi.McEta,McPhi.McPt,McPhi.McY,McKP.McPt,McKP.McEta,McKP.McPhi,McKM.McPt,McKM.McEta,McKM.McPhi,phistar); // K-
    //cout << "After MC Fill 5" << endl;
/*
    if(mEtaMode == 0)
    {
      //if( fabs(McKP.McEta) > 0.5 || fabs(McKM.McEta) > 0.5 ) continue; // eta cuts for McPhi 
      if( fabs(McKP.McEta) > 1.0 || fabs(McKM.McEta) > 1.0 ) continue; // eta cuts for McPhi 
    }
    else if(mEtaMode == 1)
    {
      if(    ! (fabs(McKP.McEta) <= 1.0 && (fabs(McKM.McEta) < 1.5 && fabs(McKM.McEta) > 1.0) ) 
          && ! (fabs(McKM.McEta) <= 1.0 && (fabs(McKP.McEta) < 1.5 && fabs(McKP.McEta) > 1.0) ) ) continue; // eta cuts for McPhi 
    }
    else if(mEtaMode == 2)
    {
      if( ! (fabs(McKM.McEta) < 1.5 && fabs(McKM.McEta) > 1.0 && fabs(McKP.McEta) < 1.5 && fabs(McKP.McEta) > 1.0 ) ) continue; // eta cuts for McPhi 
    }
    else if(mEtaMode == 3)
    {
      if( fabs(McKP.McEta) > 0.4 || fabs(McKM.McEta) > 0.4 ) continue; // eta cuts for McPhi 
    }
    else if(mEtaMode == 4)
    {
      if( fabs(McKP.McEta) > 0.6 || fabs(McKM.McEta) > 0.6 ) continue; // eta cuts for McPhi 
    }
    else if(mEtaMode == 5)
    {
      if( fabs(McKP.McEta) > 0.8 || fabs(McKM.McEta) > 0.8 ) continue; // eta cuts for McPhi 
    }



    //cout << "Right Before RC Fill " << endl;
    if( !mEffCut->passTrackCut(RcKP) ) continue; // eta and TPC cuts for RcKplus
    if( !mEffCut->passTrackCut(RcKM) ) continue; // eta and TPC cuts for RcKminus
    if( !mEffCut->passTrackCutPhi(RcPhi) ) continue;  // eta cuts for RcPhi 
    mEffHistManger->FillPhiHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcPhi,RcPhi.RcY,phistarRc,RcCosThetaStar,costhetaRc,RcKP.RcPt,RcKP.RcY);
    mEffHistManger->FillHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc,RcKP.RcPt,RcKP.RcY,costhetaRc);
    //mEffHistManger->FillKaonHistRc(McPhi.Centrality,RcPhi.RcPhi,RcPhi.RcPt,RcPhi.RcY,RcKP.RcPt,RcKP.RcY,RcKP.RcEta,RcKP.RcPhi,phistarRc,RcCosThetaStar,costhetaRc); // K+
    mEffHistManger->FillKaonHistRc(McPhi.Centrality,RcPhi.RcPhi,RcPhi.RcPt,RcPhi.RcY,RcKM.RcPt,RcKM.RcY,RcKM.RcEta,RcKM.RcPhi,phistarRc,RcCosThetaStar,costhetaRc); // K-
    mEffHistManger->FillKaonDeltaHistRc(McPhi.Centrality,RcPhi.RcPhi,RcPhi.RcEta,RcPhi.RcPt,RcPhi.RcY,RcKP.RcPt,RcKP.RcEta,RcKP.RcPhi,RcKM.RcPt,RcKM.RcEta,RcKM.RcPhi,phistarRc); // K-

    //cout << "After RC Fill " << endl;
    //cout << "Is RcPhi.RcPt = McPhi.McPt?           "; if(RcPhi.RcPt == McPhi.McPt) cout << "YES" << endl; else cout << "NO" << endl;
    //cout << "RcPhi.RcPt = " << RcPhi.RcPt << "   McPhi.McPt = " << McPhi.McPt << endl;
    //cout << "Is RcPhi.RcP = McPhi.McP?             "; if(RcPhi.RcP == McPhi.McP) cout << "YES" << endl; else cout << "NO" << endl;
    //cout << "RcPhi.RcP = " << RcPhi.RcP << "   McPhi.McP = " << McPhi.McP << endl;
    //cout << "Is RcPhi.RcEta = McPhi.McEta?         "; if(RcPhi.RcEta == McPhi.McEta) cout << "YES" << endl; else cout << "NO" << endl;
    //cout << "RcPhi.RcEta = " << RcPhi.RcEta << "   McPhi.McEta = " << McPhi.McEta << endl;
    //cout << "Is RcPhi.RcPhi = McPhi.McPhi?         "; if(RcPhi.RcPhi == McPhi.McPhi) cout << "YES" << endl; else cout << "NO" << endl;
    //cout << "RcPhi.RcPhi = " << RcPhi.RcPhi << "   McPhi.McPhi = " << McPhi.McPhi << endl;
    //cout << "Is RcPhi.RcInvMass = McPhi.McInvMass? "; if(RcPhi.RcInvMass == McPhi.McInvMass) cout << "YES" << endl; else cout << "NO" << endl;
    //cout << "RcPhi.RcInvMass = " << RcPhi.RcInvMass << "   McPhi.McInvMass = " << McPhi.McInvMass << endl;
  */
  }
  
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  cout << "Total Processed = " << totalprocessed << endl;
  cout << "Total Before MC = " << totalbeforemc << endl;

  //mEffHistManger->CalEffCosThetaStar();
  //cout << "calculated efficiency" << endl;
}

void StEffMcPhiHelicityGlobal::Finish()
{
  mFile_OutPut->cd();
  cout << "before writing hist" << endl;
  mEffHistManger->WriteHist();
  //mEffHistManger->WriteKaonHist();
  //mEffHistManger->WritePhiHist();
  cout << "after writing hist" << endl;
  mFile_OutPut->Close();
}

bool Sampling(TF1 *f_rhoPhy, float CosThetaStar, float wMax)
{
  return !(gRandom->Rndm() > f_rhoPhy->Eval(CosThetaStar)/wMax);
}
bool SamplingHelicity(TF1 *f_rhoPhy, float CosThetaStar, float wMax)
{
  return !(gRandom->Rndm() > f_rhoPhy->Eval(CosThetaStar)/wMax);
}

bool Sampling2D(TF2 *f_rhoPhy, float ThetaStar, float PhiPrime, float wMax)
{
  return !(gRandom->Rndm() > f_rhoPhy->Eval(ThetaStar,PhiPrime)/wMax);
}
bool Sampling2DWeight(TF2 *f_rhoPhy, float ThetaStar, float PhiPrime, float wMax)
{
  return f_rhoPhy->Eval(ThetaStar,PhiPrime);
}

