#include "StEffMcPhi.h"
#include "StEffHistManger.h"
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

ClassImp(StEffMcPhi)

int StEffMcPhi::mInput_flag = 1;

bool Sampling(TF1 *f_rhoPhy, float CosThetaStar, float max);
bool SamplingHelicity(TF1 *f_rhoPhy, float CosThetaStar, float max);
bool Sampling2D(TF2 *f_rhoPhy, float ThetaStar, float beta, float max);
bool Sampling2D(TF2 *f_rhoPhy, float ThetaStar, float beta, float max);

void StEffMcPhi::findFuncPID(TLorentzVector const& lKaon, int icharge, int &EtaBin, int &PhiBin)
{
  double eta = lKaon.PseudoRapidity();
  double phi = lKaon.Phi();
  while(phi < -TMath::Pi()) phi += 2.0*TMath::Pi();
  while(phi >  TMath::Pi()) phi -= 2.0*TMath::Pi();

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

bool StEffMcPhi::pass_m2_PID(int icharge, int icent, TLorentzVector const& lKaon, int EtaBin, int PhiBin)
{
  TF1 *f_m2 = NULL;

  string KEY; 
  KEY = Form("m2_PID_cent%d_charge%d_eta%d_phi%d",icent,icharge,EtaBin,PhiBin);
  f_m2 = f_m2_PID[KEY]; // only 20-60%

  double pt = lKaon.Perp();

  bool pass_m2 = false;
  double prob_m2 = 0.0;
  prob_m2 = f_m2->Eval(pt);
  pass_m2 = gRandom->Rndm() < prob_m2;
   
  return pass_m2;
}

bool StEffMcPhi::pass_nsig_PID(int icharge, int icent, TLorentzVector const& lKaon, int EtaBin, int PhiBin)
{
  TF1 *f_nsig = NULL;

  string KEY; 
  KEY = Form("nsig_PID_cent%d_charge%d_eta%d_phi%d",icent,icharge,EtaBin,PhiBin);
  f_nsig = f_nsig_PID[KEY]; // only 20-60%

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

void StEffMcPhi::readm2PID()
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

void StEffMcPhi::readnsigPID()
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

//TF1* readm2PID(int icent, int icharge, float ieta, float iphi)
//{
//  TH1F *parameterHist; 
//  string funcname = Form("m2_PID_cent%d_charge%d_eta%d_phi%d",icent,icharge,ieta,iphi);
//  TF1 *m2_PID_pT = new TF1(funcname.c_str(),m2_PID,0.1,5.0,5); 
//
//  string InPut = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/v2/m2_PID_funcparameters_19GeV.root");
//  TFile *File = TFile::Open(InPut.c_str());
//  std::cout << "m2 PID file: " << InPut << endl; 
//
//  string histname = Form("parameters_cent%d_charge%d_eta%d_phi%d",icent,icharge,ieta,iphi);
//  parameterHist = (TH1F*)File->Get(histname.c_str());
//
//  for(int i_par = 0; i_par < 5; i_par++)
//  {
//    m2_PID_pT->SetParameter(parameterHist->GetBinContent(1+i_par));
//  }
// 
//  return m2_PID_pT;
//}
void StEffMcPhi::readEfficiency(int energy)
{
  string inputPhi = Form("/star/data01/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/Efficiency/TPC/Eff_Phi_%s_Rc5Rc2.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Phi = TFile::Open(inputPhi.c_str());
  cout << "OPEN Efficiency File for phi-meson: " << inputPhi.c_str() << endl;

  h_FrameEta[0] = (TH1F*)File_Phi->Get("h_FrameEta");
  h_FramePhi[0] = (TH1F*)File_Phi->Get("h_FramePhi");

  cout << "Loaded the frames" << endl; 

  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_eta = 0; i_eta < h_FrameEta[0]->GetNbinsX()/*vmsa::BinEta*/; ++i_eta)
    {
      for(int i_phi = 0; i_phi < h_FramePhi[0]->GetNbinsX()/*vmsa::BinPhi*/; ++i_phi)
      {
	string KEY = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_EffPhi[KEY] = (TH1F*)File_Phi->Get(KEY.c_str());
        cout << "Loaded: " << KEY << endl; 
      }	
    }
  }
}

void StEffMcPhi::findHist(TLorentzVector* lPhi, int iParticleIndex, int& EtaBin, int& PhiBin)
{
  float eta = lPhi->Eta();
  EtaBin = h_FrameEta[iParticleIndex]->FindBin(eta)-1;

  float phi = lPhi->Phi();
  while(phi < 0.0) phi += TMath::Pi()*2.0;
  while(phi > TMath::Pi()*2.0) phi -= TMath::Pi()*2.0;
  PhiBin = h_FramePhi[iParticleIndex]->FindBin(phi)-1;
}

bool StEffMcPhi::tpcReconstructed(int iParticleIndex, int cent, TLorentzVector* lPhi)
{
   TH1F *h_TPC = NULL;
   int EtaBin_TPC = -1;
   int PhiBin_TPC = -1;
   findHist(lPhi,iParticleIndex,EtaBin_TPC,PhiBin_TPC);

   //cout << "EtaBin = " << EtaBin_TPC << ", PhiBin = " << PhiBin_TPC << endl;

   string KEY_TPC = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent,EtaBin_TPC,PhiBin_TPC); // get TPC eff
   h_TPC = h_EffPhi[KEY_TPC];
   
   double pt = lPhi->Perp();
   int const bin_TPC = h_TPC->FindBin(pt);
   bool is_TPC = gRandom->Rndm() < h_TPC->GetBinContent(bin_TPC);

   return is_TPC;
}

StEffMcPhi::StEffMcPhi(int Energy, long StartEvent, long StopEvent, int PID, int year, int mode, int inputpt, int startpt, int stoppt, const char* setting, int etamode, int order = 4) 
{
  mOrder = order;
  energy = Energy;
  pid = PID;
  mInputPt = inputpt;
  mStartPt = startpt;
  mStopPt = stoppt;
  mEtaMode = etamode;
  mMode = mode;

  for(int i = 0; i < 9; i++)
  {
    f_mV2[i] = readv2(energy,pid,i);
  }

  readm2PID();
  readnsigPID();
  readEfficiency(Energy);  

  std::string EP[2] = {"","2nd"};

  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_prelimv2_EPeff_randomRP/%s_%s_yabs1_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_PhiEffFinerBins_prelimv2_flatRP_yabs1p1/%s_%s_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_nov2_flatRP_yabs1_PhiEmbed/%s_%s_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_prelimv2_flatRP_yabs1_PhiEmbed/%s_%s_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_nov2_nopTspectra_flatRP_yabs1_PhiEmbed/%s_%s_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_PhiEffFinerBins_nov2_flatRP_yabs1p1_flatpt/%s_%s_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_nov2_flatRP_yabs1p1/%s_%s_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_PhiEff/%s_%s_yabs1_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_v2times3/%s_%s_yabs1_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_EP/%s_%s_eta1_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);

  SetInPutFile(InPutFile); // set input list

  SetStartEvent(StartEvent); // set start event
  SetStopEvent(StopEvent); // set stop event

  string OutPutFile = Form("Eff_%s_SingleParticle_%s_Mode%d_EtaMode%d_PtBins%d_%d.root",vmsa::mBeamEnergy[energy].c_str(),setting,mode,etamode,startpt,stoppt);
  SetOutPutFile(OutPutFile); // set output file

  mEffCut = new StEffCut();
  mEffHistManger = new StEffHistManger(energy, pid, mode, startpt, stoppt);

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

StEffMcPhi::~StEffMcPhi()
{
}

//------------------------------------------------------------
void StEffMcPhi::SetInPutFile(const string inputfile)
{
  mInPutFile = inputfile;
  cout << "Input file was set to: " << mInPutFile.c_str() << endl;
}

void StEffMcPhi::SetOutPutFile(const string outputfile)
{
  mOutPutFile = outputfile;
  cout << "Output file was set to: " << mOutPutFile.c_str() << endl;
}

void StEffMcPhi::SetStartEvent(const long StartEvent)
{
  mStartEvent = StartEvent;
  cout << "nStartEvent = " << mStartEvent << endl;
}

void StEffMcPhi::SetStopEvent(const long StopEvent)
{
  mStopEvent = StopEvent;
  cout << "nStopEvent = " << mStopEvent << endl;
}
//------------------------------------------------------------

void StEffMcPhi::Init()
{
  mEffHistManger->InitHist();
  //mEffHistManger->InitKaonHist();
  //mEffHistManger->InitPhiHist();


  // initialize the TNtuple
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
  mNtuple->SetBranchAddress("KmRcTof",&mKmRcTof);

  mNtuple->SetBranchAddress("RcPt",&mRcPt);
  mNtuple->SetBranchAddress("RcP",&mRcP);
  mNtuple->SetBranchAddress("RcEta",&mRcEta);
  mNtuple->SetBranchAddress("RcY",&mRcY);
  mNtuple->SetBranchAddress("RcPhi",&mRcPhi);
  mNtuple->SetBranchAddress("RcInvMass",&mRcInvMass);
  mNtuple->SetBranchAddress("RcTpc",&mRcTpc);

  int num_tracks = mNtuple->GetEntriesFast();
  cout << "Number of tracks in McPhiMeson = " << num_tracks<< endl;

  if(mStartEvent > num_tracks) mStartEvent = num_tracks;
  if(mStopEvent  > num_tracks) mStopEvent  = num_tracks;
  cout << "New nStartEvent = " << mStartEvent << ", new nStopEvent = " << mStopEvent << endl;

  mFile_OutPut = new TFile(mOutPutFile.c_str(),"RECREATE");
  f_y = new TF1("f_y",Form("exp(-x*x/2/%s)",mSigmay.c_str()),-1.0,1.0);

  // spin alignment information
  //std::string inputfilept     = Form("StRoot/Utility/Rho/%s/Rho_AccResSysErrors_F_0_eta1_eta1_PolySys_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
  std::string inputfilept     = Form("StRoot/Utility/Rho/%s/AccResPhiPtSys_eta1_eta1_Poly_WithRapiditySpectra_HalfSigma_FoldCos.root",vmsa::mBeamEnergy[energy].c_str());
  if(mOrder == 1) inputfilept = Form("StRoot/Utility/Rho/%s/Rho_AccResSysErrors_F_0_eta1_eta1_PolySys_FirstOrder_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
  std::string inputfilecent;
  std::string inputfiley;
  if(mOrder == 1) inputfilecent = Form("StRoot/Utility/Rho/%s/RhoCent_AccResSysErrors_eta1_eta1_PolySys_FirstOrder_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
  if(mOrder == 2) inputfilecent = Form("StRoot/Utility/Rho/%s/RhoCent_AccResSysErrors_eta1_eta1_PolySys_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
  if(mOrder == 1) inputfiley = Form("StRoot/Utility/Rho/%s/RhoEta_AccResSysErrors_eta1_eta1_PolySys_FirstOrder_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
  if(mOrder == 2) inputfiley = Form("StRoot/Utility/Rho/%s/RhoEta_AccResSysErrors_eta1_eta1_PolySys_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());

  //std::string inputfilecent = "StRoot/Utility/Rho/";
  //std::string inputfiley = "StRoot/Utility/Rho/";
  std::string inputfileptglobal = Form("StRoot/Utility/Rho/%s/AccPhiPtSys_eta1_eta1_PolySys_Global_2D_OffDiag.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *fileptglobal = TFile::Open(inputfileptglobal.c_str());
  TGraphAsymmErrors *g_rho00_global = (TGraphAsymmErrors*) fileptglobal->Get("rhoRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_real_global = (TGraphAsymmErrors*) fileptglobal->Get("realRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_imag_global = (TGraphAsymmErrors*) fileptglobal->Get("imagRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_rerho1n1_global = (TGraphAsymmErrors*) fileptglobal->Get("rerho1n1Raw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_imrho1n1_global = (TGraphAsymmErrors*) fileptglobal->Get("imrho1n1Raw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  
  std::string inputfilepthelicity = Form("StRoot/Utility/Rho/%s/AccPhiPtSys_eta1_eta1_PolySys_Helicity_2D_OffDiag.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *filepthelicity = TFile::Open(inputfilepthelicity.c_str());
  TGraphAsymmErrors *g_rho00_helicity = (TGraphAsymmErrors*) filepthelicity->Get("rhoRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_real_helicity = (TGraphAsymmErrors*) filepthelicity->Get("realRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_imag_helicity = (TGraphAsymmErrors*) filepthelicity->Get("imagRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_rerho1n1_helicity = (TGraphAsymmErrors*) filepthelicity->Get("rerho1n1Raw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_imrho1n1_helicity = (TGraphAsymmErrors*) filepthelicity->Get("imrho1n1Raw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");

  std::string yinputfilepthelicity = Form("StRoot/Utility/Rho/%s/AccPhiPtSys_eta1_eta1_PolySys_Helicity_2D_OffDiag_Rapidity.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *yfilepthelicity = TFile::Open(yinputfilepthelicity.c_str());
  TGraphAsymmErrors *g_rho00_helicity_y = (TGraphAsymmErrors*) yfilepthelicity->Get("rhoRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_real_helicity_y = (TGraphAsymmErrors*) yfilepthelicity->Get("realRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_imag_helicity_y = (TGraphAsymmErrors*) yfilepthelicity->Get("imagRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_rerho1n1_helicity_y = (TGraphAsymmErrors*) yfilepthelicity->Get("rerho1n1Raw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");
  TGraphAsymmErrors *g_imrho1n1_helicity_y = (TGraphAsymmErrors*) yfilepthelicity->Get("imrho1n1Raw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1");

  TFile *filept   = TFile::Open(inputfilept.c_str());
  TFile *filecent = TFile::Open(inputfilecent.c_str());
  TFile *filey    = TFile::Open(inputfiley.c_str());

  TGraphAsymmErrors *g_rho_pt = (TGraphAsymmErrors*) filept->Get("rhoRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Inte_Poly1_F_0_Eff_0");
  //TGraphAsymmErrors *g_rho_pt = (TGraphAsymmErrors*) filept->Get(Form("g_rho00_order%d_%s_%s_StatError",mOrder,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()));
  TGraphAsymmErrors *g_rho_y = (TGraphAsymmErrors*) filey->Get(Form("g_rho00_order%d_%s_%s_StatError",mOrder,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()));
  g_rho_y->Print();
 
  f_rhoy = new TF1("f_rhoy","[0]+[1]*x*x+[2]*x*x*x*x",0,1);
  f_rhoy->SetParameter(0,1/3);
  f_rhoy->SetParameter(1,-0.025);
  f_rhoy->SetParameter(2,0.2);
  g_rho_y->Fit(f_rhoy,"NMRI");

  f_mRhoY = new TF1("f_mRhoY",SpinDensity,-1.0,1.0,2);
  f_mRhoY->FixParameter(1,0.75);

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
    g_rho_pt->GetPoint(ipt,pt,rho);
    cout << "ipt = " << ipt << ", pt = " << pt << ", rho = " << rho << endl;
    f_mRhoPt[ipt]->FixParameter(0,rho);// set by user
    f_mRhoPt[ipt]->FixParameter(1,0.75);

    mMaxData[ipt] = f_mRhoPt[ipt]->GetMaximum(-1.0,1.0);
    cout << "Maximum of data rh00 = " << mMaxData[ipt] << endl; 

    //f_mRhoPt_Helicity[ipt]->FixParameter(0,helicityrho00[ipt]); // set by data
    f_mRhoPt_Helicity[ipt]->FixParameter(0,hrho00);       // set by user
    f_mRhoPt_Helicity[ipt]->FixParameter(1,0.75);
  
    mMaxHelicity[ipt] = f_mRhoPt_Helicity[ipt]->GetMaximum(-1.0,1.0);
    cout << "Maximum of helicity rh00 = " << mMaxHelicity[ipt] << endl; 

    double ptv;
    double rho00_data, real_data, imag_data, rerho1n1_data, imrho1n1_data;
    double rho00_data_err, real_data_err, imag_data_err, rerho1n1_data_err, imrho1n1_data_err;
    g_rho00_global->GetPoint(ipt-2,ptv,rho00_data); 
    g_real_global->GetPoint(ipt-2,ptv,real_data);
    g_imag_global->GetPoint(ipt-2,ptv,imag_data);
    g_rerho1n1_global->GetPoint(ipt-2,ptv,rerho1n1_data);
    g_imrho1n1_global->GetPoint(ipt-2,ptv,imrho1n1_data);

    //f_mRhoPt_2D[ipt]->FixParameter(0,rho);
    //f_mRhoPt_2D[ipt]->FixParameter(0,rho00_data);
    //f_mRhoPt_2D[ipt]->FixParameter(1,real_data);
    //f_mRhoPt_2D[ipt]->FixParameter(2,imag_data);
    //f_mRhoPt_2D[ipt]->FixParameter(3,rerho1n1_data);
    //f_mRhoPt_2D[ipt]->FixParameter(4,imrho1n1_data);

    //temp->FixParameter(0,rho00_data);
    //temp->FixParameter(1,real_data);
    //temp->FixParameter(2,imag_data);
    //temp->FixParameter(3,rerho1n1_data);
    //temp->FixParameter(4,imrho1n1_data);

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

    double hpt;
    double hrho00_data, hreal_data, himag_data, hrerho1n1_data, himrho1n1_data;
    double hrho00_data_err, hreal_data_err, himag_data_err, hrerho1n1_data_err, himrho1n1_data_err;
    g_rho00_helicity->GetPoint(ipt-2,hpt,hrho00_data); 
    g_real_helicity->GetPoint(ipt-2,hpt,hreal_data);
    g_imag_helicity->GetPoint(ipt-2,hpt,himag_data);
    g_rerho1n1_helicity->GetPoint(ipt-2,hpt,hrerho1n1_data);
    g_imrho1n1_helicity->GetPoint(ipt-2,hpt,himrho1n1_data);

    hrho00_data_err = g_rho00_helicity->GetErrorYhigh(ipt-2); 
    hreal_data_err = g_real_helicity->GetErrorYhigh(ipt-2);
    himag_data_err = g_imag_helicity->GetErrorYhigh(ipt-2);
    hrerho1n1_data_err = g_rerho1n1_helicity->GetErrorYhigh(ipt-2);
    himrho1n1_data_err = g_imrho1n1_helicity->GetErrorYhigh(ipt-2);

    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(0,hrho00_data);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(1,hreal_data);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(2,himag_data);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(3,hrerho1n1_data);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(4,himrho1n1_data);

    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(0,hrho00_data+hrho00_data_err);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(1,hreal_data-hreal_data_err);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(2,himag_data-himag_data_err);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(3,hrerho1n1_data-hrerho1n1_data_err);
    //f_mRhoPt_Helicity_2D[ipt]->FixParameter(4,himrho1n1_data-himrho1n1_data_err);

    cout << "ipt = " << ipt << endl;
    cout << "helicity rho00 = " << hrho00_data << endl;
    cout << "helicity real = " << hreal_data << endl;
    cout << "helicity imag = " << himag_data << endl;
    cout << "helicity rerho1n1 = " << hrerho1n1_data << endl;
    cout << "helicity imrho1n1 = " << himrho1n1_data << endl;

    //htemp->FixParameter(0,hrho00_data);
    //htemp->FixParameter(1,hreal_data);
    //htemp->FixParameter(2,himag_data);
    //htemp->FixParameter(3,hrerho1n1_data);
    //htemp->FixParameter(4,himrho1n1_data);

    f_mRhoPt_Helicity_2D[ipt]->FixParameter(0,hrho00);
    f_mRhoPt_Helicity_2D[ipt]->FixParameter(1,hreal);
    f_mRhoPt_Helicity_2D[ipt]->FixParameter(2,himag);
    f_mRhoPt_Helicity_2D[ipt]->FixParameter(3,hrerho1n1);
    f_mRhoPt_Helicity_2D[ipt]->FixParameter(4,himrho1n1);

    htemp->FixParameter(0,hrho00);
    htemp->FixParameter(1,hreal);
    htemp->FixParameter(2,himag);
    htemp->FixParameter(3,hrerho1n1);
    htemp->FixParameter(4,himrho1n1);

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

}

void StEffMcPhi::Make()
{
  long start_event_use = mStartEvent;
  long stop_event_use  = mStopEvent;
  
  double v2_7GeV[9] = {0.125818,0.125818,0.125818,0.125818,0.039951,0.039951,0.039951,0.013462,0.013462};

  gRandom->SetSeed();
  mNtuple->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry
  TF1* f_flow = new TF1("f_flow",flowSampleNorm,-TMath::Pi(),TMath::Pi(),1);

  //string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/ptyspectra_datarcratio_%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  //TFile *File_Spectra = TFile::Open(inputfile.c_str());
  //TH2D *h_mSpectraRatio = (TH2D*) ((TH2D*) File_Spectra->Get("pty_datarcratio"))->Clone();
  //double spectramax = h_mSpectraRatio->GetMaximum();
  //h_mSpectraRatio->Scale(1./spectramax);
  //cout << "spectra max = " << spectramax << ", normalzied max " << h_mSpectraRatio->GetMaximum() << endl;

  for(long i_track = start_event_use; i_track < stop_event_use; ++i_track)
  {
    if (!mNtuple->GetEntry(i_track)) // take track information
      break;  // end of data chunk
    if (floor(10.0*i_track/ static_cast<float>(stop_event_use)) > floor(10.0*(i_track-1)/ static_cast<float>(stop_event_use)))
      cout << "=> processing data: " << 100.0*i_track/ static_cast<float>(stop_event_use) << "%" << endl;

    McVecMeson McPhi; // initialize McPhi
    McPhi.Centrality = mCentrality;
    if(mMode == 0 && (( McPhi.Centrality > 5 || McPhi.Centrality < 2 ))) continue;
    if(mMode == 4 && (( McPhi.Centrality > 5 || McPhi.Centrality < 2 ))) continue;
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
    RcKM.RcTof = mKmRcTof;

    RcVecMeson RcPhi; // initialize RcPhi
    RcPhi.RcPt      = mRcPt;
    RcPhi.RcP       = mRcP;
    RcPhi.RcEta     = mRcEta;
    RcPhi.RcY       = mRcY;
    RcPhi.RcPhi     = mRcPhi;
    RcPhi.RcInvMass = mRcInvMass;
    RcPhi.RcTpc = mRcTpc;


    f_flow->ReleaseParameter(0);
    f_flow->SetParameter(0,mv2);//f_mV2[McPhi.Centrality]->Eval(McPhi.McPt));
    //f_flow->SetParameter(0,f_mV2[int(McPhi.Centrality)]->Eval(McPhi.McPt));
    //if(energy == 0) f_flow->SetParameter(0,v2_7GeV[int(McPhi.Centrality)]);//f_mV2[McPhi.Centrality]->Eval(McPhi.McPt));
    //cout << "Set the parameter" << endl;
    double maximum = f_flow->Eval(0.0);
    //cout << "Maximum = " << maximum << endl;
    double phiminusPsi = McPhi.McPhi-McPhi.Psi; 
    //cout << "grabbed maximum and phiminusPsi = " << phiminusPsi<< endl;
    while(phiminusPsi < -TMath::Pi()) phiminusPsi += 2.0*TMath::Pi();
    while(phiminusPsi >  TMath::Pi()) phiminusPsi -= 2.0*TMath::Pi();
    //cout << "Angle Wrapping = " << phiminusPsi << endl;
    double valuev2 = f_flow->Eval(phiminusPsi);
    double randomv2 = gRandom->Uniform(0,1);
    if(randomv2 > valuev2/maximum) continue; 

    //cout << "maximum = " << maximum << endl;
    //cout << "valuev2/maximum = " << valuev2/maximum << endl;    

    //cout << "mMode = " << mMode << endl;
    //cout << "Centrality = " << McPhi.Centrality << endl;

    //cout << "McPhi.McPt = " << McPhi.McPt << ", RcPhi.RcPt = " << RcPhi.RcPt << endl;
    //cout << "McPhi.McY = " << McPhi.McY << ", RcPhi.RcY = " << RcPhi.RcY << endl;

    if( mMode == 0 && (McPhi.McPt <= vmsa::pt_low[energy][mStartPt] || McPhi.McPt > vmsa::pt_up[energy][mStopPt]))    continue;
    if( mMode == 4 && (McPhi.McPt <= vmsa::pt_low[energy][mStartPt] || McPhi.McPt > vmsa::pt_up[energy][mStopPt]))    continue;
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
     

    if( mEtaMode == 0 && TMath::Abs(McPhi.McY) > 1.0 ) continue;
    //if( mEtaMode == 0 && TMath::Abs(RcPhi.RcY) > 1.0 ) continue; //Swap to RC level cut
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

    TLorentzVector lMcPhi;
    lMcPhi.SetPtEtaPhiM(McPhi.McPt,McPhi.McEta,McPhi.McPhi,vmsa::InvMass[pid]);
    TVector3 vMcPhiBeta = -1.0*lMcPhi.BoostVector();
    TVector3 phiMomentumLabUnit = lMcPhi.Vect().Unit();
     
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
    //cout << "phiMomentumLabUnit" << endl;
    //phiMomentumLabUnit.Print();
    //cout << "vBeamTotalUnit" << endl;
    //vBeamTotalUnit.Print();
    //cout << endl;

    TVector3 zaxisfrombeam = vBeamTotalUnit;
    TVector3 xaxisfrombeam = yaxisfrombeam.Cross(zaxisfrombeam).Unit();

    TLorentzVector lMcKP;
    lMcKP.SetPtEtaPhiM(McKP.McPt,McKP.McEta,McKP.McPhi,vmsa::mMassKaon);
    lMcKP.Boost(vMcPhiBeta);
    TVector3 vMcKP = lMcKP.Vect().Unit(); // direction of K+ momentum in phi-meson rest frame
 
    double costheta = vMcKP.Dot(zaxisfrombeam); 

    double xprojection = vMcKP.Dot(xaxisfrombeam);
    double yprojection = vMcKP.Dot(yaxisfrombeam);
  
    double helicityangle = TMath::ATan2(yprojection,xprojection);
    if(helicityangle < 0.0) helicityangle += 2.0*TMath::Pi();
  
    //TVector3 QVector(TMath::Sin(Psi),-1.0*TMath::Cos(Psi),0.0);
    TVector3 QVectorMc(TMath::Sin(McPhi.Psi),-1.0*TMath::Cos(McPhi.Psi),0.0);
    double phistar = vMcKP.Phi();
    TVector3 mcxprime(TMath::Cos(McPhi.Psi),TMath::Sin(McPhi.Psi),0.0);
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


    TLorentzVector lRcPhi;
    lRcPhi.SetPtEtaPhiM(RcPhi.RcPt,RcPhi.RcEta,RcPhi.RcPhi,RcPhi.RcInvMass);
    TVector3 vRcPhiBeta = -1.0*lRcPhi.BoostVector();
    TVector3 phiMomentumLabUnitRc = lRcPhi.Vect().Unit();

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
    //cout << "phiMomentumLabUnitRc" << endl;
    //phiMomentumLabUnitRc.Print();
    //cout << "vBeamTotalUnitRc" << endl;
    //vBeamTotalUnitRc.Print();
    //cout << endl;
    
    TVector3 zaxisfrombeamRc = vBeamTotalUnitRc;
    TVector3 xaxisfrombeamRc = yaxisfrombeamRc.Cross(zaxisfrombeamRc).Unit();

    TLorentzVector lRcKP;
    lRcKP.SetPtEtaPhiM(RcKP.RcPt,RcKP.RcEta,RcKP.RcPhi,vmsa::mMassKaon);
    lRcKP.Boost(vRcPhiBeta);
    TVector3 vRcKP = lRcKP.Vect().Unit(); // direction of K+ momentum in phi-meson rest frame
 
    double costhetaRc = vRcKP.Dot(zaxisfrombeamRc); 

    double xprojectionRc = vRcKP.Dot(xaxisfrombeamRc);
    double yprojectionRc = vRcKP.Dot(yaxisfrombeamRc);
  
    double helicityangleRc = TMath::ATan2(yprojectionRc,xprojectionRc);
    if(helicityangleRc < 0.0) helicityangleRc += 2.0*TMath::Pi();


    TVector3 QVectorRc(TMath::Sin(Psi),-1.0*TMath::Cos(Psi),0.0);
    double phistarRc = vRcKP.Phi();
    TVector3 xprime(TMath::Cos(Psi),TMath::Sin(Psi),0.0);
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
    
    TVector3 nQMc = QVectorMc.Unit(); // direction of QVector
    TVector3 nQRc = QVectorRc.Unit(); // direction of QVector
    //float McThetaStar = TMath::ACos(vMcKP.Dot(nQMc));
    //if(McThetaStar < 0.0) McThetaStar+=TMath::Pi();
    //if(McThetaStar > TMath::Pi()) McThetaStar-=TMath::Pi();
    float McCosThetaStar = vMcKP.Dot(nQMc);
    float RcCosThetaStar = vRcKP.Dot(nQRc);
   
   
 
    //cout << "mc costheta* - rc costheta* = " << McCosThetaStar-RcCosThetaStar << endl;
    
    int ptbin = -1; 
    int centbin = -1;
    int ybin = -1;
    if(mMode == 0)
    {
      for(int ipt = 2; ipt < 6; ipt++)
      {
        //if(RcPhi.RcPt > vmsa::pt_low[energy][ipt] && RcPhi.RcPt <= vmsa::pt_up[energy][ipt]) 

        if(McPhi.McPt > vmsa::pt_low[energy][ipt] && McPhi.McPt <= vmsa::pt_up[energy][ipt]) 
        {
          ptbin = ipt; 
          //cout << "ptbin = " << ptbin << ", pt = " << McPhi.McPt << endl;
          break;
        }
      }
      for(int iy = 0; iy < 10; iy++)
      {
        if(McPhi.McY >= float(iy-5)/5. && McPhi.McY < float(iy-4)/5.)
        {
          ybin = iy;
          //cout << "y = " << McPhi.McY << ", ybin = " << ybin << endl;
          break;
        }
        else if ( McPhi.McY ==  1.)
        {
          ybin = 9; 
        }
      }
    }
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
    if(mMode == 4)
    {
      for(int iy = 0; iy < 10; iy++)
      {
        if(McPhi.McY >= float(iy-5)/5. && McPhi.McY < float(iy-4)/5.)
        {
          ybin = iy;
          //cout << "y = " << McPhi.McY << ", ybin = " << ybin << endl;
          break;
        }
        else if ( McPhi.McY ==  1.)
        {
          ybin = 9; 
        }
      }
    }

    ////cout << "Before Sampling" << endl;
    ////cout << "ptbin = " << ptbin << ", centbin = " << centbin << ", ybin = " << ybin << endl;
    //cout << "Before Sampling" << endl;
    
    if(mMode == 0)
    {
      //double rho00 = f_rhoy->Eval(McPhi.McY);
      ////cout << "y = " << McPhi.McY << ", rho00 = " << rho00 << endl;
      //f_mRhoY->FixParameter(0,rho00);
      //double maxy = f_mRhoY->GetMaximum(-1.0,1.0);
      ////cout << "maxy = " << maxy << endl;
      // 
      //if(!Sampling(f_mRhoY,TMath::Abs(McCosThetaStar),maxy)) continue;

      //if(!Sampling(f_mRhoPt[ptbin],TMath::Abs(McCosThetaStar),mMaxData)) continue;
      //if(!Sampling(f_mRhoPt[ptbin],McCosThetaStar,mMaxData[ptbin])) continue;
      //cout << "Are we sampling? " << endl;
      //cout << "McCosThetaStar = " << McCosThetaStar << " mcphiprime = " << mcphiprime << endl;
      if(!Sampling2D(f_mRhoPt_2D[ptbin],McCosThetaStar,mcphiprime,mMax2D[ptbin])) continue;
      if(!Sampling2D(f_mRhoPt_Helicity_2D[ptbin],costheta,helicityangle,mMaxHelicity2D[ptbin])) continue;
      //if(!Sampling2D(f_mRhoY_Helicity_2D[ybin],costheta,helicityangle,mMaxHelicity2DY[ybin])) continue;
      //if(!SamplingHelicity(f_mRhoPt_Helicity[ptbin],costheta,mMaxHelicity[ptbin])) continue;
      //cout << "pass sampling" << endl;
    }
    if(mMode == 4) 
    {
      //if(!Sampling2D(f_mRhoY_Helicity_2D[ybin],costheta,helicityangle,mMaxHelicity2DY[ybin])) continue;
    }
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
    //mEffHistManger->FillHistMc(McPhi.Centrality,McPhi.McPt,McPhi.McY,McPhi.McPhi,McCosThetaStar,Psi,phistar,mcphiprime,McKP.McPt,McKP.McY,costheta,helicityangle);
    //mEffHistManger->FillHistMc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,McCosThetaStar,Psi,phistar,mcphiprime,RcKP.RcPt,RcKP.RcY,RcKM.RcPt,RcKM.RcY,costheta,helicityangle);
    mEffHistManger->FillHistMc(McPhi.Centrality,McPhi.McPt,McPhi.McY,McPhi.McPhi,McCosThetaStar,Psi,phistar,mcphiprime,McKP.McPt,McKP.McY,McKM.McPt,McKM.McY,costheta,helicityangle);
    //mEffHistManger->FillHistMc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc);
    //mEffHistManger->FillKaonHistMc(McPhi.Centrality,McPhi.McPhi,McPhi.McPt,McPhi.McY,McKP.McPt,McKP.McY,McKP.McEta,McKP.McPhi,phistar,McCosThetaStar,costheta); // K+
    //mEffHistManger->FillKaonHistMc(McPhi.Centrality,McPhi.McPhi,McPhi.McPt,McPhi.McY,McKM.McPt,McKM.McY,McKM.McEta,McKM.McPhi,phistar,McCosThetaStar,costheta); // K-
    //mEffHistManger->FillKaonDeltaHistMc(McPhi.Centrality,McPhi.McPhi,McPhi.McEta,McPhi.McPt,McPhi.McY,McKP.McPt,McKP.McEta,McKP.McPhi,McKM.McPt,McKM.McEta,McKM.McPhi,phistar); // K-
    //cout << "After MC Fill 5" << endl;
 
    if(fabs(RcPhi.RcY) > 1.0) continue;
    mEffHistManger->FillHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc,RcKP.RcPt,RcKP.RcY,RcKM.RcPt,RcKM.RcY,costhetaRc,helicityangleRc,0);

    //if( !mEffCut->passTrackCutPhi(RcPhi) ) continue;  // eta cuts for RcPhi 
    //cout << "Right before TPC Reconstructed" << endl;
    //cout << "Right after TPC Reconstructed" << endl;

    if(mEtaMode == 0)
    {
      //if( fabs(McKP.McEta) > 0.5 || fabs(McKM.McEta) > 0.5 ) continue; // eta cuts for McPhi 
      if( fabs(RcKP.RcEta) > 1.0 || fabs(RcKM.RcEta) > 1.0 ) continue; // eta cuts for McPhi 
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

    mEffHistManger->FillHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc,RcKP.RcPt,RcKP.RcY,RcKM.RcPt,RcKM.RcY,costhetaRc,helicityangleRc,1);

    if( !tpcReconstructed(0,int(McPhi.Centrality),&lRcPhi) ) continue;
    mEffHistManger->FillHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc,RcKP.RcPt,RcKP.RcY,RcKM.RcPt,RcKM.RcY,costhetaRc,helicityangleRc,2);


    //cout << "Right Before RC Fill " << endl;
    //if( !mEffCut->passTrackCutTpc(RcKP) ) continue; // eta and TPC cuts for RcKplus
    //if( !mEffCut->passTrackCutTpc(RcKM) ) continue; // eta and TPC cuts for RcKminus
    //mEffHistManger->FillHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc,RcKP.RcPt,RcKP.RcY,RcKM.RcPt,RcKM.RcY,costhetaRc,helicityangleRc,3);

    if( !mEffCut->passTrackCutTof(RcKP) ) continue; // eta and TPC cuts for RcKplus
    if( !mEffCut->passTrackCutTof(RcKM) ) continue; // eta and TPC cuts for RcKminus
    mEffHistManger->FillHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc,RcKP.RcPt,RcKP.RcY,RcKM.RcPt,RcKM.RcY,costhetaRc,helicityangleRc,3);

    TLorentzVector rl_km;
    rl_km.SetPtEtaPhiM(RcKM.RcPt,RcKM.RcEta,RcKM.RcPhi,vmsa::mMassKaon);
    int EtaBinM = -1;
    int PhiBinM = -1;
    TLorentzVector rl_kp;
    rl_kp.SetPtEtaPhiM(RcKP.RcPt,RcKP.RcEta,RcKP.RcPhi,vmsa::mMassKaon);
    int EtaBinP = -1;
    int PhiBinP = -1;

    findFuncPID(rl_km,0,EtaBinM,PhiBinM);
    findFuncPID(rl_kp,1,EtaBinP,PhiBinP);

    if( !pass_m2_PID(0,9/*McPhi.Centrality*/,rl_km,EtaBinM,PhiBinM) ) continue;
    if( !pass_m2_PID(1,9/*McPhi.Centrality*/,rl_kp,EtaBinP,PhiBinP) ) continue;
    mEffHistManger->FillHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc,RcKP.RcPt,RcKP.RcY,RcKM.RcPt,RcKM.RcY,costhetaRc,helicityangleRc,4);

    if( !pass_nsig_PID(0,9/*McPhi.Centrality*/,rl_km,EtaBinM,PhiBinM) ) continue;
    if( !pass_nsig_PID(1,9/*McPhi.Centrality*/,rl_kp,EtaBinP,PhiBinP) ) continue;
    mEffHistManger->FillHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc,RcKP.RcPt,RcKP.RcY,RcKM.RcPt,RcKM.RcY,costhetaRc,helicityangleRc,5);
    
    //mEffHistManger->FillPhiHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcPhi,RcPhi.RcY,phistarRc,RcCosThetaStar,costhetaRc,RcKP.RcPt,RcKP.RcY);
    //mEffHistManger->FillHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc,RcKP.RcPt,RcKP.RcY,costhetaRc,helicityangleRc);
    //mEffHistManger->FillKaonHistRc(McPhi.Centrality,RcPhi.RcPhi,RcPhi.RcPt,RcPhi.RcY,RcKP.RcPt,RcKP.RcY,RcKP.RcEta,RcKP.RcPhi,phistarRc,RcCosThetaStar,costhetaRc); // K+
    //mEffHistManger->FillKaonHistRc(McPhi.Centrality,RcPhi.RcPhi,RcPhi.RcPt,RcPhi.RcY,RcKM.RcPt,RcKM.RcY,RcKM.RcEta,RcKM.RcPhi,phistarRc,RcCosThetaStar,costhetaRc); // K-
    //mEffHistManger->FillKaonDeltaHistRc(McPhi.Centrality,RcPhi.RcPhi,RcPhi.RcEta,RcPhi.RcPt,RcPhi.RcY,RcKP.RcPt,RcKP.RcEta,RcKP.RcPhi,RcKM.RcPt,RcKM.RcEta,RcKM.RcPhi,phistarRc); // K-

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
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;
  mEffHistManger->CalEffCosThetaStar();
  cout << "calculated efficiency" << endl;
}

void StEffMcPhi::Finish()
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

