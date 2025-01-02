#include "StEffMcPhi.h"
#include "StEffHistManger.h"
#include "StEffCut.h"
#include <string>
#include "TFile.h"
#include "TNtuple.h"
#include "TBranch.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StRoot/Utility/functions.h"
#include "TH2D.h"
ClassImp(StEffMcPhi)

int StEffMcPhi::mInput_flag = 1;

bool Sampling(TF1 *f_rhoPhy, float CosThetaStar);
bool Sampling2D(TF2 *f_rhoPhy, float ThetaStar, float beta, float max);

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

StEffMcPhi::StEffMcPhi(int Energy, long StartEvent, long StopEvent, int PID, int year, int mode, int inputpt, int startpt, int stoppt, const char* setting, int etamode, int order = 2, float realrho1n1 = 0.0, float rho00 = 1./3., float reterms = 0.0, float imterms = 0.0, float imagrho1n1 = 0.0) 
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

  //for(int i = 0; i < 9; i++)
  //{
  //  f_mV2[i] = readv2(energy,pid,i);
  //}
  
  std::string EP[2] = {"","2nd"};

  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_prelimv2_EPeff_randomRP/%s_%s_yabs1_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  //string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_PhiEffFinerBins_prelimv2_flatRP_yabs1p1/%s_%s_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
  string InPutFile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/EffAcc_NoRapiditySpectra_nov2_flatRP_yabs1p1/%s_%s_pt%d.root",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),inputpt);
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

  //std::string inputfilept     = Form("StRoot/Utility/Rho/%s/Rho_AccResSysErrors_F_0_eta1_eta1_PolySys_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
  //if(mOrder == 1) inputfilept = Form("StRoot/Utility/Rho/%s/Rho_AccResSysErrors_F_0_eta1_eta1_PolySys_FirstOrder_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
  //std::string inputfilecent;
  //std::string inputfiley;
  //if(mOrder == 1) inputfilecent = Form("StRoot/Utility/Rho/%s/RhoCent_AccResSysErrors_eta1_eta1_PolySys_FirstOrder_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
  //if(mOrder == 2) inputfilecent = Form("StRoot/Utility/Rho/%s/RhoCent_AccResSysErrors_eta1_eta1_PolySys_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
  //if(mOrder == 1) inputfiley = Form("StRoot/Utility/Rho/%s/RhoEta_AccResSysErrors_eta1_eta1_PolySys_FirstOrder_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());
  //if(mOrder == 2) inputfiley = Form("StRoot/Utility/Rho/%s/RhoEta_AccResSysErrors_eta1_eta1_PolySys_NoRapiditySpectra_FixedFirstEP.root",vmsa::mBeamEnergy[energy].c_str());

  ////std::string inputfilecent = "StRoot/Utility/Rho/";
  ////std::string inputfiley = "StRoot/Utility/Rho/";
  //
  //TFile *filept   = TFile::Open(inputfilept.c_str());
  //TFile *filecent = TFile::Open(inputfilecent.c_str());
  //TFile *filey    = TFile::Open(inputfiley.c_str());

  //TGraphAsymmErrors *g_rho_pt = (TGraphAsymmErrors*) filept->Get(Form("g_rho00_order%d_%s_%s_StatError",mOrder,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()));
  //for(int ipt = 2; ipt < 6; ipt++)
  //{ 
  //  f_mRhoPt[ipt] = new TF1(Form("f_mRho_pt%d",ipt),SpinDensity,-1.0,1.0,2);
  //  f_mRhoPt_2D[ipt] = new TF2(Form("f_mRho2D_pt%d",ipt),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),5);
  //  TF2 *temp = new TF2("temp",SpinDensity2Dcosneg,0.0,TMath::Pi(),0.0,2.0*TMath::Pi(),5);
  //  double pt, rho; 
  //  g_rho_pt->GetPoint(ipt,pt,rho);
  //  cout << "ipt = " << ipt << ", pt = " << pt << ", rho = " << rho << endl;
  //  //cout << "ipt = " << ipt << ", pt = " << pt << ", rho = " << rho << endl;
  //  //f_mRhoPt[ipt]->FixParameter(0,0.5);
  //  //f_mRhoPt[ipt]->FixParameter(0,rho);
  //  f_mRhoPt[ipt]->FixParameter(0,mrho00);
  //  f_mRhoPt[ipt]->FixParameter(1,0.75);

  //  cout << "rho00 = " << mrho00 << ", rerho1n1 = " << rerho1n1 << endl;

  //  f_mRhoPt_2D[ipt]->FixParameter(0,mrho00);
  //  //f_mRhoPt_2D[ipt]->FixParameter(0,rho);
  //  f_mRhoPt_2D[ipt]->FixParameter(1,real);
  //  f_mRhoPt_2D[ipt]->FixParameter(2,imag);
  //  f_mRhoPt_2D[ipt]->FixParameter(3,rerho1n1);
  //  f_mRhoPt_2D[ipt]->FixParameter(4,imrho1n1);

  //  temp->FixParameter(0,mrho00);
  //  temp->FixParameter(1,real);
  //  temp->FixParameter(2,imag);
  //  temp->FixParameter(3,rerho1n1);
  //  temp->FixParameter(4,imrho1n1);

  //  //mMax = f_mRhoPt_2D[ipt]->GetMaximum();
  //  //double x[2] = {TMath::Pi()/2.,TMath::Pi()/2.};
  //  double x,y;
  //  temp->GetMinimumXY(x,y);
  //  cout << "Maximum x,y = " << x << "," << y << endl;
  //  mMax = TMath::Abs(temp->Eval(x,y));
  //  cout << "Maximum of the function is " << mMax << endl;

  //  //cout << "x,y both = TMath::Pi()/2, eval func = " << temp->Eval(1.57113,1.5708) << endl;
  //}  
 
  //for(int ipt = 0; ipt < vmsa::pt_rebin_cent; ipt++)
  //{
  //  TGraphAsymmErrors *g_rho_cent = (TGraphAsymmErrors*) filecent->Get(Form("rhoRawStat_pt_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",ipt,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str()));
  //  for(int icent = 0; icent < 9; icent++)
  //  {
  //    f_mRhoCent[ipt][icent] = new TF1(Form("f_mRho_pt%d_cent%d",ipt,icent),SpinDensity,-1.0,1.0,2);
  //    double cent, rho; 
  //    g_rho_cent->GetPoint(icent,cent,rho);
  //    cout << "ipt = " << ipt << ", cent = " << cent << ", rho = " << rho << endl;
  //    f_mRhoCent[ipt][icent]->FixParameter(0,rho);
  //    f_mRhoCent[ipt][icent]->FixParameter(1,0.75);
  //  }
  //}    

  //for(int ipt = 0; ipt < vmsa::pt_rebin_y; ipt++)
  //{
  //  for(int icent = 0; icent < vmsa::cent_rebin_total; icent++)
  //  {
  //    TGraphAsymmErrors *g_rho_y = (TGraphAsymmErrors*) filey->Get(Form("rhoRawStat_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,icent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1));
  //    for(int iy = 0; iy < 10; iy++)
  //    { 
  //      f_mRhoY[ipt][icent][iy] = new TF1(Form("f_mRho_pt%d_cent%d_y%d",ipt,icent,iy),SpinDensity,-1.0,1.0,2);
  //      double y, rho; 
  //      g_rho_y->GetPoint(iy,y,rho);
  //      cout << "ipt = " << ipt << ", icent = " << icent << ", y = " << y << ", rho = " << rho << endl;
  //      f_mRhoY[ipt][icent][iy]->FixParameter(0,rho);
  //      f_mRhoY[ipt][icent][iy]->FixParameter(1,0.75);
  //    }
  //  }
  //}    
 
  //for(int ipt = 0; ipt < vmsa::pt_rebin_cent; ipt++)
  //{ 
  //  TGraphAsymmErrors *g_rho_cent = (TGraphAsymmErrors*) filecent->GetEntry(Form(""));
  //  for(int icent = 0; icent < 9; icent++)
  //  { 
  //    f_mRhoCent[ipt][icent] = new TF1(Form("f_mRho_pt%d_cent%d",ipt,cent),SpinDensity,-1.0,1.0,2);
  //    double cent, rho; 
  //    g_rho_cent->GetPoint(icent,cent,rho);
  //    cout << "icent = " << icent << ", cent = " << cent << ", rho = " << rho << endl;
  //    f_mRhoCent[ipt][icent]->FixParameter(0,rho);
  //    f_mRhoCent[ipt][icent]->FixParameter(1,0.75);
  //}  

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
  mEffHistManger->InitKaonHist();
  mEffHistManger->InitPhiHist();


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
  f_y = new TF1("f_y",Form("exp(-x*x/2/%s)",mSigmay.c_str()),-1.0,1.0);

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


    f_flow->ReleaseParameter(0);
    //f_flow->SetParameter(0,0.2);//f_mV2[McPhi.Centrality]->Eval(McPhi.McPt));
    if(energy == 0) f_flow->SetParameter(0,v2_7GeV[int(McPhi.Centrality)]);//f_mV2[McPhi.Centrality]->Eval(McPhi.McPt));
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

    TLorentzVector lMcPhi;
    lMcPhi.SetPtEtaPhiM(McPhi.McPt,McPhi.McEta,McPhi.McPhi,vmsa::InvMass[pid]);
    TVector3 vMcPhiBeta = -1.0*lMcPhi.BoostVector();

    TLorentzVector lMcKP;
    lMcKP.SetPtEtaPhiM(McKP.McPt,McKP.McEta,McKP.McPhi,vmsa::mMassKaon);
    lMcKP.Boost(vMcPhiBeta);
    TVector3 vMcKP = lMcKP.Vect().Unit(); // direction of K+ momentum in phi-meson rest frame
    //TVector3 QVector(TMath::Sin(Psi),-1.0*TMath::Cos(Psi),0.0);
    TVector3 QVectorMc(TMath::Sin(McPhi.Psi),-1.0*TMath::Cos(McPhi.Psi),0.0);
    double phistar = vMcKP.Phi();
    TVector3 mcxprime(TMath::Cos(McPhi.Psi),TMath::Sin(McPhi.Psi),0.0);
    TVector3 mczprime(0.0,0.0,1.0);
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


    TLorentzVector lRcPhi;
    lRcPhi.SetPtEtaPhiM(RcPhi.RcPt,RcPhi.RcEta,RcPhi.RcPhi,RcPhi.RcInvMass);
    TVector3 vRcPhiBeta = -1.0*lRcPhi.BoostVector();

    TLorentzVector lRcKP;
    lRcKP.SetPtEtaPhiM(RcKP.RcPt,RcKP.RcEta,RcKP.RcPhi,vmsa::mMassKaon);
    lRcKP.Boost(vRcPhiBeta);
    TVector3 vRcKP = lRcKP.Vect().Unit(); // direction of K+ momentum in phi-meson rest frame
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
    
    //int ptbin = -1; 
    //int centbin = -1;
    //int ybin = -1;
    //if(mMode == 0)
    //{
    //  for(int ipt = 2; ipt < 6; ipt++)
    //  {
    //    if(RcPhi.RcPt > vmsa::pt_low[energy][ipt] && RcPhi.RcPt <= vmsa::pt_up[energy][ipt]) 
    //    //if(McPhi.McPt > vmsa::pt_low[energy][ipt] && McPhi.McPt <= vmsa::pt_up[energy][ipt]) 
    //    {
    //      ptbin = ipt; 
    //      //cout << "ptbin = " << ptbin << ", pt = " << McPhi.McPt << endl;
    //      break;
    //    }
    //  }
    //}
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

    //////cout << "Before Sampling" << endl;
    //////cout << "ptbin = " << ptbin << ", centbin = " << centbin << ", ybin = " << ybin << endl;

    //if(mMode == 0)
    //{
    //  //if(!Sampling(f_mRhoPt[ptbin],TMath::Abs(McCosThetaStar))) continue;
    //  //cout << "Are we sampling? " << endl;
    //  //cout << "McCosThetaStar = " << McCosThetaStar << " mcphiprime = " << mcphiprime << endl;
    //  if(!Sampling2D(f_mRhoPt_2D[ptbin],McCosThetaStar,mcphiprime,mMax)) continue;
    //  //cout << "pass sampling" << endl;
    //}
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

    //mEffHistManger->FillPhiHistMc(McPhi.Centrality,McPhi.McPt,McPhi.McPhi,McPhi.McY,phistar);
    mEffHistManger->FillHistMc(McPhi.Centrality,McPhi.McPt,McPhi.McY,McPhi.McPhi,McCosThetaStar,Psi,phistar,mcphiprime);
    //mEffHistManger->FillHistMc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc);
    //mEffHistManger->FillKaonHistMc(McPhi.Centrality,McPhi.McPhi,McPhi.McPt,McKP.McPt,McKP.McY,McKP.McEta,McKP.McPhi,phistar); // K+
    //mEffHistManger->FillKaonHistMc(McPhi.Centrality,McPhi.McPhi,McPhi.McPt,McKM.McPt,McKM.McY,McKM.McEta,McKM.McPhi,phistar); // K-
    //mEffHistManger->FillKaonDeltaHistMc(McPhi.Centrality,McPhi.McPhi,McPhi.McEta,McPhi.McPt,McKP.McPt,McKP.McEta,McKP.McPhi,McKM.McPt,McKM.McEta,McKM.McPhi,phistar); // K-

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



    if( !mEffCut->passTrackCut(RcKP) ) continue; // eta and TPC cuts for RcKplus
    if( !mEffCut->passTrackCut(RcKM) ) continue; // eta and TPC cuts for RcKminus
    if( !mEffCut->passTrackCutPhi(RcPhi) ) continue;  // eta cuts for RcPhi 
    //mEffHistManger->FillPhiHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcPhi,RcPhi.RcY,phistarRc);
    mEffHistManger->FillHistRc(McPhi.Centrality,RcPhi.RcPt,RcPhi.RcY,RcPhi.RcPhi,RcCosThetaStar,Psi,phistarRc,phiprimeRc);
    //mEffHistManger->FillKaonHistRc(McPhi.Centrality,RcPhi.RcPhi,RcPhi.RcPt,RcKP.RcPt,RcKP.RcY,RcKP.RcEta,RcKP.RcPhi,phistarRc); // K+
    //mEffHistManger->FillKaonHistRc(McPhi.Centrality,RcPhi.RcPhi,RcPhi.RcPt,RcKM.RcPt,RcKM.RcY,RcKM.RcEta,RcKM.RcPhi,phistarRc); // K-
    //mEffHistManger->FillKaonDeltaHistRc(McPhi.Centrality,RcPhi.RcPhi,RcPhi.RcEta,RcPhi.RcPt,RcKP.RcPt,RcKP.RcEta,RcKP.RcPhi,RcKM.RcPt,RcKM.RcEta,RcKM.RcPhi,phistarRc); // K-

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

bool Sampling(TF1 *f_rhoPhy, float CosThetaStar)
{
  float wMax;
  if(f_rhoPhy->GetParameter(0) <= 1.0/3.0) wMax = f_rhoPhy->Eval(0.0);
  else wMax = f_rhoPhy->Eval(1.0);
  return !(gRandom->Rndm() > f_rhoPhy->Eval(CosThetaStar)/wMax);
}

bool Sampling2D(TF2 *f_rhoPhy, float ThetaStar, float PhiPrime, float wMax)
{
  return !(gRandom->Rndm() > f_rhoPhy->Eval(ThetaStar,PhiPrime)/wMax);
}
