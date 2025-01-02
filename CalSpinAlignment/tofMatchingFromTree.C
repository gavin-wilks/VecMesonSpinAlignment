#include <THn.h>
#include <TH3F.h>
#include <TH3F.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include "../Utility/functions.h"
#include "../Utility/StSpinAlignmentCons.h"

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

TF1* readspec(int energy, int centrality)
{
 
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

  TF1 *f_spec = new TF1("f_spec",pTLevy,0.1,5.0,3);
  TGraphAsymmErrors *g_spec;
  if(energy != 3)
  {
    string InPutSpec;
    if((energy == 4 || energy == 2 || energy == 0)) InPutSpec = "pTspectra/HEPData-ins1378002-v1-root.root";
    TFile *File_Spec = TFile::Open(InPutSpec.c_str());
    cout << "Input spectra" << InPutSpec << endl;
 
    g_spec = (TGraphAsymmErrors*)File_Spec->Get("g_spec");
    if((energy == 4 || energy == 2 || energy == 0)) 
    {
      TDirectory *dir = (TDirectory*) File_Spec->Get(Form("Table %d",tableNum[energy][centrality]));
      dir->cd(); 
      g_spec = (TGraphAsymmErrors*)dir->Get("Graph1D_y1");
      for(int i = 0; i < g_spec->GetN(); i++)
      {
        double x,y;
        g_spec->GetPoint(i,x,y);
        g_spec->SetPoint(i,ptvalues[energy][centrality][i],y);
        g_spec->SetPointError(i,0.0,0.0,g_spec->GetErrorYlow(i),g_spec->GetErrorYhigh(i));
      }
    }

    TF1 *f_Levy = new TF1("f_Levy",Levy,0.1,5.0,3);
    f_Levy->SetParameter(0,1);
    if(energy == 4 && centrality >= 0 && centrality <= 1) f_Levy->SetParameter(0,0.01);
    if(energy == 0 ) f_Levy->SetParameter(0,0.01);
    f_Levy->SetParameter(1,10);
    f_Levy->SetParameter(2,0.1);
    g_spec->Print();
    cout << "Fitting full pT distribution" << endl;
    g_spec->Fit(f_Levy,"N");


    f_spec->SetParameter(0,f_Levy->GetParameter(0));
    f_spec->SetParameter(1,f_Levy->GetParameter(1));
    f_spec->SetParameter(2,f_Levy->GetParameter(2));
  }

  //double dNdy14[9] = {0.009126,0.009126,0.042299,0.042299,0.095731,0.152881,0.223975,0.329114,0.329114};    
  //double Texp14[9] = {0.211483,0.211483,0.228253,0.228253,0.241851,0.245252,0.266930,0.268719,0.268719};

  //TF1 *f_Expo = new TF1("f_Expo",specExp,0.0,5.0,2);
  //if(energy == 3)
  //{
  //  f_spec = new TF1("f_specExp",pTspecExp,0.0,5.0,2);
  //  f_spec->SetParameter(0,dNdy14[centrality]);
  //  f_spec->SetParameter(1,Texp14[centrality]);
  //  f_Expo->SetParameter(0,dNdy14[centrality]);
  //  f_Expo->SetParameter(1,Texp14[centrality]);
  //}
 
  return f_spec;
}

TF1* readv2(int energy, int centrality){
  
  string centlabel = "4080";
  if(centrality >= 4 && centrality <= 6) centlabel = "1040";
  if(centrality >= 7 && centrality <= 8) centlabel = "0010";

  TGraphAsymmErrors *g_v2;
  if(energy != 3)
  {
    string InPutV2 = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
    if((energy == 2 || energy == 0)) InPutV2 = "v2files/HEPData-ins1395151-v2-root.root";
    if(energy == 4) InPutV2 = Form("v2files/OutPhi_v2_Cent%s.root",centlabel.c_str());
    TFile *File_v2 = TFile::Open(InPutV2.c_str());
    std::cout << "v2 file: " << InPutV2 << endl;


    g_v2 = (TGraphAsymmErrors*)File_v2->Get("g_v2");

    if((energy == 2 || energy == 0)) 
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
  
  //if( energy == 3)
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


  TF1 *f_v2 = new TF1("f_v2",v2_pT_FitFunc,0.1,5.0,5);
 // if(centrality <= 3) f_v2->SetRange(0.1,2.8); 
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
  h_v2->GetYaxis()->SetRangeUser(0.0,0.50);
  h_v2->Draw("pE");
  g_v2->Draw("pE same");
  f_v2->Draw("l same");
  c_v2->SaveAs(Form("v2_%s_cent%s.pdf",vmsa::mBeamEnergy[energy].c_str(),centlabel.c_str()));

  return f_v2;
}

void tofMatchingFromTree(int energy = 4, bool datarcweight = false) 
{
    gStyle->SetHistFillColor(0);
    //gPad->SetAspectRatio(1);

    gStyle->SetOptStat(0);
 
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);    

    // INITIALIZE TTREES 
    //TFile *mKaonFileSE = TFile::Open("../data/Yields_Phi_SE_19GeV_20240830_kaontree.root");
    //TFile *mKaonFile = TFile::Open("../data/Yields_Phi_SE_19GeV_20240926_kaontree_ToFTag.root");
    TFile *mKaonFile = TFile::Open("/Volumes/File\ Backups/data/Yields_Phi_SE_19GeV_20240926_kaontree_ToFTag.root");
    //TFile *mKaonFile = TFile::Open("../data/Yields_Phi_SE_19GeV_20241219_cent2060_kaontree_ToFTag.root");
    TTree *mKaonTreeSE = (TTree*) mKaonFile->Get("kaontreeSE");
    int   mCentSE;
    float mWeightSE;
    int   mChargeSE;
    float mPtSE;
    float mRapiditySE;
    float mEtaSE;
    float mPhiSE;
    Bool_t mHasTofInfoSE;
    mKaonTreeSE->SetBranchAddress("cent", &mCentSE);
    mKaonTreeSE->SetBranchAddress("weight", &mWeightSE);
    mKaonTreeSE->SetBranchAddress("charge", &mChargeSE);
    mKaonTreeSE->SetBranchAddress("pt", &mPtSE);
    mKaonTreeSE->SetBranchAddress("rapidity", &mRapiditySE);
    mKaonTreeSE->SetBranchAddress("eta", &mEtaSE);
    mKaonTreeSE->SetBranchAddress("phi", &mPhiSE);
    mKaonTreeSE->SetBranchAddress("hastofinfo", &mHasTofInfoSE);

    TTree *mKaonTreeME = (TTree*) mKaonFile->Get("kaontreeME");
    int   mCentME;
    float mWeightME;
    int   mChargeME;
    float mPtME;
    float mRapidityME;
    float mEtaME;
    float mPhiME;
    Bool_t mHasTofInfoME;
    mKaonTreeME->SetBranchAddress("cent", &mCentME);
    mKaonTreeME->SetBranchAddress("weight", &mWeightME);
    mKaonTreeME->SetBranchAddress("charge", &mChargeME);
    mKaonTreeME->SetBranchAddress("pt", &mPtME);
    mKaonTreeME->SetBranchAddress("rapidity", &mRapidityME);
    mKaonTreeME->SetBranchAddress("eta", &mEtaME);
    mKaonTreeME->SetBranchAddress("phi", &mPhiME);
    mKaonTreeME->SetBranchAddress("hastofinfo", &mHasTofInfoME);
    // INITIALIZE TTREES 
   
     
    // INITIALIZE HISTOGRAMS
    //const int npt = 16;
    const int npt = 9;
    //const int npt = 6;
    const float ptmin = 0.0, ptmax = 4.0; 
                         //  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21
    float ptbins[npt+1] = {0.3,0.5,0.6,0.7,0.8,0.9,1.1,1.5,1.9,2.8};
    //float ptbins[npt+1] = {0.3,0.5,0.6,0.7,0.8,0.9,1.1};
                         //  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21
    //float ptbins[npt+1] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.5,3.5,5.0};
    //float ptbins[npt+1] = {0.3,0.4,.45,0.5,.55,0.6,.65,0.7,.75,0.8,.85,0.9,1.1,1.4,1.7,2.0,2.8};
    //const int ny = 30;
    //const float ymin = -1.5, ymax = 1.5;
    const int neta = 24;
    const float etamin = -1.0, etamax = 1.0;
    float etabins[neta+1] = {0.0};
    for(int ieta = 0; ieta < neta+1; ieta++)
    {
      etabins[ieta] = (float(ieta)-float(neta/2))/float(neta/2);
    }   

    const int nphi = 1;
    const float phimin = -TMath::Pi(), phimax = TMath::Pi();
    float phibins[nphi+1] = {phimin,phimax};
    //float phibins[nphi+1] = {0.0};
    //for(int iphi = 0; iphi < nphi+1; iphi++)
    //{
    //  phibins[iphi] = TMath::Pi()*(float(iphi)-float(nphi/2))/float(nphi/2);
    //}   

    TH1F*  pt[2][2][neta][nphi];
    TH1F* eta[2][2][npt][nphi];
    TH1F* phi[2][2][npt][neta];

    TH1F*  ratiopt[2][neta][nphi];
    TH1F* ratioeta[2][npt][nphi];
    TH1F* ratiophi[2][npt][neta];

    TH1F*  ptall[2][2];
    TH1F* etaall[2][2];
    TH1F* phiall[2][2];

    TH1F*  ratioptall[2];
    TH1F* ratioetaall[2];
    TH1F* ratiophiall[2];

    TH2F*  pteta[2][2];
    TH2F*  ptphi[2][2];
    TH2F* etaphi[2][2]; 

    TH2F*  ptphi_eta[2][2][neta];
    TH2F*  ratioptphi_eta[2][neta];
    TH2F*  ratio2ptphi_eta[2][neta];
    TH2F*  ratio2sigptphi_eta[2][neta];

    TH2F*  ratiopteta[2]; 
    TH2F*  ratioptphi[2];
    TH2F* ratioetaphi[2]; 

    //TH2F*  pteta[2][2][nphi];
    //TH2F*  ptphi[2][2][neta];
    //TH2F* etaphi[2][2][npt]; 

    //TH2F*  ratiopteta[2][nphi]; 
    //TH2F*  ratioptphi[2][neta];
    //TH2F* ratioetaphi[2][npt]; 

    TH3F* ptetaphi[2][3][2]; // [K+ probe, K- probe][SE, ME, SM][TPC, TPC+TOF]
    TH3F* ratioptetaphi[2];  // [K+ probe, K- probe]
 
    string charge[2]   = {"plus","minus"};
    string mixing[3]   = {"SE","ME","SM"};
    string detector[2] = {"TPC","TOF"};
 
    cout << "Before Set up histograms" << endl;

    for (int ic = 0; ic < 2; ic++)
    {
      for (int im = 0; im < 3; im++)
      {
        for(int idet = 0; idet < 2; idet++)
        {
          string hist;
          hist = Form("ptetaphi_k%s_%s_%s",charge[ic].c_str(),mixing[im].c_str(),detector[idet].c_str());
          //ptetaphi[ic][im][idet] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, neta, etamin, etamax, nphi, phimin, phimax);
          ptetaphi[ic][im][idet] = new TH3F(hist.c_str(), hist.c_str(), npt, ptbins, neta, etabins, nphi, phibins);
          ptetaphi[ic][im][idet]->Sumw2();
        }
      }
    }

    cout << "Set up histograms" << endl;

    // INITIALIZE HISTOGRAMS

    Long64_t nentries = mKaonTreeSE->GetEntries();
    cout << "nentries = " << nentries << endl;
    for (Long64_t i = 0; i < nentries; i++)
    {
      mKaonTreeSE->GetEntry(i);
      if (i%1000000 == 0)
      {
         cout << "Proccessed " << i << endl;
      }

      if(mCentSE > 5 || mCentSE < 2) continue;

      ptetaphi[mChargeSE][0][0]->Fill(mPtSE,mEtaSE,mPhiSE,mWeightSE);
      if(mHasTofInfoSE) ptetaphi[mChargeSE][0][1]->Fill(mPtSE,mEtaSE,mPhiSE,mWeightSE);

    }
      
    Long64_t nentriesME = mKaonTreeME->GetEntries();
    cout << "nentries = " << nentriesME << endl;
    for (Long64_t i = 0; i < nentriesME; i++)
    {
      mKaonTreeME->GetEntry(i);
      if (i%1000000 == 0)
      {
         cout << "Proccessed " << i << endl;
      }

      if(mCentME > 5 || mCentME < 2) continue;

      ptetaphi[mChargeME][1][0]->Fill(mPtME,mEtaME,mPhiME,mWeightME);
      if(mHasTofInfoME) ptetaphi[mChargeME][1][1]->Fill(mPtME,mEtaME,mPhiME,mWeightME);

    }
     
    for(int ic = 0; ic < 2; ic++)
    {
      for(int idet = 0; idet < 2; idet++)
      { 
        ptetaphi[ic][2][idet] = (TH3F*) ptetaphi[ic][0][idet]->Clone();
        ptetaphi[ic][2][idet]->Add(ptetaphi[ic][1][idet],-1.0);
        if(idet == 1)
        {
          ratioptetaphi[ic] = (TH3F*) ptetaphi[ic][2][1]->Clone();
          ratioptetaphi[ic]->Divide(ptetaphi[ic][2][1],ptetaphi[ic][2][0],1,1,"B");
        }

        for(int ieta = 0; ieta < neta; ieta++)
        {
          for(int iphi = 0; iphi < nphi; iphi++) 
          {
            pt[ic][idet][ieta][iphi] = (TH1F*) ptetaphi[ic][2][idet]->ProjectionX(Form("pt_%d_%d_%d_%d",ic,idet,ieta,iphi),ieta+1,ieta+1,iphi+1,iphi+1);
            if(idet == 1)
            {
              ratiopt[ic][ieta][iphi] = (TH1F*) pt[ic][1][ieta][iphi]->Clone();
              ratiopt[ic][ieta][iphi]->Divide(pt[ic][1][ieta][iphi],pt[ic][0][ieta][iphi],1,1,"B");
            }
          }
        }
        for(int ipt = 0; ipt < npt; ipt++)
        {
          for(int iphi = 0; iphi < nphi; iphi++) 
          {
            eta[ic][idet][ipt][iphi] = (TH1F*) ptetaphi[ic][2][idet]->ProjectionY(Form("eta_%d_%d_%d_%d",ic,idet,ipt,iphi),ipt+1,ipt+1,iphi+1,iphi+1);
            if(idet == 1)
            {
              ratioeta[ic][ipt][iphi] = (TH1F*) eta[ic][1][ipt][iphi]->Clone();
              ratioeta[ic][ipt][iphi]->Divide(eta[ic][1][ipt][iphi],eta[ic][0][ipt][iphi],1,1,"B");
            }
          }
        }
        for(int ipt = 0; ipt < npt; ipt++)
        {
          for(int ieta = 0; ieta < neta; ieta++) 
          {
            phi[ic][idet][ipt][ieta] = (TH1F*) ptetaphi[ic][2][idet]->ProjectionZ(Form("phi_%d_%d_%d_%d",ic,idet,ipt,ieta),ipt+1,ipt+1,ieta+1,ieta+1);
            if(idet == 1)
            {
              ratiophi[ic][ipt][ieta] = (TH1F*) phi[ic][1][ipt][ieta]->Clone();
              ratiophi[ic][ipt][ieta]->Divide(phi[ic][0][ipt][ieta],phi[ic][0][ipt][ieta],1,1,"B");
            }
          }
        }
        for(int ieta = 0; ieta < neta; ieta++) 
        {
          string nameptphi = Form("ptphi_%d_%d_%d",ic,idet,ieta);
          
          ptetaphi[ic][2][idet]->GetYaxis()->SetRange(ieta+1,ieta+1);
          ptphi_eta[ic][idet][ieta] = (TH2F*) ((TH2F*) ptetaphi[ic][2][idet]->Project3D("zx"))->Clone(nameptphi.c_str());
          //ptphi_eta[ic][idet][ieta] = (TH2F*) ptetaphi[ic][2][idet]->ProjectionXZ(nameptphi.c_str(),ieta+1,ieta+1);
          if(idet == 1)
          {
            ratioptphi_eta[ic][ieta] = (TH2F*) ptphi_eta[ic][1][ieta]->Clone();
            //ratioptphi_eta[ic][ieta]->Divide(ptphi_eta[ic][0][ieta],ptphi_eta[ic][0][ieta],1,1,"B");
            ratioptphi_eta[ic][ieta]->Divide(ptphi_eta[ic][0][ieta]);
          }
          ptetaphi[ic][2][idet]->GetYaxis()->SetRange();
        }
        for(int ieta = 0; ieta < neta; ieta++) 
        {
          if(idet == 1)
          {
            string nameptphi = Form("ratio_ptphi_%d_%d_%d",ic,idet,ieta);
            ratio2ptphi_eta[ic][ieta] = (TH2F*) ratioptphi_eta[ic][6]->Clone(nameptphi.c_str());
            ratio2ptphi_eta[ic][ieta]->Divide(ratioptphi_eta[ic][ieta]);

            nameptphi = Form("ratiosig_ptphi_%d_%d_%d",ic,idet,ieta);
            ratio2sigptphi_eta[ic][ieta] = (TH2F*) ratioptphi_eta[ic][6]->Clone(nameptphi.c_str());
            for(int ix = 1; ix <= ratio2ptphi_eta[ic][ieta]->GetNbinsX(); ix++)
            {
              for(int iy = 1; iy <= ratio2ptphi_eta[ic][ieta]->GetNbinsY(); iy++)
              {
                int global_bin = ratio2ptphi_eta[ic][ieta]->GetBin(ix,iy);   
                double val = ratio2ptphi_eta[ic][ieta]->GetBinContent(global_bin);   
                double err = ratio2ptphi_eta[ic][ieta]->GetBinError(global_bin);   
                
                double sig = 1/err;

                ratio2sigptphi_eta[ic][ieta]->SetBinContent(global_bin,sig);
              }
            }
          }
        }

        //ptall[ic][idet]  = (TH1F*) ptetaphi[ic][2][idet]->Project3D("x");
        //etaall[ic][idet] = (TH1F*) ptetaphi[ic][2][idet]->Project3D("y");
        //phiall[ic][idet] = (TH1F*) ptetaphi[ic][2][idet]->Project3D("z");

        //pteta[ic][idet]  = (TH2F*) ptetaphi[ic][2][idet]->Project3D("yx");
        //ptphi[ic][idet]  = (TH2F*) ptetaphi[ic][2][idet]->Project3D("zx");
        //etaphi[ic][idet] = (TH2F*) ptetaphi[ic][2][idet]->Project3D("zy");

      }

      //ratioptall[ic] = (TH1F*) ptall[ic][1]->Clone();
      //ratioptall[ic]->Divide(ptall[ic][1],ptall[ic][0],1,1,"B");

      //ratioetaall[ic] = (TH1F*) etaall[ic][1]->Clone();
      //ratioetaall[ic]->Divide(etaall[ic][1],etaall[ic][0],1,1,"B");

      //ratiophiall[ic] = (TH1F*) phiall[ic][1]->Clone();
      //ratiophiall[ic]->Divide(phiall[ic][1],phiall[ic][0],1,1,"B");

      //ratiopteta[ic] = (TH2F*) pteta[ic][1]->Clone();
      //ratiopteta[ic]->Divide(pteta[ic][1],pteta[ic][0],1,1,"B");

      //ratioptphi[ic] = (TH2F*) ptphi[ic][1]->Clone();
      //ratioptphi[ic]->Divide(ptphi[ic][1],ptphi[ic][0],1,1,"B");

      //ratioetaphi[ic] = (TH2F*) etaphi[ic][1]->Clone();
      //ratioetaphi[ic]->Divide(etaphi[ic][1],etaphi[ic][0],1,1,"B");
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->Divide(3,2);
    for(int i = 0; i < 6; i++)
    {
      c1->cd(i+1);
      c1->cd(i+1)->SetLeftMargin(0.12);
      c1->cd(i+1)->SetRightMargin(0.15);
      c1->cd(i+1)->SetBottomMargin(0.12);
      c1->cd(i+1)->SetTicks(1,1);
      c1->cd(i+1)->SetGrid(0,0); 
    }

    string mixinglabel[3] = {"SE","ME","SEminusME"};
    string mixingtitle[3] = {"Same Event","Mixed Event","Same - Mixed Event"};

    string outputname  = Form("figures/ToFMatching/Kaons_ToF_pt.pdf");
    string outputstart = Form("%s[",outputname.c_str());
    string outputstop  = Form("%s]",outputname.c_str());

    c1->Print(outputstart.c_str());

    for(int ieta = 0; ieta < neta; ieta++)
    {
      for(int iphi = 0; iphi < nphi; iphi++)
      {
        for(int ic = 0; ic < 2; ic++)
        {
          string histtitle;

          histtitle = Form("K%s, TPC, %1.2f<#eta<%1.2f, %1.2f<#phi<%1.2f",charge[ic].c_str(),etabins[ieta],etabins[ieta+1],phibins[iphi],phibins[iphi+1]);
          c1->cd(1+3*ic);
          c1->cd(1+3*ic)->SetLogy();
          pt[ic][0][ieta][iphi]->SetTitle(histtitle.c_str()); 
          pt[ic][0][ieta][iphi]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
          pt[ic][0][ieta][iphi]->GetYaxis()->SetTitle("Count"); 
          pt[ic][0][ieta][iphi]->Draw("pE"); 

          histtitle = Form("K%s, TPC&ToF, %1.2f<#eta<%1.2f, %1.2f<#phi<%1.2f",charge[ic].c_str(),etabins[ieta],etabins[ieta+1],phibins[iphi],phibins[iphi+1]);
          c1->cd(2+3*ic);
          c1->cd(2+3*ic)->SetLogy();
          pt[ic][1][ieta][iphi]->SetTitle(histtitle.c_str()); 
          pt[ic][1][ieta][iphi]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
          pt[ic][1][ieta][iphi]->GetYaxis()->SetTitle("Count"); 
          pt[ic][1][ieta][iphi]->Draw("pE"); 

          histtitle = Form("K%s, Eff., %1.2f<#eta<%1.2f, %1.2f<#phi<%1.2f",charge[ic].c_str(),etabins[ieta],etabins[ieta+1],phibins[iphi],phibins[iphi+1]);
          c1->cd(3+3*ic);
          ratiopt[ic][ieta][iphi]->SetTitle(histtitle.c_str()); 
          ratiopt[ic][ieta][iphi]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
          ratiopt[ic][ieta][iphi]->GetYaxis()->SetTitle("ToF Matching Efficiency"); 
          ratiopt[ic][ieta][iphi]->GetYaxis()->SetRangeUser(-2.0,3.0); 
          ratiopt[ic][ieta][iphi]->Draw("pE"); 

        }
        c1->Update();
        c1->Print(outputname.c_str());
      }
    }
    c1->Print(outputstop.c_str());

    outputname  = Form("figures/ToFMatching/Kaons_ToF_eta.pdf");
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());

    c1->Print(outputstart.c_str());

    for(int ipt = 0; ipt < npt; ipt++)
    {
      for(int iphi = 0; iphi < nphi; iphi++)
      {
        for(int ic = 0; ic < 2; ic++)
        {
          string histtitle;

          histtitle = Form("K%s, TPC, %1.2f<p_{T}<%1.2f, %1.2f<#phi<%1.2f",charge[ic].c_str(),ptbins[ipt],ptbins[ipt+1],phibins[iphi],phibins[iphi+1]);
          c1->cd(1+3*ic);
          c1->cd(1+3*ic)->SetLogy(0);
          eta[ic][0][ipt][iphi]->SetTitle(histtitle.c_str()); 
          eta[ic][0][ipt][iphi]->GetXaxis()->SetTitle("#eta"); 
          eta[ic][0][ipt][iphi]->GetYaxis()->SetTitle("Count"); 
          eta[ic][0][ipt][iphi]->Draw("pE"); 

          histtitle = Form("K%s, TPC&ToF, %1.2f<p_{T}<%1.2f, %1.2f<#phi<%1.2f",charge[ic].c_str(),ptbins[ipt],ptbins[ipt+1],phibins[iphi],phibins[iphi+1]);
          c1->cd(2+3*ic);
          c1->cd(2+3*ic)->SetLogy(0);
          eta[ic][1][ipt][iphi]->SetTitle(histtitle.c_str()); 
          eta[ic][1][ipt][iphi]->GetXaxis()->SetTitle("#eta"); 
          eta[ic][1][ipt][iphi]->GetYaxis()->SetTitle("Count"); 
          eta[ic][1][ipt][iphi]->Draw("pE"); 

          histtitle = Form("K%s, Eff., %1.2f<p_{T}<%1.2f, %1.2f<#phi<%1.2f",charge[ic].c_str(),ptbins[ipt],ptbins[ipt+1],phibins[iphi],phibins[iphi+1]);
          c1->cd(3+3*ic);
          ratioeta[ic][ipt][iphi]->SetTitle(histtitle.c_str()); 
          ratioeta[ic][ipt][iphi]->GetXaxis()->SetTitle("#eta"); 
          ratioeta[ic][ipt][iphi]->GetYaxis()->SetTitle("ToF Matching Efficiency"); 
          ratioeta[ic][ipt][iphi]->GetYaxis()->SetRangeUser(0.0,1.2); 
          ratioeta[ic][ipt][iphi]->Draw("pE"); 

        }
        c1->Update();
        c1->Print(outputname.c_str());
      }
    }
    c1->Print(outputstop.c_str());

    outputname  = Form("figures/ToFMatching/Kaons_ToF_phi.pdf");
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());

    c1->Print(outputstart.c_str());

    for(int ipt = 0; ipt < npt; ipt++)
    {
      for(int ieta = 0; ieta < neta; ieta++)
      {
        for(int ic = 0; ic < 2; ic++)
        {
          string histtitle;

          histtitle = Form("K%s, TPC, %1.2f<p_{T}<%1.2f, %1.2f<#eta<%1.2f",charge[ic].c_str(),ptbins[ipt],ptbins[ipt+1],etabins[ieta],etabins[ieta+1]);
          c1->cd(1+3*ic);
          phi[ic][0][ipt][ieta]->SetTitle(histtitle.c_str()); 
          phi[ic][0][ipt][ieta]->GetXaxis()->SetTitle("#phi"); 
          phi[ic][0][ipt][ieta]->GetYaxis()->SetTitle("Count"); 
          phi[ic][0][ipt][ieta]->Draw("pE"); 

          histtitle = Form("K%s, TPC&ToF, %1.2f<p_{T}<%1.2f, %1.2f<#eta<%1.2f",charge[ic].c_str(),ptbins[ipt],ptbins[ipt+1],etabins[ieta],etabins[ieta+1]);
          c1->cd(2+3*ic);
          phi[ic][1][ipt][ieta]->SetTitle(histtitle.c_str()); 
          phi[ic][1][ipt][ieta]->GetXaxis()->SetTitle("#phi"); 
          phi[ic][1][ipt][ieta]->GetYaxis()->SetTitle("Count"); 
          phi[ic][1][ipt][ieta]->Draw("pE"); 

          histtitle = Form("K%s, Eff., %1.2f<p_{T}<%1.2f, %1.2f<#eta<%1.2f",charge[ic].c_str(),ptbins[ipt],ptbins[ipt+1],etabins[ieta],etabins[ieta+1]);
          c1->cd(3+3*ic);
          ratiophi[ic][ipt][ieta]->SetTitle(histtitle.c_str()); 
          ratiophi[ic][ipt][ieta]->GetXaxis()->SetTitle("#phi"); 
          ratiophi[ic][ipt][ieta]->GetYaxis()->SetTitle("ToF Matching Efficiency"); 
          ratiophi[ic][ipt][ieta]->GetYaxis()->SetRangeUser(-2.0,3.0); 
          ratiophi[ic][ipt][ieta]->Draw("pE"); 

        }
        c1->Update();
        c1->Print(outputname.c_str());
      }
    }
    c1->Print(outputstop.c_str());

    outputname  = Form("figures/ToFMatching/Kaons_ToF_ptphi.pdf");
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());

    c1->Print(outputstart.c_str());

    for(int ieta = 0; ieta < neta; ieta++)
    {
      for(int ic = 0; ic < 2; ic++)
      {
        string histtitle;

        histtitle = Form("K%s, TPC, %1.2f<#eta<%1.2f",charge[ic].c_str(),etabins[ieta],etabins[ieta+1]);
        c1->cd(1+3*ic);
        ptphi_eta[ic][0][ieta]->SetTitle(histtitle.c_str()); 
        ptphi_eta[ic][0][ieta]->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
        ptphi_eta[ic][0][ieta]->GetYaxis()->SetTitle("#phi"); 
        ptphi_eta[ic][0][ieta]->Draw("colz"); 

        histtitle = Form("K%s, TPC&ToF, %1.2f<#eta<%1.2f",charge[ic].c_str(),etabins[ieta],etabins[ieta+1]);
        c1->cd(2+3*ic);
        ptphi_eta[ic][1][ieta]->SetTitle(histtitle.c_str()); 
        ptphi_eta[ic][1][ieta]->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
        ptphi_eta[ic][1][ieta]->GetYaxis()->SetTitle("#phi"); 
        ptphi_eta[ic][1][ieta]->Draw("colz"); 

        histtitle = Form("K%s, Eff., %1.2f<#eta<%1.2f",charge[ic].c_str(),etabins[ieta],etabins[ieta+1]);
        c1->cd(3+3*ic);
        ratioptphi_eta[ic][ieta]->SetTitle(histtitle.c_str()); 
        ratioptphi_eta[ic][ieta]->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
        ratioptphi_eta[ic][ieta]->GetYaxis()->SetTitle("#phi"); 
        ratioptphi_eta[ic][ieta]->Draw("colz"); 

      }
      c1->Update();
      c1->Print(outputname.c_str());
    }
    
    c1->Print(outputstop.c_str());

    
    outputname  = Form("figures/ToFMatching/Kaons_ToF_ptphi_binbybinratio.pdf");
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());

    TCanvas *c2 = new TCanvas("c2", "c2", 2400,1600);
    c2->Divide(6,4);
    for(int i = 0; i < 24; i++)
    {
      c2->cd(i+1);
      c2->cd(i+1)->SetLeftMargin(0.12);
      c2->cd(i+1)->SetRightMargin(0.15);
      c2->cd(i+1)->SetBottomMargin(0.12);
      c2->cd(i+1)->SetTicks(1,1);
      c2->cd(i+1)->SetGrid(0,0); 
    }

    c2->Print(outputstart.c_str());

    for(int ic = 0; ic < 2; ic++)
    {
      for(int ieta = 0; ieta < neta; ieta++)
      {
        string histtitle;

        histtitle = Form("K%s, Eff. Ratio, %1.2f<#eta<%1.2f",charge[ic].c_str(),etabins[ieta],etabins[ieta+1]);
        c2->cd(1+ieta);
        ratio2ptphi_eta[ic][ieta]->SetTitle(histtitle.c_str()); 
        ratio2ptphi_eta[ic][ieta]->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
        ratio2ptphi_eta[ic][ieta]->GetYaxis()->SetTitle("#phi"); 
        //ratio2ptphi_eta[ic][ieta]->GetZaxis()->SetRangeUser(0.4,5.5); 
        ratio2ptphi_eta[ic][ieta]->Draw("colz"); 

        histtitle = Form("K%s, Eff. Ratio Significance, %1.2f<#eta<%1.2f",charge[ic].c_str(),etabins[ieta],etabins[ieta+1]);
        c2->cd(13+ieta);
        ratio2sigptphi_eta[ic][ieta]->SetTitle(histtitle.c_str()); 
        ratio2sigptphi_eta[ic][ieta]->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
        ratio2sigptphi_eta[ic][ieta]->GetYaxis()->SetTitle("#phi"); 
        //ratio2sigptphi_eta[ic][ieta]->GetZaxis()->SetRangeUser(0.4,5.5); 
        ratio2sigptphi_eta[ic][ieta]->Draw("colz"); 

      }
      c2->Update();
      c2->Print(outputname.c_str());
    }
    
    c2->Print(outputstop.c_str());

    
    TFile *OutFile = new TFile("ToFMatching/ToFMatching_19GeV.root", "RECREATE");
    OutFile->cd();
    for(int ic = 0; ic < 2; ic++)
    {
      string histname = Form("k%s_ratio",charge[ic].c_str());
      cout << histname << endl;
      ratioptetaphi[ic]->SetName(histname.c_str());
      ratioptetaphi[ic]->Write();
    }
    OutFile->Close();

        
}

