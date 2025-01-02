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
   
  if(gRandom) delete gRandom;
  gRandom = new TRandom3();
  gRandom->SetSeed();
 
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

void plotKaonTTrees_Compare_SingleKaonBool(int energy = 4, bool datarcweight = false) 
{
    gStyle->SetHistFillColor(0);
    //gPad->SetAspectRatio(1);

    string folderopt = "SingleKaonBool/";
    if(datarcweight) folderopt = "ToFFromTaggingFinerPID/PhiDataRCWeighted_YRatioGaus/";

    gStyle->SetOptStat(0);
 
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);    

    TFile *inputweightsptyphi = TFile::Open(Form("phiweights/pt_y_phi_weights_v2_FunctionEmbed.root"));
    //TFile *inputweightsptyphi = TFile::Open(Form("phiweights/pt_y_phi_weights_v2_FunctionEmbed_ygaus.root"));

    TH3F *h_mpTyv2Weights[9];
    for(int i_cent = 2; i_cent < 6; ++i_cent)
    {
      std::string HistName = Form("h_ratio_cent%d",i_cent);
      h_mpTyv2Weights[i_cent] = (TH3F*) inputweightsptyphi->Get(HistName.c_str());
    }

    TFile *ToFFile = TFile::Open(Form("ToFMatching/ToFMatching_19GeV.root"));
    TH3F* ToFHist[2];
    ToFHist[0] = (TH3F*) ToFFile->Get("kplus_ratio");
    ToFHist[0]->Print();
    ToFHist[1] = (TH3F*) ToFFile->Get("kminus_ratio");
    ToFHist[1]->Print();


    // INITIALIZE TTREES 
    //TFile *mKaonFileSE = TFile::Open("../data/Yields_Phi_SE_19GeV_20240830_kaontree.root");
    TFile *mKaonFileSE = TFile::Open("../data/Yields_Phi_SE_19GeV_20240926_kaontree.root");
    TTree *mKaonTreeSE = (TTree*) mKaonFileSE->Get("kaontree");
    int   mCentSE;
    float mWeightSE;
    int   mChargeSE;
    float mPtSE;
    float mRapiditySE;
    float mEtaSE;
    float mPhiSE;
    int   mNHitsFitSE;
    int   mNHitsMaxSE;
    float mDEdxSE;
    float mDcaSE;
    mKaonTreeSE->SetBranchAddress("cent", &mCentSE);
    mKaonTreeSE->SetBranchAddress("weight", &mWeightSE);
    mKaonTreeSE->SetBranchAddress("charge", &mChargeSE);
    mKaonTreeSE->SetBranchAddress("pt", &mPtSE);
    mKaonTreeSE->SetBranchAddress("rapidity", &mRapiditySE);
    mKaonTreeSE->SetBranchAddress("eta", &mEtaSE);
    mKaonTreeSE->SetBranchAddress("phi", &mPhiSE);
    mKaonTreeSE->SetBranchAddress("nhitsfit", &mNHitsFitSE);
    mKaonTreeSE->SetBranchAddress("nhitsmax", &mNHitsMaxSE);
    mKaonTreeSE->SetBranchAddress("dedx", &mDEdxSE);
    mKaonTreeSE->SetBranchAddress("dca", &mDcaSE);

    //TFile *mKaonFileME = TFile::Open("../data/Yields_Phi_ME_19GeV_20240830_kaontree.root");
    TFile *mKaonFileME = TFile::Open("../data/Yields_Phi_ME_19GeV_20240926_kaontree.root");
    TTree *mKaonTreeME = (TTree*) mKaonFileME->Get("kaontree");
    int   mCentME;
    float mWeightME;
    int   mChargeME;
    float mPtME;
    float mRapidityME;
    float mEtaME;
    float mPhiME;
    int   mNHitsFitME;
    int   mNHitsMaxME;
    float mDEdxME;
    float mDcaME;
    mKaonTreeME->SetBranchAddress("cent", &mCentME);
    mKaonTreeME->SetBranchAddress("weight", &mWeightME);
    mKaonTreeME->SetBranchAddress("charge", &mChargeME);
    mKaonTreeME->SetBranchAddress("pt", &mPtME);
    mKaonTreeME->SetBranchAddress("rapidity", &mRapidityME);
    mKaonTreeME->SetBranchAddress("eta", &mEtaME);
    mKaonTreeME->SetBranchAddress("phi", &mPhiME);
    mKaonTreeME->SetBranchAddress("nhitsfit", &mNHitsFitME);
    mKaonTreeME->SetBranchAddress("nhitsmax", &mNHitsMaxME);
    mKaonTreeME->SetBranchAddress("dedx", &mDEdxME);
    mKaonTreeME->SetBranchAddress("dca", &mDcaME);
    // INITIALIZE TTREES 
   
     
    // INITIALIZE TTREES 
    int mCent;
    float mEpFull;
    int mNToFMatch;
    float mWeight;

    TLorentzVector *mLMeson = nullptr;
    TLorentzVector *mLKp = nullptr;
    TLorentzVector *mLKm = nullptr; 
    Bool_t mRcExists;
    TLorentzVector *mLRcKp = nullptr;
    TLorentzVector *mLRcKm = nullptr;

    int mNHitsFitP;
    int mNHitsMaxP;
    float mDEdxP;
    float mDcaP;
    int mNHitsFitM;
    int mNHitsMaxM;
    float mDEdxM;
    float mDcaM;

    float mMcPhiStar;
    float mMcCosThetaStar;
    float mMcPhiPrime;
    float mMcCosTheta;
    float mMcHelicityAngle;
    
    float mRcPhiStar;
    float mRcCosThetaStar;
    float mRcPhiPrime;
    float mRcCosTheta;
    float mRcHelicityAngle;
    
    Bool_t mPassTpcKp;
    Bool_t mPassTofKp;
    Bool_t mPassM2Kp;
    Bool_t mPassNsigKp;

    Bool_t mPassTpcKm;
    Bool_t mPassTofKm;
    Bool_t mPassM2Km;
    Bool_t mPassNsigKm;


    //TFile *mKaonFileEmbed = TFile::Open("effaccfiles/Phi/19GeV/PhiEmbedding_kaontrees/PhiEmbedding_kaontrees_fixed.root");
    //TFile *mKaonFileEmbed = TFile::Open("effaccfiles/Phi/19GeV/PhiEmbedding_kaontrees/PhiEmbedding_kaontrees_notofmatchreq_EP.root");
    //TFile *mKaonFileEmbed = TFile::Open("effaccfiles/Phi/19GeV/PhiEmbedding_kaontrees/PhiEmbedding_ToFMatching_48PhiBins.root");
    TFile *mKaonFileEmbed = TFile::Open("effaccfiles/Phi/19GeV/PhiEmbedding_kaontrees_singlakaonvariables.root");
    TTree *mKaonTreeEmbed = (TTree*) mKaonFileEmbed->Get("kaontree");
  
    mKaonTreeEmbed->Print();

    mKaonTreeEmbed->SetBranchAddress("cent", &mCent);
    mKaonTreeEmbed->SetBranchAddress("epfull", &mEpFull);
    mKaonTreeEmbed->SetBranchAddress("ntofmatch", &mNToFMatch);
    mKaonTreeEmbed->SetBranchAddress("weight", &mWeight);

    mKaonTreeEmbed->SetBranchAddress("lmeson", &mLMeson);
    mKaonTreeEmbed->SetBranchAddress("lkp", &mLKp);
    mKaonTreeEmbed->SetBranchAddress("lkm", &mLKm);
    mKaonTreeEmbed->SetBranchAddress("rcexists", &mRcExists);
    mKaonTreeEmbed->SetBranchAddress("rlkp", &mLRcKp);
    mKaonTreeEmbed->SetBranchAddress("rlkm", &mLRcKm);

    mKaonTreeEmbed->SetBranchAddress("nhitsfitp", &mNHitsFitP);
    mKaonTreeEmbed->SetBranchAddress("nhitsmaxp", &mNHitsMaxP);
    mKaonTreeEmbed->SetBranchAddress("dedxp", &mDEdxP);
    mKaonTreeEmbed->SetBranchAddress("dcap", &mDcaP);
    mKaonTreeEmbed->SetBranchAddress("nhitsfitm", &mNHitsFitM);
    mKaonTreeEmbed->SetBranchAddress("nhitsmaxm", &mNHitsMaxM);
    mKaonTreeEmbed->SetBranchAddress("dedxm", &mDEdxM);
    mKaonTreeEmbed->SetBranchAddress("dcam", &mDcaM);

    mKaonTreeEmbed->SetBranchAddress("mcphistar", &mMcPhiStar);
    mKaonTreeEmbed->SetBranchAddress("mccosthetastar", &mMcCosThetaStar);
    mKaonTreeEmbed->SetBranchAddress("mcphiprime", &mMcPhiPrime);
    mKaonTreeEmbed->SetBranchAddress("mccostheta", &mMcCosTheta);
    mKaonTreeEmbed->SetBranchAddress("mchelicityangle", &mMcHelicityAngle);
    
    mKaonTreeEmbed->SetBranchAddress("rcphistar", &mRcPhiStar);
    mKaonTreeEmbed->SetBranchAddress("rccosthetastar", &mRcCosThetaStar);
    mKaonTreeEmbed->SetBranchAddress("rcphiprime", &mRcPhiPrime);
    mKaonTreeEmbed->SetBranchAddress("rccostheta", &mRcCosTheta);
    mKaonTreeEmbed->SetBranchAddress("rchelicityangle", &mRcHelicityAngle);
    
    mKaonTreeEmbed->SetBranchAddress("passtpckp", &mPassTpcKp);
    mKaonTreeEmbed->SetBranchAddress("passtofkp", &mPassTofKp);
    mKaonTreeEmbed->SetBranchAddress("passm2kp", &mPassM2Kp);
    mKaonTreeEmbed->SetBranchAddress("passnskp", &mPassNsigKp);

    mKaonTreeEmbed->SetBranchAddress("passtpckm", &mPassTpcKm);
    mKaonTreeEmbed->SetBranchAddress("passtofkm", &mPassTofKm);
    mKaonTreeEmbed->SetBranchAddress("passm2km", &mPassM2Km);
    mKaonTreeEmbed->SetBranchAddress("passnskm", &mPassNsigKm);
    // INITIALIZE TTREES 

    TH2F* pt_correction_helicity[vmsa::pt_rebin_2D][2];
    TH2F* pt_correction_global[vmsa::pt_rebin_2D][2];
    TH2F* cent_correction_helicity[10][vmsa::pt_rebin_cent_2D][2];
    TH2F* cent_correction_global[10][vmsa::pt_rebin_cent_2D][2];
    TH2F* y_correction_helicity[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total][2];
    TH2F* y_correction_global[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total][2];
     
    for(int ipt = vmsa::pt_rebin_first_2D[energy]; ipt < vmsa::pt_rebin_last_2D[energy]; ipt++)
    {
      string histname;
      histname = Form("h_Global_MC_Cent_9_Pt_%d",ipt);
      pt_correction_global[ipt][0] = new TH2F(histname.c_str(),histname.c_str(),9,-1,1,12,0.0,2.0*TMath::Pi());
      histname = Form("h_Global_RC_Cent_9_Pt_%d",ipt);
      pt_correction_global[ipt][1] = new TH2F(histname.c_str(),histname.c_str(),9,-1,1,12,0.0,2.0*TMath::Pi());

      histname = Form("h_Helicity_MC_Cent_9_Pt_%d",ipt);
      pt_correction_helicity[ipt][0] = new TH2F(histname.c_str(),histname.c_str(),9,-1,1,12,0.0,2.0*TMath::Pi());
      histname = Form("h_Helicity_RC_Cent_9_Pt_%d",ipt);
      pt_correction_helicity[ipt][1] = new TH2F(histname.c_str(),histname.c_str(),9,-1,1,12,0.0,2.0*TMath::Pi());
    }



    // INITIALIZE HISTOGRAMS
    TH1F* pt[5][3];
    TH1F* rapidity[5][3]; 
    TH1F* eta[5][3];
    TH1F* phi[5][3];

    TH1F* spt[5][3][20];
    TH1F* srapidity[5][3][20]; 
    TH1F* seta[5][3][20];
    TH1F* sphi[5][3][20];

    //THnF* ptyetaphi[5][3];

    TH1F* ratiopt[3][2];
    TH1F* ratiorapidity[3][2]; 
    TH1F* ratioeta[3][2];
    TH1F* ratiophi[3][2];

    TH1F* sratiopt[3][2][20];
    TH1F* sratiorapidity[3][2][20]; 
    TH1F* sratioeta[3][2][20];
    TH1F* sratiophi[3][2][20];
     
    TH2F* ptrapidity[5][3];
    TH2F* pteta[5][3]; 
    TH2F* ptphi[5][3];
    TH2F* rapidityeta[5][3];
    TH2F* rapidityphi[5][3];
    TH2F* etaphi[5][3]; 

    TH2F* sptrapidity[5][3][20];
    TH2F* spteta[5][3][20]; 
    TH2F* sptphi[5][3][20];
    TH2F* srapidityeta[5][3][20];
    TH2F* srapidityphi[5][3][20];
    TH2F* setaphi[5][3][20]; 

    TH2F* ratioptrapidity[3][2];
    TH2F* ratiopteta[3][2]; 
    TH2F* ratioptphi[3][2];
    TH2F* ratiorapidityeta[3][2];
    TH2F* ratiorapidityphi[3][2];
    TH2F* ratioetaphi[3][2]; 

    TH2F* sratioptrapidity[3][2][20];
    TH2F* sratiopteta[3][2][20]; 
    TH2F* sratioptphi[3][2][20];
    TH2F* sratiorapidityeta[3][2][20];
    TH2F* sratiorapidityphi[3][2][20];
    TH2F* sratioetaphi[3][2][20]; 

    TH2F* kaonptcos[2][5]; 
    TH2F* kaonptcosratio[2][5]; 
  
    // Tracking information
    THnF* ptyetaphibasic[5][3];
    //THnF* ptyetaphia[3][2];
    //THnF* ptyetaphib[3][2];
    //THnF* ptyetaphic[3][2];
    //THnF* ptyetaphid[3][2];
    //THnF* ptyetaphie[3][2];

    //TH1F* nhitsfit[4][2]; 
    //TH2F* nhitsfitpt[4][2]; 
    //TH2F* nhitsfitrapidity[4][2]; 
    //TH2F* nhitsfiteta[4][2]; 
    //TH2F* nhitsfitphi[4][2]; 

    //TH1F* nhitsmax[4][2]; 
    //TH2F* nhitsmaxpt[4][2]; 
    //TH2F* nhitsmaxrapidity[4][2]; 
    //TH2F* nhitsmaxeta[4][2]; 
    //TH2F* nhitsmaxphi[4][2]; 

    //TH1F* nhitsratio[4][2];
    //TH2F* nhitsratiopt[4][2];
    //TH2F* nhitsratiorapidity[4][2];
    //TH2F* nhitsratioeta[4][2];
    //TH2F* nhitsratiophi[4][2];

    //TH1F* dedx[4][2];
    //TH2F* dedxpt[4][2];
    //TH2F* dedxrapidity[4][2];
    //TH2F* dedxeta[4][2];
    //TH2F* dedxphi[4][2];

    //TH1F* dca[4][2];
    //TH2F* dcapt[4][2];
    //TH2F* dcarapidity[4][2];
    //TH2F* dcaeta[4][2];
    //TH2F* dcaphi[4][2];
     
    //TH1F* rationhitsfit[2]; 
    //TH2F* rationhitsfitpt[2]; 
    //TH2F* rationhitsfitrapidity[2]; 
    //TH2F* rationhitsfiteta[2]; 
    //TH2F* rationhitsfitphi[2]; 

    //TH1F* rationhitsmax[2]; 
    //TH2F* rationhitsmaxpt[2]; 
    //TH2F* rationhitsmaxrapidity[2]; 
    //TH2F* rationhitsmaxeta[2]; 
    //TH2F* rationhitsmaxphi[2]; 

    //TH1F* rationhitsratio[2];
    //TH2F* rationhitsratiopt[2];
    //TH2F* rationhitsratiorapidity[2];
    //TH2F* rationhitsratioeta[2];
    //TH2F* rationhitsratiophi[2];

    //TH1F* ratiodedx[2];
    //TH2F* ratiodedxpt[2];
    //TH2F* ratiodedxrapidity[2];
    //TH2F* ratiodedxeta[2];
    //TH2F* ratiodedxphi[2];

    //TH1F* ratiodca[2];
    //TH2F* ratiodcapt[2];
    //TH2F* ratiodcarapidity[2];
    //TH2F* ratiodcaeta[2];
    //TH2F* ratiodcaphi[2];
     
    //TH3F* ptrapidityphi[9];
    //TH3F* ptrapidityphiDesired[9];
    //TH3F* ptrapidityphiRatio[9];
 
    const int npt = 30;
    const float ptmin = 0.0, ptmax = 3.0; 
    const int ny = 50;
    const float ymin = -1.5, ymax = 1.5;
    const int neta = 50;
    const float etamin = -1.5, etamax = 1.5;
    const int nphi = 24;
    const float phimin = -TMath::Pi(), phimax = TMath::Pi();

    Int_t bins[4] = {npt, ny, neta, nphi};
    Double_t xmin[4] = {ptmin, ymin, etamin, phimin};
    Double_t xmax[4] = {ptmax, ymax, etamax, phimax};

    const int nhf = 50;
    const float nhfmin = -0.5, nhfmax = 100.5; 
    const int nhm = 50;
    const float nhmmin = -0.5, nhmmax = 100.5; 
    const int nhr = 20;
    const float nhrmin = 0.5, nhrmax = 1.1; 
    const int nde = 75;
    const float demin = 0.0, demax = 10.0; 
    const int ndca = 75;
    const float dcamin = 0.0, dcamax = 2.1; 

    //Int_t binsb[9] = {npt, ny, neta, nphi, nhf, nhm, nhr, nde, ndca};
    //Double_t xminb[9] = {ptmin, ymin, etamin, phimin, nhfmin, nhmmin, nhrmin, demin, dcamin};
    //Double_t xmaxb[9] = {ptmax, ymax, etamax, phimax, nhfmax, nhmmax, nhrmax, demax, dcamax};

    Int_t    binsa[5] = {npt, ny, neta, nphi, nhf};
    Double_t xmina[5] = {ptmin, ymin, etamin, phimin, nhfmin};
    Double_t xmaxa[5] = {ptmax, ymax, etamax, phimax, nhfmax};

    Int_t    binsb[5] = {npt, ny, neta, nphi, nhm};
    Double_t xminb[5] = {ptmin, ymin, etamin, phimin, nhmmin};
    Double_t xmaxb[5] = {ptmax, ymax, etamax, phimax, nhmmax};

    Int_t    binsc[5] = {npt, ny, neta, nphi, nhr};
    Double_t xminc[5] = {ptmin, ymin, etamin, phimin, nhrmin};
    Double_t xmaxc[5] = {ptmax, ymax, etamax, phimax, nhrmax};

    Int_t    binsd[5] = {npt, ny, neta, nphi, nde};
    Double_t xmind[5] = {ptmin, ymin, etamin, phimin, demin};
    Double_t xmaxd[5] = {ptmax, ymax, etamax, phimax, demax};

    Int_t    binse[5] = {npt, ny, neta, nphi, ndca};
    Double_t xmine[5] = {ptmin, ymin, etamin, phimin, dcamin};
    Double_t xmaxe[5] = {ptmax, ymax, etamax, phimax, dcamax};

    string charge[3] = {"plus","minus","phi"};
    string mixing[5] = {"SE","ME","SM","RC","MC"};
  
    for(int is = 0; is < 5; is++)
    {
      for(int ic = 0; ic < 2; ic++)
      {
        string hist;
        hist = Form("kaonptcos_%d_%d",ic,is);
        kaonptcos[ic][is] = new TH2F(hist.c_str(), hist.c_str(), 25,0.0,5.0,10,-1,1);
      }
    }


    for (int im = 0; im < 5; im++)
    {
      for (int ic = 0; ic < 3; ic++)
      {
        string hist;

        // TH1F 
        hist = Form("pt_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        pt[im][ic] = new TH1F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax);
        hist = Form("rapidity_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        rapidity[im][ic] = new TH1F(hist.c_str(), hist.c_str(), ny, ymin, ymax);
        hist = Form("eta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        eta[im][ic] = new TH1F(hist.c_str(), hist.c_str(), neta, etamin, etamax);
        hist = Form("phi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        phi[im][ic] = new TH1F(hist.c_str(), hist.c_str(), nphi, phimin, phimax);
        // TH1F 

        //hist = Form("ptyetaphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        //ptyetaphi[im][ic] = new THnF(hist.c_str(), hist.c_str(), 4, bins, xmin, xmax);

        // TH2F 
        hist = Form("ptrapidity_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        ptrapidity[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, ny, ymin, ymax);
        hist = Form("pteta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        pteta[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, neta, etamin, etamax);
        hist = Form("ptphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        ptphi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nphi, phimin, phimax);
        hist = Form("rapidityeta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        rapidityeta[im][ic] = new TH2F(hist.c_str(), hist.c_str(), ny, ymin, ymax, neta, etamin, etamax);
        hist = Form("rapidityphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        rapidityphi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nphi, phimin, phimax);
        hist = Form("etaphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        etaphi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), neta, etamin, etamax, nphi, phimin, phimax);

        //hist = Form("ptrapidityphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        //ptrapidityphi[im][ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, ny, ymin, ymax, nphi, phimin, phimax);
        // TH2F 

        //if(ic < 2) 
        //{
 
           
            //hist = Form("ptyetaphibasic_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
            //ptyetaphibasic[im][ic] = new THnF(hist.c_str(), hist.c_str(), 9, binsb, xminb, xmaxb);
            hist = Form("ptyetaphibasic_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
            ptyetaphibasic[im][ic] = new THnF(hist.c_str(), hist.c_str(), 4, bins, xmin, xmax);
            ptyetaphibasic[im][ic]->Sumw2();
            //hist = Form("ptyetaphia_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
            //ptyetaphia[im][ic] = new THnF(hist.c_str(), hist.c_str(), 5, binsa, xmina, xmaxa);
            //hist = Form("ptyetaphib_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
            //ptyetaphib[im][ic] = new THnF(hist.c_str(), hist.c_str(), 5, binsb, xminb, xmaxb);
            //hist = Form("ptyetaphic_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
            //ptyetaphic[im][ic] = new THnF(hist.c_str(), hist.c_str(), 5, binsc, xminc, xmaxc);
            //hist = Form("ptyetaphid_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
            //ptyetaphid[im][ic] = new THnF(hist.c_str(), hist.c_str(), 5, binsd, xmind, xmaxd);
            //hist = Form("ptyetaphie_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
            //ptyetaphie[im][ic] = new THnF(hist.c_str(), hist.c_str(), 5, binse, xmine, xmaxe);

          //// TH1F Basic Variables
          //hist = Form("nhitsfit_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsfit[im][ic] = new TH1F(hist.c_str(), hist.c_str(), nhf, nhfmin, nhfmax);
          //hist = Form("nhitsmax_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsmax[im][ic] = new TH1F(hist.c_str(), hist.c_str(), nhm, nhmmin, nhmmax);
          //hist = Form("nhitsratio_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsratio[im][ic] = new TH1F(hist.c_str(), hist.c_str(), nhr, nhrmin, nhrmax);
          //hist = Form("dedx_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //dedx[im][ic] = new TH1F(hist.c_str(), hist.c_str(), nde, demin, demax);
          //hist = Form("dca_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //dca[im][ic] = new TH1F(hist.c_str(), hist.c_str(), ndca, dcamin, dcamax);
          //// TH1F Basic Variables

          //// TH2F Basic Vairables PT        
          //hist = Form("nhitsfitpt_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsfitpt[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nhf, nhfmin, nhfmax);
          //hist = Form("nhitsmaxpt_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsmaxpt[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nhm, nhmmin, nhmmax);
          //hist = Form("nhitsratiopt_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsratiopt[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nhr, nhrmin, nhrmax);
          //hist = Form("dedxpt_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //dedxpt[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nde, demin, demax);
          //hist = Form("dcapt_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //dcapt[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, ndca, dcamin, dcamax);

          //// TH2F Basic Vairables rapidity        
          //hist = Form("nhitsfitrapidity_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsfitrapidity[im][ic] = new TH2F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nhf, nhfmin, nhfmax);
          //hist = Form("nhitsmaxrapidity_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsmaxrapidity[im][ic] = new TH2F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nhm, nhmmin, nhmmax);
          //hist = Form("nhitsratiorapidity_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsratiorapidity[im][ic] = new TH2F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nhr, nhrmin, nhrmax);
          //hist = Form("dedxrapidity_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //dedxrapidity[im][ic] = new TH2F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nde, demin, demax);
          //hist = Form("dcarapidity_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //dcarapidity[im][ic] = new TH2F(hist.c_str(), hist.c_str(), ny, ymin, ymax, ndca, dcamin, dcamax);

          //// TH2F Basic Vairables eta        
          //hist = Form("nhitsfiteta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsfiteta[im][ic] = new TH2F(hist.c_str(), hist.c_str(), neta, etamin, etamax, nhf, nhfmin, nhfmax);
          //hist = Form("nhitsmaxeta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsmaxeta[im][ic] = new TH2F(hist.c_str(), hist.c_str(), neta, etamin, etamax, nhm, nhmmin, nhmmax);
          //hist = Form("nhitsratioeta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsratioeta[im][ic] = new TH2F(hist.c_str(), hist.c_str(), neta, etamin, etamax, nhr, nhrmin, nhrmax);
          //hist = Form("dedxeta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //dedxeta[im][ic] = new TH2F(hist.c_str(), hist.c_str(), neta, etamin, etamax, nde, demin, demax);
          //hist = Form("dcaeta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //dcaeta[im][ic] = new TH2F(hist.c_str(), hist.c_str(), neta, etamin, etamax, ndca, dcamin, dcamax);


          //// TH2F Basic Vairables phi        
          //hist = Form("nhitsfitphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsfitphi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), nphi, phimin, phimax, nhf, nhfmin, nhfmax);
          //hist = Form("nhitsmaxphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsmaxphi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), nphi, phimin, phimax, nhm, nhmmin, nhmmax);
          //hist = Form("nhitsratiophi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //nhitsratiophi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), nphi, phimin, phimax, nhr, nhrmin, nhrmax);
          //hist = Form("dedxphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //dedxphi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), nphi, phimin, phimax, nde, demin, demax);
          //hist = Form("dcaphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
          //dcaphi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), nphi, phimin, phimax, ndca, dcamin, dcamax);

        //}
      }
    }
    cout << "Set up histograms" << endl;

    // INITIALIZE HISTOGRAMS

    Long64_t nentries = mKaonTreeSE->GetEntries();
    //nentries = 1000000;// mKaonTreeSE->GetEntries();
    cout << "nentries = " << nentries << endl;
    for (Long64_t i = 0; i < nentries; i++)
    {
      mKaonTreeSE->GetEntry(i);
      if (i%1000000 == 0)
      {
         cout << "Proccessed " << i << endl;
         
         //cout << mCentSE << endl;
         //cout << mWeightSE << endl;
         //cout << mChargeSE << endl;
         //cout << mPtSE << endl;
         //cout << mRapiditySE << endl;
         //cout << mEtaSE << endl;
         //cout << mPhiSE << endl;
         //cout << mNHitsFitSE << endl;
         //cout << mNHitsMaxSE << endl;
         //cout << mDEdxSE << endl;
         //cout << mDcaSE << endl;
      }

      // TH1F 
      //pt[0][mChargeSE]->Fill(mPtSE,mWeightSE);
      //rapidity[0][mChargeSE]->Fill(mRapiditySE,mWeightSE);
      //eta[0][mChargeSE]->Fill(mEtaSE,mWeightSE);
      //phi[0][mChargeSE]->Fill(mPhiSE,mWeightSE);
      // TH1F 

      //ptyetaphi[0][mChargeSE]->Fill(mPtSE,mRapiditySE,mEtaSE,mPhiSE,mWeightSE);
      //ptyetaphibasic[0][mChargeSE]->Fill(mPtSE,mRapiditySE,mEtaSE,mPhiSE,mNHitsFitSE,mNHitsMaxSE,float(mNHitsFitSE)/float(mNHitsMaxSE),mDEdxSE,mDcaSE,mWeightSE);
      ptyetaphibasic[0][mChargeSE]->Fill(mPtSE,mRapiditySE,mEtaSE,mPhiSE,mWeightSE);
      //ptyetaphia[0][mChargeSE]->Fill(mPtSE,mRapiditySE,mEtaSE,mPhiSE,mNHitsFitSE,mWeightSE);
      //ptyetaphib[0][mChargeSE]->Fill(mPtSE,mRapiditySE,mEtaSE,mPhiSE,mNHitsMaxSE,mWeightSE);
      //ptyetaphic[0][mChargeSE]->Fill(mPtSE,mRapiditySE,mEtaSE,mPhiSE,float(mNHitsFitSE)/float(mNHitsMaxSE),mWeightSE);
      //ptyetaphid[0][mChargeSE]->Fill(mPtSE,mRapiditySE,mEtaSE,mPhiSE,mDEdxSE,mWeightSE);
      //ptyetaphie[0][mChargeSE]->Fill(mPtSE,mRapiditySE,mEtaSE,mPhiSE,mDcaSE,mWeightSE);

      // TH2F 
      //ptrapidity[0][mChargeSE]->Fill(mPtSE,mRapiditySE,mWeightSE);
      //pteta[0][mChargeSE]->Fill(mPtSE,mEtaSE,mWeightSE);
      //ptphi[0][mChargeSE]->Fill(mPtSE,mPhiSE,mWeightSE);
      //rapidityeta[0][mChargeSE]->Fill(mRapiditySE,mEtaSE,mWeightSE);
      //rapidityphi[0][mChargeSE]->Fill(mRapiditySE,mPhiSE,mWeightSE);
      //etaphi[0][mChargeSE]->Fill(mEtaSE,mPhiSE,mWeightSE);
      // TH2F 

      //nhitsfit[0][mChargeSE]  ->Fill(mNHitsFitSE,mWeightSE); 
      //nhitsmax[0][mChargeSE]  ->Fill(mNHitsMaxSE,mWeightSE); 
      //nhitsratio[0][mChargeSE]->Fill(float(mNHitsFitSE)/float(mNHitsMaxSE),mWeightSE);
      //dedx[0][mChargeSE]      ->Fill(mDEdxSE,mWeightSE);
      //dca[0][mChargeSE]       ->Fill(mDcaSE,mWeightSE);

    }
      
    Long64_t nentriesME = mKaonTreeME->GetEntries();
    //nentriesME = 1000000;//mKaonTreeME->GetEntries();
    cout << "nentriesME = " << nentriesME << endl;
    for (Long64_t i = 0; i < nentriesME; i++)
    {
      mKaonTreeME->GetEntry(i);
      if (i%1000000 == 0)
      {
         cout << "Proccessed " << i << endl;
         
         //cout << mCentME << endl;
         //cout << mWeightME << endl;
         //cout << mChargeME << endl;
         //cout << mPtME << endl;
         //cout << mRapidityME << endl;
         //cout << mEtaME << endl;
         //cout << mPhiME << endl;
         //cout << mNHitsFitME << endl;
         //cout << mNHitsMaxME << endl;
         //cout << mDEdxME << endl;
         //cout << mDcaME << endl;
      }

      // TH1F 
      //pt[1][mChargeME]->Fill(mPtME,mWeightME);
      //rapidity[1][mChargeME]->Fill(mRapidityME,mWeightME);
      //eta[1][mChargeME]->Fill(mEtaME,mWeightME);
      //phi[1][mChargeME]->Fill(mPhiME,mWeightME);
      // TH1F 

      //ptyetaphi[1][mChargeME]->Fill(mPtME,mRapidityME,mEtaME,mPhiME,mWeightME);
      //ptyetaphibasic[0][mChargeME]->Fill(mPtME,mRapidityME,mEtaME,mPhiME,mNHitsFitME,mNHitsMaxME,float(mNHitsFitME)/float(mNHitsMaxME),mDEdxME,mDcaME,mWeightME);
      ptyetaphibasic[1][mChargeME]->Fill(mPtME,mRapidityME,mEtaME,mPhiME,mWeightME);
      //ptyetaphia[0][mChargeME]->Fill(mPtME,mRapidityME,mEtaME,mPhiME,mNHitsFitME,mWeightME);
      //ptyetaphib[0][mChargeME]->Fill(mPtME,mRapidityME,mEtaME,mPhiME,mNHitsMaxME,mWeightME);
      //ptyetaphic[0][mChargeME]->Fill(mPtME,mRapidityME,mEtaME,mPhiME,float(mNHitsFitME)/float(mNHitsMaxME),mWeightME);
      //ptyetaphid[0][mChargeME]->Fill(mPtME,mRapidityME,mEtaME,mPhiME,mDEdxME,mWeightME);
      //ptyetaphie[0][mChargeME]->Fill(mPtME,mRapidityME,mEtaME,mPhiME,mDcaME,mWeightME);

      // TH2F 
      //ptrapidity[1][mChargeME]->Fill(mPtME,mRapidityME,mWeightME);
      //pteta[1][mChargeME]->Fill(mPtME,mEtaME,mWeightME);
      //ptphi[1][mChargeME]->Fill(mPtME,mPhiME,mWeightME);
      //rapidityeta[1][mChargeME]->Fill(mRapidityME,mEtaME,mWeightME);
      //rapidityphi[1][mChargeME]->Fill(mRapidityME,mPhiME,mWeightME);
      //etaphi[1][mChargeME]->Fill(mEtaME,mPhiME,mWeightME);
      // TH2F 

      //nhitsfit[1][mChargeME]->Fill(mNHitsFitME,mWeightME); 
      //nhitsmax[1][mChargeME]->Fill(mNHitsMaxME,mWeightME); 
      //nhitsratio[1][mChargeME]->Fill(float(mNHitsFitME)/float(mNHitsMaxME),mWeightME);
      //dedx[1][mChargeME]->Fill(mDEdxME,mWeightME);
      //dca[1][mChargeME]->Fill(mDcaME,mWeightME);
    }
     
    
    //Long64_t nentriesEmbed = mKaonTreeEmbed->GetEntries();
    //cout << "nentriesEmbed = " << nentriesEmbed << endl;
    //for (Long64_t i = 0; i < nentriesEmbed; i++)
    //{
    //  mKaonTreeEmbed->GetEntry(i);
    //  if (i%1000000 == 0)
    //  {
    //     cout << "Proccessed " << i << endl;
    //     
    //     //cout << mCentEmbed << endl;
    //     //cout << mWeight << endl;
    //     //cout << 0 << endl;
    //     //cout << mPtP << endl;
    //     //cout << mRapidityP << endl;
    //     //cout << mEtaP << endl;
    //     //cout << mPhiP << endl;
    //     //cout << mNHitsFitEmbed << endl;
    //     //cout << mNHitsMaxEmbed << endl;
    //     //cout << mDEdxEmbed << endl;
    //     //cout << mDcaEmbed << endl;
    //  }
    //
    //  if(mCent < 2 || mCent > 5) continue;

    //  float mMcPt = mLMeson.Pt();      
    //  float mMcRapidity = mLMeson.Rapidity();      
    //  float mMcPhi = mLMeson.Phi();      

    //  ptrapidityphi[mCent]->Fill(mMcPt,mMcRapidity,mMcPhi,mWeight);
    //} 
    Long64_t nentriesEmbed = mKaonTreeEmbed->GetEntries();
    //nentriesEmbed = 1000000;//mKaonTreeEmbed->GetEntries();

   
    string inputfile = Form("../output/AuAu%s/Phi/ptyspectra_Reweight_RawPhiPtSys_eta1_eta1_PolySys.root",vmsa::mBeamEnergy[energy].c_str());
    TFile *File_InPut = TFile::Open(inputfile.c_str());
    string KEY = Form("ptyspectra_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,"2nd",0,0,vmsa::mPID[0].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
    TH2F *ptrapidityPhiData = (TH2F*) File_InPut->Get(KEY.c_str());
    TH2F *ptrapidityPhiRc = new TH2F("ptrapidityPhiRC","ptrapidityPhiRC",vmsa::tofrebinpttotal,vmsa::tofrebinptval,vmsa::tofrebinytotal,vmsa::tofrebinyval); 
    TH2F *ptrapidityPhiRcAfter = new TH2F("ptrapidityPhiRCAfter","ptrapidityPhiRCAfter",vmsa::tofrebinpttotal,vmsa::tofrebinptval,vmsa::tofrebinytotal,vmsa::tofrebinyval); 

    int ncostheta = 14;
    float costheta[ncostheta+1];
    for(int i = 0; i < ncostheta+1; i++)
    {
      costheta[i] = (float(i)-7.)/7.;
    }
    TH3F *CosThetaStarPhiRc = new TH3F("CosThetaStarPhiRC","CosThetaStarPhiRC",vmsa::tofrebinpttotal,vmsa::tofrebinptval,vmsa::tofrebinytotal,vmsa::tofrebinyval,ncostheta,costheta); 
    TH3F *CosThetaStarPhiRcAfter = new TH3F("CosThetaStarPhiRCAfter","CosThetaStarPhiRCAfter",vmsa::tofrebinpttotal,vmsa::tofrebinptval,vmsa::tofrebinytotal,vmsa::tofrebinyval,ncostheta,costheta); 

    double norm[12] = {1.17421, 1.20569, 1.18756, 1.20748, 1.22335, 1.23579, 1.20088, 1.21209, 1.26302, 1.22473, 1.29652, 1.34468};
    double mean[12] = {-0.00329925, -0.00254986, 0.0106082, 0.00417584, 0.00726663, 0.0287173, -0.0166698, -6.69877e-05, -0.0367374, -0.0154559, 0.0763402, -0.0402667};
    double sigma[12] = {0.899472, 0.838369, 0.875892, 0.843554, 0.809015, 0.791428, 0.844824, 0.825385, 0.743799, 0.787513, 0.698393, 0.662743};

    TF1* pythiaflat[12];
    const double  pythialow[12] = {1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};
    const double pythiahigh[12] = {1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 4.2};
    
    for(int i = 0; i < 12; i++)
    { 
      pythiaflat[i]= new TF1(Form("pythiaflat_%d",i),"[0]*exp(-(x-[1])*(x-[1])/2/[2]/[2])",-1,1);
      pythiaflat[i]->SetParameter(0,norm[i]);
      pythiaflat[i]->SetParameter(1,mean[i]);
      pythiaflat[i]->SetParameter(2,sigma[i]);
    }
 
    cout << "nentriesEmbed = " << nentriesEmbed << endl;
    for (Long64_t i = 0; i < nentriesEmbed; i++)
    {
      mKaonTreeEmbed->GetEntry(i);
      if (i%1000000 == 0)
      {
         cout << "Proccessed " << i << endl;
      }

      if(mNToFMatch < 2) continue; 
      if(mCent < 2 || mCent > 5) continue;

      float mMcPt = mLMeson->Pt();      
      float mMcRapidity = mLMeson->Rapidity();      
      float mMcPhi = mLMeson->Phi();      

      float mMcPhiPsi = mMcPhi - mEpFull;
      while(mMcPhiPsi < 0.0) mMcPhiPsi += 2.0*TMath::Pi();
      while(mMcPhiPsi > 2.0*TMath::Pi()) mMcPhiPsi -= 2.0*TMath::Pi();
     
      float PhiWeight = h_mpTyv2Weights[mCent]->GetBinContent(h_mpTyv2Weights[mCent]->FindBin(mMcPt,mMcRapidity,mMcPhiPsi));
      //CosThetaStarPhiRc->Fill(mMcPt,mMcRapidity,mMcCosTheta,mWeight*PhiWeight);
      float pythiaweight = 1.0;
      for(int i = 0; i < 12; i++)
      {
        if(mMcPt >= pythialow[i] && mMcPt < pythiahigh[i]) 
        {
          pythiaweight = pythiaflat[i]->Eval(mMcRapidity);
          break;
        } 
      }

      if(mRcExists && mPassTpcKp && mPassTpcKm /*&& mPassTof*/)    
      {
        bool passToF = false;

        float mRcPtP = mLRcKp->Pt();      
        float mRcRapidityP = mLRcKp->Rapidity();      
        float mRcEtaP = mLRcKp->Eta();      
        float mRcPhiP = mLRcKp->Phi();      

        float mRcPtM = mLRcKm->Pt();      
        float mRcRapidityM = mLRcKm->Rapidity();      
        float mRcEtaM = mLRcKm->Eta();      
        float mRcPhiM = mLRcKm->Phi();      

        int GlobalBinP = -1;
        int GlobalBinM = -1;
        GlobalBinP = ToFHist[0]->FindBin(mRcPtP,mRcEtaP,mRcPhiP);
        GlobalBinM = ToFHist[1]->FindBin(mRcPtM,mRcEtaM,mRcPhiM);
        if(GlobalBinP < 0 || GlobalBinM < 0) passToF = false;   
        float valToFP = ToFHist[0]->GetBinContent(GlobalBinP);
        float valToFM = ToFHist[1]->GetBinContent(GlobalBinM);
      
        float probP = gRandom->Uniform(0,1);
        float probM = gRandom->Uniform(0,1);

        if(probP < valToFP && probM < valToFM) passToF = true;


        if(passToF && mPassM2Kp && mPassM2Km && mPassNsigKp && mPassNsigKm)    
        {
          TLorentzVector mLRcMeson = *mLRcKp + *mLRcKm;

          float mRcPt = mLRcMeson.Pt();      
          float mRcRapidity = mLRcMeson.Rapidity();      

          ptrapidityPhiRc->Fill(mRcPt,mRcRapidity,mWeight*PhiWeight*pythiaweight);
          //ptrapidityPhiRc->Fill(mMcPt,mMcRapidity,mWeight*PhiWeight);
          CosThetaStarPhiRc->Fill(mRcPt,mRcRapidity,mRcCosTheta,mWeight*PhiWeight*pythiaweight);
        }
      }
    }

    TH2F *ptrapidityPhiDataRc;

    double integralPhiData = ptrapidityPhiData->Integral(1,vmsa::tofrebinpttotal,1,vmsa::tofrebinytotal); 
    double integralPhiRc = ptrapidityPhiRc->Integral(1,vmsa::tofrebinpttotal,1,vmsa::tofrebinytotal); 
    double DataRcPhi = integralPhiData/integralPhiRc;
    
    ptrapidityPhiDataRc = (TH2F*) ptrapidityPhiData->Clone();
    ptrapidityPhiRc->Scale(DataRcPhi);
    ptrapidityPhiDataRc->Divide(ptrapidityPhiRc);

    TH1F *ptPhiData = (TH1F*) ptrapidityPhiData->ProjectionX("ptData",0,-1);
    //TH1F *ptPhiData = (TH1F*) ptrapidityPhiData->ProjectionX("ptData",4,11);
    TH1F *yPhiData  = (TH1F*) ptrapidityPhiData->ProjectionY("yData",0,-1);  
    TH1F *ptPhiRc   = (TH1F*) ptrapidityPhiRc->ProjectionX("ptRc",0,-1);  
    //TH1F *ptPhiRc   = (TH1F*) ptrapidityPhiRc->ProjectionX("ptRc",4,11);  
    TH1F *yPhiRc    = (TH1F*) ptrapidityPhiRc->ProjectionY("yRc",0,-1);  

    TH1F *ptPhiDataRc = (TH1F*) ptPhiData->Clone();
    ptPhiDataRc->Divide(ptPhiRc);
    TH1F *yPhiDataRc = (TH1F*) yPhiData->Clone();
    yPhiDataRc->Divide(yPhiRc);


    TCanvas *ct = new TCanvas("ct", "ct", 1200, 1200);
    ct->Divide(3,3);
    for(int i = 0; i < 9; i++)
    {
      ct->cd(i+1);
      ct->cd(i+1)->SetLeftMargin(0.12);
      ct->cd(i+1)->SetRightMargin(0.15);
      ct->cd(i+1)->SetBottomMargin(0.12);
      ct->cd(i+1)->SetTicks(1,1);
      ct->cd(i+1)->SetGrid(0,0); 
    }

    ct->cd(1);
    //ct->cd(1)->SetLogz();
    ptrapidityPhiData->SetTitle("Data"); 
    ptrapidityPhiData->GetXaxis()->SetTitle("p_{T} GeV/c"); 
    ptrapidityPhiData->GetYaxis()->SetTitle("y"); 
    ptrapidityPhiData->Draw("Colz"); 
    ct->cd(2);
    //ct->cd(2)->SetLogz();
    ptrapidityPhiRc->SetTitle("Embedding"); 
    ptrapidityPhiRc->GetXaxis()->SetTitle("p_{T} GeV/c"); 
    ptrapidityPhiRc->GetYaxis()->SetTitle("y"); 
    ptrapidityPhiRc->Draw("Colz"); 
    ct->cd(3);
    ct->cd(3)->SetLogy(0);
    ptrapidityPhiDataRc->SetTitle("Data/Embedding"); 
    ptrapidityPhiDataRc->GetXaxis()->SetTitle("p_{T} GeV/c"); 
    ptrapidityPhiDataRc->GetYaxis()->SetTitle("y"); 
    ptrapidityPhiDataRc->Draw("Colz"); 

    ct->cd(4);
    ct->cd(4)->SetLogy();
    ptPhiData->SetTitle("Before Data/RC Reweight #phi-meson"); 
    ptPhiData->GetXaxis()->SetTitle("p_{T} GeV/c"); 
    ptPhiData->GetYaxis()->SetTitle("Counts"); 
    ptPhiData->SetMarkerColor(kBlue); 
    ptPhiData->SetLineColor(kBlue); 
    ptPhiData->SetMarkerStyle(20); 
    ptPhiData->Draw("pE"); 
    ptPhiRc->SetMarkerColor(kOrange+7); 
    ptPhiRc->SetLineColor(kOrange+7); 
    ptPhiRc->SetMarkerStyle(20); 
    ptPhiRc->Draw("pE same");
    TLegend *legpt = new TLegend(0.6,0.6,0.8,0.8);
    legpt->AddEntry(ptPhiData,"Data","p");
    legpt->AddEntry(ptPhiRc,"RC","p");
    legpt->Draw("same");

    ct->cd(5);
    //ct->cd(5)->SetLogy();
    yPhiData->SetTitle("Before Data/RC Reweight #phi-meson"); 
    yPhiData->GetXaxis()->SetTitle("y"); 
    yPhiData->GetYaxis()->SetTitle("Counts"); 
    yPhiData->SetMarkerColor(kBlue); 
    yPhiData->SetLineColor(kBlue); 
    yPhiData->SetMarkerStyle(20); 
    yPhiData->Draw("pE"); 
    yPhiRc->SetMarkerColor(kOrange+7); 
    yPhiRc->SetLineColor(kOrange+7); 
    yPhiRc->SetMarkerStyle(20); 
    yPhiRc->Draw("pE same");
    TLegend *legy = new TLegend(0.4,0.2,0.6,0.4);
    legy->AddEntry(yPhiData,"Data","p");
    legy->AddEntry(yPhiRc,"RC","p");
    legy->Draw("same");

    ct->cd(7);
    ptPhiDataRc->SetTitle("Before Data/RC Reweight #phi-meson"); 
    ptPhiDataRc->GetXaxis()->SetTitle("p_{T} GeV/c"); 
    ptPhiDataRc->GetYaxis()->SetTitle("Data/RC"); 
    ptPhiDataRc->Draw("pE");

    ct->cd(8);
    yPhiDataRc->SetTitle("Before Data/RC Reweight #phi-meson"); 
    yPhiDataRc->GetXaxis()->SetTitle("y"); 
    yPhiDataRc->GetYaxis()->SetTitle("Data/RC"); 
    yPhiDataRc->Draw("pE");
    TF1* gausRatio = new TF1("gausRatio","[0]*exp(-x*x/2/[1]/[1])",-1.0,1.0);
    //gausRatio->SetParameter(0,1.);
    //gausRatio->SetParameter(1,1.);
    //yPhiDataRc->Fit(gausRatio,"NMR");
    //gausRatio->SetLineColor(kRed);
    //gausRatio->Draw("l same");

    ct->SaveAs(Form("figures/KaonTTrees/%sPhiMeson_DataRCAndRatio.pdf",folderopt.c_str()));
    
    TCanvas *ct2 = new TCanvas("ct2", "ct2", 1500, 600);
    ct2->Divide(5,2);

    TH1F *phi_ybins[10];   
 
    for(int i = 0; i < 10; i++)
    {
      ct2->cd(i+1);
      ct2->cd(i+1)->SetLeftMargin(0.12);
      ct2->cd(i+1)->SetRightMargin(0.15);
      ct2->cd(i+1)->SetBottomMargin(0.12);
      ct2->cd(i+1)->SetTicks(1,1);
      ct2->cd(i+1)->SetGrid(0,0); 
    }
    for(int i = 0; i < 10; i++)
    {
      ct2->cd(i+1); 
      string histname = Form("phi_datarc_y_%d",i);
      phi_ybins[i] = (TH1F*) ptrapidityPhiDataRc->ProjectionY(histname.c_str(),i+1,i+1);
      phi_ybins[i]->GetXaxis()->SetTitle("y");  
      phi_ybins[i]->GetYaxis()->SetTitle("Data/RC");  
      phi_ybins[i]->Draw("pE"); 
    }
    ct2->SaveAs(Form("figures/KaonTTrees/%sPhiMeson_DataRCAndRatio_rapiditybins.pdf",folderopt.c_str()));

    //for(int i = 0; i < 10; i++)
    //{
    //  ct2->cd(i+1); 
    //  string histname = Form("phi_datarc_y_%d",i);
    //  phi_ybins[i] = (TH1F*) ptrapidityPhiDataRc->ProjectionY(histname.c_str(),i+1,i+1);
    //  phi_ybins[i]->GetXaxis()->SetTitle("y");  
    //  phi_ybins[i]->GetYaxis()->SetTitle("Data/RC");  
    //  phi_ybins[i]->Draw("pE"); 
    //}
    //ct2->SaveAs(Form("figures/KaonTTrees/%sPhiMeson_DataRCAndRatio_rapiditybins.pdf",folderopt.c_str()));

    //Long64_t nentriesEmbed = mKaonTreeEmbed->GetEntries();
    //nentriesEmbed = 1000000;//mKaonTreeEmbed->GetEntries();
    cout << "nentriesEmbed = " << nentriesEmbed << endl;



    for (Long64_t i = 0; i < nentriesEmbed; i++)
    {
      mKaonTreeEmbed->GetEntry(i);
      if (i%1000000 == 0)
      {
         cout << "Proccessed " << i << endl;
      }

      if(mNToFMatch < 2) continue; 
      if(mCent < 2 || mCent > 5) continue;

      float mMcPt = mLMeson->Pt();      
      float mMcRapidity = mLMeson->Rapidity();      
      float mMcEta = mLMeson->Eta();      
      float mMcPhi = mLMeson->Phi();      

      float mMcPtP = mLKp->Pt();      
      float mMcRapidityP = mLKp->Rapidity();      
      float mMcEtaP = mLKp->Eta();      
      float mMcPhiP = mLKp->Phi();      

      float mMcPtM = mLKm->Pt();      
      float mMcRapidityM = mLKm->Rapidity();      
      float mMcEtaM = mLKm->Eta();      
      float mMcPhiM = mLKm->Phi();      


      //// TH1F 
      //pt[4][0]->Fill(mMcPtP,mWeight);
      //rapidity[4][0]->Fill(mMcRapidityP,mWeight);
      //eta[4][0]->Fill(mMcEtaP,mWeight);
      //phi[4][0]->Fill(mMcPhiP,mWeight);
      //pt[4][1]->Fill(mMcPtM,mWeight);
      //rapidity[4][1]->Fill(mMcRapidityM,mWeight);
      //eta[4][1]->Fill(mMcEtaM,mWeight);
      //phi[4][1]->Fill(mMcPhiM,mWeight);

      //pt[4][2]->Fill(mMcPt,mWeight);
      //rapidity[4][2]->Fill(mMcRapidity,mWeight);
      //eta[4][2]->Fill(mMcEta,mWeight);
      //phi[4][2]->Fill(mMcPhi,mWeight);
      //// TH1F 

      //// TH2F 
      //ptrapidity[4][0]->Fill(mMcPtP,mMcRapidityP,mWeight);
      //pteta[4][0]->Fill(mMcPtP,mMcEtaP,mWeight);
      //ptphi[4][0]->Fill(mMcPtP,mMcPhiP,mWeight);
      //rapidityeta[4][0]->Fill(mMcRapidityP,mMcEtaP,mWeight);
      //rapidityphi[4][0]->Fill(mMcRapidityP,mMcPhiP,mWeight);
      //etaphi[4][0]->Fill(mMcEtaP,mMcPhiP,mWeight);
      //ptrapidity[4][1]->Fill(mMcPtM,mMcRapidityM,mWeight);
      //pteta[4][1]->Fill(mMcPtM,mMcEtaM,mWeight);
      //ptphi[4][1]->Fill(mMcPtM,mMcPhiM,mWeight);
      //rapidityeta[4][1]->Fill(mMcRapidityM,mMcEtaM,mWeight);
      //rapidityphi[4][1]->Fill(mMcRapidityM,mMcPhiM,mWeight);
      //etaphi[4][1]->Fill(mMcEtaM,mMcPhiM,mWeight);

      //ptrapidity[4][2]->Fill(mMcPt,mMcRapidity,mWeight);
      //pteta[4][2]->Fill(mMcPt,mMcEta,mWeight);
      //ptphi[4][2]->Fill(mMcPt,mMcPhi,mWeight);
      //rapidityeta[4][2]->Fill(mMcRapidity,mMcEta,mWeight);
      //rapidityphi[4][2]->Fill(mMcRapidity,mMcPhi,mWeight);
      //etaphi[4][2]->Fill(mMcEta,mMcPhi,mWeight);
      //// TH2F 

      //cout << "AM I HERE? " << endl;
      //cout << "mCent = " << mCent << endl;


      //cout << " I DID GET HERE" << endl;

      //cout << "mEpFull = " << mEpFull << endl;

      float mMcPhiPsi = mMcPhi - mEpFull;
      while(mMcPhiPsi < 0.0) mMcPhiPsi += 2.0*TMath::Pi();
      while(mMcPhiPsi > 2.0*TMath::Pi()) mMcPhiPsi -= 2.0*TMath::Pi();
     
      //cout << "mMcPhiPsi = " << mMcPhiPsi << endl;

      float PhiWeight = h_mpTyv2Weights[mCent]->GetBinContent(h_mpTyv2Weights[mCent]->FindBin(mMcPt,mMcRapidity,mMcPhiPsi));
      float pythiaweight = 1.0;
      for(int i = 0; i < 12; i++)
      {
        if(mMcPt >= pythialow[i] && mMcPt < pythiahigh[i]) 
        {
          pythiaweight = pythiaflat[i]->Eval(mMcRapidity);
          break;
        } 
      }

      //float pythiaweight = pythiaflat->Eval(mRcRapidity);

      if(mMcPt < 1.2 || mMcPt > 4.2) continue;
      if(TMath::Abs(mMcRapidity) > 0.4) continue;
           
      ptyetaphibasic[4][0]->Fill(mMcPtP,mMcRapidityP,mMcEtaP,mMcPhiP,mWeight*PhiWeight*pythiaweight);
      ptyetaphibasic[4][1]->Fill(mMcPtM,mMcRapidityM,mMcEtaM,mMcPhiM,mWeight*PhiWeight*pythiaweight);
      ptyetaphibasic[4][2]->Fill(mMcPt,mMcRapidity,mMcEta,mMcPhi,mWeight*PhiWeight*pythiaweight);

      //float DataRcWeightVal = 1.0;
      //if(datarcweight) DataRcWeightVal = ptrapidityPhiDataRc->GetBinContent(ptrapidityPhiDataRc->FindBin(mRcPt,mRcRapidity));
      //if(datarcweight) DataRcWeightVal = ptrapidityPhiDataRc->GetBinContent(ptrapidityPhiDataRc->FindBin(mMcPt,mMcRapidity));

      //cout << "mPassTpc  = " << mPassTpc  << endl;
      //cout << "mPassTof  = " << mPassTof  << endl;
      //cout << "mPassM2   = " << mPassM2   << endl;
      //cout << "mPassNsig = " << mPassNsig << endl;
      //CosThetaStarPhiRcAfter->Fill(mMcPt,mMcRapidity,mMcCosTheta,mWeight*PhiWeight*DataRcWeightVal);
      //ptrapidityPhiRcAfter->Fill(mMcPt,mMcRapidity,mWeight*PhiWeight*DataRcWeightVal);

      int BinPt = -1;
      for(int ipt = vmsa::pt_rebin_first_2D[energy]; ipt < vmsa::pt_rebin_last_2D[energy]; ipt++)
      {
        if(mMcPt >= vmsa::pt_low_2D[energy][ipt] && mMcPt < vmsa::pt_up_2D[energy][ipt]) 
        {
          BinPt = ipt;
          break; 
        }
      }

      kaonptcos[0][0]->Fill(mMcPtP, mMcCosTheta,mWeight*PhiWeight*pythiaweight);
      kaonptcos[1][0]->Fill(mMcPtM,-mMcCosTheta,mWeight*PhiWeight*pythiaweight); // negative because cos(theta*) is related to the k+ and theta* for the kaon is theta*+pi
      

      
      if(BinPt >= 0) 
      {
        pt_correction_global[BinPt][0]->Fill(mMcCosThetaStar,mMcPhiPrime,mWeight*PhiWeight*pythiaweight);
        pt_correction_helicity[BinPt][0]->Fill(mMcCosTheta,mMcHelicityAngle,mWeight*PhiWeight*pythiaweight);
      }
      if(mRcExists /*&& mPassTpc && mPassTof && mPassM2 && mPassNsig*/)    
      {
        TLorentzVector mLRcMeson = *mLRcKp + *mLRcKm;

        float mRcPt = mLRcMeson.Pt();      
        float mRcRapidity = mLRcMeson.Rapidity();      
        float mRcEta = mLRcMeson.Eta();      
        float mRcPhi = mLRcMeson.Phi();      

        float mRcPtP = mLRcKp->Pt();      
        float mRcRapidityP = mLRcKp->Rapidity();      
        float mRcEtaP = mLRcKp->Eta();      
        float mRcPhiP = mLRcKp->Phi();      

        float mRcPtM = mLRcKm->Pt();      
        float mRcRapidityM = mLRcKm->Rapidity();      
        float mRcEtaM = mLRcKm->Eta();      
        float mRcPhiM = mLRcKm->Phi();      



        bool passToFKp = false;
        bool passToFKm = false;

        int GlobalBinP = -1;
        int GlobalBinM = -1;
        GlobalBinP = ToFHist[0]->FindBin(mRcPtP,mRcEtaP,mRcPhiP);
        GlobalBinM = ToFHist[1]->FindBin(mRcPtM,mRcEtaM,mRcPhiM);
        if(GlobalBinP < 0 || GlobalBinM < 0) 
        {
          passToFKp = false;   
          passToFKm = false;   
        }
        float valToFP = ToFHist[0]->GetBinContent(GlobalBinP);
        float valToFM = ToFHist[1]->GetBinContent(GlobalBinM);
      
        float probP = gRandom->Uniform(0,1);
        float probM = gRandom->Uniform(0,1);

        if(probP < valToFP ) passToFKp = true;
        if(probM < valToFM ) passToFKm = true;

        //if(passToF && mPassM2 && mPassNsig)    
        //{
          if(mRcPt < 1.2 || mRcPt > 4.2) continue;
          if(TMath::Abs(mRcRapidity) > 0.4) continue;
  
          float DataRcWeightVal = 1.0;
          //if(datarcweight) DataRcWeightVal = ptrapidityPhiDataRc->GetBinContent(ptrapidityPhiDataRc->FindBin(mRcPt,mRcRapidity));
          if(datarcweight) DataRcWeightVal = gausRatio->Eval(mRcRapidity);

          if(TMath::Abs(mRcEtaP)<1.0)
          {
            if(mPassTpcKp)
            {
              kaonptcos[0][1]->Fill(mMcPtP,mMcCosTheta, mWeight*PhiWeight*pythiaweight);
              if(passToFKp)
              {
                kaonptcos[0][2]->Fill(mMcPtP,mMcCosTheta, mWeight*PhiWeight*pythiaweight);
                if(mPassNsigKp)
                {
                  kaonptcos[0][3]->Fill(mMcPtP,mMcCosTheta, mWeight*PhiWeight*pythiaweight);
                  if(mPassM2Kp)
                  {
                    kaonptcos[0][4]->Fill(mMcPtP,mMcCosTheta, mWeight*PhiWeight*pythiaweight);
                  }
                }
              }
            }
          }
          if(TMath::Abs(mRcEtaM)<1.0)
          {
            if(mPassTpcKm)
            {
              kaonptcos[1][1]->Fill(mMcPtM,-mMcCosTheta, mWeight*PhiWeight*pythiaweight);
              if(passToFKm)
              {
                kaonptcos[1][2]->Fill(mMcPtM,-mMcCosTheta, mWeight*PhiWeight*pythiaweight);
                if(mPassNsigKm)
                {
                  kaonptcos[1][3]->Fill(mMcPtM,-mMcCosTheta, mWeight*PhiWeight*pythiaweight);
                  if(mPassM2Km)
                  {
                    kaonptcos[1][4]->Fill(mMcPtM,-mMcCosTheta, mWeight*PhiWeight*pythiaweight);
                  }
                }
              }
            }
          }

        

          //cout << "mRcPtP       = " << mRcPtP       << endl;  
          //cout << "mRcRapidityP = " << mRcRapidityP << endl;     
          //cout << "mRcEtaP      = " << mRcEtaP      << endl;
          //cout << "mRcPhiP      = " << mRcPhiP      << endl;
  
          //cout << "mRcPtM       = " << mRcPtM       << endl;  
          //cout << "mRcRapidityM = " << mRcRapidityM << endl;      
          //cout << "mRcEtaM      = " << mRcEtaM      << endl;
          //cout << "mRcPhiM      = " << mRcPhiM      << endl;
          ptrapidityPhiRcAfter->Fill(mRcPt,mRcRapidity,mWeight*PhiWeight*DataRcWeightVal*pythiaweight);
          CosThetaStarPhiRcAfter->Fill(mRcPt,mRcRapidity,mRcCosTheta,mWeight*PhiWeight*DataRcWeightVal*pythiaweight);
  
          ptyetaphibasic[3][0]->Fill(mRcPtP,mRcRapidityP,mRcEtaP,mRcPhiP,mWeight*PhiWeight*DataRcWeightVal*pythiaweight);
          ptyetaphibasic[3][1]->Fill(mRcPtM,mRcRapidityM,mRcEtaM,mRcPhiM,mWeight*PhiWeight*DataRcWeightVal*pythiaweight);
          if(BinPt >= 0) 
          {
            pt_correction_global[BinPt][1]->Fill(mRcCosThetaStar,mRcPhiPrime,mWeight*PhiWeight*pythiaweight);
            pt_correction_helicity[BinPt][1]->Fill(mRcCosTheta,mRcHelicityAngle,mWeight*PhiWeight*pythiaweight);
          }
  
          //// TH1F 
          //pt[3][0]->Fill(mRcPtP,mWeight*PhiWeight*DataRcWeightVal);
          //rapidity[3][0]->Fill(mRcRapidityP,mWeight*PhiWeight*DataRcWeightVal);
          //eta[3][0]->Fill(mRcEtaP,mWeight*PhiWeight*DataRcWeightVal);
          //phi[3][0]->Fill(mRcPhiP,mWeight*PhiWeight*DataRcWeightVal);
          //pt[3][1]->Fill(mRcPtM,mWeight*PhiWeight*DataRcWeightVal);
          //rapidity[3][1]->Fill(mRcRapidityM,mWeight*PhiWeight*DataRcWeightVal);
          //eta[3][1]->Fill(mRcEtaM,mWeight*PhiWeight*DataRcWeightVal);
          //phi[3][1]->Fill(mRcPhiM,mWeight*PhiWeight*DataRcWeightVal);
  
          //pt[3][2]->Fill(mRcPt,mWeight*PhiWeight*DataRcWeightVal);
          //rapidity[3][2]->Fill(mRcRapidity,mWeight*PhiWeight*DataRcWeightVal);
          //eta[3][2]->Fill(mRcEta,mWeight*PhiWeight*DataRcWeightVal);
          //phi[3][2]->Fill(mRcPhi,mWeight*PhiWeight*DataRcWeightVal);
          //// TH1F 
  
          //// TH2F 
          //ptrapidity[3][0]->Fill(mRcPtP,mRcRapidityP,mWeight*PhiWeight*DataRcWeightVal);
          //pteta[3][0]->Fill(mRcPtP,mRcEtaP,mWeight*PhiWeight*DataRcWeightVal);
          //ptphi[3][0]->Fill(mRcPtP,mRcPhiP,mWeight*PhiWeight*DataRcWeightVal);
          //rapidityeta[3][0]->Fill(mRcRapidityP,mRcEtaP,mWeight*PhiWeight*DataRcWeightVal);
          //rapidityphi[3][0]->Fill(mRcRapidityP,mRcPhiP,mWeight*PhiWeight*DataRcWeightVal);
          //etaphi[3][0]->Fill(mRcEtaP,mRcPhiP,mWeight*PhiWeight*DataRcWeightVal);
  
          //ptrapidity[3][1]->Fill(mRcPtM,mRcRapidityM,mWeight*PhiWeight*DataRcWeightVal);
          //pteta[3][1]->Fill(mRcPtM,mRcEtaM,mWeight*PhiWeight*DataRcWeightVal);
          //ptphi[3][1]->Fill(mRcPtM,mRcPhiM,mWeight*PhiWeight*DataRcWeightVal);
          //rapidityeta[3][1]->Fill(mRcRapidityM,mRcEtaM,mWeight*PhiWeight*DataRcWeightVal);
          //rapidityphi[3][1]->Fill(mRcRapidityM,mRcPhiM,mWeight*PhiWeight*DataRcWeightVal);
          //etaphi[3][1]->Fill(mRcEtaM,mRcPhiM,mWeight*PhiWeight*DataRcWeightVal);
  
          //ptrapidity[3][2]->Fill(mRcPt,mRcRapidity,mWeight*PhiWeight*DataRcWeightVal);
          //pteta[3][2]->Fill(mRcPt,mRcEta,mWeight*PhiWeight*DataRcWeightVal);
          //ptphi[3][2]->Fill(mRcPt,mRcPhi,mWeight*PhiWeight*DataRcWeightVal);
          //rapidityeta[3][2]->Fill(mRcRapidity,mRcEta,mWeight*PhiWeight*DataRcWeightVal);
          //rapidityphi[3][2]->Fill(mRcRapidity,mRcPhi,mWeight*PhiWeight*DataRcWeightVal);
          //etaphi[3][2]->Fill(mRcEta,mRcPhi,mWeight*PhiWeight*DataRcWeightVal);
          //// TH2F 
  
          //nhitsfit[3][0]->Fill(mNHitsFitP,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsmax[3][0]->Fill(mNHitsMaxP,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsratio[3][0]->Fill(float(mNHitsFitP)/float(mNHitsMaxP),mWeight*PhiWeight*DataRcWeightVal);
          //dedx[3][0]->Fill(mDEdxP,mWeight*PhiWeight*DataRcWeightVal);
          //dca[3][0]->Fill(mDcaP,mWeight*PhiWeight*DataRcWeightVal);
  
          //nhitsfit[3][1]->Fill(mNHitsFitM,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsmax[3][1]->Fill(mNHitsMaxM,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsratio[3][1]->Fill(float(mNHitsFitM)/float(mNHitsMaxM),mWeight*PhiWeight*DataRcWeightVal);
          //dedx[3][1]->Fill(mDEdxM,mWeight*PhiWeight*DataRcWeightVal);
          //dca[3][1]->Fill(mDcaM,mWeight*PhiWeight*DataRcWeightVal);
          
          //nhitsfitpt[3][0]  ->Fill(mRcPtP,mNHitsFitP,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsmaxpt[3][0]  ->Fill(mRcPtP,mNHitsMaxP,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsratiopt[3][0]->Fill(mRcPtP,float(mNHitsFitP)/float(mNHitsMaxP),mWeight*PhiWeight*DataRcWeightVal);
          //dedxpt[3][0]      ->Fill(mRcPtP,mDEdxP,mWeight*PhiWeight*DataRcWeightVal);
          //dcapt[3][0]       ->Fill(mRcPtP,mDcaP,mWeight*PhiWeight*DataRcWeightVal);
  
          //nhitsfitpt[3][1]  ->Fill(mRcPtM,mNHitsFitM,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsmaxpt[3][1]  ->Fill(mRcPtM,mNHitsMaxM,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsratiopt[3][1]->Fill(mRcPtM,float(mNHitsFitM)/float(mNHitsMaxM),mWeight*PhiWeight*DataRcWeightVal);
          //dedxpt[3][1]      ->Fill(mRcPtM,mDEdxM,mWeight*PhiWeight*DataRcWeightVal);
          //dcapt[3][1]       ->Fill(mRcPtM,mDcaM,mWeight*PhiWeight*DataRcWeightVal);
  
          //nhitsfitrapidity[3][0]  ->Fill(mRcRapidityP,mNHitsFitP,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsmaxrapidity[3][0]  ->Fill(mRcRapidityP,mNHitsMaxP,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsratiorapidity[3][0]->Fill(mRcRapidityP,float(mNHitsFitP)/float(mNHitsMaxP),mWeight*PhiWeight*DataRcWeightVal);
          //dedxrapidity[3][0]      ->Fill(mRcRapidityP,mDEdxP,mWeight*PhiWeight*DataRcWeightVal);
          //dcarapidity[3][0]       ->Fill(mRcRapidityP,mDcaP,mWeight*PhiWeight*DataRcWeightVal);
  
          //nhitsfitrapidity[3][1]  ->Fill(mRcRapidityM,mNHitsFitM,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsmaxrapidity[3][1]  ->Fill(mRcRapidityM,mNHitsMaxM,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsratiorapidity[3][1]->Fill(mRcRapidityM,float(mNHitsFitM)/float(mNHitsMaxM),mWeight*PhiWeight*DataRcWeightVal);
          //dedxrapidity[3][1]      ->Fill(mRcRapidityM,mDEdxM,mWeight*PhiWeight*DataRcWeightVal);
          //dcarapidity[3][1]       ->Fill(mRcRapidityM,mDcaM,mWeight*PhiWeight*DataRcWeightVal);
  
          //nhitsfiteta[3][0]  ->Fill(mRcEtaP,mNHitsFitP,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsmaxeta[3][0]  ->Fill(mRcEtaP,mNHitsMaxP,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsratioeta[3][0]->Fill(mRcEtaP,float(mNHitsFitP)/float(mNHitsMaxP),mWeight*PhiWeight*DataRcWeightVal);
          //dedxeta[3][0]      ->Fill(mRcEtaP,mDEdxP,mWeight*PhiWeight*DataRcWeightVal);
          //dcaeta[3][0]       ->Fill(mRcEtaP,mDcaP,mWeight*PhiWeight*DataRcWeightVal);
  
          //nhitsfiteta[3][1]  ->Fill(mRcEtaM,mNHitsFitM,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsmaxeta[3][1]  ->Fill(mRcEtaM,mNHitsMaxM,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsratioeta[3][1]->Fill(mRcEtaM,float(mNHitsFitM)/float(mNHitsMaxM),mWeight*PhiWeight*DataRcWeightVal);
          //dedxeta[3][1]      ->Fill(mRcEtaM,mDEdxM,mWeight*PhiWeight*DataRcWeightVal);
          //dcaeta[3][1]       ->Fill(mRcEtaM,mDcaM,mWeight*PhiWeight*DataRcWeightVal);
  
          //nhitsfitphi[3][0]  ->Fill(mRcPhiP,mNHitsFitP,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsmaxphi[3][0]  ->Fill(mRcPhiP,mNHitsMaxP,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsratiophi[3][0]->Fill(mRcPhiP,float(mNHitsFitP)/float(mNHitsMaxP),mWeight*PhiWeight*DataRcWeightVal);
          //dedxphi[3][0]      ->Fill(mRcPhiP,mDEdxP,mWeight*PhiWeight*DataRcWeightVal);
          //dcaphi[3][0]       ->Fill(mRcPhiP,mDcaP,mWeight*PhiWeight*DataRcWeightVal);
  
          //nhitsfitphi[3][1]  ->Fill(mRcPhiM,mNHitsFitM,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsmaxphi[3][1]  ->Fill(mRcPhiM,mNHitsMaxM,mWeight*PhiWeight*DataRcWeightVal); 
          //nhitsratiophi[3][1]->Fill(mRcPhiM,float(mNHitsFitM)/float(mNHitsMaxM),mWeight*PhiWeight*DataRcWeightVal);
          //dedxphi[3][1]      ->Fill(mRcPhiM,mDEdxM,mWeight*PhiWeight*DataRcWeightVal);
          //dcaphi[3][1]       ->Fill(mRcPhiM,mDcaM,mWeight*PhiWeight*DataRcWeightVal);

        //}  
      }
    }

    TH1F *CosThetaStarPhiRcDiff[vmsa::tofrebinpttotal][vmsa::tofrebinytotal];
    TH1F *CosThetaStarPhiRcDiffAfter[vmsa::tofrebinpttotal][vmsa::tofrebinytotal];
    TH1F *CosThetaStarPhiRcDiffRatio[vmsa::tofrebinpttotal][vmsa::tofrebinytotal];
    TH1F *CosThetaStarPhiRcPt[vmsa::tofrebinytotal];
    TH1F *CosThetaStarPhiRcPtAfter[vmsa::tofrebinytotal];
    TGraphAsymmErrors *gRcPt = new TGraphAsymmErrors();
    TGraphAsymmErrors *gRcPtAfter = new TGraphAsymmErrors();
    double rhoPhiRcPt[vmsa::tofrebinytotal];
    double rhoPhiRcPtAfter[vmsa::tofrebinytotal];
    double rhoErrPhiRcPt[vmsa::tofrebinytotal];
    double rhoErrPhiRcPtAfter[vmsa::tofrebinytotal];
    TH1F *CosThetaStarPhiRcPtRatio[vmsa::tofrebinytotal];
    TH1F *CosThetaStarPhiRcPtOpt[vmsa::tofrebinpttotal][3];
    TH1F *CosThetaStarPhiRcPtOptAfter[vmsa::tofrebinpttotal][3];
    TGraphAsymmErrors *gRcPtOpt[3];
    TGraphAsymmErrors *gRcPtOptAfter[3];
    for(int iy = 0; iy < 3; iy++)
    {
      gRcPtOpt[iy]      = new TGraphAsymmErrors();
      gRcPtOptAfter[iy] = new TGraphAsymmErrors();
    }
    double rhoPhiRcPtOpt[vmsa::tofrebinpttotal][3];
    double rhoPhiRcPtOptAfter[vmsa::tofrebinpttotal][3];
    double rhoErrPhiRcPtOpt[vmsa::tofrebinpttotal][3];
    double rhoErrPhiRcPtOptAfter[vmsa::tofrebinpttotal][3];
    TH1F *CosThetaStarPhiRcPtOptRatio[vmsa::tofrebinpttotal][3];
    TH1F *CosThetaStarPhiRcY[vmsa::tofrebinpttotal];
    TH1F *CosThetaStarPhiRcYAfter[vmsa::tofrebinpttotal];
    TGraphAsymmErrors *gRcY      = new TGraphAsymmErrors();
    TGraphAsymmErrors *gRcYAfter = new TGraphAsymmErrors();
    double rhoPhiRcY[vmsa::tofrebinpttotal];
    double rhoPhiRcYAfter[vmsa::tofrebinpttotal];
    double rhoErrPhiRcY[vmsa::tofrebinpttotal];
    double rhoErrPhiRcYAfter[vmsa::tofrebinpttotal];
    TH1F *CosThetaStarPhiRcYRatio[vmsa::tofrebinpttotal];

    int ylow[3] = {1,4,12};
    int yhigh[3] = {3,11,14};
    for(int ipt = 0; ipt < vmsa::tofrebinpttotal; ipt++)
    {
      for(int iy = 0; iy < 3; iy++)
      {
        string Hist = Form("costhetaopt_pt_%d_y_%d",ipt,iy);
        CosThetaStarPhiRcPtOpt[ipt][iy] = (TH1F*) CosThetaStarPhiRc->ProjectionZ(Hist.c_str(),ipt+1,ipt+1,ylow[iy],yhigh[iy],"e");
        string HistAfter = Form("costhetaoptafter_pt_%d_y_%d",ipt,iy);
        CosThetaStarPhiRcPtOptAfter[ipt][iy] = (TH1F*) CosThetaStarPhiRcAfter->ProjectionZ(HistAfter.c_str(),ipt+1,ipt+1,ylow[iy],yhigh[iy],"e");

        CosThetaStarPhiRcPtOptRatio[ipt][iy] = (TH1F*) CosThetaStarPhiRcPtOptAfter[ipt][iy]->Clone();
        CosThetaStarPhiRcPtOptRatio[ipt][iy]->Divide(CosThetaStarPhiRcPtOptAfter[ipt][iy],CosThetaStarPhiRcPtOpt[ipt][iy],1,1,"B");

        double ptmean = (vmsa::tofrebinptval[ipt]+vmsa::tofrebinptval[ipt+1])/2.0;

        TF1* spin = new TF1("rho",SpinDensity,-1,1,2);
        spin->SetParameter(0,0.33);
        spin->SetParameter(1,CosThetaStarPhiRcPtOpt[ipt][iy]->GetMaximum());
        CosThetaStarPhiRcPtOpt[ipt][iy]->Fit(spin,"NMRI");
       
        double rho00 = spin->GetParameter(0);        
        double rho00err = spin->GetParError(0);        

        gRcPtOpt[iy]->SetPoint(ipt,ptmean,rho00);
        gRcPtOpt[iy]->SetPointError(ipt,0.0,0.0,rho00err,rho00err);

        TF1* spinA = new TF1("rhoA",SpinDensity,-1,1,2);
        spinA->SetParameter(0,0.33);
        spinA->SetParameter(1,CosThetaStarPhiRcPtOptAfter[ipt][iy]->GetMaximum());
        CosThetaStarPhiRcPtOptAfter[ipt][iy]->Fit(spinA,"NMRI");
       
        double rho00A = spinA->GetParameter(0);        
        double rho00errA = spinA->GetParError(0);        

        gRcPtOptAfter[iy]->SetPoint(ipt,ptmean,rho00A);
        gRcPtOptAfter[iy]->SetPointError(ipt,0.0,0.0,rho00errA,rho00errA);
        
        
      } 
    }   

    for(int ipt = 0; ipt < vmsa::tofrebinpttotal; ipt++)
    {
      for(int iy = 0; iy < vmsa::tofrebinytotal; iy++)
      {
        string Hist = Form("costheta_pt_%d_y_%d",ipt,iy);
        CosThetaStarPhiRcDiff[ipt][iy] = (TH1F*) CosThetaStarPhiRc->ProjectionZ(Hist.c_str(),ipt+1,ipt+1,iy+1,iy+1,"e");
        string HistAfter = Form("costhetaafter_pt_%d_y_%d",ipt,iy);
        CosThetaStarPhiRcDiffAfter[ipt][iy] = (TH1F*) CosThetaStarPhiRcAfter->ProjectionZ(HistAfter.c_str(),ipt+1,ipt+1,iy+1,iy+1,"e");

        CosThetaStarPhiRcDiffRatio[ipt][iy] = (TH1F*) CosThetaStarPhiRcDiffAfter[ipt][iy]->Clone();
        CosThetaStarPhiRcDiffRatio[ipt][iy]->Divide(CosThetaStarPhiRcDiffAfter[ipt][iy],CosThetaStarPhiRcDiff[ipt][iy],1,1,"B");
      } 
    }   
    for(int ipt = 0; ipt < vmsa::tofrebinpttotal; ipt++)
    {
      string Hist = Form("costheta_pt_%d",ipt);
      CosThetaStarPhiRcY[ipt] = (TH1F*) CosThetaStarPhiRc->ProjectionZ(Hist.c_str(),ipt+1,ipt+1,1,vmsa::tofrebinytotal,"e");
      string HistAfter = Form("costhetaafter_pt_%d",ipt);
      CosThetaStarPhiRcYAfter[ipt] = (TH1F*) CosThetaStarPhiRcAfter->ProjectionZ(HistAfter.c_str(),ipt+1,ipt+1,1,vmsa::tofrebinytotal,"e");

      CosThetaStarPhiRcYRatio[ipt] = (TH1F*) CosThetaStarPhiRcYAfter[ipt]->Clone();
      CosThetaStarPhiRcYRatio[ipt]->Divide(CosThetaStarPhiRcYAfter[ipt],CosThetaStarPhiRcY[ipt],1,1,"B");

      double ptmean = (vmsa::tofrebinptval[ipt]+vmsa::tofrebinptval[ipt+1])/2.0;

      TF1* spin = new TF1("rho",SpinDensity,-1,1,2);
      spin->SetParameter(0,0.33);
      spin->SetParameter(1,CosThetaStarPhiRcY[ipt]->GetMaximum());
      CosThetaStarPhiRcY[ipt]->Fit(spin,"NMRI");
      
      double rho00 = spin->GetParameter(0);        
      double rho00err = spin->GetParError(0);        

      gRcY->SetPoint(ipt,ptmean,rho00);
      gRcY->SetPointError(ipt,0.0,0.0,rho00err,rho00err);

      TF1* spinA = new TF1("rhoA",SpinDensity,-1,1,2);
      spinA->SetParameter(0,0.33);
      spinA->SetParameter(1,CosThetaStarPhiRcYAfter[ipt]->GetMaximum());
      CosThetaStarPhiRcYAfter[ipt]->Fit(spinA,"NMRI");
      
      double rho00A = spinA->GetParameter(0);        
      double rho00errA = spinA->GetParError(0);        

      gRcYAfter->SetPoint(ipt,ptmean,rho00A);
      gRcYAfter->SetPointError(ipt,0.0,0.0,rho00errA,rho00errA);

    } 
    for(int iy = 0; iy < vmsa::tofrebinytotal; iy++)
    { 
      string Hist = Form("costheta_y_%d",iy);
      CosThetaStarPhiRcPt[iy] = (TH1F*) CosThetaStarPhiRc->ProjectionZ(Hist.c_str(),1,vmsa::tofrebinpttotal,iy+1,iy+1,"e");
      string HistAfter = Form("costhetaafter_y_%d",iy);
      CosThetaStarPhiRcPtAfter[iy] = (TH1F*) CosThetaStarPhiRcAfter->ProjectionZ(HistAfter.c_str(),1,vmsa::tofrebinpttotal,iy+1,iy+1,"e");

      CosThetaStarPhiRcPtRatio[iy] = (TH1F*) CosThetaStarPhiRcPtAfter[iy]->Clone();
      CosThetaStarPhiRcPtRatio[iy]->Divide(CosThetaStarPhiRcPtAfter[iy],CosThetaStarPhiRcPt[iy],1,1,"B");
 
       
      double ymean = (vmsa::tofrebinyval[iy]+vmsa::tofrebinyval[iy+1])/2.0;

      TF1* spin = new TF1("rho",SpinDensity,-1,1,2);
      spin->SetParameter(0,0.33);
      spin->SetParameter(1,CosThetaStarPhiRcPt[iy]->GetMaximum());
      CosThetaStarPhiRcPt[iy]->Fit(spin,"NMRI");
      
      double rho00 = spin->GetParameter(0);        
      double rho00err = spin->GetParError(0);        

      gRcPt->SetPoint(iy,ymean,rho00);
      gRcPt->SetPointError(iy,0.0,0.0,rho00err,rho00err);

      TF1* spinA = new TF1("rhoA",SpinDensity,-1,1,2);
      spinA->SetParameter(0,0.33);
      spinA->SetParameter(1,CosThetaStarPhiRcPtAfter[iy]->GetMaximum());
      CosThetaStarPhiRcPtAfter[iy]->Fit(spinA,"NMRI");
      
      double rho00A = spinA->GetParameter(0);        
      double rho00errA = spinA->GetParError(0);        

      gRcPtAfter->SetPoint(iy,ymean,rho00A);
      gRcPtAfter->SetPointError(iy,0.0,0.0,rho00errA,rho00errA);
    }
  
    TCanvas *ccos = new TCanvas("ccos", "cos", 1200, 1200);
    ccos->Divide(4,4);
    for(int i = 0; i < 16; i++)
    {
      ccos->cd(i+1);
      ccos->cd(i+1)->SetLeftMargin(0.12);
      ccos->cd(i+1)->SetRightMargin(0.15);
      ccos->cd(i+1)->SetBottomMargin(0.12);
      ccos->cd(i+1)->SetTicks(1,1);
      ccos->cd(i+1)->SetGrid(0,0); 
    } 

    string outputname = Form("figures/KaonTTrees/%sCosThetaDistributions.pdf",folderopt.c_str());
    string outputstart = Form("%s[",outputname.c_str());
    string outputstop = Form("%s]",outputname.c_str());

    ccos->Print(outputstart.c_str());

    for(int ipt = 0; ipt < vmsa::tofrebinpttotal; ipt++)
    {
      for(int iy = 0; iy < vmsa::tofrebinytotal; iy++)
      {
        ccos->cd(iy+1);

        double max = CosThetaStarPhiRcDiff[ipt][iy]->GetMaximum();
        double min = CosThetaStarPhiRcDiff[ipt][iy]->GetMinimum();

        double maxAfter = CosThetaStarPhiRcDiffAfter[ipt][iy]->GetMaximum();
        double minAfter = CosThetaStarPhiRcDiffAfter[ipt][iy]->GetMinimum();

        if(maxAfter > max) max = maxAfter;
        if(minAfter < min) min = minAfter;

        CosThetaStarPhiRcDiff[ipt][iy]->SetTitle(Form("%.2f<p_{T}<%.2f GeV/c, %.2f<y<%.2f",vmsa::tofrebinptval[ipt],vmsa::tofrebinptval[ipt+1],vmsa::tofrebinyval[iy],vmsa::tofrebinyval[iy+1]));   
        CosThetaStarPhiRcDiff[ipt][iy]->GetXaxis()->SetTitle("cos(#theta*)");
        CosThetaStarPhiRcDiff[ipt][iy]->GetYaxis()->SetTitle("Counts");
        CosThetaStarPhiRcDiff[ipt][iy]->SetMarkerStyle(20);
        CosThetaStarPhiRcDiffAfter[ipt][iy]->SetMarkerStyle(24);
        CosThetaStarPhiRcDiff[ipt][iy]->SetMarkerColor(kBlue);
        CosThetaStarPhiRcDiffAfter[ipt][iy]->SetMarkerColor(kOrange+7);
        CosThetaStarPhiRcDiff[ipt][iy]->SetLineColor(kBlue);
        CosThetaStarPhiRcDiffAfter[ipt][iy]->SetLineColor(kOrange+7);
        CosThetaStarPhiRcDiff[ipt][iy]->GetYaxis()->SetRangeUser(min*0.8,max*1.2);
        CosThetaStarPhiRcDiff[ipt][iy]->Draw("pE");
        CosThetaStarPhiRcDiffAfter[ipt][iy]->Draw("pE same");
     
        TLegend *legcos = new TLegend(0.7,0.7,0.85,0.85);
        legcos->AddEntry(CosThetaStarPhiRcDiff[ipt][iy],"No Weight","p");
        legcos->AddEntry(CosThetaStarPhiRcDiffAfter[ipt][iy],"Data/RC Weight","p");
        legcos->Draw("same");
      }
      ccos->Update();
      ccos->Print(outputname.c_str());
    }
    ccos->Print(outputstop.c_str());

    for(int i = 0; i < 16; i++)  ccos->cd(i+1)->Clear();

    TCanvas *ccosopt = new TCanvas("ccosopt", "cos", 900, 600);
    ccosopt->Divide(3,2);
    for(int i = 0; i < 6; i++)
    {
      ccosopt->cd(i+1);
      ccosopt->cd(i+1)->SetLeftMargin(0.12);
      ccosopt->cd(i+1)->SetRightMargin(0.15);
      ccosopt->cd(i+1)->SetBottomMargin(0.12);
      ccosopt->cd(i+1)->SetTicks(1,1);
      ccosopt->cd(i+1)->SetGrid(0,0); 
    } 

    outputname = Form("figures/KaonTTrees/%sCosThetaDistributions_PtYselections.pdf",folderopt.c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop = Form("%s]",outputname.c_str());

    ccosopt->Print(outputstart.c_str());

    float ylowval[3] = {-1.0,-0.5,0.5};
    float yhighval[3] = {-0.5,0.5,1.0};

    for(int ipt = 0; ipt < vmsa::tofrebinpttotal; ipt++)
    {
      for(int iy = 0; iy < 3; iy++)
      {
        ccosopt->cd(iy+1);

        double max = CosThetaStarPhiRcPtOpt[ipt][iy]->GetMaximum();
        double min = CosThetaStarPhiRcPtOpt[ipt][iy]->GetMinimum();

        double maxAfter = CosThetaStarPhiRcPtOptAfter[ipt][iy]->GetMaximum();
        double minAfter = CosThetaStarPhiRcPtOptAfter[ipt][iy]->GetMinimum();

        if(maxAfter > max) max = maxAfter;
        if(minAfter < min) min = minAfter;

        CosThetaStarPhiRcPtOpt[ipt][iy]->SetTitle(Form("%.2f<p_{T}<%.2f GeV/c, %.2f<y<%.2f",vmsa::tofrebinptval[ipt],vmsa::tofrebinptval[ipt+1],ylowval[iy],yhighval[iy]));   
        CosThetaStarPhiRcPtOpt[ipt][iy]->GetXaxis()->SetTitle("cos(#theta*)");
        CosThetaStarPhiRcPtOpt[ipt][iy]->GetYaxis()->SetTitle("Counts");
        CosThetaStarPhiRcPtOpt[ipt][iy]->SetMarkerStyle(20);
        CosThetaStarPhiRcPtOptAfter[ipt][iy]->SetMarkerStyle(24);
        CosThetaStarPhiRcPtOpt[ipt][iy]->SetMarkerColor(kBlue);
        CosThetaStarPhiRcPtOptAfter[ipt][iy]->SetMarkerColor(kOrange+7);
        CosThetaStarPhiRcPtOpt[ipt][iy]->SetLineColor(kBlue);
        CosThetaStarPhiRcPtOptAfter[ipt][iy]->SetLineColor(kOrange+7);
        CosThetaStarPhiRcPtOpt[ipt][iy]->GetYaxis()->SetRangeUser(min*0.8,max*1.2);
        CosThetaStarPhiRcPtOpt[ipt][iy]->Draw("pE");
        CosThetaStarPhiRcPtOptAfter[ipt][iy]->Draw("pE same");
     
        TLegend *legcos = new TLegend(0.7,0.7,0.85,0.85);
        legcos->AddEntry(CosThetaStarPhiRcPtOpt[ipt][iy],"No Weight","p");
        legcos->AddEntry(CosThetaStarPhiRcPtOptAfter[ipt][iy],"Data/RC Weight","p");
        legcos->Draw("same");

        ccosopt->cd(iy+4);
        CosThetaStarPhiRcPtOptRatio[ipt][iy]->SetTitle(Form("(Data/RC Reweight)/(No Reweight) %.2f<p_{T}<%.2f GeV/c, %.2f<y<%.2f",vmsa::tofrebinptval[ipt],vmsa::tofrebinptval[ipt+1],ylowval[iy],yhighval[iy]));   
        CosThetaStarPhiRcPtOptRatio[ipt][iy]->GetXaxis()->SetTitle("cos(#theta*)");
        CosThetaStarPhiRcPtOptRatio[ipt][iy]->GetYaxis()->SetTitle("(Data/RC Reweight)/(No Reweight)");
        CosThetaStarPhiRcPtOptRatio[ipt][iy]->Draw("pE");

      }
      ccosopt->Update();
      ccosopt->Print(outputname.c_str());
    }
    ccosopt->Print(outputstop.c_str());

    for(int i = 0; i < 6; i++)  ccosopt->cd(i+1)->Clear();

    outputname = Form("figures/KaonTTrees/%sCosThetaDistributions_rho00.pdf",folderopt.c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop = Form("%s]",outputname.c_str());

    ccosopt->Print(outputstart.c_str());

    //float ylowval[3] = {-1.0,-0.5,0.5};
    //float yhighval[3] = {-0.5,0.5,1.0};

    for(int iy = 0; iy < 3; iy++)
    {
      ccosopt->cd(iy+1);

      gRcPtOpt[iy]->SetTitle(Form("%.2f<y<%.2f",ylowval[iy],yhighval[iy]));   
      gRcPtOpt[iy]->GetXaxis()->SetTitle("p_{T} GeV/c");
      gRcPtOpt[iy]->GetYaxis()->SetTitle("#rho_{00}");
      gRcPtOpt[iy]->SetMarkerStyle(20);
      gRcPtOptAfter[iy]->SetMarkerStyle(24);
      gRcPtOpt[iy]->SetMarkerColor(kBlue);
      gRcPtOptAfter[iy]->SetMarkerColor(kOrange+7);
      gRcPtOpt[iy]->SetLineColor(kBlue);
      gRcPtOptAfter[iy]->SetLineColor(kOrange+7);
      gRcPtOpt[iy]->GetYaxis()->SetRangeUser(0.2,0.5);
      gRcPtOpt[iy]->Draw("ApE");
      gRcPtOptAfter[iy]->Draw("pE same");
    
      TLegend *legcos = new TLegend(0.7,0.7,0.85,0.85);
      legcos->AddEntry(gRcPtOpt[iy],"No Weight","p");
      legcos->AddEntry(gRcPtOptAfter[iy],"Data/RC Weight","p");
      legcos->Draw("same");
    }

    ccosopt->cd(4);

    gRcPt->SetTitle(Form("1.2<p_{T}<4.2"));   
    gRcPt->GetXaxis()->SetTitle("y");
    gRcPt->GetYaxis()->SetTitle("#rho_{00}");
    gRcPt->SetMarkerStyle(20);
    gRcPtAfter->SetMarkerStyle(24);
    gRcPt->SetMarkerColor(kBlue);
    gRcPtAfter->SetMarkerColor(kOrange+7);
    gRcPt->SetLineColor(kBlue);
    gRcPtAfter->SetLineColor(kOrange+7);
    gRcPt->GetYaxis()->SetRangeUser(0.2,0.5);
    gRcPt->Draw("ApE");
    gRcPtAfter->Draw("pE same");

    ccosopt->cd(5);

    gRcY->SetTitle(Form("-1.0<y<1.0"));   
    gRcY->GetXaxis()->SetTitle("p_{T} GeV/c");
    gRcY->GetYaxis()->SetTitle("#rho_{00}");
    gRcY->SetMarkerStyle(20);
    gRcYAfter->SetMarkerStyle(24);
    gRcY->SetMarkerColor(kBlue);
    gRcYAfter->SetMarkerColor(kOrange+7);
    gRcY->SetLineColor(kBlue);
    gRcYAfter->SetLineColor(kOrange+7);
    gRcY->GetYaxis()->SetRangeUser(0.2,0.5);
    gRcY->Draw("ApE");
    gRcYAfter->Draw("pE same");


    ccosopt->Update();
    ccosopt->Print(outputname.c_str());
    
    ccosopt->Print(outputstop.c_str());

//    for(int i = 0; i < 16; i++)  ccos->cd(i+1)->Clear();


    outputname = Form("figures/KaonTTrees/%sCosThetaDistributions_Pt.pdf",folderopt.c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop = Form("%s]",outputname.c_str());

    ccos->Print(outputstart.c_str());

    for(int iy = 0; iy < vmsa::tofrebinytotal; iy++)
    {
      ccos->cd(iy+1);

      double max = CosThetaStarPhiRcPt[iy]->GetMaximum();
      double min = CosThetaStarPhiRcPt[iy]->GetMinimum();

      double maxAfter = CosThetaStarPhiRcPtAfter[iy]->GetMaximum();
      double minAfter = CosThetaStarPhiRcPtAfter[iy]->GetMinimum();

      if(maxAfter > max) max = maxAfter;
      if(minAfter < min) min = minAfter;

      CosThetaStarPhiRcPt[iy]->SetTitle(Form("%.2f<p_{T}<%.2f GeV/c, %.2f<y<%.2f",vmsa::tofrebinptval[0],vmsa::tofrebinptval[vmsa::tofrebinpttotal-1],vmsa::tofrebinyval[iy],vmsa::tofrebinyval[iy+1]));   
      CosThetaStarPhiRcPt[iy]->GetXaxis()->SetTitle("cos(#theta*)");
      CosThetaStarPhiRcPt[iy]->GetYaxis()->SetTitle("Counts");
      CosThetaStarPhiRcPt[iy]->SetMarkerStyle(20);
      CosThetaStarPhiRcPtAfter[iy]->SetMarkerStyle(24);
      CosThetaStarPhiRcPt[iy]->SetMarkerColor(kBlue);
      CosThetaStarPhiRcPtAfter[iy]->SetMarkerColor(kOrange+7);
      CosThetaStarPhiRcPt[iy]->SetLineColor(kBlue);
      CosThetaStarPhiRcPtAfter[iy]->SetLineColor(kOrange+7);
      CosThetaStarPhiRcPt[iy]->GetYaxis()->SetRangeUser(min*0.8,max*1.2);
      CosThetaStarPhiRcPt[iy]->Draw("pE");
      CosThetaStarPhiRcPtAfter[iy]->Draw("pE same");
    
      TLegend *legcos = new TLegend(0.7,0.7,0.85,0.85);
      legcos->AddEntry(CosThetaStarPhiRcPt[iy],"No Weight","p");
      legcos->AddEntry(CosThetaStarPhiRcPtAfter[iy],"Data/RC Weight","p");
      legcos->Draw("same");
    }
    ccos->Update();
    ccos->Print(outputname.c_str());
    ccos->Print(outputstop.c_str());
    for(int i = 0; i < 16; i++)  ccos->cd(i+1)->Clear();

    outputname = Form("figures/KaonTTrees/%sCosThetaDistributions_Y.pdf",folderopt.c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop = Form("%s]",outputname.c_str());

    ccos->Print(outputstart.c_str());

    for(int ipt = 0; ipt < vmsa::tofrebinpttotal; ipt++)
    {
      ccos->cd(ipt+1);

      double max = CosThetaStarPhiRcY[ipt]->GetMaximum();
      double min = CosThetaStarPhiRcY[ipt]->GetMinimum();

      double maxAfter = CosThetaStarPhiRcYAfter[ipt]->GetMaximum();
      double minAfter = CosThetaStarPhiRcYAfter[ipt]->GetMinimum();

      if(maxAfter > max) max = maxAfter;
      if(minAfter < min) min = minAfter;

      CosThetaStarPhiRcY[ipt]->SetTitle(Form("%.2f<p_{T}<%.2f GeV/c, 1.0<y<1.0",vmsa::tofrebinptval[ipt],vmsa::tofrebinptval[ipt+1]));   
      CosThetaStarPhiRcY[ipt]->GetXaxis()->SetTitle("cos(#theta*)");
      CosThetaStarPhiRcY[ipt]->GetYaxis()->SetTitle("Counts");
      CosThetaStarPhiRcY[ipt]->SetMarkerStyle(20);
      CosThetaStarPhiRcYAfter[ipt]->SetMarkerStyle(24);
      CosThetaStarPhiRcY[ipt]->SetMarkerColor(kBlue);
      CosThetaStarPhiRcYAfter[ipt]->SetMarkerColor(kOrange+7);
      CosThetaStarPhiRcY[ipt]->SetLineColor(kBlue);
      CosThetaStarPhiRcYAfter[ipt]->SetLineColor(kOrange+7);
      CosThetaStarPhiRcY[ipt]->GetYaxis()->SetRangeUser(min*0.8,max*1.2);
      CosThetaStarPhiRcY[ipt]->Draw("pE");
      CosThetaStarPhiRcYAfter[ipt]->Draw("pE same");
     
      TLegend *legcos = new TLegend(0.7,0.7,0.85,0.85);
      legcos->AddEntry(CosThetaStarPhiRcY[ipt],"No Weight","p");
      legcos->AddEntry(CosThetaStarPhiRcYAfter[ipt],"Data/RC Weight","p");
      legcos->Draw("same");
    }
    ccos->Update();
    ccos->Print(outputname.c_str());
    ccos->Print(outputstop.c_str());

    for(int i = 0; i < 16; i++)  ccos->cd(i+1)->Clear();
    outputname = Form("figures/KaonTTrees/%sCosThetaRatios.pdf",folderopt.c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop = Form("%s]",outputname.c_str());

    ccos->Print(outputstart.c_str());

    for(int ipt = 0; ipt < vmsa::tofrebinpttotal; ipt++)
    {
      for(int iy = 0; iy < vmsa::tofrebinytotal; iy++)
      {
        ccos->cd(iy+1);

        CosThetaStarPhiRcDiffRatio[ipt][iy]->SetTitle(Form("(Data/RC Reweight)/(No Reweight) %.2f<p_{T}<%.2f GeV/c, %.2f<y<%.2f",vmsa::tofrebinptval[ipt],vmsa::tofrebinptval[ipt+1],vmsa::tofrebinyval[iy],vmsa::tofrebinyval[iy+1]));   
        CosThetaStarPhiRcDiffRatio[ipt][iy]->GetXaxis()->SetTitle("cos(#theta*)");
        CosThetaStarPhiRcDiffRatio[ipt][iy]->GetYaxis()->SetTitle("(Data/RC Reweight)/(No Reweight)");
        CosThetaStarPhiRcDiffRatio[ipt][iy]->Draw("pE");
      }
      ccos->Update();
      ccos->Print(outputname.c_str());
    }
    ccos->Print(outputstop.c_str());

    for(int i = 0; i < 16; i++)  ccos->cd(i+1)->Clear();
    outputname = Form("figures/KaonTTrees/%sCosThetaRatios_Y.pdf",folderopt.c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop = Form("%s]",outputname.c_str());

    ccos->Print(outputstart.c_str());

    for(int ipt = 0; ipt < vmsa::tofrebinpttotal; ipt++)
    {
      ccos->cd(ipt+1);

      CosThetaStarPhiRcYRatio[ipt]->SetTitle(Form("(Data/RC Reweight)/(No Reweight) %.2f<p_{T}<%.2f GeV/c, -1.0<y<1.0",vmsa::tofrebinptval[ipt],vmsa::tofrebinptval[ipt+1]));   
      CosThetaStarPhiRcYRatio[ipt]->GetXaxis()->SetTitle("cos(#theta*)");
      CosThetaStarPhiRcYRatio[ipt]->GetYaxis()->SetTitle("(Data/RC Reweight)/(No Reweight)");
      CosThetaStarPhiRcYRatio[ipt]->Draw("pE");
    }
    ccos->Update();
    ccos->Print(outputname.c_str());
    ccos->Print(outputstop.c_str());

    for(int i = 0; i < 16; i++)  ccos->cd(i+1)->Clear();
    outputname = Form("figures/KaonTTrees/%sCosThetaRatios_Pt.pdf",folderopt.c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop = Form("%s]",outputname.c_str());

    ccos->Print(outputstart.c_str());

    for(int iy = 0; iy < vmsa::tofrebinytotal; iy++)
    {
      ccos->cd(iy+1);

      CosThetaStarPhiRcPtRatio[iy]->SetTitle(Form("(Data/RC Reweight)/(No Reweight) %.2f<p_{T}<%.2f GeV/c, %.2f<y<%.2f",vmsa::tofrebinptval[0],vmsa::tofrebinptval[vmsa::tofrebinpttotal-1],vmsa::tofrebinyval[iy],vmsa::tofrebinyval[iy+1]));   
      CosThetaStarPhiRcPtRatio[iy]->GetXaxis()->SetTitle("cos(#theta*)");
      CosThetaStarPhiRcPtRatio[iy]->GetYaxis()->SetTitle("(Data/RC Reweight)/(No Reweight)");
      CosThetaStarPhiRcPtRatio[iy]->Draw("pE");
    }
    ccos->Update();
    ccos->Print(outputname.c_str());
    ccos->Print(outputstop.c_str());
    
    for(int i = 0; i < 16; i++)  ccos->cd(i+1)->Clear();


    TH2F *ptrapidityPhiDataRcAfter;

    //double integralPhiData = ptrapidityPhiData->Integral(1,vmsa::tofrebinpttotal,1,vmsa::tofrebinytotal); 
    double integralPhiRcAfter = ptrapidityPhiRcAfter->Integral(1,vmsa::tofrebinpttotal,1,vmsa::tofrebinytotal); 
    double DataRcPhiAfter = integralPhiData/integralPhiRcAfter;
    
    ptrapidityPhiDataRcAfter = (TH2F*) ptrapidityPhiData->Clone();
    ptrapidityPhiRcAfter->Scale(DataRcPhiAfter);
    ptrapidityPhiDataRcAfter->Divide(ptrapidityPhiRcAfter);

    TH1F *ptPhiRcAfter   = (TH1F*) ptrapidityPhiRcAfter->ProjectionX("ptRcAfter",0,-1);  
    TH1F *yPhiRcAfter    = (TH1F*) ptrapidityPhiRcAfter->ProjectionY("yRcAfter",0,-1);  

    TH1F *ptPhiDataRcAfter = (TH1F*) ptPhiData->Clone();
    ptPhiDataRcAfter->Divide(ptPhiRcAfter);
    TH1F *yPhiDataRcAfter = (TH1F*) yPhiData->Clone();
    yPhiDataRcAfter->Divide(yPhiRcAfter);

    //TCanvas *ct = new TCanvas("ct", "ct", 1200, 400);
    //ct->Divide(3,1);
    //for(int i = 0; i < 3; i++)
    //{
    //  ct->cd(i+1);
    //  ct->cd(i+1)->SetLeftMargin(0.15);
    //  ct->cd(i+1)->SetRightMargin(0.15);
    //  ct->cd(i+1)->SetBottomMargin(0.15);
    //  ct->cd(i+1)->SetTicks(1,1);
    //  ct->cd(i+1)->SetGrid(0,0); 
    //}

    ct->cd(1);
    //ct->cd(1)->SetLogz();
    ptrapidityPhiData->SetTitle("Data"); 
    ptrapidityPhiData->GetXaxis()->SetTitle("p_{T} GeV/c"); 
    ptrapidityPhiData->GetYaxis()->SetTitle("y"); 
    ptrapidityPhiData->Draw("Colz"); 
    ct->cd(2);
    //ct->cd(2)->SetLogz();
    ptrapidityPhiRcAfter->SetTitle("Embedding"); 
    ptrapidityPhiRcAfter->GetXaxis()->SetTitle("p_{T} GeV/c"); 
    ptrapidityPhiRcAfter->GetYaxis()->SetTitle("y"); 
    ptrapidityPhiRcAfter->Draw("Colz"); 
    ct->cd(3);
    ct->cd(3)->SetLogy(0);
    ptrapidityPhiDataRcAfter->SetTitle("Data/Embedding"); 
    ptrapidityPhiDataRcAfter->GetXaxis()->SetTitle("p_{T} GeV/c"); 
    ptrapidityPhiDataRcAfter->GetYaxis()->SetTitle("y"); 
    ptrapidityPhiDataRcAfter->Draw("Colz"); 

    ct->cd(4);
    ct->cd(4)->SetLogy();
    ptPhiData->SetTitle("After Data/RC Reweight"); 
    ptPhiData->GetXaxis()->SetTitle("p_{T} GeV/c"); 
    ptPhiData->GetYaxis()->SetTitle("Counts"); 
    ptPhiData->SetMarkerColor(kBlue); 
    ptPhiData->SetLineColor(kBlue); 
    ptPhiData->SetMarkerStyle(20); 
    ptPhiData->Draw("pE"); 
    ptPhiRcAfter->SetMarkerColor(kOrange+7); 
    ptPhiRcAfter->SetLineColor(kOrange+7); 
    ptPhiRcAfter->SetMarkerStyle(20); 
    ptPhiRcAfter->Draw("pE same"); 
    legpt->Draw("same"); 

    ct->cd(5);
    //ct->cd(5)->SetLogy();
    yPhiData->SetTitle("After Data/RC Reweight"); 
    yPhiData->GetXaxis()->SetTitle("y"); 
    yPhiData->GetYaxis()->SetTitle("Counts"); 
    yPhiData->SetMarkerColor(kBlue); 
    yPhiData->SetLineColor(kBlue); 
    yPhiData->SetMarkerStyle(20); 
    yPhiData->Draw("pE"); 
    yPhiRcAfter->SetMarkerColor(kOrange+7); 
    yPhiRcAfter->SetLineColor(kOrange+7); 
    yPhiRcAfter->SetMarkerStyle(20); 
    yPhiRcAfter->Draw("pE same"); 
    legy->Draw("same"); 

    ct->cd(7);
    ptPhiDataRcAfter->SetTitle("After Data/RC Reweight #phi-meson"); 
    ptPhiDataRcAfter->GetXaxis()->SetTitle("p_{T} GeV/c"); 
    ptPhiDataRcAfter->GetYaxis()->SetTitle("Data/RC"); 
    ptPhiDataRcAfter->Draw("pE");

    ct->cd(8);
    yPhiDataRcAfter->SetTitle("After Data/RC Reweight #phi-meson"); 
    yPhiDataRcAfter->GetXaxis()->SetTitle("y"); 
    yPhiDataRcAfter->GetYaxis()->SetTitle("Data/RC"); 
    yPhiDataRcAfter->Draw("pE");

    ct->SaveAs(Form("figures/KaonTTrees/%sPhiMeson_DataRCAndRatio_After.pdf",folderopt.c_str()));

    




    const int nset = 18;
    //                      y    y   y   pT
    Double_t axis[nset] = {   1,   1,  1,   1,  1,   1,  1,   0,   0,   0,   0,   0  ,0  ,0  ,0  ,0  ,0  ,0  };
    Double_t low[nset]  = {-1.5,-0.5,0.5,-1.0,0.9,-1.0,0.8, 1.0, 0.0, 2.0,   0.44,0.3,0.5,0.7,0.9,1.1,1.3,1.5};
    Double_t high[nset] = {-0.5, 0.5,1.5,-0.9,1.0,-0.8,1.0, 2.0, 0.7, ptmax, 1.5 ,0.5,0.7,0.9,1.1,1.3,1.5,1.7 };
 
    pt[4][2] = (TH1F*) ptyetaphibasic[4][2]->Projection(0);
    rapidity[4][2] = (TH1F*) ptyetaphibasic[4][2]->Projection(1);
    eta[4][2] = (TH1F*) ptyetaphibasic[4][2]->Projection(2);
    phi[4][2] = (TH1F*) ptyetaphibasic[4][2]->Projection(3);

    ptrapidity[4][2] = (TH2F*) ptyetaphibasic[4][2]->Projection(1,0);
    pteta[4][2] = (TH2F*) ptyetaphibasic[4][2]->Projection(2,0);
    ptphi[4][2] = (TH2F*) ptyetaphibasic[4][2]->Projection(3,0);
    rapidityeta[4][2] = (TH2F*) ptyetaphibasic[4][2]->Projection(2,1);
    rapidityphi[4][2] = (TH2F*) ptyetaphibasic[4][2]->Projection(3,1);
    etaphi[4][2] = (TH2F*) ptyetaphibasic[4][2]->Projection(3,2);

    for(int ic = 0; ic < 2; ic++)
    {
      ptyetaphibasic[2][ic] = (THnF*) ptyetaphibasic[0][ic]->Clone();
      ptyetaphibasic[2][ic]->Add(ptyetaphibasic[1][ic],-1.0);

      //pt[2][ic] = (TH1F*) pt[0][ic]->Clone();
      //pt[2][ic]->Add(pt[1][ic],-1.0);
      //rapidity[2][ic] = (TH1F*) rapidity[0][ic]->Clone();
      //rapidity[2][ic]->Add(rapidity[1][ic],-1.0);
      //eta[2][ic] = (TH1F*) eta[0][ic]->Clone();
      //eta[2][ic]->Add(eta[1][ic],-1.0);
      //phi[2][ic] = (TH1F*) phi[0][ic]->Clone();
      //phi[2][ic]->Add(phi[1][ic],-1.0);

      //ptrapidity[2][ic] = (TH2F*) ptrapidity[0][ic]->Clone();
      //ptrapidity[2][ic]->Add(ptrapidity[1][ic],-1.0);
      //pteta[2][ic] = (TH2F*) pteta[0][ic]->Clone();
      //pteta[2][ic]->Add(pteta[1][ic],-1.0);
      //ptphi[2][ic] = (TH2F*) ptphi[0][ic]->Clone();
      //ptphi[2][ic]->Add(ptphi[1][ic],-1.0);
      //rapidityeta[2][ic] = (TH2F*) rapidityeta[0][ic]->Clone();
      //rapidityeta[2][ic]->Add(rapidityeta[1][ic],-1.0);
      //rapidityphi[2][ic] = (TH2F*) rapidityphi[0][ic]->Clone();
      //rapidityphi[2][ic]->Add(rapidityphi[1][ic],-1.0);
      //etaphi[2][ic] = (TH2F*) etaphi[0][ic]->Clone();
      //etaphi[2][ic]->Add(etaphi[1][ic],-1.0);

      int numerator[2] = {2,3};
      int denominator[2] = {3,2};

      pt[4][ic] = (TH1F*) ptyetaphibasic[4][ic]->Projection(0);
      rapidity[4][ic] = (TH1F*) ptyetaphibasic[4][ic]->Projection(1);
      eta[4][ic] = (TH1F*) ptyetaphibasic[4][ic]->Projection(2);
      phi[4][ic] = (TH1F*) ptyetaphibasic[4][ic]->Projection(3);

      ptrapidity[4][ic] = (TH2F*) ptyetaphibasic[4][ic]->Projection(1,0);
      pteta[4][ic] = (TH2F*) ptyetaphibasic[4][ic]->Projection(2,0);
      ptphi[4][ic] = (TH2F*) ptyetaphibasic[4][ic]->Projection(3,0);
      rapidityeta[4][ic] = (TH2F*) ptyetaphibasic[4][ic]->Projection(2,1);
      rapidityphi[4][ic] = (TH2F*) ptyetaphibasic[4][ic]->Projection(3,1);
      etaphi[4][ic] = (TH2F*) ptyetaphibasic[4][ic]->Projection(3,2);
 
      for(int im = 2; im < 4; im++)
      {
        pt[im][ic] = (TH1F*) ptyetaphibasic[im][ic]->Projection(0);
        rapidity[im][ic] = (TH1F*) ptyetaphibasic[im][ic]->Projection(1);
        eta[im][ic] = (TH1F*) ptyetaphibasic[im][ic]->Projection(2);
        phi[im][ic] = (TH1F*) ptyetaphibasic[im][ic]->Projection(3);

        ptrapidity[im][ic] = (TH2F*) ptyetaphibasic[im][ic]->Projection(1,0);
        pteta[im][ic] = (TH2F*) ptyetaphibasic[im][ic]->Projection(2,0);
        ptphi[im][ic] = (TH2F*) ptyetaphibasic[im][ic]->Projection(3,0);
        rapidityeta[im][ic] = (TH2F*) ptyetaphibasic[im][ic]->Projection(2,1);
        rapidityphi[im][ic] = (TH2F*) ptyetaphibasic[im][ic]->Projection(3,1);
        etaphi[im][ic] = (TH2F*) ptyetaphibasic[im][ic]->Projection(3,2);

        
        for(int is = 0; is < nset; is++)
        {
          Int_t y_min_bin = ptyetaphibasic[im][ic]->GetAxis(axis[is])->FindBin(low[is]);
          Int_t y_max_bin = ptyetaphibasic[im][ic]->GetAxis(axis[is])->FindBin(high[is])-1;
          ptyetaphibasic[im][ic]->GetAxis(axis[is])->SetRange(y_min_bin,y_max_bin);

          spt[im][ic][is] = (TH1F*) ptyetaphibasic[im][ic]->Projection(0);
          srapidity[im][ic][is] = (TH1F*) ptyetaphibasic[im][ic]->Projection(1);
          seta[im][ic][is] = (TH1F*) ptyetaphibasic[im][ic]->Projection(2);
          sphi[im][ic][is] = (TH1F*) ptyetaphibasic[im][ic]->Projection(3);

          sptrapidity[im][ic][is] = (TH2F*) ptyetaphibasic[im][ic]->Projection(1,0);
          spteta[im][ic][is] = (TH2F*) ptyetaphibasic[im][ic]->Projection(2,0);
          sptphi[im][ic][is] = (TH2F*) ptyetaphibasic[im][ic]->Projection(3,0);
          srapidityeta[im][ic][is] = (TH2F*) ptyetaphibasic[im][ic]->Projection(2,1);
          srapidityphi[im][ic][is] = (TH2F*) ptyetaphibasic[im][ic]->Projection(3,1);
          setaphi[im][ic][is] = (TH2F*) ptyetaphibasic[im][ic]->Projection(3,2);

          ptyetaphibasic[im][ic]->GetAxis(axis[is])->SetRange();
          
        }

        int nbinsX = ptyetaphibasic[im][ic]->GetAxis(0)->GetNbins();
        int nbinsY = ptyetaphibasic[im][ic]->GetAxis(1)->GetNbins();
        int nbinsZ = ptyetaphibasic[im][ic]->GetAxis(2)->GetNbins();
        int nbinsT = ptyetaphibasic[im][ic]->GetAxis(3)->GetNbins();

        // Define the new histogram with the same binning as the original
        THnD* temp4D = new THnD("temp4D", "Copied 4D Histogram", 4, bins, xmin, xmax); 
        
        // Loop over the bins of the original histogram
        for (int binx = 1; binx <= nbinsX; ++binx) {
            for (int biny = 1; biny <= nbinsY; ++biny) {
                for (int binz = 1; binz <= nbinsZ; ++binz) {
                    for (int bint = 1; bint <= nbinsT; ++bint) {
        
                        // Get the values of the 3rd and 4th dimensions
                        double z_value = ptyetaphibasic[im][ic]->GetAxis(2)->GetBinCenter(binz);
                        double t_value = ptyetaphibasic[im][ic]->GetAxis(3)->GetBinCenter(bint);
        
                        // Condition to skip specific regions in the 3rd and 4th dimensions
                        if (ic == 0) {
                            if(z_value < 0)
                            {
                                if     (t_value > -0.5 && t_value < -0.2) continue;
                            }
                            else if(z_value > 0)
                            {
                                if     (t_value >  0.8 && t_value <  1.4) continue;
                                else if(t_value > -2.6 && t_value < -1.8) continue;
                            }
                        }
                        if (ic == 1) {
                            if(z_value < 0)
                            {
                                if     (t_value > -0.1 && t_value < 0.2) continue;
                            }
                            else if(z_value > 0)
                            {
                                if     (t_value > -2.2 && t_value < -1.4) continue;
                                else if(t_value >  1.3 && t_value <  1.6) continue;
                            }
                        }
        
                        // Copy the bin content and error from the original histogram to the new one
                        int idx[4] = {binx,biny,binz,bint};
                        int global_bin = ptyetaphibasic[im][ic]->GetBin(idx);
           
                        double content = ptyetaphibasic[im][ic]->GetBinContent(global_bin);
                        double error = ptyetaphibasic[im][ic]->GetBinError2(global_bin);
                        
                        // Set content and error in the new histogram
                        temp4D->SetBinContent(global_bin, content);
                        temp4D->SetBinError2(global_bin, error);
                    }
                }
            }
        }



        spt[im][ic][nset] = (TH1F*) temp4D->Projection(0);
        srapidity[im][ic][nset] = (TH1F*) temp4D->Projection(1);
        seta[im][ic][nset] = (TH1F*) temp4D->Projection(2);
        sphi[im][ic][nset] = (TH1F*) temp4D->Projection(3);

        sptrapidity[im][ic][nset] = (TH2F*) temp4D->Projection(1,0);
        spteta[im][ic][nset] = (TH2F*) temp4D->Projection(2,0);
        sptphi[im][ic][nset] = (TH2F*) temp4D->Projection(3,0);
        srapidityeta[im][ic][nset] = (TH2F*) temp4D->Projection(2,1);
        srapidityphi[im][ic][nset] = (TH2F*) temp4D->Projection(3,1);
        setaphi[im][ic][nset] = (TH2F*) temp4D->Projection(3,2);
      }

      for(int is = 0; is < nset+1; is++)
      {
        double integralData = spt[2][ic][is]->Integral(0,-1);
        double integralRC = spt[3][ic][is]->Integral(0,-1);
        double DataRC = integralData/integralRC;

        spt[3][ic][is]->Scale(DataRC);
        srapidity[3][ic][is]->Scale(DataRC);
        seta[3][ic][is]->Scale(DataRC);
        sphi[3][ic][is]->Scale(DataRC);

        sptrapidity[3][ic][is]->Scale(DataRC);
        spteta[3][ic][is]->Scale(DataRC);
        sptphi[3][ic][is]->Scale(DataRC);
        srapidityeta[3][ic][is]->Scale(DataRC);
        srapidityphi[3][ic][is]->Scale(DataRC);
        setaphi[3][ic][is]->Scale(DataRC);

        for(int im = 0; im < 2; im++)
        {
          sratiopt[im][ic][is] = (TH1F*) spt[numerator[im]][ic][is]->Clone();
          sratiopt[im][ic][is]->Divide(spt[denominator[im]][ic][is]);

          sratiorapidity[im][ic][is] = (TH1F*) srapidity[numerator[im]][ic][is]->Clone();
          sratiorapidity[im][ic][is]->Divide(srapidity[denominator[im]][ic][is]);

          sratioeta[im][ic][is] = (TH1F*) seta[numerator[im]][ic][is]->Clone();
          sratioeta[im][ic][is]->Divide(seta[denominator[im]][ic][is]);

          sratiophi[im][ic][is] = (TH1F*) sphi[numerator[im]][ic][is]->Clone();
          sratiophi[im][ic][is]->Divide(sphi[denominator[im]][ic][is]);

          sratioptrapidity[im][ic][is] = (TH2F*) sptrapidity[numerator[im]][ic][is]->Clone();
          sratioptrapidity[im][ic][is]->Divide(sptrapidity[denominator[im]][ic][is]);

          sratiopteta[im][ic][is] = (TH2F*) spteta[numerator[im]][ic][is]->Clone();
          sratiopteta[im][ic][is]->Divide(spteta[denominator[im]][ic][is]);

          sratioptphi[im][ic][is] = (TH2F*) sptphi[numerator[im]][ic][is]->Clone();
          sratioptphi[im][ic][is]->Divide(sptphi[denominator[im]][ic][is]);

          sratiorapidityeta[im][ic][is] = (TH2F*) srapidityeta[numerator[im]][ic][is]->Clone();
          sratiorapidityeta[im][ic][is]->Divide(srapidityeta[denominator[im]][ic][is]);

          sratiorapidityphi[im][ic][is] = (TH2F*) srapidityphi[numerator[im]][ic][is]->Clone();
          sratiorapidityphi[im][ic][is]->Divide(srapidityphi[denominator[im]][ic][is]);

          sratioetaphi[im][ic][is] = (TH2F*) setaphi[numerator[im]][ic][is]->Clone();
          sratioetaphi[im][ic][is]->Divide(setaphi[denominator[im]][ic][is]);
        }
      }   
      //nhitsfit[2][ic] = (TH1F*) nhitsfit[0][ic]->Clone();
      //nhitsfit[2][ic]->Add(nhitsfit[1][ic],-1.0);
      //nhitsmax[2][ic] = (TH1F*) nhitsmax[0][ic]->Clone();
      //nhitsmax[2][ic]->Add(nhitsmax[1][ic],-1.0);
      //nhitsratio[2][ic] = (TH1F*) nhitsratio[0][ic]->Clone();     
      //nhitsratio[2][ic]->Add(nhitsratio[1][ic],-1.0);
      //dedx[2][ic] = (TH1F*) dedx[0][ic]->Clone();
      //dedx[2][ic]->Add(dedx[1][ic],-1.0);
      //dca[2][ic] = (TH1F*) dca[0][ic]->Clone();
      //dca[2][ic]->Add(dca[1][ic],-1.0);

      //nhitsfit[2][ic]   = (TH1F*)ptyetaphibasic[2][ic]->Projection(4);
      //nhitsmax[2][ic]   = (TH1F*)ptyetaphibasic[2][ic]->Projection(5);
      //nhitsratio[2][ic] = (TH1F*)ptyetaphibasic[2][ic]->Projection(6);      
      //dedx[2][ic]       = (TH1F*)ptyetaphibasic[2][ic]->Projection(7); 
      //dca[2][ic]        = (TH1F*)ptyetaphibasic[2][ic]->Projection(8); 

      //nhitsfitpt[2][ic]   = (TH2F*)ptyetaphibasic[2][ic]->Projection(4,0);
      //nhitsmaxpt[2][ic]   = (TH2F*)ptyetaphibasic[2][ic]->Projection(5,0);
      //nhitsratiopt[2][ic] = (TH2F*)ptyetaphibasic[2][ic]->Projection(6,0);      
      //dedxpt[2][ic]       = (TH2F*)ptyetaphibasic[2][ic]->Projection(7,0); 
      //dcapt[2][ic]        = (TH2F*)ptyetaphibasic[2][ic]->Projection(8,0); 

      //nhitsfitrapidity[2][ic]   = (TH2F*)ptyetaphibasic[2][ic]->Projection(4,1);
      //nhitsmaxrapidity[2][ic]   = (TH2F*)ptyetaphibasic[2][ic]->Projection(5,1);
      //nhitsratiorapidity[2][ic] = (TH2F*)ptyetaphibasic[2][ic]->Projection(6,1);      
      //dedxrapidity[2][ic]       = (TH2F*)ptyetaphibasic[2][ic]->Projection(7,1); 
      //dcarapidity[2][ic]        = (TH2F*)ptyetaphibasic[2][ic]->Projection(8,1); 

      //nhitsfiteta[2][ic]   = (TH2F*)ptyetaphibasic[2][ic]->Projection(4,2);
      //nhitsmaxeta[2][ic]   = (TH2F*)ptyetaphibasic[2][ic]->Projection(5,2);
      //nhitsratioeta[2][ic] = (TH2F*)ptyetaphibasic[2][ic]->Projection(6,2);      
      //dedxeta[2][ic]       = (TH2F*)ptyetaphibasic[2][ic]->Projection(7,2); 
      //dcaeta[2][ic]        = (TH2F*)ptyetaphibasic[2][ic]->Projection(8,2); 

      //nhitsfitphi[2][ic]   = (TH2F*)ptyetaphibasic[2][ic]->Projection(4,3);
      //nhitsmaxphi[2][ic]   = (TH2F*)ptyetaphibasic[2][ic]->Projection(5,3);
      //nhitsratiophi[2][ic] = (TH2F*)ptyetaphibasic[2][ic]->Projection(6,3);      
      //dedxphi[2][ic]       = (TH2F*)ptyetaphibasic[2][ic]->Projection(7,3); 
      //dcaphi[2][ic]        = (TH2F*)ptyetaphibasic[2][ic]->Projection(8,3); 

      //Ratios
      double integralData = pt[2][ic]->Integral(0,-1);
      double integralRC = pt[3][ic]->Integral(0,-1);
      //double RCData = integralRC/integralData;
      double DataRC = integralData/integralRC;


      //Track Parameters

      //all
      //rationhitsfit[ic] = (TH1F*) nhitsfit[2][ic]->Clone();
      //rationhitsfit[ic]->Scale(RCData);
      //rationhitsfit[ic]->Divide(nhitsfit[3][ic]);

      //rationhitsmax[ic] = (TH1F*) nhitsmax[2][ic]->Clone();
      //rationhitsmax[ic]->Scale(RCData);
      //rationhitsmax[ic]->Divide(nhitsmax[3][ic]);

      //rationhitsratio[ic] = (TH1F*) nhitsratio[2][ic]->Clone();
      //rationhitsratio[ic]->Scale(RCData);
      //rationhitsratio[ic]->Divide(nhitsratio[3][ic]);

      //ratiodedx[ic] = (TH1F*) dedx[2][ic]->Clone();
      //ratiodedx[ic]->Scale(RCData);
      //ratiodedx[ic]->Divide(dedx[3][ic]);

      //ratiodca[ic] = (TH1F*) dca[2][ic]->Clone();
      //ratiodca[ic]->Scale(RCData);
      //ratiodca[ic]->Divide(dca[3][ic]);
      //all

      ////pt
      //rationhitsfitpt[ic] = (TH2F*) nhitsfitpt[2][ic]->Clone();
      //rationhitsfitpt[ic]->Scale(RCData);
      //rationhitsfitpt[ic]->Divide(nhitsfitpt[3][ic]);

      //rationhitsmaxpt[ic] = (TH2F*) nhitsmaxpt[2][ic]->Clone();
      //rationhitsmaxpt[ic]->Scale(RCData);
      //rationhitsmaxpt[ic]->Divide(nhitsmaxpt[3][ic]);

      //rationhitsratiopt[ic] = (TH2F*) nhitsratiopt[2][ic]->Clone();
      //rationhitsratiopt[ic]->Scale(RCData);
      //rationhitsratiopt[ic]->Divide(nhitsratiopt[3][ic]);

      //ratiodedxpt[ic] = (TH2F*) dedxpt[2][ic]->Clone();
      //ratiodedxpt[ic]->Scale(RCData);
      //ratiodedxpt[ic]->Divide(dedxpt[3][ic]);

      //ratiodcapt[ic] = (TH2F*) dcapt[2][ic]->Clone();
      //ratiodcapt[ic]->Scale(RCData);
      //ratiodcapt[ic]->Divide(dcapt[3][ic]);
      ////pt

      ////y
      //rationhitsfitrapidity[ic] = (TH2F*) nhitsfitrapidity[2][ic]->Clone();
      //rationhitsfitrapidity[ic]->Scale(RCData);
      //rationhitsfitrapidity[ic]->Divide(nhitsfitrapidity[3][ic]);

      //rationhitsmaxrapidity[ic] = (TH2F*) nhitsmaxrapidity[2][ic]->Clone();
      //rationhitsmaxrapidity[ic]->Scale(RCData);
      //rationhitsmaxrapidity[ic]->Divide(nhitsmaxrapidity[3][ic]);

      //rationhitsratiorapidity[ic] = (TH2F*) nhitsratiorapidity[2][ic]->Clone();
      //rationhitsratiorapidity[ic]->Scale(RCData);
      //rationhitsratiorapidity[ic]->Divide(nhitsratiorapidity[3][ic]);

      //ratiodedxrapidity[ic] = (TH2F*) dedxrapidity[2][ic]->Clone();
      //ratiodedxrapidity[ic]->Scale(RCData);
      //ratiodedxrapidity[ic]->Divide(dedxrapidity[3][ic]);

      //ratiodcarapidity[ic] = (TH2F*) dcarapidity[2][ic]->Clone();
      //ratiodcarapidity[ic]->Scale(RCData);
      //ratiodcarapidity[ic]->Divide(dcarapidity[3][ic]);
      ////y

      ////eta
      //rationhitsfiteta[ic] = (TH2F*) nhitsfiteta[2][ic]->Clone();
      //rationhitsfiteta[ic]->Scale(RCData);
      //rationhitsfiteta[ic]->Divide(nhitsfiteta[3][ic]);

      //rationhitsmaxeta[ic] = (TH2F*) nhitsmaxeta[2][ic]->Clone();
      //rationhitsmaxeta[ic]->Scale(RCData);
      //rationhitsmaxeta[ic]->Divide(nhitsmaxeta[3][ic]);

      //rationhitsratioeta[ic] = (TH2F*) nhitsratioeta[2][ic]->Clone();
      //rationhitsratioeta[ic]->Scale(RCData);
      //rationhitsratioeta[ic]->Divide(nhitsratioeta[3][ic]);

      //ratiodedxeta[ic] = (TH2F*) dedxeta[2][ic]->Clone();
      //ratiodedxeta[ic]->Scale(RCData);
      //ratiodedxeta[ic]->Divide(dedxeta[3][ic]);

      //ratiodcaeta[ic] = (TH2F*) dcaeta[2][ic]->Clone();
      //ratiodcaeta[ic]->Scale(RCData);
      //ratiodcaeta[ic]->Divide(dcaeta[3][ic]);
      ////eta
  
      ////phi
      //rationhitsfitphi[ic] = (TH2F*) nhitsfitphi[2][ic]->Clone();
      //rationhitsfitphi[ic]->Scale(RCData);
      //rationhitsfitphi[ic]->Divide(nhitsfitphi[3][ic]);

      //rationhitsmaxphi[ic] = (TH2F*) nhitsmaxphi[2][ic]->Clone();
      //rationhitsmaxphi[ic]->Scale(RCData);
      //rationhitsmaxphi[ic]->Divide(nhitsmaxphi[3][ic]);

      //rationhitsratiophi[ic] = (TH2F*) nhitsratiophi[2][ic]->Clone();
      //rationhitsratiophi[ic]->Scale(RCData);
      //rationhitsratiophi[ic]->Divide(nhitsratiophi[3][ic]);

      //ratiodedxphi[ic] = (TH2F*) dedxphi[2][ic]->Clone();
      //ratiodedxphi[ic]->Scale(RCData);
      //ratiodedxphi[ic]->Divide(dedxphi[3][ic]);

      //ratiodcaphi[ic] = (TH2F*) dcaphi[2][ic]->Clone();
      //ratiodcaphi[ic]->Scale(RCData);
      //ratiodcaphi[ic]->Divide(dcaphi[3][ic]);
      ////phi

      //Track Parameters


      pt[3][ic]->Scale(DataRC);
      rapidity[3][ic]->Scale(DataRC);
      eta[3][ic]->Scale(DataRC);
      phi[3][ic]->Scale(DataRC);

      ptrapidity[3][ic]->Scale(DataRC);
      pteta[3][ic]->Scale(DataRC);
      ptphi[3][ic]->Scale(DataRC);
      rapidityeta[3][ic]->Scale(DataRC);
      rapidityphi[3][ic]->Scale(DataRC);
      etaphi[3][ic]->Scale(DataRC);
      
      for(int im = 0; im < 2; im++)
      {
        ratiopt[im][ic] = (TH1F*) pt[numerator[im]][ic]->Clone();
        ratiopt[im][ic]->Divide(pt[denominator[im]][ic]);

        ratiorapidity[im][ic] = (TH1F*) rapidity[numerator[im]][ic]->Clone();
        ratiorapidity[im][ic]->Divide(rapidity[denominator[im]][ic]);

        ratioeta[im][ic] = (TH1F*) eta[numerator[im]][ic]->Clone();
        ratioeta[im][ic]->Divide(eta[denominator[im]][ic]);

        ratiophi[im][ic] = (TH1F*) phi[numerator[im]][ic]->Clone();
        ratiophi[im][ic]->Divide(phi[denominator[im]][ic]);

        ratioptrapidity[im][ic] = (TH2F*) ptrapidity[numerator[im]][ic]->Clone();
        ratioptrapidity[im][ic]->Divide(ptrapidity[denominator[im]][ic]);

        ratiopteta[im][ic] = (TH2F*) pteta[numerator[im]][ic]->Clone();
        ratiopteta[im][ic]->Divide(pteta[denominator[im]][ic]);

        ratioptphi[im][ic] = (TH2F*) ptphi[numerator[im]][ic]->Clone();
        ratioptphi[im][ic]->Divide(ptphi[denominator[im]][ic]);

        ratiorapidityeta[im][ic] = (TH2F*) rapidityeta[numerator[im]][ic]->Clone();
        ratiorapidityeta[im][ic]->Divide(rapidityeta[denominator[im]][ic]);

        ratiorapidityphi[im][ic] = (TH2F*) rapidityphi[numerator[im]][ic]->Clone();
        ratiorapidityphi[im][ic]->Divide(rapidityphi[denominator[im]][ic]);

        ratioetaphi[im][ic] = (TH2F*) etaphi[numerator[im]][ic]->Clone();
        ratioetaphi[im][ic]->Divide(etaphi[denominator[im]][ic]);
      }
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 2000, 800);
    c1->Divide(5,2);
    for(int i = 0; i < 10; i++)
    {
      c1->cd(i+1);
      c1->cd(i+1)->SetLeftMargin(0.12);
      c1->cd(i+1)->SetRightMargin(0.15);
      c1->cd(i+1)->SetBottomMargin(0.12);
      c1->cd(i+1)->SetTicks(1,1);
      c1->cd(i+1)->SetGrid(0,0); 
    }

    string mixinglabel[5] = {"SE","ME","SEminusME","RC","MC"};
    string mixingtitle[5] = {"Same Event","Mixed Event","Same - Mixed Event", "Reconstructed Embedding", "MC Embedding"};

    for(int ic = 0; ic < 3; ic++)
    {
      for(int im = 0; im < 5; im++)
      {
        string histtitle = Form("K%s, %s",charge[ic].c_str(),mixingtitle[im].c_str());

        c1->cd(1);
        c1->cd(1)->SetLogy(1);
        pt[im][ic]->SetTitle(histtitle.c_str()); 
        pt[im][ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
        pt[im][ic]->GetYaxis()->SetTitle("Count"); 
        pt[im][ic]->Draw("pE"); 

        c1->cd(2);
        rapidity[im][ic]->SetTitle(histtitle.c_str()); 
        rapidity[im][ic]->GetXaxis()->SetTitle("y"); 
        rapidity[im][ic]->GetYaxis()->SetTitle("Count"); 
        rapidity[im][ic]->Draw("pE"); 

        c1->cd(3);
        eta[im][ic]->SetTitle(histtitle.c_str()); 
        eta[im][ic]->GetXaxis()->SetTitle("#eta"); 
        eta[im][ic]->GetYaxis()->SetTitle("Count"); 
        eta[im][ic]->Draw("pE"); 

        c1->cd(4);
        phi[im][ic]->SetTitle(histtitle.c_str()); 
        phi[im][ic]->GetXaxis()->SetTitle("#phi"); 
        phi[im][ic]->GetYaxis()->SetTitle("Count"); 
        phi[im][ic]->Draw("pE"); 

        c1->cd(5);
        //c1->cd(5)->SetLogz(1);
        ptrapidity[im][ic]->SetTitle(histtitle.c_str()); 
        ptrapidity[im][ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
        ptrapidity[im][ic]->GetYaxis()->SetTitle("y"); 
        ptrapidity[im][ic]->SetContour(100); 
        ptrapidity[im][ic]->SetLineColor(0); 
        ptrapidity[im][ic]->Draw("Colz"); 

        c1->cd(6);
        //c1->cd(6)->SetLogz(1);
        pteta[im][ic]->SetTitle(histtitle.c_str()); 
        pteta[im][ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
        pteta[im][ic]->GetYaxis()->SetTitle("#eta"); 
        pteta[im][ic]->SetContour(100); 
        pteta[im][ic]->SetLineColor(0); 
        pteta[im][ic]->Draw("Colz"); 

        c1->cd(7);
        //c1->cd(7)->SetLogz(1);
        ptphi[im][ic]->SetTitle(histtitle.c_str()); 
        ptphi[im][ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
        ptphi[im][ic]->GetYaxis()->SetTitle("#phi"); 
        ptphi[im][ic]->SetContour(100); 
        ptphi[im][ic]->SetLineColor(0); 
        ptphi[im][ic]->Draw("Colz"); 

        c1->cd(8);
        //c1->cd(8)->SetLogz(1);
        rapidityeta[im][ic]->SetTitle(histtitle.c_str()); 
        rapidityeta[im][ic]->GetXaxis()->SetTitle("y"); 
        rapidityeta[im][ic]->GetYaxis()->SetTitle("#eta"); 
        rapidityeta[im][ic]->SetContour(100); 
        rapidityeta[im][ic]->SetLineColor(0); 
        rapidityeta[im][ic]->Draw("Colz"); 

        c1->cd(9);
        //c1->cd(9)->SetLogz(1);
        rapidityphi[im][ic]->SetTitle(histtitle.c_str()); 
        rapidityphi[im][ic]->GetXaxis()->SetTitle("y"); 
        rapidityphi[im][ic]->GetYaxis()->SetTitle("#phi"); 
        rapidityphi[im][ic]->SetContour(100); 
        rapidityphi[im][ic]->SetLineColor(0); 
        rapidityphi[im][ic]->Draw("Colz"); 

        c1->cd(10);
        //c1->cd(10)->SetLogz(1);
        etaphi[im][ic]->SetTitle(histtitle.c_str()); 
        etaphi[im][ic]->GetXaxis()->SetTitle("#eta"); 
        etaphi[im][ic]->GetYaxis()->SetTitle("#phi"); 
        etaphi[im][ic]->SetContour(100); 
        etaphi[im][ic]->SetLineColor(0); 
        etaphi[im][ic]->Draw("Colz"); 
   
        c1->SaveAs(Form("figures/KaonTTrees/%sK%s_%s.pdf",folderopt.c_str(),charge[ic].c_str(),mixinglabel[im].c_str()));
      }
    }

    string ratiotitle[2] = {"Data/RC","RC/Data"};
    string ratiofile[2]  = {"DataRCRatio","RCDataRatio"};

    TCanvas *cr = new TCanvas("cr", "cr", 1600, 1600);
    cr->Divide(4,4);
    for(int i = 0; i < 16; i++)
    {
      cr->cd(i+1);
      cr->cd(i+1)->SetLeftMargin(0.12);
      cr->cd(i+1)->SetRightMargin(0.15);
      cr->cd(i+1)->SetBottomMargin(0.12);
      cr->cd(i+1)->SetTicks(1,1);
      cr->cd(i+1)->SetGrid(0,0); 
    }

    for(int im = 0; im < 2; im++)
    {
      for(int ic = 0; ic < 2; ic++)
      {
        string histtitle = Form("K%s %s",charge[ic].c_str(),ratiotitle[im].c_str());

        cr->cd(1);
        cr->cd(1)->SetLogy(0);
        ratiopt[im][ic]->SetTitle(histtitle.c_str()); 
        ratiopt[im][ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
        ratiopt[im][ic]->GetYaxis()->SetTitle("Data/RC"); 
        ratiopt[im][ic]->Draw("pE"); 

        cr->cd(6);
        ratiorapidity[im][ic]->SetTitle(histtitle.c_str()); 
        ratiorapidity[im][ic]->GetXaxis()->SetTitle("y"); 
        ratiorapidity[im][ic]->GetYaxis()->SetTitle("Data/RC"); 
        ratiorapidity[im][ic]->Draw("pE"); 

        cr->cd(11);
        ratioeta[im][ic]->SetTitle(histtitle.c_str()); 
        ratioeta[im][ic]->GetXaxis()->SetTitle("#eta"); 
        ratioeta[im][ic]->GetYaxis()->SetTitle("Data/RC"); 
        ratioeta[im][ic]->Draw("pE"); 

        cr->cd(16);
        ratiophi[im][ic]->SetTitle(histtitle.c_str()); 
        ratiophi[im][ic]->GetXaxis()->SetTitle("#phi"); 
        ratiophi[im][ic]->GetYaxis()->SetTitle("Data/RC"); 
        ratiophi[im][ic]->Draw("pE"); 

        cr->cd(5);
        cr->cd(5)->SetLogz(1);
        ratioptrapidity[im][ic]->SetMinimum(0.001); 
        ratioptrapidity[im][ic]->SetMaximum(500); 
        ratioptrapidity[im][ic]->SetTitle(histtitle.c_str()); 
        ratioptrapidity[im][ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
        ratioptrapidity[im][ic]->GetYaxis()->SetTitle("y"); 
        ratioptrapidity[im][ic]->SetContour(100); 
        ratioptrapidity[im][ic]->Draw("Colz"); 

        cr->cd(9);
        cr->cd(9)->SetLogz(1);
        ratiopteta[im][ic]->SetMinimum(0.001); 
        ratiopteta[im][ic]->SetMaximum(500); 
        ratiopteta[im][ic]->SetTitle(histtitle.c_str()); 
        ratiopteta[im][ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
        ratiopteta[im][ic]->GetYaxis()->SetTitle("#eta"); 
        ratiopteta[im][ic]->SetContour(100); 
        ratiopteta[im][ic]->Draw("Colz"); 

        cr->cd(13);
        cr->cd(13)->SetLogz(1);
        ratioptphi[im][ic]->SetMinimum(0.001); 
        ratioptphi[im][ic]->SetMaximum(500); 
        ratioptphi[im][ic]->SetTitle(histtitle.c_str()); 
        ratioptphi[im][ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
        ratioptphi[im][ic]->GetYaxis()->SetTitle("#phi"); 
        ratioptphi[im][ic]->SetContour(100); 
        ratioptphi[im][ic]->Draw("Colz"); 

        cr->cd(10);
        cr->cd(10)->SetLogz(1);
        ratiorapidityeta[im][ic]->SetMinimum(0.001); 
        ratiorapidityeta[im][ic]->SetMaximum(500); 
        ratiorapidityeta[im][ic]->SetTitle(histtitle.c_str()); 
        ratiorapidityeta[im][ic]->GetXaxis()->SetTitle("y"); 
        ratiorapidityeta[im][ic]->GetYaxis()->SetTitle("#eta"); 
        ratiorapidityeta[im][ic]->SetContour(100); 
        ratiorapidityeta[im][ic]->Draw("Colz"); 

        cr->cd(14);
        cr->cd(14)->SetLogz(1);
        ratiorapidityphi[im][ic]->SetMinimum(0.001); 
        ratiorapidityphi[im][ic]->SetMaximum(500); 
        ratiorapidityphi[im][ic]->SetTitle(histtitle.c_str()); 
        ratiorapidityphi[im][ic]->GetXaxis()->SetTitle("y"); 
        ratiorapidityphi[im][ic]->GetYaxis()->SetTitle("#phi"); 
        ratiorapidityphi[im][ic]->SetContour(100); 
        ratiorapidityphi[im][ic]->Draw("Colz"); 

        cr->cd(15);
        cr->cd(15)->SetLogz(1);
        ratioetaphi[im][ic]->SetMinimum(0.001); 
        ratioetaphi[im][ic]->SetMaximum(500); 
        ratioetaphi[im][ic]->SetTitle(histtitle.c_str()); 
        ratioetaphi[im][ic]->GetXaxis()->SetTitle("#eta"); 
        ratioetaphi[im][ic]->GetYaxis()->SetTitle("#phi"); 
        ratioetaphi[im][ic]->SetContour(100); 
        ratioetaphi[im][ic]->Draw("Colz"); 
   
        cr->SaveAs(Form("figures/KaonTTrees/%sK%s_%s.pdf",folderopt.c_str(),charge[ic].c_str(),ratiofile[im].c_str()));
        
      }
    }

    
    string settitle[nset+1] = {"-1.5<y<-0.5","-0.5<y<0.5","0.5<y<1.5","-1.0<y<-0.9","0.9<y<1.0","-1.0<y<-0.8","0.8<y<1.0","1.0<p_{T}<2.0 GeV/c","0.0<p_{T}<0.7 GeV/c","p_{T}>2.0 GeV/c","0.44<p_{T}<1.5 GeV/c","0.3<p_{T}<0.5 GeV/c","0.5<p_{T}<0.7 GeV/c","0.7<p_{T}<0.9 GeV/c","0.9<p_{T}<1.1 GeV/c","1.1<p_{T}<1.3 GeV/c","1.3<p_{T}<1.5 GeV/c","1.5<p_{T}<1.7 GeV/c","Ignore Major #eta,#phi Deficiencies"};
    string setfile[nset+1]  = {"n1p5yn0p5"  ,"n0p5y0p5"  ,"0p5y1p5","n1yn0p9","0p9y1","n1yn0p8","0p8y1"  ,"1p0pt2p0",           "0pt0p7",             "pt2p0"          ,"0p44pt1p5"           ,"0p3pt0p5"           ,"0p5pt0p7"           ,"0p7pt0p9"           ,"0p9pt1p1"           ,"1p1pt1p3"           ,"1p3pt1p5"           ,"1p5pt1p7",       "IgnoreEtaPhiDef" };
    for(int is = 0; is < nset+1; is++)
    {
      for(int im = 0; im < 2; im++)
      {
        for(int ic = 0; ic < 2; ic++)
        {
          string histtitle = Form("K%s %s %s",charge[ic].c_str(),ratiotitle[im].c_str(),settitle[is].c_str());

          cr->cd(1);
          cr->cd(1)->SetLogy(0);
          sratiopt[im][ic][is]->SetTitle(histtitle.c_str()); 
          sratiopt[im][ic][is]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
          sratiopt[im][ic][is]->GetYaxis()->SetTitle("Data/RC"); 
          sratiopt[im][ic][is]->Draw("pE"); 

          cr->cd(6);
          sratiorapidity[im][ic][is]->SetTitle(histtitle.c_str()); 
          sratiorapidity[im][ic][is]->GetXaxis()->SetTitle("y"); 
          sratiorapidity[im][ic][is]->GetYaxis()->SetTitle("Data/RC"); 
          sratiorapidity[im][ic][is]->Draw("pE"); 

          cr->cd(11);
          sratioeta[im][ic][is]->SetTitle(histtitle.c_str()); 
          sratioeta[im][ic][is]->GetXaxis()->SetTitle("#eta"); 
          sratioeta[im][ic][is]->GetYaxis()->SetTitle("Data/RC"); 
          sratioeta[im][ic][is]->Draw("pE"); 

          cr->cd(16);
          sratiophi[im][ic][is]->SetTitle(histtitle.c_str()); 
          sratiophi[im][ic][is]->GetXaxis()->SetTitle("#phi"); 
          sratiophi[im][ic][is]->GetYaxis()->SetTitle("Data/RC"); 
          sratiophi[im][ic][is]->Draw("pE"); 

          cr->cd(5);
          //cr->cd(5)->SetLogz(1);
          sratioptrapidity[im][ic][is]->SetTitle(histtitle.c_str()); 
          sratioptrapidity[im][ic][is]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
          sratioptrapidity[im][ic][is]->GetYaxis()->SetTitle("y"); 
          sratioptrapidity[im][ic][is]->SetContour(100); 
          sratioptrapidity[im][ic][is]->Draw("Colz"); 

          cr->cd(9);
          //cr->cd(9)->SetLogz(1);
          sratiopteta[im][ic][is]->SetTitle(histtitle.c_str()); 
          sratiopteta[im][ic][is]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
          sratiopteta[im][ic][is]->GetYaxis()->SetTitle("#eta"); 
          sratiopteta[im][ic][is]->SetContour(100); 
          sratiopteta[im][ic][is]->Draw("Colz"); 

          cr->cd(13);
          //cr->cd(13)->SetLogz(1);
          sratioptphi[im][ic][is]->SetTitle(histtitle.c_str()); 
          sratioptphi[im][ic][is]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
          sratioptphi[im][ic][is]->GetYaxis()->SetTitle("#phi"); 
          sratioptphi[im][ic][is]->SetContour(100); 
          sratioptphi[im][ic][is]->Draw("Colz"); 

          cr->cd(10);
          //cr->cd(10)->SetLogz(1);
          sratiorapidityeta[im][ic][is]->SetTitle(histtitle.c_str()); 
          sratiorapidityeta[im][ic][is]->GetXaxis()->SetTitle("y"); 
          sratiorapidityeta[im][ic][is]->GetYaxis()->SetTitle("#eta"); 
          sratiorapidityeta[im][ic][is]->SetContour(100); 
          sratiorapidityeta[im][ic][is]->Draw("Colz"); 

          cr->cd(14);
          //cr->cd(14)->SetLogz(1);
          sratiorapidityphi[im][ic][is]->SetTitle(histtitle.c_str()); 
          sratiorapidityphi[im][ic][is]->GetXaxis()->SetTitle("y"); 
          sratiorapidityphi[im][ic][is]->GetYaxis()->SetTitle("#phi"); 
          sratiorapidityphi[im][ic][is]->SetContour(100); 
          sratiorapidityphi[im][ic][is]->Draw("Colz"); 

          cr->cd(15);
          //cr->cd(15)->SetLogz(1);
          sratioetaphi[im][ic][is]->SetTitle(histtitle.c_str()); 
          sratioetaphi[im][ic][is]->GetXaxis()->SetTitle("#eta"); 
          sratioetaphi[im][ic][is]->GetYaxis()->SetTitle("#phi"); 
          sratioetaphi[im][ic][is]->SetContour(100); 
          sratioetaphi[im][ic][is]->Draw("Colz"); 
   
          cr->SaveAs(Form("figures/KaonTTrees/%sK%s_%s_%s.pdf",folderopt.c_str(),charge[ic].c_str(),ratiofile[im].c_str(),setfile[is].c_str()));
          
        }
      }
    }

    //for(int im = 0; im < 4; im++)
    //{
    //  for(int ic = 0; ic < 2; ic++)
    //  {
    //    string histtitle = Form("K%s, %s",charge[ic].c_str(),mixingtitle[im].c_str());

    //    c1->cd(1+5*ic);
    //    c1->cd(1+5*ic)->SetLogy(0);
    //    nhitsfit[im][ic]->SetTitle(histtitle.c_str()); 
    //    nhitsfit[im][ic]->GetXaxis()->SetTitle("NHitsFit"); 
    //    nhitsfit[im][ic]->GetYaxis()->SetTitle("Count"); 
    //    nhitsfit[im][ic]->Draw("pE"); 

    //    c1->cd(2+5*ic);
    //    c1->cd(2+5*ic)->SetLogy(0);
    //    nhitsmax[im][ic]->SetTitle(histtitle.c_str()); 
    //    nhitsmax[im][ic]->GetXaxis()->SetTitle("NHitsMax"); 
    //    nhitsmax[im][ic]->GetYaxis()->SetTitle("Count"); 
    //    nhitsmax[im][ic]->Draw("pE"); 
    //    
    //    c1->cd(3+5*ic);
    //    c1->cd(3+5*ic)->SetLogy(0);
    //    nhitsratio[im][ic]->SetTitle(histtitle.c_str()); 
    //    nhitsratio[im][ic]->GetXaxis()->SetTitle("NHitsFit/NHitsMax"); 
    //    nhitsratio[im][ic]->GetYaxis()->SetTitle("Count"); 
    //    nhitsratio[im][ic]->Draw("pE"); 

    //    c1->cd(4+5*ic);
    //    c1->cd(4+5*ic)->SetLogy(1);
    //    dedx[im][ic]->SetTitle(histtitle.c_str()); 
    //    dedx[im][ic]->GetXaxis()->SetTitle("dE/dx"); 
    //    dedx[im][ic]->GetYaxis()->SetTitle("Count"); 
    //    dedx[im][ic]->Draw("pE"); 

    //    c1->cd(5+5*ic);
    //    c1->cd(5+5*ic)->SetLogy(1);
    //    dca[im][ic]->SetTitle(histtitle.c_str()); 
    //    dca[im][ic]->GetXaxis()->SetTitle("|dca| (cm)"); 
    //    dca[im][ic]->GetYaxis()->SetTitle("Count"); 
    //    dca[im][ic]->Draw("pE"); 
    //  }
    //  c1->SaveAs(Form("figures/KaonTTrees/KaonBasicVariables_%s.pdf",mixinglabel[im].c_str()));
    //}

    //for(int ic = 0; ic < 2; ic++)
    //{
    //  string histtitle = Form("K%s, Data/RC",charge[ic].c_str());

    //  c1->cd(1+5*ic);
    //  c1->cd(1+5*ic)->SetLogy(0);
    //  rationhitsfit[ic]->SetTitle(histtitle.c_str()); 
    //  rationhitsfit[ic]->GetXaxis()->SetTitle("NHitsFit"); 
    //  rationhitsfit[ic]->GetYaxis()->SetTitle("Data/RC"); 
    //  rationhitsfit[ic]->Draw("pE"); 

    //  c1->cd(2+5*ic);
    //  c1->cd(2+5*ic)->SetLogy(0);
    //  rationhitsmax[ic]->SetTitle(histtitle.c_str()); 
    //  rationhitsmax[ic]->GetXaxis()->SetTitle("NHitsMax"); 
    //  rationhitsmax[ic]->GetYaxis()->SetTitle("Data/RC"); 
    //  rationhitsmax[ic]->Draw("pE"); 
    //  
    //  c1->cd(3+5*ic);
    //  c1->cd(3+5*ic)->SetLogy(0);
    //  rationhitsratio[ic]->SetTitle(histtitle.c_str()); 
    //  rationhitsratio[ic]->GetXaxis()->SetTitle("NHitsFit/NHitsMax"); 
    //  rationhitsratio[ic]->GetYaxis()->SetTitle("Data/RC"); 
    //  rationhitsratio[ic]->Draw("pE"); 

    //  c1->cd(4+5*ic);
    //  c1->cd(4+5*ic)->SetLogy(0);
    //  ratiodedx[ic]->SetTitle(histtitle.c_str()); 
    //  ratiodedx[ic]->GetXaxis()->SetTitle("dE/dx"); 
    //  ratiodedx[ic]->GetYaxis()->SetTitle("Data/RC"); 
    //  ratiodedx[ic]->Draw("pE"); 

    //  c1->cd(5+5*ic);
    //  c1->cd(5+5*ic)->SetLogy(0);
    //  ratiodca[ic]->SetTitle(histtitle.c_str()); 
    //  ratiodca[ic]->GetXaxis()->SetTitle("|dca| (cm)"); 
    //  ratiodca[ic]->GetYaxis()->SetTitle("Data/RC"); 
    //  ratiodca[ic]->Draw("pE"); 
    //}
    //c1->SaveAs(Form("figures/KaonTTrees/KaonBasicVariables_DataRCRatio.pdf"));
    

    
    

    //// Create histograms for flat and desired pT-phi distributions
    //const int nBinsPt = 50;
    //const int nBinsY  = 30;
    //const int nBinsPhi = 50;
    //double ptMin = 0.0, ptMax = 5.0;
    //double yMin = -1.5, yMax = 1.5;
    //double phiMin = 0, phiMax = 2 * TMath::Pi();
    //
    //double resolution[9]    = {0.221173,    0.296191,    0.411608,    0.535956,    0.630288,    0.675484,    0.648544,    0.539959,    0.377757  };  
    //double resolutionerr[9] = {0.000461111, 0.000319516, 0.000202537, 0.000141295, 0.000102773, 8.85656e-05, 9.71308e-05, 0.000193934, 0.00031656};


    //string spectra = "PhiEmbedding_noweights_nomccuts_rapidity_allcent";
    //string filename = "PhiEmbedding_noweights_nomccuts_rapidity_allcent.root";

    //string inputfile = Form("effaccfiles/%s/%s/%s/%s",vmsa::mPID[0].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),filename.c_str());
    //TFile *File_Input = TFile::Open(inputfile.c_str());

    //TH3F *h_flat_pt_phi[10];
    //TH3F *h_desired_pt_phi[10];
    //TH1F *h_flat_pt[10];
    //TH1F *h_flat_y[10];
    //TH1F *h_flat_phi[10];
    //TH1F *h_desired_pt[10];
    //TH1F *h_desired_y[10];
    //TH1F *h_desired_phi[10];
    //for(int icent = 0; icent < 10; icent++)
    //{
    //  string flat = Form("h_flat_pt_phi_cent%d",icent);
    //  h_flat_pt_phi[icent] = (TH3F*) ((TH3F*)File_Input->Get(Form("h_mRc0EffPtYPhiPsi_Cent_%d",icent)))->Clone(flat.c_str());
    //  cout << "NBinsX = " << h_flat_pt_phi[icent]->GetNbinsX() << endl;
    //  cout << "NBinsY = " << h_flat_pt_phi[icent]->GetNbinsY() << endl;
    //  cout << "NBinsZ = " << h_flat_pt_phi[icent]->GetNbinsZ() << endl;
    //}

    //// Create a histogram for the ratio
    //TH2F *h_flat_pty[10];   
    //TH2F *h_flat_ptphi[10];   
    //TH2F *h_flat_yphi[10];   
    //for(int icent = 0; icent < 10; icent++)
    //{
    //  string flat = Form("h_flat_pt_cent%d",icent);
    //  h_flat_pt[icent] = (TH1F*) h_flat_pt_phi[icent]->ProjectionX(flat.c_str(),0,-1,0,-1,"e");    
    //  flat = Form("h_flat_y_cent%d",icent);
    //  h_flat_y[icent] = (TH1F*) h_flat_pt_phi[icent]->ProjectionY(flat.c_str(),0,-1,0,-1,"e");    
    //  flat = Form("h_flat_phi_cent%d",icent);
    //  h_flat_phi[icent] = (TH1F*) h_flat_pt_phi[icent]->ProjectionZ(flat.c_str(),0,-1,0,-1,"e");    
    //
    //  h_flat_pty[icent] = (TH2F*) h_flat_pt_phi[icent]->Project3D("yx");
    //  h_flat_ptphi[icent] = (TH2F*) h_flat_pt_phi[icent]->Project3D("zx");
    //  h_flat_yphi[icent] = (TH2F*) h_flat_pt_phi[icent]->Project3D("zy");
    //}

    //// Optionally, draw the histograms
    //TCanvas *c1 = new TCanvas("c1", "pT-phi Weights", 1200, 800);
    //c1->Divide(3,2);
    //for(int i = 0; i < 6; i++)
    //{
    //  c1->cd(i+1);
    //  c1->cd(i+1)->SetLeftMargin(0.15);
    //  c1->cd(i+1)->SetRightMargin(0.15);
    //  c1->cd(i+1)->SetBottomMargin(0.15);
    //  c1->cd(i+1)->SetTicks(1,1);
    //  c1->cd(i+1)->Marker    //}
    // 
    //string outputname = Form("figures/EmbeddingWeights/allcent_pt_y_phi.pdf");
    //string outputstart = Form("%s[",outputname.c_str());
    //string outputstop = Form("%s]",outputname.c_str());

    //c1->Print(outputstart.c_str());
    //for(int icent = 0; icent < 10; icent++)
    //{
    //  c1->cd(1);
    //  h_flat_pt[icent]->SetStats(0);
    //  h_flat_pt[icent]->GetXaxis()->SetTitle("p_{T} GeV/c");
    //  h_flat_pt[icent]->GetYaxis()->SetTitle("Counts");
    //  h_flat_pt[icent]->Draw("pE");

    //  c1->cd(2);
    //  h_flat_y[icent]->SetStats(0);
    //  h_flat_y[icent]->GetXaxis()->SetTitle("y");
    //  h_flat_y[icent]->GetYaxis()->SetTitle("Counts");
    //  h_flat_y[icent]->Draw("pE");

    //  c1->cd(3);
    //  h_flat_phi[icent]->SetStats(0);
    //  h_flat_phi[icent]->GetXaxis()->SetTitle("#phi-#Psi_{2}");
    //  h_flat_phi[icent]->GetYaxis()->SetTitle("Counts");
    //  h_flat_phi[icent]->Draw("pE");

    //  c1->cd(4);
    //  h_flat_pty[icent]->SetStats(0);
    //  h_flat_pty[icent]->GetXaxis()->SetTitle("p_{T} GeV/c");
    //  h_flat_pty[icent]->GetYaxis()->SetTitle("y");
    //  h_flat_pty[icent]->Draw("colz");

    //  c1->cd(5);
    //  h_flat_ptphi[icent]->SetStats(0);
    //  h_flat_ptphi[icent]->GetXaxis()->SetTitle("p_{T} GeV/c");
    //  h_flat_ptphi[icent]->GetYaxis()->SetTitle("#phi-#Psi_{2}");
    //  h_flat_ptphi[icent]->Draw("colz");

    //  c1->cd(6);
    //  h_flat_yphi[icent]->SetStats(0);
    //  h_flat_yphi[icent]->GetXaxis()->SetTitle("y");
    //  h_flat_yphi[icent]->GetYaxis()->SetTitle("#phi-#Psi_{2}");
    //  h_flat_yphi[icent]->Draw("colz");


    //  c1->Update();
    //  c1->Print(outputname.c_str());
    //}
    //c1->Print(outputstop.c_str());   

    TCanvas *c43 = new TCanvas("c43", "c43", 1600, 1200);
    c43->Divide(4,3);
    for(int i = 0; i < 12; i++)
    {
      c43->cd(i+1);
      c43->cd(i+1)->SetLeftMargin(0.12);
      c43->cd(i+1)->SetRightMargin(0.15);
      c43->cd(i+1)->SetBottomMargin(0.12);
      c43->cd(i+1)->SetTicks(1,1);
      c43->cd(i+1)->SetGrid(0,0); 
    }


    for(int ic = 0; ic < 2; ic++)
    {
      c43->cd(1);
      c43->cd(1)->SetLogy();
      pt[2][ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
      pt[2][ic]->GetYaxis()->SetTitle("Count");
      pt[2][ic]->SetMarkerStyle(20);
      pt[2][ic]->SetMarkerColor(kBlue);
      pt[2][ic]->SetMarkerColor(kBlue);
      pt[2][ic]->Draw("pE"); 
      pt[3][ic]->SetMarkerStyle(20);
      pt[3][ic]->SetMarkerColor(kOrange+7);
      pt[3][ic]->SetMarkerColor(kOrange+7);
      pt[3][ic]->Draw("pE same"); 
   
      TLegend *legpt = new TLegend(0.3,0.2,0.45,0.35);
      legpt->AddEntry(pt[2][ic],"Data","p");
      legpt->AddEntry(pt[3][ic],"RC","p");
      legpt->Draw("same");

      c43->cd(2);
      rapidity[2][ic]->GetXaxis()->SetTitle("y"); 
      rapidity[2][ic]->GetYaxis()->SetTitle("Count"); 
      rapidity[2][ic]->SetMarkerStyle(20);
      rapidity[2][ic]->SetMarkerColor(kBlue);
      rapidity[2][ic]->SetMarkerColor(kBlue);
      rapidity[2][ic]->Draw("pE"); 
      rapidity[3][ic]->SetMarkerStyle(20);
      rapidity[3][ic]->SetMarkerColor(kOrange+7);
      rapidity[3][ic]->SetMarkerColor(kOrange+7);
      rapidity[3][ic]->Draw("pE same"); 
   
      TLegend *legy = new TLegend(0.2,0.2,0.6,0.6);
      legy->AddEntry(rapidity[2][ic],"Data","p");
      legy->AddEntry(rapidity[3][ic],"RC","p");
      //legy->Draw("same");

      c43->cd(3);
      eta[2][ic]->GetXaxis()->SetTitle("#eta"); 
      eta[2][ic]->GetYaxis()->SetTitle("Count"); 
      eta[2][ic]->SetMarkerStyle(20);
      eta[2][ic]->SetMarkerColor(kBlue);
      eta[2][ic]->SetMarkerColor(kBlue);
      eta[2][ic]->Draw("pE"); 
      eta[3][ic]->SetMarkerStyle(20);
      eta[3][ic]->SetMarkerColor(kOrange+7);
      eta[3][ic]->SetMarkerColor(kOrange+7);
      eta[3][ic]->Draw("pE same"); 

      TLegend *legeta = new TLegend(0.2,0.2,0.6,0.6);
      legeta->AddEntry(eta[2][ic],"Data","p");
      legeta->AddEntry(eta[3][ic],"RC","p");
      //legeta->Draw("same");

      c43->cd(4);
      phi[2][ic]->GetXaxis()->SetTitle("#phi"); 
      phi[2][ic]->GetYaxis()->SetTitle("Count"); 
      phi[2][ic]->Draw("pE"); 
      phi[2][ic]->SetMarkerStyle(20);
      phi[2][ic]->SetMarkerColor(kBlue);
      phi[2][ic]->SetMarkerColor(kBlue);
      phi[2][ic]->Draw("pE"); 
      phi[3][ic]->SetMarkerStyle(20);
      phi[3][ic]->SetMarkerColor(kOrange+7);
      phi[3][ic]->SetMarkerColor(kOrange+7);
      phi[3][ic]->Draw("pE same"); 

      TLegend *legphi = new TLegend(0.2,0.2,0.6,0.6);
      legphi->AddEntry(phi[2][ic],"Data","p");
      legphi->AddEntry(phi[3][ic],"RC","p");
      //legphi->Draw("same");

      c43->cd(5);
      ratiopt[0][ic]->Draw("pE"); 
      ratiopt[0][ic]->GetYaxis()->SetRangeUser(0.25,1.75); 

      c43->cd(6);
      ratiorapidity[0][ic]->Draw("pE"); 
      ratiorapidity[0][ic]->GetYaxis()->SetRangeUser(0.25,1.75); 

      c43->cd(7);
      ratioeta[0][ic]->Draw("pE"); 
      ratioeta[0][ic]->GetYaxis()->SetRangeUser(0.25,1.75); 

      c43->cd(8);
      ratiophi[0][ic]->Draw("pE"); 
      ratiophi[0][ic]->GetYaxis()->SetRangeUser(0.25,1.75); 

      c43->cd(9);
      ratiopt[1][ic]->Draw("pE"); 
      ratiopt[1][ic]->GetYaxis()->SetRangeUser(0.25,1.75); 

      c43->cd(10);
      ratiorapidity[1][ic]->Draw("pE"); 
      ratiorapidity[1][ic]->GetYaxis()->SetRangeUser(0.25,1.75); 

      c43->cd(11);
      ratioeta[1][ic]->Draw("pE"); 
      ratioeta[1][ic]->GetYaxis()->SetRangeUser(0.25,1.75); 

      c43->cd(12);
      ratiophi[1][ic]->Draw("pE"); 
      ratiophi[1][ic]->GetYaxis()->SetRangeUser(0.25,1.75); 

      c43->SaveAs(Form("figures/KaonTTrees/%sK%s_1DRatios.pdf",folderopt.c_str(),charge[ic].c_str()));
    }

    for(int is = 0; is < nset+1; is++)
    {
      for(int ic = 0; ic < 2; ic++)
      {
        c43->cd(1);
        c43->cd(1)->SetLogy();
        spt[2][ic][is]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
        spt[2][ic][is]->GetYaxis()->SetTitle("Count");
        spt[2][ic][is]->SetMarkerStyle(20);
        spt[2][ic][is]->SetMarkerColor(kBlue);
        spt[2][ic][is]->SetMarkerColor(kBlue);
        spt[2][ic][is]->Draw("pE"); 
        spt[3][ic][is]->SetMarkerStyle(20);
        spt[3][ic][is]->SetMarkerColor(kOrange+7);
        spt[3][ic][is]->SetMarkerColor(kOrange+7);
        spt[3][ic][is]->Draw("pE same"); 
   
        TLegend *legpt = new TLegend(0.3,0.2,0.45,0.35);
        legpt->AddEntry(spt[2][ic][is],"Data","p");
        legpt->AddEntry(spt[3][ic][is],"RC","p");
        legpt->Draw("same");

        c43->cd(2);
        srapidity[2][ic][is]->GetXaxis()->SetTitle("y"); 
        srapidity[2][ic][is]->GetYaxis()->SetTitle("Count"); 
        srapidity[2][ic][is]->SetMarkerStyle(20);
        srapidity[2][ic][is]->SetMarkerColor(kBlue);
        srapidity[2][ic][is]->SetMarkerColor(kBlue);
        srapidity[2][ic][is]->Draw("pE"); 
        srapidity[3][ic][is]->SetMarkerStyle(20);
        srapidity[3][ic][is]->SetMarkerColor(kOrange+7);
        srapidity[3][ic][is]->SetMarkerColor(kOrange+7);
        srapidity[3][ic][is]->Draw("pE same"); 
   
        TLegend *legy = new TLegend(0.2,0.2,0.6,0.6);
        legy->AddEntry(srapidity[2][ic][is],"Data","p");
        legy->AddEntry(srapidity[3][ic][is],"RC","p");
        //legy->Draw("same");

        c43->cd(3);
        seta[2][ic][is]->GetXaxis()->SetTitle("#eta"); 
        seta[2][ic][is]->GetYaxis()->SetTitle("Count"); 
        seta[2][ic][is]->SetMarkerStyle(20);
        seta[2][ic][is]->SetMarkerColor(kBlue);
        seta[2][ic][is]->SetMarkerColor(kBlue);
        seta[2][ic][is]->Draw("pE"); 
        seta[3][ic][is]->SetMarkerStyle(20);
        seta[3][ic][is]->SetMarkerColor(kOrange+7);
        seta[3][ic][is]->SetMarkerColor(kOrange+7);
        seta[3][ic][is]->Draw("pE same"); 

        TLegend *legeta = new TLegend(0.2,0.2,0.6,0.6);
        legeta->AddEntry(seta[2][ic][is],"Data","p");
        legeta->AddEntry(seta[3][ic][is],"RC","p");
        //legeta->Draw("same");

        c43->cd(4);
        sphi[2][ic][is]->GetXaxis()->SetTitle("#phi"); 
        sphi[2][ic][is]->GetYaxis()->SetTitle("Count"); 
        sphi[2][ic][is]->Draw("pE"); 
        sphi[2][ic][is]->SetMarkerStyle(20);
        sphi[2][ic][is]->SetMarkerColor(kBlue);
        sphi[2][ic][is]->SetMarkerColor(kBlue);
        sphi[2][ic][is]->Draw("pE"); 
        sphi[3][ic][is]->SetMarkerStyle(20);
        sphi[3][ic][is]->SetMarkerColor(kOrange+7);
        sphi[3][ic][is]->SetMarkerColor(kOrange+7);
        sphi[3][ic][is]->Draw("pE same"); 

        TLegend *legphi = new TLegend(0.2,0.2,0.6,0.6);
        legphi->AddEntry(sphi[2][ic][is],"Data","p");
        legphi->AddEntry(sphi[3][ic][is],"RC","p");
        //legphi->Draw("same");

        c43->cd(5);
        sratiopt[0][ic][is]->GetYaxis()->SetRangeUser(0.25,1.75); 
        sratiopt[0][ic][is]->Draw("pE"); 

        c43->cd(6);
        sratiorapidity[0][ic][is]->GetYaxis()->SetRangeUser(0.25,1.75); 
        sratiorapidity[0][ic][is]->Draw("pE"); 

        c43->cd(7);
        sratioeta[0][ic][is]->GetYaxis()->SetRangeUser(0.25,1.75); 
        sratioeta[0][ic][is]->Draw("pE"); 

        c43->cd(8);
        sratiophi[0][ic][is]->GetYaxis()->SetRangeUser(0.25,1.75); 
        sratiophi[0][ic][is]->Draw("pE"); 

        c43->cd(9);
        sratiopt[1][ic][is]->GetYaxis()->SetRangeUser(0.25,1.75); 
        sratiopt[1][ic][is]->Draw("pE"); 

        c43->cd(10);
        sratiorapidity[1][ic][is]->GetYaxis()->SetRangeUser(0.25,1.75); 
        sratiorapidity[1][ic][is]->Draw("pE"); 

        c43->cd(11);
        sratioeta[1][ic][is]->GetYaxis()->SetRangeUser(0.25,1.75); 
        sratioeta[1][ic][is]->Draw("pE"); 

        c43->cd(12);
        sratiophi[1][ic][is]->GetYaxis()->SetRangeUser(0.25,1.75); 
        sratiophi[1][ic][is]->Draw("pE"); 

        c43->SaveAs(Form("figures/KaonTTrees/%sK%s_1DRatios_%s.pdf",folderopt.c_str(),charge[ic].c_str(),setfile[is].c_str()));
      }
    }

    TH1F *kaoncos[2][5];
    TH1F *kaoncosratio[2][5];
    for(int is = 0; is < 5; is++)
    {
      for(int ic = 0; ic < 2; ic++)
      {
        string hist;
        hist = Form("kaoncos_%d_%d",ic,is);
        kaoncos[ic][is] = (TH1F*)kaonptcos[ic][is]->ProjectionY(hist.c_str(),0,-1);
      }
    }
    for(int is = 1; is < 5; is++)
    {
      for(int ic = 0; ic < 2; ic++)
      {
        string hist;
        hist = Form("kaonptcosratio_%d_%d",ic,is);
        kaonptcosratio[ic][is] = (TH2F*)kaonptcos[ic][is]->Clone(hist.c_str());
        kaonptcosratio[ic][is]->Divide(kaonptcos[ic][is],kaonptcos[ic][0],1,1,"B");

        hist = Form("kaoncosratio_%d_%d",ic,is);
        kaoncosratio[ic][is] = (TH1F*)kaoncos[ic][is]->Clone(hist.c_str());
        kaoncosratio[ic][is]->Divide(kaoncos[ic][is],kaoncos[ic][0],1,1,"B");
      }
    }
   
    TCanvas *ckaon = new TCanvas("ckaon", "ckaon", 1600, 800);
    ckaon->Divide(4,2);
    for(int i = 0; i < 8; i++)
    {
      ckaon->cd(i+1);
      ckaon->cd(i+1)->SetLeftMargin(0.12);
      ckaon->cd(i+1)->SetRightMargin(0.15);
      ckaon->cd(i+1)->SetBottomMargin(0.12);
      ckaon->cd(i+1)->SetTicks(1,1);
      ckaon->cd(i+1)->SetGrid(0,0); 
    }

    for(int ic = 0; ic < 2; ic++) 
    { 
      ckaon->cd(1+4*ic);
      kaonptcosratio[ic][1]->SetTitle(Form("K%s PassTPC/MC",charge[ic].c_str()));
      kaonptcosratio[ic][1]->GetYaxis()->SetTitle("cos(#theta*)");
      kaonptcosratio[ic][1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      kaonptcosratio[ic][1]->Draw("Colz");

      ckaon->cd(2+4*ic);
      kaonptcosratio[ic][2]->SetTitle(Form("K%s ...+PassTOF/MC",charge[ic].c_str()));
      kaonptcosratio[ic][2]->GetYaxis()->SetTitle("cos(#theta*)");
      kaonptcosratio[ic][2]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      kaonptcosratio[ic][2]->Draw("Colz");

      ckaon->cd(3+4*ic);
      kaonptcosratio[ic][3]->SetTitle(Form("K%s ...+PassNsigKaon/MC",charge[ic].c_str()));
      kaonptcosratio[ic][3]->GetYaxis()->SetTitle("cos(#theta*)");
      kaonptcosratio[ic][3]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      kaonptcosratio[ic][3]->Draw("Colz");

      ckaon->cd(4+4*ic);
      kaonptcosratio[ic][4]->SetTitle(Form("K%s ...+PassM2/MC",charge[ic].c_str()));
      kaonptcosratio[ic][4]->GetYaxis()->SetTitle("cos(#theta*)");
      kaonptcosratio[ic][4]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      kaonptcosratio[ic][4]->Draw("Colz");
    }
    ckaon->SaveAs(Form("figures/KaonTTrees/%sKaonPtCosThetaHelicity.pdf",folderopt.c_str()));

    for(int ic = 0; ic < 2; ic++) 
    { 
      ckaon->cd(1+4*ic);
      kaoncosratio[ic][1]->SetTitle(Form("K%s PassTPC/MC",charge[ic].c_str()));
      kaoncosratio[ic][1]->GetXaxis()->SetTitle("cos(#theta*)");
      kaoncosratio[ic][1]->GetYaxis()->SetTitle("Eff");
      kaoncosratio[ic][1]->Draw("pE");

      ckaon->cd(2+4*ic);
      kaoncosratio[ic][2]->SetTitle(Form("K%s ...+PassTOF/MC",charge[ic].c_str()));
      kaoncosratio[ic][2]->GetXaxis()->SetTitle("cos(#theta*)");
      kaoncosratio[ic][2]->GetYaxis()->SetTitle("Eff");
      kaoncosratio[ic][2]->Draw("pE");

      ckaon->cd(3+4*ic);
      kaoncosratio[ic][3]->SetTitle(Form("K%s ...+PassNsigKaon/MC",charge[ic].c_str()));
      kaoncosratio[ic][3]->GetXaxis()->SetTitle("cos(#theta*)");
      kaoncosratio[ic][3]->GetYaxis()->SetTitle("Eff");
      kaoncosratio[ic][3]->Draw("pE");

      ckaon->cd(4+4*ic);
      kaoncosratio[ic][4]->SetTitle(Form("K%s ...+PassM2/MC",charge[ic].c_str()));
      kaoncosratio[ic][4]->GetXaxis()->SetTitle("cos(#theta*)");
      kaoncosratio[ic][4]->GetYaxis()->SetTitle("Eff");
      kaoncosratio[ic][4]->Draw("pE");
    }
    ckaon->SaveAs(Form("figures/KaonTTrees/%sKaon1DCosThetaHelicity.pdf",folderopt.c_str()));

    TCanvas *ckaon2 = new TCanvas("ckaon2", "ckaon2", 2000, 800);
    ckaon2->Divide(5,2);
    for(int i = 0; i < 10; i++)
    {
      ckaon2->cd(i+1);
      ckaon2->cd(i+1)->SetLeftMargin(0.12);
      ckaon2->cd(i+1)->SetRightMargin(0.15);
      ckaon2->cd(i+1)->SetBottomMargin(0.12);
      ckaon2->cd(i+1)->SetTicks(1,1);
      ckaon2->cd(i+1)->SetGrid(0,0); 
    }

    for(int ic = 0; ic < 2; ic++) 
    { 
      ckaon2->cd(1+5*ic);
      kaonptcos[ic][0]->SetTitle(Form("K%s MC",charge[ic].c_str()));
      kaonptcos[ic][0]->GetYaxis()->SetTitle("cos(#theta*)");
      kaonptcos[ic][0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      kaonptcos[ic][0]->Draw("Colz");

      ckaon2->cd(2+5*ic);
      kaonptcos[ic][1]->SetTitle(Form("K%s PassTPC",charge[ic].c_str()));
      kaonptcos[ic][1]->GetYaxis()->SetTitle("cos(#theta*)");
      kaonptcos[ic][1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      kaonptcos[ic][1]->Draw("Colz");

      ckaon2->cd(3+5*ic);
      kaonptcos[ic][2]->SetTitle(Form("K%s ...+PassTOF",charge[ic].c_str()));
      kaonptcos[ic][2]->GetYaxis()->SetTitle("cos(#theta*)");
      kaonptcos[ic][2]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      kaonptcos[ic][2]->Draw("Colz");

      ckaon2->cd(4+5*ic);
      kaonptcos[ic][3]->SetTitle(Form("K%s ...+PassNsigKaon",charge[ic].c_str()));
      kaonptcos[ic][3]->GetYaxis()->SetTitle("cos(#theta*)");
      kaonptcos[ic][3]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      kaonptcos[ic][3]->Draw("Colz");

      ckaon2->cd(5+5*ic);
      kaonptcos[ic][4]->SetTitle(Form("K%s ...+PassM2",charge[ic].c_str()));
      kaonptcos[ic][4]->GetYaxis()->SetTitle("cos(#theta*)");
      kaonptcos[ic][4]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      kaonptcos[ic][4]->Draw("Colz");
    }
    ckaon2->SaveAs(Form("figures/KaonTTrees/%sKaonPtCosThetaHelicity_Distributions.pdf",folderopt.c_str()));

    for(int ic = 0; ic < 2; ic++) 
    { 
      ckaon2->cd(1+5*ic);
      kaoncos[ic][0]->SetTitle(Form("K%s MC",charge[ic].c_str()));
      kaoncos[ic][0]->GetXaxis()->SetTitle("cos(#theta*)");
      kaoncos[ic][0]->GetYaxis()->SetTitle("Counts");
      kaoncos[ic][0]->Draw("pE");

      ckaon2->cd(2+5*ic);
      kaoncos[ic][1]->SetTitle(Form("K%s PassTPC",charge[ic].c_str()));
      kaoncos[ic][1]->GetXaxis()->SetTitle("cos(#theta*)");
      kaoncos[ic][1]->GetYaxis()->SetTitle("Counts");
      kaoncos[ic][1]->Draw("pE");

      ckaon2->cd(3+5*ic);
      kaoncos[ic][2]->SetTitle(Form("K%s ...+PassTOF",charge[ic].c_str()));
      kaoncos[ic][2]->GetXaxis()->SetTitle("cos(#theta*)");
      kaoncos[ic][2]->GetYaxis()->SetTitle("Counts");
      kaoncos[ic][2]->Draw("pE");

      ckaon2->cd(4+5*ic);
      kaoncos[ic][3]->SetTitle(Form("K%s ...+PassNsigKaon",charge[ic].c_str()));
      kaoncos[ic][3]->GetXaxis()->SetTitle("cos(#theta*)");
      kaoncos[ic][3]->GetYaxis()->SetTitle("Counts");
      kaoncos[ic][3]->Draw("pE");

      ckaon2->cd(5+5*ic);
      kaoncos[ic][4]->SetTitle(Form("K%s ...+PassM2",charge[ic].c_str()));
      kaoncos[ic][4]->GetXaxis()->SetTitle("cos(#theta*)");
      kaoncos[ic][4]->GetYaxis()->SetTitle("Counts");
      kaoncos[ic][4]->Draw("pE");
    }
    ckaon2->SaveAs(Form("figures/KaonTTrees/%sKaon1DCosThetaHelicity_Distributions.pdf",folderopt.c_str()));
 
    //TFile *Correction_Output = new TFile(Form("PhiEmbedding_%s_2DCorrections.root",vmsa::mBeamEnergy[energy].c_str()),"RECREATE");
    //Correction_Output->cd();
    //for(int ipt = vmsa::pt_rebin_first_2D[energy]; ipt < vmsa::pt_rebin_last_2D[energy]; ipt++)
    //{
    //  for(int i = 0; i < 2; i++)
    //  {
    //    pt_correction_global[ipt][i]->Write();
    //    pt_correction_helicity[ipt][i]->Write();
    //  }
    //}
    //Correction_Output->Close();
}

