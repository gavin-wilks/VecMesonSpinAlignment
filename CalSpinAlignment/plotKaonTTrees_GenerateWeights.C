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

void plotKaonTTrees_GenerateWeights(int energy = 4) 
{
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);    

    //double resolution[9]    = {0.221173,    0.296191,    0.411608,    0.535956,    0.630288,    0.675484,    0.648544,    0.539959,    0.377757  };  
    //double resolutionerr[9] = {0.000461111, 0.000319516, 0.000202537, 0.000141295, 0.000102773, 8.85656e-05, 9.71308e-05, 0.000193934, 0.00031656};


    TFile *inputweightsptyphi = TFile::Open(Form("phiweights/pt_y_phi_weights_v2_FunctionEmbed.root"));

    TH3F *h_mpTyv2Weights[9];
    for(int i_cent = 2; i_cent < 6; ++i_cent)
    {
      std::string HistName = Form("h_ratio_cent%d",i_cent);
      h_mpTyv2Weights[i_cent] = (TH3F*) inputweightsptyphi->Get(HistName.c_str());
    }

    

    // INITIALIZE TTREES 
    TFile *mKaonFileSE = TFile::Open("../data/Yields_Phi_SE_19GeV_20240830_kaontree.root");
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

    TFile *mKaonFileME = TFile::Open("../data/Yields_Phi_ME_19GeV_20240830_kaontree.root");
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
    
    Bool_t mPassTpc;
    Bool_t mPassTof;
    Bool_t mPassM2;
    Bool_t mPassNsig;


    //TFile *mKaonFileEmbed = TFile::Open("effaccfiles/Phi/19GeV/PhiEmbedding_kaontrees/PhiEmbedding_kaontrees_fixed.root");
    TFile *mKaonFileEmbed = TFile::Open("effaccfiles/Phi/19GeV/PhiEmbedding_kaontrees/PhiEmbedding_kaontrees_notofmatchreq_EP.root");
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
    
    mKaonTreeEmbed->SetBranchAddress("passtpc", &mPassTpc);
    mKaonTreeEmbed->SetBranchAddress("passtof", &mPassTof);
    mKaonTreeEmbed->SetBranchAddress("passm2", &mPassM2);
    mKaonTreeEmbed->SetBranchAddress("passns", &mPassNsig);
    // INITIALIZE TTREES 


    // INITIALIZE HISTOGRAMS
    TH1F* pt[5][3];
    TH1F* rapidity[5][3]; 
    TH1F* eta[5][3];
    TH1F* phi[5][3];

    THnF* ptyetaphi[5][3];

    TH1F* ratiopt[3];
    TH1F* ratiorapidity[3]; 
    TH1F* ratioeta[3];
    TH1F* ratiophi[3];
     
    TH2F* ptrapidity[5][3];
    TH2F* pteta[5][3]; 
    TH2F* ptphi[5][3];
    TH2F* rapidityeta[5][3];
    TH2F* rapidityphi[5][3];
    TH2F* etaphi[5][3]; 

    TH2F* ratioptrapidity[3];
    TH2F* ratiopteta[3]; 
    TH2F* ratioptphi[3];
    TH2F* ratiorapidityeta[3];
    TH2F* ratiorapidityphi[3];
    TH2F* ratioetaphi[3]; 
  
    // Tracking information
    TH1F* nhitsfit[5][3]; 
    TH2F* nhitsfitpt[5][3]; 
    TH2F* nhitsfitrapidity[5][3]; 
    TH2F* nhitsfiteta[5][3]; 
    TH2F* nhitsfitphi[5][3]; 

    TH1F* nhitsmax[5][3]; 
    TH2F* nhitsmaxpt[5][3]; 
    TH2F* nhitsmaxrapidity[5][3]; 
    TH2F* nhitsmaxeta[5][3]; 
    TH2F* nhitsmaxphi[5][3]; 

    TH1F* nhitsratio[5][3];
    TH2F* nhitsratiopt[5][3];
    TH2F* nhitsratiorapidity[5][3];
    TH2F* nhitsratioeta[5][3];
    TH2F* nhitsratiophi[5][3];

    TH1F* dedx[5][3];
    TH2F* dedxpt[5][3];
    TH2F* dedxrapidity[5][3];
    TH2F* dedxeta[5][3];
    TH2F* dedxphi[5][3];

    TH1F* dca[5][3];
    TH2F* dcapt[5][3];
    TH2F* dcarapidity[5][3];
    TH2F* dcaeta[5][3];
    TH2F* dcaphi[5][3];
     
    TH1F* rationhitsfit[2]; 
    TH2F* rationhitsfitpt[2]; 
    TH2F* rationhitsfitrapidity[2]; 
    TH2F* rationhitsfiteta[2]; 
    TH2F* rationhitsfitphi[2]; 

    TH1F* rationhitsmax[2]; 
    TH2F* rationhitsmaxpt[2]; 
    TH2F* rationhitsmaxrapidity[2]; 
    TH2F* rationhitsmaxeta[2]; 
    TH2F* rationhitsmaxphi[2]; 

    TH1F* rationhitsratio[2];
    TH2F* rationhitsratiopt[2];
    TH2F* rationhitsratiorapidity[2];
    TH2F* rationhitsratioeta[2];
    TH2F* rationhitsratiophi[2];

    TH1F* ratiodedx[2];
    TH2F* ratiodedxpt[2];
    TH2F* ratiodedxrapidity[2];
    TH2F* ratiodedxeta[2];
    TH2F* ratiodedxphi[2];

    TH1F* ratiodca[2];
    TH2F* ratiodcapt[2];
    TH2F* ratiodcarapidity[2];
    TH2F* ratiodcaeta[2];
    TH2F* ratiodcaphi[2];
     
    //TH3F* ptrapidityphi[9];
    //TH3F* ptrapidityphiDesired[9];
    //TH3F* ptrapidityphiRatio[9];
 
    const int npt = 70;
    const float ptmin = 0.0, ptmax = 3.5; 
    const int ny = 50;
    const float ymin = -1.5, ymax = 1.5;
    const int neta = 50;
    const float etamin = -1.5, etamax = 1.5;
    const int nphi = 120;
    const float phimin = -TMath::Pi(), phimax = TMath::Pi();

    Int_t bins[4] = {npt, ny, neta, nphi};
    Double_t xmin[4] = {ptmin, ymin, etamin, phimin};
    Double_t xmax[4] = {ptmax, ymax, etamax, phimax};

    const int nhf = 101;
    const float nhfmin = -0.5, nhfmax = 100.5; 
    const int nhm = 101;
    const float nhmmin = -0.5, nhmmax = 100.5; 
    const int nhr = 25;
    const float nhrmin = 0.5, nhrmax = 1.1; 
    const int nde = 250;
    const float demin = 0.0, demax = 10.0; 
    const int ndca = 150;
    const float dcamin = 0.0, dcamax = 2.1; 

    //TH3F *h_mPtYPhiPsi[9];
    //TH3F *h_mPtYPhiPsiDesired[9];
    //TH3F *h_mPtYPhiPsiWeight[9];
    //for(int i_cent = 2; i_cent < 6; ++i_cent)
    //{
    //  std::string HistName = Form("h_embed_cent%d",i_cent);
    //  TH3F *h_mPtYPhiPsi[i_cent] = new TH3F(HistName.c_str(),HistName.c_str(), 50, 1.2, 4.2, 30, -1.0, 1.0, 50, -Tmath: );
    //  HistName = Form("h_desired_cent%d",i_cent);
    //  TH3F *h_mPtYPhiPsiDesired[i_cent];
    //}


    string charge[3] = {"plus","minus","phi"};
    string mixing[5] = {"SE","ME","SM","RC","MC"};

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

        hist = Form("ptyetaphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        ptyetaphi[im][ic] = new THnF(hist.c_str(), hist.c_str(), 4, bins, xmin, xmax);

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

        // TH1F Basic Variables
        hist = Form("nhitsfit_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsfit[im][ic] = new TH1F(hist.c_str(), hist.c_str(), nhf, nhfmin, nhfmax);
        hist = Form("nhitsmax_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsmax[im][ic] = new TH1F(hist.c_str(), hist.c_str(), nhm, nhmmin, nhmmax);
        hist = Form("nhitsratio_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsratio[im][ic] = new TH1F(hist.c_str(), hist.c_str(), nhr, nhrmin, nhrmax);
        hist = Form("dedx_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        dedx[im][ic] = new TH1F(hist.c_str(), hist.c_str(), nde, demin, demax);
        hist = Form("dca_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        dca[im][ic] = new TH1F(hist.c_str(), hist.c_str(), ndca, dcamin, dcamax);
        // TH1F Basic Variables

        // TH2F Basic Vairables PT        
        hist = Form("nhitsfitpt_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsfitpt[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nhf, nhfmin, nhfmax);
        hist = Form("nhitsmaxpt_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsmaxpt[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nhm, nhmmin, nhmmax);
        hist = Form("nhitsratiopt_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsratiopt[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nhr, nhrmin, nhrmax);
        hist = Form("dedxpt_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        dedxpt[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nde, demin, demax);
        hist = Form("dcapt_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        dcapt[im][ic] = new TH2F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, ndca, dcamin, dcamax);

        // TH2F Basic Vairables rapidity        
        hist = Form("nhitsfitrapidity_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsfitrapidity[im][ic] = new TH2F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nhf, nhfmin, nhfmax);
        hist = Form("nhitsmaxrapidity_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsmaxrapidity[im][ic] = new TH2F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nhm, nhmmin, nhmmax);
        hist = Form("nhitsratiorapidity_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsratiorapidity[im][ic] = new TH2F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nhr, nhrmin, nhrmax);
        hist = Form("dedxrapidity_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        dedxrapidity[im][ic] = new TH2F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nde, demin, demax);
        hist = Form("dcarapidity_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        dcarapidity[im][ic] = new TH2F(hist.c_str(), hist.c_str(), ny, ymin, ymax, ndca, dcamin, dcamax);

        // TH2F Basic Vairables eta        
        hist = Form("nhitsfiteta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsfiteta[im][ic] = new TH2F(hist.c_str(), hist.c_str(), neta, etamin, etamax, nhf, nhfmin, nhfmax);
        hist = Form("nhitsmaxeta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsmaxeta[im][ic] = new TH2F(hist.c_str(), hist.c_str(), neta, etamin, etamax, nhm, nhmmin, nhmmax);
        hist = Form("nhitsratioeta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsratioeta[im][ic] = new TH2F(hist.c_str(), hist.c_str(), neta, etamin, etamax, nhr, nhrmin, nhrmax);
        hist = Form("dedxeta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        dedxeta[im][ic] = new TH2F(hist.c_str(), hist.c_str(), neta, etamin, etamax, nde, demin, demax);
        hist = Form("dcaeta_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        dcaeta[im][ic] = new TH2F(hist.c_str(), hist.c_str(), neta, etamin, etamax, ndca, dcamin, dcamax);


        // TH2F Basic Vairables phi        
        hist = Form("nhitsfitphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsfitphi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), nphi, phimin, phimax, nhf, nhfmin, nhfmax);
        hist = Form("nhitsmaxphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsmaxphi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), nphi, phimin, phimax, nhm, nhmmin, nhmmax);
        hist = Form("nhitsratiophi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        nhitsratiophi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), nphi, phimin, phimax, nhr, nhrmin, nhrmax);
        hist = Form("dedxphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        dedxphi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), nphi, phimin, phimax, nde, demin, demax);
        hist = Form("dcaphi_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        dcaphi[im][ic] = new TH2F(hist.c_str(), hist.c_str(), nphi, phimin, phimax, ndca, dcamin, dcamax);


      }
    }
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
      pt[0][mChargeSE]->Fill(mPtSE,mWeightSE);
      rapidity[0][mChargeSE]->Fill(mRapiditySE,mWeightSE);
      eta[0][mChargeSE]->Fill(mEtaSE,mWeightSE);
      phi[0][mChargeSE]->Fill(mPhiSE,mWeightSE);
      // TH1F 

      ptyetaphi[0][mChargeSE]->Fill(mPtSE,mRapiditySE,mEtaSE,mPhiSE,mWeightSE);

      // TH2F 
      ptrapidity[0][mChargeSE]->Fill(mPtSE,mRapiditySE,mWeightSE);
      pteta[0][mChargeSE]->Fill(mPtSE,mEtaSE,mWeightSE);
      ptphi[0][mChargeSE]->Fill(mPtSE,mPhiSE,mWeightSE);
      rapidityeta[0][mChargeSE]->Fill(mRapiditySE,mEtaSE,mWeightSE);
      rapidityphi[0][mChargeSE]->Fill(mRapiditySE,mPhiSE,mWeightSE);
      etaphi[0][mChargeSE]->Fill(mEtaSE,mPhiSE,mWeightSE);
      // TH2F 

      nhitsfit[0][mChargeSE]->Fill(mNHitsFitSE,mWeightSE); 
      nhitsmax[0][mChargeSE]->Fill(mNHitsMaxSE,mWeightSE); 
      nhitsratio[0][mChargeSE]->Fill(float(mNHitsFitSE)/float(mNHitsMaxSE),mWeightSE);
      dedx[0][mChargeSE]->Fill(mDEdxSE,mWeightSE);
      dca[0][mChargeSE]->Fill(mDcaSE,mWeightSE);

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
      pt[1][mChargeME]->Fill(mPtME,mWeightME);
      rapidity[1][mChargeME]->Fill(mRapidityME,mWeightME);
      eta[1][mChargeME]->Fill(mEtaME,mWeightME);
      phi[1][mChargeME]->Fill(mPhiME,mWeightME);
      // TH1F 

      ptyetaphi[1][mChargeME]->Fill(mPtME,mRapidityME,mEtaME,mPhiME,mWeightME);

      // TH2F 
      ptrapidity[1][mChargeME]->Fill(mPtME,mRapidityME,mWeightME);
      pteta[1][mChargeME]->Fill(mPtME,mEtaME,mWeightME);
      ptphi[1][mChargeME]->Fill(mPtME,mPhiME,mWeightME);
      rapidityeta[1][mChargeME]->Fill(mRapidityME,mEtaME,mWeightME);
      rapidityphi[1][mChargeME]->Fill(mRapidityME,mPhiME,mWeightME);
      etaphi[1][mChargeME]->Fill(mEtaME,mPhiME,mWeightME);
      // TH2F 

      nhitsfit[1][mChargeME]->Fill(mNHitsFitME,mWeightME); 
      nhitsmax[1][mChargeME]->Fill(mNHitsMaxME,mWeightME); 
      nhitsratio[1][mChargeME]->Fill(float(mNHitsFitME)/float(mNHitsMaxME),mWeightME);
      dedx[1][mChargeME]->Fill(mDEdxME,mWeightME);
      dca[1][mChargeME]->Fill(mDcaME,mWeightME);
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

      // TH1F 
      pt[4][0]->Fill(mMcPtP,mWeight);
      rapidity[4][0]->Fill(mMcRapidityP,mWeight);
      eta[4][0]->Fill(mMcEtaP,mWeight);
      phi[4][0]->Fill(mMcPhiP,mWeight);
      pt[4][1]->Fill(mMcPtM,mWeight);
      rapidity[4][1]->Fill(mMcRapidityM,mWeight);
      eta[4][1]->Fill(mMcEtaM,mWeight);
      phi[4][1]->Fill(mMcPhiM,mWeight);

      pt[4][2]->Fill(mMcPt,mWeight);
      rapidity[4][2]->Fill(mMcRapidity,mWeight);
      eta[4][2]->Fill(mMcEta,mWeight);
      phi[4][2]->Fill(mMcPhi,mWeight);
      // TH1F 

      // TH2F 
      ptrapidity[4][0]->Fill(mMcPtP,mMcRapidityP,mWeight);
      pteta[4][0]->Fill(mMcPtP,mMcEtaP,mWeight);
      ptphi[4][0]->Fill(mMcPtP,mMcPhiP,mWeight);
      rapidityeta[4][0]->Fill(mMcRapidityP,mMcEtaP,mWeight);
      rapidityphi[4][0]->Fill(mMcRapidityP,mMcPhiP,mWeight);
      etaphi[4][0]->Fill(mMcEtaP,mMcPhiP,mWeight);
      ptrapidity[4][1]->Fill(mMcPtM,mMcRapidityM,mWeight);
      pteta[4][1]->Fill(mMcPtM,mMcEtaM,mWeight);
      ptphi[4][1]->Fill(mMcPtM,mMcPhiM,mWeight);
      rapidityeta[4][1]->Fill(mMcRapidityM,mMcEtaM,mWeight);
      rapidityphi[4][1]->Fill(mMcRapidityM,mMcPhiM,mWeight);
      etaphi[4][1]->Fill(mMcEtaM,mMcPhiM,mWeight);

      ptrapidity[4][2]->Fill(mMcPt,mMcRapidity,mWeight);
      pteta[4][2]->Fill(mMcPt,mMcEta,mWeight);
      ptphi[4][2]->Fill(mMcPt,mMcPhi,mWeight);
      rapidityeta[4][2]->Fill(mMcRapidity,mMcEta,mWeight);
      rapidityphi[4][2]->Fill(mMcRapidity,mMcPhi,mWeight);
      etaphi[4][2]->Fill(mMcEta,mMcPhi,mWeight);
      // TH2F 

      //cout << "AM I HERE? " << endl;
      //cout << "mCent = " << mCent << endl;


      //cout << " I DID GET HERE" << endl;

      //cout << "mEpFull = " << mEpFull << endl;

      float mMcPhiPsi = mMcPhi - mEpFull;
      while(mMcPhiPsi < 0.0) mMcPhiPsi += 2.0*TMath::Pi();
      while(mMcPhiPsi > 2.0*TMath::Pi()) mMcPhiPsi -= 2.0*TMath::Pi();
     
      //cout << "mMcPhiPsi = " << mMcPhiPsi << endl;

      float PhiWeight = h_mpTyv2Weights[mCent]->GetBinContent(h_mpTyv2Weights[mCent]->FindBin(mMcPt,mMcRapidity,mMcPhiPsi));

      //cout << "mPassTpc  = " << mPassTpc  << endl;
      //cout << "mPassTof  = " << mPassTof  << endl;
      //cout << "mPassM2   = " << mPassM2   << endl;
      //cout << "mPassNsig = " << mPassNsig << endl;
 
      if(mRcExists && mPassTpc && mPassTof && mPassM2 && mPassNsig)    
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

        //cout << "mRcPtP       = " << mRcPtP       << endl;  
        //cout << "mRcRapidityP = " << mRcRapidityP << endl;     
        //cout << "mRcEtaP      = " << mRcEtaP      << endl;
        //cout << "mRcPhiP      = " << mRcPhiP      << endl;

        //cout << "mRcPtM       = " << mRcPtM       << endl;  
        //cout << "mRcRapidityM = " << mRcRapidityM << endl;      
        //cout << "mRcEtaM      = " << mRcEtaM      << endl;
        //cout << "mRcPhiM      = " << mRcPhiM      << endl;

        // TH1F 
        pt[3][0]->Fill(mRcPtP,mWeight*PhiWeight);
        rapidity[3][0]->Fill(mRcRapidityP,mWeight*PhiWeight);
        eta[3][0]->Fill(mRcEtaP,mWeight*PhiWeight);
        phi[3][0]->Fill(mRcPhiP,mWeight*PhiWeight);
        pt[3][1]->Fill(mRcPtM,mWeight*PhiWeight);
        rapidity[3][1]->Fill(mRcRapidityM,mWeight*PhiWeight);
        eta[3][1]->Fill(mRcEtaM,mWeight*PhiWeight);
        phi[3][1]->Fill(mRcPhiM,mWeight*PhiWeight);

        pt[3][2]->Fill(mRcPt,mWeight*PhiWeight);
        rapidity[3][2]->Fill(mRcRapidity,mWeight*PhiWeight);
        eta[3][2]->Fill(mRcEta,mWeight*PhiWeight);
        phi[3][2]->Fill(mRcPhi,mWeight*PhiWeight);
        // TH1F 

        // TH2F 
        ptrapidity[3][0]->Fill(mRcPtP,mRcRapidityP,mWeight*PhiWeight);
        pteta[3][0]->Fill(mRcPtP,mRcEtaP,mWeight*PhiWeight);
        ptphi[3][0]->Fill(mRcPtP,mRcPhiP,mWeight*PhiWeight);
        rapidityeta[3][0]->Fill(mRcRapidityP,mRcEtaP,mWeight*PhiWeight);
        rapidityphi[3][0]->Fill(mRcRapidityP,mRcPhiP,mWeight*PhiWeight);
        etaphi[3][0]->Fill(mRcEtaP,mRcPhiP,mWeight*PhiWeight);
        ptrapidity[3][1]->Fill(mRcPtM,mRcRapidityM,mWeight*PhiWeight);
        pteta[3][1]->Fill(mRcPtM,mRcEtaM,mWeight*PhiWeight);
        ptphi[3][1]->Fill(mRcPtM,mRcPhiM,mWeight*PhiWeight);
        rapidityeta[3][1]->Fill(mRcRapidityM,mRcEtaM,mWeight*PhiWeight);
        rapidityphi[3][1]->Fill(mRcRapidityM,mRcPhiM,mWeight*PhiWeight);
        etaphi[3][1]->Fill(mRcEtaM,mRcPhiM,mWeight*PhiWeight);

        ptrapidity[3][2]->Fill(mRcPt,mRcRapidity,mWeight*PhiWeight);
        pteta[3][2]->Fill(mRcPt,mRcEta,mWeight*PhiWeight);
        ptphi[3][2]->Fill(mRcPt,mRcPhi,mWeight*PhiWeight);
        rapidityeta[3][2]->Fill(mRcRapidity,mRcEta,mWeight*PhiWeight);
        rapidityphi[3][2]->Fill(mRcRapidity,mRcPhi,mWeight*PhiWeight);
        etaphi[3][2]->Fill(mRcEta,mRcPhi,mWeight*PhiWeight);
        // TH2F 

        nhitsfit[3][0]->Fill(mNHitsFitP,mWeight*PhiWeight); 
        nhitsmax[3][0]->Fill(mNHitsMaxP,mWeight*PhiWeight); 
        nhitsratio[3][0]->Fill(float(mNHitsFitP)/float(mNHitsMaxP),mWeight*PhiWeight);
        dedx[3][0]->Fill(mDEdxP,mWeight*PhiWeight);
        dca[3][0]->Fill(mDcaP,mWeight*PhiWeight);
        nhitsfit[3][1]->Fill(mNHitsFitM,mWeight*PhiWeight); 
        nhitsmax[3][1]->Fill(mNHitsMaxM,mWeight*PhiWeight); 
        nhitsratio[3][1]->Fill(float(mNHitsFitM)/float(mNHitsMaxM),mWeight*PhiWeight);
        dedx[3][1]->Fill(mDEdxM,mWeight*PhiWeight);
        dca[3][1]->Fill(mDcaM,mWeight*PhiWeight);
      }
    }

 
    for(int ic = 0; ic < 2; ic++)
    {
      ptyetaphi[2][ic] = (THnF*) ptyetaphi[0][ic]->Clone();
      ptyetaphi[2][ic]->Add(ptyetaphi[1][ic],-1.0);

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

      pt[2][ic] = (TH1F*) ptyetaphi[2][ic]->Projection(0);
      //pt[2][ic]->Add(pt[1][ic],-1.0);
      rapidity[2][ic] = (TH1F*) ptyetaphi[2][ic]->Projection(1);
      //rapidity[2][ic]->Add(rapidity[1][ic],-1.0);
      eta[2][ic] = (TH1F*) ptyetaphi[2][ic]->Projection(2);
      //eta[2][ic]->Add(eta[1][ic],-1.0);
      phi[2][ic] = (TH1F*) ptyetaphi[2][ic]->Projection(3);
      //phi[2][ic]->Add(phi[1][ic],-1.0);

      ptrapidity[2][ic] = (TH2F*) ptyetaphi[2][ic]->Projection(1,0);
      //ptrapidity[2][ic]->Add(ptrapidity[1][ic],-1.0);
      pteta[2][ic] = (TH2F*) ptyetaphi[2][ic]->Projection(2,0);
      //pteta[2][ic]->Add(pteta[1][ic],-1.0);
      ptphi[2][ic] = (TH2F*) ptyetaphi[2][ic]->Projection(3,0);
      //ptphi[2][ic]->Add(ptphi[1][ic],-1.0);
      rapidityeta[2][ic] = (TH2F*) ptyetaphi[2][ic]->Projection(2,1);
      //rapidityeta[2][ic]->Add(rapidityeta[1][ic],-1.0);
      rapidityphi[2][ic] = (TH2F*) ptyetaphi[2][ic]->Projection(3,1);
      //rapidityphi[2][ic]->Add(rapidityphi[1][ic],-1.0);
      etaphi[2][ic] = (TH2F*) ptyetaphi[2][ic]->Projection(3,2);
      //etaphi[2][ic]->Add(etaphi[1][ic],-1.0);

      nhitsfit[2][ic] = (TH1F*) nhitsfit[0][ic]->Clone();
      nhitsfit[2][ic]->Add(nhitsfit[1][ic],-1.0);
      nhitsmax[2][ic] = (TH1F*) nhitsmax[0][ic]->Clone();
      nhitsmax[2][ic]->Add(nhitsmax[1][ic],-1.0);
      nhitsratio[2][ic] = (TH1F*) nhitsratio[0][ic]->Clone();     
      nhitsratio[2][ic]->Add(nhitsratio[1][ic],-1.0);
      dedx[2][ic] = (TH1F*) dedx[0][ic]->Clone();
      dedx[2][ic]->Add(dedx[1][ic],-1.0);
      dca[2][ic] = (TH1F*) dca[0][ic]->Clone();
      dca[2][ic]->Add(dca[1][ic],-1.0);


      //Ratios
      double integralData = pt[2][ic]->Integral(0,-1);
      double integralRC = pt[3][ic]->Integral(0,-1);
      double RCData = integralRC/integralData;

      rationhitsfit[ic] = (TH1F*) nhitsfit[2][ic]->Clone();
      rationhitsfit[ic]->Scale(RCData);
      rationhitsfit[ic]->Divide(nhitsfit[3][ic]);

      rationhitsmax[ic] = (TH1F*) nhitsmax[2][ic]->Clone();
      rationhitsmax[ic]->Scale(RCData);
      rationhitsmax[ic]->Divide(nhitsmax[3][ic]);

      rationhitsratio[ic] = (TH1F*) nhitsratio[2][ic]->Clone();
      rationhitsratio[ic]->Scale(RCData);
      rationhitsratio[ic]->Divide(nhitsratio[3][ic]);

      ratiodedx[ic] = (TH1F*) dedx[2][ic]->Clone();
      ratiodedx[ic]->Scale(RCData);
      ratiodedx[ic]->Divide(dedx[3][ic]);

      ratiodca[ic] = (TH1F*) dca[2][ic]->Clone();
      ratiodca[ic]->Scale(RCData);
      ratiodca[ic]->Divide(dca[3][ic]);

      ratiopt[ic] = (TH1F*) pt[2][ic]->Clone();
      ratiopt[ic]->Scale(RCData);
      ratiopt[ic]->Divide(pt[3][ic]);

      ratiorapidity[ic] = (TH1F*) rapidity[2][ic]->Clone();
      ratiorapidity[ic]->Scale(RCData);
      ratiorapidity[ic]->Divide(rapidity[3][ic]);

      ratioeta[ic] = (TH1F*) eta[2][ic]->Clone();
      ratioeta[ic]->Scale(RCData);
      ratioeta[ic]->Divide(eta[3][ic]);

      ratiophi[ic] = (TH1F*) phi[2][ic]->Clone();
      ratiophi[ic]->Scale(RCData);
      ratiophi[ic]->Divide(phi[3][ic]);

      ratioptrapidity[ic] = (TH2F*) ptrapidity[2][ic]->Clone();
      ratioptrapidity[ic]->Scale(RCData);
      ratioptrapidity[ic]->Divide(ptrapidity[3][ic]);

      ratiopteta[ic] = (TH2F*) pteta[2][ic]->Clone();
      ratiopteta[ic]->Scale(RCData);
      ratiopteta[ic]->Divide(pteta[3][ic]);

      ratioptphi[ic] = (TH2F*) ptphi[2][ic]->Clone();
      ratioptphi[ic]->Scale(RCData);
      ratioptphi[ic]->Divide(ptphi[3][ic]);

      ratiorapidityeta[ic] = (TH2F*) rapidityeta[2][ic]->Clone();
      ratiorapidityeta[ic]->Scale(RCData);
      ratiorapidityeta[ic]->Divide(rapidityeta[3][ic]);

      ratiorapidityphi[ic] = (TH2F*) rapidityphi[2][ic]->Clone();
      ratiorapidityphi[ic]->Scale(RCData);
      ratiorapidityphi[ic]->Divide(rapidityphi[3][ic]);

      ratioetaphi[ic] = (TH2F*) etaphi[2][ic]->Clone();
      ratioetaphi[ic]->Scale(RCData);
      ratioetaphi[ic]->Divide(etaphi[3][ic]);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 2000, 800);
    c1->Divide(5,2);
    for(int i = 0; i < 10; i++)
    {
      c1->cd(i+1);
      c1->cd(i+1)->SetLeftMargin(0.15);
      c1->cd(i+1)->SetRightMargin(0.15);
      c1->cd(i+1)->SetBottomMargin(0.15);
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
        c1->cd(1)->SetLogy(0);
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
        ptrapidity[im][ic]->Draw("Colz"); 

        c1->cd(6);
        //c1->cd(6)->SetLogz(1);
        pteta[im][ic]->SetTitle(histtitle.c_str()); 
        pteta[im][ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
        pteta[im][ic]->GetYaxis()->SetTitle("#eta"); 
        pteta[im][ic]->SetContour(100); 
        pteta[im][ic]->Draw("Colz"); 

        c1->cd(7);
        //c1->cd(7)->SetLogz(1);
        ptphi[im][ic]->SetTitle(histtitle.c_str()); 
        ptphi[im][ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
        ptphi[im][ic]->GetYaxis()->SetTitle("#phi"); 
        ptphi[im][ic]->SetContour(100); 
        ptphi[im][ic]->Draw("Colz"); 

        c1->cd(8);
        rapidityeta[im][ic]->SetTitle(histtitle.c_str()); 
        rapidityeta[im][ic]->GetXaxis()->SetTitle("y"); 
        rapidityeta[im][ic]->GetYaxis()->SetTitle("#eta"); 
        rapidityeta[im][ic]->SetContour(100); 
        rapidityeta[im][ic]->Draw("Colz"); 

        c1->cd(9);
        rapidityphi[im][ic]->SetTitle(histtitle.c_str()); 
        rapidityphi[im][ic]->GetXaxis()->SetTitle("y"); 
        rapidityphi[im][ic]->GetYaxis()->SetTitle("#phi"); 
        rapidityphi[im][ic]->SetContour(100); 
        rapidityphi[im][ic]->Draw("Colz"); 

        c1->cd(10);
        etaphi[im][ic]->SetTitle(histtitle.c_str()); 
        etaphi[im][ic]->GetXaxis()->SetTitle("#eta"); 
        etaphi[im][ic]->GetYaxis()->SetTitle("#phi"); 
        etaphi[im][ic]->SetContour(100); 
        etaphi[im][ic]->Draw("Colz"); 
   
        c1->SaveAs(Form("figures/KaonTTrees/K%s_%s.pdf",charge[ic].c_str(),mixinglabel[im].c_str()));
      }
    }

    for(int ic = 0; ic < 2; ic++)
    {
      string histtitle = Form("K%s Data/RC",charge[ic].c_str());

      c1->cd(1);
      c1->cd(1)->SetLogy(0);
      ratiopt[ic]->SetTitle(histtitle.c_str()); 
      ratiopt[ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
      ratiopt[ic]->GetYaxis()->SetTitle("Data/RC"); 
      ratiopt[ic]->Draw("pE"); 

      c1->cd(2);
      ratiorapidity[ic]->SetTitle(histtitle.c_str()); 
      ratiorapidity[ic]->GetXaxis()->SetTitle("y"); 
      ratiorapidity[ic]->GetYaxis()->SetTitle("Data/RC"); 
      ratiorapidity[ic]->Draw("pE"); 

      c1->cd(3);
      ratioeta[ic]->SetTitle(histtitle.c_str()); 
      ratioeta[ic]->GetXaxis()->SetTitle("#eta"); 
      ratioeta[ic]->GetYaxis()->SetTitle("Data/RC"); 
      ratioeta[ic]->Draw("pE"); 

      c1->cd(4);
      ratiophi[ic]->SetTitle(histtitle.c_str()); 
      ratiophi[ic]->GetXaxis()->SetTitle("#phi"); 
      ratiophi[ic]->GetYaxis()->SetTitle("Data/RC"); 
      ratiophi[ic]->Draw("pE"); 

      c1->cd(5);
      //c1->cd(5)->SetLogz(1);
      ratioptrapidity[ic]->SetTitle(histtitle.c_str()); 
      ratioptrapidity[ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
      ratioptrapidity[ic]->GetYaxis()->SetTitle("y"); 
      ratioptrapidity[ic]->SetContour(100); 
      ratioptrapidity[ic]->Draw("Colz"); 

      c1->cd(6);
      //c1->cd(6)->SetLogz(1);
      ratiopteta[ic]->SetTitle(histtitle.c_str()); 
      ratiopteta[ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
      ratiopteta[ic]->GetYaxis()->SetTitle("#eta"); 
      ratiopteta[ic]->SetContour(100); 
      ratiopteta[ic]->Draw("Colz"); 

      c1->cd(7);
      //c1->cd(7)->SetLogz(1);
      ratioptphi[ic]->SetTitle(histtitle.c_str()); 
      ratioptphi[ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
      ratioptphi[ic]->GetYaxis()->SetTitle("#phi"); 
      ratioptphi[ic]->SetContour(100); 
      ratioptphi[ic]->Draw("Colz"); 

      c1->cd(8);
      ratiorapidityeta[ic]->SetTitle(histtitle.c_str()); 
      ratiorapidityeta[ic]->GetXaxis()->SetTitle("y"); 
      ratiorapidityeta[ic]->GetYaxis()->SetTitle("#eta"); 
      ratiorapidityeta[ic]->SetContour(100); 
      ratiorapidityeta[ic]->Draw("Colz"); 

      c1->cd(9);
      ratiorapidityphi[ic]->SetTitle(histtitle.c_str()); 
      ratiorapidityphi[ic]->GetXaxis()->SetTitle("y"); 
      ratiorapidityphi[ic]->GetYaxis()->SetTitle("#phi"); 
      ratiorapidityphi[ic]->SetContour(100); 
      ratiorapidityphi[ic]->Draw("Colz"); 

      c1->cd(10);
      ratioetaphi[ic]->SetTitle(histtitle.c_str()); 
      ratioetaphi[ic]->GetXaxis()->SetTitle("#eta"); 
      ratioetaphi[ic]->GetYaxis()->SetTitle("#phi"); 
      ratioetaphi[ic]->SetContour(100); 
      ratioetaphi[ic]->Draw("Colz"); 
   
      c1->SaveAs(Form("figures/KaonTTrees/K%s_DataRCRatio.pdf",charge[ic].c_str()));
      
    }

    for(int im = 0; im < 4; im++)
    {
      for(int ic = 0; ic < 2; ic++)
      {
        string histtitle = Form("K%s, %s",charge[ic].c_str(),mixingtitle[im].c_str());

        c1->cd(1+5*ic);
        c1->cd(1+5*ic)->SetLogy(0);
        nhitsfit[im][ic]->SetTitle(histtitle.c_str()); 
        nhitsfit[im][ic]->GetXaxis()->SetTitle("NHitsFit"); 
        nhitsfit[im][ic]->GetYaxis()->SetTitle("Count"); 
        nhitsfit[im][ic]->Draw("pE"); 

        c1->cd(2+5*ic);
        c1->cd(2+5*ic)->SetLogy(0);
        nhitsmax[im][ic]->SetTitle(histtitle.c_str()); 
        nhitsmax[im][ic]->GetXaxis()->SetTitle("NHitsMax"); 
        nhitsmax[im][ic]->GetYaxis()->SetTitle("Count"); 
        nhitsmax[im][ic]->Draw("pE"); 
        
        c1->cd(3+5*ic);
        c1->cd(3+5*ic)->SetLogy(0);
        nhitsratio[im][ic]->SetTitle(histtitle.c_str()); 
        nhitsratio[im][ic]->GetXaxis()->SetTitle("NHitsFit/NHitsMax"); 
        nhitsratio[im][ic]->GetYaxis()->SetTitle("Count"); 
        nhitsratio[im][ic]->Draw("pE"); 

        c1->cd(4+5*ic);
        c1->cd(4+5*ic)->SetLogy(1);
        dedx[im][ic]->SetTitle(histtitle.c_str()); 
        dedx[im][ic]->GetXaxis()->SetTitle("dE/dx"); 
        dedx[im][ic]->GetYaxis()->SetTitle("Count"); 
        dedx[im][ic]->Draw("pE"); 

        c1->cd(5+5*ic);
        c1->cd(5+5*ic)->SetLogy(1);
        dca[im][ic]->SetTitle(histtitle.c_str()); 
        dca[im][ic]->GetXaxis()->SetTitle("|dca| (cm)"); 
        dca[im][ic]->GetYaxis()->SetTitle("Count"); 
        dca[im][ic]->Draw("pE"); 
      }
      c1->SaveAs(Form("figures/KaonTTrees/KaonBasicVariables_%s.pdf",mixinglabel[im].c_str()));
    }

    for(int ic = 0; ic < 2; ic++)
    {
      string histtitle = Form("K%s, Data/RC",charge[ic].c_str());

      c1->cd(1+5*ic);
      c1->cd(1+5*ic)->SetLogy(0);
      rationhitsfit[ic]->SetTitle(histtitle.c_str()); 
      rationhitsfit[ic]->GetXaxis()->SetTitle("NHitsFit"); 
      rationhitsfit[ic]->GetYaxis()->SetTitle("Data/RC"); 
      rationhitsfit[ic]->Draw("pE"); 

      c1->cd(2+5*ic);
      c1->cd(2+5*ic)->SetLogy(0);
      rationhitsmax[ic]->SetTitle(histtitle.c_str()); 
      rationhitsmax[ic]->GetXaxis()->SetTitle("NHitsMax"); 
      rationhitsmax[ic]->GetYaxis()->SetTitle("Data/RC"); 
      rationhitsmax[ic]->Draw("pE"); 
      
      c1->cd(3+5*ic);
      c1->cd(3+5*ic)->SetLogy(0);
      rationhitsratio[ic]->SetTitle(histtitle.c_str()); 
      rationhitsratio[ic]->GetXaxis()->SetTitle("NHitsFit/NHitsMax"); 
      rationhitsratio[ic]->GetYaxis()->SetTitle("Data/RC"); 
      rationhitsratio[ic]->Draw("pE"); 

      c1->cd(4+5*ic);
      c1->cd(4+5*ic)->SetLogy(0);
      ratiodedx[ic]->SetTitle(histtitle.c_str()); 
      ratiodedx[ic]->GetXaxis()->SetTitle("dE/dx"); 
      ratiodedx[ic]->GetYaxis()->SetTitle("Data/RC"); 
      ratiodedx[ic]->Draw("pE"); 

      c1->cd(5+5*ic);
      c1->cd(5+5*ic)->SetLogy(0);
      ratiodca[ic]->SetTitle(histtitle.c_str()); 
      ratiodca[ic]->GetXaxis()->SetTitle("|dca| (cm)"); 
      ratiodca[ic]->GetYaxis()->SetTitle("Data/RC"); 
      ratiodca[ic]->Draw("pE"); 
    }
    c1->SaveAs(Form("figures/KaonTTrees/KaonBasicVariables_DataRCRatio.pdf"));
    

    
    

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
    //  c1->cd(i+1)->SetGrid(0,0); 
    //}
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

}

