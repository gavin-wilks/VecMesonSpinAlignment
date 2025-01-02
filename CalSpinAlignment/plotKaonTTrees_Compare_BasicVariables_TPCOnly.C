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

void plotKaonTTrees_Compare_BasicVariables_TPCOnly(int energy = 4, bool datarcweight = false, int variable = 0, int axis=1, double low=-0.5, double high = 0.5, string kinematic = "0") 
{

    gStyle->SetOptStat(0);

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);    

    TFile *inputweightsptyphi = TFile::Open(Form("phiweights/pt_y_phi_weights_v2_FunctionEmbed.root"));

    TH3F *h_mpTyv2Weights[9];
    for(int i_cent = 2; i_cent < 6; ++i_cent)
    {
      std::string HistName = Form("h_ratio_cent%d",i_cent);
      h_mpTyv2Weights[i_cent] = (TH3F*) inputweightsptyphi->Get(HistName.c_str());
    }

     
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
    //TH1F* pt[5][3];
    //TH1F* rapidity[5][3]; 
    //TH1F* eta[5][3];
    //TH1F* phi[5][3];

    ////THnF* ptyetaphi[5][3];

    //TH1F* ratiopt[3];
    //TH1F* ratiorapidity[3]; 
    //TH1F* ratioeta[3];
    //TH1F* ratiophi[3];
    // 
    //TH2F* ptrapidity[5][3];
    //TH2F* pteta[5][3]; 
    //TH2F* ptphi[5][3];
    //TH2F* rapidityeta[5][3];
    //TH2F* rapidityphi[5][3];
    //TH2F* etaphi[5][3]; 

    //TH2F* ratioptrapidity[3];
    //TH2F* ratiopteta[3]; 
    //TH2F* ratioptphi[3];
    //TH2F* ratiorapidityeta[3];
    //TH2F* ratiorapidityphi[3];
    //TH2F* ratioetaphi[3]; 
  
    //// Tracking information
    //THnF* ptyetaphi[3][2];
    THnF* ptyetaphi[4][2];
    //THnF* ptyetaphib[3][2];
    //THnF* ptyetaphic[3][2];
    //THnF* ptyetaphid[3][2];
    //THnF* ptyetaphie[3][2];

    TH1F* var[4][2]; 
    TH2F* varpt[4][2]; 
    TH2F* varrapidity[4][2]; 
    TH2F* vareta[4][2]; 
    TH2F* varphi[4][2]; 

    TH1F* ratiovar[2]; 
    TH2F* ratiovarpt[2]; 
    TH2F* ratiovarrapidity[2]; 
    TH2F* ratiovareta[2]; 
    TH2F* ratiovarphi[2]; 

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
    // 
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


    TH3F* pty[4][2];
    TH3F* ptphi[4][2];
    TH3F* yphi[4][2];

    const int npt = 35;
    const float ptmin = 0.0, ptmax = 3.5;
    const int ny = 48;
    const float ymin = -1.2, ymax = 1.2;
    const int neta = 30;
    const float etamin = -1.5, etamax = 1.5;
    const int nphi = 48;
    const float phimin = -TMath::Pi(), phimax = TMath::Pi();

    const int nhf = 50;
    const float nhfmin = -0.5, nhfmax = 100.5; 
    const int nhm = 50;
    const float nhmmin = -0.5, nhmmax = 100.5; 
    const int nhr = 20;
    const float nhrmin = 0.5, nhrmax = 1.1; 
    const int nde = 50;
    const float demin = 0.0, demax = 10.0; 
    const int ndca = 50;
    const float dcamin = 0.0, dcamax = 2.1; 

    Int_t    bins[4] = {npt, ny, nphi, nhf};
    Double_t xmin[4] = {ptmin, ymin, phimin, nhfmin};
    Double_t xmax[4] = {ptmax, ymax, phimax, nhfmax};

    switch(variable) {
      case 0: 
        bins[3] = nhf;
        xmin[3] = nhfmin;
        xmax[3] = nhfmax;
        break;
      case 1: 
        bins[3] = nhm;
        xmin[3] = nhmmin;
        xmax[3] = nhmmax;
        break;
      case 2: 
        bins[3] = nhr;
        xmin[3] = nhrmin;
        xmax[3] = nhrmax;
        break;
      case 3: 
        bins[3] = nde;
        xmin[3] = demin;
        xmax[3] = demax;
        break;
      case 4: 
        bins[3] = ndca;
        xmin[3] = dcamin;
        xmax[3] = dcamax;
        break;
    } 

    string charge[3] = {"plus","minus","phi"};
    string mixing[5] = {"SE","ME","SM","RC","MC"};

    for (int im = 0; im < 4; im++)
    {
      for (int ic = 0; ic < 2; ic++)
      {
        string hist;
        hist = Form("pty_k%s_%s",charge[ic].c_str(),mixing[im].c_str());
        pty[im][ic]   = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, ny, ymin, ymax, bins[3], xmin[3], xmax[3]);
        ptphi[im][ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nphi, phimin, phimax, bins[3], xmin[3], xmax[3]);
        yphi[im][ic]  = new TH3F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nphi, phimin, phimax, bins[3], xmin[3], xmax[3]);
      }
    }
    cout << "Set up histograms" << endl;

    string parameter[5] = {"nhf","nhm","nhr","de","dca"};
    TFile *mKaonFileSE = TFile::Open("../data/kaon3d_SE_19GeV_TPCOnly.root");
    TFile *mKaonFileME = TFile::Open("../data/kaon3d_ME_19GeV_TPCOnly.root");
    for (int ic = 0; ic < 2; ic++)
    {
      string hist;
      string histrename;
      hist = Form("pty%s_k%s",parameter[variable].c_str(),charge[ic].c_str());
      histrename = Form("pty%s_k%s_%s",parameter[variable].c_str(),charge[ic].c_str(),mixing[im].c_str());
      pty[0][ic] = (TH3F*) ((TH3F*)mKaonFileSE->Get(hist.c_str()))->Clone(histrename.c_str());

      hist = Form("ptphi%s_k%s",parameter[variable].c_str(),charge[ic].c_str());
      histrename = Form("ptphi%s_k%s_%s",parameter[variable].c_str(),charge[ic].c_str(),mixing[im].c_str());
      ptphi[0][ic] = (TH3F*) ((TH3F*)mKaonFileSE->Get(hist.c_str()))->Clone(histrename.c_str());

      hist = Form("yphi%s_k%s",parameter[variable].c_str(),charge[ic].c_str());
      histrename = Form("yphi%s_k%s_%s",parameter[variable].c_str(),charge[ic].c_str(),mixing[im].c_str());
      yphi[0][ic] = (TH3F*) ((TH3F*)mKaonFileSE->Get(hist.c_str()))->Clone(histrename.c_str());

      hist = Form("pty%s_k%s",parameter[variable].c_str(),charge[ic].c_str());
      histrename = Form("pty%s_k%s_%s",parameter[variable].c_str(),charge[ic].c_str(),mixing[im].c_str());
      pty[1][ic] = (TH3F*) ((TH3F*)mKaonFileME->Get(hist.c_str()))->Clone(histrename.c_str());

      hist = Form("ptphi%s_k%s",parameter[variable].c_str(),charge[ic].c_str());
      histrename = Form("ptphi%s_k%s_%s",parameter[variable].c_str(),charge[ic].c_str(),mixing[im].c_str());
      ptphi[1][ic] = (TH3F*) ((TH3F*)mKaonFileME->Get(hist.c_str()))->Clone(histrename.c_str());

      hist = Form("yphi%s_k%s",parameter[variable].c_str(),charge[ic].c_str());
      histrename = Form("yphi%s_k%s_%s",parameter[variable].c_str(),charge[ic].c_str(),mixing[im].c_str());
      yphi[1][ic] = (TH3F*) ((TH3F*)mKaonFileME->Get(hist.c_str()))->Clone(histrename.c_str());
    }
    
    



    Long64_t nentriesEmbed = mKaonTreeEmbed->GetEntries();
    //nentriesEmbed = 1000000;//mKaonTreeEmbed->GetEntries();

   
    string inputfile = Form("../output/AuAu%s/Phi/TPCOnly_ptyspectra_RawPhiPtSys_eta1_eta1_PolySys.root",vmsa::mBeamEnergy[energy].c_str());
    TFile *File_InPut = TFile::Open(inputfile.c_str());
    string KEY = Form("ptyspectra_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,"2nd",0,0,vmsa::mPID[0].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
    TH2F *ptrapidityPhiData = (TH2F*) File_InPut->Get(KEY.c_str());
    TH2F *ptrapidityPhiRc = new TH2F("ptrapidityPhiRC","ptrapidityPhiRC",vmsa::rebinpttotal,vmsa::rebinptval,vmsa::rebinytotal,vmsa::rebinyval); 
 
    if(datarcweight)
    {
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

        if(mRcExists && mPassTpc /*&& mPassTof && mPassM2*/ && mPassNsig)    
        {
          TLorentzVector mLRcMeson = *mLRcKp + *mLRcKm;

          float mRcPt = mLRcMeson.Pt();      
          float mRcRapidity = mLRcMeson.Rapidity();      

          ptrapidityPhiRc->Fill(mRcPt,mRcRapidity,mWeight*PhiWeight);
        }
      }
    }

    TH2F *ptrapidityPhiDataRc;

    if(datarcweight)
    {
      double integralPhiData = ptrapidityPhiData->Integral(1,vmsa::rebinpttotal,1,vmsa::rebinytotal); 
      double integralPhiRc = ptrapidityPhiRc->Integral(1,vmsa::rebinpttotal,1,vmsa::rebinytotal); 
      double DataRcPhi = integralPhiData/integralPhiRc;
    
      ptrapidityPhiDataRc = (TH2F*) ptrapidityPhiData->Clone();
      ptrapidityPhiRc->Scale(DataRcPhi);
      ptrapidityPhiDataRc->Divide(ptrapidityPhiRc);
    }

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

      float mMcPhiPsi = mMcPhi - mEpFull;
      while(mMcPhiPsi < 0.0) mMcPhiPsi += 2.0*TMath::Pi();
      while(mMcPhiPsi > 2.0*TMath::Pi()) mMcPhiPsi -= 2.0*TMath::Pi();
     
      float PhiWeight = h_mpTyv2Weights[mCent]->GetBinContent(h_mpTyv2Weights[mCent]->FindBin(mMcPt,mMcRapidity,mMcPhiPsi));
 
      if(mRcExists && mPassTpc /*&& mPassTof && mPassM2*/ && mPassNsig)    
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

        float DataRcWeightVal = 1.0;
        if(datarcweight) DataRcWeightVal = ptrapidityPhiDataRc->GetBinContent(ptrapidityPhiDataRc->FindBin(mRcPt,mRcRapidity));

        switch(variable)  {
          case 0:
            pty[3][0]->Fill(mRcPtP,mRcRapidityP,mNHitsFitP,mWeight*PhiWeight*DataRcWeightVal);
            pty[3][1]->Fill(mRcPtM,mRcRapidityM,mNHitsFitM,mWeight*PhiWeight*DataRcWeightVal);
            ptphi[3][0]->Fill(mRcPtP,mRcPhiP,mNHitsFitP,mWeight*PhiWeight*DataRcWeightVal);
            ptphi[3][1]->Fill(mRcPtM,mRcPhiM,mNHitsFitM,mWeight*PhiWeight*DataRcWeightVal);
            yphi[3][0]->Fill(mRcRapidityP,mRcPhiP,mNHitsFitP,mWeight*PhiWeight*DataRcWeightVal);
            yphi[3][1]->Fill(mRcRapidityM,mRcPhiM,mNHitsFitM,mWeight*PhiWeight*DataRcWeightVal);
            break;
          case 1:
            pty[3][0]->Fill(mRcPtP,mRcRapidityP,mNHitsMaxP,mWeight*PhiWeight*DataRcWeightVal);
            pty[3][1]->Fill(mRcPtM,mRcRapidityM,mNHitsMaxM,mWeight*PhiWeight*DataRcWeightVal);
            ptphi[3][0]->Fill(mRcPtP,mRcPhiP,mNHitsMaxP,mWeight*PhiWeight*DataRcWeightVal);
            ptphi[3][1]->Fill(mRcPtM,mRcPhiM,mNHitsMaxM,mWeight*PhiWeight*DataRcWeightVal);
            yphi[3][0]->Fill(mRcRapidityP,mRcPhiP,mNHitsMaxP,mWeight*PhiWeight*DataRcWeightVal);
            yphi[3][1]->Fill(mRcRapidityM,mRcPhiM,mNHitsMaxM,mWeight*PhiWeight*DataRcWeightVal);
            break;
          case 2:
            pty[3][0]->Fill(mRcPtP,mRcRapidityP,float(mNHitsFitP)/float(mNHitsMaxP),mWeight*PhiWeight*DataRcWeightVal);
            pty[3][1]->Fill(mRcPtM,mRcRapidityM,float(mNHitsFitM)/float(mNHitsMaxM),mWeight*PhiWeight*DataRcWeightVal);
            ptphi[3][0]->Fill(mRcPtP,mRcPhiP,float(mNHitsFitP)/float(mNHitsMaxP),mWeight*PhiWeight*DataRcWeightVal);
            ptphi[3][1]->Fill(mRcPtM,mRcPhiM,float(mNHitsFitM)/float(mNHitsMaxM),mWeight*PhiWeight*DataRcWeightVal);
            yphi[3][0]->Fill(mRcRapidityP,mRcPhiP,float(mNHitsFitP)/float(mNHitsMaxP),mWeight*PhiWeight*DataRcWeightVal);
            yphi[3][1]->Fill(mRcRapidityM,mRcPhiM,float(mNHitsFitM)/float(mNHitsMaxM),mWeight*PhiWeight*DataRcWeightVal);
            break;
          case 3:
            pty[3][0]->Fill(mRcPtP,mRcRapidityP,mDEdxP,mWeight*PhiWeight*DataRcWeightVal);
            pty[3][1]->Fill(mRcPtM,mRcRapidityM,mDEdxM,mWeight*PhiWeight*DataRcWeightVal);
            ptphi[3][0]->Fill(mRcPtP,mRcPhiP,mDEdxP,mWeight*PhiWeight*DataRcWeightVal);
            ptphi[3][1]->Fill(mRcPtM,mRcPhiM,mDEdxM,mWeight*PhiWeight*DataRcWeightVal);
            yphi[3][0]->Fill(mRcRapidityP,mRcPhiP,mDEdxP,mWeight*PhiWeight*DataRcWeightVal);
            yphi[3][1]->Fill(mRcRapidityM,mRcPhiM,mDEdxM,mWeight*PhiWeight*DataRcWeightVal);
            break;
          case 4:
            pty[3][0]->Fill(mRcPtP,mRcRapidityP,mDcaP,mWeight*PhiWeight*DataRcWeightVal);
            pty[3][1]->Fill(mRcPtM,mRcRapidityM,mDcaM,mWeight*PhiWeight*DataRcWeightVal);
            ptphi[3][0]->Fill(mRcPtP,mRcPhiP,mDcaP,mWeight*PhiWeight*DataRcWeightVal);
            ptphi[3][1]->Fill(mRcPtM,mRcPhiM,mDcaM,mWeight*PhiWeight*DataRcWeightVal);
            yphi[3][0]->Fill(mRcRapidityP,mRcPhiP,mDcaP,mWeight*PhiWeight*DataRcWeightVal);
            yphi[3][1]->Fill(mRcRapidityM,mRcPhiM,mDcaM,mWeight*PhiWeight*DataRcWeightVal);
            break;
        }
      }
    }

 
    for(int ic = 0; ic < 2; ic++)
    {
      pty[2][ic] = (TH3F*) pty[0][ic]->Clone();
      pty[2][ic]->Add(pty[1][ic],-1.0);
      ptphi[2][ic] = (TH3F*) ptphi[0][ic]->Clone();
      ptphi[2][ic]->Add(ptphi[1][ic],-1.0);
      yphi[2][ic] = (TH3F*) yphi[0][ic]->Clone();
      yphi[2][ic]->Add(yphi[1][ic],-1.0);

      //Int_t y_min_bin;// = ptyetaphi[2][ic]->GetAxis(axis)->FindBin(low);
      //Int_t y_max_bin;// = ptyetaphi[2][ic]->GetAxis(axis)->FindBin(high)-1;

      //switch(axis)  {
      //  case 0: // pt
      //    Int_t y_min_bin = ptyetaphi[2][ic]->GetXaxis()->FindBin(low);
      //    Int_t y_max_bin = ptyetaphi[2][ic]->GetXaxis()->FindBin(high)-1;
      //    pty[2][ic]->GetXaxis()->SetRange(y_min_bin,y_max_bin);
      //    pty[3][ic]->GetXaxis()->SetRange(y_min_bin,y_max_bin);
      //    break;
      //  case 1: // y
      //    Int_t y_min_bin = ptyetaphi[2][ic]->GetYaxis()->FindBin(low);
      //    Int_t y_max_bin = ptyetaphi[2][ic]->GetYaxis()->FindBin(high)-1;
      //    pty[2][ic]->GetYaxis()->SetRange(y_min_bin,y_max_bin);
      //    pty[3][ic]->GetYaxis()->SetRange(y_min_bin,y_max_bin);
      //    break;
      //  case 2: // phi
      //
      //Int_t y_min_bin = ptyetaphi[2][ic]->GetAxis(axis)->FindBin(low);
      //Int_t y_max_bin = ptyetaphi[2][ic]->GetAxis(axis)->FindBin(high)-1;
      //pty[2][ic]->GetAxis(axis)->SetRange(y_min_bin,y_max_bin);
      //pty[3][ic]->GetAxis(axis)->SetRange(y_min_bin,y_max_bin);

      var[2][ic]         = (TH1F*) pty[2][ic]  ->ProjectionZ();
      varpt[2][ic]       = (TH2F*) pty[2][ic]  ->Project3D("zx");
      varrapidity[2][ic] = (TH2F*) pty[2][ic]  ->Project3D("zy");
      varphi[2][ic]      = (TH2F*) ptphi[2][ic]->Project3D("zy");

      var[3][ic]         = (TH1F*) pty[3][ic]  ->ProjectionZ();
      varpt[3][ic]       = (TH2F*) pty[3][ic]  ->Project3D("zx");
      varrapidity[3][ic] = (TH2F*) pty[3][ic]  ->Project3D("zx");
      varphi[3][ic]      = (TH2F*) ptphi[3][ic]->Project3D("zy");

      double integralData = var[2][ic]->Integral(0,-1);
      double integralRC   = var[3][ic]->Integral(0,-1);
      double DataRC = integralData/integralRC;

      ratiovar[ic] = (TH1F*) var[2][ic]->Clone();
      var[3][ic]->Scale(DataRC);
      ratiovar[ic]->Divide(var[3][ic]);

      ratiovarpt[ic] = (TH2F*) varpt[2][ic]->Clone();
      varpt[3][ic]->Scale(DataRC);
      ratiovarpt[ic]->Divide(varpt[3][ic]);

      ratiovarrapidity[ic] = (TH2F*) varrapidity[2][ic]->Clone();
      varrapidity[3][ic]->Scale(DataRC);
      ratiovarrapidity[ic]->Divide(varrapidity[3][ic]);

      ratiovarphi[ic] = (TH2F*) varphi[2][ic]->Clone();
      varphi[3][ic]->Scale(DataRC);
      ratiovarphi[ic]->Divide(varphi[3][ic]);

    }

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->Divide(3,2);
    for(int i = 0; i < 6; i++)
    {
      c1->cd(i+1);
      c1->cd(i+1)->SetLeftMargin(0.15);
      c1->cd(i+1)->SetRightMargin(0.15);
      c1->cd(i+1)->SetBottomMargin(0.15);
      c1->cd(i+1)->SetTicks(1,1);
      c1->cd(i+1)->SetGrid(0,0); 
    }

    string varfile[5] = {"NHitsFit","NHitsMax","NHitsRatio","dEdx","dca"};
    string varname[5] = {"NHitsFit","NHitsMax","NHitsRatio","dE/dx","dca (cm)"};

    string setting[2] = {"Same - Mixed Event", "RC from Embedding"};
    string settingfile[2] = {"SEminusME", "RC"};

    for(int im = 2; im < 4; im++)
    {
      for(int ic = 0; ic < 2; ic++)
      {
       
        string histtitle = Form("K%s, %s %s",charge[ic].c_str(),varname[variable].c_str(),setting[im-2].c_str());

        c1->cd(1);
        c1->cd(1)->SetLogy(0);
        var[im][ic]->SetTitle(histtitle.c_str()); 
        var[im][ic]->GetXaxis()->SetTitle(varname[variable].c_str()); 
        var[im][ic]->GetYaxis()->SetTitle("Count"); 
        var[im][ic]->Draw("pE"); 

        c1->cd(2);
        c1->cd(2)->SetLogz();
        varpt[im][ic]->SetTitle(histtitle.c_str()); 
        varpt[im][ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
        varpt[im][ic]->GetYaxis()->SetTitle(varname[variable].c_str()); 
        varpt[im][ic]->Draw("Colz"); 

        c1->cd(3);
        c1->cd(3)->SetLogz();
        varrapidity[im][ic]->SetTitle(histtitle.c_str()); 
        varrapidity[im][ic]->GetXaxis()->SetTitle("y"); 
        varrapidity[im][ic]->GetYaxis()->SetTitle(varname[variable].c_str()); 
        varrapidity[im][ic]->Draw("Colz"); 

        //c1->cd(4);
        //vareta[im][ic]->SetTitle(histtitle.c_str()); 
        //vareta[im][ic]->GetXaxis()->SetTitle("#eta"); 
        //vareta[im][ic]->GetYaxis()->SetTitle(varname[variable].c_str()); 
        //vareta[im][ic]->Draw("Colz"); 

        c1->cd(4);
        c1->cd(4)->SetLogz();
        varphi[im][ic]->SetTitle(histtitle.c_str()); 
        varphi[im][ic]->GetXaxis()->SetTitle("#phi"); 
        varphi[im][ic]->GetYaxis()->SetTitle(varname[variable].c_str()); 
        varphi[im][ic]->Draw("Colz"); 

        c1->SaveAs(Form("figures/KaonTTreesTPCOnly/K%s_%s_%s_%s.pdf",charge[ic].c_str(),varfile[variable].c_str(),settingfile[im-2].c_str(),kinematic.c_str()));
      }
    }
    for(int ic = 0; ic < 2; ic++)
    {
     
      string histtitle = Form("K%s, Data/RC %s",charge[ic].c_str(),varname[variable].c_str());

      TLegend *leg = new TLegend(0.4,0.75,0.6,0.9);
   
      c1->cd(1);
      c1->cd(1)->SetLogy(0);
      var[2][ic]->SetTitle(histtitle.c_str()); 
      var[2][ic]->GetXaxis()->SetTitle(varname[variable].c_str()); 
      var[2][ic]->GetYaxis()->SetTitle("Count"); 

      double max = var[2][ic]->GetMaximum();
      double min = var[2][ic]->GetMinimum();
      if(var[3][ic]->GetMaximum() > max) max = var[3][ic]->GetMaximum();
      if(var[3][ic]->GetMinimum() < min) min = var[3][ic]->GetMinimum();
      var[2][ic]->GetYaxis()->SetRangeUser(0.9*min,1.1*max); 
      var[2][ic]->SetMarkerStyle(20);
      var[2][ic]->SetMarkerColor(kBlue);
      var[2][ic]->SetLineColor(kBlue);
      var[2][ic]->Draw("pE"); 
      var[3][ic]->SetMarkerStyle(20);
      var[3][ic]->SetMarkerColor(kOrange+7);
      var[3][ic]->SetLineColor(kOrange+7);
      var[3][ic]->Draw("pE same");
      leg->AddEntry(var[2][ic],"Data","p");
      leg->AddEntry(var[3][ic],"RC","p");
      leg->Draw("same");


      c1->cd(2);
      c1->cd(2)->SetLogy(0);
      ratiovar[ic]->SetTitle(histtitle.c_str()); 
      ratiovar[ic]->GetXaxis()->SetTitle(varname[variable].c_str()); 
      ratiovar[ic]->GetYaxis()->SetTitle("Data/RC"); 
      ratiovar[ic]->Draw("pE"); 

      c1->cd(4);
      c1->cd(4)->SetLogz();
      ratiovarpt[ic]->SetTitle(histtitle.c_str()); 
      ratiovarpt[ic]->GetXaxis()->SetTitle("p_{T} GeV/c"); 
      ratiovarpt[ic]->GetYaxis()->SetTitle(varname[variable].c_str()); 
      ratiovarpt[ic]->Draw("Colz"); 

      c1->cd(5);
      c1->cd(5)->SetLogz();
      ratiovarrapidity[ic]->SetTitle(histtitle.c_str()); 
      ratiovarrapidity[ic]->GetXaxis()->SetTitle("y"); 
      ratiovarrapidity[ic]->GetYaxis()->SetTitle(varname[variable].c_str()); 
      ratiovarrapidity[ic]->Draw("Colz"); 

      //c1->cd(4);
      //ratiovareta[ic]->SetTitle(histtitle.c_str()); 
      //ratiovareta[ic]->GetXaxis()->SetTitle("#eta"); 
      //ratiovareta[ic]->GetYaxis()->SetTitle(varname[variable].c_str()); 
      //ratiovareta[ic]->Draw("Colz"); 

      c1->cd(6);
      c1->cd(6)->SetLogz();
      ratiovarphi[ic]->SetTitle(histtitle.c_str()); 
      ratiovarphi[ic]->GetXaxis()->SetTitle("#phi"); 
      ratiovarphi[ic]->GetYaxis()->SetTitle(varname[variable].c_str()); 
      ratiovarphi[ic]->Draw("Colz"); 

      c1->SaveAs(Form("figures/KaonTTreesTPCOnly/K%s_%s_DataRCRatio_%s.pdf",charge[ic].c_str(),varfile[variable].c_str(),kinematic.c_str()));
    }
}

