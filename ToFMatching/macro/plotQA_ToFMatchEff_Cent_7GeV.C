#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/type.h"
#include "../StRoot/StToFMatchMaker/StToFMatchCons.h"

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

// TH1D* hTmp;

/*
double tof_Kaon(double* x, double* par)
{
   if ((x[0] > 0.4 && x[0] < 0.9) || (x[0] > 1.1 && x[0] < 2.5))
   {
      int bin = hTmp->FindBin(x[0]);
      return hTmp->GetBinContent(bin);
   }
   else
   {
      return par[0] / (pow(x[0] - par[1], 2) + par[2]) - par[4] / (exp(x[0] - par[3]) + par[5]) + par[6];
   }
}
*/
double tof_Kaon(double* x, double* par)
{
  return par[0]*(1.0 / (pow(x[0] - par[1], 2) + par[2]) - par[4] / (exp(x[0] - par[3]) + par[5]) + par[6]);
}

void plotQA_ToFMatchEff_Cent_7GeV(const int energy = 4, const int cent = 0, const int dpid = 0, const int charge = 0, const int pid = 0, const int mKpCent = 3, const int mKmCent = 4)
{
  int const mColor[12] = {11,11,2,2,4,4,7,7,8,8,38,38};
  int const mStyle[12] = {30,24,30,24,30,24,30,24,30,24,30,24};
  int const mLineStyle[12] = {3,2,3,2,3,2,3,2,3,2,3,2};
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);

  std::string inputfile = Form("../../output/%s/Eff_%s_ToFMatch.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  cout << "Loaded file" << endl;

  TH1DMap h_mEffCent;
  TH1DMap h_mEffCent_plot;
  TH1DMap h_mEffEta;
  TH1DMap h_mEffEta_plot;
  TH1DMap h_mEfficiency;
  TH1DMap h_mEfficiency_plot;

  string HistName = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent);
  cout << HistName << endl;
  h_mEffCent[HistName] = (TH1D*)File_InPut->Get(HistName.c_str())->Clone();
  cout << "loaded first efficiency plot" << endl;
  h_mEffCent[HistName]->SetMarkerStyle(24);
  h_mEffCent[HistName]->SetMarkerSize(1.2);
  h_mEffCent[HistName]->SetMarkerColor(2);
  h_mEffCent[HistName]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mEffCent[HistName]->GetXaxis()->CenterTitle();
  h_mEffCent[HistName]->GetYaxis()->SetTitle("Efficiency");
  h_mEffCent[HistName]->GetYaxis()->SetRangeUser(0.0,1.5);
  h_mEffCent_plot[HistName] = (TH1D*)h_mEffCent[HistName]->Clone();
  for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
  {
    cout << "Eta bin = " << i_eta << endl;
    HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta);
    h_mEffEta[HistName] = (TH1D*)File_InPut->Get(HistName.c_str())->Clone();
    h_mEffEta[HistName]->SetMarkerStyle(29);
    h_mEffEta[HistName]->SetMarkerSize(1.6);
    h_mEffEta[HistName]->SetMarkerColor(2);
    h_mEffEta[HistName]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_mEffEta[HistName]->GetXaxis()->CenterTitle();
    h_mEffEta[HistName]->GetYaxis()->SetTitle("Efficiency");
    h_mEffEta[HistName]->GetYaxis()->SetRangeUser(0.0,1.5);
    h_mEffEta_plot[HistName] = (TH1D*)h_mEffEta[HistName]->Clone();
    for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi)
    {
      cout << "Eta bin = " << i_eta << "     phi = " << i_phi << endl;
      HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta,i_phi);
      h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str())->Clone();
      h_mEfficiency[HistName]->SetMarkerStyle(mStyle[i_phi]);
      h_mEfficiency[HistName]->SetMarkerSize(1.2);
      h_mEfficiency[HistName]->SetMarkerColor(mColor[i_phi]);
      h_mEfficiency[HistName]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h_mEfficiency[HistName]->GetXaxis()->CenterTitle();
      h_mEfficiency[HistName]->GetYaxis()->SetTitle("Efficiency");
      h_mEfficiency[HistName]->GetYaxis()->SetRangeUser(0.0,1.5);
      h_mEfficiency_plot[HistName] = (TH1D*)h_mEfficiency[HistName]->Clone();
    }
  }

  cout << "Loaded all plots" << endl;

  TH1DMap h_mFitPar; // initialize histogram for fit parameters
  std::string HistName_FitPar;
  HistName_FitPar = Form("h_mFitParameters_%s%s_Cent_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent);
  h_mFitPar[HistName_FitPar] = new TH1D(HistName_FitPar.c_str(),HistName_FitPar.c_str(),7,-0.5,6.5);
  for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
  {
    HistName_FitPar = Form("h_mFitParameters_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta);
    h_mFitPar[HistName_FitPar] = new TH1D(HistName_FitPar.c_str(),HistName_FitPar.c_str(),7,-0.5,6.5);
    for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi)
    {
      HistName_FitPar = Form("h_mFitParameters_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta,i_phi);
      h_mFitPar[HistName_FitPar] = new TH1D(HistName_FitPar.c_str(),HistName_FitPar.c_str(),7,-0.5,6.5);
    }
  }

  float dip_1st[2][9] = {{0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35},
                         {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}};
  float dip_2nd[2][9] = {{0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90},
                         {0.85,0.85,1.00,0.85,0.85,1.00,0.85,0.85,0.85}};
  float dip_3rd[2][9] = {{1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00},
                         {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}};
  float dip_4th[2][9] = {{2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50},
                         {2.10,2.10,2.00,2.10,2.10,2.00,2.10,2.10,2.10}};
  double range[2][9] =  {{3.80,4.20,4.20,5.00,5.00,5.00,5.00,5.00,5.00},
                         {2.80,3.00,3.40,3.85,4.10,3.30,4.20,4.50,4.50}};

  float dip_1st_eta[2][10][9] = {{{0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //0
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //1
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //2
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //3
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //4
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //5
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //6
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //7
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //8
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}},//9
                                 {{0.35,0.35,0.35,0.35,0.30,0.35,0.30,0.30,0.30}, //0
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //1
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //2
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //3
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //4
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //5
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //6
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //7
                                  {0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35}, //8
                                  {0.35,0.35,0.35,0.35,0.30,0.30,0.30,0.30,0.30}}};//9

  float dip_2nd_eta[2][10][9] = {{{0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.60,0.70}, //0
                                  {0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90}, //1
                                  {1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00}, //2
                                  {1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00}, //3
                                  {1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00}, //4
                                  {1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00}, //5
                                  {1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00}, //6
                                  {1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00}, //7
                                  {0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90}, //8
                                  {0.70,0.70,0.70,0.70,0.60,0.70,0.70,0.60,0.70}},//9
                                 {{0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70}, //0
                                  {0.85,0.85,1.00,0.85,0.85,1.00,0.85,0.85,0.85}, //1
                                  {0.85,0.85,1.00,0.85,0.85,1.00,0.85,0.85,0.85}, //2
                                  {0.85,0.85,1.00,0.85,0.85,1.00,0.85,0.85,0.85}, //3
                                  {0.85,0.85,1.00,0.85,0.85,1.00,0.85,0.85,0.85}, //4
                                  {0.85,0.85,1.00,0.85,0.85,1.00,0.85,0.85,0.85}, //5
                                  {0.85,0.85,1.00,0.85,0.85,1.00,0.85,0.85,0.85}, //6
                                  {0.85,0.85,1.00,0.85,0.85,1.00,0.85,0.85,0.85}, //7
                                  {0.85,0.85,1.00,0.85,0.85,1.00,0.85,0.85,0.85}, //8
                                  {0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70,0.70}}};//9

  float dip_3rd_eta[2][10][9] = {{{0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80}, //0
                                  {1.10,1.10,1.10,1.10,1.10,1.10,1.10,1.10,1.10}, //1
                                  {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}, //2
                                  {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}, //3
                                  {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}, //4
                                  {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}, //5
                                  {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}, //6
                                  {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}, //7
                                  {1.10,1.10,1.10,1.10,1.10,1.10,1.10,1.10,1.10}, //8
                                  {0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80}},//9
                                 {{0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80}, //0
                                  {1.10,1.10,1.10,1.10,1.10,1.10,1.10,1.10,1.10}, //1
                                  {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}, //2
                                  {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}, //3
                                  {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}, //4
                                  {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}, //5
                                  {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}, //6
                                  {1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20,1.20}, //7
                                  {1.10,1.10,1.10,1.10,1.10,1.10,1.10,1.10,1.10}, //8
                                  {0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80}}};//9
 
  float dip_4th_eta[2][10][9] = {{{1.70,1.70,2.30,1.70,2.40,2.10,2.50,2.50,2.30}, //0
                                  {1.80,1.80,2.20,1.80,1.95,2.10,2.10,2.10,2.10}, //1
                                  {2.20,2.20,2.20,2.20,2.25,2.30,2.30,2.30,2.30}, //2
                                  {2.20,2.20,2.20,2.20,2.25,2.30,2.30,2.30,2.30}, //3
                                  {2.20,2.20,2.20,2.20,2.25,2.30,2.30,2.30,2.30}, //4
                                  {2.20,2.20,2.20,2.20,2.25,2.30,2.30,2.30,2.30}, //5
                                  {2.20,2.20,2.20,2.20,2.25,2.30,2.30,2.30,2.30}, //6
                                  {2.20,2.20,2.20,2.20,2.25,2.30,2.30,2.30,2.30}, //7
                                  {1.80,1.80,1.80,1.80,1.80,2.10,2.10,2.10,2.10}, //8
                                  {1.70,1.90,1.70,1.70,2.50,2.30,2.20,2.50,2.30}},//9
                                 {{1.40,1.50,1.50,1.50,2.30,2.00,2.10,2.10,2.10}, //0
                                  {1.40,1.80,2.00,2.10,2.30,2.00,2.20,2.10,2.10}, //1
                                  {1.40,1.80,2.00,2.10,2.10,2.00,2.20,2.10,2.10}, //2
                                  {1.40,1.80,2.00,2.10,2.10,2.00,2.20,2.10,2.10}, //3
                                  {1.40,1.80,2.00,2.10,2.10,2.00,2.20,2.10,2.10}, //4
                                  {1.40,1.80,2.00,2.10,2.10,2.00,2.20,2.10,2.10}, //5
                                  {1.40,1.80,2.00,2.10,2.10,2.00,2.20,2.10,2.10}, //6
                                  {1.40,1.80,2.00,2.10,2.10,2.00,2.20,2.10,2.10}, //7
                                  {1.40,1.80,2.00,2.10,2.30,2.00,2.20,2.10,2.10}, //8
                                  {1.40,1.50,1.50,1.50,2.30,2.00,2.10,2.10,2.10}}};//9

  float range_eta[2][10][9] =   {{{3.00,3.00,2.90,3.00,3.40,3.00,3.20,3.70,2.70}, //0
                                  {3.10,3.10,3.20,3.10,3.10,3.10,3.10,3.10,3.10}, //1
                                  {3.10,3.10,3.10,3.10,3.10,3.10,3.10,3.10,3.10}, //2
                                  {3.10,3.10,3.10,3.10,3.10,3.10,3.10,3.10,3.10}, //3
                                  {3.10,3.10,3.10,3.10,3.10,3.10,3.10,3.10,3.10}, //4
                                  {3.10,3.10,3.10,3.10,3.10,3.10,3.10,3.10,3.10}, //5
                                  {3.10,3.10,3.10,3.10,3.10,3.10,3.10,3.10,3.10}, //6
                                  {3.10,3.10,3.10,3.10,3.10,3.10,3.10,3.10,3.10}, //7
                                  {3.10,3.30,3.10,3.10,3.10,3.10,3.10,3.10,3.10}, //8
                                  {3.00,3.30,3.00,3.00,3.70,3.20,3.30,3.85,2.70}},//9
                                 {{2.50,2.20,2.50,3.00,3.20,3.20,2.80,3.20,3.20}, //0
                                  {1.90,2.60,2.70,3.85,4.10,3.30,3.10,4.50,4.50}, //1
                                  {1.90,2.60,2.70,3.85,4.10,3.30,3.10,4.50,4.50}, //2
                                  {1.90,2.60,2.70,3.85,4.10,3.30,3.10,4.50,4.50}, //3
                                  {1.90,2.60,2.70,3.85,4.10,3.30,3.10,4.50,4.50}, //4
                                  {1.90,2.60,2.70,3.85,4.10,3.30,3.10,4.50,4.50}, //5
                                  {1.90,2.60,2.70,3.85,4.10,3.30,3.10,4.50,4.50}, //6
                                  {1.90,2.60,2.70,3.85,4.10,3.30,3.10,4.50,4.50}, //7
                                  {1.90,2.60,2.70,3.85,4.10,3.30,3.10,4.50,4.50}, //8
                                  {2.50,2.40,2.50,3.00,3.20,3.20,2.80,3.20,3.20}}};//9

  ////double range[9] = {4.15,4.25,4.5,5.0,5.0,5.0,5.0,5.0,5.0}; // Kplus
  //if(charge == 0 && cent == 0) range[energy] = range[energy] - 0.75;
  //if(charge == 1 && cent >= 0 && cent <= 3) range[energy] = range[energy] - 0.75;
  //if(charge == 1 && cent >= 7) range[energy] = range[energy] + 0.25;
  //if(charge == 1 && cent == 6) range[energy] = range[energy] + 0.25;
  //if(charge == 1 && cent == 0) range[energy] = 2.75;
  //if(charge == 1 && cent == 1) range[energy] = 3.5;
  // double range[7] = {4.0,3.5,4.5,4.0,4.7,4.2,5.0}; // Kminus

  int cent_kaon = -1;
  if(charge == 0) cent_kaon = mKpCent;
  if(charge == 1) cent_kaon = mKmCent;
  //if(charge == 0) cent_kaon = cent;
  //if(charge == 1) cent_kaon = cent;
  double par0[9] = {-0.018979, -0.0322611, -0.0680754, -0.0698575, -0.0315267, -0.00589929, -0.00226724, -0.00212137, -0.00389514 };
  double par1[9] = {0.0308943, -0.0939411, -0.14377, -0.19003, -0.116323, 0.180593, 0.207874, 0.208863, 0.194876};
  double par2[9] = {-0.00592033, -0.0600635, -0.0515391, -0.0708703, -0.0756912, 0.00912449, 0.00500487, 0.00497987, 0.00824164};
  double par3[9] = {1.28883, 1.53952, 1.52213, 1.01707, 1.5415, 2.75657, -0.326349, 1.11207, 1.37717};
  double par4[9] = {1.58923e-06, -0.00130657, -0.00973403, -0.0163526, -0.00162583, -2.20034e-05, 0.773984, 0.119933, -1.14531e-06};
  double par5[9] = {-0.340008, -0.261115, -0.236246, -0.345869, -0.260416, -0.0777638, 6.00519, 1.05048, -0.311133};
  double par6[9] = {0.685642, 0.687831, 0.682733, 0.674683, 0.659559, 0.639223, 0.658294, 0.626411, 0.578902};

  HistName = Form("h_mEfficiency_K%s_Cent_%d",tof::mCharge[charge].c_str(),cent); // ToF efficiency at cent 9
  for(int i_bin = 0; i_bin < h_mEffCent[HistName]->GetNbinsX(); ++i_bin) 
  {
    float pt = h_mEffCent[HistName]->GetBinCenter(i_bin+1); 
    if ((pt > dip_1st[charge][cent] && pt < dip_2nd[charge][cent]) || (pt > dip_3rd[charge][cent] && pt < dip_4th[charge][cent])) h_mEffCent[HistName]->SetBinError(i_bin+1,0.0); 
    if ((pt > dip_1st[charge][cent] && pt < dip_2nd[charge][cent]) || (pt > dip_3rd[charge][cent] && pt < dip_4th[charge][cent])) h_mEffCent[HistName]->SetBinContent(i_bin+1,0.0); 
  }

  TF1 *f_kaon_cent = new TF1("f_kaon_cent",tof_Kaon,0.1,10,7);
  f_kaon_cent->SetParameter(0,par0[cent_kaon]);
  f_kaon_cent->SetParameter(1,par1[cent_kaon]);
  f_kaon_cent->SetParameter(2,par2[cent_kaon]);
  f_kaon_cent->SetParameter(3,par3[cent_kaon]);
  f_kaon_cent->SetParameter(4,par4[cent_kaon]);
  f_kaon_cent->SetParameter(5,par5[cent_kaon]);
  f_kaon_cent->SetParameter(6,par6[cent_kaon]);
  for(int i_par = 0; i_par < 7; ++i_par)
  {
    //f_kaon_cent->SetParLimits(i_par,-10.0,10.0);
  }
  f_kaon_cent->SetRange(0.1,range[charge][cent]);

  h_mEffCent[HistName]->Fit(f_kaon_cent, "MWNR");

  double par_fit[7];
  HistName_FitPar = Form("h_mFitParameters_%s%s_Cent_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent);
  for(int i_par = 0; i_par < 7; ++i_par)
  {
    par_fit[i_par] = f_kaon_cent->GetParameter(i_par);
    h_mFitPar[HistName_FitPar]->SetBinContent(i_par+1,f_kaon_cent->GetParameter(i_par));
  }

  // fit eta distribution
  TF1 *f_kaon_eta[tof::BinEta];
  TF1 *f_kaon_phi[tof::BinEta][tof::BinPhi];
  for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta) // ToF efficiency at different eta and phi bin at cent 9
  {
    HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta);
    for(int i_bin = 0; i_bin < h_mEffEta[HistName]->GetNbinsX(); ++i_bin) 
    {
      float pt = h_mEffEta[HistName]->GetBinCenter(i_bin+1);
      //if((i_eta == 0 || i_eta == 9) || (cent == 0))
      ////if(i_eta == 0 || i_eta == 9)
      //{
      //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.15) || (pt > dip_3rd[energy]-0.15 && pt < dip_4th[energy]-0.5)) h_mEffEta[HistName]->SetBinError(i_bin+1,0.0); 
      //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.15) || (pt > dip_3rd[energy]-0.15 && pt < dip_4th[energy]-0.5)) h_mEffEta[HistName]->SetBinContent(i_bin+1,0.0); 
      //} 
      ////if(cent >= 0 && cent <= 3 && charge == 1)
      //////if(i_eta == 0 || i_eta == 9)
      ////{
      ////  if ((pt > dip_1st[energy]+0.05 && pt < dip_2nd[energy]+0.1) || (pt > dip_3rd[energy]+0.1 && pt < dip_4th[energy]-0.1)) h_mEffEta[HistName]->SetBinError(i_bin+1,0.0); 
      ////  if ((pt > dip_1st[energy]+0.05 && pt < dip_2nd[energy]+0.1) || (pt > dip_3rd[energy]+0.1 && pt < dip_4th[energy]-0.1)) h_mEffEta[HistName]->SetBinContent(i_bin+1,0.0); 
      ////}
      //if(cent == 6 && charge == 1 && i_eta == 2)
      ////if(i_eta == 0 || i_eta == 9)
      //{
      //  if ((pt > dip_1st[energy]+0.05 && pt < dip_2nd[energy]-0.05) || (pt > dip_3rd[energy]-0.1 && pt < dip_4th[energy]-0.25)) h_mEffEta[HistName]->SetBinError(i_bin+1,0.0); 
      //  if ((pt > dip_1st[energy]+0.05 && pt < dip_2nd[energy]-0.05) || (pt > dip_3rd[energy]-0.1 && pt < dip_4th[energy]-0.25)) h_mEffEta[HistName]->SetBinContent(i_bin+1,0.0); 
      //}
      //if(cent == 6 && charge == 1 && i_eta >= 3)
      ////if(i_eta == 0 || i_eta == 9)
      //{
      //  if ((pt > dip_1st[energy]+0.05 && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.25)) h_mEffEta[HistName]->SetBinError(i_bin+1,0.0); 
      //  if ((pt > dip_1st[energy]+0.05 && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.25)) h_mEffEta[HistName]->SetBinContent(i_bin+1,0.0); 
      //}
      //if(cent == 0 && charge == 1 && (i_eta <= 1 || i_eta >= 8))
      ////if(i_eta == 0 || i_eta == 9)
      //{
      //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.1) || (pt > dip_3rd[energy]-0.15 && pt < dip_4th[energy]-0.75)) h_mEffEta[HistName]->SetBinError(i_bin+1,0.0); 
      //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.1) || (pt > dip_3rd[energy]-0.15 && pt < dip_4th[energy]-0.75)) h_mEffEta[HistName]->SetBinContent(i_bin+1,0.0); 
      //}
       
      //else if(i_eta == 2 || i_eta == 12)
      //{ 
      //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.2) || (pt > dip_3rd[energy]-0.2 && pt < dip_4th[energy]-0.5)) h_mEffEta[HistName]->SetBinError(i_bin+1,0.0); 
      //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.2) || (pt > dip_3rd[energy]-0.2 && pt < dip_4th[energy]-0.5)) h_mEffEta[HistName]->SetBinContent(i_bin+1,0.0); 
      //}
      //else if(i_eta == 3 || i_eta == 11)
      //{
      //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.1) || (pt > dip_3rd[energy]-0.1 && pt < dip_4th[energy]-0.2)) h_mEffEta[HistName]->SetBinError(i_bin+1,0.0); 
      //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.1) || (pt > dip_3rd[energy]-0.1 && pt < dip_4th[energy]-0.2)) h_mEffEta[HistName]->SetBinContent(i_bin+1,0.0); 
      //
      //}
      //else
      //{
      //  if ((pt > dip_1st[cent] && pt < dip_2nd[cent]) || (pt > dip_3rd[cent] && pt < dip_4th[cent])) h_mEffEta[HistName]->SetBinError(i_bin+1,0.0); 
      //  if ((pt > dip_1st[cent] && pt < dip_2nd[cent]) || (pt > dip_3rd[cent] && pt < dip_4th[cent])) h_mEffEta[HistName]->SetBinContent(i_bin+1,0.0);
      //}
      {
        if ((pt > dip_1st_eta[charge][i_eta][cent] && pt < dip_2nd_eta[charge][i_eta][cent]) || (pt > dip_3rd_eta[charge][i_eta][cent] && pt < dip_4th_eta[charge][i_eta][cent])) h_mEffEta[HistName]->SetBinError(i_bin+1,0.0); 
        if ((pt > dip_1st_eta[charge][i_eta][cent] && pt < dip_2nd_eta[charge][i_eta][cent]) || (pt > dip_3rd_eta[charge][i_eta][cent] && pt < dip_4th_eta[charge][i_eta][cent])) h_mEffEta[HistName]->SetBinContent(i_bin+1,0.0);
      }
    }
    double par_fite[7] = {0.0};
    std::string FuncName = Form("f_kaon_cent_%d_eta_%d",cent,i_eta);
    f_kaon_eta[i_eta] = new TF1(FuncName.c_str(),tof_Kaon,0.2,10,7);
    f_kaon_eta[i_eta]->SetParameter(0,par_fit[0]);
    //if(i_eta == 1 || i_eta == 13) f_kaon_eta[i_eta]->SetParameter(0,0.003);
    //if(i_eta == 0 || i_eta == 14) f_kaon_eta[i_eta]->FixParameter(0,0.0);
    for(int i_par = 1; i_par < 7; ++i_par)
    {
      f_kaon_eta[i_eta]->SetParameter(i_par,par_fit[i_par]);
    }
    f_kaon_eta[i_eta]->SetRange(0.1,range_eta[charge][i_eta][cent]);
    //if(i_eta < 4 || i_eta > 10) f_kaon_eta[i_eta]->SetRange(0.2,range[energy]);

    h_mEffEta[HistName]->Fit(f_kaon_eta[i_eta], "NR");

    HistName_FitPar = Form("h_mFitParameters_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta);
    for(int i_bin = 0; i_bin < 7; ++i_bin)
    {
      par_fite[i_bin] = f_kaon_eta[i_eta]->GetParameter(i_bin);
      h_mFitPar[HistName_FitPar]->SetBinContent(i_bin+1,f_kaon_eta[i_eta]->GetParameter(i_bin));
    }

    for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi) // eta & phi dependence
    {
      HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta,i_phi);
      for(int i_bin = 0; i_bin < h_mEfficiency[HistName]->GetNbinsX(); ++i_bin) 
      {
	float pt = h_mEfficiency[HistName]->GetBinCenter(i_bin+1); 
	//if(i_phi == 3 || i_phi == 5 || i_phi == 6) // 200 GeV second
	// if(i_phi == 4 || i_phi == 5 || i_phi == 6) // 62 GeV second
	// if(i_phi == 5 || i_phi == 6) // 39 GeV second
	// if(i_phi == 3 || i_phi == 4 || i_phi == 5 || i_phi == 6) // 27 GeV second
	// if(i_phi == 3 || i_phi == 5 || i_phi == 6) // 19 GeV second
	// if(i_phi == 5 || i_phi == 6) // 11 GeV second
	//{
	//  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
	//  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0); 
	//}
	//else
	//{
        //if((i_eta == 0 || i_eta == 9) || (cent == 0))
        //{
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.15) || (pt > dip_3rd[energy]-0.15 && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.15) || (pt > dip_3rd[energy]-0.15 && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0);  
        //}
        //if(cent == 0 && charge == 1)
        //{
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0);  
        //}
        //if(cent == 6 && charge == 1 && i_eta == 2)
        //{
        //  if ((pt > dip_1st[energy]+0.05 && pt < dip_2nd[energy]-0.05) || (pt > dip_3rd[energy]-0.1 && pt < dip_4th[energy]-0.25)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
        //  if ((pt > dip_1st[energy]+0.05 && pt < dip_2nd[energy]-0.05) || (pt > dip_3rd[energy]-0.1 && pt < dip_4th[energy]-0.25)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0);  
        //}
        //if(cent == 6 && charge == 1 && i_eta >= 3)
        //{
        //  if ((pt > dip_1st[energy]+0.05 && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.25)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
        //  if ((pt > dip_1st[energy]+0.05 && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.25)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0);  
        //}
        //if(cent == 4 && charge == 1 && i_eta >= 1 && i_eta <= 8)
        //{
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0);  
        //}
        //if(cent == 3 && charge == 1)
        //{
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0);  
        //}
        //if(cent == 2 && charge == 1 && i_eta >= 1 && i_eta <= 8)
        //{
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0);  
        //}
        //if(cent == 1 && charge == 1)
        //{
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0);  
        //}
        //if(cent == 0 && charge == 1)
        //{
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.75)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]-0.75)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0);  
        //}
        //else if(i_eta == 2 || i_eta == 12)
        //{ 
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.2) || (pt > dip_3rd[energy]-0.2 && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.2) || (pt > dip_3rd[energy]-0.2 && pt < dip_4th[energy]-0.5)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0); 
        //}
        //else if(i_eta == 3 || i_eta == 11)
        //{
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.1) || (pt > dip_3rd[energy]-0.1 && pt < dip_4th[energy]-0.2)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
        //  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]-0.1) || (pt > dip_3rd[energy]-0.1 && pt < dip_4th[energy]-0.2)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0); 
        //
        //}
        //else
        //{
        //  if ((pt > dip_1st[cent] && pt < dip_2nd[cent]) || (pt > dip_3rd[cent] && pt < dip_4th[cent])) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
        //  if ((pt > dip_1st[cent] && pt < dip_2nd[cent]) || (pt > dip_3rd[cent] && pt < dip_4th[cent])) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0);
        //}
        {
          if ((pt > dip_1st_eta[charge][i_eta][cent] && pt < dip_2nd_eta[charge][i_eta][cent]) || (pt > dip_3rd_eta[charge][i_eta][cent] && pt < dip_4th_eta[charge][i_eta][cent])) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
          if ((pt > dip_1st_eta[charge][i_eta][cent] && pt < dip_2nd_eta[charge][i_eta][cent]) || (pt > dip_3rd_eta[charge][i_eta][cent] && pt < dip_4th_eta[charge][i_eta][cent])) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0);
        }


	//  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]+0.5)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
	//  if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]+0.5)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0); 
	//}
      }

      std::string FuncName = Form("f_kaon_cent_%d_eta_%d_phi_%d",cent,i_eta,i_phi);
      f_kaon_phi[i_eta][i_phi] = new TF1(FuncName.c_str(),tof_Kaon,0.1,10,7);
      f_kaon_phi[i_eta][i_phi]->SetParameter(0,par_fit[0]);
      //if(i_eta == 1 || i_eta == 13) f_kaon_phi[i_eta][i_phi]->SetParameter(0,0.003);
      //if(i_eta == 0 || i_eta == 14) f_kaon_phi[i_eta][i_phi]->FixParameter(0,0.0);
      for(int i_par = 1; i_par < 7; ++i_par)
      {
	f_kaon_phi[i_eta][i_phi]->FixParameter(i_par,par_fite[i_par]);
      }
      //f_kaon_phi[i_eta][i_phi]->SetRange(0.2,range[energy]);
      f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th_eta[charge][i_eta][cent],range_eta[charge][i_eta][cent]);
      //if(i_eta == 0 || i_eta == 9 ) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.5,range[energy]);
      //if(cent == 0) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.5,range[energy]);
      //if(cent >= 0 && cent <= 3 && charge == 1) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.1,range[energy]);
      //if(cent == 6 && charge == 1 && i_eta == 2) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.25,range[energy]);
      //if(cent == 6 && charge == 1 && i_eta >= 3) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.25,range[energy]);
      //if(cent == 6 && charge == 1 && i_eta == 1) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.25,range[energy]-0.5);
      //if(cent == 4 && charge == 1 && i_eta >= 1 && i_eta <= 8) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.5,range[energy]-0.75);
      //if(cent == 0 && charge == 1) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.75,range[energy]);
      //if(cent == 3 && charge == 1) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.5,range[energy]-0.75);
      //if(cent == 2 && charge == 1 && i_eta >= 1 && i_eta <= 8) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.5,range[energy]-0.75);
      //if(cent == 1 && charge == 1) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.5,range[energy]);

      //if(i_eta == 2 || i_eta == 12) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.5,range[energy]);
      //if(i_eta == 3 || i_eta == 11) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]-0.2,range[energy]);
      //if(i_phi == 3 || i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(0.2,dip_3rd[energy]); // 200 GeV first
      // if(i_phi == 3 || i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy],range[energy]); // 200 GeV second
      // if(i_phi == 4 || i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(0.2,dip_3rd[energy]); // 62 GeV first
      // if(i_phi == 4 || i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy],range[energy]); // 62 GeV second
      // if(i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(0.2,dip_3rd[energy]); // 39 GeV first
      // if(i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy],range[energy]); // 39 GeV second
      // if(i_phi == 3 || i_phi == 4 || i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(0.2,dip_3rd[energy]); // 27 GeV first
      // if(i_phi == 3 || i_phi == 4 || i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy],range[energy]); // 27 GeV second
      // if(i_phi == 3 || i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(0.2,dip_3rd[energy]); // 19 GeV first
      // if(i_phi == 3 || i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy],range[energy]); // 19 GeV second
      // if(i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(0.2,dip_3rd[energy]); // 11 GeV first
      // if(i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy],range[energy]); // 11 GeV second

      h_mEfficiency[HistName]->Fit(f_kaon_phi[i_eta][i_phi], "NR");

      HistName_FitPar = Form("h_mFitParameters_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta,i_phi);
      for(int i_bin = 0; i_bin < 7; ++i_bin)
      {
	h_mFitPar[HistName_FitPar]->SetBinContent(i_bin+1,f_kaon_phi[i_eta][i_phi]->GetParameter(i_bin));
      }
    }
  }

#if _PlotQA_

  TCanvas *c_Efficiency = new TCanvas("c_Efficiency","c_Efficiency",10,10,800,800);
  c_Efficiency->cd()->SetLeftMargin(0.1);
  c_Efficiency->cd()->SetBottomMargin(0.1);
  c_Efficiency->cd()->SetGrid(0,0);
  c_Efficiency->cd()->SetTicks(1,1);

  string output_start; 
  if(charge == 0) output_start = Form("figures/Efficiency/ToF/QAs/%s/EffFit_AuAu%s_%s%s_QA_cent%d_KpDefault%d.pdf[",vmsa::mPID[pid].c_str(),tof::mPID_ToF[dpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str(),cent,mKpCent);
  if(charge == 1) output_start = Form("figures/Efficiency/ToF/QAs/%s/EffFit_AuAu%s_%s%s_QA_cent%d_KmDefault%d.pdf[",vmsa::mPID[pid].c_str(),tof::mPID_ToF[dpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str(),cent,mKmCent);
  c_Efficiency->Print(output_start.c_str());

  string outputname;
  if(charge == 0) outputname = Form("figures/Efficiency/ToF/QAs/%s/EffFit_AuAu%s_%s%s_QA_cent%d_KpDefault%d.pdf",vmsa::mPID[pid].c_str(),tof::mPID_ToF[dpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str(),cent,mKpCent);
  if(charge == 1) outputname = Form("figures/Efficiency/ToF/QAs/%s/EffFit_AuAu%s_%s%s_QA_cent%d_KmDefault%d.pdf",vmsa::mPID[pid].c_str(),tof::mPID_ToF[dpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str(),cent,mKmCent);

  HistName = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent); // ToF efficiency at cent 9
  h_mEffCent[HistName]->Draw("pE");
  h_mEffCent_plot[HistName]->Draw("pE same");
  f_kaon_cent->SetLineColor(4);
  f_kaon_cent->SetLineWidth(2);
  f_kaon_cent->SetLineStyle(2);
  f_kaon_cent->SetNpx(1000);
  f_kaon_cent->SetRange(0.2,8.0);
  f_kaon_cent->Draw("l same");
  c_Efficiency->Update();
  c_Efficiency->Print(outputname.c_str());

  for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta) // ToF efficiency at different eta and phi bin at cent 9
  {
    TLegend *leg = new TLegend(0.5,0.7,0.8,0.9);
    leg->SetFillColor(10);
    leg->SetBorderSize(0);
    HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta);
    // h_mEffEta[HistName]->Draw("pE");
    //if(i_eta == 1 || i_eta == 13) h_mEffEta_plot[HistName]->GetYaxis()->SetRangeUser(0.0,0.03);
    h_mEffEta_plot[HistName]->Draw("pE");
    f_kaon_eta[i_eta]->SetLineColor(2);
    f_kaon_eta[i_eta]->SetLineWidth(4);
    f_kaon_eta[i_eta]->SetLineStyle(2);
    f_kaon_eta[i_eta]->SetNpx(1000);
    f_kaon_eta[i_eta]->SetRange(0.2,8.0);
    f_kaon_eta[i_eta]->Draw("l same");
    leg->AddEntry(h_mEffEta[HistName],"eta & phi integrated","p");
    for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi) 
    {
      HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta,i_phi);
      // h_mEfficiency[HistName]->Draw("pE same");
      h_mEfficiency_plot[HistName]->Draw("pE same");
      f_kaon_phi[i_eta][i_phi]->SetLineColor(mColor[i_phi]);
      f_kaon_phi[i_eta][i_phi]->SetLineWidth(2);
      f_kaon_phi[i_eta][i_phi]->SetLineStyle(mLineStyle[i_phi]);
      f_kaon_phi[i_eta][i_phi]->SetNpx(1000);
      // f_kaon_phi[i_eta][i_phi]->SetRange(0.2,range[energy]);
      f_kaon_phi[i_eta][i_phi]->SetRange(0.2,8.0);
      f_kaon_phi[i_eta][i_phi]->Draw("l same");
      string Leg_phi = Form("%d bin",i_phi);
      leg->AddEntry(h_mEfficiency[HistName],Leg_phi.c_str(),"p");
    }
    HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta);
    // h_mEffEta[HistName]->Draw("pE same");
    h_mEffEta_plot[HistName]->Draw("pE same");

    leg->Draw("same");
    c_Efficiency->Update();
    c_Efficiency->Print(outputname.c_str());
  }
  string output_stop; 
  if(charge == 0) output_stop = Form("figures/Efficiency/ToF/QAs/%s/EffFit_AuAu%s_%s%s_QA_cent%d_KpDefault%d.pdf]",vmsa::mPID[pid].c_str(),tof::mPID_ToF[dpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str(),cent,mKpCent);
  if(charge == 1) output_stop = Form("figures/Efficiency/ToF/QAs/%s/EffFit_AuAu%s_%s%s_QA_cent%d_KmDefault%d.pdf]",vmsa::mPID[pid].c_str(),tof::mPID_ToF[dpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str(),cent,mKmCent);
  c_Efficiency->Print(output_stop.c_str());


  TCanvas *c_Efficiency2 = new TCanvas("c_Efficiency2","c_Efficiency2",10,10,1600,1200);
  c_Efficiency2->Divide(4,3);
  for(int i = 0; i < tof::BinPhi; i++)
  {
    c_Efficiency2->cd(i+1)->SetLeftMargin(0.1);
    c_Efficiency2->cd(i+1)->SetBottomMargin(0.1);
    c_Efficiency2->cd(i+1)->SetGrid(0,0);
    c_Efficiency2->cd(i+1)->SetTicks(1,1);
  }

  TCanvas *c_Efficiency2s = new TCanvas("c_Efficiency2s","c_Efficiency2s",10,10,800,800);
  c_Efficiency2s->cd()->SetLeftMargin(0.1);
  c_Efficiency2s->cd()->SetBottomMargin(0.1);
  c_Efficiency2s->cd()->SetGrid(0,0);
  c_Efficiency2s->cd()->SetTicks(1,1);
  

  output_start; 
  if(charge == 0) output_start = Form("figures/Efficiency/ToF/QAs/%s/EffFit_AuAu%s_%s%s_QA_cent%d_IndividualPlots.pdf[",vmsa::mPID[pid].c_str(),tof::mPID_ToF[dpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str(),cent);
  if(charge == 1) output_start = Form("figures/Efficiency/ToF/QAs/%s/EffFit_AuAu%s_%s%s_QA_cent%d_IndividualPlots.pdf[",vmsa::mPID[pid].c_str(),tof::mPID_ToF[dpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str(),cent);
  c_Efficiency->Print(output_start.c_str());

  outputname;
  if(charge == 0) outputname = Form("figures/Efficiency/ToF/QAs/%s/EffFit_AuAu%s_%s%s_QA_cent%d_IndividualPlots.pdf",vmsa::mPID[pid].c_str(),tof::mPID_ToF[dpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str(),cent);
  if(charge == 1) outputname = Form("figures/Efficiency/ToF/QAs/%s/EffFit_AuAu%s_%s%s_QA_cent%d_IndividualPlots.pdf",vmsa::mPID[pid].c_str(),tof::mPID_ToF[dpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str(),cent);

  HistName = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent); // ToF efficiency at cent 9

  for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta) // ToF efficiency at different eta and phi bin at cent 9
  {
    TLegend *leg = new TLegend(0.5,0.7,0.8,0.9);
    leg->SetFillColor(10);
    leg->SetBorderSize(0);
    HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta);
    // h_mEffEta[HistName]->Draw("pE");
    //if(i_eta == 1 || i_eta == 13) h_mEffEta_plot[HistName]->GetYaxis()->SetRangeUser(0.0,0.03);
    c_Efficiency2s->cd();
    h_mEffEta_plot[HistName]->Draw("pE");
    f_kaon_eta[i_eta]->SetLineColor(2);
    f_kaon_eta[i_eta]->SetLineWidth(4);
    f_kaon_eta[i_eta]->SetLineStyle(2);
    f_kaon_eta[i_eta]->SetNpx(1000);
    f_kaon_eta[i_eta]->SetRange(0.3,8.0);
    f_kaon_eta[i_eta]->Draw("l same");
    //leg->AddEntry(h_mEffEta[HistName],"eta & phi integrated","p");
    c_Efficiency2s->Update();
    c_Efficiency2s->Print(outputname.c_str());
    for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi) 
    {
      c_Efficiency2->cd(i_phi+1);
      HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta,i_phi);
      // h_mEfficiency[HistName]->Draw("pE same");
      h_mEfficiency_plot[HistName]->Draw("pE");
      f_kaon_phi[i_eta][i_phi]->SetLineColor(mColor[i_phi]);
      f_kaon_phi[i_eta][i_phi]->SetLineWidth(2);
      f_kaon_phi[i_eta][i_phi]->SetLineStyle(mLineStyle[i_phi]);
      f_kaon_phi[i_eta][i_phi]->SetNpx(1000);
      // f_kaon_phi[i_eta][i_phi]->SetRange(0.2,range[energy]);
      f_kaon_phi[i_eta][i_phi]->SetRange(0.3,8.0);
      f_kaon_phi[i_eta][i_phi]->Draw("l same");
      //string Leg_phi = Form("%d bin",i_phi);
      //leg->AddEntry(h_mEfficiency[HistName],Leg_phi.c_str(),"p");
    }
    c_Efficiency2->Update();
    c_Efficiency2->Print(outputname.c_str());
    //HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta);
    // h_mEffEta[HistName]->Draw("pE same");
    //h_mEffEta_plot[HistName]->Draw("pE same");

    //leg->Draw("same");
    //c_Efficiency->Update();
    //c_Efficiency->Print(outputname.c_str());
  }
  output_stop; 
  if(charge == 0) output_stop = Form("figures/Efficiency/ToF/QAs/%s/EffFit_AuAu%s_%s%s_QA_cent%d_IndividualPlots.pdf]",vmsa::mPID[pid].c_str(),tof::mPID_ToF[dpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str(),cent);
  if(charge == 1) output_stop = Form("figures/Efficiency/ToF/QAs/%s/EffFit_AuAu%s_%s%s_QA_cent%d_IndividualPlots.pdf]",vmsa::mPID[pid].c_str(),tof::mPID_ToF[dpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str(),cent);
  c_Efficiency2->Print(output_stop.c_str());

#endif

  string outputfile = Form("../../output/%s/FitPar_AuAu%s_%s%s_cent%d.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent);
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  HistName_FitPar = Form("h_mFitParameters_%s%s_Cent_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent);
  h_mFitPar[HistName_FitPar]->Write();
  for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
  {
    HistName_FitPar = Form("h_mFitParameters_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta);
    h_mFitPar[HistName_FitPar]->Write();
    for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi) 
    {
      HistName_FitPar = Form("h_mFitParameters_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[dpid].c_str(),tof::mCharge[charge].c_str(),cent,i_eta,i_phi);
      h_mFitPar[HistName_FitPar]->Write();
    }
  }
}
