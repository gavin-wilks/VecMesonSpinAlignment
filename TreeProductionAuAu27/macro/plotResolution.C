#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "../StRoot/Utility/StSpinAlignmentCons.h"

using namespace std;

double ResolutionFull(double *x_val, double *par)
{
  double y;
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

void plotResolution(int beamEnergy = 6)
{
  //string inputfile = Form("./file_%s_Resolution.root",vmsa::mBeamEnergy[beamEnergy].c_str());
  //string inputfile = Form("../UtilityFilesEta1p5OfficialCent/file_%s_Resolution.root",vmsa::mBeamEnergy[beamEnergy].c_str());
  //string inputfileXu = Form("/star/u/gwilks3/Workspace/FileTransfers/file_%s_Resolution.root",vmsa::mBeamEnergy[beamEnergy].c_str());
  TGraphAsymmErrors *g_mTpcSubRes2[7];
  TFile *File_InPut[7];
  for(int i  = 6; i  < 7; i++)
  {
    string inputfile = Form("../StRoot/Utility/Resolution/file_%s_Resolution.root",vmsa::mBeamEnergy[i].c_str());
    File_InPut[i] = TFile::Open(inputfile.c_str());
    g_mTpcSubRes2[i] = new TGraphAsymmErrors();
  }
  //TFile *File_InPutXu = TFile::Open(inputfileXu.c_str());
  TF1 *f_res = new TF1("f_res",ResolutionFull,0,10,0);
  double Centrality_start[9] = {0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05, 0.0};
  double Centrality_stop[9]  = {0.8,0.7,0.6,0.5,0.4,0.3,0.2, 0.1,0.05};


  //string paperfile = "../ComparisonFiles/HEPData-ins1119620-v1-root.root";
  //TFile *File_Paper = TFile::Open(paperfile.c_str());
  //TDirectory *dir = (TDirectory*) File_Paper->Get("Table 2;1");
  //dir->cd();
  //TGraphAsymmErrors *g_mRes2Off_Raw = dir->Get("Graph1D_y3;1");
  TGraphAsymmErrors *g_mRes2Off = new TGraphAsymmErrors();
  //TGraphAsymmErrors *g_mRes2OffSub = new TGraphAsymmErrors();

  double CentVal[9] = {2.5,7.5,15,25,35,45,55,65,75};
  double LiKeVal[9] = {0.29528, 0.41235, 0.49976, 0.519075, 0.472266, 0.38668, 0.283204, 0.193584, 0.137737};
  double LiKeErr[9] = {0.000271561, 0.000192375, 0.000108864, 0.000103988, 0.000115842, 0.000148325, 0.000201894, 0.000315366, 0.000469525}; //from centrality 0-5% to 70-80%.
 
  //double val_19GeV[9]={0.271658,0.372193,0.453476,0.472727,0.429947,0.350802,0.256684,0.177754,0.126203};
  //double err_19GeV[9]={0.003,0.001,0.001,0.001,0.001,0.001,0.002,0.002,0.003}; 
  //for(int i_cent = 0; i_cent < 9; ++i_cent)
  //{
  //  double x = 0;
  //  double y = 0;
  //  g_mRes2Off_Raw->GetPoint(i_cent, x, y);
  //  g_mRes2Off->SetPoint(i_cent,CentVal[i_cent],y);
  //  g_mRes2Off->SetPointError(i_cent,0.0,0.0,g_mRes2Off_Raw->GetErrorYlow(i_cent),g_mRes2Off_Raw->GetErrorYhigh(i_cent));
  //}

  /*double mTpcSubRes1Val[9];
  double mTpcSubRes1Err[9];
  double mTpcFullRes1Val[9];
  double mTpcFullRes1Err[9];
  TGraphAsymmErrors *g_mTpcSubRes1 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_mTpcFullRes1 = new TGraphAsymmErrors();

  // calculate sub event plane resolution
  cout << "TPC Event Plane Resolution:" << endl;
  cout << "Sub Event Plane:" << endl;
  TProfile *p_mTpcSubRes1 = (TProfile*)File_InPut->Get("p_mTpcSubRes1");
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    const double resRaw = p_mTpcSubRes1->GetBinContent(p_mTpcSubRes1->FindBin(i_cent));
    const double errRaw = p_mTpcSubRes1->GetBinError(p_mTpcSubRes1->FindBin(i_cent));
    if(resRaw > 0)
    {
      mTpcSubRes1Val[i_cent] = TMath::Sqrt(resRaw);
      mTpcSubRes1Err[i_cent] = errRaw/(2.0*TMath::Sqrt(resRaw));
    }
    cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resSub = " << mTpcSubRes1Val[i_cent] << " +/- " << mTpcSubRes1Err[i_cent] << endl;
    g_mTpcSubRes1->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mTpcSubRes1Val[i_cent]*100.0);
    g_mTpcSubRes1->SetPointError(i_cent,0.0,0.0,mTpcSubRes1Err[i_cent]*100.0,mTpcSubRes1Err[i_cent]*100.0);
  }

  // calculate full event plane resolution
  cout << "Full Event Plane:" << endl;
  TProfile *p_mTpcRanRes1 = (TProfile*)File_InPut->Get("p_mTpcRanRes1");
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    const double resRaw = p_mTpcRanRes1->GetBinContent(p_mTpcRanRes1->FindBin(i_cent));
    const double errRaw = p_mTpcRanRes1->GetBinError(p_mTpcRanRes1->FindBin(i_cent));
    if(resRaw > 0)
    {
      const double resSub = TMath::Sqrt(resRaw);
      const double errSub = errRaw/(2.0*TMath::Sqrt(resRaw));

      const double chiSub = f_res->GetX(resSub);
      const double errChiSub = errSub/f_res->Derivative(chiSub);
      const double chiFull = chiSub*TMath::Sqrt(2.0);
      mTpcFullRes1Val[i_cent] = f_res->Eval(chiFull);
      mTpcFullRes1Err[i_cent] = f_res->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0);
    }
    cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resFull = " << mTpcFullRes1Val[i_cent] << " +/- " << mTpcFullRes1Err[i_cent] << endl;
    g_mTpcFullRes1->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mTpcFullRes1Val[i_cent]*100.0);
    g_mTpcFullRes1->SetPointError(i_cent,0.0,0.0,mTpcFullRes1Err[i_cent]*100.0,mTpcFullRes1Err[i_cent]*100.0);
  }
  */
  // TPC event plane resolution
  double mTpcSubRes2Val[7][9];
  double mTpcSubRes2Err[7][9];
  //double mTpcSubRes2ValXu[9];
  //double mTpcSubRes2ErrXu[9];
  //double mTpcFullRes2Val[9];
  //double mTpcFullRes2Err[9];
 // TGraphAsymmErrors *g_mTpcSubRes2Xu = new TGraphAsymmErrors();
 // TGraphAsymmErrors *g_mTpcFullRes2 = new TGraphAsymmErrors();

  // calculate sub event plane resolution
  cout << "TPC Event Plane Resolution:" << endl;
  cout << "Sub Event Plane:" << endl;
  TProfile *p_mTpcSubRes2[7];
  for(int i  = 6; i  < 7; i++)
  {
    p_mTpcSubRes2[i] = (TProfile*)File_InPut[i]->Get("p_mRes2_Sub");
   }
   //TProfile *p_mTpcSubRes2Xu = (TProfile*)File_InPutXu->Get("p_mRes2_Sub");

  for(int i  = 6; i  < 7; i++)
  {
    for(int i_cent = 0; i_cent < 9; ++i_cent)
    {
      const double resRaw = p_mTpcSubRes2[i]->GetBinContent(p_mTpcSubRes2[i]->FindBin(i_cent));
      const double errRaw = p_mTpcSubRes2[i]->GetBinError(p_mTpcSubRes2[i]->FindBin(i_cent));
      //const double resRawXu = p_mTpcSubRes2Xu->GetBinContent(p_mTpcSubRes2Xu->FindBin(i_cent));
      //const double errRawXu = p_mTpcSubRes2Xu->GetBinError(p_mTpcSubRes2Xu->FindBin(i_cent));
      if(resRaw > 0)
      {
        mTpcSubRes2Val[i][i_cent] = TMath::Sqrt(resRaw);
        mTpcSubRes2Err[i][i_cent] = errRaw/(2.0*TMath::Sqrt(resRaw));
        //mTpcSubRes2ValXu[i_cent] = TMath::Sqrt(resRawXu);
        //mTpcSubRes2ErrXu[i_cent] = errRaw/(2.0*TMath::Sqrt(resRawXu));
      }
      cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resSub = " << mTpcSubRes2Val[i][i_cent] << " +/- " << mTpcSubRes2Err[i][i_cent] << endl;
      g_mTpcSubRes2[i]->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mTpcSubRes2Val[i][i_cent]);
      g_mTpcSubRes2[i]->SetPointError(i_cent,0.0,0.0,mTpcSubRes2Err[i][i_cent],mTpcSubRes2Err[i][i_cent]);
      //g_mTpcSubRes2Xu->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mTpcSubRes2ValXu[i_cent]);
      //g_mTpcSubRes2Xu->SetPointError(i_cent,0.0,0.0,mTpcSubRes2ErrXu[i_cent],mTpcSubRes2ErrXu[i_cent]);
    }
  }
  // calculate full event plane resolution
  //cout << "Full Event Plane:" << endl;
  //TProfile *p_mTpcRanRes2 = (TProfile*)File_InPut->Get("p_mRes2_Ran");
  //for(int i_cent = 0; i_cent < 9; ++i_cent)
  //{
  //  const double resRaw = p_mTpcRanRes2->GetBinContent(p_mTpcRanRes2->FindBin(i_cent));
  //  const double errRaw = p_mTpcRanRes2->GetBinError(p_mTpcRanRes2->FindBin(i_cent));
  //  if(resRaw > 0)
  //  {
  //    const double resSub = TMath::Sqrt(resRaw);
  //    const double errSub = errRaw/(2.0*TMath::Sqrt(resRaw));

  //    const double chiSub = f_res->GetX(resSub);
  //    const double errChiSub = errSub/f_res->Derivative(chiSub);
  //    const double chiFull = chiSub*TMath::Sqrt(2.0);
  //    mTpcFullRes2Val[i_cent] = f_res->Eval(chiFull);
  //    mTpcFullRes2Err[i_cent] = f_res->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0);
  //  }
  //  cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resFull = " << mTpcFullRes2Val[i_cent] << " +/- " << mTpcFullRes2Err[i_cent] << endl;
  //  g_mTpcFullRes2->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mTpcFullRes2Val[i_cent]);
  //  g_mTpcFullRes2->SetPointError(i_cent,0.0,0.0,mTpcFullRes2Err[i_cent],mTpcFullRes2Err[i_cent]);
  //}

  //TGraphAsymmErrors *g_mRes2Ratio = new TGraphAsymmErrors();
  //for(int i_cent = 0; i_cent < 9; ++i_cent)
  //{
  //  g_mRes2Off->SetPoint(i_cent,CentVal[i_cent],LiKeVal[i_cent]);
  //  g_mRes2Off->SetPointError(i_cent,0.0,0.0,LiKeErr[i_cent],LiKeErr[i_cent]);
  //  g_mRes2Ratio->SetPoint(i_cent,CentVal[i_cent],mTpcSubRes2Val[8-i_cent]/LiKeVal[i_cent]);
  //  double error1 = mTpcSubRes2Err[8-i_cent]/mTpcSubRes2Val[8-i_cent];
  //  double error2 = LiKeErr[i_cent]/LiKeVal[i_cent];
  //  double error = mTpcSubRes2Val[8-i_cent]/LiKeVal[i_cent]*TMath::Sqrt(error1*error1+error2*error2);

  //  g_mRes2Ratio->SetPointError(i_cent,0.0,0.0,error,error);
  //}
  //TGraphAsymmErrors *g_mRes2SubOff = new TGraphAsymmErrors();
  //for(int i_cent = 0; i_cent < 9; ++i_cent)
 // {
   // g_mRes2SubOff->SetPoint(i_cent,CentVal[i_cent],val_19GeV[i_cent]);
  //  g_mRes2SubOff->SetPointError(i_cent,0.0,0.0,err_19GeV[i_cent],err_19GeV[i_cent]);
  //}
/* 
  double mTpcSubRes3Val[9];
  double mTpcSubRes3Err[9];
  double mTpcFullRes3Val[9];
  double mTpcFullRes3Err[9];
  TGraphAsymmErrors *g_mTpcSubRes3 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_mTpcFullRes3 = new TGraphAsymmErrors();

  // calculate sub event plane resolution
  cout << "TPC Event Plane Resolution:" << endl;
  cout << "Sub Event Plane:" << endl;
  TProfile *p_mTpcSubRes3 = (TProfile*)File_InPut->Get("p_mTpcSubRes3");
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    const double resRaw = p_mTpcSubRes3->GetBinContent(p_mTpcSubRes3->FindBin(i_cent));
    const double errRaw = p_mTpcSubRes3->GetBinError(p_mTpcSubRes3->FindBin(i_cent));
    if(resRaw > 0)
    {
      mTpcSubRes3Val[i_cent] = TMath::Sqrt(resRaw);
      mTpcSubRes3Err[i_cent] = errRaw/(2.0*TMath::Sqrt(resRaw));
    }
    cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resSub = " << mTpcSubRes3Val[i_cent] << " +/- " << mTpcSubRes3Err[i_cent] << endl;
    g_mTpcSubRes3->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mTpcSubRes3Val[i_cent]*100.0);
    g_mTpcSubRes3->SetPointError(i_cent,0.0,0.0,mTpcSubRes3Err[i_cent]*100.0,mTpcSubRes3Err[i_cent]*100.0);
  }

  // calculate full event plane resolution
  cout << "Full Event Plane:" << endl;
  TProfile *p_mTpcRanRes3 = (TProfile*)File_InPut->Get("p_mTpcRanRes3");
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    const double resRaw = p_mTpcRanRes3->GetBinContent(p_mTpcRanRes3->FindBin(i_cent));
    const double errRaw = p_mTpcRanRes3->GetBinError(p_mTpcRanRes3->FindBin(i_cent));
    if(resRaw > 0)
    {
      const double resSub = TMath::Sqrt(resRaw);
      const double errSub = errRaw/(2.0*TMath::Sqrt(resRaw));

      const double chiSub = f_res->GetX(resSub);
      const double errChiSub = errSub/f_res->Derivative(chiSub);
      const double chiFull = chiSub*TMath::Sqrt(2.0);
      mTpcFullRes3Val[i_cent] = f_res->Eval(chiFull);
      mTpcFullRes3Err[i_cent] = f_res->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0);
    }
    cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resFull = " << mTpcFullRes3Val[i_cent] << " +/- " << mTpcFullRes3Err[i_cent] << endl;
    g_mTpcFullRes3->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mTpcFullRes3Val[i_cent]*100.0);
    g_mTpcFullRes3->SetPointError(i_cent,0.0,0.0,mTpcFullRes3Err[i_cent]*100.0,mTpcFullRes3Err[i_cent]*100.0);
  }
*/
  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.15);
  c_play->SetGrid(0,0);
  c_play->SetTicks(1,1);
  c_play->cd();

  TH1F *h_play = new TH1F("h_play","h_play",100,0,100);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("Centrality (%)");
  h_play->GetYaxis()->SetTitle("R_{2}");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(0,80);
  h_play->GetYaxis()->SetRangeUser(-0.01,1.0);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");
  h_play->Draw("pE");

  //g_mTpcSubRes1->SetMarkerStyle(24);
  //g_mTpcSubRes1->SetMarkerColor(kGreen+2);
  //g_mTpcSubRes1->SetMarkerSize(1.5);
  //g_mTpcSubRes1->Draw("pE Same");

  //g_mTpcFullRes1->SetMarkerStyle(20);
  //g_mTpcFullRes1->SetMarkerColor(kGreen+2);
  //g_mTpcFullRes1->SetMarkerSize(1.5);
  //g_mTpcFullRes1->Draw("pE Same");

  //g_mTpcSubRes2Xu->SetMarkerStyle(20);
  //g_mTpcSubRes2Xu->SetMarkerColor(kBlack);
  //g_mTpcSubRes2Xu->SetMarkerSize(1.5);
  //g_mTpcSubRes2Xu->Draw("pE Same");
 
  //g_mTpcSubRes2[3]->SetMarkerStyle(20);
  //g_mTpcSubRes2[3]->SetMarkerColor(kOrange+7);
  //g_mTpcSubRes2[3]->SetMarkerSize(2);
  //g_mTpcSubRes2[3]->Draw("pE Same");

  //g_mTpcSubRes2[4]->SetMarkerStyle(20);
  //g_mTpcSubRes2[4]->SetMarkerColor(kBlue);
  //g_mTpcSubRes2[4]->SetMarkerSize(2);
  //g_mTpcSubRes2[4]->Draw("pE Same");

  //g_mTpcSubRes2[5]->SetMarkerStyle(20);
  //g_mTpcSubRes2[5]->SetMarkerColor(kBlue);
  //g_mTpcSubRes2[5]->SetMarkerSize(2);
  //g_mTpcSubRes2[5]->Draw("pE Same");
  g_mTpcSubRes2[6]->SetMarkerStyle(20);
  g_mTpcSubRes2[6]->SetMarkerColor(kBlue);
  g_mTpcSubRes2[6]->SetMarkerSize(2);
  g_mTpcSubRes2[6]->Draw("pE Same");

  //g_mTpcFullRes2->SetMarkerStyle(20);
  //g_mTpcFullRes2->SetMarkerColor(kAzure+2);
  //g_mTpcFullRes2->SetMarkerSize(1.5);
  //g_mTpcFullRes2->Draw("pE Same");

  //g_mRes2Off->SetMarkerStyle(20);
  //g_mRes2Off->SetMarkerColor(kBlue);
  //g_mRes2Off->SetMarkerSize(1.5);
  //g_mRes2Off->Draw("pE Same");

  //g_mRes2SubOff->SetMarkerStyle(20);
  //g_mRes2SubOff->SetMarkerColor(kBlue);
  //g_mRes2SubOff->SetMarkerSize(1.5);
  //g_mRes2SubOff->Draw("pE Same");

  //g_mTpcSubRes3->SetMarkerStyle(24);
  //g_mTpcSubRes3->SetMarkerColor(kGray+2);
  //g_mTpcSubRes3->SetMarkerSize(1.5);
  //g_mTpcSubRes3->Draw("pE Same");

  //g_mTpcFullRes3->SetMarkerStyle(20);
  //g_mTpcFullRes3->SetMarkerColor(kGray+2);
  //g_mTpcFullRes3->SetMarkerSize(1.5);
  //g_mTpcFullRes3->Draw("pE Same");

  TLegend *leg = new TLegend(0.55,0.65,0.85,0.85);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  //leg->AddEntry(g_mTpcFullRes1,"1^{st} Full TPC EP","p");
  //leg->AddEntry(g_mTpcSubRes1,"1^{st} #eta_{sub} TPC EP","p");
  //leg->AddEntry(g_mTpcFullRes2,"2^{nd} Full TPC EP","p");
  //leg->AddEntry(g_mTpcSubRes2[3],"14.6 GeV #eta_{sub} TPC","p");
  //leg->AddEntry(g_mTpcSubRes2[4],"19.6 GeV #eta_{sub} TPC","p");
  //leg->AddEntry(g_mTpcSubRes2[5],"27 GeV #eta_{sub} TPC","p");
  leg->AddEntry(g_mTpcSubRes2[6],"200 GeV Run14 #eta_{sub} TPC","p");
  //leg->AddEntry(g_mTpcSubRes2Xu,"#eta_{sub} TPC Xu Sun's BES-I","p");
  //leg->AddEntry(g_mTpcFullRes2,"2^{nd} Full Ran TPC","p");
  //leg->AddEntry(g_mRes2Off,"2^{nd} Full Ran TPC, BES-I","p");
  //leg->AddEntry(g_mRes2Off,"2^{nd} #eta_{sub} TPC Li-Ke (|#eta|<1.5)","p");
  //leg->AddEntry(g_mRes2SubOff,"2^{nd} #eta_{sub} TPC BES-I Official","p");
  //leg->AddEntry(g_mTpcFullRes3,"3^{rd} Full TPC EP","p");
  //leg->AddEntry(g_mTpcSubRes3,"3^{rd} #eta_{sub} TPC EP","p");//leg->AddEntry(g_mZdcFullRes1,"1^{st} ZDC-SMD Full EP","p");
  //leg->AddEntry(g_mZdcFullRes2,"2^{nd} ZDC-SMD Full EP","p");
  leg->Draw("same");

  string FigureName = Form("./figures/c_mEpResolution.pdf");
  c_play->SaveAs(FigureName.c_str());



  //TCanvas *c_play2 = new TCanvas("c_play2","c_play2",10,10,800,800);
  //c_play2->SetLeftMargin(0.15);
  //c_play2->SetBottomMargin(0.15);
  //c_play2->SetGrid(0,0);
  //c_play2->SetTicks(1,1);
  //c_play2->cd();

  //TH1F *h_play = new TH1F("h_play","h_play",100,0,100);
  //for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  //{
  //  h_play->SetBinContent(i_bin+1,-10.0);
  //  h_play->SetBinError(i_bin+1,1.0);
  //}
  //h_play->SetTitle("");
  //h_play->SetStats(0);
  //h_play->GetXaxis()->SetTitle("Centrality (%)");
  //h_play->GetYaxis()->SetTitle("Gavin/LiKe");
  //h_play->GetXaxis()->CenterTitle();
  //h_play->GetYaxis()->CenterTitle();
  //h_play->GetXaxis()->SetTitleSize(0.06);
  //h_play->GetYaxis()->SetTitleSize(0.06);
  //h_play->GetXaxis()->SetRangeUser(0,80);
  //h_play->GetYaxis()->SetRangeUser(0.95,1.05);
  //h_play->GetXaxis()->SetLabelSize(0.04);
  //h_play->GetYaxis()->SetLabelSize(0.04);
  //h_play->SetNdivisions(505,"X");
  //h_play->SetNdivisions(505,"Y");
  //h_play->Draw("pE");

  //g_mRes2Ratio->SetMarkerStyle(20);
  //g_mRes2Ratio->SetMarkerColor(kBlue);
  //g_mRes2Ratio->SetMarkerSize(1.5);
  //g_mRes2Ratio->Draw("pE Same");
  //g_mRes2Ratio->Draw("pE same");
}
