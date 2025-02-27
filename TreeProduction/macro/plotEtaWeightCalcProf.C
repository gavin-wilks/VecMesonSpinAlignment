#include "../StRoot/Utility/StSpinAlignmentCons.h"
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include <fstream>
#include <iostream>
#include "TGraphAsymmErrors.h"
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

void plotEtaWeightCalcProf(int beamEnergy = 0)
{
  //string inputfile = Form("../StRoot/StEventPlaneUtility/ChargedFlow/file_%s_ChargedFlow.root",vmsa::mBeamEnergy[beamEnergy].c_str());
  //string inputfileEpd = Form("../StRoot/Utility/EpdResolution/file_%s_EpdResEta_2.root",vmsa::mBeamEnergy[beamEnergy].c_str());
  string inputfileEpd = Form("../StRoot/Utility/EpdResolution/file_%s_EpdFlow_4.root",vmsa::mBeamEnergy[beamEnergy].c_str());
 // string inputfileEpd = Form("../StRoot/Utility/EpdResolution/file_19p6GeV_2019_EpdFlow_6.root");
  //string inputfileTpc = Form("../StRoot/StEventPlaneUtility/Resolution/file_%s_Resolution.root",vmsa::mBeamEnergy[beamEnergy].c_str());
  //string inputfileEpd = Form("../StRoot/StEventPlaneUtility/EpdFlow/file_%s_EpdFlow.root",vmsa::mBeamEnergy[beamEnergy].c_str());

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TFile *File_InPut_Epd = TFile::Open(inputfileEpd.c_str());
  //TFile *File_InPut_Tpc = TFile::Open(inputfileTpc.c_str());

  std::string centrality[10] = {"70%-80%","60%-70%","50%-60%","40%-50%","30%-40%","20%-30%","10%-20%","5%-10%","0%-5%","All"};

  TProfile *p_mTpcVEta[3][10]; // <v*Ep> vs. eta
  TProfile *p_mTpcFlowEta[3][10];
  TProfile *p_mTpcSubRes[3]; 
  double mTpcSubResVal[3][10];
  double mTpcSubResErr[3][10];

  TProfile *p_mEpdVEta[3][10];
  TProfile *p_mEpdFlowEta[3][10];
  TProfile *p_mEpdFlowEtaWeights[3][10];
  TProfile *p_mEpdSubRes[3];
  double mEpdSubResVal[3][10];
  double mEpdSubResErr[3][10];

  TProfile *p_mFlowEta[3][10];
  TGraphAsymmErrors *g_mFlowEta[3][10];  

  for(int order = 1; order <= 1; ++order)
  {
    //TPC
    //string ProName = Form("p_mTpcSubRes%d",order);
    //p_mTpcSubRes[order-1] = (TProfile*) File_InPut_Tpc->Get(ProName.c_str());
  
    //for(int i_cent = 0; i_cent < 9; ++i_cent) 
    //{ 
    //  const double resRaw = p_mTpcSubRes[order-1]->GetBinContent(p_mTpcSubRes[order-1]->FindBin(i_cent)); 
    //  const double errRaw = p_mTpcSubRes[order-1]->GetBinError(p_mTpcSubRes[order-1]->FindBin(i_cent)); 
    //  cout << "TPC: order " << order << "  cent " << i_cent << "   Res Raw = " << resRaw << "   Err = " << errRaw << endl;
    //  //if(order == 1) resRaw*=-1;
    //  if(resRaw > 0) 
    //  { 
    //    mTpcSubResVal[order-1][i_cent] = TMath::Sqrt(resRaw); 
    //    mTpcSubResErr[order-1][i_cent] = errRaw/(2.0*TMath::Sqrt(resRaw)); 
    //  } 
    //}

    //EPD
    std::string ProName = Form("p_mEpdSubRes%d",order);
    p_mEpdSubRes[order-1] = (TProfile*) File_InPut_Epd->Get(ProName.c_str());

    for(int i_cent = 0; i_cent < 9; ++i_cent) 
    { 
      const double resRaw = p_mEpdSubRes[order-1]->GetBinContent(p_mEpdSubRes[order-1]->FindBin(i_cent)); 
      const double errRaw = p_mEpdSubRes[order-1]->GetBinError(p_mEpdSubRes[order-1]->FindBin(i_cent)); 
      if(resRaw > 0) 
      { 
        mEpdSubResVal[order-1][i_cent] = TMath::Sqrt(resRaw); 
        mEpdSubResErr[order-1][i_cent] = errRaw/(2.0*TMath::Sqrt(resRaw));    
      cout << "EPD: order " << order << "  cent " << i_cent << "   Res Raw = " << TMath::Sqrt(resRaw) << "   Err = " << endl;
      } 
    }

    for(int i_cent = 0; i_cent < 9; ++i_cent)
    {
      //TPC
      //ProName = Form("p_mTpcFlowEta%d_Cent%d",order,i_cent);
      //p_mTpcFlowEta[order-1][i_cent] = (TProfile*)File_InPut_Tpc->Get(ProName.c_str());

      //ProName = Form("p_mTpcV%dEta_Cent%d",order,i_cent);
      //p_mTpcVEta[order-1][i_cent] = (TProfile*)p_mTpcFlowEta[order-1][i_cent]->Clone(ProName.c_str()); 
      //p_mTpcVEta[order-1][i_cent]->Scale(1.0/mTpcSubResVal[order-1][i_cent]);

      //EPD
      ProName = Form("p_mEpdFlowEta%d_Cent%d",order,i_cent);
      p_mEpdFlowEta[order-1][i_cent] = (TProfile*)File_InPut_Epd->Get(ProName.c_str());  

      ProName = Form("p_mEpdFlowEtaWeights%d_Cent%d",order,i_cent);
      p_mEpdFlowEtaWeights[order-1][i_cent] = (TProfile*)File_InPut_Epd->Get(ProName.c_str());  
      
      ProName = Form("p_mEpdV%dEta_Cent%d",order,i_cent);
      p_mEpdVEta[order-1][i_cent] = (TProfile*)p_mEpdFlowEta[order-1][i_cent]->Clone(ProName.c_str()); 
      TH1F* hWeightsX = (TH1F*) p_mEpdFlowEtaWeights[order-1][i_cent]->ProjectionX("");
      p_mEpdVEta[order-1][i_cent]->Divide(hWeightsX);
      p_mEpdVEta[order-1][i_cent]->Scale(1.0/mEpdSubResVal[order-1][i_cent]); 
      cout << "Order: " << order << "   Cent: " << i_cent << endl;
      //TPC and EPD
      ///ProName = Form("p_mV%dEta_Cent%d",order,i_cent);
      ///p_mFlowEta[order-1][i_cent] = (TProfile*) p_mTpcVEta[order-1][i_cent]->Clone(ProName.c_str());
      //p_mFlowEta[order-1][i_cent]->Add(p_mEpdVEta[order-1][i_cent]); 
      p_mFlowEta[order-1][i_cent] = (TProfile*) p_mEpdVEta[order-1][i_cent]->Clone(ProName.c_str());
      
      //g_mFlowEta[order-1][i_cent] = new TGraphAsymmErrors();
      //for(int ibin = 1; ibin <= p_mFlowEta[order-1][i_cent]->GetNbinsX(); ibin++)
      //{
      //  double binC = p_mFlowEta[order-1][i_cent]->GetBinCenter(ibin);
      //  double val  = p_mFlowEta[order-1][i_cent]->GetBinContent(ibin);
      //  double err  = p_mFlowEta[order-1][i_cent]->GetBinError(ibin);
      //  g_mFlowEta[order-1][i_cent]->SetPoint(ibin,binC,val);
      //  g_mFlowEta[order-1][i_cent]->SetPointError(ibin,0.0,0.0,err,err);
      ////  cout << "BinC: " << binC << "    val: " << val << "    err: " << err << endl;
      //}
    }
  }

  /*TH1F *h_play = new TH1F("h_play","h_play",100,-5.5,5.5);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("#eta");
  h_play->GetYaxis()->SetTitle("v(%)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(-5.5,5.5);
  h_play->GetYaxis()->SetRangeUser(-1.0,1.0);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");
*/
  std::string outputname = Form("./figures/EtaWeightsVerson3_%s.pdf",vmsa::mBeamEnergy[beamEnergy].c_str());

  TCanvas *c_v = new TCanvas("c_v","c_v2",2400,1800);
  c_v->Divide(3,3);
  for(int ipad = 0; ipad < 9; ipad++)
  {
    c_v->cd(ipad+1)->SetLeftMargin(0.15);
    c_v->cd(ipad+1)->SetBottomMargin(0.15);
    c_v->cd(ipad+1)->SetGrid(0,0);
    c_v->cd(ipad+1)->SetTicks(1,1);
  }

  std::string output_start = Form("%s[",outputname.c_str());
  c_v->Print(output_start.c_str(),"pdf");

  std::ofstream outfile;
  outfile.open("epdweights.txt");

  TLegend *leg[3][10]; 
  TF1 *fit[3][10];
  for(int order = 1; order <= 1; order++)
  {
    for(int i_cent = 0; i_cent < 9; i_cent++)
    {
      c_v->cd(i_cent+1)->Clear();
      c_v->cd(i_cent+1);
      //cout << "Start plotting" << endl;
      //h_play->DrawCopy("pE");

      fit[order-1][i_cent] = new TF1(Form("f%d%d",order,i_cent),"[0]*x+[1]*x*x*x+[2]*x*x*x*x*x",0,5.4);
      if(i_cent < 3) fit[order-1][i_cent] = new TF1(Form("f%d%d",order,i_cent),"[0]*x+[1]*x*x*x+[2]*x*x*x*x*x",0,5.5);
      fit[order-1][i_cent]->SetLineWidth(1);
      fit[order-1][i_cent]->SetLineColor(2);
      /*p_mFlowEta[order-1][i_cent]->Fit(fit[order-1][i_cent]);
      p_mTpcVEta[order-1][i_cent]->SetTitle(centrality[i_cent].c_str());
      p_mTpcVEta[order-1][i_cent]->GetXaxis()->SetTitle("#eta");
      p_mTpcVEta[order-1][i_cent]->GetYaxis()->SetTitle(Form("v_{%d}",order));
      p_mTpcVEta[order-1][i_cent]->SetMarkerStyle(24);
      p_mTpcVEta[order-1][i_cent]->SetMarkerColor(kGray+2);
      p_mTpcVEta[order-1][i_cent]->SetMarkerSize(1.5);
      p_mTpcVEta[order-1][i_cent]->SetLineColor(kGray+2);
      p_mTpcVEta[order-1][i_cent]->Draw("pE same");
*/
      p_mFlowEta[order-1][i_cent]->Fit(fit[order-1][i_cent],"MIER");
      p_mFlowEta[order-1][i_cent]->SetTitle(centrality[i_cent].c_str());
      p_mFlowEta[order-1][i_cent]->GetXaxis()->SetTitle("#eta");
      p_mFlowEta[order-1][i_cent]->GetYaxis()->SetTitle(Form("v_{%d}",order));
      p_mFlowEta[order-1][i_cent]->SetMarkerStyle(24);
      p_mFlowEta[order-1][i_cent]->SetMarkerColor(kGray+2);
      p_mFlowEta[order-1][i_cent]->SetMarkerSize(1.5);
      p_mFlowEta[order-1][i_cent]->SetLineColor(kGray+2);
      p_mFlowEta[order-1][i_cent]->Draw("pE");
 

  /*    p_mEpdVEta[order-1][i_cent]->SetMarkerStyle(20);
      p_mEpdVEta[order-1][i_cent]->SetMarkerColor(kGray+2);
      p_mEpdVEta[order-1][i_cent]->SetMarkerSize(1.5);
      p_mEpdVEta[order-1][i_cent]->SetLineColor(kGray+2);
      p_mEpdVEta[order-1][i_cent]->Draw("pE same");
*/
      //fit[order-1][i_cent]->Draw("C same");
    
      leg[order-1][i_cent] = new TLegend(0.20,0.70,0.45,0.85);
      leg[order-1][i_cent]->SetFillColor(10);
      leg[order-1][i_cent]->SetBorderSize(0);
      //leg[order-1][i_cent]->AddEntry(p_mTpcVEta[order-1][i_cent],"TPC","p");
      //leg[order-1][i_cent]->AddEntry(p_mEpdVEta[order-1][i_cent],"EPD","p");
      leg[order-1][i_cent]->AddEntry(fit[order-1][i_cent],"p0*x+p1*x^{3}+p2*x^{5}","l");
      leg[order-1][i_cent]->Draw("same");

      outfile << fit[order-1][i_cent]->GetParameter(0) << " " << fit[order-1][i_cent]->GetParameter(1) << " " << fit[order-1][i_cent]->GetParameter(2) << "\n";
      if(i_cent == 8) outfile << "\n";
      //cout << "One plot loop" << endl;*/
      //c_v->Update();
      //c_v->Print(outputname.c_str(),"pdf");

    } 

    c_v->Update();
    c_v->Print(outputname.c_str(),"pdf");
  }
  outfile.close();
  std::string output_stop = Form("%s]",outputname.c_str());
  c_v->Print(output_stop.c_str(),"pdf");
}
