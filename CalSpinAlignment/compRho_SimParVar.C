#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
//#include "../Utility/draw.h"
#include "../../../draw.h"
#include "../Utility/functions.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void compRho_SimParVar(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "AccRes")
{

  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  string inputfile = Form("../output/AuAu%s/%s/Poly/%sPhiPtSys_eta1_eta1_Poly_SimParVar.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());

  TFile *File_InPut = TFile::Open(inputfile.c_str());
  
  const int nfiles = 38;
  const int defaultfile = 2;

  TGraphAsymmErrors *g_Stat[nfiles+1];
  for(int ifile = 1; ifile <= nfiles; ifile++)
  {
    TString StatErrorRho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1_Setting%d",9,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),ifile-1);
    cout << StatErrorRho << endl;
    g_Stat[ifile] = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut->Get(StatErrorRho.Data()))->Clone();
  }

  cout << "Loaded all the files " << endl;

  TCanvas *c_rho = new TCanvas("c_rho_SysError","c_rho_SysError",10,10,1200,600);
  c_rho->Divide(4,2);
  for(int i = 0; i < 9; i++)
  {
    c_rho->cd(i+1);
    c_rho->cd(i+1)->SetLeftMargin(0.15);
    c_rho->cd(i+1)->SetBottomMargin(0.15);
    c_rho->cd(i+1)->SetTicks(1,1);
    c_rho->cd(i+1)->SetGrid(0,0);
  }
 
  int color[7] = {kBlack,kOrange+7,kBlue,kGray+2,kRed,kGreen,kViolet} ;

  {  
    // global rho00 input 
    c_rho->cd(1);
    TLegend *leg_rhog1 = new TLegend(0.2,0.2,0.45,0.45);
    leg_rhog1->SetHeader("#rho_{00}^{global}","C");
    const int nglobal = 4; 
    std::string rhogtext[nglobal] = {"0.3333 (Default)","0.3000","0.3667","0.4000"};
    int s_rhog[nglobal] = {3,2,4,5};  

    for(int i = 0; i < nglobal; i++)
    {
      int ifile = s_rhog[i]; 
      g_Stat[ifile]->SetMarkerStyle(20);
      g_Stat[ifile]->SetMarkerColor(color[i]);
      g_Stat[ifile]->SetLineColor(color[i]);
      g_Stat[ifile]->GetXaxis()->SetTitle("p_{T}");
      g_Stat[ifile]->GetYaxis()->SetTitle("global #rho_{00}");
      g_Stat[ifile]->GetYaxis()->SetRangeUser(0.29,0.40);
      if(i != 0)
      {
        for(int ip = 0; ip < g_Stat[ifile]->GetN(); ip++)
        {
          double pt, rho;
          g_Stat[ifile]->GetPoint(ip,pt,rho);
          g_Stat[ifile]->SetPoint(ip,pt+0.04*i,rho); 
        }
      }
      if(i == 0) g_Stat[ifile]->Draw("APE");
      else       g_Stat[ifile]->Draw("PE same");
      leg_rhog1->AddEntry(g_Stat[ifile],rhogtext[i].c_str(),"p");
    }
    leg_rhog1->Draw("same");
  }
  cout << "Study 1 " << endl;
  {  
    // helicity rho00 input 
    c_rho->cd(2);
    TLegend *leg_rhog2 = new TLegend(0.5,0.6,0.75,0.85);
    leg_rhog2->SetHeader("#rho_{00}^{helicity}","C");
    const int nglobal = 5;
    std::string rhogtext[nglobal] = {"1/3 (Default)","0.20","0.25","0.30","0.35"};
    int s_rhog[nglobal] = {3,6,7,8,9};  

    for(int i = 0; i < nglobal; i++)
    {
      int ifile = s_rhog[i]; 
      g_Stat[ifile]->SetMarkerStyle(20);
      g_Stat[ifile]->SetMarkerColor(color[i]);
      g_Stat[ifile]->SetLineColor(color[i]);
      g_Stat[ifile]->GetXaxis()->SetTitle("p_{T}");
      g_Stat[ifile]->GetYaxis()->SetTitle("global #rho_{00}");
      g_Stat[ifile]->GetYaxis()->SetRangeUser(0.29,0.40);
      if(i != 0)
      {
        for(int ip = 0; ip < g_Stat[ifile]->GetN(); ip++)
        {
          double pt, rho;
          g_Stat[ifile]->GetPoint(ip,pt,rho);
          g_Stat[ifile]->SetPoint(ip,pt+0.04*i,rho); 
        }
      }
      if(i == 0) g_Stat[ifile]->Draw("APE");
      else       g_Stat[ifile]->Draw("PE same");
      leg_rhog2->AddEntry(g_Stat[ifile],rhogtext[i].c_str(),"p");
    }
    leg_rhog2->Draw("same");
  }
  cout << "Study 2 " << endl;
  //{  
  //  // global and helicity cos(theta*) weight 
  //  c_rho->cd(3);
  //  TLegend *leg_rhog3 = new TLegend(0.2,0.2,0.45,0.45);
  //  leg_rhog3->SetHeader("Weighting","C");
  //  const int nglobal = 2; 
  //  std::string rhogtext[nglobal] = {"No weight (Default)","Weight 2D cos(#theta*)"};
  //  int s_rhog[nglobal] = {3,10};  

  //  for(int i = 0; i < nglobal; i++)
  //  {
  //    int ifile = s_rhog[i]; 
  //    g_Stat[ifile]->SetMarkerStyle(20);
  //    g_Stat[ifile]->SetMarkerColor(color[i]);
  //    g_Stat[ifile]->SetLineColor(color[i]);
  //    g_Stat[ifile]->GetXaxis()->SetTitle("p_{T}");
  //    g_Stat[ifile]->GetYaxis()->SetTitle("global #rho_{00}");
  //    g_Stat[ifile]->GetYaxis()->SetRangeUser(0.29,0.40);
  //    if(i == 0) g_Stat[ifile]->Draw("APE");
  //    else       g_Stat[ifile]->Draw("PE same");
  //    leg_rhog3->AddEntry(g_Stat[ifile],rhogtext[i].c_str(),"p");
  //  }
  //  leg_rhog3->Draw("same");
  //}
  cout << "Study 3 " << endl;
  {  
    // spectra reweighting 
    c_rho->cd(4);
    TLegend *leg_rhog4 = new TLegend(0.2,0.2,0.45,0.45);
    leg_rhog4->SetHeader("Reweighting","C");
    const int nglobal = 3; 
    std::string rhogtext[nglobal] = {"no weight (Default)","#phi reweight (pt,y)","kaon reweight (pt,y)"};
    int s_rhog[nglobal] = {3,37,38};  

    for(int i = 0; i < nglobal; i++)
    {
      int ifile = s_rhog[i]; 
      g_Stat[ifile]->SetMarkerStyle(20);
      g_Stat[ifile]->SetMarkerColor(color[i]);
      g_Stat[ifile]->SetLineColor(color[i]);
      g_Stat[ifile]->GetXaxis()->SetTitle("p_{T}");
      g_Stat[ifile]->GetYaxis()->SetTitle("global #rho_{00}");
      g_Stat[ifile]->GetYaxis()->SetRangeUser(0.29,0.40);
      if(i != 0)
      {
        for(int ip = 0; ip < g_Stat[ifile]->GetN(); ip++)
        {
          double pt, rho;
          g_Stat[ifile]->GetPoint(ip,pt,rho);
          g_Stat[ifile]->SetPoint(ip,pt+0.04*i,rho); 
        }
      }
      if(i == 0) g_Stat[ifile]->Draw("APE");
      else       g_Stat[ifile]->Draw("PE same");
      leg_rhog4->AddEntry(g_Stat[ifile],rhogtext[i].c_str(),"p");
    }
    leg_rhog4->Draw("same");
  }
  cout << "Study 4 " << endl;
  {  
    // Re(rho1-1) input 
    c_rho->cd(5);
    TLegend *leg_rhog5 = new TLegend(0.2,0.2,0.45,0.45);
    leg_rhog5->SetHeader("Re(#rho_{1,-1})","C");
    const int nglobal = 7; 
    std::string rhogtext[nglobal] = {"0.0 (Default)","-0.3","-0.2","-0.1","0.1","0.2","0.3"};
    int s_rhog[nglobal] = {3,13,14,15,16,17,18};  

    for(int i = 0; i < nglobal; i++)
    {
      int ifile = s_rhog[i]; 
      g_Stat[ifile]->SetMarkerStyle(20);
      g_Stat[ifile]->SetMarkerColor(color[i]);
      g_Stat[ifile]->SetLineColor(color[i]);
      g_Stat[ifile]->GetXaxis()->SetTitle("p_{T}");
      g_Stat[ifile]->GetYaxis()->SetTitle("global #rho_{00}");
      g_Stat[ifile]->GetYaxis()->SetRangeUser(0.29,0.40);
      if(i != 0)
      {
        for(int ip = 0; ip < g_Stat[ifile]->GetN(); ip++)
        {
          double pt, rho;
          g_Stat[ifile]->GetPoint(ip,pt,rho);
          g_Stat[ifile]->SetPoint(ip,pt+0.04*i,rho); 
        }
      }
      if(i == 0) g_Stat[ifile]->Draw("APE");
      else       g_Stat[ifile]->Draw("PE same");
      leg_rhog5->AddEntry(g_Stat[ifile],rhogtext[i].c_str(),"p");
    }
    leg_rhog5->Draw("same");
  }
  cout << "Study 5 " << endl;
  {  
    // Re(rho1-1) input 
    c_rho->cd(6);
    TLegend *leg_rhog6 = new TLegend(0.2,0.2,0.45,0.45);
    leg_rhog6->SetHeader("Im(#rho_{1,-1})","C");
    const int nglobal = 7; 
    std::string rhogtext[nglobal] = {"0.0 (Default)","-0.3","-0.2","-0.1","0.1","0.2","0.3"};
    int s_rhog[nglobal] = {3,31,32,33,34,35,36};  

    for(int i = 0; i < nglobal; i++)
    {
      int ifile = s_rhog[i]; 
      g_Stat[ifile]->SetMarkerStyle(20);
      g_Stat[ifile]->SetMarkerColor(color[i]);
      g_Stat[ifile]->SetLineColor(color[i]);
      g_Stat[ifile]->GetXaxis()->SetTitle("p_{T}");
      g_Stat[ifile]->GetYaxis()->SetTitle("global #rho_{00}");
      g_Stat[ifile]->GetYaxis()->SetRangeUser(0.29,0.40);
      if(i != 0)
      {
        for(int ip = 0; ip < g_Stat[ifile]->GetN(); ip++)
        {
          double pt, rho;
          g_Stat[ifile]->GetPoint(ip,pt,rho);
          g_Stat[ifile]->SetPoint(ip,pt+0.04*i,rho); 
        }
      }
      if(i == 0) g_Stat[ifile]->Draw("APE");
      else       g_Stat[ifile]->Draw("PE same");
      leg_rhog6->AddEntry(g_Stat[ifile],rhogtext[i].c_str(),"p");
    }
    leg_rhog6->Draw("same");
  }
  cout << "Study 6 " << endl;
  {  
    // Re(rho1-1) input 
    c_rho->cd(7);
    TLegend *leg_rhog7 = new TLegend(0.2,0.2,0.45,0.45);
    leg_rhog7->SetHeader("Re(#rho_{1,0})-Re(#rho_{0,-1})","C");
    const int nglobal = 7; 
    std::string rhogtext[nglobal] = {"0.0 (Default)","-0.3","-0.2","-0.1","0.1","0.2","0.3"};
    int s_rhog[nglobal] = {3,19,20,21,22,23,24};  

    for(int i = 0; i < nglobal; i++)
    {
      int ifile = s_rhog[i]; 
      g_Stat[ifile]->SetMarkerStyle(20);
      g_Stat[ifile]->SetMarkerColor(color[i]);
      g_Stat[ifile]->SetLineColor(color[i]);
      g_Stat[ifile]->GetXaxis()->SetTitle("p_{T}");
      g_Stat[ifile]->GetYaxis()->SetTitle("global #rho_{00}");
      g_Stat[ifile]->GetYaxis()->SetRangeUser(0.29,0.40);
      if(i != 0)
      {
        for(int ip = 0; ip < g_Stat[ifile]->GetN(); ip++)
        {
          double pt, rho;
          g_Stat[ifile]->GetPoint(ip,pt,rho);
          g_Stat[ifile]->SetPoint(ip,pt+0.04*i,rho); 
        }
      }
      if(i == 0) g_Stat[ifile]->Draw("APE");
      else       g_Stat[ifile]->Draw("PE same");
      leg_rhog7->AddEntry(g_Stat[ifile],rhogtext[i].c_str(),"p");
    }
    leg_rhog7->Draw("same");
  }
  cout << "Study 7 " << endl;
  {  
    // Re(rho1-1) input 
    c_rho->cd(8);
    TLegend *leg_rhog8 = new TLegend(0.2,0.2,0.45,0.45);
    leg_rhog8->SetHeader("Im(#rho_{1,0})-Im(#rho_{0,-1})","C");
    const int nglobal = 7; 
    std::string rhogtext[nglobal] = {"0.0 (Default)","-0.3","-0.2","-0.1","0.1","0.2","0.3"};
    int s_rhog[nglobal] = {3,25,26,27,28,29,30};  

    for(int i = 0; i < nglobal; i++)
    {
      int ifile = s_rhog[i]; 
      g_Stat[ifile]->SetMarkerStyle(20);
      g_Stat[ifile]->SetMarkerColor(color[i]);
      g_Stat[ifile]->SetLineColor(color[i]);
      g_Stat[ifile]->GetXaxis()->SetTitle("p_{T}");
      g_Stat[ifile]->GetYaxis()->SetTitle("global #rho_{00}");
      g_Stat[ifile]->GetYaxis()->SetRangeUser(0.29,0.40);
      if(i != 0)
      {
        for(int ip = 0; ip < g_Stat[ifile]->GetN(); ip++)
        {
          double pt, rho;
          g_Stat[ifile]->GetPoint(ip,pt,rho);
          g_Stat[ifile]->SetPoint(ip,pt+0.04*i,rho); 
        }
      }
      if(i == 0) g_Stat[ifile]->Draw("APE");
      else       g_Stat[ifile]->Draw("PE same");
      leg_rhog8->AddEntry(g_Stat[ifile],rhogtext[i].c_str(),"p");
    }
    leg_rhog8->Draw("same");
  }
  cout << "Study 8 " << endl;
  
  c_rho->SaveAs("figures/globalrho00_SimParVar.pdf");
}


void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE)
{
  const int nEnergy = g_rho->GetN();
  TBox *bSys[nEnergy];
  for(int i_energy = 0; i_energy < g_rho->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double err = g_rho->GetErrorYhigh(i_energy);
    
    //bSys[i_energy] = new TBox(energy/1.08,rho-err,energy*1.08,rho+err);
    bSys[i_energy] = new TBox(vmsa::pt_low[beamE][i_energy],rho-err,vmsa::pt_up[beamE][i_energy],rho+err);
    bSys[i_energy]->SetFillColor(0);
    bSys[i_energy]->SetFillStyle(0);
    bSys[i_energy]->SetLineStyle(1);
    bSys[i_energy]->SetLineWidth(1);
    bSys[i_energy]->SetLineColor(plot_color);
    bSys[i_energy]->Draw("l Same");
  }
}


void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color)
{
  const int nPt = g_rho->GetN();
  TBox *bSys[nPt];
  for(int i_pt = 0; i_pt < g_rho->GetN(); ++i_pt) // plot sys errors
  {
    double pt, rho;
    g_rho->GetPoint(i_pt,pt,rho);
    double err = g_rho->GetErrorYhigh(i_pt);

    bSys[i_pt] = new TBox(pt-0.08,rho-err,pt+0.08,rho+err);
    bSys[i_pt]->SetFillColor(0);
    bSys[i_pt]->SetFillStyle(0);
    bSys[i_pt]->SetLineStyle(1);
    bSys[i_pt]->SetLineWidth(1);
    bSys[i_pt]->SetLineColor(plot_color);
    bSys[i_pt]->Draw("l Same");
  }
}
