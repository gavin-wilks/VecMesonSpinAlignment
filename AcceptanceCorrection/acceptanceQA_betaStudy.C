#include "Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"

void acceptanceQA_betaStudy(const int energy = 4, const int pid = 0, bool doall = true, std::string res = "0p4", double resval = 0.4, int cut = 1, int method = 1, int v2 = 1) {


  std::string cutoption = "";
  if(cut == 1) cutoption = "_EtaCut";

  gROOT->Reset();
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(22);
  gStyle->SetLabelFont(22);
  gStyle->SetTitleFont(22);
 
  //std::string ytextfile = "0p2y0p4";
  //std::string ytext = "0.2<y<0.4";
  std::string ytextfile = "yabs1";
  std::string ytext = "-1<y<1";

  Double_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};
  Double_t pt_set[8] = {0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};

 
  TH1F*  h_cosbeta[6][3];
  TH1F* h_cos2beta[6][3];
  TH1F* h_cos4beta[6][3];
  TH1F*  h_cosbetaP[6][3];
  TH1F* h_cos2betaP[6][3];
  TH1F* h_cos4betaP[6][3];
 
  TProfile*  p_cosbeta[6][3];
  TProfile* p_cos2beta[6][3];
  TProfile* p_cos4beta[6][3];
  TProfile*  p_cosbetaP[6][3];
  TProfile* p_cos2betaP[6][3];
  TProfile* p_cos4betaP[6][3];

  TH1F* h_beta[6][3];
  TH1F* h_tstar[6][3];

  TH1F* h_betaP[6][3];
  TH1F* h_tstarP[6][3];

  TH2F* h_betatstar[6][3];
  TH2F* h_betatstarP[6][3];
 
 
  TFile *MCFiles[3]; 
  std::string etatext[6] = {"","1.0","0.8","0.6","0.4","0.2"};
  double rho00in[3] = {0.2500,0.3333,0.4000};
  int rho00val[3] = {2500,3333,4000};

  for(int irho = 0; irho < 3; irho++)
  {
    MCFiles[irho] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/ResolutionTesting/ResolutionTest_0p4_%s_v2%d_withBeta_MANYPLOTS/McAcceptanceOutput_pt1_energy%d_pid%d_cent4_EtaMode_0_nrho%d.root",vmsa::mPID[pid].c_str(),ytextfile.c_str(),v2,energy,pid,rho00val[irho]),"READ");

    for(int ieta = 0; ieta < 6; ieta++)
    {
      std::string histname;
      histname = Form("h_cosbeta_%d",ieta);
      h_cosbeta[ieta][irho]      =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_cosbeta_%d_%d",ieta,irho));
      histname = Form("h_cos2beta_%d",ieta);
      h_cos2beta[ieta][irho]     =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_cos2beta_%d_%d",ieta,irho));
      histname = Form("h_cos4beta_%d",ieta);
      h_cos4beta[ieta][irho]     =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_cos4beta_%d_%d",ieta,irho));
      
      histname = Form("h_cosbetaP_%d",ieta);
      h_cosbetaP[ieta][irho]     =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_cosbetaP_%d_%d",ieta,irho));
      histname = Form("h_cos2betaP_%d",ieta);
      h_cos2betaP[ieta][irho]    =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_cos2betaP_%d_%d",ieta,irho));
      histname = Form("h_cos4betaP_%d",ieta);
      h_cos4betaP[ieta][irho]    =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_cos4betaP_%d_%d",ieta,irho));
      
      histname = Form("p_cosbeta_%d",ieta);
      p_cosbeta[ieta][irho]      = (TProfile*)((TProfile*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("p_cosbeta_%d_%d",ieta,irho));
      histname = Form("p_cos2beta_%d",ieta);
      p_cos2beta[ieta][irho]     = (TProfile*)((TProfile*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("p_cos2beta_%d_%d",ieta,irho));
      histname = Form("p_cos4beta_%d",ieta);
      p_cos4beta[ieta][irho]     = (TProfile*)((TProfile*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("p_cos4beta_%d_%d",ieta,irho));
      
      histname = Form("p_cosbetaP_%d",ieta);
      p_cosbetaP[ieta][irho]     = (TProfile*)((TProfile*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("p_cosbetaP_%d_%d",ieta,irho));
      histname = Form("p_cos2betaP_%d",ieta);
      p_cos2betaP[ieta][irho]    = (TProfile*)((TProfile*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("p_cos2betaP_%d_%d",ieta,irho));
      histname = Form("p_cos4betaP_%d",ieta);
      p_cos4betaP[ieta][irho]    = (TProfile*)((TProfile*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("p_cos4betaP_%d_%d",ieta,irho));
      
      histname = Form("h_beta_%d",ieta);
      h_beta[ieta][irho]         =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_beta_%d_%d",ieta,irho));
      histname = Form("h_tstar_%d",ieta);
      h_tstar[ieta][irho]        =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_tstar_%d_%d",ieta,irho));
      
      histname = Form("h_betaP_%d",ieta);
      h_betaP[ieta][irho]        =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_betaP_%d_%d",ieta,irho));
      histname = Form("h_tstarP_%d",ieta);
      h_tstarP[ieta][irho]       =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_tstarP_%d_%d",ieta,irho));
      
      histname = Form("h_betatstar_%d",ieta);
      h_betatstar[ieta][irho]    =         (TH2F*)((TH2F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_betatstar_%d_%d",ieta,irho));
      histname = Form("h_betatstarP_%d",ieta);
      h_betatstarP[ieta][irho]   =         (TH2F*)((TH2F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_betatstarP_%d_%d",ieta,irho));
 
    }
  }
 

  TH1F *h_play = new TH1F("h_play","h_play",100,0,1);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetYaxis()->SetTitle("");
  h_play->GetXaxis()->SetTitle("|cos#theta*|");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetLimits(0.0,1.0);
  h_play->GetYaxis()->SetRangeUser(-0.15,0.15);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");

  std::string FigureName;

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.15);
  c_play->SetGrid(0,0);
  c_play->SetTicks(1,1);
  c_play->cd();

  //int styleb[6] = {29,21,34,20};
  //int stylea[6] = {30,25,28,24};
  //int color[6] = {3,4,6,2};
  int styleb[6] = {20,20,34,21,29,33};
  int stylea[6] = {20,24,28,25,30,27};
  int color[6] = {1,2,6,4,3,9};

  std::string outputname = Form("./figures/AcceptanceQA/BetaStudies/cosbeta_rhocomp_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  std::string output_start = Form("%s[",outputname.c_str());
  std::string output_stop = Form("%s]",outputname.c_str());
  c_play->Print(output_start.c_str());

  // <cos(beta)> vs |cos(theta*)|
  for(int ieta = 0; ieta < 6; ieta++)
  {
    TLegend *leg2 = new TLegend(0.6,0.70,0.8,0.85);
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);

    double min = 999;
    double max = -999;

    for(int irho = 0; irho < 3; irho++)
    {
      p_cosbeta[ieta][irho]->SetMarkerStyle(styleb[irho]);
      p_cosbeta[ieta][irho]->SetMarkerColor(color[irho]);
      p_cosbeta[ieta][irho]->SetLineColor(color[irho]);
      p_cosbeta[ieta][irho]->SetMarkerSize(2.0);
      p_cosbeta[ieta][irho]->SetLineWidth(2);
      
      leg2->AddEntry(p_cosbeta[ieta][irho],Form("#rho_{00}=%1.4f",rho00in[irho]));

      double tmin = p_cosbeta[ieta][irho]->GetMinimum();
      double tmax = p_cosbeta[ieta][irho]->GetMaximum();

      if(tmin < min) min = tmin;
      if(tmax > max) max = tmax;
     
      cout << "min = " << tmin << endl;
      cout << "max = " << tmax << endl;
    }

    h_play->GetYaxis()->SetRangeUser(min-0.02,max+0.02);

    if(ieta == 0) h_play->SetTitle(Form("Before Cut, %s",ytext.c_str()));
    else          h_play->SetTitle(Form("After |#eta|<%s Cut, %s",etatext[ieta].c_str(),ytext.c_str()));
    h_play->GetYaxis()->SetTitle("#LTcos#beta#GT");
    h_play->GetXaxis()->SetTitle("|cos#theta^{*}|");
    h_play->DrawCopy("pE");
    TF1 *zero = new TF1("zero","0",0.0,1.0);
    zero->SetLineStyle(2);
    zero->Draw("same");
    for(int irho = 0; irho < 3; irho++)
    {
      p_cosbeta[ieta][irho]->Draw("pE Same");
    }
    leg2->Draw("same");

    c_play->Update();
    c_play->Print(outputname.c_str());

  }
  c_play->Print(output_stop.c_str());


  outputname = Form("./figures/AcceptanceQA/BetaStudies/cosbetaP_rhocomp_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  output_stop = Form("%s]",outputname.c_str());
  c_play->Print(output_start.c_str());

  // <cos(beta')> vs |cos(theta*')|
  for(int ieta = 0; ieta < 6; ieta++)
  {
    TLegend *leg2 = new TLegend(0.6,0.70,0.8,0.85);
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);

    double min = 999;
    double max = -999;

    for(int irho = 0; irho < 3; irho++)
    {
      p_cosbetaP[ieta][irho]->SetMarkerStyle(styleb[irho]);
      p_cosbetaP[ieta][irho]->SetMarkerColor(color[irho]);
      p_cosbetaP[ieta][irho]->SetLineColor(color[irho]);
      p_cosbetaP[ieta][irho]->SetMarkerSize(2.0);
      p_cosbetaP[ieta][irho]->SetLineWidth(2);
      
      leg2->AddEntry(p_cosbetaP[ieta][irho],Form("#rho_{00}=%1.4f",rho00in[irho]));

      double tmin = p_cosbetaP[ieta][irho]->GetMinimum();
      double tmax = p_cosbetaP[ieta][irho]->GetMaximum();

      if(tmin < min) min = tmin;
      if(tmax > max) max = tmax;
     
      cout << "min = " << tmin << endl;
      cout << "max = " << tmax << endl;
    }

    h_play->GetYaxis()->SetRangeUser(min-0.02,max+0.02);

    if(ieta == 0) h_play->SetTitle(Form("Before Cut, %s",ytext.c_str()));
    else          h_play->SetTitle(Form("After |#eta|<%s Cut, %s",etatext[ieta].c_str(),ytext.c_str()));
    h_play->GetYaxis()->SetTitle("#LTcos#beta'#GT");
    h_play->GetXaxis()->SetTitle("|cos#theta^{*}'|");
    h_play->DrawCopy("pE");
    TF1 *zero = new TF1("zero","0",0.0,1.0);
    zero->SetLineStyle(2);
    zero->Draw("same");
    for(int irho = 0; irho < 3; irho++)
    {
      p_cosbetaP[ieta][irho]->Draw("pE Same");
    }
    leg2->Draw("same");

    c_play->Update();
    c_play->Print(outputname.c_str());

  }
  c_play->Print(output_stop.c_str());

  outputname = Form("./figures/AcceptanceQA/BetaStudies/cos2beta_rhocomp_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  output_stop = Form("%s]",outputname.c_str());
  c_play->Print(output_start.c_str());

  // <cos(2beta)> vs |cos(theta*)|
  for(int ieta = 0; ieta < 6; ieta++)
  {
    TLegend *leg2 = new TLegend(0.6,0.70,0.8,0.85);
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);

    double min = 999;
    double max = -999;

    for(int irho = 0; irho < 3; irho++)
    {
      p_cos2beta[ieta][irho]->SetMarkerStyle(styleb[irho]);
      p_cos2beta[ieta][irho]->SetMarkerColor(color[irho]);
      p_cos2beta[ieta][irho]->SetLineColor(color[irho]);
      p_cos2beta[ieta][irho]->SetMarkerSize(2.0);
      p_cos2beta[ieta][irho]->SetLineWidth(2);
      
      leg2->AddEntry(p_cos2beta[ieta][irho],Form("#rho_{00}=%1.4f",rho00in[irho]));

      double tmin = p_cos2beta[ieta][irho]->GetMinimum();
      double tmax = p_cos2beta[ieta][irho]->GetMaximum();

      if(tmin < min) min = tmin;
      if(tmax > max) max = tmax;
     
      cout << "min = " << tmin << endl;
      cout << "max = " << tmax << endl;
    }

    h_play->GetYaxis()->SetRangeUser(min-0.02,max+0.02);

    if(ieta == 0) h_play->SetTitle(Form("Before Cut, %s",ytext.c_str()));
    else          h_play->SetTitle(Form("After |#eta|<%s Cut, %s",etatext[ieta].c_str(),ytext.c_str()));
    h_play->GetYaxis()->SetTitle("#LTcos2#beta#GT");
    h_play->GetXaxis()->SetTitle("|cos#theta^{*}|");
    h_play->DrawCopy("pE");
    TF1 *zero = new TF1("zero","0",0.0,1.0);
    zero->SetLineStyle(2);
    zero->Draw("same");
    for(int irho = 0; irho < 3; irho++)
    {
      p_cos2beta[ieta][irho]->Draw("pE Same");
    }
    leg2->Draw("same");

    c_play->Update();
    c_play->Print(outputname.c_str());

  }
  c_play->Print(output_stop.c_str());


  outputname = Form("./figures/AcceptanceQA/BetaStudies/cos2betaP_rhocomp_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  output_stop = Form("%s]",outputname.c_str());
  c_play->Print(output_start.c_str());

  // <cos(2beta')> vs |cos(theta*')|
  for(int ieta = 0; ieta < 6; ieta++)
  {
    TLegend *leg2 = new TLegend(0.6,0.70,0.8,0.85);
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);

    double min = 999;
    double max = -999;

    for(int irho = 0; irho < 3; irho++)
    {
      p_cos2betaP[ieta][irho]->SetMarkerStyle(styleb[irho]);
      p_cos2betaP[ieta][irho]->SetMarkerColor(color[irho]);
      p_cos2betaP[ieta][irho]->SetLineColor(color[irho]);
      p_cos2betaP[ieta][irho]->SetMarkerSize(2.0);
      p_cos2betaP[ieta][irho]->SetLineWidth(2);
      
      leg2->AddEntry(p_cos2betaP[ieta][irho],Form("#rho_{00}=%1.4f",rho00in[irho]));

      double tmin = p_cos2betaP[ieta][irho]->GetMinimum();
      double tmax = p_cos2betaP[ieta][irho]->GetMaximum();

      if(tmin < min) min = tmin;
      if(tmax > max) max = tmax;
     
      cout << "min = " << tmin << endl;
      cout << "max = " << tmax << endl;
    }

    h_play->GetYaxis()->SetRangeUser(min-0.02,max+0.02);

    if(ieta == 0) h_play->SetTitle(Form("Before Cut, %s",ytext.c_str()));
    else          h_play->SetTitle(Form("After |#eta|<%s Cut, %s",etatext[ieta].c_str(),ytext.c_str()));
    h_play->GetYaxis()->SetTitle("#LTcos2#beta'#GT");
    h_play->GetXaxis()->SetTitle("|cos#theta^{*}'|");
    h_play->DrawCopy("pE");
    TF1 *zero = new TF1("zero","0",0.0,1.0);
    zero->SetLineStyle(2);
    zero->Draw("same");
    for(int irho = 0; irho < 3; irho++)
    {
      p_cos2betaP[ieta][irho]->Draw("pE Same");
    }
    leg2->Draw("same");

    c_play->Update();
    c_play->Print(outputname.c_str());

  }
  c_play->Print(output_stop.c_str());

  outputname = Form("./figures/AcceptanceQA/BetaStudies/cos4beta_rhocomp_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  output_stop = Form("%s]",outputname.c_str());
  c_play->Print(output_start.c_str());

  // <cos(4beta)> vs |cos(theta*)|
  for(int ieta = 0; ieta < 6; ieta++)
  {
    TLegend *leg2 = new TLegend(0.6,0.70,0.8,0.85);
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);

    double min = 999;
    double max = -999;

    for(int irho = 0; irho < 3; irho++)
    {
      p_cos4beta[ieta][irho]->SetMarkerStyle(styleb[irho]);
      p_cos4beta[ieta][irho]->SetMarkerColor(color[irho]);
      p_cos4beta[ieta][irho]->SetLineColor(color[irho]);
      p_cos4beta[ieta][irho]->SetMarkerSize(2.0);
      p_cos4beta[ieta][irho]->SetLineWidth(2);
      
      leg2->AddEntry(p_cos4beta[ieta][irho],Form("#rho_{00}=%1.4f",rho00in[irho]));

      double tmin = p_cos4beta[ieta][irho]->GetMinimum();
      double tmax = p_cos4beta[ieta][irho]->GetMaximum();

      if(tmin < min) min = tmin;
      if(tmax > max) max = tmax;
     
      cout << "min = " << tmin << endl;
      cout << "max = " << tmax << endl;
    }

    h_play->GetYaxis()->SetRangeUser(min-0.02,max+0.02);

    if(ieta == 0) h_play->SetTitle(Form("Before Cut, %s",ytext.c_str()));
    else          h_play->SetTitle(Form("After |#eta|<%s Cut, %s",etatext[ieta].c_str(),ytext.c_str()));
    h_play->GetYaxis()->SetTitle("#LTcos4#beta#GT");
    h_play->GetXaxis()->SetTitle("|cos#theta^{*}|");
    h_play->DrawCopy("pE");
    TF1 *zero = new TF1("zero","0",0.0,1.0);
    zero->SetLineStyle(2);
    zero->Draw("same");
    for(int irho = 0; irho < 3; irho++)
    {
      p_cos4beta[ieta][irho]->Draw("pE Same");
    }
    leg2->Draw("same");

    c_play->Update();
    c_play->Print(outputname.c_str());

  }
  c_play->Print(output_stop.c_str());


  outputname = Form("./figures/AcceptanceQA/BetaStudies/cos4betaP_rhocomp_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  output_stop = Form("%s]",outputname.c_str());
  c_play->Print(output_start.c_str());

  // <cos(4beta')> vs |cos(theta*')|
  for(int ieta = 0; ieta < 6; ieta++)
  {
    TLegend *leg2 = new TLegend(0.6,0.70,0.8,0.85);
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);

    double min = 999;
    double max = -999;

    for(int irho = 0; irho < 3; irho++)
    {
      p_cos4betaP[ieta][irho]->SetMarkerStyle(styleb[irho]);
      p_cos4betaP[ieta][irho]->SetMarkerColor(color[irho]);
      p_cos4betaP[ieta][irho]->SetLineColor(color[irho]);
      p_cos4betaP[ieta][irho]->SetMarkerSize(2.0);
      p_cos4betaP[ieta][irho]->SetLineWidth(2);
      
      leg2->AddEntry(p_cos4betaP[ieta][irho],Form("#rho_{00}=%1.4f",rho00in[irho]));

      double tmin = p_cos4betaP[ieta][irho]->GetMinimum();
      double tmax = p_cos4betaP[ieta][irho]->GetMaximum();

      if(tmin < min) min = tmin;
      if(tmax > max) max = tmax;
     
      cout << "min = " << tmin << endl;
      cout << "max = " << tmax << endl;
    }

    h_play->GetYaxis()->SetRangeUser(min-0.02,max+0.02);

    if(ieta == 0) h_play->SetTitle(Form("Before Cut, %s",ytext.c_str()));
    else          h_play->SetTitle(Form("After |#eta|<%s Cut, %s",etatext[ieta].c_str(),ytext.c_str()));
    h_play->GetYaxis()->SetTitle("#LTcos4#beta'#GT");
    h_play->GetXaxis()->SetTitle("|cos#theta^{*}'|");
    h_play->DrawCopy("pE");
    TF1 *zero = new TF1("zero","0",0.0,1.0);
    zero->SetLineStyle(2);
    zero->Draw("same");
    for(int irho = 0; irho < 3; irho++)
    {
      p_cos4betaP[ieta][irho]->Draw("pE Same");
    }
    leg2->Draw("same");

    c_play->Update();
    c_play->Print(outputname.c_str());

  }
  c_play->Print(output_stop.c_str());

  // <cos(beta)> vs |cos(theta*)|
  for(int irho = 0; irho < 3; irho++)
  {
    TLegend *leg2 = new TLegend(0.6,0.65,0.9,0.85);
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);

    double min = 999;
    double max = -999;

    for(int ieta = 0; ieta < 6; ieta++)
    {
      p_cosbeta[ieta][irho]->SetMarkerStyle(styleb[ieta]);
      p_cosbeta[ieta][irho]->SetMarkerColor(color[ieta]);
      p_cosbeta[ieta][irho]->SetLineColor(color[ieta]);
      p_cosbeta[ieta][irho]->SetMarkerSize(2.0);
      p_cosbeta[ieta][irho]->SetLineWidth(2);

      if(ieta == 0) leg2->AddEntry(p_cosbeta[ieta][irho],Form("Before Cut"),"p");
      else leg2->AddEntry(p_cosbeta[ieta][irho],Form("After |#eta|<%s Cut",etatext[ieta].c_str()),"p");

      double tmin = p_cosbeta[ieta][irho]->GetMinimum();
      double tmax = p_cosbeta[ieta][irho]->GetMaximum();

      if(tmin < min) min = tmin;
      if(tmax > max) max = tmax;
     
      cout << "min = " << tmin << endl;
      cout << "max = " << tmax << endl;
    }

    h_play->GetYaxis()->SetRangeUser(min-0.02,max+0.02);
    //if(tmin <= 0.0 && tmax > 0.0) h_play->GetYaxis()->SetRangeUser(min*1.5,max+0.02);
    //if(tmin <= 0.0 && tmax < 0.0) h_play->GetYaxis()->SetRangeUser(min*1.5,max+0.02);
    //if(tmin >  0.0 && tmax > 0.0) h_play->GetYaxis()->SetRangeUser(min*0.5,max+0.02);

    h_play->SetTitle(Form("#rho_{00}=%1.4f, %s",rho00in[irho],ytext.c_str()));
    h_play->GetYaxis()->SetTitle("#LTcos#beta#GT");
    h_play->GetXaxis()->SetTitle("|cos#theta^{*}|");
    h_play->DrawCopy("pE");
    TF1 *zero = new TF1("zero","0",0.0,1.0);
    zero->SetLineStyle(2);
    zero->Draw("same");
    for(int ieta = 0; ieta < 6; ieta++)
    {
      p_cosbeta[ieta][irho]->Draw("pE Same");
    }
    leg2->Draw("same");

    FigureName = Form("./figures/AcceptanceQA/BetaStudies/cosbeta_%s_%s_irho%d.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str(),irho);
    c_play->SaveAs(FigureName.c_str());
  }


  // <cos(beta')> vs |cos(theta*')|
  for(int irho = 0; irho < 3; irho++)
  {
    TLegend *leg2 = new TLegend(0.6,0.65,0.9,0.85);
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);

    double min = 999;
    double max = -999;

    for(int ieta = 0; ieta < 6; ieta++)
    {
      p_cosbetaP[ieta][irho]->SetMarkerStyle(styleb[ieta]);
      p_cosbetaP[ieta][irho]->SetMarkerColor(color[ieta]);
      p_cosbetaP[ieta][irho]->SetLineColor(color[ieta]);
      p_cosbetaP[ieta][irho]->SetMarkerSize(2.0);
      p_cosbetaP[ieta][irho]->SetLineWidth(2);

      if(ieta == 0) leg2->AddEntry(p_cosbetaP[ieta][irho],Form("Before Cut"),"p");
      else leg2->AddEntry(p_cosbetaP[ieta][irho],Form("After |#eta|<%s Cut",etatext[ieta].c_str()),"p");

      double tmin = p_cosbetaP[ieta][irho]->GetMinimum();
      double tmax = p_cosbetaP[ieta][irho]->GetMaximum();

      if(tmin < min) min = tmin;
      if(tmax > max) max = tmax;
     
      cout << "min = " << tmin << endl;
      cout << "max = " << tmax << endl;
    }

    h_play->GetYaxis()->SetRangeUser(min-0.02,max+0.02);
    //if(tmin <= 0.0 && tmax > 0.0) h_play->GetYaxis()->SetRangeUser(min*1.5,max+0.02);
    //if(tmin <= 0.0 && tmax < 0.0) h_play->GetYaxis()->SetRangeUser(min*1.5,max+0.02);
    //if(tmin >  0.0 && tmax > 0.0) h_play->GetYaxis()->SetRangeUser(min*0.5,max+0.02);

    h_play->SetTitle(Form("#rho_{00}=%1.4f, %s",rho00in[irho],ytext.c_str()));
    h_play->GetYaxis()->SetTitle("#LTcos#beta'#GT");
    h_play->GetXaxis()->SetTitle("|cos#theta^{*}'|");
    h_play->DrawCopy("pE");
    TF1 *zero = new TF1("zero","0",0.0,1.0);
    zero->SetLineStyle(2);
    zero->Draw("same");
    for(int ieta = 0; ieta < 6; ieta++)
    {
      p_cosbetaP[ieta][irho]->Draw("pE Same");
    }
    leg2->Draw("same");

    FigureName = Form("./figures/AcceptanceQA/BetaStudies/cosbetaP_%s_%s_irho%d.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str(),irho);
    c_play->SaveAs(FigureName.c_str());
  }

  // <cos(2beta)> vs |cos(theta*)|
  for(int irho = 0; irho < 3; irho++)
  {
    TLegend *leg2 = new TLegend(0.6,0.65,0.9,0.85);
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);

    double min = 999;
    double max = -999;

    for(int ieta = 0; ieta < 6; ieta++)
    {
      p_cos2beta[ieta][irho]->SetMarkerStyle(styleb[ieta]);
      p_cos2beta[ieta][irho]->SetMarkerColor(color[ieta]);
      p_cos2beta[ieta][irho]->SetLineColor(color[ieta]);
      p_cos2beta[ieta][irho]->SetMarkerSize(2.0);
      p_cos2beta[ieta][irho]->SetLineWidth(2);

      if(ieta == 0) leg2->AddEntry(p_cos2beta[ieta][irho],Form("Before Cut"),"p");
      else leg2->AddEntry(p_cos2beta[ieta][irho],Form("After |#eta|<%s Cut",etatext[ieta].c_str()),"p");

      double tmin = p_cos2beta[ieta][irho]->GetMinimum();
      double tmax = p_cos2beta[ieta][irho]->GetMaximum();

      if(tmin < min) min = tmin;
      if(tmax > max) max = tmax;
     
      cout << "min = " << tmin << endl;
      cout << "max = " << tmax << endl;
    }

    h_play->GetYaxis()->SetRangeUser(min-0.02,max+0.02);
    //if(tmin <= 0.0 && tmax > 0.0) h_play->GetYaxis()->SetRangeUser(min*1.5,max+0.02);
    //if(tmin <= 0.0 && tmax < 0.0) h_play->GetYaxis()->SetRangeUser(min*1.5,max+0.02);
    //if(tmin >  0.0 && tmax > 0.0) h_play->GetYaxis()->SetRangeUser(min*0.5,max+0.02);

    h_play->SetTitle(Form("#rho_{00}=%1.4f, %s",rho00in[irho],ytext.c_str()));
    h_play->GetYaxis()->SetTitle("#LTcos2#beta#GT");
    h_play->GetXaxis()->SetTitle("|cos#theta^{*}|");
    h_play->DrawCopy("pE");
    TF1 *zero = new TF1("zero","0",0.0,1.0);
    zero->SetLineStyle(2);
    zero->Draw("same");
    for(int ieta = 0; ieta < 6; ieta++)
    {
      p_cos2beta[ieta][irho]->Draw("pE Same");
    }
    leg2->Draw("same");

    FigureName = Form("./figures/AcceptanceQA/BetaStudies/cos2beta_%s_%s_irho%d.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str(),irho);
    c_play->SaveAs(FigureName.c_str());
  }


  // <cos(2beta')> vs |cos(theta*')|
  for(int irho = 0; irho < 3; irho++)
  {
    TLegend *leg2 = new TLegend(0.6,0.65,0.9,0.85);
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);

    double min = 999;
    double max = -999;

    for(int ieta = 0; ieta < 6; ieta++)
    {
      p_cos2betaP[ieta][irho]->SetMarkerStyle(styleb[ieta]);
      p_cos2betaP[ieta][irho]->SetMarkerColor(color[ieta]);
      p_cos2betaP[ieta][irho]->SetLineColor(color[ieta]);
      p_cos2betaP[ieta][irho]->SetMarkerSize(2.0);
      p_cos2betaP[ieta][irho]->SetLineWidth(2);

      if(ieta == 0) leg2->AddEntry(p_cos2betaP[ieta][irho],Form("Before Cut"),"p");
      else leg2->AddEntry(p_cos2betaP[ieta][irho],Form("After |#eta|<%s Cut",etatext[ieta].c_str()),"p");

      double tmin = p_cos2betaP[ieta][irho]->GetMinimum();
      double tmax = p_cos2betaP[ieta][irho]->GetMaximum();

      if(tmin < min) min = tmin;
      if(tmax > max) max = tmax;
     
      cout << "min = " << tmin << endl;
      cout << "max = " << tmax << endl;
    }

    h_play->GetYaxis()->SetRangeUser(min-0.02,max+0.02);
    //if(tmin <= 0.0 && tmax > 0.0) h_play->GetYaxis()->SetRangeUser(min*1.5,max+0.02);
    //if(tmin <= 0.0 && tmax < 0.0) h_play->GetYaxis()->SetRangeUser(min*1.5,max+0.02);
    //if(tmin >  0.0 && tmax > 0.0) h_play->GetYaxis()->SetRangeUser(min*0.5,max+0.02);

    h_play->SetTitle(Form("#rho_{00}=%1.4f, %s",rho00in[irho],ytext.c_str()));
    h_play->GetYaxis()->SetTitle("#LTcos2#beta'#GT");
    h_play->GetXaxis()->SetTitle("|cos#theta^{*}'|");
    h_play->DrawCopy("pE");
    TF1 *zero = new TF1("zero","0",0.0,1.0);
    zero->SetLineStyle(2);
    zero->Draw("same");
    for(int ieta = 0; ieta < 6; ieta++)
    {
      p_cos2betaP[ieta][irho]->Draw("pE Same");
    }
    leg2->Draw("same");

    FigureName = Form("./figures/AcceptanceQA/BetaStudies/cos2betaP_%s_%s_irho%d.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str(),irho);
    c_play->SaveAs(FigureName.c_str());
  }


  // <cos(4beta)> vs |cos(theta*)|
  for(int irho = 0; irho < 3; irho++)
  {
    TLegend *leg2 = new TLegend(0.6,0.65,0.9,0.85);
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);

    double min = 999;
    double max = -999;

    for(int ieta = 0; ieta < 6; ieta++)
    {
      p_cos4beta[ieta][irho]->SetMarkerStyle(styleb[ieta]);
      p_cos4beta[ieta][irho]->SetMarkerColor(color[ieta]);
      p_cos4beta[ieta][irho]->SetLineColor(color[ieta]);
      p_cos4beta[ieta][irho]->SetMarkerSize(2.0);
      p_cos4beta[ieta][irho]->SetLineWidth(2);

      if(ieta == 0) leg2->AddEntry(p_cos4beta[ieta][irho],Form("Before Cut"),"p");
      else leg2->AddEntry(p_cos4beta[ieta][irho],Form("After |#eta|<%s Cut",etatext[ieta].c_str()),"p");

      double tmin = p_cos4beta[ieta][irho]->GetMinimum();
      double tmax = p_cos4beta[ieta][irho]->GetMaximum();

      if(tmin < min) min = tmin;
      if(tmax > max) max = tmax;
     
      cout << "min = " << tmin << endl;
      cout << "max = " << tmax << endl;
    }

    h_play->GetYaxis()->SetRangeUser(min-0.02,max+0.02);
    //if(tmin <= 0.0 && tmax > 0.0) h_play->GetYaxis()->SetRangeUser(min*1.5,max+0.02);
    //if(tmin <= 0.0 && tmax < 0.0) h_play->GetYaxis()->SetRangeUser(min*1.5,max+0.02);
    //if(tmin >  0.0 && tmax > 0.0) h_play->GetYaxis()->SetRangeUser(min*0.5,max+0.02);

    h_play->SetTitle(Form("#rho_{00}=%1.4f, %s",rho00in[irho],ytext.c_str()));
    h_play->GetYaxis()->SetTitle("#LTcos4#beta#GT");
    h_play->GetXaxis()->SetTitle("|cos#theta^{*}|");
    h_play->DrawCopy("pE");
    TF1 *zero = new TF1("zero","0",0.0,1.0);
    zero->SetLineStyle(2);
    zero->Draw("same");
    for(int ieta = 0; ieta < 6; ieta++)
    {
      p_cos4beta[ieta][irho]->Draw("pE Same");
    }
    leg2->Draw("same");

    FigureName = Form("./figures/AcceptanceQA/BetaStudies/cos4beta_%s_%s_irho%d.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str(),irho);
    c_play->SaveAs(FigureName.c_str());
  }


  // <cos(4beta')> vs |cos(theta*')|
  for(int irho = 0; irho < 3; irho++)
  {
    TLegend *leg2 = new TLegend(0.6,0.65,0.9,0.85);
    leg2->SetFillColor(10);
    leg2->SetBorderSize(0);

    double min = 999;
    double max = -999;

    for(int ieta = 0; ieta < 6; ieta++)
    {
      p_cos4betaP[ieta][irho]->SetMarkerStyle(styleb[ieta]);
      p_cos4betaP[ieta][irho]->SetMarkerColor(color[ieta]);
      p_cos4betaP[ieta][irho]->SetLineColor(color[ieta]);
      p_cos4betaP[ieta][irho]->SetMarkerSize(2.0);
      p_cos4betaP[ieta][irho]->SetLineWidth(2);

      if(ieta == 0) leg2->AddEntry(p_cos4betaP[ieta][irho],Form("Before Cut"),"p");
      else leg2->AddEntry(p_cos4betaP[ieta][irho],Form("After |#eta|<%s Cut",etatext[ieta].c_str()),"p");

      double tmin = p_cos4betaP[ieta][irho]->GetMinimum();
      double tmax = p_cos4betaP[ieta][irho]->GetMaximum();

      if(tmin < min) min = tmin;
      if(tmax > max) max = tmax;
     
      cout << "min = " << tmin << endl;
      cout << "max = " << tmax << endl;
    }

    h_play->GetYaxis()->SetRangeUser(min-0.02,max+0.02);
    //if(tmin <= 0.0 && tmax > 0.0) h_play->GetYaxis()->SetRangeUser(min*1.5,max+0.02);
    //if(tmin <= 0.0 && tmax < 0.0) h_play->GetYaxis()->SetRangeUser(min*1.5,max+0.02);
    //if(tmin >  0.0 && tmax > 0.0) h_play->GetYaxis()->SetRangeUser(min*0.5,max+0.02);

    h_play->SetTitle(Form("#rho_{00}=%1.4f, %s",rho00in[irho],ytext.c_str()));
    h_play->GetYaxis()->SetTitle("#LTcos4#beta'#GT");
    h_play->GetXaxis()->SetTitle("|cos#theta^{*}'|");
    h_play->DrawCopy("pE");
    TF1 *zero = new TF1("zero","0",0.0,1.0);
    zero->SetLineStyle(2);
    zero->Draw("same");
    for(int ieta = 0; ieta < 6; ieta++)
    {
      p_cos4betaP[ieta][irho]->Draw("pE Same");
    }
    leg2->Draw("same");

    FigureName = Form("./figures/AcceptanceQA/BetaStudies/cos4betaP_%s_%s_irho%d.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str(),irho);
    c_play->SaveAs(FigureName.c_str());
  }


  TCanvas *c_play2 = new TCanvas("c_play2","c_play2",10,10,1200,400);
  c_play2->Divide(3,1);
  for(int irho = 0; irho < 3; irho++)
  {
    c_play2->cd(irho+1)->SetLeftMargin(0.15);
    c_play2->cd(irho+1)->SetBottomMargin(0.15);
    c_play2->cd(irho+1)->SetGrid(0,0);
    c_play2->cd(irho+1)->SetTicks(1,1);
  }


  gStyle->SetOptStat("eM");

  outputname = Form("./figures/AcceptanceQA/BetaStudies/cosbeta_histograms_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());

  for(int ieta = 0; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
      c_play2->cd(irho+1);
      h_cosbeta[ieta][irho]->GetXaxis()->SetTitle("cos#beta");
      if(ieta == 0) h_cosbeta[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, Before Cut",rho00in[irho],ytext.c_str()));
      else          h_cosbeta[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, After |#eta|<%s Cut",rho00in[irho],ytext.c_str(),etatext[ieta].c_str()));
      h_cosbeta[ieta][irho]->Draw();
    }

    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }

  output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file


  outputname = Form("./figures/AcceptanceQA/BetaStudies/cosbetaP_histograms_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());

  for(int ieta = 0; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
      c_play2->cd(irho+1);
      h_cosbetaP[ieta][irho]->GetXaxis()->SetTitle("cos#beta'");
      if(ieta == 0) h_cosbetaP[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, Before Cut",rho00in[irho],ytext.c_str()));
      else          h_cosbetaP[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, After |#eta|<%s Cut",rho00in[irho],ytext.c_str(),etatext[ieta].c_str()));
      h_cosbetaP[ieta][irho]->Draw();
    }

    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }

  output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file


  outputname = Form("./figures/AcceptanceQA/BetaStudies/cos2beta_histograms_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());

  for(int ieta = 0; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
      c_play2->cd(irho+1);
      h_cos2beta[ieta][irho]->GetXaxis()->SetTitle("cos2#beta");
      if(ieta == 0) h_cos2beta[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, Before Cut",rho00in[irho],ytext.c_str()));
      else          h_cos2beta[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, After |#eta|<%s Cut",rho00in[irho],ytext.c_str(),etatext[ieta].c_str()));
      h_cos2beta[ieta][irho]->Draw();
    }

    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }

  output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file


  outputname = Form("./figures/AcceptanceQA/BetaStudies/cos2betaP_histograms_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());

  for(int ieta = 0; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
      c_play2->cd(irho+1);
      h_cos2betaP[ieta][irho]->GetXaxis()->SetTitle("cos2#beta'");
      if(ieta == 0) h_cos2betaP[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, Before Cut",rho00in[irho],ytext.c_str()));
      else          h_cos2betaP[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, After |#eta|<%s Cut",rho00in[irho],ytext.c_str(),etatext[ieta].c_str()));
      h_cos2betaP[ieta][irho]->Draw();
    }

    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }

  output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file


  outputname = Form("./figures/AcceptanceQA/BetaStudies/cos4beta_histograms_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());

  for(int ieta = 0; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
      c_play2->cd(irho+1);
      h_cos4beta[ieta][irho]->GetXaxis()->SetTitle("cos4#beta");
      if(ieta == 0) h_cos4beta[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, Before Cut",rho00in[irho],ytext.c_str()));
      else          h_cos4beta[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, After |#eta|<%s Cut",rho00in[irho],ytext.c_str(),etatext[ieta].c_str()));
      h_cos4beta[ieta][irho]->Draw();
    }

    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }

  output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file


  outputname = Form("./figures/AcceptanceQA/BetaStudies/cos4betaP_histograms_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());

  for(int ieta = 0; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
      c_play2->cd(irho+1);
      h_cos4betaP[ieta][irho]->GetXaxis()->SetTitle("cos4#beta'");
      if(ieta == 0) h_cos4betaP[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, Before Cut",rho00in[irho],ytext.c_str()));
      else          h_cos4betaP[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, After |#eta|<%s Cut",rho00in[irho],ytext.c_str(),etatext[ieta].c_str()));
      h_cos4betaP[ieta][irho]->Draw();
    }

    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }

  output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file


  outputname = Form("./figures/AcceptanceQA/BetaStudies/beta_histograms_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());

  for(int ieta = 0; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
      c_play2->cd(irho+1);
      h_beta[ieta][irho]->GetXaxis()->SetTitle("#beta");
      if(ieta == 0) h_beta[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, Before Cut",rho00in[irho],ytext.c_str()));
      else          h_beta[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, After |#eta|<%s Cut",rho00in[irho],ytext.c_str(),etatext[ieta].c_str()));
      h_beta[ieta][irho]->Draw();
    }

    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }

  output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file

  outputname = Form("./figures/AcceptanceQA/BetaStudies/betaP_histograms_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());

  for(int ieta = 0; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
      c_play2->cd(irho+1);
      h_betaP[ieta][irho]->GetXaxis()->SetTitle("#beta'");
      if(ieta == 0) h_betaP[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, Before Cut",rho00in[irho],ytext.c_str()));
      else          h_betaP[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, After |#eta|<%s Cut",rho00in[irho],ytext.c_str(),etatext[ieta].c_str()));
      h_betaP[ieta][irho]->Draw();
    }

    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }

  output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file

  outputname = Form("./figures/AcceptanceQA/BetaStudies/tstar_histograms_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());

  for(int ieta = 0; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
      c_play2->cd(irho+1);
      h_tstar[ieta][irho]->GetXaxis()->SetTitle("#theta^{*}");
      if(ieta == 0) h_tstar[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, Before Cut",rho00in[irho],ytext.c_str()));
      else          h_tstar[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, After |#eta|<%s Cut",rho00in[irho],ytext.c_str(),etatext[ieta].c_str()));
      h_tstar[ieta][irho]->Draw();
    }

    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }

  output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file

  outputname = Form("./figures/AcceptanceQA/BetaStudies/tstarP_histograms_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());

  for(int ieta = 0; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
      c_play2->cd(irho+1);
      h_tstarP[ieta][irho]->GetXaxis()->SetTitle("#theta^{*}'");
      if(ieta == 0) h_tstarP[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, Before Cut",rho00in[irho],ytext.c_str()));
      else          h_tstarP[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, After |#eta|<%s Cut",rho00in[irho],ytext.c_str(),etatext[ieta].c_str()));
      h_tstarP[ieta][irho]->Draw();
    }

    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }

  output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file

  outputname = Form("./figures/AcceptanceQA/BetaStudies/betatstar_histograms_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());

  for(int ieta = 0; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
      c_play2->cd(irho+1);
      h_betatstar[ieta][irho]->GetXaxis()->SetTitle("#theta^{*}");
      h_betatstar[ieta][irho]->GetYaxis()->SetTitle("#beta");
      if(ieta == 0) h_betatstar[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, Before Cut",rho00in[irho],ytext.c_str()));
      else          h_betatstar[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, After |#eta|<%s Cut",rho00in[irho],ytext.c_str(),etatext[ieta].c_str()));
      h_betatstar[ieta][irho]->Draw("colz");
    }

    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }

  output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file

  outputname = Form("./figures/AcceptanceQA/BetaStudies/betatstarP_histograms_%s_%s.pdf",ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str());
  output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());

  for(int ieta = 0; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
      c_play2->cd(irho+1);
      h_betatstarP[ieta][irho]->GetXaxis()->SetTitle("#theta^{*}'");
      h_betatstarP[ieta][irho]->GetYaxis()->SetTitle("#beta'");
      if(ieta == 0) h_betatstarP[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, Before Cut",rho00in[irho],ytext.c_str()));
      else          h_betatstarP[ieta][irho]->SetTitle(Form("#rho_{00}=%1.4f, %s, After |#eta|<%s Cut",rho00in[irho],ytext.c_str(),etatext[ieta].c_str()));
      h_betatstarP[ieta][irho]->Draw("colz");
    }

    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }

  output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file


}

void Correction(Double_t Res, Double_t Res_error, Double_t obv_rho, Double_t obv_rho_error, Double_t &real_rho, Double_t &real_rho_error) {

  real_rho = (Res-1+4*obv_rho)/(1+3*Res);

  Double_t real_rho_error_1 = obv_rho_error*4/(1+3*Res);
  Double_t real_rho_error_2 = Res_error*(4-12*obv_rho)/(1+3*Res)/(1+3*Res);
  real_rho_error = TMath::Sqrt(real_rho_error_1*real_rho_error_1+real_rho_error_2*real_rho_error_2);

}

Double_t chi(Double_t res) {

  Double_t chi = 2.0;
  Double_t delta = 1.0;

  for(Int_t i=0; i<15; i++) {
    chi = (res1(chi) < res) ? chi + delta : chi - delta;
    delta = delta/2.;
  }

  return chi;
}

Double_t res1(Double_t chi) {

  Double_t con = 0.626657;                   // TMath::Sqrt(pi/2)/2
  Double_t arg = chi * chi / 4.;

  Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

  return res;
}

Double_t resEventPlane(Double_t chi) {

  Double_t con = 0.626657;                   // TMath::Sqrt(pi/2)/2
  Double_t arg = chi * chi / 4.;
  Double_t halfpi = TMath::Pi()/2.;

  Double_t besselOneHalf = TMath::Sqrt(arg/halfpi) * TMath::SinH(arg)/arg;
  Double_t besselThreeHalfs = TMath::Sqrt(arg/halfpi) * (TMath::CosH(arg)/arg - TMath::SinH(arg)/(arg*arg));
  Double_t res = con * chi * exp(-arg) * (besselOneHalf + besselThreeHalfs);

  return res;
}

double FuncAD(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double D = par[2];
  double R = par[3];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double result = (1.+Bs*D/2.) + (As+D-Bs*D)*CosTheta*CosTheta + (As*D+Bs*D/2.)*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}

double FuncAF_Updated(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double R = par[3];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double result = (2. + F - Bs*F/2.) + (2.*As - F + As*F + Bs*F)*CosTheta*CosTheta - (As*F + Bs*F/2.)*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}

double Func4th(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double F = par[1];
  double G = par[2];

  double result = 1. + (4.*F+3.*G)/8. - (2.*F+3.*G)/4.*CosTheta*CosTheta + 3.*G/8.*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}

/*double Func4th(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double F = par[1];
  double G = par[2];

  double result = 2. + F + 3.*G/4. - F*CosTheta*CosTheta - 3.*G/2.*CosTheta*CosTheta + 3.*G/4.*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}*/

double FuncAFG(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double G = par[3];
  double R = par[4];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double order0 = 2. + F - Bs*F/2. + 3.*G/4. - Bs*G/2.;    
  double order2 = (2.*As - F + As*F + Bs*F - 3.*G/2. + 3.*As*G/4. + 3.*Bs*G/2.)*CosTheta*CosTheta;
  double order4 = (-1.*As*F - Bs*F/2. + 3.*G/4. - 3.*As*G/2. - 3.*Bs*G/2.)*CosTheta*CosTheta*CosTheta*CosTheta;
  double order6 = (3.*As*G/4. + Bs*G/2.)*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta;

  //double result = order0 + order2 + order4 + order6;

  return N*(order0 + order2 + order4 + order6);

}

double FuncAFGH(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double G = par[3];
  double H = par[4];
  double R = par[5];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double order0 = 64. - 16.*(-2. + Bs)*F + 8.*(3. - 2.*Bs)*G + 5.*(4. - 3.*Bs)*H;    
  double order2 = 4.*(As*(16. + 8.*F + 6.*G + 5.*H) + (-1. + Bs)*(8.*F + 12.*G + 15.*H))*CosTheta*CosTheta;
  double order4 = -2.*(8.*(2.*As + Bs)*F - 12.*G + 24.*(As + Bs)*G + 15.*(-2. + 2.*As + 3.*Bs)*H)*CosTheta*CosTheta*CosTheta*CosTheta;
  double order6 = 4.*(6.*As*G + 4.*Bs*G - 5.*H + 15.*(As + Bs)*H)*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta;
  double order8 = -5.*(4.*As + 3.*Bs)*H*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta;

  //double result = order0 + order2 + order4 + order6;

  return N*(order0 + order2 + order4 + order6 + order8);

}
