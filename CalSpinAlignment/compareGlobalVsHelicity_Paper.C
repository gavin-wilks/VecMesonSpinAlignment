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
#include "resolution_pt.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void compareGlobalVsHelicity_Paper(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 23, std::string simmode = "Mc", int inputpar = 0)//defaultF = 0 is BESII, defaultF = 1 is BESI
{

  const int start = 0;
  const int stop  = 5;

  const int nvar = 5;
  float sigmay[nvar] = {1./3.-0.2,1./3.-0.1,0.0,1./3.+0.1,1./3.+0.2};

  const int defaultfile = 3;

  const int nkin = 6;
  const int ptkin[nkin] = {1,1,2,2,5,5};
  const int ykin[nkin]  = {0,1,0,1,0,1};


  //const int nopts = 2;
  //std::string option = "v2_On_Off";
  //std::string date[nopts] = {"20240809g","20240809g"}; // nov2, prelimv2
  //std::string opts[nopts] = {"nov2","prelimv2"};
  //int color[nopts] = {kBlue, kOrange+7}; 
  //std::string label[nopts] = {"v_{2} OFF","v_{2} ON (Prelim)"};

  const int nopts = 11;
  std::string option = "v2_Variance";
  //std::string date[nopts] = {"20240809_2","20240815he","20240815he","20240815he","20240815he","20240815he","20240815he","20240815he","20240815he","20240815he"}; // nov2, prelimv2
  //std::string opts[nopts] = {"nov2","v20.0333","v20.0667","v20.1000","v20.1333","v20.1667","v20.2000","v20.2333","v20.2667","v20.3000"};
  float vals[nopts] = {0.0,1./30.,2./30,3./30.,4./30.,5./30,6./30.,7./30.,8./30,9./30.,10./30.};
  //int color[nvar] = {kBlue, kOrange+7, kBlack, kGray+2, kViolet, kGreen+2, kCyan+1}; 
  int color[nvar] = {kBlue, kOrange+7, kBlack, kGray+2, kViolet}; 
  //std::string label[nopts] = {"","v_{2} ON (Prelim)"};

  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  TFile *File_InPutRC[nkin]; 
  TH2D *h_m2D_RC[nopts][nvar][nkin];  
  TH2D *h_m2D_RC_Global[nopts][nvar][nkin];  
  for(int ikin = 0; ikin < nkin; ikin++)
  {
    string inputfileRC = Form("effaccfiles/%s/%s/Paper/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0_pt%d_y%d.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),ptkin[ikin],ykin[ikin]);
    File_InPutRC[ikin] = TFile::Open(inputfileRC.c_str());
    for(int iopt = 0; iopt < nopts; iopt++)
    {
      for(int i = 0; i < nvar; i++)
      {
        string KEY_2D_RC = Form("h_m%sEffCosPhiPrimeH_v2_%d_rho_%d",simmode.c_str(),iopt,i);
        string KEY_2D_RC_Global = Form("h_m%sEffCosPhiPrime_v2_%d_rho_%d",simmode.c_str(),iopt,i);
        string KEY_2D = Form("h_m%sEffCosPhiPrimeH_v2_%d_rho_%d_kin_%d",simmode.c_str(),iopt,i,ikin);
        string KEY_2D_Global = Form("h_m%sEffCosPhiPrime_v2_%d_rho_%d_%d",simmode.c_str(),iopt,i,ikin);
        cout << "Load " << KEY_2D_RC << endl;
        h_m2D_RC[iopt][i][ikin] = (TH2D*) ((TH2D*) File_InPutRC[ikin]->Get(KEY_2D_RC.c_str()))->Clone();
        cout << "Load " << KEY_2D_RC_Global << endl;
        h_m2D_RC_Global[iopt][i][ikin] = (TH2D*) ((TH2D*) File_InPutRC[ikin]->Get(KEY_2D_RC_Global.c_str()))->Clone();

        h_m2D_RC[iopt][i][ikin]->SetDirectory(0);
        h_m2D_RC_Global[iopt][i][ikin]->SetDirectory(0);

      }
    }
    File_InPutRC[ikin]->Close();
    cout << inputfileRC << endl;
  }
  TGraphAsymmErrors *g2D[nvar][nkin][5]; 
  TGraphAsymmErrors *g2D_Global[nvar][nkin][5]; 
  TF1 *fits2D[nopts][nvar][nkin];
  TF1 *fits2D_Global[nopts][nvar][nkin];

  for(int i = 0; i < nvar; i++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int par = 0; par < 5; par++) 
      {
        g2D[i][ikin][par] = new TGraphAsymmErrors();
        g2D_Global[i][ikin][par] = new TGraphAsymmErrors();
      }   
      for(int iopt = 0; iopt < nopts; iopt++)
      {
        fits2D[iopt][i][ikin] = new TF2(Form("fit2D_%d_%d_%d",iopt,i,ikin),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),6);
        fits2D[iopt][i][ikin]->SetParameter(0,1./3.);
        fits2D[iopt][i][ikin]->SetParLimits(0,0.0,1.0);
        fits2D[iopt][i][ikin]->SetParameter(5,h_m2D_RC[iopt][i][ikin]->GetMaximum());
        h_m2D_RC[iopt][i][ikin]->Fit(fits2D[iopt][i][ikin],"NMRI");
        for(int par = 0; par < 5; par++)
        {
          //if(par == 0) g2D[i][ikin][par]->SetPoint(iopt,vals[iopt]+0.002*(float(i)-float(nvar)/2.), fits2D[iopt][i][ikin]->GetParameter(0)-1./3.);
          //if(par != 0) g2D[i][ikin][par]->SetPoint(iopt,vals[iopt]+0.002*(float(i)-float(nvar)/2.), fits2D[iopt][i][ikin]->GetParameter(par));
          if(par == 0) g2D[i][ikin][par]->SetPoint(iopt,vals[iopt], fits2D[iopt][i][ikin]->GetParameter(0)-1./3.);
          if(par != 0) g2D[i][ikin][par]->SetPoint(iopt,vals[iopt], fits2D[iopt][i][ikin]->GetParameter(par));
          g2D[i][ikin][par]->SetPointError(iopt,0.0,0.0,fits2D[iopt][i][ikin]->GetParError(par),fits2D[iopt][i][ikin]->GetParError(par));
        }

        fits2D_Global[iopt][i][ikin] = new TF2(Form("fit2D_Global_%d_%d_%d",iopt,i,ikin),SpinDensity2Dcos,-1.0,1.0,0,2.0*TMath::Pi(),6);
        fits2D_Global[iopt][i][ikin]->SetParameter(0,1./3.);
        fits2D_Global[iopt][i][ikin]->SetParLimits(0,0.0,1.0);
        fits2D_Global[iopt][i][ikin]->SetParameter(5,h_m2D_RC_Global[iopt][i][ikin]->GetMaximum());
        h_m2D_RC_Global[iopt][i][ikin]->Fit(fits2D_Global[iopt][i][ikin],"NMI");
        for(int par = 0; par < 5; par++)
        {
          if(par == 0) g2D_Global[i][ikin][par]->SetPoint(iopt,vals[iopt]+0.002*(float(i)-float(nvar)/2.), fits2D_Global[iopt][i][ikin]->GetParameter(0)-1./3.);
          if(par != 0) g2D_Global[i][ikin][par]->SetPoint(iopt,vals[iopt]+0.002*(float(i)-float(nvar)/2.), fits2D_Global[iopt][i][ikin]->GetParameter(par));
          g2D_Global[i][ikin][par]->SetPointError(iopt,0.0,0.0,fits2D_Global[iopt][i][ikin]->GetParError(par),fits2D_Global[iopt][i][ikin]->GetParError(par));
        }
      }
    }
  }


  TCanvas *cfit2D = new TCanvas("cfit2D","cfit2D",10,10,2000,400);
  cfit2D->Divide(5,1);
  
  for(int i = 0; i < 5; i++)
  {
    cfit2D->cd(i+1);
    cfit2D->cd(i+1)->SetLeftMargin(0.20);  
    cfit2D->cd(i+1)->SetBottomMargin(0.15);
    cfit2D->cd(i+1)->SetTicks(1,1);
    cfit2D->cd(i+1)->SetGrid(0,0);
  }
 
  string param[5] = {"#Delta#rho_{00}","Re(#rho_{10}-#rho_{0-1})","Im(#rho_{10}-#rho_{0-1})","Re(#rho_{1-1})","Im(#rho_{1-1})"};

  std::string outputname = Form("figures/%s/%s/pTstudy/GlobalvsHelicity_PAPER_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  std::string outputstart = Form("%s[",outputname.c_str()); 
  std::string outputstop = Form("%s]",outputname.c_str()); 

  cfit2D->Print(outputstart.c_str());

  TF1 *rhoGv2[nvar][nkin];

  TLegend *leg;// = new TLegend(0.4,0.75,0.7,0.9);
  for(int ikin = 0; ikin < nkin; ikin++)
  {
    if(ikin == 0) leg = new TLegend(0.4,0.75,0.7,0.9);
    for(int par = 0; par < 5; par++)
    {
      cfit2D->cd(1+par);
      //TLegend *leg = new TLegend(0.4,0.6,0.6,0.8);
      double max = -999999;
      double min = 999999;
      for(int i = 0; i < nvar; i++)
      {
        double tmax = g2D_Global[i][ikin][par]->GetHistogram()->GetMaximum();   
        double tmin = g2D_Global[i][ikin][par]->GetHistogram()->GetMinimum(); 

        if(tmax > max) max = tmax;
        if(tmin < min) min = tmin;

        cout << "max = " << max << ", min = " << min << endl;
      }
      for(int i = 0; i < nvar; i++)
      {
        if(par == 0) 
        {
          rhoGv2[i][ikin] = new TF1(Form("rhoGv2_%d_%d",i,ikin),rhoGlobalHelicity,0.0,1./3.,3);
          rhoGv2[i][ikin]->SetParameter(0,double(ptkin[ikin]));      
          rhoGv2[i][ikin]->SetParameter(1,double(ykin[ikin]));      
          rhoGv2[i][ikin]->SetParameter(2,(double(i)-2)*0.1+1./3.);      
          rhoGv2[i][ikin]->SetLineColor(color[i]);
        }

        g2D_Global[i][ikin][par]->SetMarkerStyle(20);
        g2D_Global[i][ikin][par]->SetMarkerSize(1.2);
        g2D_Global[i][ikin][par]->SetMarkerColor(color[i]);
        g2D_Global[i][ikin][par]->SetLineColor(color[i]);
        g2D_Global[i][ikin][par]->SetTitle(Form("Input %s^{Helicity}, p_{T}<=%d GeV/c, y=%d",param[inputpar].c_str(),ptkin[ikin],ykin[ikin]));
        g2D_Global[i][ikin][par]->GetXaxis()->SetTitle(Form("v_{2} Input"));
        g2D_Global[i][ikin][par]->GetYaxis()->SetTitle(Form("%s^{Global} from 2D Fit",param[par].c_str()));
        if(min < 0.0 && max < 0.0) g2D_Global[i][ikin][par]->GetYaxis()->SetRangeUser(1.5*min,0.5*max);
        if(min < 0.0 && max > 0.0) g2D_Global[i][ikin][par]->GetYaxis()->SetRangeUser(1.5*min,1.5*max);
        if(min > 0.0 && max > 0.0) g2D_Global[i][ikin][par]->GetYaxis()->SetRangeUser(0.5*min,1.5*max);

        if(i == 0) g2D_Global[i][ikin][par]->Draw("APE"); 
        else       g2D_Global[i][ikin][par]->Draw("PE same"); 

        if(par == 0) rhoGv2[i][ikin]->Draw("l same");
          
	if(par == 0) leg->AddEntry(g2D_Global[i][ikin][par],Form("%s^{Helicity} = %1.1f",param[inputpar].c_str(),(float(i)-2)*0.1),"p");
      }
      leg->Draw("same");
    }
    
    cfit2D->Print(outputname.c_str());
    cfit2D->Update();
  }
  cfit2D->Print(outputstop.c_str());

  outputname = Form("figures/%s/%s/pTstudy/HelicityvsHelicity_PAPER_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  cfit2D->Print(outputstart.c_str());


  //TLegend *leg;// = new TLegend(0.4,0.75,0.7,0.9);
  for(int ikin = 0; ikin < nkin; ikin++)
  {
    if(ikin == 0) leg = new TLegend(0.4,0.75,0.7,0.9);
    for(int par = 0; par < 5; par++)
    {
      cfit2D->cd(1+par);
      //TLegend *leg = new TLegend(0.4,0.6,0.6,0.8);
      double max = -999999;
      double min = 999999;
      for(int i = 0; i < nvar; i++)
      {
        double tmax = g2D[i][ikin][par]->GetHistogram()->GetMaximum();   
        double tmin = g2D[i][ikin][par]->GetHistogram()->GetMinimum(); 

        if(tmax > max) max = tmax;
        if(tmin < min) min = tmin;

        cout << "max = " << max << ", min = " << min << endl;
      }
      for(int i = 0; i < nvar; i++)
      {
        g2D[i][ikin][par]->SetMarkerStyle(20);
        g2D[i][ikin][par]->SetMarkerSize(1.2);
        g2D[i][ikin][par]->SetMarkerColor(color[i]);
        g2D[i][ikin][par]->SetLineColor(color[i]);
        g2D[i][ikin][par]->SetTitle(Form("Input %s^{Helicity}, p_{T}<=%d GeV/c, y=%d",param[inputpar].c_str(),ptkin[ikin],ykin[ikin]));
        g2D[i][ikin][par]->GetXaxis()->SetTitle(Form("v_{2} Input"));
        g2D[i][ikin][par]->GetYaxis()->SetTitle(Form("%s^{Helicity} from 2D Fit",param[par].c_str()));
        if(min < 0.0 && max < 0.0) g2D[i][ikin][par]->GetYaxis()->SetRangeUser(1.5*min,0.5*max);
        if(min < 0.0 && max > 0.0) g2D[i][ikin][par]->GetYaxis()->SetRangeUser(1.5*min,1.5*max);
        if(min > 0.0 && max > 0.0) g2D[i][ikin][par]->GetYaxis()->SetRangeUser(0.5*min,1.5*max);

        if(i == 0) g2D[i][ikin][par]->Draw("APE"); 
        else       g2D[i][ikin][par]->Draw("PE same"); 
          
	//if(par == 0) leg->AddEntry(g2D[i][ikin][par],Form("%s^{Helicity} = %1.1f",param[inputpar].c_str(),(float(i)-2)*0.1),"p");
      }
      leg->Draw("same");
    }
    
    cfit2D->Print(outputname.c_str());
    cfit2D->Update();
  }
  cfit2D->Print(outputstop.c_str());

  //outputname = Form("figures/%s/%s/pTstudy/HelicityvsHelicity_PAPER_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  //outputstart = Form("%s[",outputname.c_str()); 
  //outputstop = Form("%s]",outputname.c_str()); 

  //cfit2D->Print(outputstart.c_str());

  ////TLegend *leg = new TLegend(0.4,0.75,0.7,0.9);
  //TLegend *leg2[5];// = new TLegend(0.4,0.75,0.7,0.9);
  //for(int i = 0; i < 4; i++)
  //{
  //  for(int irho = 0; irho < 5; irho++)
  //  {
  //    if(i == 0) leg2[irho] = new TLegend(0.4,0.75,0.7,0.9);
  //    for(int par = 0; par < 5; par++)
  //    {
  //      cfit2D->cd(1+par*5+irho);
  //      //TLegend *leg = new TLegend(0.4,0.6,0.6,0.8);
  //      double max = -999999;
  //      double min = 999999;
  //      for(int iopt = 0; iopt < nvar; iopt++)
  //      {
  //        double tmax = g2D[iopt][irho][par][i]->GetHistogram()->GetMaximum();   
  //        double tmin = g2D[iopt][irho][par][i]->GetHistogram()->GetMinimum(); 

  //        if(tmax > max) max = tmax;
  //        if(tmin < min) min = tmin;

  //        cout << "max = " << max << ", min = " << min << endl;
  //      }
  //      for(int iopt = 0; iopt < nvar; iopt++)
  //      {

  //        //cfit->cd(i+5)->SetLogx();
  //        //g2D[irho][par][i]->Print();
  //        g2D[iopt][irho][par][i]->SetMarkerStyle(20);
  //        g2D[iopt][irho][par][i]->SetMarkerSize(1.2);
  //        g2D[iopt][irho][par][i]->SetMarkerColor(color[iopt]);
  //        g2D[iopt][irho][par][i]->SetLineColor(color[iopt]);
  //        g2D[iopt][irho][par][i]->SetTitle(Form("Input %s^{Helicity}, %1.1f<p_{T}<%1.1f GeV/c",param[irho].c_str(),ptedges[i],ptedges[i+1]));
  //        //if(irho == 0) g2D[iopt][irho][par][i]->SetTitle(Form("Input %s^{Global} = %1.4f, %1.1f<p_{T}<%1.1f GeV/c",param[irho].c_str(),sigmay[i],ptedges[i],ptedges[i+1]));
  //        //if(irho != 0) g2D[iopt][irho][par][i]->SetTitle(Form("Input %s^{Global} = %1.1f, %1.1f<p_{T}<%1.1f GeV/c",param[irho].c_str(),sigmay_offdiag[i],ptedges[i],ptedges[i+1]));
  //        g2D[iopt][irho][par][i]->GetXaxis()->SetTitle(Form("v_{2} Input"));
  //        g2D[iopt][irho][par][i]->GetYaxis()->SetTitle(Form("%s^{Helicity} from 2D Fit",param[par].c_str()));
  //        if(min < 0.0 && max < 0.0) g2D[iopt][irho][par][i]->GetYaxis()->SetRangeUser(1.5*min,0.5*max);
  //        if(min < 0.0 && max > 0.0) g2D[iopt][irho][par][i]->GetYaxis()->SetRangeUser(1.5*min,1.5*max);
  //        if(min > 0.0 && max > 0.0) g2D[iopt][irho][par][i]->GetYaxis()->SetRangeUser(0.5*min,1.5*max);

  //        if(iopt == 0) g2D[iopt][irho][par][i]->Draw("APE"); 
  //        else          g2D[iopt][irho][par][i]->Draw("PE same"); 
  //          
  //        if(i == 0 && irho == 0 && par == 0) leg2[irho]->AddEntry(g2D[iopt][irho][par][i],Form("%s^{Helicity} = %1.4f",param[irho].c_str(),sigmay[iopt]),"p");
  //        if(i == 0 && irho != 0 && par == 0) leg2[irho]->AddEntry(g2D[iopt][irho][par][i],Form("%s^{Helicity} = %1.4f",param[irho].c_str(),sigmay_offdiag[iopt]),"p");
  //        //if(i == 0 && irho == 0 && par == 0) leg->AddEntry(g2D[iopt][irho][par][i],label[iopt].c_str(),"p");
  //      }
  //      leg2[irho]->Draw("same");
  //    }
  //  }
  //  cfit2D->Print(outputname.c_str());
  //  cfit2D->Update();
  //}
  //cfit2D->Print(outputstop.c_str());

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

