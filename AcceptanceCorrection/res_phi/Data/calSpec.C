#include <string>
#include <iostream>
#include <fstream>
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TFile.h"

using namespace std;
string const Cent[4] = {"2030","3040","4050","5060"};

double Levy(double *var, double *par)
{
  double const m0 = 1.01940; // phi-meson mass
  double pT   = var[0];
  double mT   = sqrt(pT*pT+m0*m0);
  double dNdy = par[0];
  double n    = par[1];
  double T    = par[2];

  double numer = dNdy*(n-1)*(n-2); 
  double denom = n*T*(n*T+m0*(n-2));
  double power = pow(1+(mT-m0)/(n*T),-1.0*n);

  double y = numer*power/denom;

  return y;
}

void calSpec()
{
  TGraphAsymmErrors *g_pt = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_pt_spectra[4];
  for(int i_cent = 0; i_cent < 4; ++i_cent)
  {
    g_pt_spectra[i_cent] = new TGraphAsymmErrors();
    string inputfile = Form("Phi_Spec_%s.txt",Cent[i_cent].c_str());
    cout << "inputfile = " << inputfile.c_str() << endl;
    FILE *f_spec = fopen(inputfile.c_str(),"r");
    if(f_spec == NULL)
    {
      perror("Error opening file");
    }
    else
    {
      float pT, yield, error;
      char line[80];

      int n_points = 0;
      while(fgets(line,80,f_spec))
      {
	sscanf(&line[0],"%f %f %f", &pT, &yield, &error);
	//cout << "pT = " << pT << ", yield = " << yield << ", error = " << error << endl;

	g_pt_spectra[i_cent]->SetPoint(n_points,pT,yield);
	g_pt_spectra[i_cent]->SetPointError(n_points,0.0,0.0,error,error);
	n_points++;
      }
    }
  }

  for(int i_point = 0; i_point < g_pt_spectra[0]->GetN(); ++i_point)
  {
    double pt;
    double yields = 0.0; 
    double err_yields = 0.0;
    for(int i_cent = 0; i_cent < 4; ++i_cent)
    {
      double y, err_y;
      g_pt_spectra[i_cent]->GetPoint(i_point,pt,y);
      err_y = g_pt_spectra[i_cent]->GetErrorYlow(i_point);
      yields += y;
      err_yields += err_y*err_y;
    }
    g_pt->SetPoint(i_point,pt,yields);
    g_pt->SetPointError(i_point,0.0,0.0,sqrt(err_yields),sqrt(err_yields));
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.15);
  c_play->SetGrid(0,0);
  c_play->SetTicks(1,1);
  c_play->SetLogy(1);
  TH1F *h_play = new TH1F("h_play","h_play",100,0,10.0);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1);
  } 
  h_play->SetTitle("");
  h_play->SetStats(0);

  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetXaxis()->SetNdivisions(505);
  h_play->GetXaxis()->SetRangeUser(0.0,5.0);

  h_play->GetYaxis()->SetTitle("dN/p_{T}dp_{T}dy");
  h_play->GetYaxis()->CenterTitle();
  // h_play->GetYaxis()->SetNdivisions(505);
  h_play->GetYaxis()->SetRangeUser(1E-5,10);
  h_play->GetYaxis()->SetLabelSize(0.03);
  h_play->Draw("pE");
  for(int i_cent = 0; i_cent < 4; ++i_cent)
  {
    g_pt_spectra[i_cent]->SetMarkerColor(i_cent+1);
    g_pt_spectra[i_cent]->SetMarkerSize(1.1);
    g_pt_spectra[i_cent]->SetMarkerStyle(20+i_cent);
    g_pt_spectra[i_cent]->Draw("pE same");
  }
  g_pt->SetMarkerColor(kGray+2);
  g_pt->SetMarkerSize(1.4);
  g_pt->SetMarkerStyle(24);
  g_pt->SetName("g_spec");
  g_pt->Draw("pE same");

  TF1 *f_Levy = new TF1("f_Levy",Levy,0,5,3);
  f_Levy->SetParameter(0,1);
  f_Levy->SetParameter(1,10);
  f_Levy->SetParameter(2,0.1);
  g_pt->Fit(f_Levy,"N");
  f_Levy->SetLineStyle(2);
  f_Levy->SetLineColor(2);
  f_Levy->SetLineWidth(2);
  f_Levy->Draw("l same");

  TFile *File_OutPut = new TFile("Phi_Spec.root","RECREATE");
  File_OutPut->cd();
  g_pt->Write();
  // f_Levy->Write();
  File_OutPut->Close();
}
