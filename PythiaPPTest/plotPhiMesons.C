#include "TFile.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"

#include <string>

using namespace std;

void plotPhiMesons()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TFile *file = TFile::Open("phi_mesons.root");

  const int nopts = 14;     //  1    2    3    4    5    6    7    8    9   10   11   12   13    14
  const int axis[nopts]    = {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,    0}; 
  const double  low[nopts] = {1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 1.2,  0.0};
  const double high[nopts] = {1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 4.2, 4.2, 10.0};
  string out[nopts];
  for(int iopt = 0; iopt < nopts; iopt++) out[iopt] = Form("%1.1fpT%1.1f",low[iopt],high[iopt]);

  TH2F *pty[nopts]; 
  TH2F *ptphi[nopts];  
  TH2F *yphi[nopts]; 
             
  TH1F *pt[nopts]; 
  TH1F *y[nopts];  
  TH1F *phi[nopts];

  TH1F *yratiobinned[nopts];

  TH3F *ptyphi = (TH3F*) file->Get("h_phi_ptyphi")->Clone();

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    if(axis[iopt] == 0) ptyphi->GetXaxis()->SetRangeUser(low[iopt],high[iopt]);
    if(axis[iopt] == 1) ptyphi->GetYaxis()->SetRangeUser(low[iopt],high[iopt]);
    if(axis[iopt] == 2) ptyphi->GetZaxis()->SetRangeUser(low[iopt],high[iopt]);

    string hist;
   
    hist = Form("pty_%d",iopt);
    pty[iopt]   = (TH2F*) ptyphi->Project3D("yx")->Clone(hist.c_str());

    hist = Form("ptphi_%d",iopt);
    ptphi[iopt] = (TH2F*) ptyphi->Project3D("zx")->Clone(hist.c_str());

    hist = Form("yphi_%d",iopt);
    yphi[iopt]  = (TH2F*) ptyphi->Project3D("zy")->Clone(hist.c_str());

    hist = Form("pt_%d",iopt);
    pt[iopt]  = (TH1F*) ptyphi->Project3D("x")->Clone(hist.c_str());

    hist = Form("y_%d",iopt);
    y[iopt]   = (TH1F*) ptyphi->Project3D("y")->Clone(hist.c_str());

    hist = Form("phi_%d",iopt);
    phi[iopt] = (TH1F*) ptyphi->Project3D("z")->Clone(hist.c_str()); 
  }

  TCanvas *c = new TCanvas("c","c",10,10,1200,800);
  c->Divide(3,2);
  for(int i = 0; i < 6; i++)
  {
    c->cd(i+1);
    c->cd(i+1)->SetLeftMargin(0.15);
    c->cd(i+1)->SetRightMargin(0.15);
    c->cd(i+1)->SetBottomMargin(0.15);
    c->cd(i+1)->SetGrid(0,0);
  }

  TF1* yfunc[nopts];
 
  string outputname = "figures/pp_19p6GeV_phi_mesons.pdf";
  string outputstart = Form("%s[",outputname.c_str());
  string outputstop =  Form("%s]",outputname.c_str());

  c->Print(outputstart.c_str());

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    string histname = Form("%1.1f < p_{T} < %1.1f GeV/c",low[iopt],high[iopt]);
 
    c->cd(1);
    c->cd(1)->SetLogz();
    pty[iopt]->SetTitle(histname.c_str());
    pty[iopt]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    pty[iopt]->GetYaxis()->SetTitle("y");
    pty[iopt]->Draw("Colz");

    c->cd(2);
    c->cd(2)->SetLogz();
    ptphi[iopt]->SetTitle(histname.c_str());
    ptphi[iopt]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    ptphi[iopt]->GetYaxis()->SetTitle("#phi");
    ptphi[iopt]->Draw("Colz");

    c->cd(3);
    yphi[iopt]->SetTitle(histname.c_str());
    yphi[iopt]->GetXaxis()->SetTitle("y");
    yphi[iopt]->GetYaxis()->SetTitle("#phi");
    yphi[iopt]->Draw("Colz");

    c->cd(4);
    c->cd(4)->SetLogy();
    pt[iopt]->SetTitle(histname.c_str());
    pt[iopt]->GetXaxis()->SetTitle("p_{T} GeV/c");
    pt[iopt]->GetYaxis()->SetTitle("Counts");
    pt[iopt]->Draw("PE");

    c->cd(5);

    yfunc[iopt] = new TF1(Form("yfunc%d",iopt),"[0]*(exp(-(x-[1])*(x-[1])/2/[2]/[2]))",-1,1);
    //yfunc[iopt] = new TF1(Form("yfunc%d",iopt),"[0]*(exp(-(x-[1])*(x-[1])/2/[2]/[2]) + exp(-(x-[3])*(x-[3])/2/[2]/[2]))",-5,5);
    yfunc[iopt]->SetParameter(0,y[iopt]->GetMaximum());
    yfunc[iopt]->SetParameter(2,0.75);
    //yfunc[iopt]->SetParLimits(3,-0.000001,5);
    //yfunc[iopt] = new TF1(Form("yfunc%d",iopt),"[0]*TMath::Student(x,[1])",-5,5);
    //yfunc[iopt]->SetParameter(0,y[iopt]->GetMaximum());
    //yfunc[iopt]->SetParameter(1,20);
    //yfunc[iopt]->SetParLimits(1,-5,0.000001);
    //yfunc[iopt]->SetParameter(2,0.75);
    //yfunc[iopt]->SetParLimits(3,-0.000001,5);

    y[iopt]->SetTitle(histname.c_str());
    y[iopt]->GetXaxis()->SetTitle("y");
    y[iopt]->GetYaxis()->SetTitle("Counts");
    if(iopt < nopts-2) y[iopt]->Fit(yfunc[iopt],"MRI");
    y[iopt]->Draw("PE"); 

    c->cd(6);
    phi[iopt]->SetTitle(histname.c_str());
    phi[iopt]->GetXaxis()->SetTitle("#phi");
    phi[iopt]->GetYaxis()->SetTitle("Counts");
    phi[iopt]->Draw("PE");

    //c->SaveAs(Form("figures/pp19p6_phimesons_%s.pdf",out[iopt].c_str()));
    c->Update();
    c->Print(outputname.c_str());
  }
  c->Print(outputstop.c_str());
  

  TCanvas *c2 = new TCanvas("c2","c2",10,10,1600,1600);
  c2->Divide(4,4);
  for(int i = 0; i < 16; i++)
  {
    c2->cd(i+1);
    c2->cd(i+1)->SetLeftMargin(0.15);
    c2->cd(i+1)->SetRightMargin(0.15);
    c2->cd(i+1)->SetBottomMargin(0.15);
    c2->cd(i+1)->SetGrid(0,0);
  }
  
  TH1F* ygaus_sum = (TH1F*) y[0]->Clone("ygaus_sum");
  //const int bins = ygaus_sum->GetNbinsX();
  double sum[1000000] = {0.0};

  for(int iopt = 0; iopt < nopts-2; iopt++)
  {  
    c2->cd(iopt+1);
    y[iopt]->GetXaxis()->SetRangeUser(-1.0,1.0);
    y[iopt]->Draw("PE");

    int min = y[iopt]->GetXaxis()->FindBin(-1.0);
    int max = y[iopt]->GetXaxis()->FindBin(1.0);
    double integralbin = y[iopt]->Integral(min,max);
    double binwidth = y[iopt]->GetBinWidth(1);

    double scalingfactor = integralbin/2.*binwidth;
    yratiobinned[iopt] = (TH1F*) y[iopt]->Clone(Form("yratio_%d",iopt));
    yratiobinned[iopt]->Scale(1./scalingfactor);
 
    for(int i = 1; i <= ygaus_sum->GetNbinsX(); i++)
    {
      double center = y[iopt]->GetBinCenter(i);
      sum[i-1] += yfunc[iopt]->Eval(center);
    }     
  }

  for(int i = 1; i <= ygaus_sum->GetNbinsX(); i++)
  {
    ygaus_sum->SetBinContent(i,sum[i-1]);
  }
  
  c2->cd(nopts-1);
  y[nopts-2]->Draw("PE");
  int bin_min = y[nopts-2]->GetXaxis()->FindBin(-1.0);
  int bin_max = y[nopts-2]->GetXaxis()->FindBin(1.0);
  double integral = y[nopts-2]->Integral(bin_min,bin_max);
  double binwidth = y[nopts-2]->GetBinWidth(1);
  cout << "integral from -1 y < 1 = " << integral << endl;
  cout << "We want to fill a uniform histogram with the same number of phi-meson" << endl;
  cout << "First we divide the integral by the range, which will give us a starting point for the scale of the histogram." << endl;
  cout << "      integral/2 = " << integral/2. << endl;
  cout << "Then we multiply by the bin width to calculate the actual scale of the individual histogram bins that we need." << endl;
  cout << "      integral/2*binwidth = " << integral/2.*binwidth << endl; 
  cout << "Now let's generate the plot with all bins from -1 to 1 set equal to this value. " << endl;

  TH1F *yflat = (TH1F*) y[nopts-1]->Clone("yflat");
  for(int i = 1; i <= yflat->GetNbinsX(); i++)
  {
    if( i < bin_min || i >= bin_max ) yflat->SetBinContent(i,0.0);
    else {
      yflat->SetBinContent(i,integral/2.*binwidth);
      yflat->SetBinError(i,0.0);
    }
  }
  ygaus_sum->SetLineColor(kRed);
  ygaus_sum->SetMarkerColor(kRed);
  //ygaus_sum->Draw("C same");
  
  c2->Update();
  
  double flatintegral = yflat->Integral(bin_min,bin_max);
  cout << "The integral of the flat function is = " << flatintegral << endl;

  c2->cd(nopts);
  y[nopts-2]->GetXaxis()->SetRangeUser(-1.0,1.0);
  y[nopts-2]->Draw("PE");
  yflat->GetXaxis()->SetRangeUser(-1.0,1.0);
  yflat->SetLineColor(kOrange+7);
  //yflat->Scale(integral/flatintegral);
  yflat->Draw("h same");


  double flatintegral2 = yflat->Integral(bin_min,bin_max);
  cout << "The integral of the flat function is = " << flatintegral2 << endl;
  

  
  c2->cd(nopts+1);
  TH1F *yratio = (TH1F*) y[nopts-2]->Clone("yratio");
  yratio->Divide(yflat);
  yratio->Fit(yfunc[nopts-1],"MRI","C",-1,1);
  yratio->GetYaxis()->SetTitle("Pythia/Flat");
  yratio->Draw("PE");
 

  c2->SaveAs("figures/pp_19GeV_rapiditydistributions_withFits.pdf");

  TF1* yfuncratio[nopts];

  double norm[nopts] = {0.0};
  double mean[nopts] = {0.0};
  double sigma[nopts] = {0.0};


  for(int iopt = 0; iopt < nopts-2; iopt++)
  {  
    c2->cd(iopt+1);
    yfuncratio[iopt] = new TF1(Form("yfuncratio%d",iopt),"[0]*(exp(-(x-[1])*(x-[1])/2/[2]/[2]))",-1,1);
    yfuncratio[iopt]->SetParameter(0,yratiobinned[iopt]->GetMaximum());
    yfuncratio[iopt]->SetParameter(2,0.75);
    yratiobinned[iopt]->GetXaxis()->SetRangeUser(-1.0,1.0);
    yratiobinned[iopt]->Fit(yfuncratio[iopt],"MRI","C",-1,1);
     
    norm[iopt]  = yfuncratio[iopt]->GetParameter(0);
    mean[iopt]  = yfuncratio[iopt]->GetParameter(1);
    sigma[iopt] = yfuncratio[iopt]->GetParameter(2);

    yratiobinned[iopt]->Draw("PE");
  }
  c2->SaveAs("figures/pp_19GeV_rapiditydistributionsratios_withFits.pdf");

  cout << "double norm[" << nopts-2 << "] = {";
  for(int iopt = 0; iopt < nopts-2; iopt++) 
  {
    cout << norm[iopt];
    if(iopt < nopts-3) cout << ", ";
  }
  cout << "};"<< endl;

  cout << "double mean[" << nopts-2 << "] = {";
  for(int iopt = 0; iopt < nopts-2; iopt++) 
  {
    cout << mean[iopt];
    if(iopt < nopts-3) cout << ", ";
  }
  cout << "};"<< endl;

  cout << "double sigma[" << nopts-2 << "] = {";
  for(int iopt = 0; iopt < nopts-2; iopt++) 
  {
    cout << sigma[iopt];
    if(iopt < nopts-3) cout << ", ";
  }
  cout << "};"<< endl;

}
