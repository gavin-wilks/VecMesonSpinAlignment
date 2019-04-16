#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/type.h"
#include "../StRoot/StToFMatchMaker/StToFMatchCons.h"

/*
double tof_Kaon(double* x, double* par)
{
  // return par[0] / (pow(x[0] - par[1], 2) + par[2]) - par[4] / (exp(x[0] - par[3]) + par[5]) + par[6] + par[7]*x[0];
  return par[0] / (pow(x[0] - par[1], 2) + par[2]) - par[4] / (exp(x[0] - par[3]) + par[5]) + par[6];
}
*/

double tof_Kaon(double* x, double* par)
{
  return par[0]*(1.0 / (pow(x[0] - par[1], 2) + par[2]) - par[4] / (exp(x[0] - par[3]) + par[5]) + par[6]);
}

/*
double Cheb(double *x_val, double *par)
{
  double x = x_val[0];
  double t0 = 1.0;
  double t1 = x;
  double t2 = 2.0*x*t1-1.0*t0;
  double t3 = 2.0*x*t2-1.0*t1;
  double t4 = 2.0*x*t3-1.0*t2;
  double t5 = 2.0*x*t4-1.0*t3;
  double t6 = 2.0*x*t5-1.0*t4;
  double y = par[0]*t0 + par[1]*t1 + par[2]*t2 + par[3]*t3 + par[4]*t4 + par[5]*t5 +par[6]*t6;

  return y;
}
*/

void fitToFMatchEffEnergy()
{
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);

  std::string const Energy[6] = {"11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  TFile *File_InPut[6];
  TH1F *h_mKplus[6];
  TH1F *h_mKminus[6];
  float dip_1st[6] = {0.4,0.4,0.4,0.4,0.4,0.4};
  float dip_2nd[6] = {0.8,0.8,0.8,0.8,0.8,0.8};
  float dip_3rd[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
  float dip_4th[6] = {2.5,2.5,2.5,2.5,2.5,2.5}; // Kplus
  // float dip_4th[6] = {2.0,2.5,2.3,2.5,2.5,2.5}; // Kminus
  for(int i_energy = 0; i_energy < 6; ++i_energy)
  {
    std::string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/ToFMatch/Eff_%s_ToFMatch_2060.root",Energy[i_energy].c_str(),Energy[i_energy].c_str());
    File_InPut[i_energy] = TFile::Open(inputfile.c_str());
    h_mKplus[i_energy] = (TH1F*)File_InPut[i_energy]->Get("h_mEfficiency_Kplus_Cent_9");
    h_mKplus[i_energy]->SetTitle(Energy[i_energy].c_str());
    h_mKplus[i_energy]->SetMarkerStyle(24);
    h_mKplus[i_energy]->SetMarkerSize(1.2);
    h_mKplus[i_energy]->SetMarkerColor(2);
    h_mKplus[i_energy]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_mKplus[i_energy]->GetXaxis()->CenterTitle();
    h_mKplus[i_energy]->GetYaxis()->SetTitle("Efficiency");
    h_mKplus[i_energy]->GetYaxis()->SetRangeUser(0.0,1.2);
    for(int i_bin = 0; i_bin < h_mKplus[i_energy]->GetNbinsX(); ++i_bin) // increase errors to make it inrelavant to fit
    {
      float pt = h_mKplus[i_energy]->GetBinCenter(i_bin+1); 
      // if ((pt > dip_1st[i_energy] && pt < dip_2nd[i_energy]) || (pt > dip_3rd[i_energy] && pt < dip_4th[i_energy])) h_mKplus[i_energy]->SetBinError(i_bin+1,0.5); 
      if ((pt > dip_1st[i_energy] && pt < dip_2nd[i_energy]) || (pt > dip_3rd[i_energy] && pt < dip_4th[i_energy])) h_mKplus[i_energy]->SetBinError(i_bin+1,0.0); 
      if ((pt > dip_1st[i_energy] && pt < dip_2nd[i_energy]) || (pt > dip_3rd[i_energy] && pt < dip_4th[i_energy])) h_mKplus[i_energy]->SetBinContent(i_bin+1,0.0); 
    }

    h_mKminus[i_energy] = (TH1F*)File_InPut[i_energy]->Get("h_mEfficiency_Kminus_Cent_9");
    h_mKminus[i_energy]->SetTitle(Energy[i_energy].c_str());
    h_mKminus[i_energy]->SetMarkerStyle(24);
    h_mKminus[i_energy]->SetMarkerSize(1.2);
    h_mKminus[i_energy]->SetMarkerColor(4);
    h_mKminus[i_energy]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_mKminus[i_energy]->GetXaxis()->CenterTitle();
    h_mKminus[i_energy]->GetYaxis()->SetTitle("Efficiency");
    h_mKminus[i_energy]->GetYaxis()->SetRangeUser(0.0,1.2);
    for(int i_bin = 0; i_bin < h_mKminus[i_energy]->GetNbinsX(); ++i_bin) // increase errors to make it inrelavant to fit
    {
      float pt = h_mKminus[i_energy]->GetBinCenter(i_bin+1); 
      // if ((pt > dip_1st[i_energy] && pt < dip_2nd[i_energy]) || (pt > dip_3rd[i_energy] && pt < dip_4th[i_energy])) h_mKminus[i_energy]->SetBinError(i_bin+1,0.5); 
      if ((pt > dip_1st[i_energy] && pt < dip_2nd[i_energy]) || (pt > dip_3rd[i_energy] && pt < dip_4th[i_energy])) h_mKminus[i_energy]->SetBinError(i_bin+1,0.0); 
      if ((pt > dip_1st[i_energy] && pt < dip_2nd[i_energy]) || (pt > dip_3rd[i_energy] && pt < dip_4th[i_energy])) h_mKminus[i_energy]->SetBinContent(i_bin+1,0.0); 
    }
  }

  // fit parameter from guannan
  TF1 *f_kplus[6];
  TF1 *f_kplus_plot[6];
  TF1 *f_kminus[6];
  TF1 *f_kminus_plot[6];
  int cent_kplus = 2; 
  int cent_kminus = 4; 
  double range[6] = {4.0,5.0,4.0,4.7,4.8,5.0}; // Kplus
  // double range[6] = {4.0,4.5,4.0,4.7,4.2,5.0}; // Kminus

  double par0[9] = {-0.018979, -0.0322611, -0.0680754, -0.0698575, -0.0315267, -0.00589929, -0.00226724, -0.00212137, -0.00389514 };
  double par1[9] = {0.0308943, -0.0939411, -0.14377, -0.19003, -0.116323, 0.180593, 0.207874, 0.208863, 0.194876};
  double par2[9] = {-0.00592033, -0.0600635, -0.0515391, -0.0708703, -0.0756912, 0.00912449, 0.00500487, 0.00497987, 0.00824164};
  double par3[9] = {1.28883, 1.53952, 1.52213, 1.01707, 1.5415, 2.75657, -0.326349, 1.11207, 1.37717};
  double par4[9] = {1.58923e-06, -0.00130657, -0.00973403, -0.0163526, -0.00162583, -2.20034e-05, 0.773984, 0.119933, -1.14531e-06};
  double par5[9] = {-0.340008, -0.261115, -0.236246, -0.345869, -0.260416, -0.0777638, 6.00519, 1.05048, -0.311133};
  double par6[9] = {0.685642, 0.687831, 0.682733, 0.674683, 0.659559, 0.639223, 0.658294, 0.626411, 0.578902};

  TCanvas *c_ToF_Kplus = new TCanvas("c_ToF_Kplus","c_ToF_Kplus",1200,800);
  c_ToF_Kplus->Divide(3,2);
  for(int i_energy = 0; i_energy < 6; ++i_energy)
  {
    c_ToF_Kplus->cd(i_energy+1);
    c_ToF_Kplus->cd(i_energy+1)->SetLeftMargin(0.15);
    c_ToF_Kplus->cd(i_energy+1)->SetBottomMargin(0.15);
    c_ToF_Kplus->cd(i_energy+1)->SetTicks(1,1);
    c_ToF_Kplus->cd(i_energy+1)->SetGrid(0,0);
    h_mKplus[i_energy]->Draw("pE");

    std::string FuncName = Form("f_kplus%d",i_energy);
    f_kplus[i_energy] = new TF1(FuncName.c_str(),tof_Kaon,0.20,10.0,7);
    f_kplus[i_energy]->SetParameter(0,par0[cent_kplus]);
    f_kplus[i_energy]->SetParameter(1,par1[cent_kplus]);
    f_kplus[i_energy]->SetParameter(2,par2[cent_kplus]);
    f_kplus[i_energy]->SetParameter(3,par3[cent_kplus]);
    f_kplus[i_energy]->SetParameter(4,par4[cent_kplus]);
    f_kplus[i_energy]->SetParameter(5,par5[cent_kplus]);
    f_kplus[i_energy]->SetParameter(6,par6[cent_kplus]);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_kplus[i_energy]->SetParLimits(i_par,-10.0,10.0);
    }
    f_kplus[i_energy]->SetRange(0.2,range[i_energy]);

    h_mKplus[i_energy]->Fit(f_kplus[i_energy],"MWNR");
    f_kplus[i_energy]->SetLineColor(1);
    f_kplus[i_energy]->SetLineStyle(2);
    f_kplus[i_energy]->SetLineWidth(4);
    f_kplus[i_energy]->Draw("l same");

    FuncName = Form("f_kplus_plot_%d",i_energy);
    f_kplus_plot[i_energy] = new TF1(FuncName.c_str(),tof_Kaon,0.20,10.0,7);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_kplus_plot[i_energy]->SetParameter(i_par,f_kplus[i_energy]->GetParameter(i_par));
    }
    f_kplus_plot[i_energy]->SetLineColor(2);
    f_kplus_plot[i_energy]->SetLineStyle(2);
    f_kplus_plot[i_energy]->SetLineWidth(4);
    f_kplus_plot[i_energy]->SetLineWidth(4);
    f_kplus_plot[i_energy]->Draw("l same");

    TLegend *leg = new TLegend(0.2,0.65,0.4,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->AddEntry(h_mKplus[i_energy],"K^{+}","P");
    leg->Draw("Same");
    cout << "K+ " << Energy[i_energy].c_str() << endl;
  }

  TCanvas *c_ToF_Kminus = new TCanvas("c_ToF_Kminus","c_ToF_Kminus",1200,800);
  c_ToF_Kminus->Divide(3,2);
  for(int i_energy = 0 ;i_energy < 6; ++i_energy)
  {
    c_ToF_Kminus->cd(i_energy+1);
    c_ToF_Kminus->cd(i_energy+1)->SetLeftMargin(0.15);
    c_ToF_Kminus->cd(i_energy+1)->SetBottomMargin(0.15);
    c_ToF_Kminus->cd(i_energy+1)->SetTicks(1,1);
    c_ToF_Kminus->cd(i_energy+1)->SetGrid(0,0);
    h_mKminus[i_energy]->Draw("pE");

    std::string FuncName = Form("f_kminu%d",i_energy);
    f_kminus[i_energy] = new TF1(FuncName.c_str(),tof_Kaon,0.20,10.0,7);
    f_kminus[i_energy]->SetParameter(0,par0[cent_kminus]);
    f_kminus[i_energy]->SetParameter(1,par1[cent_kminus]);
    f_kminus[i_energy]->SetParameter(2,par2[cent_kminus]);
    f_kminus[i_energy]->SetParameter(3,par3[cent_kminus]);
    f_kminus[i_energy]->SetParameter(4,par4[cent_kminus]);
    f_kminus[i_energy]->SetParameter(5,par5[cent_kminus]);
    f_kminus[i_energy]->SetParameter(6,par6[cent_kminus]);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_kminus[i_energy]->SetParLimits(i_par,-10.0,10.0);
    }
    f_kminus[i_energy]->SetRange(0.2,range[i_energy]);

    h_mKminus[i_energy]->Fit(f_kminus[i_energy],"MWNR");
    f_kminus[i_energy]->SetLineColor(1);
    f_kminus[i_energy]->SetLineStyle(2);
    f_kminus[i_energy]->SetLineWidth(4);
    f_kminus[i_energy]->Draw("l same");

    FuncName = Form("f_kminus_plot_%d",i_energy);
    f_kminus_plot[i_energy] = new TF1(FuncName.c_str(),tof_Kaon,0.20,10.0,7);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_kminus_plot[i_energy]->SetParameter(i_par,f_kminus[i_energy]->GetParameter(i_par));
    }
    f_kminus_plot[i_energy]->SetLineColor(4);
    f_kminus_plot[i_energy]->SetLineStyle(2);
    f_kminus_plot[i_energy]->SetLineWidth(4);
    f_kminus_plot[i_energy]->SetLineWidth(4);
    f_kminus_plot[i_energy]->Draw("l same");

    TLegend *leg = new TLegend(0.2,0.65,0.4,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->AddEntry(h_mKminus[i_energy],"K^{-}","P");
    leg->Draw("Same");

    cout << "K- " << Energy[i_energy].c_str() << endl;
  }

  /*
  TCanvas *c_ToF_test = new TCanvas("c_ToF_test","c_ToF_test",800,800);
  c_ToF_test->SetLeftMargin(0.15);
  c_ToF_test->SetBottomMargin(0.15);
  c_ToF_test->SetTicks(1,1);
  c_ToF_test->SetGrid(0,0);
  h_mKminus[energy]->Draw("pE");
  TF1 *f_kminus_test = new TF1("f_kminus_test",tof_Kaon,0.2,10.0,8);
  f_kminus_test->SetParameter(0,par0[cent_kminus]);
  f_kminus_test->SetParameter(1,par1[cent_kminus]);
  f_kminus_test->SetParameter(2,par2[cent_kminus]);
  f_kminus_test->SetParameter(3,par3[cent_kminus]);
  f_kminus_test->SetParameter(4,par4[cent_kminus]);
  f_kminus_test->SetParameter(5,par5[cent_kminus]);
  f_kminus_test->SetParameter(6,par6[cent_kminus]);
  // f_kminus_test->SetParameter(0,0.01);
  // f_kminus_test->SetParameter(1,0.01);
  // f_kminus_test->SetParameter(2,0.01);
  // f_kminus_test->SetParameter(3,0.01);
  // f_kminus_test->SetParameter(4,0.01);
  // f_kminus_test->SetParameter(5,0.01);
  // f_kminus_test->SetParameter(6,0.01);
  f_kminus_test->SetParameter(7,0.01);
  for(int i_par = 0; i_par < 7; ++i_par)
  {
    f_kminus_test->SetParLimits(i_par,-10.0,10.0);
  }
  f_kminus_test->SetRange(0.2,range[energy]);
  h_mKminus[energy]->Fit(f_kminus_test,"NR");
  f_kminus_test->SetLineColor(2);
  f_kminus_test->SetLineStyle(2);
  f_kminus_test->SetLineWidth(2);
  f_kminus_test->Draw("l same");
  */

}
