#include "Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"

void acceptanceQA(const int energy = 4, const int pid = 0, bool doall = true) {

  gROOT->Reset();
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(22);
  gStyle->SetLabelFont(22);
  gStyle->SetTitleFont(22);

  Double_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};
  Double_t pt_set[8] = {0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};
//  Double_t pt_set[8] = {3.0, 3.3, 3.6, 3.9, 4.3, 4.6, 4.9, 5.2};
//  Double_t pt_set[8] = {0.6, 1.4, 2.2, 3.0, 3.8, 4.6, 5.4, 7.2};

  TH1D *h_theta_star_before[5];
  TH1D *h_theta[5];
  TH1D *h_theta_star[5];
  TH1D *h_out1[5];
  TH1D *h_out2[5];

  TFile *MCFiles[5]; 
  double rho00in[5] = {0.2667,0.30,0.3333,0.3667,0.40};
  TGraphAsymmErrors *g_rho00outBefore = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rho00outAfter  = new TGraphAsymmErrors();
  double Fval[5] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Ferr[5] = {0.0};
  double Fpval[5] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Fperr[5] = {0.0};

  for(int irho = 0; irho < 5; irho++)
  {
    MCFiles[irho] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/VaryRho/McAcceptanceOutput_pt3_energy%d_pid%d_cent9_rho%d.root",vmsa::mPID[pid].c_str(),energy,pid,irho),"READ");
    //MCFiles[irho] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/VaryRho_0p8y1/McAcceptanceOutput_pt1_energy%d_pid%d_cent5_rho%d.root",vmsa::mPID[pid].c_str(),energy,pid,irho),"READ");

    h_theta_star_before[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_theta_star_before"))->Clone(Form("h_theta_star_before_%d",irho));
    h_theta[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_theta"))->Clone(Form("h_theta_%d",irho));
    h_theta_star[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_theta_star"))->Clone(Form("h_theta_star_%d",irho));
    //h_out1[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_out1"))->Clone(Form("h_out1_%d",irho));
    //h_out2[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_out2"))->Clone(Form("h_out2_%d",irho));
 
    TF1 *Func_rho = new TF1("Func_rho","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
    TF1 *Func_A = new TF1("Func_A","[0]*(1.+[1]*(x*x))",0,1);

    Func_rho->SetParameter(0,h_theta_star[irho]->GetBinContent(1));
    Func_rho->SetParameter(1,1./3.);
    h_theta_star[irho]->Fit(Func_rho,"ERQ");
    g_rho00outBefore->SetPoint(irho,rho00in[irho],Func_rho->GetParameter(1));
    g_rho00outBefore->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));

    TH1D *h_theta_clone = (TH1D*)h_theta[irho]->Clone("h_theta_clone");
    h_theta_clone->Sumw2();
    h_theta[irho]->Sumw2();
    h_theta_clone->Divide(h_theta_before[irho]);
    Func_A->SetParameter(0,h_theta_clone->GetBinContent(1));
    Func_A->SetParameter(1,0);
    h_theta_clone->Fit(Func_A,"ER");
    Fpval[irho] = Func_A->GetParameter(1);
    Fperr[irho] = Func_A->GetParError(1);

    TH1D *h_theta_star_clone = (TH1D*)h_theta_star[irho]->Clone("h_theta_star_clone");
    h_theta_star_clone->Sumw2();
    h_theta_star_before[irho]->Sumw2();
    h_theta_star_clone->Divide(h_theta_star_before[irho]); 
    Func_A->SetParameter(0,h_theta_star_clone->GetBinContent(1));
    Func_A->SetParameter(1,0);
    h_theta_star_clone->Fit(Func_A,"ER");
    Fval[irho] = Func_A->GetParameter(1);
    Ferr[irho] = Func_A->GetParError(1);
    
    
    TCanvas *c1 = new TCanvas();
    c1->SetFillColor(0);
    c1->SetGrid(0,0);
    c1->SetTitle(0);
    c1->SetBottomMargin(0.15);
    c1->SetLeftMargin(0.15);
    h_theta_star_clone->Draw("pE");
    Func_A->Draw("same"); 
    c1->SaveAs(Form("FValues/0p8y1/2ndorder_rho%d.pdf",irho));
    delete c1;

    delete h_theta_clone;
    delete h_theta_star_clone;
    delete Func_rho;
    delete Func_A;
    //MCFiles[irho]->Close();
  }
 
  double FfromFp[5] = {0.0};
  for(int irho = 0; irho < 5; irho++)
  {
    FfromFp[irho] = -Fpval[irho]/(2.+Fpval[irho]);
    cout << "F = " << Fval[irho] << " +/- " << Ferr[irho] << "     F* = " << Fpval[irho] << " +/- " << Fperr[irho] << "     F from F* = " << FfromFp[irho] << endl;
  }


  Char_t Char[9][10] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
  Double_t centCent[9] = {75.0,65.0,55.0,45.0,35.0,25.0,15.0,7.5,2.5};

  //double eta_D[6] = {0.0855382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};//19GeV


  TF1 *line = new TF1("line","1/3",-0.5,8.5);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
     
  TFile *fres1 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/Resolution_file_19GeV_EpdCorrections_4.root","READ");
  TFile *fres2 = new TFile("../TreeProduction/StRoot/Utility/Resolution/file_19GeV_Resolution.root","READ");
  TFile  *fd12 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/file_19GeV_EpdFlow.root","READ");
  //TFile event("all_event_output.root");
  
  TH1D *DeltaPsi1  = (TH1D*) fres1->Get("AveCosDeltaPsi1");
  TH1D *DeltaPsi2  = (TH1D*) fres2->Get("p_mRes2_Sub");
  TH1D *DeltaPsi12 = (TH1D*)  fd12->Get("p_mD12");

  TH1D *resolution_1 = new TH1D("resolution_1","resolution_1", 9, -0.5, 8.5);
  TH1D *resolution_12 = new TH1D("resolution_12","resolution_12", 9, -0.5, 8.5);
  TH1D *resolution_2 = new TH1D("resolution_2","resolution_2", 9, -0.5, 8.5);

  TGraphAsymmErrors *g_res1  = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_res2  = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_res12 = new TGraphAsymmErrors();

  for(Int_t cent=0; cent<9; cent++) {

    Double_t CosMean = DeltaPsi1->GetBinContent(cent+1); // load from resolution file
    Double_t CosError = DeltaPsi1->GetBinError(cent+1);
    Double_t ZDCSMD_resSub = CosMean>0.? (TMath::Sqrt(CosMean)) : 0.;
    Double_t ZDCSMD_resSubErr = CosMean>0.? CosError/2./CosMean : 0.;
    Double_t ZDCSMD_chiSub = chi(ZDCSMD_resSub);
    Double_t ZDCSMD_chiSubDelta = chi(ZDCSMD_resSub+0.005);
    Double_t ResZ = resEventPlane(TMath::Sqrt(2.)*ZDCSMD_chiSub);
    Double_t ResZDelta = resEventPlane(TMath::Sqrt(2.)*ZDCSMD_chiSubDelta);
    Double_t ResErrZ = ZDCSMD_resSubErr * TMath::Abs((ResZ-ResZDelta)/0.005);

    resolution_1->SetBinContent(cent+1, ResZ);
    resolution_1->SetBinError(cent+1, ResErrZ);

    Double_t Cos12 = DeltaPsi12->GetBinContent(cent+1); // /TMath::Sqrt(2); // load file for D12 
    Double_t CosError12 = DeltaPsi12->GetBinError(cent+1); // /TMath::Sqrt(2);
    Double_t mean = Cos12/ResZ;
    Double_t error = ResErrZ*ResErrZ/ResZ/ResZ + CosError12*CosError12/Cos12/Cos12;

    resolution_12->SetBinContent(cent+1,mean);
    resolution_12->SetBinError(cent+1, mean*TMath::Sqrt(error));

    Double_t CosMean_2 = DeltaPsi2->GetBinContent(cent+1);
    Double_t CosError_2 = DeltaPsi2->GetBinError(cent+1);
    Double_t ZDCSMD_resSub_2 = CosMean_2>0.? (TMath::Sqrt(CosMean_2)) : 0.;
    Double_t ZDCSMD_resSubErr_2 = CosMean_2>0.? CosError_2/2./CosMean_2 : 0.;
    Double_t ZDCSMD_chiSub_2 = chi(ZDCSMD_resSub_2);
    Double_t ZDCSMD_chiSubDelta_2 = chi(ZDCSMD_resSub_2+0.005);
    Double_t ResZ_2 = res1(TMath::Sqrt(2.)*ZDCSMD_chiSub_2);
    Double_t ResZDelta_2 = res1(TMath::Sqrt(2.)*ZDCSMD_chiSubDelta_2);
    Double_t ResErrZ_2 = ZDCSMD_resSubErr_2 * TMath::Abs((ResZ_2-ResZDelta_2)/0.005);

    resolution_2->SetBinContent(cent+1,ResZ_2);
    resolution_2->SetBinError(cent+1,ResErrZ_2);

    cout<<"Cent"<<cent_set[cent+1]<<"_"<<cent_set[cent]<<": "<<ResZ<<" +/- "<<ResErrZ<<" "<<ResZ_2<<" +/- "<<ResErrZ_2<<" "<<mean<<" +/- "<<mean*TMath::Sqrt(error)<<endl;

    g_res1->SetPoint(cent,centCent[cent],ResZ);
    g_res1->SetPointError(cent,0.0,0.0,ResErrZ,ResErrZ);
    g_res2->SetPoint(cent,centCent[cent],ResZ_2);
    g_res2->SetPointError(cent,0.0,0.0,ResErrZ_2,ResErrZ_2);
    g_res12->SetPoint(cent,centCent[cent],mean);
    g_res12->SetPointError(cent,0.0,0.0,mean*TMath::Sqrt(error),mean*TMath::Sqrt(error));
  }
  
  double Res_12 = 0;
  double Res_12_weight = 0;
  for(int cent=2; cent<=5; cent++) {
    Res_12 += resolution_12->GetBinContent(cent+1)/resolution_12->GetBinError(cent+1)/resolution_12->GetBinError(cent+1);
    Res_12_weight += 1./resolution_12->GetBinError(cent+1)/resolution_12->GetBinError(cent+1);
  }
  Res_12 = Res_12/Res_12_weight;


  TString *Title;

  TCanvas *c1[5];

  for(int irho = 0; irho < 5; irho++) 
  {
    TF1 *Func_rho = new TF1("Func_rho",FuncAD,0,1,4);
    //TF1 *Func_rdl = new TF1("Func_rdl","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);

    h_theta_star[irho]->GetXaxis()->SetTitleOffset(1.2);
    h_theta_star[irho]->GetXaxis()->SetTitle("cos#theta*");
    h_theta_star[irho]->GetYaxis()->SetTitle("yield");
    h_theta_star[irho]->GetYaxis()->SetTitleOffset(1.0);
    h_theta_star[irho]->SetMarkerColor(2);
    h_theta_star[irho]->SetMarkerSize(1.8);
    h_theta_star[irho]->SetMarkerStyle(21);
    //h_theta_star[irho]->Draw();
    Func_rho->SetParameter(0, h_theta_star[irho]->GetBinContent(5));
    Func_rho->SetParameter(1, 0.3);
    Func_rho->FixParameter(2, Fval[irho]);
    Func_rho->FixParameter(3, 1.0);//Res_12);

    h_theta_star[irho]->Fit(Func_rho, "NMI");

    cout << "rho after = " << Func_rho->GetParameter(1) << " +/- " << Func_rho->GetParError(1) << endl;

    TString *Title = new TString(Form("fit/acceptanceQA/yield_pt_1_cent_9_rho%d_4th.pdf",irho));
  
    TCanvas *c1 = new TCanvas(Form("c%d",irho),Form("c%d",irho),10,10,800,800);
    c1->cd();
    c1->SetFillColor(0);
    c1->SetGrid(0,0);
    c1->SetTitle(0);
    c1->SetBottomMargin(0.15);
    c1->SetLeftMargin(0.15);
    h_theta_star[irho]->Draw("pE");
    h_theta_star[irho]->Fit(Func_rho, "NMI");
    //Func_rho->Draw("same");
    c1->SaveAs(Title->Data());

    g_rho00outAfter->SetPoint(irho,rho00in[irho],Func_rho->GetParameter(1));
    g_rho00outAfter->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));
    //Func_rho->Draw("same");
    delete Func_rho;
    //delete Func_rdl;
    //Title = new TString(Form("fit/acceptanceQA/yield_pt_3_cent_9_rho%d.pdf",irho));
  
    //c1[irho] = new TCanvas(Form("c%d",irho));
    //c1[irho]->cd();
    //c1[irho]->SetFillColor(0);
    //c1[irho]->SetGrid(0,0);
    //c1[irho]->SetTitle(0);
    //c1[irho]->SetBottomMargin(0.15);
    //c1[irho]->SetLeftMargin(0.15);
    //h_theta_star[irho]->Draw("pE");
    //Func_rho->Draw("same");
    ////c1->SaveAs(Title->Data());
  
    //delete Title;
  }

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
  h_play->GetXaxis()->SetTitle("Input #rho_{00}");
  h_play->GetYaxis()->SetTitle("Output #rho_{00}");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetLimits(0.24,0.42);
  h_play->GetYaxis()->SetRangeUser(0.24,0.55);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");
  h_play->Draw("pE");

  TF1 *linear = new TF1("linear","x",0.24,0.42);
  linear->Draw("same");

  g_rho00outBefore->SetMarkerStyle(20);
  g_rho00outBefore->SetMarkerColor(kBlack);
  g_rho00outBefore->SetMarkerSize(1.5);
  g_rho00outBefore->Draw("pE Same");

  g_rho00outAfter->SetMarkerStyle(20);
  g_rho00outAfter->SetMarkerColor(kRed);
  g_rho00outAfter->SetMarkerSize(1.5);
  g_rho00outAfter->Draw("pE Same");

  TLegend *leg1 = new TLegend(0.2,0.70,0.5,0.85);
  leg1->SetFillColor(10);
  leg1->SetBorderSize(0);
  leg1->AddEntry(g_rho00outBefore,"Before Correction","p");
  leg1->AddEntry(g_rho00outAfter,"After Correction","p");
  leg1->Draw("same");

  string FigureName = Form("./figures/AcceptanceQA/AcceptanceQA_%s.pdf",vmsa::mBeamEnergy[energy].c_str());
  c_play->SaveAs(FigureName.c_str());

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

