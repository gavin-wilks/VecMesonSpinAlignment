#include "Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "Utility/type.h"
#include <string>
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"


void calculateFPt_EPRes_EffAcc_Poly_PhiPsi(const int energy = 4, const int pid = 0, bool doall = true, bool isBesI = false, bool random3D = false, int etamode = 0, int order = 2, int yspectra = 2, int EP = 0, int v2 = 1) {


  std::string spectra;
  if(yspectra == 2) spectra = "_WithRapiditySpectra_HalfSigma";  
  if(yspectra == 1) spectra = "_WithRapiditySpectra";  
  if(yspectra == 0) spectra = "_NoRapiditySpectra";  
  if(yspectra == -1) spectra = "_NoRapiditySpectra_FixedFirstEP";  
  if(yspectra == -2) spectra = "_NoRapiditySpectra_FixedFirstEP_PhiPsi";  
  if(yspectra == -3) spectra = "_NoRapiditySpectra_PhiPsi";  

  std::string eptext = "";
  if(EP == 1) eptext = "_EP";
 
  std::string v2text = "";
  if(v2 == 0) v2text = "_noV2";

  std::string ordertext[2] = {"","2nd"};
  
  std::string etastring;
  if(etamode == 0) etastring = "eta1_eta1";
  if(etamode == 3) etastring = "eta0p4";
  if(etamode == 4) etastring = "eta0p6";
  if(etamode == 5) etastring = "eta0p8";

  //gROOT->Reset();
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(22);
  gStyle->SetLabelFont(22);
  gStyle->SetTitleFont(22);

  Double_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};
  Double_t pt_set[8] = {0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};
//  Double_t pt_set[8] = {3.0, 3.3, 3.6, 3.9, 4.3, 4.6, 4.9, 5.2};
//  Double_t pt_set[8] = {0.6, 1.4, 2.2, 3.0, 3.8, 4.6, 5.4, 7.2};

  TH1D *h_theta_star_before[6];
  TH1D *h_theta[6];
  TH1D *h_theta_before[6];
  TH1D *h_theta_star[6];
  TH1D *h_out1[6];
  TH1D *h_out2[6];

  TFile *MCFiles[6];
  double Fval[6] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Ferr[6] = {0.0};
  double Gval[6] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Gerr[6] = {0.0};
  double Fpval[6]={0.0};
  double Fperr[6]={0.0};
  double Gpval[6]={0.0};
  double Gperr[6]={0.0};
  double FfromFp[6]={0.0};
  double FfromFperr[6]={0.0};
  double FvalBESI[6] = {0.0,-0.00989006, -0.0102287, -0.0102549, -0.00800524, -0.00767692};

  TH1D *h_theta_star_clone[6];
  TF1 *Func_A[6];

  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,800);
  c1->Divide(2,2);

  /*for(int ipt = 2; ipt < 6; ipt++)
  {
    if(etamode == 0) MCFiles[ipt] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/50MEvents_OnlyEtaCut/McAcceptanceOutput_pt%d_energy%d_pid%d_cent9.root",vmsa::mPID[pid].c_str(),ipt+1,energy,pid),"READ");
    if(etamode != 0) MCFiles[ipt] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/EtaModePtStudy/McAcceptanceOutput_pt%d_energy%d_pid%d_cent9_EtaMode_%d.root",vmsa::mPID[pid].c_str(),ipt+1,energy,pid,etamode),"READ");
    //if(ipt == 5) MCFiles[ipt] = new TFile("McAcceptanceOutput_pt6_energy4_pid0_cent9_EtaMode_3_1.root","READ");
    //if(ipt == 2) MCFiles[ipt] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/McAcceptanceOutput_pt%d_energy%d_pid%d_cent9_00004999.root",vmsa::mPID[pid].c_str(),ipt+1,energy,pid),"READ");

    h_theta_before[ipt] = (TH1D*) MCFiles[ipt]->Get("h_theta_before");
    h_theta_star_before[ipt] = (TH1D*) MCFiles[ipt]->Get("h_theta_star_before");
    h_theta[ipt] = (TH1D*) MCFiles[ipt]->Get("h_theta");
    h_theta_star[ipt] = (TH1D*) MCFiles[ipt]->Get("h_theta_star");
    //h_out1[ipt] = (TH1D*) MCFiles[ipt]->Get("h_out1");
    //h_out2[ipt] = (TH1D*) MCFiles[ipt]->Get("h_out2");
 
    TF1 *Func_rho = new TF1("Func_rho","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
    //TF1 *Func_A = new TF1("Func_A","[0]*(1.+[1]*(x*x))",0,1);
    Func_A[ipt] = new TF1("Func4th",Func4th,0,1,3);
    TF1 *Func_Theta = new TF1("Func_Theta","[0]*(1.+[1]*(x*x)+[2]*(x*x*x*x))",0,1);

    //Func_rho->SetParameter(0,h_theta_star[ipt]->GetBinContent(1));
    //Func_rho->SetParameter(1,1./3.);
    TH1D *h_theta_clone = (TH1D*)h_theta[ipt]->Clone("h_theta_clone");
    h_theta_clone->Sumw2();
    h_theta_before[ipt]->Sumw2();
    h_theta_clone->Divide(h_theta_before[ipt]); 
    Func_Theta->SetParameter(0,h_theta_clone->GetBinContent(1));
    Func_Theta->SetParameter(1,0);
    h_theta_clone->Fit(Func_Theta,"ERQ0");
    //cout << "rho00 = " << Func_rho->GetParameter(1) << endl;;    


    //h_theta[ipt]->Fit(Func_Theta,"ER");
    //Fpval[ipt] = Func_Theta->GetParameter(1);
    //Fperr[ipt] = Func_Theta->GetParError(1);
    //Gpval[ipt] = Func_Theta->GetParameter(2);
    //Gperr[ipt] = Func_Theta->GetParError(2);


    h_theta_star_clone[ipt] = (TH1D*)h_theta_star[ipt]->Clone(Form("h_theta_star_clone_%d",ipt));
    h_theta_star_clone[ipt]->Sumw2();
    h_theta_star_before[ipt]->Sumw2();
    h_theta_star_clone[ipt]->Divide(h_theta_star_before[ipt]); 
    Func_A[ipt]->SetParameter(0,h_theta_star_clone[ipt]->GetBinContent(1));
    Func_A[ipt]->SetParameter(1,0);
    h_theta_star_clone[ipt]->Fit(Func_A[ipt],"ER0");

    Fval[ipt] = Func_A[ipt]->GetParameter(1);
    Ferr[ipt] = Func_A[ipt]->GetParError(1);
    Gval[ipt] = Func_A[ipt]->GetParameter(2);
    Gerr[ipt] = Func_A[ipt]->GetParError(2);

    
    //TCanvas *c_play2 = new TCanvas("c_play2","c_play2",10,10,800,800);
    //c_play2->SetLeftMargin(0.15);
    //c_play2->SetBottomMargin(0.15);
    //c_play2->SetGrid(0,0);
    //c_play2->SetTicks(1,1);
    //c_play2->cd();
    //h_theta_clone->Draw("pE");
    //Func_Theta->Draw("same");
    //c_play2->SaveAs(Form("FValues/pt%d_theta_EtaMode%d.pdf",ipt,etamode));

    //FfromFp[ipt] = -Fpval[ipt]/(2.+Fpval[ipt]);
    //FfromFperr[ipt] = TMath::Abs((Fpval[ipt]/(2.-Fpval[ipt])/(2.-Fpval[ipt])-1./(2.-Fpval[ipt]))*Fperr[ipt]);

    //TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
    cout << "ipt - 1 = " << ipt-1 << endl;
    c1->cd(ipt-1);
    c1->cd(ipt-1)->SetLeftMargin(0.15);
    c1->cd(ipt-1)->SetBottomMargin(0.15);
    c1->cd(ipt-1)->SetGrid(0,0);
    c1->cd(ipt-1)->SetTicks(1,1);
    h_theta_star_clone[ipt]->GetXaxis()->SetTitle("|cos#theta^{*}|");
    h_theta_star_clone[ipt]->SetTitle(Form("%1.1f < p_{T} < %1.1f GeV/c",vmsa::pt_low[energy][ipt],vmsa::pt_up[energy][ipt]));
    h_theta_star_clone[ipt]->Draw("pE");
    Func_A[ipt]->Draw("same");

   // delete c_play;
   // delete c_play2;
    //delete h_theta_star_clone;
    //delete h_theta_clone;
    //delete Func_rho;
    //delete Func_A;
  }*/
  //c1->cd(0);
  //c1->SaveAs(Form("FValues/pt_thetaStar_EtaMode%d.pdf",etamode));
 
  //if(isBesI)
  //{ 
  //  Fval[1] = -0.00989006;
  //  Fval[2] = -0.0102287;
  //  Fval[3] = -0.0102549;
  //  Fval[4] = -0.00800524;
  //  Fval[5] = -0.000652552;
  //}
  ////double FfromF[6] ={0.0};
  //for(int ipt = 2; ipt < 6; ipt++)
  //{
  //  //FfromF[ipt] = -Fpval[ipt]/(2.+Fpval[ipt]);
  //  cout << std::fixed << std::setprecision(7) << "pt bin: " <<  pt_set[ipt] << "-" << pt_set[ipt+1] << "GeV/c    F = " << Fval[ipt] << " +/- " << Ferr[ipt] << "    F from Theta = " << Fpval[ipt] << " +/- " << Fperr[ipt] << endl;   
  //  cout << std::fixed << std::setprecision(7) << "pt bin: " <<  pt_set[ipt] << "-" << pt_set[ipt+1] << "GeV/c    G = " << Gval[ipt] << " +/- " << Gerr[ipt] << "    G from Theta = " << Gpval[ipt] << " +/- " << Gperr[ipt] << endl;   
  //  //cout << "pt bin: " << std::fixed << std::setprecision(5) << pt_set[ipt] << "-" << pt_set[ipt+1] << "GeV/c    F = " << Fval[ipt] << " +/- " << Ferr[ipt] << "    F' = " << Fpval[ipt] << " +/- " << Fperr[ipt] << "    F from F' (swapping them)= " << -Fval[ipt]/(2.+Fval[ipt]) << " +/- " << TMath::Abs((Fval[ipt]/(2.-Fval[ipt])/(2.-Fval[ipt])-1./(2.-Fval[ipt]))*Ferr[ipt]) << endl;   
  ////out<<"-D'/(2+D'): "<<-D_theta/(2.+D_theta)<<endl;
  //}

  Char_t Char[9][10] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
  Double_t centCent[9] = {75.0,65.0,55.0,45.0,35.0,25.0,15.0,7.5,2.5};

  //double eta_D[6] = {0.0855382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};//19GeV


  TF1 *line = new TF1("line","1/3",-0.5,8.5);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
     
  TFile *fres1 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/Resolution_file_19GeV_EpdCorrections_4.root","READ");
  if(energy == 3) fres1 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/Resolution_file_14GeV_EpdCorrections_6.root","READ");
  TFile *fres2 = new TFile("../TreeProduction/UtilityFilesEta1p5OfficialCent/file_19GeV_Resolution.root","READ");
  if(energy == 3) fres2 = new TFile("../TreeProduction/StRoot/Utility/Resolution/file_14GeV_Resolution.root","READ");
  TFile  *fd12 = new TFile("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/Phi/file_19GeV_EpdFlow_4.root","READ");
  if(energy == 3) fd12 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/file_14GeV_EpdFlow_6.root","READ");
  
  TH1D *DeltaPsi1  = (TH1D*) fres1->Get("AveCosDeltaPsi1");
  TH1D *DeltaPsi2  = (TH1D*) fres2->Get("p_mRes2_Sub");
  TH1D *DeltaPsi12 = (TH1D*)  fd12->Get("p_mD12");

  TH1D *resolution_1 = new TH1D("resolution_1","resolution_1", 9, -0.5, 8.5);
  TH1D *resolution_12 = new TH1D("resolution_12","resolution_12", 9, -0.5, 8.5);
  TH1D *resolution_2 = new TH1D("resolution_2","resolution_2", 9, -0.5, 8.5);
  TH1D *resolution_2sub = new TH1D("resolution_2sub","resolution_2sub", 9, -0.5, 8.5);

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

    Double_t Cos12 = DeltaPsi12->GetBinContent(cent+1)/TMath::Sqrt(2.); // load file for D12 
    Double_t CosError12 = DeltaPsi12->GetBinError(cent+1)/TMath::Sqrt(2.);
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
    resolution_2sub->SetBinContent(cent+1,ZDCSMD_resSub_2);
    resolution_2sub->SetBinError(cent+1,ZDCSMD_resSubErr_2);

    cout<<"Cent"<<cent_set[cent+1]<<"_"<<cent_set[cent]<<": "<<ResZ<<" +/- "<<ResErrZ<<" "<<ResZ_2<<" +/- "<<ResErrZ_2<<" "<<mean<<" +/- "<<mean*TMath::Sqrt(error)<<endl;

    g_res1->SetPoint(cent,centCent[cent],ResZ);
    g_res1->SetPointError(cent,0.0,0.0,ResErrZ,ResErrZ);
    g_res2->SetPoint(cent,centCent[cent],ResZ_2);
    g_res2->SetPointError(cent,0.0,0.0,ResErrZ_2,ResErrZ_2);
    g_res12->SetPoint(cent,centCent[cent],mean);
    g_res12->SetPointError(cent,0.0,0.0,mean*TMath::Sqrt(error),mean*TMath::Sqrt(error));
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
  h_play->GetXaxis()->SetTitle("Centrality (%)");
  h_play->GetYaxis()->SetTitle("Resolution");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(0,80);
  h_play->GetYaxis()->SetRangeUser(0.0,0.8);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");
  h_play->Draw("pE");

  g_res1->SetMarkerStyle(24);
  g_res1->SetMarkerColor(kGreen+2);
  g_res1->SetMarkerSize(1.5);
  g_res1->Draw("pE Same");

  g_res2->SetMarkerStyle(20);
  g_res2->SetMarkerColor(kGreen+2);
  g_res2->SetMarkerSize(1.5);
  g_res2->Draw("pE Same");

  g_res12->SetMarkerStyle(20);
  g_res12->SetMarkerColor(kRed);
  g_res12->SetMarkerSize(1.5);
  g_res12->Draw("pE Same");

  //p_mEpdSubRes->Draw("pE Same");
  //g_mTpcSubRes2->SetMarkerStyle(24);
  //g_mTpcSubRes2->SetMarkerColor(kAzure+2);
  //g_mTpcSubRes2->SetMarkerSize(1.5);
  //g_mTpcSubRes2->Draw("pE Same");

  //g_mTpcFullRes2->SetMarkerStyle(20);
  //g_mTpcFullRes2->SetMarkerColor(kAzure+2);
  //g_mTpcFullRes2->SetMarkerSize(1.5);
  //g_mTpcFullRes2->Draw("pE Same");

  //g_mTpcSubRes3->SetMarkerStyle(24);
  //g_mTpcSubRes3->SetMarkerColor(kGray+2);
  //g_mTpcSubRes3->SetMarkerSize(1.5);
  //g_mTpcSubRes3->Draw("pE Same");

  //g_mTpcFullRes3->SetMarkerStyle(20);
  //g_mTpcFullRes3->SetMarkerColor(kGray+2);
  //g_mTpcFullRes3->SetMarkerSize(1.5);
  //g_mTpcFullRes3->Draw("pE Same");

  TLegend *leg1 = new TLegend(0.60,0.70,0.85,0.85);
  leg1->SetFillColor(10);
  leg1->SetBorderSize(0);
  leg1->AddEntry(g_res1,"R_{1} EPD","p");
  leg1->AddEntry(g_res2,"R_{2} TPC","p");
  leg1->AddEntry(g_res12,"R_{21}^{Sub} = D_{12}/(R_{1}#sqrt{2})","p");
  //leg->AddEntry(g_mTpcFullRes2,"2^{nd} Full EPD EP","p");
  //leg->AddEntry(g_mTpcSubRes2,"2^{nd} #eta_{sub} EPD EP","p");
  //leg->AddEntry(g_mTpcFullRes3,"3^{rd} RanFull EPD","p");
  //leg->AddEntry(g_mTpcSubRes3,"3^{rd} #eta_{sub}#sqrt{2} EPD","p");//leg->AddEntry(g_mZdcFullRes1,"1^{st} ZDC-SMD Full EP","p");
  //leg->AddEntry(g_mZdcFullRes2,"2^{nd} ZDC-SMD Full EP","p");
  leg1->Draw("same");

  string FigureName = Form("./figures/Resolution/Resolutions_%s.pdf",vmsa::mBeamEnergy[energy].c_str());
  c_play->SaveAs(FigureName.c_str());

  
  double Res_12 = 0;
  double Res_12_weight = 0;
  double Res_12_err = 0;

  if(!isBesI)
  {
    if(order == 1)
    {
      for(int cent=2; cent<=5; cent++) {
        Res_12 += resolution_1->GetBinContent(cent+1)/resolution_1->GetBinError(cent+1)/resolution_1->GetBinError(cent+1);
        Res_12_weight += 1./resolution_1->GetBinError(cent+1)/resolution_1->GetBinError(cent+1);
      }
      Res_12 = Res_12/Res_12_weight;
      Res_12_err = TMath::Sqrt(1./Res_12_weight); 
      cout << Res_12 << " +/- " << Res_12_err;
    }
    if(order == 2)
    {
      for(int cent=2; cent<=5; cent++) {
        Res_12 += resolution_12->GetBinContent(cent+1)/resolution_12->GetBinError(cent+1)/resolution_12->GetBinError(cent+1);
        Res_12_weight += 1./resolution_12->GetBinError(cent+1)/resolution_12->GetBinError(cent+1);
      }
      Res_12 = Res_12/Res_12_weight;
      Res_12_err = TMath::Sqrt(1./Res_12_weight); 
      cout << Res_12 << " +/- " << Res_12_err;
    } 
  }

  double R21_BESI[4] = {0.204623, 0.296309, 0.360465, 0.396092};
  double R21_BESIe[4] = {0.00703685, 0.00541934, 0.00497943, 0.00549284}; 

  if(isBesI)
  {
  Res_12 = 0.0;
  Res_12_weight = 0.0;
  for(int cent=0; cent<=3; cent++) {
    Res_12 += R21_BESI[cent]/R21_BESIe[cent]/R21_BESIe[cent];
    Res_12_weight += 1./R21_BESIe[cent]/R21_BESIe[cent];
  }
  Res_12 = Res_12/Res_12_weight;
  cout << Res_12;
  }


  if(doall)
  {
  double eff_plateau[10][7][7][7];
  double eff_plateau_error[10][7][7][7];
  double eff_eta[10][7][7][7];
  double eff_eta_error[10][7][7][7];
  //if(eff_mode==0)
  //TFile *eff_file_plateau = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/PlateauEtaStudy/Eff_%s_SingleParticle_noToF_Mode0_EtaMode%d.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
  
  TFile *eff_file_plateau;
    if(order == 1) eff_file_plateau = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order1%s%s%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),eptext.c_str(),v2text.c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
    if(order == 2) eff_file_plateau = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2%s%s%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),eptext.c_str(),v2text.c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
  //TFile *eff_file_eta = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/Eff_%s_SingleParticle_2060_fittoeta.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()),"READ");
  if(isBesI) eff_file = new TFile(Form("../../../FileTransfers/Eff_%s_SingleKaon_second.root",vmsa::mBeamEnergy[energy].c_str()),"READ");
  eff_file_plateau->Print();
  //eff_file_eta->Print();
  //else
  //  TFile *eff_file = new TFile("Eff_19GeV_SingleKaon_eta.root","READ");

  for(int i=9; i<10;i++)
  {
    for(int j=3; j<=6; j++)
    { 
      for(int k=1; k<=7; k++)
      {
        for(int l=1; l<=7; l++)
        {
          //if(j==7) {
          //  eff[i][j-1][k-1] = 1;
          //  eff_error[i][j-1][k-1] = 0;
          //  continue;
          //}

          Title = new TString("h_mEffCosEP_Cent_");
          *Title += 9;
          *Title += "_Pt_";
          *Title += (j-1);
          //Title->Print();
          TH2D *eff_hist = (TH2D*)eff_file_plateau->Get(Title->Data());
          eff_plateau[i][j-1][k-1][l-1] = eff_hist->GetBinContent(k,l);
          eff_plateau_error[i][j-1][k-1][l-1] = eff_hist->GetBinError(k,l)/eff_plateau[i][j-1][k-1][l-1];
          cout << "efficiency = " << eff_plateau[i][j-1][k-1][l-1] << " +/- " << eff_plateau_error[i][j-1][k-1][l-1] << endl;
          delete Title;
          delete eff_hist;
        }
      }
    }
  }
  eff_file_plateau->Close();

  /*for(int i=0; i<10;i++)
  {
    for(int j=1; j<=6; j++)
    { 
     for(int k=1; k<=7; k++)
      {
        //if(j==7) {
        //  eff[i][j-1][k-1] = 1;
        //  eff_error[i][j-1][k-1] = 0;
        //  continue;
        //}

        Title = new TString("h_mEffCos_Cent_");
        *Title += 9;
        *Title += "_Pt_";
        *Title += (j-1);
        TH1D *eff_hist = (TH1D*)eff_file_eta->Get(Title->Data());
        eff_eta[i][j-1][k-1] = eff_hist->GetBinContent(k);
        eff_eta_error[i][j-1][k-1] = eff_hist->GetBinError(k)/eff_eta[i][j-1][k-1];
        delete Title;
        delete eff_hist;

      }
    }
  }
  eff_file_eta->Close();
  */
  TFile *input;
  if(order == 1) input = new TFile(Form("rho00/%s/%s/Poly/PhiPsi_Raw%sPtSys_%s_PolySys_FirstOrder.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str()),"READ");
  if(order == 2) input = new TFile(Form("rho00/%s/%s/Poly/PhiPsi_Raw%sPtSys_%s_PolySys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str()),"READ");
  if(random3D) *input = new TFile(Form("rho00/%s/%s/3DRandom/Raw%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()),"READ");
  if(isBesI) input = new TFile(Form("rho00/%s/%s/BESI/Raw%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()),"READ");

  TF1 *Func_rho = new TF1("Func_rho",FuncAFG,0,1,5);
  TF1 *Func_rdl = new TF1("Func_rdl","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);

  TH1D *Big_yield = new TH1D("Big_yield","Big_yield", 7, 0, 1);
  double big_cos[9][2] = {0};

  TString *Title;

  Double_t weight_rho00 = 0;
  Double_t weight_error_rho00 = 0;
  Double_t weight_all = 0;

  Double_t pt_rho00[6] = {0};
  Double_t pt_error_rho00[6] = {0};
  Double_t pt_all[6] = {0};

  TH1D *rho00_hist[6][6][6][6][6];
  TGraphAsymmErrors *g_rho00[6][6][6][6][6];
  //TH1DMap PtCos;
  TH1F *clonePt; 
  
  string outputname;
  if(order == 1) outputname = Form("output/%s/%s/PhiPsi_AccRes%sPtSys_%s_Poly_FirstOrder%s%s%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str(),spectra.c_str(),eptext.c_str(),v2text.c_str());
  if(order == 2) outputname = Form("output/%s/%s/PhiPsi_AccRes%sPtSys_%s_Poly%s%s%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str(),spectra.c_str(),eptext.c_str(),v2text.c_str());
  if(isBesI) outputname = Form("output/%s/%s/BESI/AccRes%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *output = new TFile(outputname.c_str(),"RECREATE");
  output->cd();
  Double_t ptbin_rho00[7] = {0.4,0.8,1.2, 1.8, 2.4, 3.0, 4.2};

  for(int i_cent = 9/*vmsa::Cent_start*/; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
        for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
        {
          for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
          {
            for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
            {
              for(int i_poly = 0; i_poly < 3; i_poly++)
              {
                for(int i_F = 0; i_F < 1; i_F++)
                {
                  for(int i_eff = 0; i_eff < 1; i_eff++)
                  {
  
                    Double_t weight_rho00 = 0;
                    Double_t weight_error_rho00 = 0;
                    Double_t weight_all = 0;

                    Double_t pt_rho00[6] = {0};
                    Double_t pt_error_rho00[6] = {0};
                    Double_t pt_all[6] = {0};
         
                    
                    TCanvas *c_fit = new TCanvas("c_fit","c_fit",10,10,600,600);
                    c_fit->Divide(2,2);
                    for(int i = 0; i < 4; i++)
                    {
                      c_fit->cd(i+1)->SetLeftMargin(0.15);
                      c_fit->cd(i+1)->SetBottomMargin(0.15);
                      c_fit->cd(i+1)->SetGrid(0,0);
                      c_fit->cd(i+1)->SetTicks(1,1);
                    }
                    Title = new TString(Form("fit/Phi/%s/Order%d/pTstudy/yield_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_EtaMode%d_Divided.pdf",vmsa::mBeamEnergy[energy].c_str(),order,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1,etamode));

                    for(Int_t PtBin=3; PtBin<=6; PtBin++) 
                    {
                      string key = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",PtBin-1,i_cent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);

                      TH2F *PtCos_raw = (TH2F*)input->Get(key.c_str())->Clone("PtCos_raw");
                      TH2F *PtCos_2D = new TH2F(key.c_str(),key.c_str(), 7, 0, 1, 7, -TMath::Pi()/2.0, TMath::Pi()/2.0);
                      TH1F *PtCos;
                      //delete Title;
                      for(int i=1; i<=7; i++) 
                      {
                        for(int j=1; j<=7; j++) 
                        {
                          float inte_mean = PtCos_raw->GetBinContent(i,j);
                          float inte_mean_error = PtCos_raw->GetBinError(i,j);
                          if(i_eff == 0)
                          {
                            PtCos_2D->SetBinContent(i, j, inte_mean/eff_plateau[i_cent][PtBin-1][i-1][j-1]);
                            PtCos_2D->SetBinError(i, j, inte_mean/eff_plateau[i_cent][PtBin-1][i-1][j-1]*TMath::Sqrt(inte_mean_error*inte_mean_error/inte_mean/inte_mean+eff_plateau_error[i_cent][PtBin-1][i-1][j-1]*eff_plateau_error[i_cent][PtBin-1][i-1][j-1]));
                          }
                          //if(i_eff == 1)
                          //{
                          //  PtCos->SetBinContent(i, inte_mean/eff_eta[i_cent][PtBin-1][i-1]);
                          //  PtCos->SetBinError(i, inte_mean/eff_eta[i_cent][PtBin-1][i-1]*TMath::Sqrt(inte_mean_error*inte_mean_error/inte_mean/inte_mean+eff_eta_error[i][PtBin-1][i-1]*eff_eta_error[i][PtBin-1][i-1]));
                          //}
                        }
                      }
                      PtCos = (TH1F*) PtCos_2D->ProjectionX();

                      c_fit->cd(PtBin-2);
                      PtCos->GetXaxis()->SetTitleOffset(1.2);
                      PtCos->SetTitle(Form("%1.1f<p_{T}<%1.1fGeV/c,Cent %d-%d",vmsa::pt_low[energy][PtBin-1],vmsa::pt_up[energy][PtBin-1],20,60));
                      PtCos->GetXaxis()->SetTitle("|cos#theta*|");
                      PtCos->GetYaxis()->SetTitle("EffxAcc Corrected Yield");
                      PtCos->GetYaxis()->SetTitleOffset(1.0);
                      PtCos->SetMarkerColor(kBlack);
                      PtCos->SetMarkerSize(1.8);
                      PtCos->SetMarkerStyle(20);
                      PtCos->DrawCopy("pE");

                      //PtCos->GetXaxis()->SetTitleOffset(1.2);
                      //PtCos->GetXaxis()->SetTitle("cos#theta*");
                      //PtCos->GetYaxis()->SetTitle("yield");
                      //PtCos->GetYaxis()->SetTitleOffset(1.0);
                      //PtCos->SetMarkerColor(2);
                      //PtCos->SetMarkerSize(1.8);
                      //PtCos->SetMarkerStyle(21);
                      //PtCos->Draw();
                      //Func_rho->SetParameter(0, PtCos->GetBinContent(5));
                      //Func_rho->SetParameter(1, 0.3333);
                      //if(i_F == 0) Func_rho->FixParameter(2, Fval[PtBin-1]);
                      //if(i_F == 0) Func_rho->FixParameter(3, Gval[PtBin-1]);
                      ////if(i_F == 1) Func_rho->FixParameter(2, Fval[PtBin-1]+Ferr[PtBin-1]);
                      ////if(i_F == 2) Func_rho->FixParameter(2, Fval[PtBin-1]-Ferr[PtBin-1]);
                      ////if(i_F == 1) Func_rho->FixParameter(2, FvalBESI[PtBin-1]);
                      //Func_rho->FixParameter(4, Res_12);
                      //if(random3D) Func_rho->FixParameter(3, 0.);
                      //cout << "Resolution = " << Func_rho->GetParameter(3);
                      ////cout<<"Fit with real EP:"<<endl;
                      //PtCos->Fit(Func_rho, "NMI");
                      //Func_rho->Draw("same");

                      //Title = new TString(Form("fit/yield_pt_%d_cent_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s.pdf",PtBin-1,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str()));
           
                      //TCanvas *cy = new TCanvas();
                      //cy->SetFillColor(0);
                      //cy->SetGrid(0,0);
                      //cy->SetTitle(0);
                      //cy->SetBottomMargin(0.15);
                      //cy->SetLeftMargin(0.15);
                      //PtCos->Draw();
                      //Func_rho->Draw("same");
                      //cy->SaveAs(Title->Data());
                      //delete cy;

                      //Float_t real_rho = Func_rho->GetParameter(1);
                      //Float_t real_rho_error = Func_rho->GetParError(1);
                      
                      TF1 *Func_obs = new TF1("Func_obs","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);

                      PtCos->Fit(Func_obs,"NMI"); // fit corrected distribution for observerved rho00
                      double rho_obs = Func_obs->GetParameter(1);
                      double rho_obs_err = Func_obs->GetParError(1);
                      Func_obs->SetLineColor(kRed);      
                      Func_obs->DrawCopy("same");
       
                      // Error calculation
                      double drhodobs = 4./(1.+3.*Res_12);  // calculate d(rho)/d(rho_obs)
                      double drhodR = -12.*(rho_obs - 1./3.)/(1.+3.*Res_12)/(1.+3.*Res_12); // calculate d(rho)/d(R)

                      double real_rho = 1./3. + 4./(1.+3.*Res_12)*(rho_obs - 1./3.);
                      double real_rho_error = TMath::Sqrt((rho_obs_err*rho_obs_err)*(drhodobs*drhodobs) + (Res_12_err*Res_12_err)*(drhodR*drhodR));

                      Float_t weight = PtCos->Integral(1,7);

                      if(PtBin>=3) {
                        weight_rho00 += real_rho*weight;
                        weight_error_rho00 += real_rho_error*real_rho_error*weight*weight;
                        weight_all += weight;
                      }

                      pt_rho00[PtBin-1] += real_rho*weight;
                      pt_error_rho00[PtBin-1] += real_rho_error*real_rho_error*weight*weight;
                      pt_all[PtBin-1] += weight;

                    
                      //if(i_cent == 9 && PtBin == 3 && i_dca == 0 && i_sig == 0 && i_norm == 0 && i_sigma == 0 && i_method == 1 && i_F == 0 && i_eff == 0)
                      //{ 
                      //   string name = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_%d",PtBin-1,i_cent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F,i_eff);                             
                      //   cout << name.c_str() << endl; 
                      //   clonePt = new TH1F(name.c_str(),name.c_str(),7,0,1);
                      //   for(int i=1; i<=7; i++) {
                      //     float inte_mean = PtCos->GetBinContent(i);
                      //     float inte_mean_error = PtCos->GetBinError(i);
                      //     clonePt->SetBinContent(i, inte_mean);
                      //     clonePt->SetBinError(i, inte_mean_error);
                      //   }
                      //  
                      //}
                      delete PtCos;
                      delete PtCos_raw;
                    }
                    c_fit->SaveAs(Title->Data());

                    weight_rho00 = weight_rho00/weight_all;
                    weight_error_rho00 = TMath::Sqrt(weight_error_rho00)/weight_all;

                    for(int PtBin=3; PtBin<=6; PtBin++) {
                      pt_rho00[PtBin-1] = pt_rho00[PtBin-1]/pt_all[PtBin-1];
                      pt_error_rho00[PtBin-1] = TMath::Sqrt(pt_error_rho00[PtBin-1])/pt_all[PtBin-1];
                    }

                    cout<<"rho00"<<endl;
                    cout<<weight_rho00<<" +/- "<<weight_error_rho00<<endl;

                    for(int PtBin=3; PtBin<=6; PtBin++) {
                      cout<<pt_rho00[PtBin-1];
                      if(PtBin==6) cout<<" "<<weight_rho00<<endl;
                      else cout<<" ";
                    }

                    for(int PtBin=3; PtBin<=6; PtBin++) {
                      cout<<pt_error_rho00[PtBin-1];
                      if(PtBin==6) cout<<" "<<weight_error_rho00<<endl;
                      else cout<<" ";
                    }
 
                    Title = new TString(Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_F_%d_Eff_%d",i_cent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1,i_F,i_eff));                             

                    g_rho00[i_dca][i_sig][i_norm][i_sigma][i_method] = new TGraphAsymmErrors();
                    g_rho00[i_dca][i_sig][i_norm][i_sigma][i_method]->SetName(Title->Data());
                    for(Int_t PtBin=3; PtBin<=6; PtBin++) {
                      double ptmean = (ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                      g_rho00[i_dca][i_sig][i_norm][i_sigma][i_method]->SetPoint(PtBin-1,ptmean,pt_rho00[PtBin-1]);
                      g_rho00[i_dca][i_sig][i_norm][i_sigma][i_method]->SetPointError(PtBin-1,0.0,0.0,pt_error_rho00[PtBin-1],pt_error_rho00[PtBin-1]);
                    }
                    delete Title;

                    Title = new TString(Form("rhoFinalWeighted_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_F_%d_Eff_%d",i_cent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1,i_F,i_eff));                             

                    rho00_hist[i_dca][i_sig][i_norm][i_sigma][i_method] = new TH1D(Title->Data(), Title->Data(), 7, 0.5, 7.5);
                    rho00_hist[i_dca][i_sig][i_norm][i_sigma][i_method]->SetBinContent(7,weight_rho00);
                    rho00_hist[i_dca][i_sig][i_norm][i_sigma][i_method]->SetBinError(7,weight_error_rho00);
                    for(Int_t PtBin=3; PtBin<=6; PtBin++) {
                      rho00_hist[i_dca][i_sig][i_norm][i_sigma][i_method]->SetBinContent(PtBin,pt_rho00[PtBin-1]);
                      rho00_hist[i_dca][i_sig][i_norm][i_sigma][i_method]->SetBinError(PtBin,pt_error_rho00[PtBin-1]);
                    }

                    g_rho00[i_dca][i_sig][i_norm][i_sigma][i_method]->Write();
                    rho00_hist[i_dca][i_sig][i_norm][i_sigma][i_method]->Write();
                    cout << g_rho00[i_dca][i_sig][i_norm][i_sigma][i_method]->GetN() << endl;;
                    delete Title;
                    //delete rho00_hist;
                    //delete g_rho00;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  //clonePt->Write();
  output->Close();
  }
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

  double result = (1.+Bs*D/2.) + (As+D)*CosTheta*CosTheta + (As*D-Bs*D/2.)*CosTheta*CosTheta*CosTheta*CosTheta;

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

//double FuncAFG(double *x_val, double *par) {
//
//  double CosTheta = x_val[0];
//  double N = par[0];
//  double rho = par[1];
//  double F = par[2];
//  double G = par[2];
//  double R = par[3];
//
//  double A = (3.*rho-1.)/(1.-rho);
//  double As = A*(1.+3.*R)/(4.+A*(1.-R));
//  double Bs = A*(1.-R)/(4.+A*(1.-R));
//
//  double order0 = 16. + 8.*F + 6.*G + As*(8.+2.*F+G) - Bs*(8.+6.*F+5.*G);    
//  double order2 = 2.*(-1.*((4.+As)*F) - (6.+As)*G + Bs*(8.+9.*F+10.*G))*CosTheta*CosTheta;
//  double order4 = (-12.*Bs*F + (6.+As-25.*Bs)*G)*CosTheta*CosTheta*CosTheta*CosTheta;
//  double order6 = 10.*Bs*G*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta;
//
//  double result = order0 + order2 + order4 + order6;
//
//  return N*result;
//
//}

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

  double result = order0 + order2 + order4 + order6;

  return N*result;

}

