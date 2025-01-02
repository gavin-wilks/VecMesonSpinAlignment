#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "Utility/type.h"
#include <string>

using namespace std;

void calculateFCent_EPRes_EffAcc_Poly(const int energy = 4, const int pid = 0, int ptbin = 0, bool doall = false, bool isBesI = false, bool random3D = false, int mode = 0, int etamode = 0, std::string ep = "sub", int order = 1, int yspectra = 0, int EP = 1, int v2 = 0) {

  std::string spectra;
  if(yspectra == 2) spectra = "_WithRapiditySpectra_HalfSigma";  
  if(yspectra == 1) spectra = "_WithRapiditySpectra";  
  if(yspectra == 0) spectra = "_NoRapiditySpectra";  
  if(yspectra == -1) spectra = "_NoRapiditySpectra_FixedFirstEP";  

  std::string eptext = "";
  if(EP == 1) eptext = "_EP";
 
  std::string v2text = "";
  if(v2 == 0) v2text = "_noV2";

  int centVal[10] = {80,70,60,50,40,30,20,10,5,0};
 
  std::string ordertext[2] = {"1st","2nd"};
 
  std::string xtraopt = "_noDelta";
 
  std::string etastring;
  if(etamode == 0) etastring = "eta1_eta1";
  if(etamode == 3) etastring = "eta0p4";
  if(etamode == 4) etastring = "eta0p6";
  if(etamode == 5) etastring = "eta0p8";

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

  TH1D *h_theta_star_ratio[10][10];
  TH1D *h_theta_ratio[10][10];
  //TH1D *h_theta_star_before[6];
  //TH1D *h_theta[6];
  //TH1D *h_theta_star[6];
  //TH1D *h_out1[6];
  //TH1D *h_out2[6];

  TFile *MCFile;
  double Fval[10][10] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Ferr[10][10] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Gval[10][10] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Gerr[10][10] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};

  //MCFile = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/AcceptanceProccessTupleEPRes_%s%s/Acceptance_%s_19GeV_Mode%d_EtaMode%d.root",ep.c_str(),xtraopt.c_str(),vmsa::mPID[pid].c_str(),mode,etamode),"READ");

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,900,900);
  c_play->Divide(3,3);
  for(int i = 0; i < 9; i++)
  {
    //c_play->cd(i+1)->SetTopMargin(0.10);
    c_play->cd(i+1)->SetLeftMargin(0.15);
    c_play->cd(i+1)->SetBottomMargin(0.15);
    c_play->cd(i+1)->SetGrid(0,0);
    c_play->cd(i+1)->SetTicks(1,1);
  }

  //string outputname = Form("figures/Fplots_%s_%s_Mode%d_EtaMode%d_rapidity.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),mode,etamode);
  //string output_start = Form("%s[",outputname.c_str());
 // c_play->Print(output_start.c_str());
  //cout << "c_play->Print()" << endl;
  
  //for(int ipt = vmsa::pt_rebin_first_y[energy]; ipt <= vmsa::pt_rebin_last_y[energy]; ipt++)
  //{
  //  cout << "ipt = " << ipt << endl;
  //  for(int icent = 0; icent < vmsa::cent_rebin_total; icent++)
  //  {
  //    cout << "icent = " << icent << endl;
  //    for(int iy = 0; iy < vmsa::y_total; iy++)
  //    {
  //      c_play->cd(iy+1);
  //      //if(ipt == 2) MCFiles[ipt] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/McAcceptanceOutput_pt%d_energy%d_pid%d_cent9_00004999.root",vmsa::mPID[pid].c_str(),ipt+1,energy,pid),"READ");

  //      //h_theta_star_before[ipt] = (TH1D*) MCFiles->Get("h_theta_star_before");
  //      h_theta_star_ratio[ipt][icent][iy] = (TH1D*) MCFile->Get(Form("h_mEffCosS_Cent_%d_Pt_%d_Y_%d",icent,ipt,iy))->Clone();
  //      //h_theta_star_ratio[ipt][icent]->GetYaxis()->SetTitle(Form("Cent %d  %1.1f<pT<%1.1f",icent,vmsa::pt_low_cent[energy][ipt],vmsa::pt_up_cent[energy][ipt]));
  //      cout << Form("Cent %d  %1.1f<pT<%1.1f",icent,vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt]) << endl;
  //      h_theta_star_ratio[ipt][icent][iy]->GetXaxis()->SetTitle(Form("Cent %d  %1.1f <y<%1.1f  %1.1f<pT<%1.1f    cos#theta",icent,vmsa::y_bin[iy],vmsa::y_bin[iy+1],vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt]));
  //      //h_theta[ipt] = (TH1D*) MCFiles->Get("h_theta");
  //      //h_out1[ipt] = (TH1D*) MCFiles->Get("h_out1");
  //      //h_out2[ipt] = (TH1D*) MCFiles->Get("h_out2");
 
  //      //TF1 *Func_rho = new TF1("Func_rho","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
  //      //TF1 *Func_A = new TF1("Func_A","[0]*(1.+[1]*(x*x))",0,1);
  //      TF1 *Func_A = new TF1("Func_A",Func4th,0,1,3);

  //      //Func_rho->SetParameter(0,h_theta_star[ipt]->GetBinContent(1));
  //      //Func_rho->SetParameter(1,1./3.);
  //      //h_theta_star[ipt]->Fit(Func_rho,"ERQ");
  //      //cout << "rho00 = " << Func_rho->GetParameter(1) << endl;;    


  //      //Func_A->SetParameter(0,h_theta[ipt]->GetBinContent(1));
  //      //Func_A->SetParameter(1,0);
  //      //h_theta[ipt]->Fit(Func_A,"ER");
  //      //Fpval[ipt] = Func_A->GetParameter(1);
  //      //Fperr[ipt] = Func_A->GetParError(1);


  //      //TH1D *h_theta_star_clone = (TH1D*)h_theta_star[ipt]->Clone("h_theta_star_clone");
  //      //h_theta_star_clone->Sumw2();
  //      //h_theta_star_before[ipt]->Sumw2();
  //      //h_theta_star_clone->Divide(h_theta_star_before[ipt]); 
  //      Func_A->SetParameter(0,h_theta_star_ratio[ipt][icent][iy]->GetBinContent(1));
  //      Func_A->SetParameter(1,0);
  //      h_theta_star_ratio[ipt][icent][iy]->Fit(Func_A,"ER");

  //      Fval2[ipt][icent][iy] = Func_A->GetParameter(1);
  //      Ferr2[ipt][icent][iy] = Func_A->GetParError(1);
  //      Fval4[ipt][icent][iy] = Func_A->GetParameter(2);
  //      Ferr4[ipt][icent][iy] = Func_A->GetParError(2);

  //      //FfromFp[ipt] = -Fpval[ipt]/(2.+Fpval[ipt]);
  //      //FfromFperr[ipt] = TMath::Abs((Fpval[ipt]/(2.-Fpval[ipt])/(2.-Fpval[ipt])-1./(2.-Fpval[ipt]))*Fperr[ipt]);

  //      h_theta_star_ratio[ipt][icent][iy]->Draw("pE");
  //      Func_A->SetLineColor(kRed);
  //      Func_A->SetLineWidth(1);
  //      Func_A->Draw("same");
  //      //c_play->SaveAs(Form("FValues/pt%d.pdf",ipt));

  //      //delete c_play;
  //      //delete h_theta_star_clone;
  //      //delete Func_rho;
  //      //delete Func_A;
  //    }
  //    c_play->Update();
  //    c_play->Print(outputname.c_str());
  //  }
  //}
  //string output_stop = Form("%s]",outputname.c_str());
  //c_play->Print(output_stop.c_str()); // close pdf file
 

/*
  std::string outputname = Form("figures/Fplots_%s_%s_Mode%d_EtaMode%d_cent_thetastar_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),mode,etamode,ep.c_str(),xtraopt.c_str());
  std::string output_start = Form("%s[",outputname.c_str());
  c_play->Print(output_start.c_str());
  cout << "c_play->Print()" << endl;
  
  for(int ipt = vmsa::pt_rebin_first_cent[energy]; ipt <= vmsa::pt_rebin_last_cent[energy]; ipt++)
  {
    cout << "ipt = " << ipt << endl;
    for(int icent = 0; icent < 9; icent++)
    {
      cout << "icent = " << icent << endl;
      c_play->cd(icent+1);
      //if(ipt == 2) MCFiles[ipt] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/McAcceptanceOutput_pt%d_energy%d_pid%d_cent9_00004999.root",vmsa::mPID[pid].c_str(),ipt+1,energy,pid),"READ");

      //h_theta_star_before[ipt] = (TH1D*) MCFiles->Get("h_theta_star_before");
      h_theta_star_ratio[ipt][icent] = (TH1D*) MCFile->Get(Form("h_mEffCosS_Cent_%d_Pt_%d",icent,ipt))->Clone();
      //h_theta_star_ratio[ipt][icent]->GetYaxis()->SetTitle(Form("Cent %d  %1.1f<pT<%1.1f",icent,vmsa::pt_low_cent[energy][ipt],vmsa::pt_up_cent[energy][ipt]));
      cout << Form("Cent %d  %1.1f<pT<%1.1f",icent,vmsa::pt_low_cent[energy][ipt],vmsa::pt_up_cent[energy][ipt]) << endl;
      h_theta_star_ratio[ipt][icent]->GetXaxis()->SetTitle(Form("Cent %d    %1.1f<pT<%1.1f    cos#theta*",icent,vmsa::pt_low_cent[energy][ipt],vmsa::pt_up_cent[energy][ipt]));
      //h_theta[ipt] = (TH1D*) MCFiles->Get("h_theta");
      //h_out1[ipt] = (TH1D*) MCFiles->Get("h_out1");
      //h_out2[ipt] = (TH1D*) MCFiles->Get("h_out2");
 
      //TF1 *Func_rho = new TF1("Func_rho","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
      //TF1 *Func_A = new TF1("Func_A","[0]*(1.+[1]*(x*x))",0,1);
      TF1 *Func_A = new TF1("Func_A",Func4th,0,1,3);

      //Func_rho->SetParameter(0,h_theta_star[ipt]->GetBinContent(1));
      //Func_rho->SetParameter(1,1./3.);
      //h_theta_star[ipt]->Fit(Func_rho,"ERQ");
      //cout << "rho00 = " << Func_rho->GetParameter(1) << endl;;    


      //Func_A->SetParameter(0,h_theta[ipt]->GetBinContent(1));
      //Func_A->SetParameter(1,0);
      //h_theta[ipt]->Fit(Func_A,"ER");
      //Fpval[ipt] = Func_A->GetParameter(1);
      //Fperr[ipt] = Func_A->GetParError(1);


      //TH1D *h_theta_star_clone = (TH1D*)h_theta_star[ipt]->Clone("h_theta_star_clone");
      //h_theta_star_clone->Sumw2();
      //h_theta_star_before[ipt]->Sumw2();
      //h_theta_star_clone->Divide(h_theta_star_before[ipt]); 
      Func_A->SetParameter(0,h_theta_star_ratio[ipt][icent]->GetBinContent(1));
      Func_A->SetParameter(1,0);
      h_theta_star_ratio[ipt][icent]->Fit(Func_A,"ER");

      Fval[ipt][icent] = Func_A->GetParameter(1);
      Ferr[ipt][icent] = Func_A->GetParError(1);
      Gval[ipt][icent] = Func_A->GetParameter(2);
      Gerr[ipt][icent] = Func_A->GetParError(2);

      //Fval[ipt] = Func_A->GetParameter(1);
      //Ferr[ipt] = Func_A->GetParError(1);

      //FfromFp[ipt] = -Fpval[ipt]/(2.+Fpval[ipt]);
      //FfromFperr[ipt] = TMath::Abs((Fpval[ipt]/(2.-Fpval[ipt])/(2.-Fpval[ipt])-1./(2.-Fpval[ipt]))*Fperr[ipt]);

      h_theta_star_ratio[ipt][icent]->Draw("pE");
      Func_A->SetLineColor(kRed);
      Func_A->SetLineWidth(1);
      Func_A->Draw("same");
      //c_play->SaveAs(Form("FValues/pt%d.pdf",ipt));

      //delete c_play;
      //delete h_theta_star_clone;
      //delete Func_rho;
      //delete Func_A;
      
      c_play->Update();
      c_play->Print(outputname.c_str());
    }
  }
  string output_stop = Form("%s]",outputname.c_str());
  c_play->Print(output_stop.c_str()); // close pdf file

    
  TGraphAsymmErrors *gFvalues[10][2];

  TCanvas *c_play2 = new TCanvas("c_play2","c_play2",10,10,800,800);
  c_play2->Divide(2,2);
  for(int i = 0; i < 4; i++)
  {
    //c_play2->cd(i+1)->SetTopMargin(0.10);
    c_play2->cd(i+1)->SetLeftMargin(0.15);
    c_play2->cd(i+1)->SetBottomMargin(0.15);
    c_play2->cd(i+1)->SetGrid(0,0);
    c_play2->cd(i+1)->SetTicks(1,1);
  }

  string outputname = Form("figures/Fplots_%s_%s_Mode%d_EtaMode%d_centrality_Fvalues_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),mode,etamode,ep.c_str(),xtraopt.c_str());
  string output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());
  cout << "c_play2->Print()" << endl;

  TH1F *h_play = new TH1F("h_play","h_play",100,-0.5,99.5);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-1000.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetYaxis()->SetTitle("");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(0.0,80.0);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");

  for(int iF = 0; iF < 2; iF++)
  {
    for(int ipt = vmsa::pt_rebin_first_cent[energy]; ipt <= vmsa::pt_rebin_last_cent[energy]; ipt++)
    {
      c_play2->cd(ipt+1);
      cout << "ipt = " << ipt << endl;
      gFvalues[ipt][iF] = new TGraphAsymmErrors();
      for(int icent = 0; icent < 9; icent++)
      {
        cout << "icent = " << icent << endl;
        double centval = (cent_set[icent]+cent_set[icent+1])/2.;
        if(iF == 0) 
        {
          gFvalues[ipt][iF]->SetPoint(icent,centval,Fval[ipt][icent]);
          gFvalues[ipt][iF]->SetPointError(icent,0.0,0.0,Ferr[ipt][icent],Ferr[ipt][icent]);
        }
        if(iF == 1) 
        {
          gFvalues[ipt][iF]->SetPoint(icent,centval,Gval[ipt][icent]);
          gFvalues[ipt][iF]->SetPointError(icent,0.0,0.0,Gerr[ipt][icent],Gerr[ipt][icent]);
        }
        
        double min = TMath::MinElement(gFvalues[ipt][iF]->GetN(),gFvalues[ipt][iF]->GetY());
        double max = TMath::MaxElement(gFvalues[ipt][iF]->GetN(),gFvalues[ipt][iF]->GetY());
        if(min < 0.0)  min =  1.2*min;
        if(min > 0.0)  min =  0.6*min;
        if(min == 0.0) min = -0.2*max;
        if(max < 0.0)  max =  0.8*max;
        if(max > 0.0)  max =  1.2*max;
        if(max == 0.0) max = -0.2*min;
        h_play->GetYaxis()->SetRangeUser(min,max);
        h_play->GetXaxis()->SetTitle(Form("%1.1f<pT<%1.1f     Centrality (%)",vmsa::pt_low_cent[energy][ipt],vmsa::pt_up_cent[energy][ipt]));
        //h_play->SetTitle(Form("Cent %d  %1.1f<pT<%1.1f",icent,vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt]));
        h_play->DrawCopy("pE");
        gFvalues[ipt][iF]->SetMarkerStyle(20); 
        gFvalues[ipt][iF]->SetLineColor(kBlue); 
        gFvalues[ipt][iF]->SetMarkerColor(kBlue); 
        gFvalues[ipt][iF]->Draw("pE same"); 
      }
    }
    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }
  string output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file
  */  

  //if(isBesI)
  //{ 
  //  Fval[1] = -0.00989006;
  //  Fval[2] = -0.0102287;
  //  Fval[3] = -0.0102549;
  //  Fval[4] = -0.00800524;
  //  Fval[5] = -0.000652552;
  //}
  ////double FfromF[6] ={0.0};
  //for(int ipt = 0; ipt < 6; ipt++)
  //{
  //  //FfromF[ipt] = -Fpval[ipt]/(2.+Fpval[ipt]);
  //  cout << "pt bin: " <<  pt_set[ipt] << "-" << pt_set[ipt+1] << std::fixed << std::setprecision(7) << "GeV/c    F = " << Fval[ipt] << " +/- " << Ferr[ipt] << "    F* = " << Fpval[ipt] << " +/- " << Fperr[ipt] << "    F from F* = " << FfromFp[ipt] << " +/- " << FfromFperr[ipt] << endl;   
  //  //cout << "pt bin: " << std::fixed << std::setprecision(5) << pt_set[ipt] << "-" << pt_set[ipt+1] << "GeV/c    F = " << Fval[ipt] << " +/- " << Ferr[ipt] << "    F' = " << Fpval[ipt] << " +/- " << Fperr[ipt] << "    F from F' (swapping them)= " << -Fval[ipt]/(2.+Fval[ipt]) << " +/- " << TMath::Abs((Fval[ipt]/(2.-Fval[ipt])/(2.-Fval[ipt])-1./(2.-Fval[ipt]))*Ferr[ipt]) << endl;   
  ////out<<"-D'/(2+D'): "<<-D_theta/(2.+D_theta)<<endl;
  //}

  Char_t Char[9][10] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
  Double_t centCent[9] = {75.0,65.0,55.0,45.0,35.0,25.0,15.0,7.5,2.5};

  ////double eta_D[6] = {0.0855382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};//19GeV


  TF1 *line = new TF1("line","1/3",-0.5,8.5);
  line->SetLineStyle(2);
  line->SetLineWidth(1);

  TFile *fres1 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/Resolution_file_19GeV_EpdCorrections_4.root","READ");
  if(energy == 3) fres1 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/Resolution_file_14GeV_EpdCorrections_6.root","READ");
  if(energy == 0) fres1 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/Resolution_file_7GeV_EpdCorrections_6.root","READ");
  TFile *fres2 = new TFile("../TreeProduction/UtilityFilesEta1p5OfficialCent/file_19GeV_Resolution.root","READ");
  if(energy == 3) fres2 = new TFile("../TreeProduction/StRoot/Utility/Resolution/file_14GeV_Resolution.root","READ");
  if(energy == 0) fres2 = new TFile("../TreeProduction/StRoot/Utility/Resolution/file_7GeV_Resolution.root","READ");
  TFile  *fd12 = new TFile("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/Phi/file_19GeV_EpdFlow_4.root","READ");
  if(energy == 3) fd12 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/file_14GeV_EpdFlow_6.root","READ");
  if(energy == 0) fd12 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/file_7GeV_EpdFlow_6.root","READ");
     
  //TFile event("all_event_output.root");
  
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

  std::ofstream outfile(Form("resolutioncent_%s_%d.h",vmsa::mBeamEnergy[energy].c_str(),order));

  outfile << "const double r" << vmsa::mBeamEnergy[energy] << "ep" << order << "[3][2] = {" << std::endl;

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

  
  double Res_12[9] = {0.0};
  double Res_12_err[9] = {0.0};

  if(order == 1)
  {
    for(int icent = 0; icent < 9; icent++) {
      Res_12[icent] = resolution_1->GetBinContent(icent+1);
      Res_12_err[icent] = resolution_1->GetBinError(icent+1);
      outfile << "{" << Res_12[icent] << "," << Res_12_err[icent] << "}," << std::endl;
    }
  }
  if(order == 2)
  {
    for(int icent = 0; icent < 9; icent++) {
      Res_12[icent] = resolution_12->GetBinContent(icent+1);
      Res_12_err[icent] = resolution_12->GetBinError(icent+1);
      outfile << "{" << Res_12[icent] << "," << Res_12_err[icent] << "}," << std::endl;
    }
  }
  //for(int icent = 0; icent < 9; icent++) {
  //  Res_12[icent] = resolution_2sub->GetBinContent(icent+1);
  //}
  outfile << "};" << std::endl;
  outfile.close();

  if(doall)
  {
    double eff[10][10][7];
    double eff_error[10][10][7];

    //TFile *eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/Eff_%s_SingleParticle_2060.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()),"READ");
    TFile *eff_file;
    if(energy == 3 || energy == 4)
    {
      if(order == 1) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order1%s%s%s/Eff_%s_SingleParticle_noToF_Mode1_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),eptext.c_str(),v2text.c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
      if(order == 2) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2%s%s%s/Eff_%s_SingleParticle_noToF_Mode1_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),eptext.c_str(),v2text.c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
    }
    if(energy == 0)
    {
      if(order == 1) eff_file = new TFile(Form("../RcPhiEffCorr/order1/Eff_%s_SingleParticle_noToF_Mode1_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
      if(order == 2) eff_file = new TFile(Form("../RcPhiEffCorr/order2/Eff_%s_SingleParticle_noToF_Mode1_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
    }
    //if(order == 1) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order1/Eff_19GeV_SingleParticle_noToF_Mode1_EtaMode%d.root",etamode),"READ");
    //if(order == 2) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_FixedRes/Eff_19GeV_SingleParticle_noToF_Mode1_EtaMode%d.root",etamode),"READ");
    eff_file->Print();

    for(int icent = 0; icent < 9; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_cent; ipt++)
      {  
        for(int itheta = 0; itheta < 7; itheta++)
        {
          TString *Title = new TString(Form("h_mEffCos_Cent_%d_Pt_%d",icent,ipt));
          TH1D *eff_hist = (TH1D*)eff_file->Get(Title->Data());
          eff[icent][ipt][itheta] = eff_hist->GetBinContent(itheta+1);
          eff_error[icent][ipt][itheta] = eff_hist->GetBinError(itheta+1)/eff[icent][ipt][itheta];
          cout << "icent = " << icent << " ipt = " << ipt << " itheta = " << itheta << " eff = " << eff[icent][ipt][itheta] << " +/- " << eff_error[icent][ipt][itheta] << endl;
          delete Title;
          delete eff_hist;
        }
        
      }
    }
    eff_file->Close();

    TFile *input;
    if(order == 1) input = new TFile(Form("rho00/%s/%s/Poly/RawRhoCentSys_%s_Poly_FirstOrder.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str()),"READ");
    if(order == 2) input = new TFile(Form("rho00/%s/%s/Poly/RawRhoCentSys_%s_Poly.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str()),"READ");
    //if(random3D) *input = new TFile(Form("rho00/%s/%s/3DRandom/Raw%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()),"READ");
    
    TString *Title;

    //TH1F *clonePt; 
    
    string outputname;
    if(order == 1) outputname = Form("output/%s/%s/AccResRhoCentSys_%s_Poly_FirstOrder_pt_%d%s%s%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str(),ptbin,spectra.c_str(),eptext.c_str(),v2text.c_str());
    if(order == 2) outputname = Form("output/%s/%s/AccResRhoCentSys_%s_Poly_pt_%d%s%s%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str(),ptbin,spectra.c_str(),eptext.c_str(),v2text.c_str());
    TFile *output = new TFile(outputname.c_str(),"RECREATE");
    output->cd();

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
                //if(i_dca != 0 || i_sig != 0 || i_norm != 0 || i_sigma != 0 || i_method != 0) continue;
    
                TF1 *Func_rho = new TF1("Func_rho",FuncAFG,0,1,5);
                TF1 *Func_obs = new TF1("Func_obs","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
      
                TH1D *rho00_hist;
                TGraphAsymmErrors *g_rho00 = new TGraphAsymmErrors();
    
                double weight_rho00 = 0.0;
                double weight_error_rho00 = 0.0;
                double weight_all = 0.0;

                double cent_rho00[10] = {0.0};
                double cent_error_rho00[10] = {0.0};
                double cent_all[10] = {0.0};

                TCanvas *c_fit = new TCanvas("c_fit","c_fit",10,10,900,900);
                c_fit->Divide(3,3);
                for(int i = 0; i < 9; i++)
                {
                  c_fit->cd(i+1)->SetLeftMargin(0.15);
                  c_fit->cd(i+1)->SetBottomMargin(0.15);
                  c_fit->cd(i+1)->SetGrid(0,0);
                  c_fit->cd(i+1)->SetTicks(1,1);
                }
                Title = new TString(Form("fit/Phi/%s/Order%d/centrality/yield_pt_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_EtaMode%d_Divided.pdf",vmsa::mBeamEnergy[energy].c_str(),order,ptbin,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1,etamode));


                for(int icent = 0; icent < 9; ++icent) // Centrality loop
                {
                  string key = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ptbin,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                  cout << key << endl;

                  //string newkey = Form("EffCorrected_eta_%d_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",iy,ptbin,icent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());

                  TH1F *PtCos_raw = (TH1F*)input->Get(key.c_str())->Clone("PtCos_raw");
                  TH1F *PtCos = new TH1F(key.c_str(),key.c_str(), 7, 0, 1);
                  PtCos->Sumw2();
                  //delete Title;
                  for(int itheta = 0; itheta < 7; itheta++) 
                  {
                    float inte_mean = PtCos_raw->GetBinContent(itheta+1);
                    float inte_mean_error = PtCos_raw->GetBinError(itheta+1);
                    cout << "inte_mean  = " << inte_mean << " +/- " << inte_mean_error << endl;
                    PtCos->SetBinContent(itheta+1, inte_mean/eff[icent][ptbin][itheta]);
                    PtCos->SetBinError(itheta+1, inte_mean/eff[icent][ptbin][itheta]*TMath::Sqrt(inte_mean_error*inte_mean_error/inte_mean/inte_mean+eff_error[icent][ptbin][itheta]*eff_error[icent][ptbin][itheta]));
                  }
                  //PtCos->Divide(h_theta_star_ratio[ptbin][icent]);
                  PtCos->Write();
                  c_fit->cd(icent+1);
                  PtCos->GetXaxis()->SetTitleOffset(1.2);
                  PtCos->GetXaxis()->SetTitle(Form("|cos#theta*|     %1.1f<p_{T}<%1.1fGeV/c,Cent %d-%d",vmsa::pt_low_cent[energy][ptbin],vmsa::pt_up_cent[energy][ptbin],centVal[icent+1],centVal[icent]));
                  std::string title = Form("%1.1f<p_{T}<%1.1fGeV/c,Cent %d-%d",vmsa::pt_low_cent[energy][ptbin],vmsa::pt_up_cent[energy][ptbin],centVal[icent+1],centVal[icent]);
                  std::cout << title.c_str() << std::endl;
                  //PtCos->GetXaxis()->SetTitle("|cos#theta*|");
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
                  //Func_rho->SetParameter(1, 0.333);
                  //Func_rho->FixParameter(2, 0.0);//Fval[ptbin][icent]);
                  //Func_rho->FixParameter(3, 0.0);//Gval[ptbin][icent]);
                  //Func_rho->FixParameter(4, Res_12[icent]);
                  //if(random3D) Func_rho->FixParameter(3, 0.);
                  //cout << "Resolution = " << Func_rho->GetParameter(3);
                  //cout<<"Fit with real EP:"<<endl;
                  //PtCos->Fit(Func_rho, "NMI");
                  //Func_rho->Draw("same");

                  //Title = new TString(Form("fit/yield_pt_%d_cent_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s.pdf",PtBin-1,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str()));
           
                  //TCanvas *c1 = new TCanvas();
                  //c1->SetFillColor(0);
                  //c1->SetGrid(0,0);
                  //c1->SetTitle(0);
                  //c1->SetBottomMargin(0.15);
                  //c1->SetLeftMargin(0.15);
                  //PtCos->Draw();
                  //Func_rho->Draw("same");
                  //c1->SaveAs(Title->Data());
                  //delete c1;
                  //
                   
                  PtCos->Fit(Func_obs,"NMI"); // fit corrected distribution for observerved rho00
                  double rho_obs = Func_obs->GetParameter(1);
                  double rho_obs_err = Func_obs->GetParError(1);
                  Func_obs->SetLineColor(kRed);      
                  Func_obs->DrawCopy("same");

                  // Error calculation
                  double drhodobs = 4./(1.+3.*Res_12[icent]);  // calculate d(rho)/d(rho_obs)
                  double drhodR = -12.*(rho_obs - 1./3.)/(1.+3.*Res_12[icent])/(1.+3.*Res_12[icent]); // calculate d(rho)/d(R)

                  double real_rho = 1./3. + 4./(1.+3.*Res_12[icent])*(rho_obs - 1./3.);
                  double real_rho_error = TMath::Sqrt((rho_obs_err*rho_obs_err)*(drhodobs*drhodobs) + (Res_12_err[icent]*Res_12_err[icent])*(drhodR*drhodR));

                  //float real_rho = Func_rho->GetParameter(1);
                  //float real_rho_error = Func_rho->GetParError(1);
                  float weight = PtCos->Integral(1,7);

                  //if(PtBin>=3) {
                  //  weight_rho00 += real_rho*weight;
                  //  weight_error_rho00 += real_rho_error*real_rho_error*weight*weight;
                  //  weight_all += weight;
                  //}

                  weight_rho00 += real_rho*weight;
                  weight_error_rho00 += real_rho_error*real_rho_error*weight*weight;
                  weight_all += weight;

                  cent_rho00[icent] += real_rho*weight;
                  cent_error_rho00[icent] += real_rho_error*real_rho_error*weight*weight;
                  cent_all[icent] += weight;

                
                  //if(i_cent == 9 && PtBin == 3 && i_dca == 0 && i_sig == 0 && i_norm == 0 && i_sigma == 0 && i_method == 1)
                  //{ 
                  //   string name = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d",PtBin-1,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F);                             
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

                for(int icent = 0; icent < 9; icent++) 
                {
                  cent_rho00[icent] = cent_rho00[icent]/cent_all[icent];
                  cent_error_rho00[icent] = TMath::Sqrt(cent_error_rho00[icent])/cent_all[icent];
                }

                cout << "rho00 = " << weight_rho00<< " +/- " << weight_error_rho00 << endl;

                for(int icent = 0; icent < 9; icent++) 
                {
                  cout << cent_rho00[icent];
                  if(icent == 8) cout << " " << weight_rho00 << endl;
                  else cout << " ";
                }

                //for(int PtBin=1; PtBin<=6; PtBin++) {
                //  cout<<pt_error_rho00[PtBin-1];
                //  if(PtBin==6) cout<<" "<<weight_error_rho00<<endl;
                //  else cout<<" ";
                // }
 
                Title = new TString(Form("rhoRaw_pt_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ptbin,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             

                g_rho00 = new TGraphAsymmErrors();
                g_rho00->SetName(Title->Data());
                for(int icent = 0; icent < 9; icent++)
                {
                  double centmean = (cent_set[icent] + cent_set[icent+1])/2.0;
                  g_rho00->SetPoint(icent,centmean,cent_rho00[icent]);
                  g_rho00->SetPointError(icent,0.0,0.0,cent_error_rho00[icent],cent_error_rho00[icent]);
                  cout << "icent = " << icent << " centmean = " << centmean << " cent_rho00[icent] = " << cent_rho00[icent] << " +/- " << cent_error_rho00[icent] << endl;
                }
                delete Title;

                Title = new TString(Form("rhoFinalWeighted_pt_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ptbin,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             

                rho00_hist = new TH1D(Title->Data(), Title->Data(), 10, 0.5, 10.5);
                rho00_hist->SetBinContent(10,weight_rho00);
                rho00_hist->SetBinError(10,weight_error_rho00);
                for(int icent = 0; icent < 9; icent++) 
                {
                  rho00_hist->SetBinContent(icent+1,cent_rho00[icent]);
                  rho00_hist->SetBinError(icent+1,cent_error_rho00[icent]);
                }

                g_rho00->Write();
                rho00_hist->Write();

                delete Title;
                delete rho00_hist;
                delete g_rho00;
                delete Func_rho;
              }
            }
          }
        }
      }
    }    
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
