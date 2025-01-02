#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "Utility/type.h"
#include <string>

using namespace std;

void calculateFY_BinbyBin_2ndOrder(const int energy = 4, const int pid = 0, bool doall = true, bool isBesI = false, bool random3D = false, int mode = 1, int etamode = 0, int ypadding = 2, int ipt = 0, int icent = 1, std::string etastring = "eta1_eta1") {

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

  TH1D *h_theta_star_ratio[10][10][15];
  TH1D *h_theta_ratio[10][10][15];
  //TH1D *h_theta_star_before[6];
  //TH1D *h_theta[6];
  //TH1D *h_theta_star[6];
  //TH1D *h_out1[6];
  //TH1D *h_out2[6];

  TFile *MCFile;
  double Fval[10][10][15] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Ferr[10][10][15] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Gval[10][10][15] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Gerr[10][10][15] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};

  MCFile = new TFile(Form("./Rebinned/Acceptance_%s_19GeV_Mode%d_EtaMode%d.root",vmsa::mPID[pid].c_str(),mode,etamode),"READ");

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,1500,900);
  c_play->Divide(5,3);
  for(int i = 0; i < 15; i++)
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
 


  std::string outputname = Form("figures/Fplots_2ndOrder_%s_%s_Mode%d_EtaMode%d_rapidity_thetastar_cent%d_pt%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),mode,etamode,icent,ipt);
  std::string output_start = Form("%s[",outputname.c_str());
  c_play->Print(output_start.c_str());
  cout << "c_play->Print()" << endl;
  
  for(int iy = 0; iy < vmsa::y_total; iy++)
  {
    c_play->cd(iy+1);
    //if(ipt == 2) MCFiles[ipt] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/McAcceptanceOutput_pt%d_energy%d_pid%d_cent9_00004999.root",vmsa::mPID[pid].c_str(),ipt+1,energy,pid),"READ");

    //h_theta_star_before[ipt] = (TH1D*) MCFiles->Get("h_theta_star_before");
    h_theta_ratio[ipt][icent][iy] = (TH1D*) MCFile->Get(Form("h_mEffCosS_Cent_%d_Pt_%d_Y_%d",icent,ipt,iy))->Clone();
    //h_theta_star_ratio[ipt][icent]->GetYaxis()->SetTitle(Form("Cent %d  %1.1f<pT<%1.1f",icent,vmsa::pt_low_cent[energy][ipt],vmsa::pt_up_cent[energy][ipt]));
    cout << Form("Cent %d  %1.1f<pT<%1.1f",icent,vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt]) << endl;
    h_theta_ratio[ipt][icent][iy]->GetXaxis()->SetTitle(Form("Cent %d  %1.1f <y<%1.1f  %1.1f<pT<%1.1f    cos#theta*",icent,vmsa::y_bin[iy],vmsa::y_bin[iy+1],vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt]));
    //h_theta[ipt] = (TH1D*) MCFiles->Get("h_theta");
    //h_out1[ipt] = (TH1D*) MCFiles->Get("h_out1");
    //h_out2[ipt] = (TH1D*) MCFiles->Get("h_out2");
 
    //TF1 *Func_rho = new TF1("Func_rho","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
    TF1 *Func_A = new TF1("Func_A","[0]*(1.+[1]*(x*x))",0,1);
    //TF1 *Func_A = new TF1("Func_A",Func4th,0,1,3);

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
    Func_A->SetParameter(0,h_theta_ratio[ipt][icent][iy]->GetBinContent(1));
    Func_A->SetParameter(1,0);
    h_theta_ratio[ipt][icent][iy]->Fit(Func_A,"ER");

    Fval[ipt][icent][iy] = Func_A->GetParameter(1);
    Ferr[ipt][icent][iy] = Func_A->GetParError(1);
    //Gval[ipt][icent][iy] = Func_A->GetParameter(2);
    //Gerr[ipt][icent][iy] = Func_A->GetParError(2);

    //Fval[ipt] = Func_A->GetParameter(1);
    //Ferr[ipt] = Func_A->GetParError(1);

    //FfromFp[ipt] = -Fpval[ipt]/(2.+Fpval[ipt]);
    //FfromFperr[ipt] = TMath::Abs((Fpval[ipt]/(2.-Fpval[ipt])/(2.-Fpval[ipt])-1./(2.-Fpval[ipt]))*Fperr[ipt]);

    h_theta_ratio[ipt][icent][iy]->Draw("pE");
    Func_A->SetLineColor(kRed);
    Func_A->SetLineWidth(1);
    Func_A->Draw("same");
    //c_play->SaveAs(Form("FValues/pt%d.pdf",ipt));

    //delete c_play;
    //delete h_theta_star_clone;
    //delete Func_rho;
    //delete Func_A;
  }
  c_play->Update();
  c_play->Print(outputname.c_str());

  string output_stop = Form("%s]",outputname.c_str());
  c_play->Print(output_stop.c_str()); // close pdf file

    
  TGraphAsymmErrors *gFvalues[10][10][2];

  TCanvas *c_play2 = new TCanvas("c_play2","c_play2",10,10,800,800);
  c_play2->cd()->SetLeftMargin(0.15);
  c_play2->cd()->SetBottomMargin(0.15);
  c_play2->cd()->SetGrid(0,0);
  c_play2->cd()->SetTicks(1,1);

  string outputname = Form("figures/Fplots_2ndOrder_%s_%s_Mode%d_EtaMode%d_rapidity_Fvalues_pt%d_cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),mode,etamode,ipt,icent);
  string output_start = Form("%s[",outputname.c_str());
  c_play2->Print(output_start.c_str());
  cout << "c_play2->Print()" << endl;

  TH1F *h_play = new TH1F("h_play","h_play",100,-2.0,2.0);
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
  h_play->GetXaxis()->SetRangeUser(-1.5,1.5);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");

  for(int iF = 0; iF < 1; iF++)
  {
    c_play2->cd();
    gFvalues[ipt][icent][iF] = new TGraphAsymmErrors();
    cout << "icent = " << icent << endl;
    for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++)
    {
      double yval = (vmsa::y_bin[iy]+vmsa::y_bin[iy+1])/2.;
      if(iF == 0) 
      {
        gFvalues[ipt][icent][iF]->SetPoint(iy-ypadding,yval,Fval[ipt][icent][iy]);
        gFvalues[ipt][icent][iF]->SetPointError(iy-ypadding,0.0,0.0,Ferr[ipt][icent][iy],Ferr[ipt][icent][iy]);
      }
      if(iF == 1) 
      {
        //gFvalues[ipt][icent][iF]->SetPoint(iy-ypadding,yval,Gval[ipt][icent][iy]);
        //gFvalues[ipt][icent][iF]->SetPointError(iy-ypadding,0.0,0.0,Gerr[ipt][icent][iy],Gerr[ipt][icent][iy]);
      }
    }
    double min = TMath::MinElement(gFvalues[ipt][icent][iF]->GetN(),gFvalues[ipt][icent][iF]->GetY());
    double max = TMath::MaxElement(gFvalues[ipt][icent][iF]->GetN(),gFvalues[ipt][icent][iF]->GetY());
    if(min < 0.0)  min =  1.2*min;
    if(min > 0.0)  min =  0.6*min;
    if(min == 0.0) min = -0.2*max;
    if(max < 0.0)  max =  0.8*max;
    if(max > 0.0)  max =  1.2*max;
    if(max == 0.0) max = -0.2*min;
    h_play->GetYaxis()->SetRangeUser(min,max);
    h_play->GetXaxis()->SetTitle(Form("Cent %d  %1.1f<pT<%1.1f      y",icent,vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt]));
    //h_play->SetTitle(Form("Cent %d  %1.1f<pT<%1.1f",icent,vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt]));
    h_play->DrawCopy("pE");
    gFvalues[ipt][icent][iF]->SetMarkerStyle(20); 
    gFvalues[ipt][icent][iF]->SetLineColor(kBlue); 
    gFvalues[ipt][icent][iF]->SetMarkerColor(kBlue); 
    gFvalues[ipt][icent][iF]->Draw("pE same"); 
     
    c_play2->Update();
    c_play2->Print(outputname.c_str());
  }
  string output_stop = Form("%s]",outputname.c_str());
  c_play2->Print(output_stop.c_str()); // close pdf file
    

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
  TFile *fres2 = new TFile("../TreeProduction/UtilityFilesEta1p5OfficialCent/file_19GeV_Resolution.root","READ");
  //TFile *fres2 = new TFile("../TreeProduction/OldUtilityFilesEta1OfficialCent/file_19GeV_Resolution.root","READ");
  //TFile  *fd12 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/file_19GeV_EpdFlow.root","READ");
  TFile  *fd12 = new TFile("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/Phi/file_19GeV_EpdFlow_4.root","READ");
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

    cout<<"Cent"<<cent_set[cent+1]<<"_"<<cent_set[cent]<<": "<<ResZ<<" +/- "<<ResErrZ<<" "<<ResZ_2<<" +/- "<<ResErrZ_2<<" "<<mean<<" +/- "<<mean*TMath::Sqrt(error)<<endl;

    g_res1->SetPoint(cent,centCent[cent],ResZ);
    g_res1->SetPointError(cent,0.0,0.0,ResErrZ,ResErrZ);
    g_res2->SetPoint(cent,centCent[cent],ResZ_2);
    g_res2->SetPointError(cent,0.0,0.0,ResErrZ_2,ResErrZ_2);
    g_res12->SetPoint(cent,centCent[cent],mean);
    g_res12->SetPointError(cent,0.0,0.0,mean*TMath::Sqrt(error),mean*TMath::Sqrt(error));
  }

  
  double Res_12[vmsa::cent_rebin_total] = {0.0};
  double Res_12_weight[vmsa::cent_rebin_total] = {0.0};
  for(int icent_rebin = 0; icent_rebin < vmsa::cent_rebin_total; icent_rebin++)
  {
    for(int i_cent = vmsa::cent_rebin[icent_rebin]; i_cent < vmsa::cent_rebin[icent_rebin+1]; i_cent++) {
      Res_12[icent_rebin] += resolution_12->GetBinContent(i_cent+1)/resolution_12->GetBinError(i_cent+1)/resolution_12->GetBinError(i_cent+1);
      Res_12_weight[icent_rebin] += 1./resolution_12->GetBinError(i_cent+1)/resolution_12->GetBinError(i_cent+1);
    }
    Res_12[icent_rebin] = Res_12[icent_rebin]/Res_12_weight[icent_rebin];
    cout << "CentBin = " << icent_rebin << "    for cent range = [" << vmsa::cent_rebin[icent_rebin] << "," << vmsa::cent_rebin[icent_rebin+1]-1 << "]       with Res12 = " << Res_12[icent_rebin] << endl;
  }
  
  TString *Title;
  if(doall)
  {
    double eff[vmsa::y_total][7];
    double eff_error[vmsa::y_total][7];

    //TFile *eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/Eff_%s_SingleParticle_2060.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()),"READ");
    TFile *eff_file = new TFile(Form("../RcPhiEffCorr/Rebinned/Eff_19GeV_SingleParticle_noToF_Mode2_EtaMode%d.root",etamode),"READ");
    eff_file->Print();

    for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++)
    { 
      Title = new TString(Form("h_mEffCos_Cent_%d_Pt_%d_Y_%d",icent,ipt,iy));
      //TString *Title = new TString(Form("h_mEffCos_Cent_%d_Pt_%d_Y_%d",icent,ipt,iy));
      cout << Title->Data() << endl;
      TH1D *eff_hist = (TH1D*)eff_file->Get(Title->Data());
      eff_hist->Print();
      for(int itheta = 0; itheta < 7; itheta++)
      {
        eff[iy][itheta] = eff_hist->GetBinContent(itheta+1);
        eff_error[iy][itheta] = eff_hist->GetBinError(itheta+1)/eff[iy][itheta];
        cout << "icent = " << icent << " ipt = " << ipt << " iy = " << iy << " itheta = " << itheta << " eff = " << eff[iy][itheta] << " +/- " << eff_error[iy][itheta] << endl;
      }
      delete Title;
      delete eff_hist;
    }
    eff_file->Close();

    TFile *input = new TFile(Form("rho00/%s/%s/RawRhoEtaSys_%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str()),"READ");
    //if(random3D) *input = new TFile(Form("rho00/%s/%s/3DRandom/Raw%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()),"READ");
    
    //TString *Title;

    //TH1F *clonePt; 
    
    string outputname = Form("output/%s/%s/AccResRhoEtaSys_2ndOrder_%s_cent_%d_pt_%d.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str(),icent,ipt);
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
              //if(i_dca != 0 || i_sig != 0 || i_norm != 0 || i_sigma != 0 || i_method != 0) continue;
    

              TH1D *rho00_hist;
              TGraphAsymmErrors *g_rho00 = new TGraphAsymmErrors();
    
              double weight_rho00 = 0.0;
              double weight_error_rho00 = 0.0;
              double weight_all = 0.0;

              double y_rho00[vmsa::y_total] = {0.0};
              double y_error_rho00[vmsa::y_total] = {0.0};
              double y_all[vmsa::y_total] = {0.0};
            
              TCanvas *c_fit = new TCanvas("c_fit","c_fit",10,10,1500,900);
              c_fit->Divide(5,3);
              for(int i = 0; i < 15; i++)
              {
                //c_play->cd(i+1)->SetTopMargin(0.10);
                c_fit->cd(i+1)->SetLeftMargin(0.15);
                c_fit->cd(i+1)->SetBottomMargin(0.15);
                c_fit->cd(i+1)->SetGrid(0,0);
                c_fit->cd(i+1)->SetTicks(1,1);
              }
              Title = new TString(Form("fit/yield_pt_%d_cent_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s.pdf",ipt,icent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str()));
              
              TH1F *PtCos[vmsa::y_total];
              TF1 *Func_rho[vmsa::y_total];
              for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
              {
                string key = Form("eta_%d_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",iy,ipt,icent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
                cout << key << endl;

                //string newkey = Form("EffCorrected_eta_%d_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",iy,ipt,icent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());

                TH1F *PtCos_raw = (TH1F*)input->Get(key.c_str())->Clone("PtCos_raw");
                PtCos[iy] = new TH1F(key.c_str(),key.c_str(), 7, 0, 1);
                //delete Title;
                for(int itheta = 0; itheta < 7; itheta++) 
                {
                  float inte_mean = PtCos_raw->GetBinContent(itheta+1);
                  float inte_mean_error = PtCos_raw->GetBinError(itheta+1);
                  cout << "inte_mean  = " << inte_mean << " +/- " << inte_mean_error << endl;
                  PtCos[iy]->SetBinContent(itheta+1, inte_mean/eff[iy][itheta]);
                  PtCos[iy]->SetBinError(itheta+1, inte_mean/eff[iy][itheta]*TMath::Sqrt(inte_mean_error*inte_mean_error/inte_mean/inte_mean+eff_error[iy][itheta]*eff_error[iy][itheta]));
                }
                PtCos[iy]->Write();

                c_fit->cd(iy+1);
                PtCos[iy]->GetXaxis()->SetTitleOffset(1.2);
                PtCos[iy]->GetXaxis()->SetTitle("cos#theta*");
                PtCos[iy]->GetYaxis()->SetTitle("Efficiency Corrected Yield");
                PtCos[iy]->GetYaxis()->SetTitleOffset(1.0);
                PtCos[iy]->SetMarkerColor(kBlack);
                PtCos[iy]->SetMarkerSize(1.8);
                PtCos[iy]->SetMarkerStyle(20);
                PtCos[iy]->Draw("pE");
              
                //Func_rho[iy] = new TF1(Form("Func_rho_%d",iy),FuncAFG,0,1,6);
                Func_rho[iy] = new TF1(Form("Func_rho_%d",iy),FuncAD,0,1,4);
                Func_rho[iy]->SetParameter(0, PtCos[iy]->GetBinContent(5));
                Func_rho[iy]->SetParameter(1, 0.333);
                Func_rho[iy]->SetParLimits(1, 0.0, 1.0);
                Func_rho[iy]->FixParameter(2, Fval[icent][ipt][iy]);
               // Func_rho[iy]->FixParameter(3, Gval[icent][ipt][iy]);
                Func_rho[iy]->FixParameter(3, Res_12[icent]);
                //Func_rho[iy]->FixParameter(4, 1.0);
                //if(random3D) Func_rho->FixParameter(3, 0.);
                //cout << "Resolution = " << Func_rho->GetParameter(3);
                //cout<<"Fit with real EP:"<<endl;
                PtCos[iy]->Fit(Func_rho[iy], "NMRI");
                Func_rho[iy]->SetLineColor(kRed);
                Func_rho[iy]->Draw("same");

           
                float real_rho = Func_rho[iy]->GetParameter(1);
                float real_rho_error = Func_rho[iy]->GetParError(1);
                float weight = PtCos[iy]->Integral(1,7);

                //if(PtBin>=3) {
                //  weight_rho00 += real_rho*weight;
                //  weight_error_rho00 += real_rho_error*real_rho_error*weight*weight;
                //  weight_all += weight;
                //}

                weight_rho00 += real_rho*weight;
                weight_error_rho00 += real_rho_error*real_rho_error*weight*weight;
                weight_all += weight;

                y_rho00[iy] += real_rho*weight;
                y_error_rho00[iy] += real_rho_error*real_rho_error*weight*weight;
                y_all[iy] += weight;

              
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
                //delete PtCos;
                delete PtCos_raw;
              }
              c_fit->SaveAs(Title->Data());

              weight_rho00 = weight_rho00/weight_all;
              weight_error_rho00 = TMath::Sqrt(weight_error_rho00)/weight_all;

              for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
              {
                y_rho00[iy] = y_rho00[iy]/y_all[iy];
                y_error_rho00[iy] = TMath::Sqrt(y_error_rho00[iy])/y_all[iy];
              }

              cout << "rho00 = " << weight_rho00<< " +/- " << weight_error_rho00 << endl;

              for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
              {
                cout << y_rho00[iy];
                if(iy == vmsa::y_total-ypadding-1) cout << " " << weight_rho00 << endl;
                else cout << " ";
              }

              //for(int PtBin=1; PtBin<=6; PtBin++) {
              //  cout<<pt_error_rho00[PtBin-1];
              //  if(PtBin==6) cout<<" "<<weight_error_rho00<<endl;
              //  else cout<<" ";
             // }
 
              Title = new TString(Form("rhoRaw_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",ipt,icent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str()));                             

              g_rho00 = new TGraphAsymmErrors();
              g_rho00->SetName(Title->Data());
              for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++)
              {
                double ymean = (vmsa::y_bin[iy] + vmsa::y_bin[iy+1])/2.0;
                g_rho00->SetPoint(iy-ypadding,ymean,y_rho00[iy]);
                g_rho00->SetPointError(iy-ypadding,0.0,0.0,y_error_rho00[iy],y_error_rho00[iy]);
                cout << "iy = " << iy << " ymean = " << ymean << " y_rho00[iy] = " << y_rho00[iy] << " +/- " << y_error_rho00[iy] << endl;
              }
              delete Title;

              Title = new TString(Form("rhoFinalWeighted_pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",ipt,icent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str()));                             

              rho00_hist = new TH1D(Title->Data(), Title->Data(), vmsa::y_total+1, 0.5, (double)vmsa::y_total+1.5);
              rho00_hist->SetBinContent(vmsa::y_total+1,weight_rho00);
              rho00_hist->SetBinError(vmsa::y_total+1,weight_error_rho00);
              for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
              {
                rho00_hist->SetBinContent(iy+1,y_rho00[iy]);
                rho00_hist->SetBinError(iy+1,y_error_rho00[iy]);
              }

              g_rho00->Write();
              rho00_hist->Write();

              delete Title;
              delete rho00_hist;
              delete g_rho00;
              //delete Func_rho;
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

  //double result = 1. + (4.*F+3.*G)/8. - (2.*F+3.*G)/4.*CosTheta*CosTheta + 3.*G/8.*CosTheta*CosTheta*CosTheta*CosTheta;
   
  double order0 = 2. + F + 3.*G/4.;
  double order2 = (-1.*F - 3.*G/2.)*CosTheta*CosTheta;
  double order4 = (3.*G/4.)*CosTheta*CosTheta*CosTheta*CosTheta;
 
  double result = order0 + order2 + order4;

  return N*result;

}

double FuncAFG(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double G = par[3];
  double R = par[4];
  double rho10 = par[5];

  double A = (3.*rho-1.)/(1.-rho);
  double B = rho10/(1.-rho);
  //double As = A*(1.+3.*R)/(4.+A*(1.-R));
  //double Bs = A*(1.-R)/(4.+A*(1.-R));
  double denom = 4. + A*(1.-R) + B*(-1.+R);
  double As = (A*(1.+3.*R) + B*(3.-3.*R))/denom;
  double Bs = (A*(1.-R)    + B*(3.+R)   )/denom;


  double order0 = 2. + F - Bs*F/2. + 3.*G/4. - Bs*G/2.;    
  double order2 = (2.*As - F + As*F + Bs*F - 3.*G/2. + 3.*As*G/4. + 3.*Bs*G/2.)*CosTheta*CosTheta;
  double order4 = (-1.*As*F - Bs*F/2. + 3.*G/4. - 3.*As*G/2. - 3.*Bs*G/2.)*CosTheta*CosTheta*CosTheta*CosTheta;
  double order6 = (3.*As*G/4. + Bs*G/2.)*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta;

  double result = order0 + order2 + order4 + order6;

  return N*result;

}
