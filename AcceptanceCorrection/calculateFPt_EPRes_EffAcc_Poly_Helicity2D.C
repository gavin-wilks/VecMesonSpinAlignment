#include "Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "Utility/type.h"
#include "Utility/functions.h"
#include <string>
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH2.h"
#include "TH1.h"

void calculateFPt_EPRes_EffAcc_Poly_Helicity2D(const int energy = 4, const int pid = 0, bool doall = true, bool isBesI = false, bool random3D = false, int etamode = 0, int order = 2, int yspectra = 2, int EP = 0, int v2 = 1) {


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
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.4);
  gStyle->SetLegendFont(22);
  gStyle->SetLabelFont(22);
  gStyle->SetTitleFont(22);

  Double_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};
  Double_t pt_set[8] = {0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};
  Char_t Char[9][10] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
  Double_t centCent[9] = {75.0,65.0,55.0,45.0,35.0,25.0,15.0,7.5,2.5};



  if(doall)
  {
    double eff_plateau[10][10][9][12];
    double eff_plateau_error[10][10][9][12];
    //if(eff_mode==0)
    //TFile *eff_file_plateau = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/PlateauEtaStudy/Eff_%s_SingleParticle_noToF_Mode0_EtaMode%d.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
    
    TFile *eff_file_plateau;
      if(order == 1) eff_file_plateau = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order1%s%s%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),eptext.c_str(),v2text.c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
      //if(order == 2) eff_file_plateau = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut//Eff_%s_SingleParticle_noToF_Mode0_EtaMode%d_PtBins2_5.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
      //if(order == 2) eff_file_plateau = new TFile(Form("../RcPhiEffCorr_Cos2PhiStarPhi_PhiEmbed/HelicityCorrectionsMcDataRcWeight/Eff_%s_SingleParticle_noToF_Mode4_EtaMode%d_PtBins2_5.root",vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
      if(order == 2) eff_file_plateau = new TFile(Form("../RcPhiEffCorr_Cos2PhiStarPhi_PhiEmbed/HelicityCorrectionsFixedBins/Eff_%s_SingleParticle_noToF_Mode0_EtaMode%d_PtBins2_5.root",vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
   
    TH2D *eff_hist_plot[10];
    TH2D *histmc[10];
    TH2D *histrc[10];
    
    TH1D *eff_hist_plot_cos[10];
    TH1D *histmc_cos[10];
    TH1D *histrc_cos[10];

    TH1D *eff_hist_plot_phi[10];
    TH1D *histmc_phi[10];
    TH1D *histrc_phi[10];

    TH2D *yield_hist_plot[10];
    TH1D *yield_hist_plot_cos[10];
    TH1D *yield_hist_plot_phi[10];

    TH2F *rawyield_hist_plot[10];
    TH1F *rawyield_hist_plot_cos[10];
    TH1F *rawyield_hist_plot_phi[10];
    
    for(int i=9; i<10;i++)
    {
      for(int j=3; j<=6; j++)
      {
        string histnamemc = Form("h_mMcEffCosHPhiPrime_Cent_9_Pt_%d",j-1);
        string histnamerc = Form("h_mRc5EffCosHPhiPrime_Cent_9_Pt_%d",j-1);
        histmc[j-1] = (TH2D*) ((TH2D*) eff_file_plateau->Get(histnamemc.c_str()))->Clone(Form("mc%d",j-1));    
        histrc[j-1] = (TH2D*) ((TH2D*) eff_file_plateau->Get(histnamerc.c_str()))->Clone(Form("rc%d",j-1));    

        histmc_cos[j-1] = (TH1D*) histmc[j-1]->ProjectionX(Form("mccos%d",j-1),0,-1,"e");
        histrc_cos[j-1] = (TH1D*) histrc[j-1]->ProjectionX(Form("rccos%d",j-1),0,-1,"e");

        histmc_phi[j-1] = (TH1D*) histmc[j-1]->ProjectionY(Form("mcphi%d",j-1),0,-1,"e");
        histrc_phi[j-1] = (TH1D*) histrc[j-1]->ProjectionY(Form("rcphi%d",j-1),0,-1,"e");

        eff_hist_plot[j-1] = (TH2D*) histrc[j-1]->Clone(Form("eff%d",j-1));
        eff_hist_plot[j-1]->Divide(histrc[j-1],histmc[j-1],1,1,"B");
        eff_hist_plot[j-1]->Print();

        eff_hist_plot_cos[j-1] = (TH1D*) histrc_cos[j-1]->Clone(Form("effcos%d",j-1));
        eff_hist_plot_cos[j-1]->Divide(histrc_cos[j-1],histmc_cos[j-1],1,1,"B");
        eff_hist_plot_cos[j-1]->Print();

        eff_hist_plot_phi[j-1] = (TH1D*) histrc_phi[j-1]->Clone(Form("effphi%d",j-1));
        eff_hist_plot_phi[j-1]->Divide(histrc_phi[j-1],histmc_phi[j-1],1,1,"B");
        eff_hist_plot_phi[j-1]->Print();

        for(int k=1; k<=9; k++)
        {
          for(int l=1; l<=12; l++)
          {
            eff_plateau[i][j-1][k-1][l-1] = eff_hist_plot[j-1]->GetBinContent(k,l);
            eff_plateau_error[i][j-1][k-1][l-1] = eff_hist_plot[j-1]->GetBinError(k,l)/eff_plateau[i][j-1][k-1][l-1];
            cout << "efficiency = " << eff_plateau[i][j-1][k-1][l-1] << " +/- " << eff_plateau_error[i][j-1][k-1][l-1] << endl;
          }
        }
      }
    }
    //eff_file_plateau->Close();
    //delete eff_file_plateau ;
 
    cout << "Closed the efficiency file" << endl;

    TFile *input;
    if(order == 1) input = new TFile(Form("rho00/%s/%s/Poly/PhiPsi_Raw%sPtSys_%s_PolySys_FirstOrder.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str()),"READ");
    if(order == 2) input = new TFile(Form("rho00/%s/%s/Poly/RawPhiPtSys_%s_PolySys_Helicity_2D_OffDiag.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str()),"READ");
    if(random3D) *input = new TFile(Form("rho00/%s/%s/3DRandom/Raw%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()),"READ");
    if(isBesI) input = new TFile(Form("rho00/%s/%s/BESI/Raw%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()),"READ");

    TF1 *Func_rho = new TF1("Func_rho",FuncAFG,0,1,5);
    TF1 *Func_rdl = new TF1("Func_rdl","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);

    TH1D *Big_yield = new TH1D("Big_yield","Big_yield", 7, 0, 1);
    double big_cos[9][2] = {0};

    TString *Title;

    TGraphAsymmErrors *g_rho00;
    TGraphAsymmErrors *g_real;
    TGraphAsymmErrors *g_imag;
    TGraphAsymmErrors *g_rerho;
    TGraphAsymmErrors *g_imrho;
    
    Double_t ptbin_rho00[7] = {0.4,0.8,1.2, 1.8, 2.4, 3.0, 4.2};

    cout << "Before the loop" << endl;
    double rho00[10];
    double real[10];
    double imag[10];
    double rerho[10];
    double imrho[10];
    double rho00err[10];
    double realerr[10];
    double imagerr[10];
    double rerhoerr[10];
    double imrhoerr[10];
    TF2 *Func_obs[10];
    for(int i_cent = 9/*vmsa::Cent_start*/; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
      {
        for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
        {
          if( i_dca != 0 && i_sig != 0 ) continue;
          for(int i_norm = vmsa::Norm_start; i_norm < 1/*vmsa::Norm_stop*/; ++i_norm)
          {
            for(int i_sigma = vmsa::Sig_start; i_sigma < 1/*vmsa::Sig_stop*/; ++i_sigma)
            {
              for(int i_method = 1/*vmsa::Method_start*/; i_method < vmsa::Method_stop; ++i_method)
              {
                for(int i_poly = 0; i_poly < 1/*3*/; i_poly++)
                {
                  //cout << "Within the loop" << endl;
                  for(Int_t PtBin=3; PtBin<=6; PtBin++) 
                  {
                    //eff_hist_plot[PtBin-1]->Print();
                    string key = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",PtBin-1,i_cent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                    cout << key << endl;
                    //TH2D *PtCos_raw = (TH2D*)input->Get(key.c_str())->Clone("PtCos_raw");
                    rawyield_hist_plot[PtBin-1] = (TH2F*)((TH2F*)input->Get(key.c_str()))->Clone(Form("raweta%d",PtBin-1));// (TH2D*) PtCos_raw->Clone();
                     
                    //rawyield_hist_plot[PtBin-1] = (TH2D*)input->Get(key.c_str());
                    rawyield_hist_plot[PtBin-1]->Print();
                    yield_hist_plot[PtBin-1] = new TH2D(Form("correta%d",PtBin-1),Form("correta%d",PtBin-1), 9, -1, 1, 12, 0.0, TMath::Pi()*2.0);
                    for(int i=1; i<=9; i++) 
                    {
                      for(int j=1; j<=12; j++) 
                      {
                        float inte_mean = rawyield_hist_plot[PtBin-1]->GetBinContent(i,j);
                        float inte_mean_error = rawyield_hist_plot[PtBin-1]->GetBinError(i,j);
                        //cout << "Inte = " << inte_mean << " +/- " << inte_mean_error << endl;

                        double value = inte_mean/eff_plateau[i_cent][PtBin-1][i-1][j-1]; 
                        double value_err = inte_mean/eff_plateau[i_cent][PtBin-1][i-1][j-1]*TMath::Sqrt(inte_mean_error*inte_mean_error/inte_mean/inte_mean+eff_plateau_error[i_cent][PtBin-1][i-1][j-1]*eff_plateau_error[i_cent][PtBin-1][i-1][j-1]);
                        //cout << "corrected yield = " << value << " +/- " << value_err << endl;

                        yield_hist_plot[PtBin-1]->SetBinContent(i, j, value);
                        yield_hist_plot[PtBin-1]->SetBinError(i, j, value_err);
                   
                        //cout << "EfficiencyCalculation" << endl;
                      }
                    }
                    cout << "Right before 1D projections" << endl;
                    rawyield_hist_plot[PtBin-1]->Print();
                    rawyield_hist_plot_cos[PtBin-1] = (TH1F*) rawyield_hist_plot[PtBin-1]->ProjectionX(Form("rawyieldcos_%d",PtBin-1),0,-1,"e");
                    rawyield_hist_plot_phi[PtBin-1] = (TH1F*) rawyield_hist_plot[PtBin-1]->ProjectionY(Form("rawyieldphi_%d",PtBin-1),0,-1,"e");
                    yield_hist_plot_cos[PtBin-1] = (TH1D*) yield_hist_plot[PtBin-1]->ProjectionX(Form("yieldcos_%d",PtBin-1),0,-1,"e");
                    yield_hist_plot_phi[PtBin-1] = (TH1D*) yield_hist_plot[PtBin-1]->ProjectionY(Form("yieldphi_%d",PtBin-1),0,-1,"e");
                    cout << "After 1D projections" << endl;
                    
                    Func_obs[PtBin-1] = new TF2(Form("Func_obs%d",PtBin-1),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),6);
                  
                    cout << "PtBin - 1 = " << PtBin-1 << endl;

                    cout << "Before fit" << endl;
                    Func_obs[PtBin-1]->SetParameter(0,1./3.);
                    Func_obs[PtBin-1]->SetParameter(1,0.0);
                    Func_obs[PtBin-1]->SetParameter(2,0.0);
                    Func_obs[PtBin-1]->SetParameter(3,0.0);
                    Func_obs[PtBin-1]->SetParameter(4,0.0);
                    Func_obs[PtBin-1]->SetParameter(5,yield_hist_plot[PtBin-1]->GetMaximum());
                    yield_hist_plot[PtBin-1]->Print();
                    yield_hist_plot[PtBin-1]->Fit(Func_obs[PtBin-1],"NMRI"); // fit corrected distribution for observerved rho00
                    cout << "After fit" << endl;
                  
                    cout << "PtBin - 1 = " << PtBin-1 << endl;
                    rho00[PtBin-1] = Func_obs[PtBin-1]->GetParameter(0);
                    rho00err[PtBin-1] = Func_obs[PtBin-1]->GetParError(0);
                    real[PtBin-1] = Func_obs[PtBin-1]->GetParameter(1);
                    realerr[PtBin-1] = Func_obs[PtBin-1]->GetParError(1);
                    imag[PtBin-1] = Func_obs[PtBin-1]->GetParameter(2);
                    imagerr[PtBin-1] = Func_obs[PtBin-1]->GetParError(2);
                    rerho[PtBin-1] = Func_obs[PtBin-1]->GetParameter(3);
                    rerhoerr[PtBin-1] = Func_obs[PtBin-1]->GetParError(3);
                    imrho[PtBin-1] = Func_obs[PtBin-1]->GetParameter(4);
                    imrhoerr[PtBin-1] = Func_obs[PtBin-1]->GetParError(4);
                    
                    cout << "rho00 = " << rho00[PtBin-1] << " +/- " << rho00err[PtBin-1] << endl;
                    cout << "real = " << real[PtBin-1] << " +/- " << realerr[PtBin-1] << endl;
                    cout << "imag = " << imag[PtBin-1] << " +/- " << imagerr[PtBin-1] << endl;
                    cout << "rerho1n1 = " << rerho[PtBin-1] << " +/- " << rerhoerr[PtBin-1] << endl;
                    cout << "imrho1n1 = " << imrho[PtBin-1] << " +/- " << imrhoerr[PtBin-1] << endl;

                    //Func_obs->SetLineColor(kRed);      
                    //Func_obs->DrawCopy("same");
         
                  }

                  TCanvas *c = new TCanvas("c","c",10,10,800,800);
                  c->Divide(2,2);
                  for(int i = 0; i < 4; i++)
                  {   
                    c->cd(i+1);
                    c->cd(i+1)->SetLeftMargin(0.075);
                    c->cd(i+1)->SetRightMargin(0.2);
                    c->cd(i+1)->SetBottomMargin(0.15);
                    c->cd(i+1)->SetTicks(1,1);
                    c->cd(i+1)->SetGrid(0,0);

                    eff_hist_plot[i+2]->SetTitle(Form("Efficiency %1.1f<p_{T}<%1.1f",ptbin_rho00[i+2],ptbin_rho00[i+3]));
                    eff_hist_plot[i+2]->GetXaxis()->SetTitle("cos(#theta*)");
                    eff_hist_plot[i+2]->GetYaxis()->SetTitle("#phi'");
                    eff_hist_plot[i+2]->Draw("Colz");
                  }
                  c->SaveAs(Form("figures/%s/%s/Pt2DHelicityEfficiency.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));
                  for(int i = 0; i < 4; i++)
                  {   
                    c->cd(i+1);
                      
                    yield_hist_plot[i+2]->SetTitle(Form("Corrected Yields %1.1f<p_{T}<%1.1f",ptbin_rho00[i+2],ptbin_rho00[i+3]));
                    yield_hist_plot[i+2]->GetXaxis()->SetTitle("cos(#theta*)");
                    yield_hist_plot[i+2]->GetYaxis()->SetTitle("#phi'");
                    yield_hist_plot[i+2]->Draw("Colz");
                    //Func_obs[i+2]->Draw("same");
                  }
                  c->SaveAs(Form("figures/%s/%s/Pt2DHelicityEfficiencyCorrectedYields.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));
                  for(int i = 0; i < 4; i++)
                  {   
                    c->cd(i+1);
                      
                    yield_hist_plot[i+2]->SetTitle(Form("Corrected Yields %1.1f<p_{T}<%1.1f",ptbin_rho00[i+2],ptbin_rho00[i+3]));
                    yield_hist_plot[i+2]->GetXaxis()->SetTitle("cos(#theta*)");
                    yield_hist_plot[i+2]->GetYaxis()->SetTitle("#phi'");
                    yield_hist_plot[i+4]->Draw("Colz");
                    Func_obs[i+2]->Draw("Colz");
                  }
                  c->SaveAs(Form("figures/%s/%s/Pt2DHelicity2DFunc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));
                  for(int i = 0; i < 4; i++)
                  {   
                    c->cd(i+1);
              
                    rawyield_hist_plot[i+2]->Print();
                    rawyield_hist_plot[i+2]->SetTitle(Form("Raw Yields %1.1f<p_{T}<%1.1f",ptbin_rho00[i+2],ptbin_rho00[i+3]));
                    rawyield_hist_plot[i+2]->GetXaxis()->SetTitle("cos(#theta*)");
                    rawyield_hist_plot[i+2]->GetYaxis()->SetTitle("#phi'");
                    rawyield_hist_plot[i+2]->Draw("Colz");
                  }
                  c->SaveAs(Form("figures/%s/%s/Pt2DHelicityRawYields.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));
                  for(int i = 0; i < 4; i++)
                  { 
                    c->cd(i+1);           
   
                    eff_hist_plot_cos[i+2]->SetTitle(Form("Efficiency %1.1f<p_{T}<%1.1f",ptbin_rho00[i+2],ptbin_rho00[i+3]));
                    eff_hist_plot_cos[i+2]->GetXaxis()->SetTitle("cos(#theta*)");
                    eff_hist_plot_cos[i+2]->GetYaxis()->SetTitle("Efficiency");
                    eff_hist_plot_cos[i+2]->Draw("pE");
                  }
                  c->SaveAs(Form("figures/%s/%s/Pt1DCosHelicityEfficiency.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));
                  for(int i = 0; i < 4; i++)
                  {   
                    c->cd(i+1);
                      
                    yield_hist_plot_cos[i+2]->SetTitle(Form("Corrected Yields %1.1f<p_{T}<%1.1f",ptbin_rho00[i+2],ptbin_rho00[i+3]));
                    yield_hist_plot_cos[i+2]->GetXaxis()->SetTitle("cos(#theta*)");
                    yield_hist_plot_cos[i+2]->GetYaxis()->SetTitle("count");
                    yield_hist_plot_cos[i+2]->Draw("pE");
                  }
                  c->SaveAs(Form("figures/%s/%s/Pt1DCosHelicityEfficiencyCorrectedYields.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));
                  for(int i = 0; i < 4; i++)
                  {   
                    c->cd(i+1);
                      
                    rawyield_hist_plot_cos[i+2]->SetTitle(Form("Raw Yields %1.1f<p_{T}<%1.1f",ptbin_rho00[i+2],ptbin_rho00[i+3]));
                    rawyield_hist_plot_cos[i+2]->GetXaxis()->SetTitle("cos(#theta*)");
                    rawyield_hist_plot_cos[i+2]->GetYaxis()->SetTitle("count");
                    rawyield_hist_plot_cos[i+2]->Draw("pE");
                  }
                  c->SaveAs(Form("figures/%s/%s/Pt1DCosHelicityRawYields.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));
                  for(int i = 0; i < 4; i++)
                  {   
                    c->cd(i+1);           

                    eff_hist_plot_phi[i+2]->SetTitle(Form("Efficiency %1.1f<p_{T}<%1.1f",ptbin_rho00[i+2],ptbin_rho00[i+3]));
                    eff_hist_plot_phi[i+2]->GetXaxis()->SetTitle("phi'");
                    eff_hist_plot_phi[i+2]->GetYaxis()->SetTitle("Efficiency");
                    eff_hist_plot_phi[i+2]->Draw("pE");
                  }
                  c->SaveAs(Form("figures/%s/%s/Pt1DPhiHelicityEfficiency.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));
                  for(int i = 0; i < 4; i++)
                  {   
                    c->cd(i+1);
                      
                    yield_hist_plot_phi[i+2]->SetTitle(Form("Corrected Yields %1.1f<p_{T}<%1.1f",ptbin_rho00[i+2],ptbin_rho00[i+3]));
                    yield_hist_plot_phi[i+2]->GetXaxis()->SetTitle("phi'");
                    yield_hist_plot_phi[i+2]->GetYaxis()->SetTitle("count");
                    yield_hist_plot_phi[i+2]->Draw("pE");
                  }
                  c->SaveAs(Form("figures/%s/%s/Pt1DPhiHelicityEfficiencyCorrectedYields.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));
                  for(int i = 0; i < 4; i++)
                  {   
                    c->cd(i+1);
                      
                    rawyield_hist_plot_phi[i+2]->SetTitle(Form("Raw Yields %1.1f<p_{T}<%1.1f",ptbin_rho00[i+2],ptbin_rho00[i+3]));
                    rawyield_hist_plot_phi[i+2]->GetXaxis()->SetTitle("phi'");
                    rawyield_hist_plot_phi[i+2]->GetYaxis()->SetTitle("count");
                    rawyield_hist_plot_phi[i+2]->Draw("pE");
                  }
                  c->SaveAs(Form("figures/%s/%s/Pt1DPhiHelicityRawYields.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));
                }
              }
            }
          }
        }
      }
    }

    string outputname;
    if(order == 1) outputname = Form("output/%s/%s/PhiPsi_AccRes%sPtSys_%s_Poly_FirstOrder%s%s%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str(),spectra.c_str(),eptext.c_str(),v2text.c_str());
    if(order == 2) outputname = Form("output/%s/%s/AccPhiPtSys_%s_PolySys_Helicity_2D_OffDiag.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str(),spectra.c_str());
    if(isBesI) outputname = Form("output/%s/%s/BESI/AccRes%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
    TFile *output = new TFile(outputname.c_str(),"RECREATE");
    output->cd();

    for(int i_cent = 9/*vmsa::Cent_start*/; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
      {
        for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
        {
          if( i_dca != 0 && i_sig != 0 ) continue;
          for(int i_norm = vmsa::Norm_start; i_norm < 1/*vmsa::Norm_stop*/; ++i_norm)
          {
            for(int i_sigma = vmsa::Sig_start; i_sigma < 1/*vmsa::Sig_stop*/; ++i_sigma)
            {
              for(int i_method = 1/*vmsa::Method_start*/; i_method < vmsa::Method_stop; ++i_method)
              {
                for(int i_poly = 0; i_poly < 1/*3*/; i_poly++)
                {

                  Title = new TString(Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             
                  g_rho00 = new TGraphAsymmErrors();
                  g_rho00->SetName(Title->Data());
                  for(Int_t PtBin=3; PtBin<=6; PtBin++) {
                    double ptmean = (ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                    g_rho00->SetPoint(PtBin-3,ptmean,rho00[PtBin-1]);
                    g_rho00->SetPointError(PtBin-3,0.0,0.0,rho00err[PtBin-1],rho00err[PtBin-1]);
                  }
                  g_rho00->Write();
                  delete Title;

                  Title = new TString(Form("realRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             
                  g_real = new TGraphAsymmErrors();
                  g_real->SetName(Title->Data());
                  for(Int_t PtBin=3; PtBin<=6; PtBin++) {
                    //double ptmean = (float(PtBin)-5.5)/5.; //(ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                    double ptmean = (ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                    g_real->SetPoint(PtBin-3,ptmean,real[PtBin-1]);
                    g_real->SetPointError(PtBin-3,0.0,0.0,realerr[PtBin-1],realerr[PtBin-1]);
                  }
                  g_real->Write();
                  delete Title;

                  Title = new TString(Form("imagRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             
                  g_imag = new TGraphAsymmErrors();
                  g_imag->SetName(Title->Data());
                  for(Int_t PtBin=3; PtBin<=6; PtBin++) {
                    double ptmean = (ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                    //double ptmean = (float(PtBin)-5.5)/5.; //(ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                    g_imag->SetPoint(PtBin-3,ptmean,imag[PtBin-1]);
                    g_imag->SetPointError(PtBin-3,0.0,0.0,imagerr[PtBin-1],imagerr[PtBin-1]);
                  }
                  g_imag->Write();
                  delete Title;

                  Title = new TString(Form("rerho1n1Raw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             
                  g_rerho = new TGraphAsymmErrors();
                  g_rerho->SetName(Title->Data());
                  for(Int_t PtBin=3; PtBin<=6; PtBin++) {
                    double ptmean = (ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                    //double ptmean = (float(PtBin)-5.5)/5.; //(ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                    g_rerho->SetPoint(PtBin-3,ptmean,rerho[PtBin-1]);
                    g_rerho->SetPointError(PtBin-3,0.0,0.0,rerhoerr[PtBin-1],rerhoerr[PtBin-1]);
                  }
                  g_rerho->Write();
                  delete Title;

                  Title = new TString(Form("imrho1n1Raw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             
                  g_imrho = new TGraphAsymmErrors();
                  g_imrho->SetName(Title->Data());
                  for(Int_t PtBin=3; PtBin<=6; PtBin++) {
                    double ptmean = (ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                    //double ptmean = (float(PtBin)-5.5)/5.; //(ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                    g_imrho->SetPoint(PtBin-3,ptmean,imrho[PtBin-1]);
                    g_imrho->SetPointError(PtBin-3,0.0,0.0,imrhoerr[PtBin-1],imrhoerr[PtBin-1]);
                  }
                  g_imrho->Write();
                  delete Title;

                  g_rho00->Print();
                  g_real->Print();
                  g_imag->Print();
                  g_rerho->Print();
                  g_imrho->Print();
                  
                  output->Close();
                }
              }
            }
          }
        }
      }
    }
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

