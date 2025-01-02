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
#include <iostream>
#include <iomanip>

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void compareHelicityVsGlobal_Paper_All_Y_Eta_PtCut(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 23, std::string simmode = "Mc", int inputpar = 0)//defaultF = 0 is BESII, defaultF = 1 is BESI
{
  gStyle->SetOptStat(0);

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(200000);

  float etacuts[11] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

  const int start = 0;
  const int stop  = 5;

  const int nvar = 5;
  float sigmay[nvar] = {1./3.-0.2,1./3.-0.1,0.0,1./3.+0.1,1./3.+0.2};

  const int defaultfile = 3;

  //const int nkin = 6;
  //const int ptkin[nkin] = {1,1,2,2,5,5};
  //const int ykin[nkin]  = {0,1,0,1,0,1};

  const int nkin = 5;
  const float ptkin[nkin] = {0.1,0.3,0.5,0.7,0.9};
  const float ykin[nkin]  = {9.0,9.0,9.0,9.0,9.0};


  //const int nopts = 2;
  //std::string option = "v2_On_Off";
  //std::string date[nopts] = {"20240809g","20240809g"}; // nov2, prelimv2
  //std::string opts[nopts] = {"nov2","prelimv2"};
  //int color[nopts] = {kBlue, kOrange+7}; 
  //std::string label[nopts] = {"v_{2} OFF","v_{2} ON (Prelim)"};

  const int nopts = 1;

  //std::string option = "y_Variance_RC_EtaCuts_GlobalInput_pt2_individualy";
  std::string option = "y_Variance_RC_EtaCuts_GlobalInput_pt1_individualy_ptcuts";
  //std::string date[nopts] = {"20240809_2","20240815he","2040815he","20240815he","20240815he","20240815he","20240815he","20240815he","20240815he","20240815he"}; // nov2, prelimv2
  //std::string opts[nopts] = {"nov2","v20.0333","v20.0667","v20.1000","v20.1333","v20.1667","v20.2000","v20.2333","v20.2667","v20.3000"};
  //float vals[nopts] = {0.0,1./30.,2./30,3./30.,4./30.,5./30,6./30.,7./30.,8./30,9./30.,10./30.};
  float vals[nopts] = {9.0};
  //for(int i = 0; i < nopts; i++) vals[i] = float(i-5)/5.;


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

  std::string outputname; 
  std::string outputstart;
  std::string outputstop; 
  string param[5] = {"#Delta#rho_{00}","Re(#rho_{10}-#rho_{0-1})","Im(#rho_{10}-#rho_{0-1})","Re(#rho_{1-1})","Im(#rho_{1-1})"};

  TFile *File_InPutRC[nopts][nkin]; 
  TH2D *h_m2D_RC[nopts][5][nvar][nkin][11][9];  
  TH2D *h_m2D_RC_Global[nopts][5][nvar][nkin][11][9]; 
  //TH2D *h_m2D_RC[nopts][5][nvar][nkin][11];  
  TH1D *h_m1D_RC_Global_Cos[nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Global_Cos_All[nopts][5][nvar][nkin][11][9]; 
  TH1D *h_m1D_RC_Global_Ratio_Cos_All[nopts][5][nvar][nkin][11][9]; 
  TH1D *h_m1D_RC_Global_Corrected_Cos_All[nopts][5][nvar][nkin][11][9]; 
  TH1D *h_m1D_RC_Global_Beta[nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Global_Ratio_Cos[nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Global_Ratio_Beta[nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Cos[nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Beta[nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Ratio_Cos[nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Ratio_Beta[nopts][5][nvar][nkin][11][9][20]; 

  TH2D *h_m2D_RC_Global_Ratio[nopts][5][nvar][nkin][11][9];
  TH2D *h_m2D_RC_Ratio[nopts][5][nvar][nkin][11][9];
  TH2D *h_m2D_RC_Global_Corrected[nopts][5][nvar][nkin][11][9];
  TH2D *h_m2D_RC_Global_Corrected_BackToMC[nopts][5][nvar][nkin][11][9];
  TH2D *h_m2D_RC_Global_Corrected_BackToMC_E[nopts][5][nvar][nkin][11][9];
  TH2D *h_m2D_RC_Corrected[nopts][5][nvar][nkin][11][9];
 
  TH2D *h_m2D_Weight[nopts][nkin];

  TGraphAsymmErrors *g_m1D_RC_Global_Ratio_Cos[nopts][5][nvar][nkin][11][9][20]; 
  TGraphAsymmErrors *g_m1D_RC_Global_Ratio_Beta[nopts][5][nvar][nkin][11][9][20]; 
  TGraphAsymmErrors *g_m1D_RC_Ratio_Cos[nopts][5][nvar][nkin][11][9][20]; 
  TGraphAsymmErrors *g_m1D_RC_Ratio_Beta[nopts][5][nvar][nkin][11][9][20]; 

  for(int ikin = 0; ikin < nkin; ikin++)
  {
    for(int iopt = 0; iopt < nopts; iopt++)
    {
      string inputfileRC = Form("effaccfiles/%s/%s/PaperAllYRC_EtaFixed_GlobalInput_nov2_ptcuts/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0_pt%1.1f_y%1.1f.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),ptkin[ikin],vals[iopt]);
      File_InPutRC[iopt][ikin] = TFile::Open(inputfileRC.c_str());
      //for(int irho = 0; irho < 1; irho++)
      //{
      //  for(int i = 2; i < 3; i++)
      //  //for(int i = 0; i < nvar; i++)
      //  {
      //    for(int ieta = 0; ieta < 1; ieta+=10)
      //    {
      //      string KEY_2D_RC = Form("h_m%sEffCosPhiPrimeH_v2_0_rhoinput_%d_rho_%d_eta_%d",simmode.c_str(),irho,i,ieta);
      //      string KEY_2D_RC_Global = Form("h_m%sEffCosPhiPrime_v2_0_rhoinput_%d_rho_%d_eta_%d",simmode.c_str(),irho,i,ieta);
      //      string KEY_2D = Form("h_m%sEffCosPhiPrimeH_v2_0_rhoinput_%d_rho_%d_eta_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,iopt,ikin);
      //      string KEY_2D_Global = Form("h_m%sEffCosPhiPrime_v2_0_rhoinput_%d_rho_%d_eta_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,iopt,ikin);

      //      h_m2D_RC[iopt][irho][i][ikin][ieta] = (TH2D*) ((TH2D*) File_InPutRC[iopt][ikin]->Get(KEY_2D_RC.c_str()))->Clone(KEY_2D.c_str());
      //      h_m2D_RC_Global[iopt][irho][i][ikin][ieta] = (TH2D*) ((TH2D*) File_InPutRC[iopt][ikin]->Get(KEY_2D_RC_Global.c_str()))->Clone(KEY_2D_Global.c_str());
      //      if(irho == 0 && i == 2 && ieta == 0)
      //      {
      //        string KEY_Weight = Form("Weight_v2_0_%d_pt_%d",iopt,ikin);
      //        h_m2D_Weight[iopt][ikin] = (TH2D*) h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->Clone(KEY_Weight.c_str());

      //        double integral = h_m2D_Weight[iopt][ikin]->Integral();           
      //        int totalbins = h_m2D_Weight[iopt][ikin]->GetNbinsX() * h_m2D_Weight[iopt][ikin]->GetNbinsY();
 
      //        cout << "integral = " << integral << ", totalbins = " << totalbins << endl;
      //   
      //        for(int ix = 1; ix <= h_m2D_Weight[iopt][ikin]->GetNbinsX(); ix++)
      //        {
      //          for(int iy = 1; iy <= h_m2D_Weight[iopt][ikin]->GetNbinsY(); iy++)
      //          {
      //            double val = h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->GetBinContent(ix,iy);
      //            double weight = integral/double(totalbins)/val;
      // 
      //            h_m2D_Weight[iopt][ikin]->SetBinContent(ix,iy,weight);
      //            h_m2D_Weight[iopt][ikin]->SetBinError(ix,iy,0.0);
      //            cout << "weight = " << weight << endl;
      //          }
      //        }           
      //        //h_m2D_Weight[iopt][ikin]->SetDirectory(0);
      //      }
      //    }
      //  }
      //}
    }
  }

  for(int ikin = 0; ikin < nkin; ikin++)
  {
    for(int iopt = 0; iopt < nopts; iopt++)
    {
      //string inputfileRC = Form("effaccfiles/%s/%s/PaperAllYRC_EtaFixed_GlobalInput_nov2/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0_pt%1.1f_y%1.1f.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),ptkin[ikin],vals[iopt]);
      //File_InPutRC[iopt][ikin] = TFile::Open(inputfileRC.c_str());
      for(int irho = 0; irho < 5; irho++)
      {
        for(int i = 0; i < nvar; i++)
        //for(int i = 0; i < nvar; i++)
        {
          for(int ieta = 0; ieta < 11; ieta+=10)
          {
            for(int ipt = 0; ipt < 9; ipt++)
            {
              string KEY_2D_RC = Form("h_m%sEffCosPhiPrimeH_v2_0_rhoinput_%d_rho_%d_eta_%d_pt%d",simmode.c_str(),irho,i,ieta,ipt);
              string KEY_2D_RC_Global = Form("h_m%sEffCosPhiPrime_v2_0_rhoinput_%d_rho_%d_eta_%d_pt%d",simmode.c_str(),irho,i,ieta,ipt);
              string KEY_2D = Form("h_m%sEffCosPhiPrimeH_v2_0_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin);
              string KEY_2D_Global = Form("h_m%sEffCosPhiPrime_v2_0_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin);
              cout << KEY_2D_RC << endl;
              h_m2D_RC[iopt][irho][i][ikin][ieta][ipt] = (TH2D*) ((TH2D*) File_InPutRC[iopt][ikin]->Get(KEY_2D_RC.c_str()))->Clone(KEY_2D.c_str());
              cout << KEY_2D_RC_Global << endl;
              h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt] = (TH2D*) ((TH2D*) File_InPutRC[iopt][ikin]->Get(KEY_2D_RC_Global.c_str()))->Clone(KEY_2D_Global.c_str());

              string KEY_1D_RC_Global = Form("h_m%sEffCosPhiPrimeCos_v2_0_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin);
              h_m1D_RC_Global_Cos_All[iopt][irho][i][ikin][ieta][ipt] = (TH1D*) h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->ProjectionX(KEY_1D_RC_Global.c_str());

              //if(irho == 0 && i == 2 && ieta == 0)
              //{
              //  string KEY_Weight = Form("Weight_v2_0_%d_pt_%d",iopt,ptkin[ikin]);
              //  h_m2D_Weight[iopt][ikin] = (TH2D*) h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->Clone(KEY_Weight.c_str());

              //  double integral = h_m2D_Weight[iopt][ikin]->Integral();           
              //  int totalbins = h_m2D_Weight[iopt][ikin]->GetNbinsX() * h_m2D_Weight[iopt][ikin]->GetNbinsY();
 
              //  cout << "integral = " << integral << ", totalbins = " << totalbins << endl;
         
              //  for(int ix = 1; ix <= h_m2D_Weight[iopt][ikin]->GetNbinsX(); ix++)
              //  {
              //    for(int iy = 1; iy <= h_m2D_Weight[iopt][ikin]->GetNbinsY(); iy++)
              //    {
              //      double val = h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->GetBinContent(ix,iy);
              //      double weight = integral/double(totalbins)/val;
       
              //      h_m2D_Weight[iopt][ikin]->SetBinContent(ix,iy,weight);
              //      //cout << "weight = " << weight << endl;
              //    }
              //  }           
              //  //h_m2D_Weight[iopt][ikin]->SetDirectory(0);
              //}
              //h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->Multiply(h_m2D_Weight[iopt][ikin]); 
              //h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->Multiply(h_m2D_Weight[iopt][ikin]); 

              //for(int ibin = 0; ibin < 20; ibin++)
              //{
              //  string KEY = Form("h_m%sEffCosPhiPrimeCos_v2_6_rhoinput_%d_rho_%d_eta_%d_y_%d_kin_%d_bin_%d",simmode.c_str(),irho,i,ieta,iopt,ikin,ibin);
              //  h_m1D_RC_Global_Cos[iopt][irho][i][ikin][ieta][ipt][ibin] = (TH1D*) h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->ProjectionX(KEY.c_str(),ibin+1,ibin+1);
              //  //h_m1D_RC_Global_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->SetDirectory(0);
              //  KEY = Form("h_m%sEffCosPhiPrimeCosH_v2_6_rhoinput_%d_rho_%d_eta_%d_y_%d_kin_%d_bin_%d",simmode.c_str(),irho,i,ieta,iopt,ikin,ibin);
              //  h_m1D_RC_Cos[iopt][irho][i][ikin][ieta][ipt][ibin] = (TH1D*) h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->ProjectionX(KEY.c_str(),ibin+1,ibin+1);
              //  //h_m1D_RC_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->SetDirectory(0);
              //}
              //for(int ibin = 0; ibin < 20; ibin++)
              //{
              //  string KEY = Form("h_m%sEffCosPhiPrimeBeta_v2_6_rhoinput_%d_rho_%d_eta_%d_y_%d_kin_%d_bin_%d",simmode.c_str(),irho,i,ieta,iopt,ikin,ibin);
              //  h_m1D_RC_Global_Beta[iopt][irho][i][ikin][ieta][ipt][ibin] = (TH1D*) h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->ProjectionY(KEY.c_str(),ibin+1,ibin+1);
              //  //h_m1D_RC_Global_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->SetDirectory(0);
              //  KEY = Form("h_m%sEffCosPhiPrimeBetaH_v2_6_rhoinput_%d_rho_%d_eta_%d_y_%d_kin_%d_bin_%d",simmode.c_str(),irho,i,ieta,iopt,ikin,ibin);
              //  h_m1D_RC_Beta[iopt][irho][i][ikin][ieta][ipt][ibin] = (TH1D*) h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->ProjectionY(KEY.c_str(),ibin+1,ibin+1);
              //  //h_m1D_RC_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->SetDirectory(0);
              //}


              

              //if(ieta >= 1)
              //{
                string KEY = Form("h_m%s2DGLOBALCOSRATIO_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin);
                h_m1D_RC_Global_Ratio_Cos_All[iopt][irho][i][ikin][ieta][ipt] = (TH1D*) h_m1D_RC_Global_Cos_All[iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                h_m1D_RC_Global_Ratio_Cos_All[iopt][irho][i][ikin][ieta][ipt]->Divide(h_m1D_RC_Global_Cos_All[iopt][irho][i][ikin][ieta][ipt],h_m1D_RC_Global_Cos_All[iopt][irho][i][ikin][0][0],1,1,"B");

                KEY = Form("h_m%s2DGLOBALRATIO_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin);
                h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt],h_m2D_RC_Global[iopt][irho][i][ikin][0][0],1,1,"B");
                KEY = Form("h_m%s2DHELICITYRATIO_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin);
                h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC[iopt][irho][i][ikin][ieta][ipt],h_m2D_RC[iopt][irho][i][ikin][0][0],1,1,"B");


                KEY = Form("h_m%s2DGLOBALCOSCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin);
                h_m1D_RC_Global_Corrected_Cos_All[iopt][irho][i][ikin][ieta][ipt] = (TH1D*) h_m1D_RC_Global_Cos_All[iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                h_m1D_RC_Global_Corrected_Cos_All[iopt][irho][i][ikin][ieta][ipt]->Divide(h_m1D_RC_Global_Ratio_Cos_All[iopt][irho][i][ikin][ieta][ipt]);

                KEY = Form("h_m%s2DGLOBALCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin);
                h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]);
                KEY = Form("h_m%s2DGLOBALCORRECTEDBTMC_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin);
                h_m2D_RC_Global_Corrected_BackToMC[iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                h_m2D_RC_Global_Corrected_BackToMC[iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Global[iopt][irho][i][ikin][0][ipt]);
                KEY = Form("h_m%s2DHELICITYCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin);
                h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt]);

       //         h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta-1]->SetDirectory(0);
       //         h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta-1]->SetDirectory(0);
       //         h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta-1]->SetDirectory(0);
       //         h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta-1]->SetDirectory(0);
              //}
       //       h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->SetDirectory(0);
       //       h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->SetDirectory(0);
            }
          }
        }
      }
      //File_InPutRC[iopt][ikin]->Close();
      //cout << inputfileRC << endl;
    }
  }


  for(int ikin = 0; ikin < nkin; ikin++)
  {
    for(int iopt = 0; iopt < nopts; iopt++)
    {
      for(int irho = 0; irho < 5; irho++)
      {
        for(int i = 0; i < nvar; i++)
        //for(int i = 0; i < nvar; i++)
        {
          for(int ieta = 0; ieta < 11; ieta+=10)
          {
            for(int ipt = 0; ipt < 9; ipt++) 
            {
	      string KEY = Form("h_m%s2DGLOBALCOSCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin);
	      h_m1D_RC_Global_Corrected_Cos_All[iopt][irho][i][ikin][ieta][ipt] = (TH1D*) h_m1D_RC_Global_Cos_All[iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
	      h_m1D_RC_Global_Corrected_Cos_All[iopt][irho][i][ikin][ieta][ipt]->Divide(h_m1D_RC_Global_Ratio_Cos_All[iopt][irho][2][ikin][ieta][ipt]);

              KEY = Form("h_m%s2DGLOBALCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin);
              h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
              h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Global_Ratio[iopt][irho][2][ikin][ieta][ipt]);

              //for(int ibin = 0; ibin < 20; ibin++)
              //{
              //  h_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin] = (TH1D*) h_m1D_RC_Global_Cos[iopt][irho][i][ikin][ieta+1][ipt][ibin];
              //  h_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin] = (TH1D*) h_m1D_RC_Global_Beta[iopt][irho][i][ikin][ieta+1][ipt][ibin];

              //  h_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->Divide(h_m1D_RC_Global_Cos[iopt][irho][i][ikin][ieta+1][ipt][ibin],h_m1D_RC_Global_Cos[iopt][irho][i][ikin][0][0][ibin],1,1,"B");
              //  h_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->Divide(h_m1D_RC_Global_Beta[iopt][irho][i][ikin][ieta+1][ipt][ibin],h_m1D_RC_Global_Beta[iopt][irho][i][ikin][0][0][ibin],1,1,"B");

              //  h_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin] = (TH1D*) h_m1D_RC_Cos[iopt][irho][i][ikin][ieta+1][ipt][ibin];
              //  h_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin] = (TH1D*) h_m1D_RC_Beta[iopt][irho][i][ikin][ieta+1][ipt][ibin];

              //  h_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->Divide(h_m1D_RC_Cos[iopt][irho][i][ikin][ieta+1][ipt][ibin],h_m1D_RC_Cos[iopt][irho][i][ikin][0][0][ibin],1,1,"B");
              //  h_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->Divide(h_m1D_RC_Beta[iopt][irho][i][ikin][ieta+1][ipt][ibin],h_m1D_RC_Beta[iopt][irho][i][ikin][0][0][ibin],1,1,"B");
              //}
            }
          }
        }
      }
    }
  }





  TGraphAsymmErrors *g2D[5][nvar][nkin][5][nopts][11]; 
  TGraphAsymmErrors *g2D_Global[5][nvar][nkin][5][nopts][11]; 
  TGraphAsymmErrors *g2D_Global_Corrected[5][nvar][nkin][5][nopts][11]; 
  TGraphAsymmErrors *g2D_Corrected[5][nvar][nkin][5][nopts][11]; 
  TF2 *fits2D[nopts][5][nvar][nkin][11][9];
  TF2 *fits2D_Global[nopts][5][nvar][nkin][11][9];
  TF2 *fits2D_Global_Corrected[nopts][5][nvar][nkin][11][9];
  TF2 *fits2D_Corrected[nopts][5][nvar][nkin][11][9];


  TGraphAsymmErrors *g1D_Global[5][nvar][nkin][nopts][11]; 
  TGraphAsymmErrors *g1D_Global_Corrected[5][nvar][nkin][nopts][11]; 
  TF1 *fits1D_Global[nopts][5][nvar][nkin][11][9];
  TF1 *fits1D_Global_Corrected[nopts][5][nvar][nkin][11][9];

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int i = 0; i < nvar; i++)
    {
      for(int ikin = 0; ikin < nkin; ikin++)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          for(int ieta = 0; ieta < 11; ieta+=10)
          {
            g1D_Global[irho][i][ikin][iopt][ieta] = new TGraphAsymmErrors();
            g1D_Global_Corrected[irho][i][ikin][iopt][ieta] = new TGraphAsymmErrors();
            for(int ipt = 0; ipt < 9; ipt++) 
            {
              cout << "ipt = " << ipt << ", ieta = " << ieta << endl;
             
              fits1D_Global[iopt][irho][i][ikin][ieta][ipt] = new TF1(Form("fit1D_GlobalCos_%d_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta,ipt),SpinDensity,-1.0,1.0,2);
              fits1D_Global[iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,1./3.);
              fits1D_Global[iopt][irho][i][ikin][ieta][ipt]->SetParLimits(0,0.0,1.0);
              fits1D_Global[iopt][irho][i][ikin][ieta][ipt]->SetParameter(1,h_m1D_RC_Global_Cos_All[iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
              h_m1D_RC_Global_Cos_All[iopt][irho][i][ikin][ieta][ipt]->Fit(fits1D_Global[iopt][irho][i][ikin][ieta][ipt],"NMRI");
              if(irho == 0) g1D_Global[irho][i][ikin][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025,fits1D_Global[iopt][irho][i][ikin][ieta][ipt]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
              if(irho != 0) g1D_Global[irho][i][ikin][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025,fits1D_Global[iopt][irho][i][ikin][ieta][ipt]->GetParameter(0)-1./3.);
              g1D_Global[irho][i][ikin][iopt][ieta]->SetPointError(ipt,0.0,0.0,fits1D_Global[iopt][irho][i][ikin][ieta][ipt]->GetParError(0),fits1D_Global[iopt][irho][i][ikin][ieta][ipt]->GetParError(0));

              //if(ieta >= 1) 
             // {
                fits1D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt] = new TF1(Form("fit1D_Global_CorrectedCos_%d_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta,ipt),SpinDensity,-1.0,1.0,2);
                fits1D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,1./3.);
                fits1D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetParLimits(0,0.0,1.0);
                fits1D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetParameter(1,h_m1D_RC_Global_Corrected_Cos_All[iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
                h_m1D_RC_Global_Corrected_Cos_All[iopt][irho][i][ikin][ieta][ipt]->Fit(fits1D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt],"NMRI");
                if(irho == 0) g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025,fits1D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
                if(irho != 0) g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025,fits1D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetParameter(0)-1./3.);
                g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->SetPointError(ipt,0.0,0.0,fits1D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetParError(0),fits1D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetParError(0));
             // }
              cout << "FIT GLOBAL CORRECTED" << endl;
            }
          }
        }
      }
    }
  }


  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int i = 0; i < nvar; i++)
    {
      for(int ikin = 0; ikin < nkin; ikin++)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          for(int ieta = 0; ieta < 11; ieta+=10)
          {
            for(int par = 0; par < 5; par++) 
            {
              g2D[irho][i][ikin][par][iopt][ieta] = new TGraphAsymmErrors();
              g2D_Global[irho][i][ikin][par][iopt][ieta] = new TGraphAsymmErrors();
            }   
            for(int ipt = 0; ipt < 9; ipt++)
            {
              //fits2D[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),6);
              //fits2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,1./3.);
              //fits2D[iopt][irho][i][ikin][ieta][ipt]->SetParLimits(0,0.0,1.0);
              //fits2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(1,0.0);
              //fits2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(2,0.0);
              //fits2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(3,0.0);
              //fits2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(4,0.0);
              //fits2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(5,h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
              //h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->Fit(fits2D[iopt][irho][i][ikin][ieta][ipt],"NMR");
              //for(int par = 0; par < 5; par++)
              //{
              //  if(par == 0 && par == irho) g2D[irho][i][ikin][par][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025, fits2D[iopt][irho][i][ikin][ieta][ipt]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
              //  if(par != 0 && par == irho) g2D[irho][i][ikin][par][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025, fits2D[iopt][irho][i][ikin][ieta][ipt]->GetParameter(par)-(float(i)-2)*0.1);
              //  if(par == 0 && par != irho) g2D[irho][i][ikin][par][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025, fits2D[iopt][irho][i][ikin][ieta][ipt]->GetParameter(0)-1./3.);
              //  if(par != 0 && par != irho) g2D[irho][i][ikin][par][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025, fits2D[iopt][irho][i][ikin][ieta][ipt]->GetParameter(par));
              //  g2D[irho][i][ikin][par][iopt][ieta]->SetPointError(ipt,0.0,0.0,fits2D[iopt][irho][i][ikin][ieta][ipt]->GetParError(par),fits2D[iopt][irho][i][ikin][ieta][ipt]->GetParError(par));
              //}
              //cout << "Fit helicity " << endl;

              fits2D_Global[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_Global_%d_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta,ipt),SpinDensity2Dcos,-1.0,1.0,0,2.0*TMath::Pi(),6);
              fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,1./3.);
              fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->SetParLimits(0,0.0,1.0);
              fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->SetParameter(1,0.0);
              fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->SetParameter(2,0.0);
              fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->SetParameter(3,0.0);
              fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->SetParameter(4,0.0);
              fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->SetParameter(5,h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
              //h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->Fit(fits2D_Global[iopt][irho][i][ikin][ieta][ipt],"NMRI");
              for(int par = 0; par < 5; par++)
              {
                if(par == 0 && par == irho) g2D_Global[irho][i][ikin][par][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025, fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
                if(par != 0 && par == irho) g2D_Global[irho][i][ikin][par][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025, fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->GetParameter(par)-(float(i)-2)*0.1);
                if(par == 0 && par != irho) g2D_Global[irho][i][ikin][par][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025, fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->GetParameter(0)-1./3.);
                if(par != 0 && par != irho) g2D_Global[irho][i][ikin][par][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025, fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->GetParameter(par));
                g2D_Global[irho][i][ikin][par][iopt][ieta]->SetPointError(ipt,0.0,0.0,fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->GetParError(par),fits2D_Global[iopt][irho][i][ikin][ieta][ipt]->GetParError(par));
              }
            } 
          }
        }
      }
    }
  }

  TCanvas *cfit2D = new TCanvas("cfit2D","cfit2D",10,10,2000,2000);
  cfit2D->Divide(5,5);
  
  for(int i = 0; i < 25; i++)
  {
    cfit2D->cd(i+1);
    cfit2D->cd(i+1)->SetLeftMargin(0.175);  
    cfit2D->cd(i+1)->SetBottomMargin(0.175);
    cfit2D->cd(i+1)->SetTicks(1,1);
    cfit2D->cd(i+1)->SetGrid(0,0);
  }
 

  outputname = Form("figures/%s/%s/pTstudy/GlobalvsGlobal_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  cfit2D->Print(outputstart.c_str());

  //TF1 *rhoGv2[5][5][nvar][nkin][11][9];

  TLegend *leg[5];// = new TLegend(0.4,0.75,0.7,0.9);

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int ieta = 0; ieta < 11; ieta+=10)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          if(ikin == 0 && iopt == 0 && ieta == 0) 
          {
            leg[irho] = new TLegend(0.25,0.775,0.85,0.865);
            leg[irho]->SetNColumns(2);
            leg[irho]->SetBorderSize(0);
          }
          //double max = -999999;
          //double min = 999999;
          //for(int i = 0; i < nvar; i++)
          //{
          //  for(int par = 0; par < 5; par++)
          //  {
          //    double tmax = g2D_Global[irho][i][ikin][par][iopt]->GetHistogram()->GetMaximum();   
          //    double tmin = g2D_Global[irho][i][ikin][par][iopt]->GetHistogram()->GetMinimum(); 

          //    if(tmax > max) max = tmax;
          //    if(tmin < min) min = tmin;
          //  }
          //}
          for(int par = 0; par < 5; par++)
          {
            //cfit2D->cd(1+irho+5*par);
            cfit2D->cd(1+irho+5*par);
            double max = -999999;
            double min = 999999;
            for(int i = 0; i < nvar; i++)
            {
                double tmax = g2D_Global[irho][i][ikin][par][iopt][ieta]->GetHistogram()->GetMaximum();   
                double tmin = g2D_Global[irho][i][ikin][par][iopt][ieta]->GetHistogram()->GetMinimum(); 

                if(tmax > max) max = tmax;
                if(tmin < min) min = tmin;
            }
            
            //double max = -999999;
            //double min = 999999;
            //for(int i = 0; i < nvar; i++)
            //{
            //  for(int ir = 0; ir < 5; ir++)
            //  {
            //    double tmax = g2D_Global[ir][i][ikin][par]->GetHistogram()->GetMaximum();   
            //    double tmin = g2D_Global[ir][i][ikin][par]->GetHistogram()->GetMinimum(); 

            //    if(tmax > max) max = tmax;
            //    if(tmin < min) min = tmin;
            //  }
            //}
            for(int i = 0; i < nvar; i++)
            {
            //  switch(par)
            //  {
            //    case 0: 
            //      rhoGv2[par][irho][i][ikin][ieta][ipt] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),rhoGfromHy,-1.0,1.0,7);
            //      break;
            //    case 1: 
            //      rhoGv2[par][irho][i][ikin][ieta][ipt] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),realGfromHy,-1.0,1.0,7);
            //      break;
            //    case 2: 
            //      rhoGv2[par][irho][i][ikin][ieta][ipt] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),imagGfromHy,-1.0,1.0,7);
            //      break;
            //    case 3: 
            //      rhoGv2[par][irho][i][ikin][ieta][ipt] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),rerho1n1GfromHy,-1.0,1.0,7);
            //      break;
            //    case 4: 
            //      rhoGv2[par][irho][i][ikin][ieta][ipt] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),imrho1n1GfromHy,-1.0,1.0,7);
            //      break;
            //  }

            //  rhoGv2[par][irho][i][ikin][ieta][ipt]->SetParameter(0,0.2);      
            //  rhoGv2[par][irho][i][ikin][ieta][ipt]->SetParameter(1,double(ptkin[ikin]));      
            //  for(int ip = 0; ip < 5; ip++)
            //  {
            //    if(ip == irho)
            //    {
            //      if(ip == 0)
            //      {
            //        rhoGv2[par][irho][i][ikin][ieta][ipt]->SetParameter(2+ip,(double(i)-2)*0.1+1./3.);      
            //      }
            //      else
            //      {
            //        rhoGv2[par][irho][i][ikin][ieta][ipt]->SetParameter(2+ip,(double(i)-2)*0.1);       
            //      }
            //    }
            //    else
            //    {
            //      if(ip == 0)
            //      {
            //        rhoGv2[par][irho][i][ikin][ieta][ipt]->SetParameter(2+ip,1./3.);       
            //      }
            //      else
            //      {
            //        rhoGv2[par][irho][i][ikin][ieta][ipt]->SetParameter(2+ip,0.0);       
            //      }
            //    }
            //  }          
            //  rhoGv2[par][irho][i][ikin][ieta][ipt]->SetLineColor(color[i]);
            //  rhoGv2[par][irho][i][ikin][ieta][ipt]->SetLineStyle(1);
            //  rhoGv2[par][irho][i][ikin][ieta][ipt]->SetLineWidth(1);

              g2D_Global[irho][i][ikin][par][iopt][ieta]->SetMarkerStyle(20);
              g2D_Global[irho][i][ikin][par][iopt][ieta]->SetMarkerSize(1.0);
              g2D_Global[irho][i][ikin][par][iopt][ieta]->SetMarkerColor(color[i]);
              g2D_Global[irho][i][ikin][par][iopt][ieta]->SetLineColor(color[i]);
              //g2D_Global[irho][i][ikin][par][iet[ieta]a]->SetTitle(Form("Input %s^{h}, v_{2}=0.2, p_{T}=%1.1f GeV/c, |#eta|#leq%1.1f",param[irho].c_str(),ptkin[ikin],etacuts[ieta]));
              g2D_Global[irho][i][ikin][par][iopt][ieta]->SetTitle(Form("Input %s^{g}, v_{2}=0.2, p_{T}=%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),ptkin[ikin],etacuts[ieta]));
              g2D_Global[irho][i][ikin][par][iopt][ieta]->GetXaxis()->SetTitle(Form("(GeV/c) p_{T,K^{+/-}}>"));
              g2D_Global[irho][i][ikin][par][iopt][ieta]->GetYaxis()->SetTitle(Form("%s^{g} from 2D Fit - Truth",param[par].c_str()));
              g2D_Global[irho][i][ikin][par][iopt][ieta]->GetXaxis()->SetRangeUser(-0.025,0.225);
              if(min < 0.0 && max < 0.0) g2D_Global[irho][i][ikin][par][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,0.4*max);
              if(min < 0.0 && max > 0.0) g2D_Global[irho][i][ikin][par][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,1.6*max);
              if(min > 0.0 && max > 0.0) g2D_Global[irho][i][ikin][par][iopt][ieta]->GetYaxis()->SetRangeUser(0.4*min,1.6*max);
              //g2D_Global[irho][i][ikin][par]->GetXaxis()->SetLimits(-1.0,1.0);
              //g2D_Global[irho][i][ikin][par]->GetXaxis()->SetRangeUser(-1.0,1.0);

              if(i == 0) g2D_Global[irho][i][ikin][par][iopt][ieta]->Draw("APE"); 
              else       g2D_Global[irho][i][ikin][par][iopt][ieta]->Draw("PE same"); 

              if(ikin == 0 && iopt == 0 && ieta == 0 && par == 0) leg[irho]->AddEntry(g2D_Global[irho][i][ikin][par][iopt][ieta],Form("%s^{g} = %1.1f",param[irho].c_str(),(float(i)-2)*0.1),"p");
            }
            leg[irho]->Draw("same");
          }
        }
        cfit2D->Print(outputname.c_str());
        cfit2D->Update();
      }   
    }
  }
   
  cfit2D->Print(outputstop.c_str());



  TCanvas *cfit1D = new TCanvas("cfit1D","cfit1D",10,10,2000,400);
  cfit1D->Divide(5,1);
  
  for(int i = 0; i < 5; i++)
  {
    cfit1D->cd(i+1);
    cfit1D->cd(i+1)->SetLeftMargin(0.175);  
    cfit1D->cd(i+1)->SetBottomMargin(0.175);
    cfit1D->cd(i+1)->SetTicks(1,1);
    cfit1D->cd(i+1)->SetGrid(0,0);
  }
 

  outputname = Form("figures/%s/%s/pTstudy/Global1DvsGlobal_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  cfit1D->Print(outputstart.c_str());

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int ieta = 0; ieta < 11; ieta+=10)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          cfit1D->cd(1+irho);
          double max = -999999;
          double min = 999999;
          for(int i = 0; i < nvar; i++)
          {
            double tmax = g1D_Global[irho][i][ikin][iopt][ieta]->GetHistogram()->GetMaximum();   
            double tmin = g1D_Global[irho][i][ikin][iopt][ieta]->GetHistogram()->GetMinimum(); 

            if(tmax > max) max = tmax;
            if(tmin < min) min = tmin;
          }
          
          for(int i = 0; i < nvar; i++)
          {
            g1D_Global[irho][i][ikin][iopt][ieta]->SetMarkerStyle(20);
            g1D_Global[irho][i][ikin][iopt][ieta]->SetMarkerSize(1.0);
            g1D_Global[irho][i][ikin][iopt][ieta]->SetMarkerColor(color[i]);
            g1D_Global[irho][i][ikin][iopt][ieta]->SetLineColor(color[i]);
            g1D_Global[irho][i][ikin][iopt][ieta]->SetTitle(Form("Input %s^{g}, v_{2}=0.0, p_{T}=%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),ptkin[ikin],etacuts[ieta]));
            g1D_Global[irho][i][ikin][iopt][ieta]->GetXaxis()->SetTitle(Form("(GeV/c) p_{T,K^{+/-}}>"));
            g1D_Global[irho][i][ikin][iopt][ieta]->GetYaxis()->SetTitle(Form("%s^{g} from 1D Fit - Truth",param[0].c_str()));
            g1D_Global[irho][i][ikin][iopt][ieta]->GetXaxis()->SetRangeUser(-0.025,0.225);
            if(min < 0.0 && max < 0.0) g1D_Global[irho][i][ikin][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,0.4*max);
            if(min < 0.0 && max > 0.0) g1D_Global[irho][i][ikin][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,1.6*max);
            if(min > 0.0 && max > 0.0) g1D_Global[irho][i][ikin][iopt][ieta]->GetYaxis()->SetRangeUser(0.4*min,1.6*max);

            if(i == 0) g1D_Global[irho][i][ikin][iopt][ieta]->Draw("APE"); 
            else       g1D_Global[irho][i][ikin][iopt][ieta]->Draw("PE same"); 

          }
          leg[irho]->Draw("same");
        }   
        cfit1D->Print(outputname.c_str());
        cfit1D->Update();
      }
    }
  }
   
  cfit1D->Print(outputstop.c_str());

  outputname = Form("figures/%s/%s/pTstudy/Global1DCorrectedvsGlobal_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  cfit1D->Print(outputstart.c_str());

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int ieta = 0; ieta < 11; ieta+=10)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          cfit1D->cd(1+irho);
          double max = -999999;
          double min = 999999;
          for(int i = 0; i < nvar; i++)
          {
            double tmax = g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->GetHistogram()->GetMaximum();   
            double tmin = g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->GetHistogram()->GetMinimum(); 

            if(tmax > max) max = tmax;
            if(tmin < min) min = tmin;
          }
          
          for(int i = 0; i < nvar; i++)
          {
            g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->SetMarkerStyle(20);
            g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->SetMarkerSize(1.0);
            g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->SetMarkerColor(color[i]);
            g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->SetLineColor(color[i]);
            g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->SetTitle(Form("Input %s^{g}, v_{2}=0.0, p_{T}=%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),ptkin[ikin],etacuts[ieta]));
            g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->GetXaxis()->SetTitle(Form("(GeV/c) p_{T,K^{+/-}}>"));
            g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->GetYaxis()->SetTitle(Form("%s^{g} from 1D Fit - Truth",param[0].c_str()));
            g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->GetXaxis()->SetRangeUser(-0.05,0.225);
            if(min < 0.0 && max < 0.0) g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,0.4*max);
            if(min < 0.0 && max > 0.0) g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,1.6*max);
            if(min > 0.0 && max > 0.0) g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->GetYaxis()->SetRangeUser(0.4*min,1.6*max);

            if(i == 0) g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->Draw("APE"); 
            else       g1D_Global_Corrected[irho][i][ikin][iopt][ieta]->Draw("PE same"); 

          }
          leg[irho]->Draw("same");
        }   
        cfit1D->Print(outputname.c_str());
        cfit1D->Update();
      }
    }
  }
   
  cfit1D->Print(outputstop.c_str());








  outputname = Form("figures/%s/%s/pTstudy/HelicityvsGlobal_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  cfit2D->Print(outputstart.c_str());


  //TLegend *leg;// = new TLegend(0.4,0.75,0.7,0.9);
  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int ieta = 0; ieta < 11; ieta+=10)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          //double max = -999999;
          //double min = 999999;
          //for(int par = 0; par < 5; par++)
          //{
          //  for(int i = 0; i < nvar; i++)
          //  {
          //    double tmax = g2D[irho][i][ikin][par][iopt]->GetHistogram()->GetMaximum();   
          //    double tmin = g2D[irho][i][ikin][par][iopt]->GetHistogram()->GetMinimum(); 

          //    if(tmax > max) max = tmax;
          //    if(tmin < min) min = tmin;

          //    //cout << "max = " << max << ", min = " << min << endl;
          //  }
          //}
          for(int par = 0; par < 5; par++)
          {
            cfit2D->cd(1+irho+5*par);
            //cfit2D->cd(1+irho+5*par);
            //TLegend *leg = new TLegend(0.4,0.6,0.6,0.8);
            double max = -999999;
            double min = 999999;
            for(int i = 0; i < nvar; i++)
            {
              double tmax = g2D[irho][i][ikin][par][iopt][ieta]->GetHistogram()->GetMaximum();   
              double tmin = g2D[irho][i][ikin][par][iopt][ieta]->GetHistogram()->GetMinimum(); 

              if(tmax > max) max = tmax;
              if(tmin < min) min = tmin;

              //cout << "max = " << max << ", min = " << min << endl;
            }
            for(int i = 0; i < nvar; i++)
            {
              g2D[irho][i][ikin][par][iopt][ieta]->SetMarkerStyle(20);
              g2D[irho][i][ikin][par][iopt][ieta]->SetMarkerSize(1.2);
              g2D[irho][i][ikin][par][iopt][ieta]->SetMarkerColor(color[i]);
              g2D[irho][i][ikin][par][iopt][ieta]->SetLineColor(color[i]);
              g2D[irho][i][ikin][par][iopt][ieta]->SetTitle(Form("Input %s^{g}, v_{2}=0.2, p_{T}=%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),ptkin[ikin],etacuts[ieta]));
              g2D[irho][i][ikin][par][iopt][ieta]->GetXaxis()->SetTitle(Form("(GeV/c) p_{T,K^{+/-}}>"));
              g2D[irho][i][ikin][par][iopt][ieta]->GetYaxis()->SetTitle(Form("%s^{Helicity} from 2D Fit",param[par].c_str()));
              if(min < 0.0 && max < 0.0) g2D[irho][i][ikin][par][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,0.4*max);
              if(min < 0.0 && max > 0.0) g2D[irho][i][ikin][par][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,1.6*max);
              if(min > 0.0 && max > 0.0) g2D[irho][i][ikin][par][iopt][ieta]->GetYaxis()->SetRangeUser(0.4*min,1.6*max);
              g2D[irho][i][ikin][par][iopt][ieta]->GetXaxis()->SetRangeUser(0.025,0.225);

              if(i == 0) g2D[irho][i][ikin][par][iopt][ieta]->Draw("APE"); 
              else       g2D[irho][i][ikin][par][iopt][ieta]->Draw("PE same"); 
                
              //if(par == 0) leg->AddEntry(g2D[irho][i][ikin][par],Form("%s^{Helicity} = %1.1f",param[inputpar].c_str(),(float(i)-2)*0.1),"p");
            }
            leg[irho]->Draw("same");
          }
        }
        cfit2D->Print(outputname.c_str());
        cfit2D->Update();
      }
    }
  }
  cfit2D->Print(outputstop.c_str());

  TCanvas *c2D = new TCanvas("c2D","c2D",10,10,1600,1200);
  c2D->Divide(4,3);
  
  for(int i = 0; i < 12; i++)
  {
    c2D->cd(i+1);
    c2D->cd(i+1)->SetLeftMargin(0.15);  
    c2D->cd(i+1)->SetRightMargin(0.15);  
    c2D->cd(i+1)->SetBottomMargin(0.15);
    c2D->cd(i+1)->SetTicks(1,1);
    c2D->cd(i+1)->SetGrid(0,0);
  }
 
  outputname = Form("figures/%s/%s/pTstudy/Global2DDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  c2D->Print(outputstart.c_str());

  for(int ipt = 0; ipt < 9; ipt++)
  {
    for(int ieta = 0; ieta < 11; ieta+=10)
    {
      for(int ikin = 0; ikin < nkin; ikin++)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          for(int i = 0; i < nvar; i++)
          //for(int i = 0; i < nvar; i++)
          {
            for(int iopt = 0; iopt < nopts; iopt++)
            {
              c2D->cd(1+iopt);
              h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{h}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta]));
              h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
              h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{g}");
              h_m2D_RC_Global[iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
            }
            c2D->Print(outputname.c_str());
            c2D->Update();
          }
        }
      }
    }
  }
  c2D->Print(outputstop.c_str());

  outputname = Form("figures/%s/%s/pTstudy/Helicity2DDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  c2D->Print(outputstart.c_str());

  for(int ipt = 0; ipt < 9; ipt++)
  {
    for(int ieta = 0; ieta < 11; ieta+=10)
    {
      for(int ikin = 0; ikin < nkin; ikin++)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          for(int i = 0; i < nvar; i++)
          {
            for(int iopt = 0; iopt < nopts; iopt++)
            {
              c2D->cd(1+iopt);
              h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{h}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta]));
              h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{h}");
              h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{h}");
              h_m2D_RC[iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
            }
            c2D->Print(outputname.c_str());
            c2D->Update();
          }
        }
      }
    }
  }
  c2D->Print(outputstop.c_str());


  //TF2* fitacc2D[nopts][5][nvar][nkin][10][9];
  //TGraphAsymmErrors *g_m2D_chi2ndf[5][nvar][nkin][10][9];
  //for(int ieta = 9; ieta >= 0; ieta--)
  //{
  //  for(int ikin = 0; ikin < nkin; ikin++)
  //  {
  //    for(int irho = 0; irho < 5; irho++)
  //    {
  //      for(int i = 0; i < nvar; i++)
  //      //for(int i = 0; i < nvar; i++)
  //      {
  //        g_m2D_chi2ndf[irho][i][ikin][ieta][ipt] = new TGraphAsymmErrors();
  //        for(int iopt = 0; iopt < nopts; iopt++)
  //        {
  //          cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;

  //          fitacc2D[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
  //          //fitacc2D[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple7,-1.0,1.0,0.0,2.0*TMath::Pi(),12);
  //          //fitacc2D[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple8,-1.0,1.0,0.0,2.0*TMath::Pi(),16);
  //          fitacc2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
  //          
  //          if(ieta < 9)
  //          {
  //            for(int ipar = 1; ipar < 9; ipar++) fitacc2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(ipar,fitacc2D[iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
  //            //for(int ipar = 1; ipar < 12; ipar++) fitacc2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(ipar,fitacc2D[iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
  //            //for(int ipar = 1; ipar < 16; ipar++) fitacc2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(ipar,fitacc2D[iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
  //          }
  //          //h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[iopt][irho][i][ikin][ieta][ipt],"NMR");
 
  //          double chi2 = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->GetChisquare();
  //          double ndf  = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->GetNDF();

  //          cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
  //          g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->SetPoint(iopt,vals[iopt],chi2/ndf);

  //          for(int iy = 1; iy <= 20; iy++)
  //          { 
  //            g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][iy-1] = new TGraphAsymmErrors();
  //          }
  //          for(int ix = 1; ix <= 20; ix++)
  //          {
  //            g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ix-1] = new TGraphAsymmErrors();
  //            for(int iy = 1; iy <= 20; iy++)
  //            { 
  //              double widthx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinWidth(ix);
  //              double startx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinLowEdge(ix);
  //              double stopx  = startx+widthx;
  //              double widthy = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinWidth(iy);
  //              double starty = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinLowEdge(iy);
  //              double stopy  = starty+widthy;
  //              double value = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
  //              
  //              double bincenterx = (startx+stopx)/2.0;
  //              double bincentery = (starty+stopy)/2.0;
  //  
  //              if(h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetBinContent(ix,iy) <= 0.0) 
  //              cout << "RATIO IS "<<h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetBinContent(ix,iy)<<"+/-"<<h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetBinError(ix,iy)<<"!?   iopt="<<iopt<<",irho="<<irho<<",i="<<i<<",ikin"<<ikin<<",ieta"<<ieta<<endl;

  //              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][iy-1]->SetPoint(ix-1,bincenterx,value); 
  //              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ix-1]->SetPoint(iy-1,bincentery,value); 

  //            }
  //          }
  //        }
  //      }
  //    }
  //  }
  //}

  //TF2* fitacc2DH[nopts][5][nvar][nkin][10][9];
  //TGraphAsymmErrors *g_m2DH_chi2ndf[5][nvar][nkin][10][9];
  //for(int ieta = 9; ieta >= 0; ieta--)
  //{
  //  for(int ikin = 0; ikin < nkin; ikin++)
  //  {
  //    for(int irho = 0; irho < 5; irho++)
  //    {
  //      for(int i = 0; i < nvar; i++)
  //      //for(int i = 0; i < nvar; i++)
  //      {
  //        g_m2DH_chi2ndf[irho][i][ikin][ieta][ipt] = new TGraphAsymmErrors();
  //        for(int iopt = 0; iopt < nopts; iopt++)
  //        {
  //          cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;

  //          fitacc2DH[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
  //          //fitacc2D[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple7,-1.0,1.0,0.0,2.0*TMath::Pi(),12);
  //          fitacc2DH[iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
  //          
  //          if(ieta < 9)
  //          {
  //            for(int ipar = 1; ipar < 9; ipar++) fitacc2DH[iopt][irho][i][ikin][ieta][ipt]->SetParameter(ipar,fitacc2DH[iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
  //            //for(int ipar = 1; ipar < 12; ipar++) fitacc2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(ipar,fitacc2D[iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
  //          }
  //          //h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2DH[iopt][irho][i][ikin][ieta][ipt],"NMR");

  //          double chi2 = fitacc2DH[iopt][irho][i][ikin][ieta][ipt]->GetChisquare();
  //          double ndf  = fitacc2DH[iopt][irho][i][ikin][ieta][ipt]->GetNDF();

  //          cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
  //          g_m2DH_chi2ndf[irho][i][ikin][ieta][ipt]->SetPoint(iopt,vals[iopt],chi2/ndf);

  //          for(int iy = 1; iy <= 20; iy++)
  //          { 
  //            g_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][iy-1] = new TGraphAsymmErrors();
  //          }
  //          for(int ix = 1; ix <= 20; ix++)
  //          {
  //            g_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ix-1] = new TGraphAsymmErrors();
  //            for(int iy = 1; iy <= 20; iy++)
  //            { 
  //              double widthx = h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinWidth(ix);
  //              double startx = h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinLowEdge(ix);
  //              double stopx  = startx+widthx;
  //              double widthy = h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinWidth(iy);
  //              double starty = h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinLowEdge(iy);
  //              double stopy  = starty+widthy;
  //              double value = fitacc2DH[iopt][irho][i][ikin][ieta][ipt]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                

  //              double bincenterx = (startx+stopx)/2.0;
  //              double bincentery = (starty+stopy)/2.0;

  //              g_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][iy-1]->SetPoint(ix-1,bincenterx,value); 
  //              g_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ix-1]->SetPoint(iy-1,bincentery,value); 

  //            }
  //          }
  //        }
  //      }
  //    }
  //  }
  //}

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int i = 0; i < nvar; i++)
    {
      for(int ikin = 0; ikin < nkin; ikin++)
      {
        for(int ieta = 0; ieta < 11; ieta+=10) 
        {
          for(int irho = 0; irho < 5; irho++)
          {
            for(int par = 0; par < 5; par++) 
            {
              g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta] = new TGraphAsymmErrors();
            }   
            for(int ipt = 0; ipt < 9; ipt++)
            {
              
              fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta,ipt),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),6);
              //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2DcosAcc,-1.0,1.0,0.0,2.0*TMath::Pi(),14);
              //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2DcosAccLonger,-1.0,1.0,0.0,2.0*TMath::Pi(),17);
              //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2DcosAccLongerer,-1.0,1.0,0.0,2.0*TMath::Pi(),21);
              fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,1./3.);
              fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetParLimits(0,0.0,1.0);
              fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetParameter(1,0.0);
              fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetParameter(2,0.0);
              fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetParameter(3,0.0);
              fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetParameter(4,0.0);
              fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetParameter(5,h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
              //for(int ifit = 1; ifit < 9; ifit++)
              ////for(int ifit = 1; ifit < 12; ifit++)
              ////for(int ifit = 1; ifit < 16; ifit++)
              //{
              //  double value = fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(ifit);
              //  double error = fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParError(ifit);
              //  fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->FixParameter(5+ifit,value);
              //  //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetParLimits(5+ifit,value-error,value+error);
              //}
              //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->FixParameter(7,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(2));
              //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->FixParameter(8,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(3));
              //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->FixParameter(9,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(4));
              //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->FixParameter(10,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(5));
              //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->FixParameter(11,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(6));
              //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->FixParameter(12,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(7));
              //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->FixParameter(13,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(8));
             
              h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->Fit(fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt],"NMRI");
              for(int par = 0; par < 5; par++)
              {
                if(par == 0 && par == irho) g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025,fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
                if(par != 0 && par == irho) g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025,fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetParameter(par)-(float(i)-2)*0.1);
                if(par == 0 && par != irho) g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025,fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetParameter(0)-1./3.);
                if(par != 0 && par != irho) g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025,fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetParameter(par));
                g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->SetPointError(ipt,0.0,0.0,fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetParError(par),fits2D_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetParError(par));
              }
            }
          }
        }
      }
    }
  }

  outputname = Form("figures/%s/%s/pTstudy/GlobalCorrectedvsGlobal_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  cfit2D->Print(outputstart.c_str());

  //TF1 *rhoGv2[5][5][nvar][nkin][11][9];

  //TLegend *leg[5];// = new TLegend(0.4,0.75,0.7,0.9);

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int ieta = 0; ieta < 11; ieta+=10)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          //double max = -999999;
          //double min = 999999;
          //for(int i = 0; i < nvar; i++)
          //{
          //  for(int par = 0; par < 5; par++)
          //  {
          //    double tmax = g2D_Global[irho][i][ikin][par][iopt]->GetHistogram()->GetMaximum();   
          //    double tmin = g2D_Global[irho][i][ikin][par][iopt]->GetHistogram()->GetMinimum(); 

          //    if(tmax > max) max = tmax;
          //    if(tmin < min) min = tmin;
          //  }
          //}
          for(int par = 0; par < 5; par++)
          {
            //cfit2D->cd(1+irho+5*par);
            cfit2D->cd(1+irho+5*par);
            double max = -999999;
            double min = 999999;
            for(int i = 0; i < nvar; i++)
            {
                double tmax = g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->GetHistogram()->GetMaximum();   
                double tmin = g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->GetHistogram()->GetMinimum(); 

                if(tmax > max) max = tmax;
                if(tmin < min) min = tmin;
            }
            
            //double max = -999999;
            //double min = 999999;
            //for(int i = 0; i < nvar; i++)
            //{
            //  for(int ir = 0; ir < 5; ir++)
            //  {
            //    double tmax = g2D_Global[ir][i][ikin][par]->GetHistogram()->GetMaximum();   
            //    double tmin = g2D_Global[ir][i][ikin][par]->GetHistogram()->GetMinimum(); 

            //    if(tmax > max) max = tmax;
            //    if(tmin < min) min = tmin;
            //  }
            //}
            for(int i = 0; i < nvar; i++)
            {
            //  switch(par)
            //  {
            //    case 0: 
            //      rhoGv2[par][irho][i][ikin][ieta][ipt] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),rhoGfromHy,-1.0,1.0,7);
            //      break;
            //    case 1: 
            //      rhoGv2[par][irho][i][ikin][ieta][ipt] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),realGfromHy,-1.0,1.0,7);
            //      break;
            //    case 2: 
            //      rhoGv2[par][irho][i][ikin][ieta][ipt] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),imagGfromHy,-1.0,1.0,7);
            //      break;
            //    case 3: 
            //      rhoGv2[par][irho][i][ikin][ieta][ipt] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),rerho1n1GfromHy,-1.0,1.0,7);
            //      break;
            //    case 4: 
            //      rhoGv2[par][irho][i][ikin][ieta][ipt] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),imrho1n1GfromHy,-1.0,1.0,7);
            //      break;
            //  }

            //  rhoGv2[par][irho][i][ikin][ieta][ipt]->SetParameter(0,0.2);      
            //  rhoGv2[par][irho][i][ikin][ieta][ipt]->SetParameter(1,double(ptkin[ikin]));      
            //  for(int ip = 0; ip < 5; ip++)
            //  {
            //    if(ip == irho)
            //    {
            //      if(ip == 0)
            //      {
            //        rhoGv2[par][irho][i][ikin][ieta][ipt]->SetParameter(2+ip,(double(i)-2)*0.1+1./3.);      
            //      }
            //      else
            //      {
            //        rhoGv2[par][irho][i][ikin][ieta][ipt]->SetParameter(2+ip,(double(i)-2)*0.1);       
            //      }
            //    }
            //    else
            //    {
            //      if(ip == 0)
            //      {
            //        rhoGv2[par][irho][i][ikin][ieta][ipt]->SetParameter(2+ip,1./3.);       
            //      }
            //      else
            //      {
            //        rhoGv2[par][irho][i][ikin][ieta][ipt]->SetParameter(2+ip,0.0);       
            //      }
            //    }
            //  }          
            //  rhoGv2[par][irho][i][ikin][ieta][ipt]->SetLineColor(color[i]);
            //  rhoGv2[par][irho][i][ikin][ieta][ipt]->SetLineStyle(1);
            //  rhoGv2[par][irho][i][ikin][ieta][ipt]->SetLineWidth(1);

              g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->SetMarkerStyle(20);
              g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->SetMarkerSize(1.0);
              g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->SetMarkerColor(color[i]);
              g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->SetLineColor(color[i]);
              g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->SetTitle(Form("Input %s^{g}, v_{2}=0.2, p_{T}=%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),ptkin[ikin],etacuts[ieta]));
              g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->GetXaxis()->SetTitle(Form("(GeV/c) p_{T,K^{+/-}}<"));
              g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->GetYaxis()->SetTitle(Form("Corrected %s^{g} from 2D Fit - Truth",param[par].c_str()));
              g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->GetXaxis()->SetRangeUser(-0.025,0.225);
              if(min < 0.0 && max < 0.0) g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,0.4*max);
              if(min < 0.0 && max > 0.0) g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,1.6*max);
              if(min > 0.0 && max > 0.0) g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->GetYaxis()->SetRangeUser(0.4*min,1.6*max);
              //g2D_Global[irho][i][ikin][par]->GetXaxis()->SetLimits(-1.0,1.0);
              //g2D_Global[irho][i][ikin][par]->GetXaxis()->SetRangeUser(-1.0,1.0);

              if(i == 0) g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->Draw("APE"); 
              else       g2D_Global_Corrected[irho][i][ikin][par][iopt][ieta]->Draw("PE same"); 

              //rhoGv2[par][irho][i][ikin][ieta][ipt]->Draw("l same");
                
              //if(par == 0 && ikin == 0 && iopt == 0) leg[irho]->AddEntry(g2D_Global[irho][i][ikin][par][iopt],Form("%s^{h} = %1.1f",param[irho].c_str(),(float(i)-2)*0.1),"p");
            }
            leg[irho]->Draw("same");
          }
        }   
        cfit2D->Print(outputname.c_str());
        cfit2D->Update();
      }
    }
  }
   
  cfit2D->Print(outputstop.c_str());

//  outputname = Form("figures/%s/%s/pTstudy/Global2DRatioDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c2D->Print(outputstart.c_str());
//
//  for(int ieta = 0; ieta < 10; ieta+=10)
//  {
//    for(int ikin = 0; ikin < nkin; ikin++)
//    {
//      for(int irho = 0; irho < 5; irho++)
//      {
//        for(int i = 0; i < nvar; i++)
//        //for(int i = 0; i < nvar; i++)
//        {
//          //g_m2D_chi2ndf[irho][i][ikin][ieta][ipt] = new TGraphAsymmErrors();
//          for(int iopt = 0; iopt < nopts; iopt++)
//          {
//            cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;
//            c2D->cd(1+iopt);
//
//            //fitacc2D[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
//            //fitacc2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
//
//            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1]));
//            cout << "Set title" << endl;
//            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
//            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{g}");
//            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
//
////            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[iopt][irho][i][ikin][ieta][ipt],"NMR");
////            //h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[iopt][irho][i][ikin][ieta][ipt],"NMR");
////            double chi2 = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->GetChisquare();
////            double ndf  = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->GetNDF();
////            cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
////            g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->SetPoint(iopt,vals[iopt],chi2/ndf);
////
////            for(int iy = 1; iy <= 20; iy++)
////            { 
////              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][iy-1] = new TGraphAsymmErrors();
////            }
////            for(int ix = 1; ix <= 20; ix++)
////            {
////              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ix-1] = new TGraphAsymmErrors();
////              for(int iy = 1; iy <= 20; iy++)
////              { 
////                double widthx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinWidth(ix);
////                double startx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinLowEdge(ix);
////                double stopx  = startx+widthx;
////                double widthy = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinWidth(iy);
////                double starty = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinLowEdge(iy);
////                double stopy  = starty+widthy;
////                double value = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
////
////                double bincenterx = (startx+stopx)/2.0;
////                double bincentery = (starty+stopy)/2.0;
////
////                g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][iy-1]->SetPoint(ix-1,bincenterx,value); 
////                g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ix-1]->SetPoint(iy-1,bincentery,value); 
////
////              }
////            }
//            
//          }
//          c2D->Print(outputname.c_str());
//          c2D->Update();
//        }
//      }
//    }
//  }
//  c2D->Print(outputstop.c_str());
//
//  outputname = Form("figures/%s/%s/pTstudy/Global2DCorrectedDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c2D->Print(outputstart.c_str());
//
//  for(int ieta = 0; ieta < 10; ieta+=10)
//  {
//    for(int ikin = 0; ikin < nkin; ikin++)
//    {
//      for(int irho = 0; irho < 5; irho++)
//      {
//        for(int i = 0; i < nvar; i++)
//        //for(int i = 0; i < nvar; i++)
//        {
//          //g_m2D_chi2ndf[irho][i][ikin][ieta][ipt] = new TGraphAsymmErrors();
//          for(int iopt = 0; iopt < nopts; iopt++)
//          {
//            cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;
//            c2D->cd(1+iopt);
//
//            //fitacc2D[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
//            //fitacc2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
//
//            h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1]));
//            h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
//            h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{g}");
//            h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
//
////            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[iopt][irho][i][ikin][ieta][ipt],"NMR");
////            //h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[iopt][irho][i][ikin][ieta][ipt],"NMR");
////            double chi2 = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->GetChisquare();
////            double ndf  = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->GetNDF();
////            cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
////            g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->SetPoint(iopt,vals[iopt],chi2/ndf);
////
////            for(int iy = 1; iy <= 20; iy++)
////            { 
////              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][iy-1] = new TGraphAsymmErrors();
////            }
////            for(int ix = 1; ix <= 20; ix++)
////            {
////              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ix-1] = new TGraphAsymmErrors();
////              for(int iy = 1; iy <= 20; iy++)
////              { 
////                double widthx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinWidth(ix);
////                double startx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinLowEdge(ix);
////                double stopx  = startx+widthx;
////                double widthy = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinWidth(iy);
////                double starty = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinLowEdge(iy);
////                double stopy  = starty+widthy;
////                double value = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
////
////                double bincenterx = (startx+stopx)/2.0;
////                double bincentery = (starty+stopy)/2.0;
////
////                g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][iy-1]->SetPoint(ix-1,bincenterx,value); 
////                g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ix-1]->SetPoint(iy-1,bincentery,value); 
////
////              }
////            }
//            
//          }
//          c2D->Print(outputname.c_str());
//          c2D->Update();
//        }
//      }
//    }
//  }
//  c2D->Print(outputstop.c_str());
//  outputname = Form("figures/%s/%s/pTstudy/Global2DCorrectedDistributionsOverMC_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c2D->Print(outputstart.c_str());
//
//  for(int ieta = 0; ieta < 10; ieta+=10)
//  {
//    for(int ikin = 0; ikin < nkin; ikin++)
//    {
//      for(int irho = 0; irho < 5; irho++)
//      {
//        for(int i = 0; i < nvar; i++)
//        //for(int i = 0; i < nvar; i++)
//        {
//          //g_m2D_chi2ndf[irho][i][ikin][ieta][ipt] = new TGraphAsymmErrors();
//          for(int iopt = 0; iopt < nopts; iopt++)
//          {
//            cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;
//            c2D->cd(1+iopt);
//
//            //fitacc2D[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
//            //fitacc2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
//
//            h_m2D_RC_Global_Corrected_BackToMC[iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1]));
//            h_m2D_RC_Global_Corrected_BackToMC[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
//            h_m2D_RC_Global_Corrected_BackToMC[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{g}");
//            h_m2D_RC_Global_Corrected_BackToMC[iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
//
//            h_m2D_RC_Global_Corrected_BackToMC_E[iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC_Global_Corrected_BackToMC[iopt][irho][i][ikin][ieta][ipt]->Clone();
//            for(int ix = 1; ix <= 20; ix++)
//            {
//              for(int iy = 1; iy <= 20; iy++)
//              { 
//                double value = h_m2D_RC_Global_Corrected_BackToMC[iopt][irho][i][ikin][ieta][ipt]->GetBinError(ix,iy);                
//                cout << "CORRECTED/MC = " << std::fixed << std::setprecision(15) << h_m2D_RC_Global_Corrected_BackToMC[iopt][irho][i][ikin][ieta][ipt]->GetBinContent(ix,iy) << endl;
//                h_m2D_RC_Global_Corrected_BackToMC_E[iopt][irho][i][ikin][ieta][ipt]->SetBinContent(ix,iy,value);
//              }
//            }
//
////            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[iopt][irho][i][ikin][ieta][ipt],"NMR");
////            //h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[iopt][irho][i][ikin][ieta][ipt],"NMR");
////            double chi2 = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->GetChisquare();
////            double ndf  = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->GetNDF();
////            cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
////            g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->SetPoint(iopt,vals[iopt],chi2/ndf);
////
////            for(int iy = 1; iy <= 20; iy++)
////            { 
////              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][iy-1] = new TGraphAsymmErrors();
////            }
////            for(int ix = 1; ix <= 20; ix++)
////            {
////              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ix-1] = new TGraphAsymmErrors();
////              for(int iy = 1; iy <= 20; iy++)
////              { 
////                double widthx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinWidth(ix);
////                double startx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinLowEdge(ix);
////                double stopx  = startx+widthx;
////                double widthy = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinWidth(iy);
////                double starty = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinLowEdge(iy);
////                double stopy  = starty+widthy;
////                double value = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
////
////                double bincenterx = (startx+stopx)/2.0;
////                double bincentery = (starty+stopy)/2.0;
////
////                g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][iy-1]->SetPoint(ix-1,bincenterx,value); 
////                g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ix-1]->SetPoint(iy-1,bincentery,value); 
////
////              }
////            }
//            
//          }
//          c2D->Print(outputname.c_str());
//          c2D->Update();
//        }
//      }
//    }
//  }
//  c2D->Print(outputstop.c_str());
//  outputname = Form("figures/%s/%s/pTstudy/Global2DCorrectedDistributionsOverMCError_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c2D->Print(outputstart.c_str());
//
//  for(int ieta = 0; ieta < 10; ieta+=10)
//  {
//    for(int ikin = 0; ikin < nkin; ikin++)
//    {
//      for(int irho = 0; irho < 5; irho++)
//      {
//        for(int i = 0; i < nvar; i++)
//        //for(int i = 0; i < nvar; i++)
//        {
//          //g_m2D_chi2ndf[irho][i][ikin][ieta][ipt] = new TGraphAsymmErrors();
//          for(int iopt = 0; iopt < nopts; iopt++)
//          {
//            cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;
//            c2D->cd(1+iopt);
//
//            //fitacc2D[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
//            //fitacc2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
//
//            h_m2D_RC_Global_Corrected_BackToMC_E[iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1]));
//            h_m2D_RC_Global_Corrected_BackToMC_E[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
//            h_m2D_RC_Global_Corrected_BackToMC_E[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{g}");
//            h_m2D_RC_Global_Corrected_BackToMC_E[iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
//
//          }
//          c2D->Print(outputname.c_str());
//          c2D->Update();
//        }
//      }
//    }
//  }
//  c2D->Print(outputstop.c_str());
//
//  outputname = Form("figures/%s/%s/pTstudy/Helicity2DCorrectedDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c2D->Print(outputstart.c_str());
//
//  for(int ieta = 0; ieta < 10; ieta+=10)
//  {
//    for(int ikin = 0; ikin < nkin; ikin++)
//    {
//      for(int irho = 0; irho < 5; irho++)
//      {
//        for(int i = 0; i < nvar; i++)
//        //for(int i = 0; i < nvar; i++)
//        {
//          //g_m2D_chi2ndf[irho][i][ikin][ieta][ipt] = new TGraphAsymmErrors();
//          for(int iopt = 0; iopt < nopts; iopt++)
//          {
//            cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;
//            c2D->cd(1+iopt);
//
//            //fitacc2D[iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
//            //fitacc2D[iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
//
//            h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1]));
//            h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{h}");
//            h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{h}");
//            h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
//
////            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[iopt][irho][i][ikin][ieta][ipt],"NMR");
////            //h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[iopt][irho][i][ikin][ieta][ipt],"NMR");
////            double chi2 = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->GetChisquare();
////            double ndf  = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->GetNDF();
////            cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
////            g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->SetPoint(iopt,vals[iopt],chi2/ndf);
////
////            for(int iy = 1; iy <= 20; iy++)
////            { 
////              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][iy-1] = new TGraphAsymmErrors();
////            }
////            for(int ix = 1; ix <= 20; ix++)
////            {
////              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ix-1] = new TGraphAsymmErrors();
////              for(int iy = 1; iy <= 20; iy++)
////              { 
////                double widthx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinWidth(ix);
////                double startx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinLowEdge(ix);
////                double stopx  = startx+widthx;
////                double widthy = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinWidth(iy);
////                double starty = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinLowEdge(iy);
////                double stopy  = starty+widthy;
////                double value = fitacc2D[iopt][irho][i][ikin][ieta][ipt]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
////
////                double bincenterx = (startx+stopx)/2.0;
////                double bincentery = (starty+stopy)/2.0;
////
////                g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][iy-1]->SetPoint(ix-1,bincenterx,value); 
////                g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ix-1]->SetPoint(iy-1,bincentery,value); 
////
////              }
////            }
//            
//          }
//          c2D->Print(outputname.c_str());
//          c2D->Update();
//        }
//      }
//    }
//  }
//  c2D->Print(outputstop.c_str());
//
//
//  outputname = Form("figures/%s/%s/pTstudy/Global2DRatioDistributionsCHI2NDF_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c2D->Print(outputstart.c_str());
//
//  for(int ikin = 0; ikin < nkin; ikin++)
//  {
//    for(int irho = 0; irho < 5; irho++)
//    {
//      for(int i = 0; i < nvar; i++)
//      //for(int i = 0; i < nvar; i++)
//      {
//        for(int ieta = 0; ieta < 12; ieta+=10)
//        {
//          //cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;
//          c2D->cd(1+ieta);
//          c2D->cd(1+ieta)->Clear();
//          if(ieta > 9/* || ieta < 2*/) continue;
//
//          g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,etacuts[ieta+1]));
//          g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("y");
//          g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#chi^2/NDF");
//          g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->Draw("APE");
//        }
//        c2D->Print(outputname.c_str());
//        c2D->Update();
//      }
//    }
//  }
//  c2D->Print(outputstop.c_str());
//
//
//  outputname = Form("figures/%s/%s/pTstudy/Helicity2DRatioDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c2D->Print(outputstart.c_str());
//
//  for(int ieta = 0; ieta < 10; ieta+=10)
//  {
//    for(int ikin = 0; ikin < nkin; ikin++)
//    {
//      for(int irho = 0; irho < 5; irho++)
//      {
//        //for(int i = 0; i < nvar; i++)
//        for(int i = 0; i < nvar; i++)
//        {
//          for(int iopt = 0; iopt < nopts; iopt++)
//          {
//            c2D->cd(1+iopt);
//            h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{h}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1]));
//            h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{h}");
//            h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{h}");
//            h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
//          }
//          c2D->Print(outputname.c_str());
//          c2D->Update();
//        }
//      }
//    }
//  }
//  c2D->Print(outputstop.c_str());
//
//  TCanvas *c1D = new TCanvas("c2D","c2D",10,10,2000,1600);
//  c1D->Divide(5,4);
//  
//  for(int i = 0; i < 20; i++)
//  {
//    c1D->cd(i+1);
//    c1D->cd(i+1)->SetLeftMargin(0.15);  
//    c1D->cd(i+1)->SetBottomMargin(0.15);
//    c1D->cd(i+1)->SetTicks(1,1);
//    c1D->cd(i+1)->SetGrid(0,0);
//  }

//  outputname = Form("figures/%s/%s/pTstudy/Global1DCosRatioDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c1D->Print(outputstart.c_str());
//
//  for(int ieta = 0; ieta < 10; ieta+=10)
//  {
//    for(int ikin = 0; ikin < nkin; ikin++)
//    {
//      for(int irho = 0; irho < 5; irho++)
//      {
//        for(int i = 0; i < nvar; i++)
//        //for(int i = 0; i < nvar; i++)
//        {
//          for(int iopt = 0; iopt < nopts; iopt++)
//          {
//            for(int ibin = 0; ibin < 20; ibin++)
//            {
//              c1D->cd(1+ibin);
//              c1D->cd(1+ibin)->SetLogy();
//              h_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f, #beta%d",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1],ibin));
//              h_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
//              h_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->GetYaxis()->SetTitle("Acceptance");
//              h_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("Colz");
//              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->SetLineColor(kRed);
//              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("C");
//            }
//            c1D->Print(outputname.c_str());
//            c1D->Update();
//          }
//        }
//      }
//    }
//  }
//  c1D->Print(outputstop.c_str());
//
//  outputname = Form("figures/%s/%s/pTstudy/Global1DBetaRatioDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c1D->Print(outputstart.c_str());
//
//  for(int ieta = 0; ieta < 10; ieta+=10)
//  {
//    for(int ikin = 0; ikin < nkin; ikin++)
//    {
//      for(int irho = 0; irho < 5; irho++)
//      {
//        for(int i = 0; i < nvar; i++)
//        //for(int i = 0; i < nvar; i++)
//        {
//          for(int iopt = 0; iopt < nopts; iopt++)
//          {
//            for(int ibin = 0; ibin < 20; ibin++)
//            {
//              c1D->cd(1+ibin);
//              h_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f, cos%d",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1],ibin));
//              h_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->GetXaxis()->SetTitle("#beta_{g}");
//              h_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->GetYaxis()->SetTitle("Acceptance");
//              h_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("Colz");
//              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->SetLineColor(kRed);
//              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("C");
//            }
//            c1D->Print(outputname.c_str());
//            c1D->Update();
//          }
//        }
//      }
//    }
//  }
//  c1D->Print(outputstop.c_str());
//
//
//  outputname = Form("figures/%s/%s/pTstudy/Helicity1DCosRatioDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c1D->Print(outputstart.c_str());
//
//  for(int ieta = 0; ieta < 10; ieta+=10)
//  {
//    for(int ikin = 0; ikin < nkin; ikin++)
//    {
//      for(int irho = 0; irho < 5; irho++)
//      {
//        for(int i = 0; i < nvar; i++)
//        //for(int i = 0; i < nvar; i++)
//        {
//          for(int iopt = 0; iopt < nopts; iopt++)
//          {
//            for(int ibin = 0; ibin < 20; ibin++)
//            {
//              c1D->cd(1+ibin);
//              c1D->cd(1+ibin)->SetLogy();
//              h_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f, #beta%d",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1],ibin));
//              h_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->GetXaxis()->SetTitle("cos(#theta^{*}_{h}");
//              h_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->GetYaxis()->SetTitle("Acceptance");
//              h_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("Colz");
//              g_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->SetLineColor(kRed);
//              g_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("C");
//            }
//            c1D->Print(outputname.c_str());
//            c1D->Update();
//          }
//        }
//      }
//    }
//  }
//  c1D->Print(outputstop.c_str());
//
//  outputname = Form("figures/%s/%s/pTstudy/Helicity1DBetaRatioDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c1D->Print(outputstart.c_str());
//
//  for(int ieta = 0; ieta < 10; ieta+=10)
//  {
//    for(int ikin = 0; ikin < nkin; ikin++)
//    {
//      for(int irho = 0; irho < 5; irho++)
//      {
//        for(int i = 0; i < nvar; i++)
//        //for(int i = 0; i < nvar; i++)
//        {
//          for(int iopt = 0; iopt < nopts; iopt++)
//          {
//            for(int ibin = 0; ibin < 20; ibin++)
//            {
//              c1D->cd(1+ibin);
//              h_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f, cos%d",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1],ibin));
//              h_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->GetXaxis()->SetTitle("#beta_{h}");
//              h_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->GetYaxis()->SetTitle("Acceptance");
//              h_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("Colz");
//              g_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->SetLineColor(kRed);
//              g_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("C");
//            }
//            c1D->Print(outputname.c_str());
//            c1D->Update();
//          }
//        }
//      }
//    }
//  }
//  c1D->Print(outputstop.c_str());
  //String outputfile = "v2weighting_19GeV.root";
  //TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  //File_OutPut->cd();
  //For(int iopt = 0; iopt < nopts; iopt++)
  //{
  //  for(int ikin = 0; ikin < nkin; ikin++)
  //  {
  //    h_m2D_Weight[iopt][ikin]->Write(); 
  //  } 
  //}
  //File_OutPut->Close();  
 
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

