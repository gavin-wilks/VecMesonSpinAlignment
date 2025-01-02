#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
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

void compareHelicityVsGlobal_PtCut_PtRes_PtDep(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 23, std::string simmode = "Mc", int inputpar = 0)//defaultF = 0 is BESII, defaultF = 1 is BESI
{
  cout << "Did we even start the macro?" << endl;
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

  const int nkin = 1;
  const float ptkin[nkin] = {10.0};
  const float ykin[nkin]  = {10.0};

  const int npt = 20;
  float ptbins[21] = {0.0};
  for(int ipt = 0; ipt < npt+1; ipt++)
  {
    ptbins[ipt] = 0.1*ipt;
  }
  //const int nopts = 2;
  //std::string option = "v2_On_Off";
  //std::string date[nopts] = {"20240809g","20240809g"}; // nov2, prelimv2
  //std::string opts[nopts] = {"nov2","prelimv2"};
  //int color[nopts] = {kBlue, kOrange+7}; 
  //std::string label[nopts] = {"v_{2} OFF","v_{2} ON (Prelim)"};

  const int nopts = 1;

  //std::string option = "y_Variance_RC_EtaCuts_GlobalInput_pt2_individualy";
  std::string option = "GlobalInput_ptcuts_ptres_ptdep_corrected";
  //std::string date[nopts] = {"20240809_2","20240815he","2040815he","20240815he","20240815he","20240815he","20240815he","20240815he","20240815he","20240815he"}; // nov2, prelimv2
  //std::string opts[nopts] = {"nov2","v20.0333","v20.0667","v20.1000","v20.1333","v20.1667","v20.2000","v20.2333","v20.2667","v20.3000"};
  //float vals[nopts] = {0.0,1./30.,2./30,3./30.,4./30.,5./30,6./30.,7./30.,8./30,9./30.,10./30.};
  float vals[nopts] = {10.0};
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

  TFile *File_InPutRC[2][nopts][nkin]; 
  TH3D *h_m3D_RC[2][nopts][5][nvar][nkin][11][9];  
  TH3D *h_m3D_RC_Global[2][nopts][5][nvar][nkin][11][9]; 
  TH2D *h_m2D_RC[2][nopts][5][nvar][nkin][11][9][20];  
  TH2D *h_m2D_RC_Global[2][nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Global_Cos[2][nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Global_Cos_All[2][nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Global_Ratio_Cos_All[2][nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Global_Corrected_Cos_All[2][nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Global_Beta[2][nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Global_Ratio_Cos[2][nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Global_Ratio_Beta[2][nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Cos[2][nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Beta[2][nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Ratio_Cos[2][nopts][5][nvar][nkin][11][9][20]; 
  TH1D *h_m1D_RC_Ratio_Beta[2][nopts][5][nvar][nkin][11][9][20]; 

  TH2D *h_m2D_RC_Global_Ratio[2][nopts][5][nvar][nkin][11][9][npt];
  TH2D *h_m2D_RC_Ratio[2][nopts][5][nvar][nkin][11][9][npt];
  TH2D *h_m2D_RC_Global_Corrected[2][nopts][5][nvar][nkin][11][9][npt];
  TH2D *h_m2D_RC_Global_Corrected_BackToMC[2][nopts][5][nvar][nkin][11][9][npt];
  TH2D *h_m2D_RC_Global_Corrected_BackToMC_E[2][nopts][5][nvar][nkin][11][9][npt];
  TH2D *h_m2D_RC_Corrected[2][nopts][5][nvar][nkin][11][9][npt];
 
  TH2D *h_m2D_Weight[2][nopts][nkin];

  TGraphAsymmErrors *g_m1D_RC_Global_Ratio_Cos[2][nopts][5][nvar][nkin][11][9][20]; 
  TGraphAsymmErrors *g_m1D_RC_Global_Ratio_Beta[2][nopts][5][nvar][nkin][11][9][20]; 
  TGraphAsymmErrors *g_m1D_RC_Ratio_Cos[2][nopts][5][nvar][nkin][11][9][20]; 
  TGraphAsymmErrors *g_m1D_RC_Ratio_Beta[2][nopts][5][nvar][nkin][11][9][20]; 

  

  for(int ires = 0; ires < 2; ires++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int iopt = 0; iopt < nopts; iopt++)
      {
        string inputfileRC = Form("effaccfiles/%s/%s/PaperAllYRC_EtaFixed_GlobalInput_nov2_ptcuts_ptres/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0_pt%1.1f_y%1.1f_noptres_2.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),ptkin[ikin],vals[iopt]);
        if(ires == 1) inputfileRC = Form("effaccfiles/%s/%s/PaperAllYRC_EtaFixed_GlobalInput_nov2_ptcuts_ptres/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0_pt%1.1f_y%1.1f_ptres_2.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),ptkin[ikin],vals[iopt]);
        File_InPutRC[ires][iopt][ikin] = TFile::Open(inputfileRC.c_str());
        cout << "Opened a file" << endl;
        //for(int irho = 0; irho < 5; irho++)
        //{
        //  for(int i = 2; i < 3; i++)
        //  //for(int i = 0; i < nvar; i++)
        //  {
        //    for(int ieta = 0; ieta < 1; ieta++)
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
  }

  for(int ires = 0; ires < 2; ires++)
  {
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
            //for(int ieta = 0; ieta < 1; ieta++)
            for(int ieta = 0; ieta < 1; ieta+=10)
            {
              for(int ipt = 0; ipt < 9; ipt+=4)
              {
                string KEY_3D_RC = Form("h3_m%sEffCosPhiPrimeH_v2_0_rhoinput_%d_rho_%d_eta_%d_pt%d",simmode.c_str(),irho,i,ieta,ipt);
                string KEY_3D_RC_Global = Form("h3_m%sEffCosPhiPrime_v2_0_rhoinput_%d_rho_%d_eta_%d_pt%d",simmode.c_str(),irho,i,ieta,ipt);
                string KEY_3D = Form("h3_m%sEffCosPhiPrimeH_v2_0_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires);
                string KEY_3D_Global = Form("h3_m%sEffCosPhiPrime_v2_0_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires);

                h_m3D_RC[ires][iopt][irho][i][ikin][ieta][ipt] = (TH3D*) ((TH3D*) File_InPutRC[ires][iopt][ikin]->Get(KEY_3D_RC.c_str()))->Clone(KEY_3D.c_str());
                h_m3D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt] = (TH3D*) ((TH3D*) File_InPutRC[ires][iopt][ikin]->Get(KEY_3D_RC_Global.c_str()))->Clone(KEY_3D_Global.c_str());

                for(int iphipt = 0; iphipt < npt; iphipt++)
                {

                  string KEY_2D_Global = Form("h_m%sEffCosPhiPrime_v2_0_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d_pt_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires,iphipt);
                  h_m3D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetRange(iphipt+1,iphipt+1);
                  h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = (TH2D*) h_m3D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt]->Project3D("zy");
                  h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetName(KEY_2D_Global.c_str());
                  h_m3D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetRange();
                   
                  string KEY_2D = Form("h_m%sEffCosPhiPrimeH_v2_0_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d_pt_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires,iphipt);
                  h_m3D_RC[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetRange(iphipt+1,iphipt+1);
                  h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = (TH2D*) h_m3D_RC[ires][iopt][irho][i][ikin][ieta][ipt]->Project3D("zy");
                  h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetName(KEY_2D.c_str());
                  h_m3D_RC[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetRange();

                  string KEY_1D_RC_Global = Form("h_m%sEffCosPhiPrimeCos_v2_0_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d_pt_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires,iphipt);
                  h_m1D_RC_Global_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = (TH1D*) h_m3D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt]->ProjectionY(KEY_1D_RC_Global.c_str(),iphipt+1,iphipt+1,1,20);


                  string KEY = Form("h_m%s2DGLOBALRATIO_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d_pt_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires,iphipt);
                  h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = (TH2D*) h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Clone(KEY.c_str());
                  h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Divide(h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt],h_m2D_RC_Global[ires][iopt][irho][i][ikin][0][0][iphipt],1,1,"B");
                  KEY = Form("h_m%s2DHELICITYRATIO_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d_pt_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires,iphipt);
                  h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = (TH2D*) h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Clone(KEY.c_str());
                  h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Divide(h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt][iphipt],h_m2D_RC[ires][iopt][irho][i][ikin][0][0][iphipt],1,1,"B");

                  KEY = Form("h_m%s1DGLOBALCOSRATIO_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d_pt_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires,iphipt);
                  h_m1D_RC_Global_Ratio_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = (TH1D*) h_m1D_RC_Global_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Clone(KEY.c_str());
                  h_m1D_RC_Global_Ratio_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Divide(h_m1D_RC_Global_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt][iphipt],h_m1D_RC_Global_Cos_All[ires][iopt][irho][i][ikin][0][0][iphipt],1,1,"B");

                  //KEY = Form("h_m%s2DGLOBALCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires);
                  //h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                  //h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]);
                  //KEY = Form("h_m%s2DGLOBALCORRECTEDBTMC_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires);
                  //h_m2D_RC_Global_Corrected_BackToMC[ires][iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                  //h_m2D_RC_Global_Corrected_BackToMC[ires][iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Global[ires][iopt][irho][i][ikin][0][ipt]);
                  //KEY = Form("h_m%s2DHELICITYCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires);
                  //h_m2D_RC_Corrected[ires][iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                  //h_m2D_RC_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]);
                }


                //KEY = Form("h_m%s2DGLOBALRATIO_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires);
                //h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                //h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt],h_m2D_RC_Global[ires][iopt][irho][i][ikin][0][0],1,1,"B");
                //KEY = Form("h_m%s2DHELICITYRATIO_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires);
                //h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                //h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt],h_m2D_RC[ires][iopt][irho][i][ikin][0][0],1,1,"B");


                //KEY = Form("h_m%s2DGLOBALCOSCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires);
                //h_m1D_RC_Global_Corrected_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt] = (TH1D*) h_m1D_RC_Global_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                //h_m1D_RC_Global_Corrected_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt]->Divide(h_m1D_RC_Global_Ratio_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt]);

                //KEY = Form("h_m%s2DGLOBALCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires);
                //h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                //h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]);
                //KEY = Form("h_m%s2DGLOBALCORRECTEDBTMC_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires);
                //h_m2D_RC_Global_Corrected_BackToMC[ires][iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                //h_m2D_RC_Global_Corrected_BackToMC[ires][iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Global[ires][iopt][irho][i][ikin][0][ipt]);
                //KEY = Form("h_m%s2DHELICITYCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires);
                //h_m2D_RC_Corrected[ires][iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                //h_m2D_RC_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]);

              }
            }
          }
        }
      }
    }
  }


  for(int ires = 0; ires < 2; ires++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int iopt = 0; iopt < nopts; iopt++)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          for(int i = 0; i < nvar; i++)
          //for(int i = 0; i < nvar; i++)
          {
            for(int ieta = 0; ieta < 1; ieta+=10)
            {
              for(int ipt = 0; ipt < 9; ipt+=4) 
              {
                for(int iphipt = 0; iphipt < npt; iphipt++)
                {
                  string KEY = Form("h_m%s1DGLOBALCOSCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d_pt_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires,iphipt);
                  h_m1D_RC_Global_Corrected_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = (TH1D*) h_m1D_RC_Global_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Clone(KEY.c_str());
                  h_m1D_RC_Global_Corrected_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Divide(h_m1D_RC_Global_Ratio_Cos_All[0][iopt][irho][2][ikin][ieta][ipt][iphipt]);

                  KEY = Form("h_m%s2DGLOBALCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d_pt_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires,iphipt);
                  h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = (TH2D*) h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Clone(KEY.c_str());
                  h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Divide(h_m2D_RC_Global_Ratio[0][iopt][irho][2][ikin][ieta][ipt][iphipt]);
                }
                //KEY = Form("h_m%s2DGLOBALCORRECTED_v2_6_rhoinput_%d_rho_%d_eta_%d_pt_%d_y_%d_kin_%d_res_%d",simmode.c_str(),irho,i,ieta,ipt,iopt,ikin,ires);
                //h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt]->Clone(KEY.c_str());
                //h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->Divide(h_m2D_RC_Global_Ratio[0][iopt][irho][2][ikin][ieta][ipt]);

                //for(int ibin = 0; ibin < 20; ibin++)
                //{
                //  h_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin] = (TH1D*) h_m1D_RC_Global_Cos[ires][iopt][irho][i][ikin][ieta+1][ipt][ibin];
                //  h_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin] = (TH1D*) h_m1D_RC_Global_Beta[ires][iopt][irho][i][ikin][ieta+1][ipt][ibin];

                //  h_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->Divide(h_m1D_RC_Global_Cos[ires][iopt][irho][i][ikin][ieta+1][ipt][ibin],h_m1D_RC_Global_Cos[ires][iopt][irho][i][ikin][0][0][ibin],1,1,"B");
                //  h_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->Divide(h_m1D_RC_Global_Beta[ires][iopt][irho][i][ikin][ieta+1][ipt][ibin],h_m1D_RC_Global_Beta[ires][iopt][irho][i][ikin][0][0][ibin],1,1,"B");

                //  h_m1D_RC_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin] = (TH1D*) h_m1D_RC_Cos[ires][iopt][irho][i][ikin][ieta+1][ipt][ibin];
                //  h_m1D_RC_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin] = (TH1D*) h_m1D_RC_Beta[ires][iopt][irho][i][ikin][ieta+1][ipt][ibin];

                //  h_m1D_RC_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->Divide(h_m1D_RC_Cos[ires][iopt][irho][i][ikin][ieta+1][ipt][ibin],h_m1D_RC_Cos[ires][iopt][irho][i][ikin][0][0][ibin],1,1,"B");
                //  h_m1D_RC_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->Divide(h_m1D_RC_Beta[ires][iopt][irho][i][ikin][ieta+1][ipt][ibin],h_m1D_RC_Beta[ires][iopt][irho][i][ikin][0][0][ibin],1,1,"B");
                //}
              }
            }
          }
        }
      }
    }
  }





  TGraphAsymmErrors *g2D[5][nvar][nkin][5][2][nopts][11][9]; 
  TGraphAsymmErrors *g2D_Global[5][nvar][nkin][5][2][nopts][11][9]; 
  TGraphAsymmErrors *g2D_Global_Corrected[5][nvar][nkin][5][2][nopts][11][9]; 
  TGraphAsymmErrors *g2D_Corrected[5][nvar][nkin][5][2][nopts][11][9]; 
  TF2 *fits2D[2][nopts][5][nvar][nkin][11][9][20];
  TF2 *fits2D_Global[2][nopts][5][nvar][nkin][11][9][20];
  TF2 *fits2D_Global_Corrected[2][nopts][5][nvar][nkin][11][9][20];
  TF2 *fits2D_Corrected[2][nopts][5][nvar][nkin][11][9][20];


  TGraphAsymmErrors *g1D_Global[5][nvar][nkin][2][nopts][11][9]; 
  TGraphAsymmErrors *g1D_GlobalDiff[5][nvar][nkin][nopts][11][9]; 
  TGraphAsymmErrors *g1D_Global_Corrected[5][nvar][nkin][2][nopts][11][9]; 
  TF1 *fits1D_Global[2][nopts][5][nvar][nkin][11][9][npt];
  TF1 *fits1D_Global_Corrected[2][nopts][5][nvar][nkin][11][9][npt];

  for(int ires = 0; ires < 2; ires++)
  {
    for(int iopt = 0; iopt < nopts; iopt++)
    {
      for(int i = 0; i < nvar; i++)
      {
        for(int ikin = 0; ikin < nkin; ikin++)
        {
          for(int irho = 0; irho < 5; irho++)
          {
            for(int ieta = 0; ieta < 1; ieta+=10)
            {
              for(int ipt = 0; ipt < 9; ipt+=4) 
              {
                g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt] = new TGraphAsymmErrors();
                g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt] = new TGraphAsymmErrors();
                for(int iphipt = 0; iphipt < npt; iphipt++)
                {
                  fits1D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = new TF1(Form("fit1D_GlobalCos_%d_%d_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta,ipt,iphipt),SpinDensity,-1.0,1.0,2);
                  fits1D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(0,1./3.);
                  fits1D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParLimits(0,0.0,1.0);
                  fits1D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(1,h_m1D_RC_Global_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetMaximum());
                  //h_m1D_RC_Global_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Fit(fits1D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt],"NMRI");
                  if(irho == 0) g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires,fits1D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
                  if(irho != 0) g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires,fits1D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(0)-1./3.);
                  g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->SetPointError(iphipt,0.0,0.0,fits1D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(0),fits1D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(0));

                  fits1D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = new TF1(Form("fit1D_Global_CorrectedCos_%d_%d_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta,ipt,iphipt),SpinDensity,-1.0,1.0,2);
                  fits1D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(0,1./3.);
                  fits1D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParLimits(0,0.0,1.0);
                  fits1D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(1,h_m1D_RC_Global_Corrected_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetMaximum());
                  h_m1D_RC_Global_Corrected_Cos_All[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Fit(fits1D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt],"NMRI");
                  if(irho == 0) g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires,fits1D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
                  if(irho != 0) g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires,fits1D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(0)-1./3.);
                  g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->SetPointError(iphipt,0.0,0.0,fits1D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(0),fits1D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(0));
                }
              }
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
          for(int ieta = 0; ieta < 1; ieta+=10)
          {
            for(int ipt = 0; ipt < 9; ipt+=4) 
            {
              g1D_GlobalDiff[irho][i][ikin][iopt][ieta][ipt] = new TGraphAsymmErrors();
              for(int iphipt = 0; iphipt < npt; iphipt++)
              {
                double val1, pt1; 
                double val2, pt2; 
                g1D_Global[irho][i][ikin][0][iopt][ieta][ipt]->GetPoint(iphipt,pt1,val1);
                g1D_Global[irho][i][ikin][1][iopt][ieta][ipt]->GetPoint(iphipt,pt2,val2);
                g1D_GlobalDiff[irho][i][ikin][iopt][ieta][ipt]->SetPoint(iphipt,pt1,val2-val1);
                //g1D_Global[irho][i][ikin][iopt][ieta][ipt]->SetPointError(iphipt,0.0,0.0,fits1D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(0),fits1D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(0));
              }
            }
          }
        }
      }
    }
  }

  for(int ires = 0; ires < 2; ires++)
  {
    for(int iopt = 0; iopt < nopts; iopt++)
    {
      for(int i = 0; i < nvar; i++)
      {
        for(int ikin = 0; ikin < nkin; ikin++)
        {
          for(int irho = 0; irho < 5; irho++)
          {
            for(int ieta = 0; ieta < 1; ieta+=10)
            {
              for(int ipt = 0; ipt < 9; ipt+=4)
              {
                for(int par = 0; par < 5; par++) 
                {
                  g2D[irho][i][ikin][par][ires][iopt][ieta][ipt] = new TGraphAsymmErrors();
                  g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt] = new TGraphAsymmErrors();
                }
                for(int iphipt = 0; iphipt < npt; iphipt++)
                {
                  fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d_%d_%d_%d",ires,iopt,irho,i,ikin,ieta,ipt,iphipt),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),6);
                  fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(0,1./3.);
                  fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParLimits(0,0.0,1.0);
                  fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(1,0.0);
                  fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(2,0.0);
                  fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(3,0.0);
                  fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(4,0.0);
                  fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(5,h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetMaximum());
                  //h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Fit(fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt],"NMRI");
                  for(int par = 0; par < 5; par++)
                  {
                    if(par == 0 && par == irho) g2D[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
                    if(par != 0 && par == irho) g2D[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(par)-(float(i)-2)*0.1);
                    if(par == 0 && par != irho) g2D[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(0)-1./3.);
                    if(par != 0 && par != irho) g2D[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(par));
                    g2D[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPointError(iphipt,0.0,0.0,fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(par),fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(par));
                  }

                  fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = new TF2(Form("fit2D_Global_%d_%d_%d_%d_%d_%d_%d_%d",ires,iopt,irho,i,ikin,ieta,ipt,iphipt),SpinDensity2Dcos,-1.0,1.0,0,2.0*TMath::Pi(),6);
                  fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(0,1./3.);
                  fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParLimits(0,0.0,1.0);
                  fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(1,0.0);
                  fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(2,0.0);
                  fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(3,0.0);
                  fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(4,0.0);
                  fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(5,h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetMaximum());
                  h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Fit(fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt],"NMRI");
                  for(int par = 0; par < 5; par++)
                  {
                    if(par == 0 && par == irho) g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
                    if(par != 0 && par == irho) g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(par)-(float(i)-2)*0.1);
                    if(par == 0 && par != irho) g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(0)-1./3.);
                    if(par != 0 && par != irho) g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(par));
                    g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPointError(iphipt,0.0,0.0,fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(par),fits2D_Global[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(par));
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  for(int ires = 0; ires < 2; ires++)
  {
    for(int iopt = 0; iopt < nopts; iopt++)
    {
      for(int i = 0; i < nvar; i++)
      {
        for(int ikin = 0; ikin < nkin; ikin++)
        {
          for(int irho = 0; irho < 5; irho++)
          {
            for(int ieta = 0; ieta < 1; ieta+=10)
            {
              for(int ipt = 0; ipt < 9; ipt+=4)
              {
                for(int par = 0; par < 5; par++) 
                {
                  g2D_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt] = new TGraphAsymmErrors();
                  g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt] = new TGraphAsymmErrors();
                }
                for(int iphipt = 0; iphipt < npt; iphipt++)
                {
                  fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d_%d_%d_%d",ires,iopt,irho,i,ikin,ieta,ipt,iphipt),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),6);
                  fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(0,1./3.);
                  fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParLimits(0,0.0,1.0);
                  fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(1,0.0);
                  fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(2,0.0);
                  fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(3,0.0);
                  fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(4,0.0);
                  fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(5,h_m2D_RC_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetMaximum());
                  //h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Fit(fits2D[ires][iopt][irho][i][ikin][ieta][ipt][iphipt],"NMRI");
                  for(int par = 0; par < 5; par++)
                  {
                    if(par == 0 && par == irho) g2D_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
                    if(par != 0 && par == irho) g2D_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(par)-(float(i)-2)*0.1);
                    if(par == 0 && par != irho) g2D_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(0)-1./3.);
                    if(par != 0 && par != irho) g2D_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(par));
                    g2D_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPointError(iphipt,0.0,0.0,fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(par),fits2D_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(par));
                  }
//                  cout << "Fit helicity " << endl;

                  fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt] = new TF2(Form("fit2D_Global_Corrected_%d_%d_%d_%d_%d_%d_%d_%d",ires,iopt,irho,i,ikin,ieta,ipt,iphipt),SpinDensity2Dcos,-1.0,1.0,0,2.0*TMath::Pi(),6);
                  fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(0,1./3.);
                  fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParLimits(0,0.0,1.0);
                  fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(1,0.0);
                  fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(2,0.0);
                  fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(3,0.0);
                  fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(4,0.0);
                  fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->SetParameter(5,h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetMaximum());
                  h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->Fit(fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt],"NMRI");
                  for(int par = 0; par < 5; par++)
                  {
                    if(par == 0 && par == irho) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
                    if(par != 0 && par == irho) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(par)-(float(i)-2)*0.1);
                    if(par == 0 && par != irho) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(0)-1./3.);
                    if(par != 0 && par != irho) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPoint(iphipt,(ptbins[iphipt]+ptbins[iphipt+1])/2.+0.025*ires, fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParameter(par));
                    g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetPointError(iphipt,0.0,0.0,fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(par),fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt][iphipt]->GetParError(par));
                  }
                }
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

  TLegend *leg;// = new TLegend(0.4,0.75,0.7,0.9);
  leg = new TLegend(0.25,0.775,0.85,0.865);
  leg->SetNColumns(2);
  leg->SetBorderSize(0);

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int ieta = 0; ieta < 1; ieta+=10)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          for(int ipt = 0; ipt < 9; ipt+=4) 
          {
            for(int par = 0; par < 5; par++)
            {
              for(int i = 0; i < nvar; i++)
              {
                double max = -999999;
                double min = 999999;

                cfit2D->cd(1+i+5*par);
                for(int ires = 0; ires < 2; ires++)
                {
                  double tmax = g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetHistogram()->GetMaximum();   
                  double tmin = g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetHistogram()->GetMinimum(); 

                  if(tmax > max) max = tmax;
                  if(tmin < min) min = tmin;
                }
                for(int ires = 0; ires < 2; ires++)
                {
                  g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetMarkerStyle(20);
                  g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetMarkerSize(1.0);
                  g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetMarkerColor(color[ires]);
                  g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetLineColor(color[ires]);
                  g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetTitle(Form("Input %s^{g}=%1.1f, v_{2}=0.0, p_{T,K^{+/-}}>%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),(float(i)-2.)*0.1,float(ipt)*0.025,etacuts[ieta]));
                  g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetXaxis()->SetTitle(Form("(GeV/c) p_{T,#phi}"));
                  g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetYaxis()->SetTitle(Form("%s^{g} from 2D Fit - Truth",param[par].c_str()));
                  g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetXaxis()->SetRangeUser(-0.1,2.1);
                  if(min < 0.0 && max < 0.0) g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(1.5*min,0.5*max);
                  if(min < 0.0 && max > 0.0) g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(1.5*min,1.5*max);
                  if(min > 0.0 && max > 0.0) g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(0.5*min,1.5*max);

                  if(ires == 0) g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->Draw("APE"); 
                  else       g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt]->Draw("PE same"); 

                  if(ikin == 0 && iopt == 0 && ieta == 0 && par == 0 && i == 0 && irho == 0 && par == 0 && ires == 0 && ipt == 0) leg->AddEntry(g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt],"MC p_{T}","p");
                  if(ikin == 0 && iopt == 0 && ieta == 0 && par == 0 && i == 0 && irho == 0 && par == 0 && ires == 1 && ipt == 0) leg->AddEntry(g2D_Global[irho][i][ikin][par][ires][iopt][ieta][ipt],"p_{T} res. = 5%","p");
                }
                leg->Draw("same");
              }
            }
            cfit2D->Print(outputname.c_str());
            cfit2D->Update();
          }
        }
      }   
    }
  }
   
  cfit2D->Print(outputstop.c_str());

  outputname = Form("figures/%s/%s/pTstudy/GlobalCorrrectedvsGlobal_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  cfit2D->Print(outputstart.c_str());


  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int ieta = 0; ieta < 1; ieta+=10)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          for(int ipt = 0; ipt < 9; ipt+=4) 
          {
            for(int par = 0; par < 5; par++)
            {
              for(int i = 0; i < nvar; i++)
              {
                double max = -999999;
                double min = 999999;

                cfit2D->cd(1+i+5*par);
                for(int ires = 0; ires < 2; ires++)
                {
                  double tmax = g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetHistogram()->GetMaximum();   
                  double tmin = g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetHistogram()->GetMinimum(); 

                  if(tmax > max) max = tmax;
                  if(tmin < min) min = tmin;
                }
                for(int ires = 0; ires < 2; ires++)
                {
                  g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetMarkerStyle(20);
                  g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetMarkerSize(1.0);
                  g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetMarkerColor(color[ires]);
                  g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetLineColor(color[ires]);
                  g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->SetTitle(Form("Input %s^{g}=%1.1f, v_{2}=0.0, p_{T,K^{+/-}}>%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),(float(i)-2.)*0.1,float(ipt)*0.025,etacuts[ieta]));
                  g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetXaxis()->SetTitle(Form("(GeV/c) p_{T,#phi}"));
                  g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetYaxis()->SetTitle(Form("%s^{g} from 2D Fit - Truth",param[par].c_str()));
                  g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetXaxis()->SetRangeUser(-0.1,2.1);
                  if(min < 0.0 && max < 0.0) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(1.5*min,0.5*max);
                  if(min < 0.0 && max > 0.0) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(1.5*min,1.5*max);
                  if(min > 0.0 && max > 0.0) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(0.5*min,1.5*max);

                  if(ires == 0) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->Draw("APE"); 
                  else       g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta][ipt]->Draw("PE same"); 

                }
                leg->Draw("same");
              }
            }
            cfit2D->Print(outputname.c_str());
            cfit2D->Update();
          }
        }
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
      for(int ieta = 0; ieta < 1; ieta+=10)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          //cfit2D->cd(1+irho+5*par);
          for(int ipt = 0; ipt < 9; ipt+=4) 
          {
            for(int i = 0; i < nvar; i++)
            {
              double max = -999999;
              double min = 999999;

              cfit1D->cd(1+i);
              for(int ires = 0; ires < 2; ires++)
              {
                double tmax = g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetHistogram()->GetMaximum();   
                double tmin = g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetHistogram()->GetMinimum(); 

                if(tmax > max) max = tmax;
                if(tmin < min) min = tmin;
              }
              for(int ires = 0; ires < 2; ires++)
              {
                g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->SetMarkerStyle(20);
                g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->SetMarkerSize(1.0);
                g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->SetMarkerColor(color[ires]);
                g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->SetLineColor(color[ires]);
                g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->SetTitle(Form("Input %s^{g}=%1.1f, v_{2}=0.0, p_{T,K^{+/-}}>%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),(float(i)-2.)*0.1,float(ipt)*0.025,etacuts[ieta]));
                g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetXaxis()->SetTitle(Form("p_{T,#phi} (GeV/c)"));
                g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetYaxis()->SetTitle(Form("%s^{g} from 1D Fit - Truth",param[0].c_str()));
                g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetXaxis()->SetRangeUser(-0.1,2.1);
                if(min < 0.0 && max < 0.0) g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(1.5*min,0.5*max);
                if(min < 0.0 && max > 0.0) g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(1.5*min,1.5*max);
                if(min > 0.0 && max > 0.0) g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(0.5*min,1.5*max);

                if(ires == 0) g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->Draw("APE"); 
                else       g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->Draw("PE same"); 

                //if(ikin == 0 && iopt == 0 && ieta == 0 && i == 0 && irho == 0 && ipt == 0 &&  ires == 0) leg->AddEntry(g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt],"MC p_{T}","p");
                //if(ikin == 0 && iopt == 0 && ieta == 0 && i == 0 && irho == 0 && ipt == 0 &&  ires == 1) leg->AddEntry(g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt],"p_{T} res. = 5%","p");

              }
              leg->Draw("same");
            }
            cfit1D->Print(outputname.c_str());
            cfit1D->Update();
          }
        }
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
      for(int ieta = 0; ieta < 1; ieta+=10)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          //cfit2D->cd(1+irho+5*par);
          for(int ipt = 0; ipt < 9; ipt+=4) 
          {
            for(int i = 0; i < nvar; i++)
            {
              double max = -999999;
              double min = 999999;

              cfit1D->cd(1+i);
              for(int ires = 0; ires < 2; ires++)
              {
                double tmax = g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->GetHistogram()->GetMaximum();   
                double tmin = g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->GetHistogram()->GetMinimum(); 

                if(tmax > max) max = tmax;
                if(tmin < min) min = tmin;
              }
              for(int ires = 0; ires < 2; ires++)
              {
                g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->SetMarkerStyle(20);
                g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->SetMarkerSize(1.0);
                g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->SetMarkerColor(color[ires]);
                g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->SetLineColor(color[ires]);
                g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->SetTitle(Form("Input %s^{g}=%1.1f, v_{2}=0.0, p_{T,K^{+/-}}>%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),(float(i)-2.)*0.1,float(ipt)*0.025,etacuts[ieta]));
                g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->GetXaxis()->SetTitle(Form("p_{T,#phi} (GeV/c)"));
                g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->GetYaxis()->SetTitle(Form("%s^{g} from 1D Fit - Truth",param[0].c_str()));
                g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->GetXaxis()->SetRangeUser(-0.1,2.1);
                if(min < 0.0 && max < 0.0) g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(1.5*min,0.5*max);
                if(min < 0.0 && max > 0.0) g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(1.5*min,1.5*max);
                if(min > 0.0 && max > 0.0) g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(0.5*min,1.5*max);

                if(ires == 0) g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->Draw("APE"); 
                else       g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta][ipt]->Draw("PE same"); 

              }
              leg->Draw("same");
            }
            cfit1D->Print(outputname.c_str());
            cfit1D->Update();
          }
        }
      }
    }
  }
   
  cfit1D->Print(outputstop.c_str());

  //outputname = Form("figures/%s/%s/pTstudy/Global1DResMNoRes_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  //outputstart = Form("%s[",outputname.c_str()); 
  //outputstop = Form("%s]",outputname.c_str()); 

  //cfit1D->Print(outputstart.c_str());

  //for(int iopt = 0; iopt < nopts; iopt++)
  //{
  //  for(int ikin = 0; ikin < nkin; ikin++)
  //  {
  //    for(int ieta = 0; ieta < 1; ieta+=10)
  //    {
  //      for(int irho = 0; irho < 5; irho++)
  //      {
  //        //cfit2D->cd(1+irho+5*par);
  //        for(int ipt = 0; ipt < 9; ipt+=4) 
  //        {
  //          for(int i = 0; i < nvar; i++)
  //          {
  //            double max = -999999;
  //            double min = 999999;

  //            cfit1D->cd(1+i);
  //            //for(int ires = 0; ires < 2; ires++)
  //            //{
  //            //  double tmax = g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetHistogram()->GetMaximum();   
  //            //  double tmin = g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetHistogram()->GetMinimum(); 

  //            //  if(tmax > max) max = tmax;
  //            //  if(tmin < min) min = tmin;
  //            //}
  //            //for(int ires = 0; ires < 2; ires++)
  //            //{
  //              g1D_GlobalDiff[irho][i][ikin][iopt][ieta][ipt]->SetMarkerStyle(20);
  //              g1D_GlobalDiff[irho][i][ikin][iopt][ieta][ipt]->SetMarkerSize(1.0);
  //              //g1D_GlobalDiff[irho][i][ikin][iopt][ieta][ipt]->SetMarkerColor(color[ires]);
  //              //g1D_GlobalDiff[irho][i][ikin][iopt][ieta][ipt]->SetLineColor(color[ires]);
  //              g1D_GlobalDiff[irho][i][ikin][iopt][ieta][ipt]->SetTitle(Form("Input %s^{g}=%1.1f, v_{2}=0.0, p_{T,K^{+/-}}>%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),(float(i)-2.)*0.1,float(ipt)*0.025,etacuts[ieta]));
  //              g1D_GlobalDiff[irho][i][ikin][iopt][ieta][ipt]->GetXaxis()->SetTitle(Form("p_{T,#phi} (GeV/c)"));
  //              g1D_GlobalDiff[irho][i][ikin][iopt][ieta][ipt]->GetYaxis()->SetTitle(Form("%s^{g} (5%% p_{T} Res. - MC)",param[0].c_str()));
  //              g1D_GlobalDiff[irho][i][ikin][iopt][ieta][ipt]->GetXaxis()->SetRangeUser(-0.1,2.1);
  //              //if(min < 0.0 && max < 0.0) g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(1.5*min,0.5*max);
  //              //if(min < 0.0 && max > 0.0) g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(1.5*min,1.5*max);
  //              //if(min > 0.0 && max > 0.0) g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->GetYaxis()->SetRangeUser(0.5*min,1.5*max);

  //              g1D_GlobalDiff[irho][i][ikin][iopt][ieta][ipt]->Draw("APE"); 
  //              //else       g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt]->Draw("PE same"); 

  //              //if(ikin == 0 && iopt == 0 && ieta == 0 && i == 0 && irho == 0 && ipt == 0 &&  ires == 0) leg->AddEntry(g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt],"MC p_{T}","p");
  //              //if(ikin == 0 && iopt == 0 && ieta == 0 && i == 0 && irho == 0 && ipt == 0 &&  ires == 1) leg->AddEntry(g1D_Global[irho][i][ikin][ires][iopt][ieta][ipt],"p_{T} res. = 5%","p");

  //            //}
  //            //leg->Draw("same");
  //          }
  //          cfit1D->Print(outputname.c_str());
  //          cfit1D->Update();
  //        }
  //      }
  //    }
  //  }
  //}
  // 
  //cfit1D->Print(outputstop.c_str());

  //outputname = Form("figures/%s/%s/pTstudy/Global1DCorrectedvsGlobal_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  //outputstart = Form("%s[",outputname.c_str()); 
  //outputstop = Form("%s]",outputname.c_str()); 

  //cfit1D->Print(outputstart.c_str());

  //for(int iopt = 0; iopt < nopts; iopt++)
  //{
  //  for(int ikin = 0; ikin < nkin; ikin++)
  //  {
  //    for(int ieta = 0; ieta < 1; ieta+=10)
  //    {
  //      for(int irho = 0; irho < 5; irho++)
  //      {
  //        //cfit2D->cd(1+irho+5*par);
  //        for(int i = 0; i < nvar; i++)
  //        {
  //          double max = -999999;
  //          double min = 999999;

  //          cfit1D->cd(1+i);
  //          for(int ires = 0; ires < 2; ires++)
  //          {
  //            double tmax = g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->GetHistogram()->GetMaximum();   
  //            double tmin = g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->GetHistogram()->GetMinimum(); 

  //            if(tmax > max) max = tmax;
  //            if(tmin < min) min = tmin;
  //          }
  //          for(int ires = 0; ires < 2; ires++)
  //          {
  //            g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->SetMarkerStyle(20);
  //            g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->SetMarkerSize(1.0);
  //            g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->SetMarkerColor(color[ires]);
  //            g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->SetLineColor(color[ires]);
  //            g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->SetTitle(Form("Input %s^{g}=%1.1f, v_{2}=0.0, p_{T}=%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),(float(i)-2.)*0.1,ptkin[ikin],etacuts[ieta]));
  //            g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->GetXaxis()->SetTitle(Form("(GeV/c) p_{T,K^{+/-}}<"));
  //            g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->GetYaxis()->SetTitle(Form("%s^{g} from 2D Fit - Truth",param[0].c_str()));
  //            g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->GetXaxis()->SetRangeUser(-0.025,0.225);
  //            if(min < 0.0 && max < 0.0) g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,0.4*max);
  //            if(min < 0.0 && max > 0.0) g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,1.6*max);
  //            if(min > 0.0 && max > 0.0) g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->GetYaxis()->SetRangeUser(0.4*min,1.6*max);

  //            if(ires == 0) g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->Draw("APE"); 
  //            else       g1D_Global_Corrected[irho][i][ikin][ires][iopt][ieta]->Draw("PE same"); 

  //          }
  //          leg->Draw("same");
  //        }
  //        cfit1D->Print(outputname.c_str());
  //        cfit1D->Update();
  //      }
  //    }
  //  }
  //}
  // 
  //cfit1D->Print(outputstop.c_str());





  //outputname = Form("figures/%s/%s/pTstudy/HelicityvsGlobal_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  //outputstart = Form("%s[",outputname.c_str()); 
  //outputstop = Form("%s]",outputname.c_str()); 

  //cfit2D->Print(outputstart.c_str());


  ////TLegend *leg;// = new TLegend(0.4,0.75,0.7,0.9);
  //for(int iopt = 0; iopt < nopts; iopt++)
  //{
  //  for(int ikin = 0; ikin < nkin; ikin++)
  //  {
  //    for(int ieta = 0; ieta < 1; ieta++)
  //    {
  //      for(int irho = 0; irho < 5; irho++)
  //      {
  //        //double max = -999999;
  //        //double min = 999999;
  //        //for(int par = 0; par < 5; par++)
  //        //{
  //        //  for(int i = 0; i < nvar; i++)
  //        //  {
  //        //    double tmax = g2D[irho][i][ikin][par][ires][iopt]->GetHistogram()->GetMaximum();   
  //        //    double tmin = g2D[irho][i][ikin][par][ires][iopt]->GetHistogram()->GetMinimum(); 

  //        //    if(tmax > max) max = tmax;
  //        //    if(tmin < min) min = tmin;

  //        //    //cout << "max = " << max << ", min = " << min << endl;
  //        //  }
  //        //}
  //        for(int par = 0; par < 5; par++)
  //        {
  //          cfit2D->cd(1+par+5*irho);
  //          //cfit2D->cd(1+irho+5*par);
  //          //TLegend *leg = new TLegend(0.4,0.6,0.6,0.8);
  //          double max = -999999;
  //          double min = 999999;
  //          for(int i = 0; i < nvar; i++)
  //          {
  //            double tmax = g2D[irho][i][ikin][par][ires][iopt][ieta]->GetHistogram()->GetMaximum();   
  //            double tmin = g2D[irho][i][ikin][par][ires][iopt][ieta]->GetHistogram()->GetMinimum(); 

  //            if(tmax > max) max = tmax;
  //            if(tmin < min) min = tmin;

  //            //cout << "max = " << max << ", min = " << min << endl;
  //          }
  //          for(int i = 0; i < nvar; i++)
  //          {
  //            g2D[irho][i][ikin][par][ires][iopt][ieta]->SetMarkerStyle(20);
  //            g2D[irho][i][ikin][par][ires][iopt][ieta]->SetMarkerSize(1.2);
  //            g2D[irho][i][ikin][par][ires][iopt][ieta]->SetMarkerColor(color[i]);
  //            g2D[irho][i][ikin][par][ires][iopt][ieta]->SetLineColor(color[i]);
  //            g2D[irho][i][ikin][par][ires][iopt][ieta]->SetTitle(Form("Input %s^{g}, v_{2}=0.2, p_{T}=%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),ptkin[ikin],etacuts[ieta]));
  //            g2D[irho][i][ikin][par][ires][iopt][ieta]->GetXaxis()->SetTitle(Form("(GeV/c) p_{T,K^{+/-}}<"));
  //            g2D[irho][i][ikin][par][ires][iopt][ieta]->GetYaxis()->SetTitle(Form("%s^{Helicity} from 2D Fit",param[par].c_str()));
  //            if(min < 0.0 && max < 0.0) g2D[irho][i][ikin][par][ires][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,0.4*max);
  //            if(min < 0.0 && max > 0.0) g2D[irho][i][ikin][par][ires][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,1.6*max);
  //            if(min > 0.0 && max > 0.0) g2D[irho][i][ikin][par][ires][iopt][ieta]->GetYaxis()->SetRangeUser(0.4*min,1.6*max);
  //            g2D[irho][i][ikin][par][ires][iopt][ieta]->GetXaxis()->SetRangeUser(0.025,0.225);

  //            if(i == 0) g2D[irho][i][ikin][par][ires][iopt][ieta]->Draw("APE"); 
  //            else       g2D[irho][i][ikin][par][ires][iopt][ieta]->Draw("PE same"); 
  //              
  //            //if(par == 0) leg->AddEntry(g2D[irho][i][ikin][par],Form("%s^{Helicity} = %1.1f",param[inputpar].c_str(),(float(i)-2)*0.1),"p");
  //          }
  //          leg[irho]->Draw("same");
  //        }
  //      }
  //      cfit2D->Print(outputname.c_str());
  //      cfit2D->Update();
  //    }
  //  }
  //}
  //cfit2D->Print(outputstop.c_str());

  //TCanvas *c2D = new TCanvas("c2D","c2D",10,10,1600,1200);
  //c2D->Divide(4,3);
  //
  //for(int i = 0; i < 12; i++)
  //{
  //  c2D->cd(i+1);
  //  c2D->cd(i+1)->SetLeftMargin(0.15);  
  //  c2D->cd(i+1)->SetRightMargin(0.15);  
  //  c2D->cd(i+1)->SetBottomMargin(0.15);
  //  c2D->cd(i+1)->SetTicks(1,1);
  //  c2D->cd(i+1)->SetGrid(0,0);
  //}
 
  //outputname = Form("figures/%s/%s/pTstudy/Global2DDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  //outputstart = Form("%s[",outputname.c_str()); 
  //outputstop = Form("%s]",outputname.c_str()); 

  //c2D->Print(outputstart.c_str());

  //for(int ipt = 0; ipt < 9; ipt+=4)
  //{
  //  for(int ieta = 0; ieta < 1; ieta++)
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
  //            c2D->cd(1+iopt);
  //            h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{h}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[ires][iopt],etacuts[ieta]));
  //            h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
  //            h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{g}");
  //            h_m2D_RC_Global[ires][iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
  //          }
  //          c2D->Print(outputname.c_str());
  //          c2D->Update();
  //        }
  //      }
  //    }
  //  }
  //}
  //c2D->Print(outputstop.c_str());

  //outputname = Form("figures/%s/%s/pTstudy/Helicity2DDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  //outputstart = Form("%s[",outputname.c_str()); 
  //outputstop = Form("%s]",outputname.c_str()); 

  //c2D->Print(outputstart.c_str());

  //for(int ipt = 0; ipt < 9; ipt+=4)
  //{
  //  for(int ieta = 1; ieta < 1; ieta++)
  //  {
  //    for(int ikin = 0; ikin < nkin; ikin++)
  //    {
  //      for(int irho = 0; irho < 5; irho++)
  //      {
  //        for(int i = 0; i < nvar; i++)
  //        {
  //          for(int iopt = 0; iopt < nopts; iopt++)
  //          {
  //            c2D->cd(1+iopt);
  //            h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{h}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[ires][iopt],etacuts[ieta]));
  //            h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{h}");
  //            h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{h}");
  //            h_m2D_RC[ires][iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
  //          }
  //          c2D->Print(outputname.c_str());
  //          c2D->Update();
  //        }
  //      }
  //    }
  //  }
  //}
  //c2D->Print(outputstop.c_str());


  //TF2* fitacc2D[2][nopts][5][nvar][nkin][10][9];
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

  //          fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
  //          //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple7,-1.0,1.0,0.0,2.0*TMath::Pi(),12);
  //          //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple8,-1.0,1.0,0.0,2.0*TMath::Pi(),16);
  //          fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
  //          
  //          if(ieta < 9)
  //          {
  //            for(int ipar = 1; ipar < 9; ipar++) fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(ipar,fitacc2D[ires][iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
  //            //for(int ipar = 1; ipar < 12; ipar++) fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(ipar,fitacc2D[ires][iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
  //            //for(int ipar = 1; ipar < 16; ipar++) fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(ipar,fitacc2D[ires][iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
  //          }
  //          //h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt],"NMR");
 
  //          double chi2 = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->GetChisquare();
  //          double ndf  = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->GetNDF();

  //          cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
  //          g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->SetPoint(iopt,vals[ires][iopt],chi2/ndf);

  //          for(int iy = 1; iy <= 20; iy++)
  //          { 
  //            g_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][iy-1] = new TGraphAsymmErrors();
  //          }
  //          for(int ix = 1; ix <= 20; ix++)
  //          {
  //            g_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ix-1] = new TGraphAsymmErrors();
  //            for(int iy = 1; iy <= 20; iy++)
  //            { 
  //              double widthx = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinWidth(ix);
  //              double startx = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinLowEdge(ix);
  //              double stopx  = startx+widthx;
  //              double widthy = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinWidth(iy);
  //              double starty = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinLowEdge(iy);
  //              double stopy  = starty+widthy;
  //              double value = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
  //              
  //              double bincenterx = (startx+stopx)/2.0;
  //              double bincentery = (starty+stopy)/2.0;
  //  
  //              if(h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetBinContent(ix,iy) <= 0.0) 
  //              cout << "RATIO IS "<<h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetBinContent(ix,iy)<<"+/-"<<h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetBinError(ix,iy)<<"!?   iopt="<<iopt<<",irho="<<irho<<",i="<<i<<",ikin"<<ikin<<",ieta"<<ieta<<endl;

  //              g_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][iy-1]->SetPoint(ix-1,bincenterx,value); 
  //              g_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ix-1]->SetPoint(iy-1,bincentery,value); 

  //            }
  //          }
  //        }
  //      }
  //    }
  //  }
  //}

  //TF2* fitacc2DH[2][nopts][5][nvar][nkin][10][9];
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

  //          fitacc2DH[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
  //          //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple7,-1.0,1.0,0.0,2.0*TMath::Pi(),12);
  //          fitacc2DH[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
  //          
  //          if(ieta < 9)
  //          {
  //            for(int ipar = 1; ipar < 9; ipar++) fitacc2DH[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(ipar,fitacc2DH[ires][iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
  //            //for(int ipar = 1; ipar < 12; ipar++) fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(ipar,fitacc2D[ires][iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
  //          }
  //          //h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2DH[ires][iopt][irho][i][ikin][ieta][ipt],"NMR");

  //          double chi2 = fitacc2DH[ires][iopt][irho][i][ikin][ieta][ipt]->GetChisquare();
  //          double ndf  = fitacc2DH[ires][iopt][irho][i][ikin][ieta][ipt]->GetNDF();

  //          cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
  //          g_m2DH_chi2ndf[irho][i][ikin][ieta][ipt]->SetPoint(iopt,vals[ires][iopt],chi2/ndf);

  //          for(int iy = 1; iy <= 20; iy++)
  //          { 
  //            g_m1D_RC_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][iy-1] = new TGraphAsymmErrors();
  //          }
  //          for(int ix = 1; ix <= 20; ix++)
  //          {
  //            g_m1D_RC_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ix-1] = new TGraphAsymmErrors();
  //            for(int iy = 1; iy <= 20; iy++)
  //            { 
  //              double widthx = h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinWidth(ix);
  //              double startx = h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinLowEdge(ix);
  //              double stopx  = startx+widthx;
  //              double widthy = h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinWidth(iy);
  //              double starty = h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinLowEdge(iy);
  //              double stopy  = starty+widthy;
  //              double value = fitacc2DH[ires][iopt][irho][i][ikin][ieta][ipt]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                

  //              double bincenterx = (startx+stopx)/2.0;
  //              double bincentery = (starty+stopy)/2.0;

  //              g_m1D_RC_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][iy-1]->SetPoint(ix-1,bincenterx,value); 
  //              g_m1D_RC_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ix-1]->SetPoint(iy-1,bincentery,value); 

  //            }
  //          }
  //        }
  //      }
  //    }
  //  }
  //}

  //for(int ires = 0; ires < 2; ires++)
  //{
  //  for(int iopt = 0; iopt < nopts; iopt++)
  //  {
  //    for(int i = 0; i < nvar; i++)
  //    {
  //      for(int ikin = 0; ikin < nkin; ikin++)
  //      {
  //        for(int ieta = 0; ieta < 1; ieta+=10) 
  //        {
  //          for(int irho = 0; irho < 5; irho++)
  //          {
  //            for(int par = 0; par < 5; par++) 
  //            {
  //              g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta] = new TGraphAsymmErrors();
  //            }   
  //            for(int ipt = 0; ipt < 9; ipt+=4)
  //            {
  //              
  //              fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),6);
  //              //fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2DcosAcc,-1.0,1.0,0.0,2.0*TMath::Pi(),14);
  //              //fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2DcosAccLonger,-1.0,1.0,0.0,2.0*TMath::Pi(),17);
  //              //fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2DcosAccLongerer,-1.0,1.0,0.0,2.0*TMath::Pi(),21);
  //              fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,1./3.);
  //              fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->SetParLimits(0,0.0,1.0);
  //              fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(1,0.0);
  //              fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(2,0.0);
  //              fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(3,0.0);
  //              fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(4,0.0);
  //              fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(5,h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
  //              //for(int ifit = 1; ifit < 9; ifit++)
  //              ////for(int ifit = 1; ifit < 12; ifit++)
  //              ////for(int ifit = 1; ifit < 16; ifit++)
  //              //{
  //              //  double value = fitacc2D[ires][iopt][irho][i][ikin][ieta-1]->GetParameter(ifit);
  //              //  double error = fitacc2D[ires][iopt][irho][i][ikin][ieta-1]->GetParError(ifit);
  //              //  fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->FixParameter(5+ifit,value);
  //              //  //fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->SetParLimits(5+ifit,value-error,value+error);
  //              //}
  //              //fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->FixParameter(7,fitacc2D[ires][iopt][irho][i][ikin][ieta-1]->GetParameter(2));
  //              //fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->FixParameter(8,fitacc2D[ires][iopt][irho][i][ikin][ieta-1]->GetParameter(3));
  //              //fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->FixParameter(9,fitacc2D[ires][iopt][irho][i][ikin][ieta-1]->GetParameter(4));
  //              //fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->FixParameter(10,fitacc2D[ires][iopt][irho][i][ikin][ieta-1]->GetParameter(5));
  //              //fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->FixParameter(11,fitacc2D[ires][iopt][irho][i][ikin][ieta-1]->GetParameter(6));
  //              //fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->FixParameter(12,fitacc2D[ires][iopt][irho][i][ikin][ieta-1]->GetParameter(7));
  //              //fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->FixParameter(13,fitacc2D[ires][iopt][irho][i][ikin][ieta-1]->GetParameter(8));
  //             
  //              h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->Fit(fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt],"NMRI");
  //              for(int par = 0; par < 5; par++)
  //              {
  //                if(par == 0 && par == irho) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025,fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
  //                if(par != 0 && par == irho) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025,fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->GetParameter(par)-(float(i)-2)*0.1);
  //                if(par == 0 && par != irho) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025,fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->GetParameter(0)-1./3.);
  //                if(par != 0 && par != irho) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->SetPoint(ipt,float(ipt)*0.025,fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->GetParameter(par));
  //                g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->SetPointError(ipt,0.0,0.0,fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->GetParError(par),fits2D_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->GetParError(par));
  //              }
  //            }
  //          }
  //        }
  //      }
  //    }
  //  }
  //}

  //outputname = Form("figures/%s/%s/pTstudy/GlobalCorrectedvsGlobal_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  //outputstart = Form("%s[",outputname.c_str()); 
  //outputstop = Form("%s]",outputname.c_str()); 

  //cfit2D->Print(outputstart.c_str());


  //for(int iopt = 0; iopt < nopts; iopt++)
  //{
  //  for(int ikin = 0; ikin < nkin; ikin++)
  //  {
  //    for(int ieta = 0; ieta < 1; ieta+=10)
  //    {
  //      for(int irho = 0; irho < 5; irho++)
  //      {
  //        for(int par = 0; par < 5; par++)
  //        {
  //          //cfit2D->cd(1+irho+5*par);
  //          for(int i = 0; i < nvar; i++)
  //          {
  //            double max = -999999;
  //            double min = 999999;

  //            cfit2D->cd(1+i+5*par);
  //            for(int ires = 0; ires < 2; ires++)
  //            {
  //              double tmax = g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->GetHistogram()->GetMaximum();   
  //              double tmin = g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->GetHistogram()->GetMinimum(); 

  //              if(tmax > max) max = tmax;
  //              if(tmin < min) min = tmin;
  //            }
  //            for(int ires = 0; ires < 2; ires++)
  //            {
  //              g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->SetMarkerStyle(20);
  //              g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->SetMarkerSize(1.0);
  //              g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->SetMarkerColor(color[ires]);
  //              g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->SetLineColor(color[ires]);
  //              g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->SetTitle(Form("Input %s^{g}=%1.1f, v_{2}=0.0, p_{T}=%1.1f GeV/c, |#eta|<%1.1f",param[irho].c_str(),(float(i)-2.)*0.1,ptkin[ikin],etacuts[ieta]));
  //              g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->GetXaxis()->SetTitle(Form("(GeV/c) p_{T,K^{+/-}}<"));
  //              g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->GetYaxis()->SetTitle(Form("%s^{g} from 2D Fit - Truth",param[par].c_str()));
  //              g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->GetXaxis()->SetRangeUser(-0.025,0.225);
  //              if(min < 0.0 && max < 0.0) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,0.4*max);
  //              if(min < 0.0 && max > 0.0) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->GetYaxis()->SetRangeUser(1.6*min,1.6*max);
  //              if(min > 0.0 && max > 0.0) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->GetYaxis()->SetRangeUser(0.4*min,1.6*max);

  //              if(ires == 0) g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->Draw("APE"); 
  //              else       g2D_Global_Corrected[irho][i][ikin][par][ires][iopt][ieta]->Draw("PE same"); 

  //            }
  //            leg->Draw("same");
  //          }
  //        }
  //        cfit2D->Print(outputname.c_str());
  //        cfit2D->Update();
  //      }
  //    }
  //  }
  //}
  // 
  //cfit2D->Print(outputstop.c_str());

//  outputname = Form("figures/%s/%s/pTstudy/Global2DRatioDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c2D->Print(outputstart.c_str());
//
//  for(int ieta = 0; ieta < 10; ieta++)
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
//            //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
//            //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
//
//            h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[ires][iopt],etacuts[ieta+1]));
//            cout << "Set title" << endl;
//            h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
//            h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{g}");
//            h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
//
////            h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt],"NMR");
////            //h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt],"NMR");
////            double chi2 = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->GetChisquare();
////            double ndf  = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->GetNDF();
////            cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
////            g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->SetPoint(iopt,vals[ires][iopt],chi2/ndf);
////
////            for(int iy = 1; iy <= 20; iy++)
////            { 
////              g_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][iy-1] = new TGraphAsymmErrors();
////            }
////            for(int ix = 1; ix <= 20; ix++)
////            {
////              g_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ix-1] = new TGraphAsymmErrors();
////              for(int iy = 1; iy <= 20; iy++)
////              { 
////                double widthx = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinWidth(ix);
////                double startx = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinLowEdge(ix);
////                double stopx  = startx+widthx;
////                double widthy = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinWidth(iy);
////                double starty = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinLowEdge(iy);
////                double stopy  = starty+widthy;
////                double value = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
////
////                double bincenterx = (startx+stopx)/2.0;
////                double bincentery = (starty+stopy)/2.0;
////
////                g_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][iy-1]->SetPoint(ix-1,bincenterx,value); 
////                g_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ix-1]->SetPoint(iy-1,bincentery,value); 
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
//  for(int ieta = 0; ieta < 10; ieta++)
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
//            //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
//            //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
//
//            h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[ires][iopt],etacuts[ieta+1]));
//            h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
//            h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{g}");
//            h_m2D_RC_Global_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
//
////            h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt],"NMR");
////            //h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt],"NMR");
////            double chi2 = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->GetChisquare();
////            double ndf  = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->GetNDF();
////            cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
////            g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->SetPoint(iopt,vals[ires][iopt],chi2/ndf);
////
////            for(int iy = 1; iy <= 20; iy++)
////            { 
////              g_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][iy-1] = new TGraphAsymmErrors();
////            }
////            for(int ix = 1; ix <= 20; ix++)
////            {
////              g_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ix-1] = new TGraphAsymmErrors();
////              for(int iy = 1; iy <= 20; iy++)
////              { 
////                double widthx = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinWidth(ix);
////                double startx = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinLowEdge(ix);
////                double stopx  = startx+widthx;
////                double widthy = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinWidth(iy);
////                double starty = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinLowEdge(iy);
////                double stopy  = starty+widthy;
////                double value = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
////
////                double bincenterx = (startx+stopx)/2.0;
////                double bincentery = (starty+stopy)/2.0;
////
////                g_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][iy-1]->SetPoint(ix-1,bincenterx,value); 
////                g_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ix-1]->SetPoint(iy-1,bincentery,value); 
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
//  for(int ieta = 0; ieta < 10; ieta++)
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
//            //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
//            //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
//
//            h_m2D_RC_Global_Corrected_BackToMC[ires][iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[ires][iopt],etacuts[ieta+1]));
//            h_m2D_RC_Global_Corrected_BackToMC[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
//            h_m2D_RC_Global_Corrected_BackToMC[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{g}");
//            h_m2D_RC_Global_Corrected_BackToMC[ires][iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
//
//            h_m2D_RC_Global_Corrected_BackToMC_E[ires][iopt][irho][i][ikin][ieta][ipt] = (TH2D*) h_m2D_RC_Global_Corrected_BackToMC[ires][iopt][irho][i][ikin][ieta][ipt]->Clone();
//            for(int ix = 1; ix <= 20; ix++)
//            {
//              for(int iy = 1; iy <= 20; iy++)
//              { 
//                double value = h_m2D_RC_Global_Corrected_BackToMC[ires][iopt][irho][i][ikin][ieta][ipt]->GetBinError(ix,iy);                
//                cout << "CORRECTED/MC = " << std::fixed << std::setprecision(15) << h_m2D_RC_Global_Corrected_BackToMC[ires][iopt][irho][i][ikin][ieta][ipt]->GetBinContent(ix,iy) << endl;
//                h_m2D_RC_Global_Corrected_BackToMC_E[ires][iopt][irho][i][ikin][ieta][ipt]->SetBinContent(ix,iy,value);
//              }
//            }
//
////            h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt],"NMR");
////            //h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt],"NMR");
////            double chi2 = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->GetChisquare();
////            double ndf  = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->GetNDF();
////            cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
////            g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->SetPoint(iopt,vals[ires][iopt],chi2/ndf);
////
////            for(int iy = 1; iy <= 20; iy++)
////            { 
////              g_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][iy-1] = new TGraphAsymmErrors();
////            }
////            for(int ix = 1; ix <= 20; ix++)
////            {
////              g_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ix-1] = new TGraphAsymmErrors();
////              for(int iy = 1; iy <= 20; iy++)
////              { 
////                double widthx = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinWidth(ix);
////                double startx = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinLowEdge(ix);
////                double stopx  = startx+widthx;
////                double widthy = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinWidth(iy);
////                double starty = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinLowEdge(iy);
////                double stopy  = starty+widthy;
////                double value = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
////
////                double bincenterx = (startx+stopx)/2.0;
////                double bincentery = (starty+stopy)/2.0;
////
////                g_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][iy-1]->SetPoint(ix-1,bincenterx,value); 
////                g_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ix-1]->SetPoint(iy-1,bincentery,value); 
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
//  for(int ieta = 0; ieta < 10; ieta++)
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
//            //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
//            //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
//
//            h_m2D_RC_Global_Corrected_BackToMC_E[ires][iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[ires][iopt],etacuts[ieta+1]));
//            h_m2D_RC_Global_Corrected_BackToMC_E[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
//            h_m2D_RC_Global_Corrected_BackToMC_E[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{g}");
//            h_m2D_RC_Global_Corrected_BackToMC_E[ires][iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
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
//  for(int ieta = 0; ieta < 10; ieta++)
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
//            //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
//            //fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->SetParameter(0,h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetMaximum());
//
//            h_m2D_RC_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[ires][iopt],etacuts[ieta+1]));
//            h_m2D_RC_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{h}");
//            h_m2D_RC_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{h}");
//            h_m2D_RC_Corrected[ires][iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
//
////            h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt],"NMR");
////            //h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Fit(fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt],"NMR");
////            double chi2 = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->GetChisquare();
////            double ndf  = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->GetNDF();
////            cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
////            g_m2D_chi2ndf[irho][i][ikin][ieta][ipt]->SetPoint(iopt,vals[ires][iopt],chi2/ndf);
////
////            for(int iy = 1; iy <= 20; iy++)
////            { 
////              g_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][iy-1] = new TGraphAsymmErrors();
////            }
////            for(int ix = 1; ix <= 20; ix++)
////            {
////              g_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ix-1] = new TGraphAsymmErrors();
////              for(int iy = 1; iy <= 20; iy++)
////              { 
////                double widthx = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinWidth(ix);
////                double startx = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->GetBinLowEdge(ix);
////                double stopx  = startx+widthx;
////                double widthy = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinWidth(iy);
////                double starty = h_m2D_RC_Global_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->GetBinLowEdge(iy);
////                double stopy  = starty+widthy;
////                double value = fitacc2D[ires][iopt][irho][i][ikin][ieta][ipt]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
////
////                double bincenterx = (startx+stopx)/2.0;
////                double bincentery = (starty+stopy)/2.0;
////
////                g_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][iy-1]->SetPoint(ix-1,bincenterx,value); 
////                g_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ix-1]->SetPoint(iy-1,bincentery,value); 
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
//        for(int ieta = 0; ieta < 12; ieta++)
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
//  for(int ieta = 0; ieta < 10; ieta++)
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
//            h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{h}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[ires][iopt],etacuts[ieta+1]));
//            h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetXaxis()->SetTitle("cos(#theta^{*}_{h}");
//            h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->GetYaxis()->SetTitle("#beta_{h}");
//            h_m2D_RC_Ratio[ires][iopt][irho][i][ikin][ieta][ipt]->Draw("Colz");
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
//  for(int ieta = 0; ieta < 10; ieta++)
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
//              h_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f, #beta%d",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[ires][iopt],etacuts[ieta+1],ibin));
//              h_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
//              h_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->GetYaxis()->SetTitle("Acceptance");
//              h_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("Colz");
//              g_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->SetLineColor(kRed);
//              g_m1D_RC_Global_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("C");
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
//  for(int ieta = 0; ieta < 10; ieta++)
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
//              h_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f, cos%d",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[ires][iopt],etacuts[ieta+1],ibin));
//              h_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->GetXaxis()->SetTitle("#beta_{g}");
//              h_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->GetYaxis()->SetTitle("Acceptance");
//              h_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("Colz");
//              g_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->SetLineColor(kRed);
//              g_m1D_RC_Global_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("C");
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
//  for(int ieta = 0; ieta < 10; ieta++)
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
//              h_m1D_RC_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f, #beta%d",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[ires][iopt],etacuts[ieta+1],ibin));
//              h_m1D_RC_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->GetXaxis()->SetTitle("cos(#theta^{*}_{h}");
//              h_m1D_RC_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->GetYaxis()->SetTitle("Acceptance");
//              h_m1D_RC_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("Colz");
//              g_m1D_RC_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->SetLineColor(kRed);
//              g_m1D_RC_Ratio_Cos[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("C");
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
//  for(int ieta = 0; ieta < 10; ieta++)
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
//              h_m1D_RC_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f, cos%d",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[ires][iopt],etacuts[ieta+1],ibin));
//              h_m1D_RC_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->GetXaxis()->SetTitle("#beta_{h}");
//              h_m1D_RC_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->GetYaxis()->SetTitle("Acceptance");
//              h_m1D_RC_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("Colz");
//              g_m1D_RC_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->SetLineColor(kRed);
//              g_m1D_RC_Ratio_Beta[ires][iopt][irho][i][ikin][ieta][ipt][ibin]->Draw("C");
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
  //    h_m2D_Weight[ires][iopt][ikin]->Write(); 
  //  } 
  //}
  //File_OutPut->Close();  
 
}

