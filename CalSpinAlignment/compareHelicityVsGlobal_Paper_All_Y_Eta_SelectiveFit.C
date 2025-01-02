#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TF2.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
//#include "../Utility/draw.h"
#include "../../../draw.h"
#include "../Utility/functions.h"
#include "resolution_pt.h"
#include "TMinuit.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

int globalbinidx[1000]  = {0};
double xlows[1000]  = {0.0};
double ylows[1000]  = {0.0};
double xhighs[1000] = {0.0};
double yhighs[1000] = {0.0};
double xvalues[1000] = {0.0};
double yvalues[1000] = {0.0};
double zvalues[1000] = {0.0};
double zerrors[1000] = {0.0};

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

double SpinDensity2DCosBeta(double x, double y, double *p)
{
  return p[5]*((1.-p[0])+(3.*p[0] -1.)*x*x
               -sqrt(2)*p[1]*2*x*TMath::Sqrt(1-x*x)*cos(y)
               +sqrt(2)*p[2]*2*x*TMath::Sqrt(1-x*x)*sin(y)
               -2.*p[3]*(1-x*x)*cos(2.*y)
               +2.*p[4]*(1-x*x)*sin(2.*y));
}


void chi2func(int &npar, double *gin, double &f, double *par, int iflag)
{
   const int nbins = 400;
   TF2* spindensity2D = new TF2("spindensity",SpinDensity2Dcos,-1.0,1.0,0.0,2*TMath::Pi(),6);

   //if (spindensity2D != nullptr) {
   //  //cout << "spindensity2D is initialized!" << endl;
   //  
   //} else {
   //  cout << "spindensity2D is not initialized!" << endl;
   //}

   for(int i = 0; i < 6; i++)
   {
     //cout << "Before spindensity"<<endl; 
     spindensity2D->SetParameter(i,par[i]);
     //cout << "After spindensity"<<endl;  
   }
   //calculate chisquare
   double chisq = 0.;
   double delta = 0.;
   for (int i = 0; i < nbins; i++) 
   {
     if(zvalues[i] <= 0.0000001 || zerrors[i] <= 0.0000001) continue;  

     // Integrate the model over the bin range
     //cout << "Before Integral"<<endl; 
     double binwidth = (xhighs[i]-xlows[i])*(yhighs[i]-ylows[i]); 
     double lambda_ij = spindensity2D->Integral(xlows[i], xhighs[i], ylows[i], yhighs[i],0.01)/binwidth;
     //cout << "After Integral"<<endl;  
     
     delta  = (zvalues[i]-lambda_ij)/zerrors[i];
     //delta  = (zvalues[i]-SpinDensity2DCosBeta(xvalues[i],yvalues[i],par))/zerrors[i];
     chisq += delta*delta;
   }
   f = chisq;
}

void compareHelicityVsGlobal_Paper_All_Y_Eta_SelectiveFit(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 23, std::string simmode = "Mc", int inputpar = 0)//defaultF = 0 is BESII, defaultF = 1 is BESI
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

  const int nkin = 1;
  const float ptkin[nkin] = {1.0};
  const float ykin[nkin]  = {0};


  //const int nopts = 2;
  //std::string option = "v2_On_Off";
  //std::string date[nopts] = {"20240809g","20240809g"}; // nov2, prelimv2
  //std::string opts[nopts] = {"nov2","prelimv2"};
  //int color[nopts] = {kBlue, kOrange+7}; 
  //std::string label[nopts] = {"v_{2} OFF","v_{2} ON (Prelim)"};

  const int nopts = 1;

  //std::string option = "y_Variance_RC_EtaCuts_GlobalInput_pt2_individualy";
  std::string option = "y_Variance_RC_EtaCuts_GlobalInput_pt1_selectivefit";
  //std::string date[nopts] = {"20240809_2","20240815he","20240815he","20240815he","20240815he","20240815he","20240815he","20240815he","20240815he","20240815he"}; // nov2, prelimv2
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

  TFile *File_InPutRC[nkin]; 
  TH2D *h_m2D_RC[nopts][5][nvar][nkin][11];  
  TH2D *h_m2D_RC_Global[nopts][5][nvar][nkin][11]; 
  //TH2D *h_m2D_RC[nopts][5][nvar][nkin][11];  
  TH1D *h_m1D_RC_Global_Cos[nopts][5][nvar][nkin][11][20]; 
  TH1D *h_m1D_RC_Global_Beta[nopts][5][nvar][nkin][11][20]; 
  TH1D *h_m1D_RC_Global_Ratio_Cos[nopts][5][nvar][nkin][11][20]; 
  TH1D *h_m1D_RC_Global_Ratio_Beta[nopts][5][nvar][nkin][11][20]; 
  TH1D *h_m1D_RC_Cos[nopts][5][nvar][nkin][11][20]; 
  TH1D *h_m1D_RC_Beta[nopts][5][nvar][nkin][11][20]; 
  TH1D *h_m1D_RC_Ratio_Cos[nopts][5][nvar][nkin][11][20]; 
  TH1D *h_m1D_RC_Ratio_Beta[nopts][5][nvar][nkin][11][20]; 

  TH2D *h_m2D_RC_Global_Ratio[nopts][5][nvar][nkin][10];
  TH2D *h_m2D_RC_Ratio[nopts][5][nvar][nkin][10];
  TH2D *h_m2D_RC_Global_Corrected[nopts][5][nvar][nkin][10];
  TH2D *h_m2D_RC_Corrected[nopts][5][nvar][nkin][10];
 
  TH2D *h_m2D_Weight[nopts][nkin];

  TGraphAsymmErrors *g_m1D_RC_Global_Ratio_Cos[nopts][5][nvar][nkin][11][20]; 
  TGraphAsymmErrors *g_m1D_RC_Global_Ratio_Beta[nopts][5][nvar][nkin][11][20]; 
  TGraphAsymmErrors *g_m1D_RC_Ratio_Cos[nopts][5][nvar][nkin][11][20]; 
  TGraphAsymmErrors *g_m1D_RC_Ratio_Beta[nopts][5][nvar][nkin][11][20]; 

  for(int ikin = 0; ikin < nkin; ikin++)
  {
    for(int iopt = 0; iopt < nopts; iopt++)
    {
      string inputfileRC = Form("effaccfiles/%s/%s/PaperAllYRC_EtaFixed_GlobalInput/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0_pt%1.1f_y%1.1f.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),ptkin[ikin],vals[iopt]);
      File_InPutRC[ikin] = TFile::Open(inputfileRC.c_str());
      for(int irho = 0; irho < 5; irho++)
      {
        for(int i = 0; i < nvar; i++)
        //for(int i = 0; i < nvar; i++)
        {
          for(int ieta = 0; ieta < 11; ieta++)
          {
            string KEY_2D_RC = Form("h_m%sEffCosPhiPrimeH_v2_6_rhoinput_%d_rho_%d_eta_%d",simmode.c_str(),irho,i,ieta);
            string KEY_2D_RC_Global = Form("h_m%sEffCosPhiPrime_v2_6_rhoinput_%d_rho_%d_eta_%d",simmode.c_str(),irho,i,ieta);
            string KEY_2D = Form("h_m%sEffCosPhiPrimeH_v2_6_rhoinput_%d_rho_%d_eta_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,iopt,ikin);
            string KEY_2D_Global = Form("h_m%sEffCosPhiPrime_v2_6_rhoinput_%d_rho_%d_eta_%d_y_%d_kin_%d",simmode.c_str(),irho,i,ieta,iopt,ikin);

            h_m2D_RC[iopt][irho][i][ikin][ieta] = (TH2D*) ((TH2D*) File_InPutRC[ikin]->Get(KEY_2D_RC.c_str()))->Clone(KEY_2D.c_str());
            //h_m2D_RC[iopt][irho][i][ikin][ieta]->Sumw2();
            h_m2D_RC_Global[iopt][irho][i][ikin][ieta] = (TH2D*) ((TH2D*) File_InPutRC[ikin]->Get(KEY_2D_RC_Global.c_str()))->Clone(KEY_2D_Global.c_str());
            //h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->Sumw2();
            for(int ibin = 0; ibin < 20; ibin++)
            {
              string KEY = Form("h_m%sEffCosPhiPrimeCos_v2_6_rhoinput_%d_rho_%d_eta_%d_y_%d_kin_%d_bin_%d",simmode.c_str(),irho,i,ieta,iopt,ikin,ibin);
              h_m1D_RC_Global_Cos[iopt][irho][i][ikin][ieta][ibin] = (TH1D*) h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->ProjectionX(KEY.c_str(),ibin+1,ibin+1);
              h_m1D_RC_Global_Cos[iopt][irho][i][ikin][ieta][ibin]->SetDirectory(0);
              KEY = Form("h_m%sEffCosPhiPrimeCosH_v2_6_rhoinput_%d_rho_%d_eta_%d_y_%d_kin_%d_bin_%d",simmode.c_str(),irho,i,ieta,iopt,ikin,ibin);
              h_m1D_RC_Cos[iopt][irho][i][ikin][ieta][ibin] = (TH1D*) h_m2D_RC[iopt][irho][i][ikin][ieta]->ProjectionX(KEY.c_str(),ibin+1,ibin+1);
              h_m1D_RC_Cos[iopt][irho][i][ikin][ieta][ibin]->SetDirectory(0);
            }
            for(int ibin = 0; ibin < 20; ibin++)
            {
              string KEY = Form("h_m%sEffCosPhiPrimeBeta_v2_6_rhoinput_%d_rho_%d_eta_%d_y_%d_kin_%d_bin_%d",simmode.c_str(),irho,i,ieta,iopt,ikin,ibin);
              h_m1D_RC_Global_Beta[iopt][irho][i][ikin][ieta][ibin] = (TH1D*) h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->ProjectionY(KEY.c_str(),ibin+1,ibin+1);
              h_m1D_RC_Global_Beta[iopt][irho][i][ikin][ieta][ibin]->SetDirectory(0);
              KEY = Form("h_m%sEffCosPhiPrimeBetaH_v2_6_rhoinput_%d_rho_%d_eta_%d_y_%d_kin_%d_bin_%d",simmode.c_str(),irho,i,ieta,iopt,ikin,ibin);
              h_m1D_RC_Beta[iopt][irho][i][ikin][ieta][ibin] = (TH1D*) h_m2D_RC[iopt][irho][i][ikin][ieta]->ProjectionY(KEY.c_str(),ibin+1,ibin+1);
              h_m1D_RC_Beta[iopt][irho][i][ikin][ieta][ibin]->SetDirectory(0);
            }

            if(ieta >= 1)
            {
              h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta-1] = (TH2D*) h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->Clone();
              //h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta-1]->Sumw2();
              h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta-1]->Divide(h_m2D_RC_Global[iopt][irho][i][ikin][ieta],h_m2D_RC_Global[iopt][irho][i][ikin][0],1,1,"B");
              h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta-1] = (TH2D*) h_m2D_RC[iopt][irho][i][ikin][ieta]->Clone();
              //h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta-1]->Sumw2();
              h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta-1]->Divide(h_m2D_RC[iopt][irho][i][ikin][ieta],h_m2D_RC[iopt][irho][i][ikin][0],1,1,"B");

              h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta-1]->SetDirectory(0);
              h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta-1]->SetDirectory(0);

              h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta-1] = (TH2D*) h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->Clone();
              //h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta-1]->Sumw2();
              h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta-1]->Divide(h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta-1]);
              h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta-1] = (TH2D*) h_m2D_RC[iopt][irho][i][ikin][ieta]->Clone();
              //h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta-1]->Sumw2();
              h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta-1]->Divide(h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta-1]);

              h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta-1]->SetDirectory(0);
              h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta-1]->SetDirectory(0);
            }

            h_m2D_RC[iopt][irho][i][ikin][ieta]->SetDirectory(0);
            h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->SetDirectory(0);

          }
        }
      }
      File_InPutRC[ikin]->Close();
      cout << inputfileRC << endl;
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
          for(int ieta = 0; ieta < 10; ieta++)
          {
            for(int ibin = 0; ibin < 20; ibin++)
            {
              h_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin] = (TH1D*) h_m1D_RC_Global_Cos[iopt][irho][i][ikin][ieta+1][ibin];
              h_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin] = (TH1D*) h_m1D_RC_Global_Beta[iopt][irho][i][ikin][ieta+1][ibin];

              h_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->Divide(h_m1D_RC_Global_Cos[iopt][irho][i][ikin][ieta+1][ibin],h_m1D_RC_Global_Cos[iopt][irho][i][ikin][0][ibin],1,1,"B");
              h_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->Divide(h_m1D_RC_Global_Beta[iopt][irho][i][ikin][ieta+1][ibin],h_m1D_RC_Global_Beta[iopt][irho][i][ikin][0][ibin],1,1,"B");

              h_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin] = (TH1D*) h_m1D_RC_Cos[iopt][irho][i][ikin][ieta+1][ibin];
              h_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin] = (TH1D*) h_m1D_RC_Beta[iopt][irho][i][ikin][ieta+1][ibin];

              h_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->Divide(h_m1D_RC_Cos[iopt][irho][i][ikin][ieta+1][ibin],h_m1D_RC_Cos[iopt][irho][i][ikin][0][ibin],1,1,"B");
              h_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->Divide(h_m1D_RC_Beta[iopt][irho][i][ikin][ieta+1][ibin],h_m1D_RC_Beta[iopt][irho][i][ikin][0][ibin],1,1,"B");
            }
          }
        }
      }
    }
  }





  TGraphAsymmErrors *g2D[5][nvar][nkin][5][nopts]; 
  TGraphAsymmErrors *g2D_Global[5][nvar][nkin][5][nopts]; 
  TGraphAsymmErrors *g2D_Global_Corrected[5][nvar][nkin][5][nopts]; 
  TGraphAsymmErrors *g2D_Corrected[5][nvar][nkin][5][nopts]; 
  TF2 *fits2D[nopts][5][nvar][nkin][11];
  TF2 *fits2D_Global[nopts][5][nvar][nkin][11];
  TF2 *fits2D_Global_Corrected[nopts][5][nvar][nkin][11];
  TF2 *fits2D_Corrected[nopts][5][nvar][nkin][11];

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int i = 0; i < nvar; i++)
    {
      for(int ikin = 0; ikin < nkin; ikin++)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          for(int par = 0; par < 5; par++) 
          {
            g2D[irho][i][ikin][par][iopt] = new TGraphAsymmErrors();
            g2D_Global[irho][i][ikin][par][iopt] = new TGraphAsymmErrors();
          }   
          for(int ieta = 1; ieta < 11; ieta++)
          {
            fits2D[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),6);
            fits2D[iopt][irho][i][ikin][ieta]->SetParameter(0,1./3.);
            fits2D[iopt][irho][i][ikin][ieta]->SetParLimits(0,0.0,1.0);
            fits2D[iopt][irho][i][ikin][ieta]->SetParameter(1,0.0);
            fits2D[iopt][irho][i][ikin][ieta]->SetParameter(2,0.0);
            fits2D[iopt][irho][i][ikin][ieta]->SetParameter(3,0.0);
            fits2D[iopt][irho][i][ikin][ieta]->SetParameter(4,0.0);
            fits2D[iopt][irho][i][ikin][ieta]->SetParameter(5,h_m2D_RC[iopt][irho][i][ikin][ieta]->GetMaximum());
            //h_m2D_RC[iopt][irho][i][ikin][ieta]->Fit(fits2D[iopt][irho][i][ikin][ieta],"NMR");
            for(int par = 0; par < 5; par++)
            {
              if(par == 0 && par == irho) g2D[irho][i][ikin][par][iopt]->SetPoint(ieta-1,etacuts[ieta], fits2D[iopt][irho][i][ikin][ieta]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
              if(par != 0 && par == irho) g2D[irho][i][ikin][par][iopt]->SetPoint(ieta-1,etacuts[ieta], fits2D[iopt][irho][i][ikin][ieta]->GetParameter(par)-(float(i)-2)*0.1);
              if(par == 0 && par != irho) g2D[irho][i][ikin][par][iopt]->SetPoint(ieta-1,etacuts[ieta], fits2D[iopt][irho][i][ikin][ieta]->GetParameter(0)-1./3.);
              if(par != 0 && par != irho) g2D[irho][i][ikin][par][iopt]->SetPoint(ieta-1,etacuts[ieta], fits2D[iopt][irho][i][ikin][ieta]->GetParameter(par));
              g2D[irho][i][ikin][par][iopt]->SetPointError(ieta-1,0.0,0.0,fits2D[iopt][irho][i][ikin][ieta]->GetParError(par),fits2D[iopt][irho][i][ikin][ieta]->GetParError(par));
            }
            cout << "Fit helicity " << endl;

            fits2D_Global[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_Global_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2Dcos,-1.0,1.0,0,2.0*TMath::Pi(),6);
            fits2D_Global[iopt][irho][i][ikin][ieta]->SetParameter(0,1./3.);
            fits2D_Global[iopt][irho][i][ikin][ieta]->SetParLimits(0,0.0,1.0);
            fits2D_Global[iopt][irho][i][ikin][ieta]->SetParameter(5,h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->GetMaximum());
            //h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->Fit(fits2D_Global[iopt][irho][i][ikin][ieta],"NMR");
            for(int par = 0; par < 5; par++)
            {
              if(par == 0 && par == irho) g2D_Global[irho][i][ikin][par][iopt]->SetPoint(ieta-1,etacuts[ieta]/*+0.002*(float(i)-float(nvar)/2.)*/, fits2D_Global[iopt][irho][i][ikin][ieta]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
              if(par != 0 && par == irho) g2D_Global[irho][i][ikin][par][iopt]->SetPoint(ieta-1,etacuts[ieta]/*+0.002*(float(i)-float(nvar)/2.)*/, fits2D_Global[iopt][irho][i][ikin][ieta]->GetParameter(par)-(float(i)-2)*0.1);
              if(par == 0 && par != irho) g2D_Global[irho][i][ikin][par][iopt]->SetPoint(ieta-1,etacuts[ieta]/*+0.002*(float(i)-float(nvar)/2.)*/, fits2D_Global[iopt][irho][i][ikin][ieta]->GetParameter(0)-1./3.);
              if(par != 0 && par != irho) g2D_Global[irho][i][ikin][par][iopt]->SetPoint(ieta-1,etacuts[ieta]/*+0.002*(float(i)-float(nvar)/2.)*/, fits2D_Global[iopt][irho][i][ikin][ieta]->GetParameter(par));
              g2D_Global[irho][i][ikin][par][iopt]->SetPointError(ieta-1,0.0,0.0,fits2D_Global[iopt][irho][i][ikin][ieta]->GetParError(par),fits2D_Global[iopt][irho][i][ikin][ieta]->GetParError(par));
            }
            cout << "Fit global " << endl;
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

  //TF1 *rhoGv2[5][5][nvar][nkin][11];

  TLegend *leg[5];// = new TLegend(0.4,0.75,0.7,0.9);

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int irho = 0; irho < 5; irho++)
      {
        if(ikin == 0 && iopt == 0) 
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
          cfit2D->cd(1+par+5*irho);
          double max = -999999;
          double min = 999999;
          for(int i = 0; i < nvar; i++)
          {
              double tmax = g2D_Global[irho][i][ikin][par][iopt]->GetHistogram()->GetMaximum();   
              double tmin = g2D_Global[irho][i][ikin][par][iopt]->GetHistogram()->GetMinimum(); 

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
          //      rhoGv2[par][irho][i][ikin][ieta] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),rhoGfromHy,-1.0,1.0,7);
          //      break;
          //    case 1: 
          //      rhoGv2[par][irho][i][ikin][ieta] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),realGfromHy,-1.0,1.0,7);
          //      break;
          //    case 2: 
          //      rhoGv2[par][irho][i][ikin][ieta] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),imagGfromHy,-1.0,1.0,7);
          //      break;
          //    case 3: 
          //      rhoGv2[par][irho][i][ikin][ieta] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),rerho1n1GfromHy,-1.0,1.0,7);
          //      break;
          //    case 4: 
          //      rhoGv2[par][irho][i][ikin][ieta] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),imrho1n1GfromHy,-1.0,1.0,7);
          //      break;
          //  }

          //  rhoGv2[par][irho][i][ikin][ieta]->SetParameter(0,0.2);      
          //  rhoGv2[par][irho][i][ikin][ieta]->SetParameter(1,double(ptkin[ikin]));      
          //  for(int ip = 0; ip < 5; ip++)
          //  {
          //    if(ip == irho)
          //    {
          //      if(ip == 0)
          //      {
          //        rhoGv2[par][irho][i][ikin][ieta]->SetParameter(2+ip,(double(i)-2)*0.1+1./3.);      
          //      }
          //      else
          //      {
          //        rhoGv2[par][irho][i][ikin][ieta]->SetParameter(2+ip,(double(i)-2)*0.1);       
          //      }
          //    }
          //    else
          //    {
          //      if(ip == 0)
          //      {
          //        rhoGv2[par][irho][i][ikin][ieta]->SetParameter(2+ip,1./3.);       
          //      }
          //      else
          //      {
          //        rhoGv2[par][irho][i][ikin][ieta]->SetParameter(2+ip,0.0);       
          //      }
          //    }
          //  }          
          //  rhoGv2[par][irho][i][ikin][ieta]->SetLineColor(color[i]);
          //  rhoGv2[par][irho][i][ikin][ieta]->SetLineStyle(1);
          //  rhoGv2[par][irho][i][ikin][ieta]->SetLineWidth(1);

            g2D_Global[irho][i][ikin][par][iopt]->SetMarkerStyle(20);
            g2D_Global[irho][i][ikin][par][iopt]->SetMarkerSize(1.0);
            g2D_Global[irho][i][ikin][par][iopt]->SetMarkerColor(color[i]);
            g2D_Global[irho][i][ikin][par][iopt]->SetLineColor(color[i]);
            //g2D_Global[irho][i][ikin][par][ieta]->SetTitle(Form("Input %s^{h}, v_{2}=0.2, p_{T}=%1.1f GeV/c, |#eta|#leq%1.1f",param[irho].c_str(),ptkin[ikin],etacuts[ieta]));
            g2D_Global[irho][i][ikin][par][iopt]->SetTitle(Form("Input %s^{g}, v_{2}=0.2, p_{T}=%1.1f GeV/c",param[irho].c_str(),ptkin[ikin]));
            g2D_Global[irho][i][ikin][par][iopt]->GetXaxis()->SetTitle(Form("Kaon |#eta|<"));
            g2D_Global[irho][i][ikin][par][iopt]->GetYaxis()->SetTitle(Form("%s^{g} from 2D Fit - Truth",param[par].c_str()));
            g2D_Global[irho][i][ikin][par][iopt]->GetXaxis()->SetRangeUser(0.25,1.0+0.05);
            if(min < 0.0 && max < 0.0) g2D_Global[irho][i][ikin][par][iopt]->GetYaxis()->SetRangeUser(1.6*min,0.4*max);
            if(min < 0.0 && max > 0.0) g2D_Global[irho][i][ikin][par][iopt]->GetYaxis()->SetRangeUser(1.6*min,1.6*max);
            if(min > 0.0 && max > 0.0) g2D_Global[irho][i][ikin][par][iopt]->GetYaxis()->SetRangeUser(0.4*min,1.6*max);
            //g2D_Global[irho][i][ikin][par]->GetXaxis()->SetLimits(-1.0,1.0);
            //g2D_Global[irho][i][ikin][par]->GetXaxis()->SetRangeUser(-1.0,1.0);

            if(i == 0) g2D_Global[irho][i][ikin][par][iopt]->Draw("APE"); 
            else       g2D_Global[irho][i][ikin][par][iopt]->Draw("PE same"); 

            //rhoGv2[par][irho][i][ikin][ieta]->Draw("l same");
              
            if(par == 0 && ikin == 0 && iopt == 0) leg[irho]->AddEntry(g2D_Global[irho][i][ikin][par][iopt],Form("%s^{h} = %1.1f",param[irho].c_str(),(float(i)-2)*0.1),"p");
          }
          leg[irho]->Draw("same");
        }
      }   
      cfit2D->Print(outputname.c_str());
      cfit2D->Update();
    }
  }
   
  cfit2D->Print(outputstop.c_str());

  outputname = Form("figures/%s/%s/pTstudy/HelicityvsGlobal_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  cfit2D->Print(outputstart.c_str());


  //TLegend *leg;// = new TLegend(0.4,0.75,0.7,0.9);
  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
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
          cfit2D->cd(1+par+5*irho);
          //cfit2D->cd(1+irho+5*par);
          //TLegend *leg = new TLegend(0.4,0.6,0.6,0.8);
          double max = -999999;
          double min = 999999;
          for(int i = 0; i < nvar; i++)
          {
            double tmax = g2D[irho][i][ikin][par][iopt]->GetHistogram()->GetMaximum();   
            double tmin = g2D[irho][i][ikin][par][iopt]->GetHistogram()->GetMinimum(); 

            if(tmax > max) max = tmax;
            if(tmin < min) min = tmin;

            //cout << "max = " << max << ", min = " << min << endl;
          }
          for(int i = 0; i < nvar; i++)
          {
            g2D[irho][i][ikin][par][iopt]->SetMarkerStyle(20);
            g2D[irho][i][ikin][par][iopt]->SetMarkerSize(1.2);
            g2D[irho][i][ikin][par][iopt]->SetMarkerColor(color[i]);
            g2D[irho][i][ikin][par][iopt]->SetLineColor(color[i]);
//            g2D[irho][i][ikin][par][iopt]->SetTitle(Form("Input %s^{Helicity}, v_{2}=0.2, p_{T}=%1.1f GeV/c, |#eta|#leq%1.1f",param[irho].c_str(),ptkin[ikin],etacuts[ieta]));
            g2D[irho][i][ikin][par][iopt]->SetTitle(Form("Input %s^{g}, v_{2}=0.2, p_{T}=%1.1f GeV/c",param[irho].c_str(),ptkin[ikin]));
            g2D[irho][i][ikin][par][iopt]->GetXaxis()->SetTitle(Form("p_{T} (GeV/c)"));
            g2D[irho][i][ikin][par][iopt]->GetYaxis()->SetTitle(Form("%s^{Helicity} from 2D Fit",param[par].c_str()));
            if(min < 0.0 && max < 0.0) g2D[irho][i][ikin][par][iopt]->GetYaxis()->SetRangeUser(1.6*min,0.4*max);
            if(min < 0.0 && max > 0.0) g2D[irho][i][ikin][par][iopt]->GetYaxis()->SetRangeUser(1.6*min,1.6*max);
            if(min > 0.0 && max > 0.0) g2D[irho][i][ikin][par][iopt]->GetYaxis()->SetRangeUser(0.4*min,1.6*max);
            g2D[irho][i][ikin][par][iopt]->GetXaxis()->SetRangeUser(0.05,1.0+0.05);

            if(i == 0) g2D[irho][i][ikin][par][iopt]->Draw("APE"); 
            else       g2D[irho][i][ikin][par][iopt]->Draw("PE same"); 
              
            //if(par == 0) leg->AddEntry(g2D[irho][i][ikin][par],Form("%s^{Helicity} = %1.1f",param[inputpar].c_str(),(float(i)-2)*0.1),"p");
          }
          leg[irho]->Draw("same");
        }
      }
      cfit2D->Print(outputname.c_str());
      cfit2D->Update();
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
 
//  outputname = Form("figures/%s/%s/pTstudy/Global2DDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c2D->Print(outputstart.c_str());
//
//  for(int ieta = 0; ieta < 11; ieta++)
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
//            h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{h}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta]));
//            h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
//            h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->GetYaxis()->SetTitle("#beta_{g}");
//            h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->Draw("Colz");
//          }
//          c2D->Print(outputname.c_str());
//          c2D->Update();
//        }
//      }
//    }
//  }
//  c2D->Print(outputstop.c_str());
//
//  outputname = Form("figures/%s/%s/pTstudy/Helicity2DDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
//  outputstart = Form("%s[",outputname.c_str()); 
//  outputstop = Form("%s]",outputname.c_str()); 
//
//  c2D->Print(outputstart.c_str());
//
//  for(int ieta = 1; ieta < 11; ieta++)
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
//            h_m2D_RC[iopt][irho][i][ikin][ieta]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{h}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta]));
//            h_m2D_RC[iopt][irho][i][ikin][ieta]->GetXaxis()->SetTitle("cos(#theta^{*}_{h}");
//            h_m2D_RC[iopt][irho][i][ikin][ieta]->GetYaxis()->SetTitle("#beta_{h}");
//            h_m2D_RC[iopt][irho][i][ikin][ieta]->Draw("Colz");
//          }
//          c2D->Print(outputname.c_str());
//          c2D->Update();
//        }
//      }
//    }
//  }
//  c2D->Print(outputstop.c_str());


  TF2* fitacc2D[nopts][5][nvar][nkin][10];
  TGraphAsymmErrors *g_m2D_chi2ndf[5][nvar][nkin][10];
  for(int ieta = 9; ieta >= 0; ieta--)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int irho = 0; irho < 5; irho++)
      {
        for(int i = 0; i < nvar; i++)
        //for(int i = 0; i < nvar; i++)
        {
          g_m2D_chi2ndf[irho][i][ikin][ieta] = new TGraphAsymmErrors();
          for(int iopt = 0; iopt < nopts; iopt++)
          {
            //cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;

            fitacc2D[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
            //fitacc2D[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple7,-1.0,1.0,0.0,2.0*TMath::Pi(),12);
            //fitacc2D[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple8,-1.0,1.0,0.0,2.0*TMath::Pi(),16);
            fitacc2D[iopt][irho][i][ikin][ieta]->SetParameter(0,h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetMaximum());
            
            if(ieta < 9)
            {
              for(int ipar = 1; ipar < 9; ipar++) fitacc2D[iopt][irho][i][ikin][ieta]->SetParameter(ipar,fitacc2D[iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
              //for(int ipar = 1; ipar < 12; ipar++) fitacc2D[iopt][irho][i][ikin][ieta]->SetParameter(ipar,fitacc2D[iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
              //for(int ipar = 1; ipar < 16; ipar++) fitacc2D[iopt][irho][i][ikin][ieta]->SetParameter(ipar,fitacc2D[iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
            }
            //h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->Fit(fitacc2D[iopt][irho][i][ikin][ieta],"NMR");
 
            double chi2 = fitacc2D[iopt][irho][i][ikin][ieta]->GetChisquare();
            double ndf  = fitacc2D[iopt][irho][i][ikin][ieta]->GetNDF();

            //cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
            g_m2D_chi2ndf[irho][i][ikin][ieta]->SetPoint(iopt,vals[iopt],chi2/ndf);

            for(int iy = 1; iy <= 20; iy++)
            { 
              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][iy-1] = new TGraphAsymmErrors();
            }
            for(int ix = 1; ix <= 20; ix++)
            {
              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ix-1] = new TGraphAsymmErrors();
              for(int iy = 1; iy <= 20; iy++)
              { 
                double widthx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinWidth(ix);
                double startx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinLowEdge(ix);
                double stopx  = startx+widthx;
                double widthy = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinWidth(iy);
                double starty = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinLowEdge(iy);
                double stopy  = starty+widthy;
                double value = fitacc2D[iopt][irho][i][ikin][ieta]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
                
                double bincenterx = (startx+stopx)/2.0;
                double bincentery = (starty+stopy)/2.0;
    
                if(h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetBinContent(ix,iy) <= 0.0) 
                //cout << "RATIO IS "<<h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetBinContent(ix,iy)<<"+/-"<<h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetBinError(ix,iy)<<"!?   iopt="<<iopt<<",irho="<<irho<<",i="<<i<<",ikin"<<ikin<<",ieta"<<ieta<<endl;

                g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][iy-1]->SetPoint(ix-1,bincenterx,value); 
                g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ix-1]->SetPoint(iy-1,bincentery,value); 

              }
            }
          }
        }
      }
    }
  }

  TF2* fitacc2DH[nopts][5][nvar][nkin][10];
  TGraphAsymmErrors *g_m2DH_chi2ndf[5][nvar][nkin][10];
  for(int ieta = 9; ieta >= 0; ieta--)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int irho = 0; irho < 5; irho++)
      {
        for(int i = 0; i < nvar; i++)
        //for(int i = 0; i < nvar; i++)
        {
          g_m2DH_chi2ndf[irho][i][ikin][ieta] = new TGraphAsymmErrors();
          for(int iopt = 0; iopt < nopts; iopt++)
          {
            //cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;

            fitacc2DH[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
            //fitacc2D[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple7,-1.0,1.0,0.0,2.0*TMath::Pi(),12);
            fitacc2DH[iopt][irho][i][ikin][ieta]->SetParameter(0,h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta]->GetMaximum());
            
            if(ieta < 9)
            {
              for(int ipar = 1; ipar < 9; ipar++) fitacc2DH[iopt][irho][i][ikin][ieta]->SetParameter(ipar,fitacc2DH[iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
              //for(int ipar = 1; ipar < 12; ipar++) fitacc2D[iopt][irho][i][ikin][ieta]->SetParameter(ipar,fitacc2D[iopt][irho][i][ikin][ieta+1]->GetParameter(ipar));
            }
            //h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta]->Fit(fitacc2DH[iopt][irho][i][ikin][ieta],"NMR");

            double chi2 = fitacc2DH[iopt][irho][i][ikin][ieta]->GetChisquare();
            double ndf  = fitacc2DH[iopt][irho][i][ikin][ieta]->GetNDF();

            //cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
            g_m2DH_chi2ndf[irho][i][ikin][ieta]->SetPoint(iopt,vals[iopt],chi2/ndf);

            for(int iy = 1; iy <= 20; iy++)
            { 
              g_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][iy-1] = new TGraphAsymmErrors();
            }
            for(int ix = 1; ix <= 20; ix++)
            {
              g_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ix-1] = new TGraphAsymmErrors();
              for(int iy = 1; iy <= 20; iy++)
              { 
                double widthx = h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinWidth(ix);
                double startx = h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinLowEdge(ix);
                double stopx  = startx+widthx;
                double widthy = h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinWidth(iy);
                double starty = h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinLowEdge(iy);
                double stopy  = starty+widthy;
                double value = fitacc2DH[iopt][irho][i][ikin][ieta]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                

                double bincenterx = (startx+stopx)/2.0;
                double bincentery = (starty+stopy)/2.0;

                g_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][iy-1]->SetPoint(ix-1,bincenterx,value); 
                g_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ix-1]->SetPoint(iy-1,bincentery,value); 

              }
            }
          }
        }
      }
    }
  }

  //for(int iopt = 0; iopt < nopts; iopt++)
  //{
  //  for(int i = 0; i < nvar; i++)
  //  {
  //    for(int ikin = 0; ikin < nkin; ikin++)
  //    {
  //      for(int irho = 0; irho < 5; irho++)
  //      {
  //        for(int par = 0; par < 5; par++) 
  //        {
  //          g2D_Global_Corrected[irho][i][ikin][par][iopt] = new TGraphAsymmErrors();
  //        }   
  //        for(int ieta = 1; ieta < 11; ieta++)
  //        {
  //          fits2D_Global_Corrected[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2DcosAcc,-1.0,1.0,0.0,2.0*TMath::Pi(),14);
  //          //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2DcosAccLonger,-1.0,1.0,0.0,2.0*TMath::Pi(),17);
  //          //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2DcosAccLongerer,-1.0,1.0,0.0,2.0*TMath::Pi(),21);
  //          fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParameter(0,1./3.);
  //          fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParLimits(0,0.0,1.0);
  //          fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParameter(1,0.0);
  //          fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParameter(2,0.0);
  //          fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParameter(3,0.0);
  //          fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParameter(4,0.0);
  //          fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParameter(5,h_m2D_RC[iopt][irho][i][ikin][ieta]->GetMaximum());
  //          for(int ifit = 1; ifit < 9; ifit++)
  //          //for(int ifit = 1; ifit < 12; ifit++)
  //          //for(int ifit = 1; ifit < 16; ifit++)
  //          {
  //            double value = fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(ifit);
  //            double error = fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParError(ifit);
  //            fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->FixParameter(5+ifit,value);
  //            //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParLimits(5+ifit,value-error,value+error);
  //          }
  //          //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->FixParameter(7,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(2));
  //          //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->FixParameter(8,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(3));
  //          //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->FixParameter(9,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(4));
  //          //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->FixParameter(10,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(5));
  //          //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->FixParameter(11,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(6));
  //          //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->FixParameter(12,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(7));
  //          //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->FixParameter(13,fitacc2D[iopt][irho][i][ikin][ieta-1]->GetParameter(8));
  //         
  //          //h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->Fit(fits2D_Global_Corrected[iopt][irho][i][ikin][ieta],"NMR");
  //          for(int par = 0; par < 5; par++)
  //          {
  //            if(par == 0 && par == irho) g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPoint(ieta,etacuts[ieta]/*+0.002*(float(i)-float(nvar)/2.)*/, fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->GetParameter(0)-1./3.-(float(i)-2)*0.1);
  //            if(par != 0 && par == irho) g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPoint(ieta,etacuts[ieta]/*+0.002*(float(i)-float(nvar)/2.)*/, fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->GetParameter(par)-(float(i)-2)*0.1);
  //            if(par == 0 && par != irho) g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPoint(ieta,etacuts[ieta]/*+0.002*(float(i)-float(nvar)/2.)*/, fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->GetParameter(0)-1./3.);
  //            if(par != 0 && par != irho) g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPoint(ieta,etacuts[ieta]/*+0.002*(float(i)-float(nvar)/2.)*/, fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->GetParameter(par));
  //            g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPointError(ieta,0.0,0.0,fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->GetParError(par),fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->GetParError(par));
  //          }
  //          cout << "Fit global " << endl;
  //        }
  //      }
  //    }
  //  }
  //}


  outputname = Form("figures/%s/%s/pTstudy/Global2DRatioDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  c2D->Print(outputstart.c_str());

  for(int ieta = 0; ieta < 10; ieta++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int irho = 0; irho < 5; irho++)
      {
        for(int i = 0; i < nvar; i++)
        //for(int i = 0; i < nvar; i++)
        {
          //g_m2D_chi2ndf[irho][i][ikin][ieta] = new TGraphAsymmErrors();
          for(int iopt = 0; iopt < nopts; iopt++)
          {
            //cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;
            c2D->cd(1+iopt);

            //fitacc2D[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
            //fitacc2D[iopt][irho][i][ikin][ieta]->SetParameter(0,h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetMaximum());

            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1]));
            //cout << "Set title" << endl;
            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetYaxis()->SetTitle("#beta_{g}");
            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->Draw("Colz");

//            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->Fit(fitacc2D[iopt][irho][i][ikin][ieta],"NMR");
//            //h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->Fit(fitacc2D[iopt][irho][i][ikin][ieta],"NMR");
//            double chi2 = fitacc2D[iopt][irho][i][ikin][ieta]->GetChisquare();
//            double ndf  = fitacc2D[iopt][irho][i][ikin][ieta]->GetNDF();
//            cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
//            g_m2D_chi2ndf[irho][i][ikin][ieta]->SetPoint(iopt,vals[iopt],chi2/ndf);
//
//            for(int iy = 1; iy <= 20; iy++)
//            { 
//              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][iy-1] = new TGraphAsymmErrors();
//            }
//            for(int ix = 1; ix <= 20; ix++)
//            {
//              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ix-1] = new TGraphAsymmErrors();
//              for(int iy = 1; iy <= 20; iy++)
//              { 
//                double widthx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinWidth(ix);
//                double startx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinLowEdge(ix);
//                double stopx  = startx+widthx;
//                double widthy = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinWidth(iy);
//                double starty = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinLowEdge(iy);
//                double stopy  = starty+widthy;
//                double value = fitacc2D[iopt][irho][i][ikin][ieta]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
//
//                double bincenterx = (startx+stopx)/2.0;
//                double bincentery = (starty+stopy)/2.0;
//
//                g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][iy-1]->SetPoint(ix-1,bincenterx,value); 
//                g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ix-1]->SetPoint(iy-1,bincentery,value); 
//
//              }
//            }
            
          }
          c2D->Print(outputname.c_str());
          c2D->Update();
        }
      }
    }
  }
  c2D->Print(outputstop.c_str());

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int i = 0; i < nvar; i++)
    {
      for(int ikin = 0; ikin < nkin; ikin++)
      {
        for(int irho = 0; irho < 5; irho++)
        {
          for(int par = 0; par < 5; par++) 
          {
            g2D_Global_Corrected[irho][i][ikin][par][iopt] = new TGraphAsymmErrors();
          }   
          for(int ieta = 0; ieta < 10; ieta++)
          {
            //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_Corrected_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),6);
            //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParameter(0,1./3.);
            //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParLimits(0,0.0,1.0);
            //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParameter(1,0.0);
            //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParameter(2,0.0);
            //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParameter(3,0.0);
            //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParameter(4,0.0);
            //fits2D_Global_Corrected[iopt][irho][i][ikin][ieta]->SetParameter(5,h_m2D_RC_Global[iopt][irho][i][ikin][ieta]->GetMaximum());
 
            // 2D Fit with certain selections missing
            // Define the number of bins
            const int nbinsx = h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetNbinsX();
            const int nbinsy = h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetNbinsY();
            //const int nbinstot = nbinsx*nbinsy;

            int idx = 0;
 
            for(int ix = 1; ix <= nbinsx; ix++)
            {
              for(int iy = 1; iy <= nbinsy; iy++)
              {
                int globalbin = h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetBin(ix,iy);
                //double xvalue = h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinCenter(ix);
                //double yvalue = h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinCenter(iy);
                double xlow = h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinLowEdge(ix);
                double ylow = h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinLowEdge(iy);
                double xhigh = h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinWidth(ix)+xlow;
                double yhigh = h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinWidth(iy)+ylow;
                double zvalue = h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetBinContent(ix,iy);   
                double zerror = h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetBinError(ix,iy);   
     
                //xlows[globalbin-1] = xlow;
                //ylows[globalbin-1] = ylow;
                //xhighs[globalbin-1] = xhigh;
                //yhighs[globalbin-1] = yhigh;
                //zvalues[globalbin-1] = zvalue;
                //zerrors[globalbin-1] = zerror;

                //cout << "GLobal bin = " << globalbin << endl;
                //cout << " xlows[globalbin-1]  = " << xlows[globalbin-1]  << endl;  
                //cout << " ylows[globalbin-1]  = " << ylows[globalbin-1]  << endl;
                //cout << " xhighs[globalbin-1] = " << xhighs[globalbin-1] << endl;
                //cout << " yhighs[globalbin-1] = " << yhighs[globalbin-1] << endl;

                xlows[idx] = xlow;
                ylows[idx] = ylow;
                xhighs[idx] = xhigh;
                yhighs[idx] = yhigh;
                zvalues[idx] = zvalue;
                zerrors[idx] = zerror;

                //cout << "GLobal bin = " << globalbin << endl;
                //cout << " xlows[idx]  = " << xlows[idx]  << endl;  
                //cout << " ylows[idx]  = " << ylows[idx]  << endl;
                //cout << " xhighs[idx] = " << xhighs[idx] << endl;
                //cout << " yhighs[idx] = " << yhighs[idx] << endl;
                idx++;               

              }
            }     
            cout << "Final idx = " << idx << endl;
            //Create a TMinuit with a maximum of 6 parameters         
            TMinuit *gMinuit = new TMinuit(6);
            gMinuit->SetFCN(chi2func); // This is the function we minimize

            double arglist[10];
            int ierflg = 0;

            arglist[0] = 1;
            gMinuit->mnexcm("SET ERR",arglist,1,ierflg);

            cout << "Maximum should be set to " << h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetMaximum() << endl;

            //Set Starting values and step sizes for parameters
            static double startpar[6] = {1./3.,0.,0.,0.,0.,h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetMaximum()};
            static double step[6] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.001};
            for(int iparm = 0; iparm < 6; iparm++)
            {
              if(iparm == 0) gMinuit->mnparm(iparm,Form("p%d",iparm),startpar[iparm],step[iparm],0,1,ierflg);
              gMinuit->mnparm(iparm,Form("p%d",iparm),startpar[iparm],step[iparm],0,0,ierflg);
            }


            //Now we minimize
            ierflg = 0;
            arglist[0] = 10000;
            arglist[1] = 0.1;
            cout << "About to minimize" << endl;
            gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
               

            //Print results
            double fmin, fedm, errdef;
            int npari, nparx, istat;
            gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);
            std::cout << "Minimum function value (fmin): " << fmin << std::endl;
            std::cout << "Estimated distance to minimum (fedm): " << fedm << std::endl;
            std::cout << "Error definition (errdef): " << errdef << std::endl;
            std::cout << "Number of free parameters (npari): " << npari << std::endl;
            std::cout << "Total number of parameters (nparx): " << nparx << std::endl;
            std::cout << "Covariance matrix status (istat): " << istat << std::endl;

            if(istat == 3)
            {
              ierflg = 0;
              arglist[0] = 2; // strategy level (2 for robust)
              gMinuit->mnexcm("SET STR",arglist,1,ierflg);

              ierflg = 0;
              arglist[0] = 10000;
              arglist[1] = 0.00001;
              gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);

              gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);

              if(istat == 3)
              {
                ierflg = 0;
                arglist[0] = 3; // strategy level (2 for robust)
                gMinuit->mnexcm("SET STR",arglist,1,ierflg);

                ierflg = 0;
                arglist[0] = 10000;
                arglist[1] = 0.000000001;
                gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);

                gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);
              }
            }

            //if(istat == 4) gMinuit->mnexcm("MINOS",arglist,0,ierflg);
            //if(istat == 4) gMinuit->mnexcm("HESSE",arglist,0,ierflg);
 
            // Second pass
            static double startpar2[6] = {1./3.,0.,0.,0.,0.,h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetMaximum()};
            static double errorpar2[6] = {0.0};
            static double step2[6] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};
            for(int iparm = 0; iparm < 6; iparm++)
            {
              gMinuit->GetParameter(iparm,startpar2[iparm],errorpar2[iparm]);
              cout << "p" << iparm << " = " << startpar2[iparm] << " +/- " << errorpar2[iparm] << endl;
              //gMinuit->mnparm(iparm,Form("p%d",iparm),startpar2[iparm],step2[iparm],0,0,ierflg);
            }


            ////Now we minimize
            //arglist[0] = 1000;
            //arglist[1] = 1.;
            //gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
         
            ////Print results
            //gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);
            //std::cout << "Minimum function value (fmin): " << fmin << std::endl;
            //std::cout << "Estimated distance to minimum (fedm): " << fedm << std::endl;
            //std::cout << "Error definition (errdef): " << errdef << std::endl;
            //std::cout << "Number of free parameters (npari): " << npari << std::endl;
            //std::cout << "Total number of parameters (nparx): " << nparx << std::endl;
            //std::cout << "Covariance matrix status (istat): " << istat << std::endl;

            int etastart = 0;//ieta;
            for(int par = 0; par < 5; par++)
            {
              if(errorpar2[par] < TMath::Sqrt(2.)+0.00001 && errorpar2[par] > TMath::Sqrt(2.)-0.00001)
              {
                errorpar2[par] = 0.0;
                etastart = ieta;
                //startpar2[par] = 0.0;
               // if(par == 0 && par == irho) g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPoint(ieta,etacuts[ieta+1],startpar2[0]-1./3.);
               // if(par != 0 && par == irho) g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPoint(ieta,etacuts[ieta+1],startpar2[par]);
               // if(par == 0 && par != irho) g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPoint(ieta,etacuts[ieta+1],startpar2[0]-1./3.);
               // if(par != 0 && par != irho) g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPoint(ieta,etacuts[ieta+1],startpar2[par]);
              }
              else
              {
                if(par == 0 && par == irho) g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPoint(ieta-etastart,etacuts[ieta+1],startpar2[0]-1./3.-(float(i)-2)*0.1);
                if(par != 0 && par == irho) g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPoint(ieta-etastart,etacuts[ieta+1],startpar2[par]-(float(i)-2)*0.1);
                if(par == 0 && par != irho) g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPoint(ieta-etastart,etacuts[ieta+1],startpar2[0]-1./3.);
                if(par != 0 && par != irho) g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPoint(ieta-etastart,etacuts[ieta+1],startpar2[par]);
              }
              g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetPointError(ieta,0.0,0.0,errorpar2[par],errorpar2[par]);
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

  //TF1 *rhoGv2[5][5][nvar][nkin][11];

  //TLegend *leg2[5];// = new TLegend(0.4,0.75,0.7,0.9);

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int irho = 0; irho < 5; irho++)
      {
        //if(ikin == 0) 
        //{
        //  leg2[irho] = new TLegend(0.25,0.775,0.85,0.865);
        //  leg2[irho]->SetNColumns(2);
        //  leg2[irho]->SetBorderSize(0);
        //}
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
          cfit2D->cd(1+par+5*irho);
          double max = -999999;
          double min = 999999;
          for(int i = 0; i < nvar; i++)
          {
              double tmax = g2D_Global_Corrected[irho][i][ikin][par][iopt]->GetHistogram()->GetMaximum();   
              double tmin = g2D_Global_Corrected[irho][i][ikin][par][iopt]->GetHistogram()->GetMinimum(); 

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
          //      rhoGv2[par][irho][i][ikin][ieta] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),rhoGfromHy,-1.0,1.0,7);
          //      break;
          //    case 1: 
          //      rhoGv2[par][irho][i][ikin][ieta] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),realGfromHy,-1.0,1.0,7);
          //      break;
          //    case 2: 
          //      rhoGv2[par][irho][i][ikin][ieta] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),imagGfromHy,-1.0,1.0,7);
          //      break;
          //    case 3: 
          //      rhoGv2[par][irho][i][ikin][ieta] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),rerho1n1GfromHy,-1.0,1.0,7);
          //      break;
          //    case 4: 
          //      rhoGv2[par][irho][i][ikin][ieta] = new TF1(Form("rhoGv2_%d_%d_%d_%d_%d",par,irho,i,ikin,ieta),imrho1n1GfromHy,-1.0,1.0,7);
          //      break;
          //  }

          //  rhoGv2[par][irho][i][ikin][ieta]->SetParameter(0,0.2);      
          //  rhoGv2[par][irho][i][ikin][ieta]->SetParameter(1,double(ptkin[ikin]));      
          //  for(int ip = 0; ip < 5; ip++)
          //  {
          //    if(ip == irho)
          //    {
          //      if(ip == 0)
          //      {
          //        rhoGv2[par][irho][i][ikin][ieta]->SetParameter(2+ip,(double(i)-2)*0.1+1./3.);      
          //      }
          //      else
          //      {
          //        rhoGv2[par][irho][i][ikin][ieta]->SetParameter(2+ip,(double(i)-2)*0.1);       
          //      }
          //    }
          //    else
          //    {
          //      if(ip == 0)
          //      {
          //        rhoGv2[par][irho][i][ikin][ieta]->SetParameter(2+ip,1./3.);       
          //      }
          //      else
          //      {
          //        rhoGv2[par][irho][i][ikin][ieta]->SetParameter(2+ip,0.0);       
          //      }
          //    }
          //  }          
          //  rhoGv2[par][irho][i][ikin][ieta]->SetLineColor(color[i]);
          //  rhoGv2[par][irho][i][ikin][ieta]->SetLineStyle(1);
          //  rhoGv2[par][irho][i][ikin][ieta]->SetLineWidth(1);

            g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetMarkerStyle(20);
            g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetMarkerSize(1.0);
            g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetMarkerColor(color[i]);
            g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetLineColor(color[i]);
            //g2D_Global[irho][i][ikin][par][ieta]->SetTitle(Form("Input %s^{h}, v_{2}=0.2, p_{T}=%1.1f GeV/c, |#eta|#leq%1.1f",param[irho].c_str(),ptkin[ikin],etacuts[ieta]));
            g2D_Global_Corrected[irho][i][ikin][par][iopt]->SetTitle(Form("Input %s^{g}, v_{2}=0.2, p_{T}=%1.1fGeV/,c, y=%1.1f",param[irho].c_str(),ptkin[ikin],vals[iopt]));
            g2D_Global_Corrected[irho][i][ikin][par][iopt]->GetXaxis()->SetTitle(Form("Kaon |#eta|<"));
            g2D_Global_Corrected[irho][i][ikin][par][iopt]->GetYaxis()->SetTitle(Form("Corrected %s^{g} from 2D Fit - Truth",param[par].c_str()));
            g2D_Global_Corrected[irho][i][ikin][par][iopt]->GetXaxis()->SetRangeUser(0.05,1.0+0.05);
            if(min < 0.0 && max < 0.0) g2D_Global_Corrected[irho][i][ikin][par][iopt]->GetYaxis()->SetRangeUser(1.6*min,0.4*max);
            if(min < 0.0 && max > 0.0) g2D_Global_Corrected[irho][i][ikin][par][iopt]->GetYaxis()->SetRangeUser(1.6*min,1.6*max);
            if(min > 0.0 && max > 0.0) g2D_Global_Corrected[irho][i][ikin][par][iopt]->GetYaxis()->SetRangeUser(0.4*min,1.6*max);
            //g2D_Global_Corrected[irho][i][ikin][par][iopt]->GetYaxis()->SetRangeUser(-0.3,0.3);
            //g2D_Global[irho][i][ikin][par]->GetXaxis()->SetLimits(-1.0,1.0);
            //g2D_Global[irho][i][ikin][par]->GetXaxis()->SetRangeUser(-1.0,1.0);

            if(i == 0) g2D_Global_Corrected[irho][i][ikin][par][iopt]->Draw("APE"); 
            else       g2D_Global_Corrected[irho][i][ikin][par][iopt]->Draw("PE same"); 

            //rhoGv2[par][irho][i][ikin][ieta]->Draw("l same");
              
            //if(par == 0 && ikin == 0 && iopt == 0) leg2[irho]->AddEntry(g2D_Global_Corrected[irho][i][ikin][par][iopt],Form("%s^{g} = %1.1f",param[irho].c_str(),(float(i)-2)*0.1),"p");
          }
          leg[irho]->Draw("same");
        }
      }   
      cfit2D->Print(outputname.c_str());
      cfit2D->Update();
    }
  }
   
  cfit2D->Print(outputstop.c_str());

  outputname = Form("figures/%s/%s/pTstudy/Global2DCorrectedDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  c2D->Print(outputstart.c_str());

  for(int ieta = 0; ieta < 10; ieta++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int irho = 0; irho < 5; irho++)
      {
        for(int i = 0; i < nvar; i++)
        //for(int i = 0; i < nvar; i++)
        {
          //g_m2D_chi2ndf[irho][i][ikin][ieta] = new TGraphAsymmErrors();
          for(int iopt = 0; iopt < nopts; iopt++)
          {
            cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;
            c2D->cd(1+iopt);

            //fitacc2D[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
            //fitacc2D[iopt][irho][i][ikin][ieta]->SetParameter(0,h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetMaximum());

            h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1]));
            h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
            h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->GetYaxis()->SetTitle("#beta_{g}");
            h_m2D_RC_Global_Corrected[iopt][irho][i][ikin][ieta]->Draw("Colz");

//            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->Fit(fitacc2D[iopt][irho][i][ikin][ieta],"NMR");
//            //h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->Fit(fitacc2D[iopt][irho][i][ikin][ieta],"NMR");
//            double chi2 = fitacc2D[iopt][irho][i][ikin][ieta]->GetChisquare();
//            double ndf  = fitacc2D[iopt][irho][i][ikin][ieta]->GetNDF();
//            cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
//            g_m2D_chi2ndf[irho][i][ikin][ieta]->SetPoint(iopt,vals[iopt],chi2/ndf);
//
//            for(int iy = 1; iy <= 20; iy++)
//            { 
//              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][iy-1] = new TGraphAsymmErrors();
//            }
//            for(int ix = 1; ix <= 20; ix++)
//            {
//              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ix-1] = new TGraphAsymmErrors();
//              for(int iy = 1; iy <= 20; iy++)
//              { 
//                double widthx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinWidth(ix);
//                double startx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinLowEdge(ix);
//                double stopx  = startx+widthx;
//                double widthy = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinWidth(iy);
//                double starty = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinLowEdge(iy);
//                double stopy  = starty+widthy;
//                double value = fitacc2D[iopt][irho][i][ikin][ieta]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
//
//                double bincenterx = (startx+stopx)/2.0;
//                double bincentery = (starty+stopy)/2.0;
//
//                g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][iy-1]->SetPoint(ix-1,bincenterx,value); 
//                g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ix-1]->SetPoint(iy-1,bincentery,value); 
//
//              }
//            }
            
          }
          c2D->Print(outputname.c_str());
          c2D->Update();
        }
      }
    }
  }
  c2D->Print(outputstop.c_str());

  outputname = Form("figures/%s/%s/pTstudy/Helicity2DCorrectedDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  c2D->Print(outputstart.c_str());

  for(int ieta = 0; ieta < 10; ieta++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int irho = 0; irho < 5; irho++)
      {
        for(int i = 0; i < nvar; i++)
        //for(int i = 0; i < nvar; i++)
        {
          //g_m2D_chi2ndf[irho][i][ikin][ieta] = new TGraphAsymmErrors();
          for(int iopt = 0; iopt < nopts; iopt++)
          {
            cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;
            c2D->cd(1+iopt);

            //fitacc2D[iopt][irho][i][ikin][ieta] = new TF2(Form("fit2D_%d_%d_%d_%d_%d",iopt,irho,i,ikin,ieta),acceptance2Dsimple4,-1.0,1.0,0.0,2.0*TMath::Pi(),9);
            //fitacc2D[iopt][irho][i][ikin][ieta]->SetParameter(0,h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetMaximum());

            h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1]));
            h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta]->GetXaxis()->SetTitle("cos(#theta^{*}_{h}");
            h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta]->GetYaxis()->SetTitle("#beta_{h}");
            h_m2D_RC_Corrected[iopt][irho][i][ikin][ieta]->Draw("Colz");

//            h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->Fit(fitacc2D[iopt][irho][i][ikin][ieta],"NMR");
//            //h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->Fit(fitacc2D[iopt][irho][i][ikin][ieta],"NMR");
//            double chi2 = fitacc2D[iopt][irho][i][ikin][ieta]->GetChisquare();
//            double ndf  = fitacc2D[iopt][irho][i][ikin][ieta]->GetNDF();
//            cout << "chi2 = " << chi2 << ", ndf = " << ndf << endl;
//            g_m2D_chi2ndf[irho][i][ikin][ieta]->SetPoint(iopt,vals[iopt],chi2/ndf);
//
//            for(int iy = 1; iy <= 20; iy++)
//            { 
//              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][iy-1] = new TGraphAsymmErrors();
//            }
//            for(int ix = 1; ix <= 20; ix++)
//            {
//              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ix-1] = new TGraphAsymmErrors();
//              for(int iy = 1; iy <= 20; iy++)
//              { 
//                double widthx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinWidth(ix);
//                double startx = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetXaxis()->GetBinLowEdge(ix);
//                double stopx  = startx+widthx;
//                double widthy = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinWidth(iy);
//                double starty = h_m2D_RC_Global_Ratio[iopt][irho][i][ikin][ieta]->GetYaxis()->GetBinLowEdge(iy);
//                double stopy  = starty+widthy;
//                double value = fitacc2D[iopt][irho][i][ikin][ieta]->Integral(startx,stopx,starty,stopy)/widthx/widthy;                
//
//                double bincenterx = (startx+stopx)/2.0;
//                double bincentery = (starty+stopy)/2.0;
//
//                g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][iy-1]->SetPoint(ix-1,bincenterx,value); 
//                g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ix-1]->SetPoint(iy-1,bincentery,value); 
//
//              }
//            }
            
          }
          c2D->Print(outputname.c_str());
          c2D->Update();
        }
      }
    }
  }
  c2D->Print(outputstop.c_str());


  outputname = Form("figures/%s/%s/pTstudy/Global2DRatioDistributionsCHI2NDF_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  c2D->Print(outputstart.c_str());

  for(int ikin = 0; ikin < nkin; ikin++)
  {
    for(int irho = 0; irho < 5; irho++)
    {
      for(int i = 0; i < nvar; i++)
      //for(int i = 0; i < nvar; i++)
      {
        for(int ieta = 0; ieta < 12; ieta++)
        {
          //cout << "ieta = " << ieta << ", ikin = " << ikin << ", irho = " << irho << ", i = " << i << ", iopt = " << iopt << endl;
          c2D->cd(1+ieta);
          c2D->cd(1+ieta)->Clear();
          if(ieta > 9/* || ieta < 2*/) continue;

          g_m2D_chi2ndf[irho][i][ikin][ieta]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,etacuts[ieta+1]));
          g_m2D_chi2ndf[irho][i][ikin][ieta]->GetXaxis()->SetTitle("y");
          g_m2D_chi2ndf[irho][i][ikin][ieta]->GetYaxis()->SetTitle("#chi^2/NDF");
          g_m2D_chi2ndf[irho][i][ikin][ieta]->Draw("APE");
        }
        c2D->Print(outputname.c_str());
        c2D->Update();
      }
    }
  }
  c2D->Print(outputstop.c_str());


  outputname = Form("figures/%s/%s/pTstudy/Helicity2DRatioDistributions_PAPERALL_parametercomparison_%s_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str(),option.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  c2D->Print(outputstart.c_str());

  for(int ieta = 0; ieta < 10; ieta++)
  {
    for(int ikin = 0; ikin < nkin; ikin++)
    {
      for(int irho = 0; irho < 5; irho++)
      {
        //for(int i = 0; i < nvar; i++)
        for(int i = 0; i < nvar; i++)
        {
          for(int iopt = 0; iopt < nopts; iopt++)
          {
            c2D->cd(1+iopt);
            h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{h}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1]));
            h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta]->GetXaxis()->SetTitle("cos(#theta^{*}_{h}");
            h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta]->GetYaxis()->SetTitle("#beta_{h}");
            h_m2D_RC_Ratio[iopt][irho][i][ikin][ieta]->Draw("Colz");
          }
          c2D->Print(outputname.c_str());
          c2D->Update();
        }
      }
    }
  }
  c2D->Print(outputstop.c_str());

  TCanvas *c1D = new TCanvas("c2D","c2D",10,10,2000,1600);
  c1D->Divide(5,4);
  
  for(int i = 0; i < 20; i++)
  {
    c1D->cd(i+1);
    c1D->cd(i+1)->SetLeftMargin(0.15);  
    c1D->cd(i+1)->SetBottomMargin(0.15);
    c1D->cd(i+1)->SetTicks(1,1);
    c1D->cd(i+1)->SetGrid(0,0);
  }

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
//              h_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f, #beta%d",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1],ibin));
//              h_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->GetXaxis()->SetTitle("cos(#theta^{*}_{g}");
//              h_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->GetYaxis()->SetTitle("Acceptance");
//              h_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->Draw("Colz");
//              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->SetLineColor(kRed);
//              g_m1D_RC_Global_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->Draw("C");
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
//              h_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f, cos%d",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1],ibin));
//              h_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->GetXaxis()->SetTitle("#beta_{g}");
//              h_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->GetYaxis()->SetTitle("Acceptance");
//              h_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->Draw("Colz");
//              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->SetLineColor(kRed);
//              g_m1D_RC_Global_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->Draw("C");
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
//              h_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f, #beta%d",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1],ibin));
//              h_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->GetXaxis()->SetTitle("cos(#theta^{*}_{h}");
//              h_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->GetYaxis()->SetTitle("Acceptance");
//              h_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->Draw("Colz");
//              g_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->SetLineColor(kRed);
//              g_m1D_RC_Ratio_Cos[iopt][irho][i][ikin][ieta][ibin]->Draw("C");
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
//              h_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->SetTitle(Form("p_{T}=%1.1f GeV/c, Input %s^{g}=%1.1f, v_{2}=0.2, y=%1.1f, |#eta|#leq%1.1f, cos%d",ptkin[ikin],param[irho].c_str(),(float(i)-2)*0.1,vals[iopt],etacuts[ieta+1],ibin));
//              h_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->GetXaxis()->SetTitle("#beta_{h}");
//              h_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->GetYaxis()->SetTitle("Acceptance");
//              h_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->Draw("Colz");
//              g_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->SetLineColor(kRed);
//              g_m1D_RC_Ratio_Beta[iopt][irho][i][ikin][ieta][ibin]->Draw("C");
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

