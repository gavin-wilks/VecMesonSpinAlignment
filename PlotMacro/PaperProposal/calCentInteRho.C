#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include <TString.h>

void calCentInteRho()
{
  int const NumBeamEnergy = 4;
  std::string const mBeamEnergy[NumBeamEnergy] = {"27GeV_Run18","39GeV","62GeV","200GeV_Run14"};
  float const mEnergyValue[NumBeamEnergy] = {27.0,39.0,62.4,200.0};
  int cent_start[4] = {3,2,2,3}; // not the defination in StRefMultCorr
  int cent_stop[4]  = {6,5,5,6};

  TFile *File_Input[4];
  TGraphAsymmErrors *g_rhoCent_1st_stat[4];
  TGraphAsymmErrors *g_rhoCent_2nd_stat[4];

  for(int i_energy = 0; i_energy < 4; ++i_energy)
  {
    string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/CentDependence/rhoCent_%s_LXaxis.root",mBeamEnergy[i_energy].c_str());
    File_Input[i_energy] = TFile::Open(inputfile.c_str());
    string GrapName; 
    GrapName = Form("g_rhoCent_%s_1st_stat",mBeamEnergy[i_energy].c_str());
    g_rhoCent_1st_stat[i_energy] = (TGraphAsymmErrors*)File_Input[i_energy]->Get(GrapName.c_str());
    GrapName = Form("g_rhoCent_%s_2nd_stat",mBeamEnergy[i_energy].c_str());
    g_rhoCent_2nd_stat[i_energy] = (TGraphAsymmErrors*)File_Input[i_energy]->Get(GrapName.c_str());
  }

  TGraphAsymmErrors *g_rhoCentInte_1st = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rhoCentInte_2nd = new TGraphAsymmErrors();

  // double rhoInte_1st[4] = {0,0,0,0};
  // double errInte_1st[4] = {0,0,0,0};
  // double rhoInte_2nd[4] = {0,0,0,0};
  // double errInte_2nd[4] = {0,0,0,0};
  for(int i_energy = 0; i_energy < 4; ++i_energy)
  {
    double rhoSum_1st = 0.0; 
    double errSum_1st = 0.0;
    double wgtSum_1st = 0.0;
    double rhoSum_2nd = 0.0; 
    double errSum_2nd = 0.0;
    double wgtSum_2nd = 0.0;

    for(int i_cent = 3; i_cent <= 6; ++i_cent)
    { // integrated from 20-60%
      double cent_1st, rho_1st, err_1st;
      g_rhoCent_1st_stat[i_energy]->GetPoint(i_cent,cent_1st,rho_1st);
      err_1st = g_rhoCent_1st_stat[i_energy]->GetErrorYhigh(i_cent);
      double weight_1st = 1.0/(err_1st*err_1st);
      rhoSum_1st += rho_1st*weight_1st;
      errSum_1st += err_1st*err_1st*weight_1st*weight_1st;
      wgtSum_1st += weight_1st;
      cout << "Energy: " << mBeamEnergy[i_energy].c_str() << ", cent = " << cent_1st << ", rho_1st = " << rho_1st << " +/- " << err_1st << endl;
    }
    double rhoInte_1st = rhoSum_1st/wgtSum_1st;
    double errInte_1st = sqrt(errSum_1st)/wgtSum_1st;
    g_rhoCentInte_1st->SetPoint(i_energy,mEnergyValue[i_energy],rhoInte_1st);
    g_rhoCentInte_1st->SetPointError(i_energy,0.0,0.0,errInte_1st,errInte_1st);
    cout << "Energy: " << mBeamEnergy[i_energy].c_str() << ", rho_1st = " << rhoInte_1st << " +/- " << errInte_1st << endl;
  }

  for(int i_energy = 0; i_energy < 4; ++i_energy)
  {
    double rhoSum_1st = 0.0; 
    double errSum_1st = 0.0;
    double wgtSum_1st = 0.0;
    double rhoSum_2nd = 0.0; 
    double errSum_2nd = 0.0;
    double wgtSum_2nd = 0.0;

    for(int i_cent = cent_start[i_energy]; i_cent <= cent_stop[i_energy]; ++i_cent)
    { // integrated from 20-60%
      double cent_2nd, rho_2nd, err_2nd;
      g_rhoCent_2nd_stat[i_energy]->GetPoint(i_cent,cent_2nd,rho_2nd);
      err_2nd = g_rhoCent_2nd_stat[i_energy]->GetErrorYhigh(i_cent);
      double weight_2nd = 1.0/(err_2nd*err_2nd);
      rhoSum_2nd += rho_2nd*weight_2nd;
      errSum_2nd += err_2nd*err_2nd*weight_2nd*weight_2nd;
      wgtSum_2nd += weight_2nd;
      cout << "Energy: " << mBeamEnergy[i_energy].c_str() << ", cent = " << cent_2nd << ", rho_2nd = " << rho_2nd << " +/- " << err_2nd << endl;
    }
    double rhoInte_2nd = rhoSum_2nd/wgtSum_2nd;
    double errInte_2nd = sqrt(errSum_2nd)/wgtSum_2nd;
    g_rhoCentInte_2nd->SetPoint(i_energy,mEnergyValue[i_energy],rhoInte_2nd);
    g_rhoCentInte_2nd->SetPointError(i_energy,0.0,0.0,errInte_2nd,errInte_2nd);
    cout << "Energy: " << mBeamEnergy[i_energy].c_str() << ", rho_2nd = " << rhoInte_2nd << " +/- " << errInte_2nd << endl;
  }

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/CentDependence/rhoCentEnergyDependence_LXaxis.root";
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  g_rhoCentInte_1st->SetName("g_rhoCentInte_1st");
  g_rhoCentInte_1st->Write();
  g_rhoCentInte_2nd->SetName("g_rhoCentInte_2nd");
  g_rhoCentInte_2nd->Write();
  File_OutPut->Close();
}
