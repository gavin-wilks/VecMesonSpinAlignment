#include <string>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <map>

typedef std::map<std::string,TH1D*> TH1DMap;
typedef std::map<std::string,TF1*> TF1Map;

using namespace std;

// tof matching efficiency
double tof_Kaon(double* x, double* par)
{
  return par[0]*(1.0 / (pow(x[0] - par[1], 2) + par[2]) - par[4] / (exp(x[0] - par[3]) + par[5]) + par[6]);
}

void plotQA_TofFit(int mEnergy = 6, int mPID = 0)
{
  string const mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  string const mParType[2] = {"Kplus","Kminus"};
  string const mParTex[2] = {"K^{+}","K^{-}"};
  string const Centrality[10] = {"70%-80%","60%-70%","50%-60%","40%-50%","30%-40%","20%-30%","10%-20%","5%-10%","0%-5%","20%-60%"}; // Centrality bin
  int const mColor[12] = {11,11,2,2,4,4,7,7,8,8,38,38};
  int const mStyle[12] = {30,24,30,24,30,24,30,24,30,24,30,24};
  int const mLineStyle[12] = {3,2,3,2,3,2,3,2,3,2,3,2};

  string inputEff = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/ToFMatch/Eff_%s_ToFMatch_2060.root",mBeamEnergy[mEnergy].c_str(),mBeamEnergy[mEnergy].c_str());
  cout << "open input file: " << inputEff.c_str() << endl;

  TFile *File_TofEff = TFile::Open(inputEff.c_str());

  TH1DMap h_mEfficiency;
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    string HistName;
    HistName = Form("h_mEfficiency_%s_Cent_%d",mParType[mPID].c_str(),i_cent);
    h_mEfficiency[HistName] = (TH1D*)File_TofEff->Get(HistName.c_str());
    cout << "read in => " << HistName.c_str() << endl;
    for(int i_eta = 0; i_eta < 12; ++i_eta)
    {
      HistName = Form("h_mEfficiency_%s_Cent_%d_Eta_%d",mParType[mPID].c_str(),i_cent,i_eta);
      h_mEfficiency[HistName] = (TH1D*)File_TofEff->Get(HistName.c_str());
      cout << "read in => " << HistName.c_str() << endl;
      for(int i_phi = 0; i_phi < 12; ++i_phi)
      {
	HistName = Form("h_mEfficiency_%s_Cent_%d_Eta_%d_Phi_%d",mParType[mPID].c_str(),i_cent,i_eta,i_phi);
	h_mEfficiency[HistName] = (TH1D*)File_TofEff->Get(HistName.c_str());
	cout << "read in => " << HistName.c_str() << endl;
      }
    }
  }

  string inputfirst = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/ToFMatch/FitPar_AuAu%s_%s_first_2060.root",mBeamEnergy[mEnergy].c_str(),mBeamEnergy[mEnergy].c_str(),mParType[mPID].c_str());
  TFile *File_First = TFile::Open(inputfirst.c_str());
  cout << "OPEN ToF Matching Efficiency Fit File for First: " << inputfirst.c_str() << endl;

  string inputsecond = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/ToFMatch/FitPar_AuAu%s_%s_second_2060.root",mBeamEnergy[mEnergy].c_str(),mBeamEnergy[mEnergy].c_str(),mParType[mPID].c_str());
  TFile *File_Second = TFile::Open(inputsecond.c_str());
  cout << "OPEN ToF Matching Efficiency Fit File for First: " << inputsecond.c_str() << endl;

  TF1Map f_TofFirst;
  TF1Map f_TofSecond;

  for(int i_eta = 0; i_eta < 12; ++i_eta) // eta
  {
    TH1D *h_first = NULL;
    string KEY;
    KEY = Form("h_mFitParameters_%s_Cent_9_Eta_%d",mParType[mPID].c_str(),i_eta);
    h_first = (TH1D*)File_First->Get(KEY.c_str());
    KEY = Form("f_mToFMatch_%s_Cent_9_Eta_%d",mParType[mPID].c_str(),i_eta);
    f_TofFirst[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.2,10,7);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_TofFirst[KEY]->FixParameter(i_par,h_first->GetBinContent(i_par+1));
    }

    TH1D *h_second = NULL;
    KEY = Form("h_mFitParameters_%s_Cent_9_Eta_%d",mParType[mPID].c_str(),i_eta);
    h_second = (TH1D*)File_Second->Get(KEY.c_str());
    KEY = Form("f_mToFMatch_%s_Cent_9_Eta_%d",mParType[mPID].c_str(),i_eta);
    f_TofSecond[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.2,10,7);
    for(int i_par = 0; i_par < 7; ++i_par)
    {
      f_TofSecond[KEY]->FixParameter(i_par,h_second->GetBinContent(i_par+1));
    }
  }

  for(int i_eta = 0; i_eta < 12; ++i_eta) // eta & phi
  {
    for(int i_phi = 0; i_phi < 12; ++i_phi)
    {
      TH1D *h_first = NULL;
      string KEY;
      KEY = Form("h_mFitParameters_%s_Cent_9_Eta_%d_Phi_%d",mParType[mPID].c_str(),i_eta,i_phi);
      h_first = (TH1D*)File_First->Get(KEY.c_str());
      KEY = Form("f_mToFMatch_%s_Cent_9_Eta_%d_Phi_%d",mParType[mPID].c_str(),i_eta,i_phi);
      f_TofFirst[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.2,10,7);
      for(int i_par = 0; i_par < 7; ++i_par)
      {
	f_TofFirst[KEY]->FixParameter(i_par,h_first->GetBinContent(i_par+1));
      }

      TH1D *h_second = NULL;
      KEY = Form("h_mFitParameters_%s_Cent_9_Eta_%d_Phi_%d",mParType[mPID].c_str(),i_eta,i_phi);
      h_second = (TH1D*)File_Second->Get(KEY.c_str());
      KEY = Form("f_mToFMatch_%s_Cent_9_Eta_%d_Phi_%d",mParType[mPID].c_str(),i_eta,i_phi);
      f_TofSecond[KEY] = new TF1(KEY.c_str(),tof_Kaon,0.2,10,7);
      for(int i_par = 0; i_par < 7; ++i_par)
      {
	f_TofSecond[KEY]->FixParameter(i_par,h_second->GetBinContent(i_par+1));
      }
    }	
  }

  TCanvas *c_Efficiency = new TCanvas("c_Efficiency","c_Efficiency",10,10,1600,800);
  c_Efficiency->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_Efficiency->cd(i_pad+1)->SetLeftMargin(0.1);
    c_Efficiency->cd(i_pad+1)->SetBottomMargin(0.1);
    c_Efficiency->cd(i_pad+1)->SetGrid(0,0);
    c_Efficiency->cd(i_pad+1)->SetTicks(1,1);
  }

  int const eta_bin = 4;

  string title;
  c_Efficiency->cd(1);
  string HistName = Form("h_mEfficiency_%s_Cent_9_Eta_%d",mParType[mPID].c_str(),eta_bin);
  TLegend *leg = new TLegend(0.5,0.7,0.8,0.9);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  title = Form("%s fit to peak",mParTex[mPID].c_str());
  h_mEfficiency[HistName]->SetTitle(title.c_str());
  h_mEfficiency[HistName]->SetStats(0);
  h_mEfficiency[HistName]->SetMarkerStyle(29);
  h_mEfficiency[HistName]->SetMarkerSize(1.6);
  h_mEfficiency[HistName]->SetMarkerColor(2);
  h_mEfficiency[HistName]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mEfficiency[HistName]->GetXaxis()->CenterTitle();
  h_mEfficiency[HistName]->GetYaxis()->SetTitle("Efficiency");
  h_mEfficiency[HistName]->GetYaxis()->SetRangeUser(0.0,1.2);
  h_mEfficiency[HistName]->DrawCopy("pE");
  leg->AddEntry(h_mEfficiency[HistName],"eta & phi integrated","p");
  string KEY_first = Form("f_mToFMatch_%s_Cent_9_Eta_%d",mParType[mPID].c_str(),eta_bin);
  f_TofFirst[KEY_first]->SetLineColor(2);
  f_TofFirst[KEY_first]->SetLineWidth(4);
  f_TofFirst[KEY_first]->SetLineStyle(2);
  f_TofFirst[KEY_first]->SetNpx(1000);
  f_TofFirst[KEY_first]->SetRange(0.2,8.0);
  f_TofFirst[KEY_first]->Draw("l same");

  for(int i_phi = 0; i_phi < 12; ++i_phi) 
  {
    HistName = Form("h_mEfficiency_%s_Cent_9_Eta_%d_Phi_%d",mParType[mPID].c_str(),eta_bin,i_phi);
    h_mEfficiency[HistName]->SetMarkerColor(mColor[i_phi]);
    h_mEfficiency[HistName]->SetMarkerStyle(mStyle[i_phi]);
    h_mEfficiency[HistName]->SetMarkerSize(1.2);
    h_mEfficiency[HistName]->DrawCopy("pE same");
    string Leg_phi = Form("%d bin",i_phi);
    leg->AddEntry(h_mEfficiency[HistName],Leg_phi.c_str(),"p");
    KEY_first = Form("f_mToFMatch_%s_Cent_9_Eta_%d_Phi_%d",mParType[mPID].c_str(),eta_bin,i_phi);
    f_TofFirst[KEY_first]->SetLineColor(mColor[i_phi]);
    f_TofFirst[KEY_first]->SetLineWidth(4);
    f_TofFirst[KEY_first]->SetLineStyle(mLineStyle[i_phi]);
    f_TofFirst[KEY_first]->SetNpx(1000);
    f_TofFirst[KEY_first]->SetRange(0.2,8.0);
    f_TofFirst[KEY_first]->Draw("l same");
  }
  HistName = Form("h_mEfficiency_%s_Cent_9_Eta_%d",mParType[mPID].c_str(),eta_bin);
  h_mEfficiency[HistName]->DrawCopy("pE same");
  KEY_first = Form("f_mToFMatch_%s_Cent_9_Eta_%d",mParType[mPID].c_str(),eta_bin);
  f_TofFirst[KEY_first]->Draw("l same");
  leg->Draw("same");

  c_Efficiency->cd(2);
  HistName = Form("h_mEfficiency_%s_Cent_9_Eta_%d",mParType[mPID].c_str(),eta_bin);
  title = Form("%s fit to plateau",mParTex[mPID].c_str());
  h_mEfficiency[HistName]->SetTitle(title.c_str());
  h_mEfficiency[HistName]->SetStats(0);
  h_mEfficiency[HistName]->DrawCopy("pE");
  string KEY_second = Form("f_mToFMatch_%s_Cent_9_Eta_%d",mParType[mPID].c_str(),eta_bin);
  f_TofSecond[KEY_second]->SetLineColor(2);
  f_TofSecond[KEY_second]->SetLineWidth(4);
  f_TofSecond[KEY_second]->SetLineStyle(2);
  f_TofSecond[KEY_second]->SetNpx(1000);
  f_TofSecond[KEY_second]->SetRange(0.2,8.0);
  f_TofSecond[KEY_second]->Draw("l same");

  for(int i_phi = 0; i_phi < 12; ++i_phi) 
  {
    HistName = Form("h_mEfficiency_%s_Cent_9_Eta_%d_Phi_%d",mParType[mPID].c_str(),eta_bin,i_phi);
    h_mEfficiency[HistName]->DrawCopy("pE same");
    KEY_second = Form("f_mToFMatch_%s_Cent_9_Eta_%d_Phi_%d",mParType[mPID].c_str(),eta_bin,i_phi);
    f_TofSecond[KEY_second]->SetLineColor(mColor[i_phi]);
    f_TofSecond[KEY_second]->SetLineWidth(4);
    f_TofSecond[KEY_second]->SetLineStyle(mLineStyle[i_phi]);
    f_TofSecond[KEY_second]->SetNpx(1000);
    f_TofSecond[KEY_second]->SetRange(0.2,8.0);
    f_TofSecond[KEY_second]->Draw("l same");
  }
  HistName = Form("h_mEfficiency_%s_Cent_9_Eta_%d",mParType[mPID].c_str(),eta_bin);
  h_mEfficiency[HistName]->DrawCopy("pE same");
  KEY_second = Form("f_mToFMatch_%s_Cent_9_Eta_%d",mParType[mPID].c_str(),eta_bin);
  f_TofSecond[KEY_second]->Draw("l same");

  string FigName = Form("/Users/xusun/WorkSpace/Papers/VecMesonSpinAlignment/figures/Efficiency/ToF/c_TofEffDifFitEta%d_%s%s.eps",eta_bin,mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  c_Efficiency->SaveAs(FigName.c_str());
}
