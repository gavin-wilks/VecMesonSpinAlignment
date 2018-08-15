#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

void plotQA_ToFCom()
{
  TString inputfile = "/global/homes/x/xusun/AuAu200GeV/SpinAlignment/ToFMatch/Eff_200GeV_ToFMatch.root";
  TFile *File_InPut = TFile::Open(inputfile.Data());
  TH1D *h_mEffKPlus = (TH1D*)File_InPut->Get("h_mEfficiency_Kplus_Cent_9_Eta_0_Phi_0")->Clone();
  TH1D *h_mEffKMinus = (TH1D*)File_InPut->Get("h_mEfficiency_Kminus_Cent_9_Eta_0_Phi_0")->Clone();
  TH1D *h_mEffPiPlus = (TH1D*)File_InPut->Get("h_mEfficiency_Piplus_Cent_9_Eta_0_Phi_0")->Clone();
  TH1D *h_mEffPiMinus = (TH1D*)File_InPut->Get("h_mEfficiency_Piminus_Cent_9_Eta_0_Phi_0")->Clone();

  TString inputShuai = "/global/homes/x/xusun/AuAu200GeV/SpinAlignment/ToFMatch/allTofEffForXu.root";
  TFile *File_Shuai = TFile::Open(inputShuai.Data());
  TH1D *h_mEffKPlus_Shuai = (TH1D*)File_Shuai->Get("hKPlusEff1D")->Clone();
  TH1D *h_mEffKMinus_Shuai = (TH1D*)File_Shuai->Get("hKMinusEff1D")->Clone();
  TH1D *h_mEffPiPlus_Shuai = (TH1D*)File_Shuai->Get("hPiPlusEff1D")->Clone();
  TH1D *h_mEffPiMinus_Shuai = (TH1D*)File_Shuai->Get("hPiMinusEff1D")->Clone();

  TCanvas *c_Efficiency = new TCanvas("c_Efficiency","c_Efficiency",10,10,1200,1200);
  c_Efficiency->Divide(2,2);
  for(int i_pad = 0; i_pad < 4; ++i_pad)
  {
    c_Efficiency->cd(i_pad+1);
    c_Efficiency->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Efficiency->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Efficiency->cd(i_pad+1)->SetTicks(1,1);
    c_Efficiency->cd(i_pad+1)->SetGrid(0,0);
  }

  c_Efficiency->cd(1); // Kplus
  h_mEffKPlus_Shuai->Draw("pE");
  h_mEffKPlus->SetMarkerStyle(24);
  h_mEffKPlus->SetMarkerSize(1.2);
  h_mEffKPlus->SetMarkerColor(4);
  h_mEffKPlus->Draw("pE same");

  c_Efficiency->cd(2); // Kplus
  h_mEffKMinus_Shuai->Draw("pE");
  h_mEffKMinus->SetMarkerStyle(24);
  h_mEffKMinus->SetMarkerSize(1.2);
  h_mEffKMinus->SetMarkerColor(4);
  h_mEffKMinus->Draw("pE same");

  c_Efficiency->cd(3); // Kplus
  h_mEffPiPlus_Shuai->Draw("pE");
  h_mEffPiPlus->SetMarkerStyle(24);
  h_mEffPiPlus->SetMarkerSize(1.2);
  h_mEffPiPlus->SetMarkerColor(4);
  h_mEffPiPlus->Draw("pE same");

  c_Efficiency->cd(4); // Kplus
  h_mEffPiMinus_Shuai->Draw("pE");
  h_mEffPiMinus->SetMarkerStyle(24);
  h_mEffPiMinus->SetMarkerSize(1.2);
  h_mEffPiMinus->SetMarkerColor(4);
  h_mEffPiMinus->Draw("pE same");

  c_Efficiency->SaveAs("../../figures/c_ToFMatchEff_ComToShuai.eps");
  c_Efficiency->SaveAs("../../figures/c_ToFMatchEff_ComToShuai.pdf");
  c_Efficiency->SaveAs("../../figures/c_ToFMatchEff_ComToShuai.png");
}
