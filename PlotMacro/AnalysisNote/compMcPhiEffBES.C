#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/draw.h"

using namespace std;

float pos_x[8] = {0.05,0.54,0.05,0.54,0.05,0.54,0.05,0.54};
float pos_y[8] = {1.40,1.40,1.3,1.3,1.2,1.2,1.1,1.1};

void compMcPhiEffBES(int energy = 4, int cent = 9)
{
  gStyle->SetOptDate(0);
  string inputfile = Form("/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/RcPhiEffCorr/Eff_%s_SingleParticle_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputfileBESI = Form("/star/u/gwilks3/Workspace/FileTransfers/Eff_%s_SingleKaon_second.root",vmsa::mBeamEnergy[energy].c_str());
  // string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/Efficiency/Eff_%s_SingleKaon.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1D *h_mEff[vmsa::pt_rebin];
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  {
    string HistName = Form("h_mEffCos_Cent_%d_Pt_%d",cent,i_pt);
    h_mEff[i_pt] = (TH1D*)File_InPut->Get(HistName.c_str());
  }

  TFile *File_InPutBESI = TFile::Open(inputfileBESI.c_str());
  TH1D *h_mEffBESI[vmsa::pt_rebin];
  TH1D *h_mEffRatio[vmsa::pt_rebin];
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  {
    string HistName = Form("h_mEffCos_Cent_%d_Pt_%d",cent,i_pt);
    h_mEffBESI[i_pt] = (TH1D*)File_InPutBESI->Get(HistName.c_str());
    h_mEffRatio[i_pt] = (TH1D*)h_mEff[i_pt]->Clone();
    h_mEffRatio[i_pt]->Divide(h_mEff[i_pt],h_mEffBESI[i_pt],1,1,"B");
  }

  TCanvas *c_eff = new TCanvas("c_eff","c_eff",10,10,800,800);
  c_eff->SetLeftMargin(0.15);
  c_eff->SetBottomMargin(0.15);
  c_eff->SetGrid(0,0);
  c_eff->SetTicks(1,1);
  TH1D *h_frame = new TH1D("h_frame","h_frame",100,0.0,1.0);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10);
    h_frame->SetBinError(i_bin+1,1);
  }
  string legEnergy = Form("AuAu %s 20%%-60%%",vmsa::mBeamEnergy[energy].c_str());
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetTitle("cos(#theta*)");
  h_frame->GetXaxis()->CenterTitle();
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetNdivisions(505);

  h_frame->GetYaxis()->SetTitle("Efficiency Ratio BES-II/BES-I");
  h_frame->GetYaxis()->SetTitleSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->SetNdivisions(505);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetRangeUser(1.0,2.2);
  h_frame->Draw("pE");

  double intVal[6];
  double intErr[6];
  double slopeVal[6];
  double slopeErr[6];

  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  {
    h_mEffRatio[i_pt]->SetMarkerStyle(vmsa::Style[i_pt]);
    h_mEffRatio[i_pt]->SetMarkerColor(vmsa::Color[i_pt]);
    h_mEffRatio[i_pt]->SetMarkerSize(1.4);
    h_mEffRatio[i_pt]->Draw("pEX0 same");
    TF1 *f_poly = new TF1("f_poly","pol1",0.0,1.0);
    h_mEffRatio[i_pt]->Fit(f_poly,"N");

    intVal[i_pt] = f_poly->GetParameter(0);
    intErr[i_pt] = f_poly->GetParError(0);
    slopeVal[i_pt] = f_poly->GetParameter(1);
    slopeErr[i_pt] = f_poly->GetParError(1);

    f_poly->SetLineColor(vmsa::Color[i_pt]);
    f_poly->SetLineStyle(2);
    f_poly->SetLineWidth(2);
    f_poly->Draw("l same");
    string pt_range = Form("p_{T} = %1.1f-%1.1f GeV/c",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    Draw_TGAE_Point_new_Symbol(pos_x[i_pt],pos_y[i_pt]-0.05,0.0,0.0,0.0,0.0,vmsa::Style[i_pt],vmsa::Color[i_pt],1.2);
    plotTopLegend((char*)pt_range.c_str(),pos_x[i_pt]+0.03,pos_y[i_pt]-0.055,0.04,1,0.0,42,0,1);
  }
  plotTopLegend((char*)legEnergy.c_str(),0.17,1.025,0.05,1,0.0,42,0,1);

  string ptbin[6] = {"0.4-0.8 GeV/c","0.8-1.2 GeV/c","1.2-1.8 GeV/c","1.8-2.4 GeV/c","2.4-3.0 GeV/c","3.0-4.2 GeV/c"};
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  {
    string consistent = "Yes";
    if(slopeVal[i_pt] > 0.0) {if(slopeVal[i_pt]-slopeErr[i_pt] > 0) consistent = "No, BES-II is higher";}
    if(slopeVal[i_pt] <= 0.0) {if(slopeVal[i_pt]+slopeErr[i_pt] < 0) consistent = "No, BES-II is lower";}
    cout << "pT = " << ptbin[i_pt] << std::fixed << std::setprecision(4) << "    y-int = " << intVal[i_pt] << " +/- " << intErr[i_pt] << "    slope = ";
    if(slopeVal[i_pt] > 0) cout << " ";
    cout << slopeVal[i_pt] << " +/- " << slopeErr[i_pt] << "    Consistent? " << consistent << endl;
  }

  string FigName = Form("figures/AnalysisNote/Efficiency/c_effAuAu%s_comp_BES.pdf",vmsa::mBeamEnergy[energy].c_str());
  c_eff->SaveAs(FigName.c_str());
  // c_eff->SaveAs("../figures/effPt.png");
}
