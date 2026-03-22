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
float pos_y[8] = {0.56,0.56,0.45,0.45,0.40,0.40,0.35,0.35};

int centval[10] = {80,70,60,50,40,30,20,10,5,0};

void plotMcPhiEff_Centrality(int energy = 4, int pid = 0, int etamode = 0, int order = 1)
{
  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(1);
  //string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/Eff_%s_SingleParticle_2060_fitto%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),setting.c_str());
  //string inputfile = "../../RcPhiEffCorr/Eff_19GeV_SingleParticle_noToF_Mode2_EtaMode0.root";
  //string inputfile = Form("../../RcPhiEffCorr/Eff_19GeV_SingleParticle_noToF_Mode2_EtaMode%d_PtBins0_1.root",etamode);
  string inputfile;
  if(energy == 3) inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu14GeV_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d/Eff_14GeV_SingleParticle_noToF_Mode1_EtaMode0.root",order);
  if(energy == 4 && order == 2) inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d_FixedRes/Eff_19GeV_SingleParticle_noToF_Mode1_EtaMode0.root",order);
  if(energy == 4 && order == 1) inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d/Eff_19GeV_SingleParticle_noToF_Mode1_EtaMode0.root",order);
  //string inputfile = Form("/star/u/gwilks3/Workspace/FileTransfers/Eff_%s_SingleKaon_second.root",vmsa::mBeamEnergy[energy].c_str());
  // string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/Efficiency/Eff_%s_SingleKaon.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  int mpt_rebin;
  int mpt_rebin_last;
  int mpt_rebin_first; 

  TH1D *h_mEff[10][10];
  for(int i_cent = 0; i_cent < 9; i_cent++)
  {
    for(int i_pt = 0/*vmsa::pt_rebin_first_y[energy]*/; i_pt <= 2/*vmsa::pt_rebin_last_y[energy]*/; ++i_pt)
    {
        string HistName = Form("h_mEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEff[i_cent][i_pt] = (TH1D*)File_InPut->Get(HistName.c_str());
    }
  }

  TCanvas *c_eff = new TCanvas("c_eff","c_eff",10,10,1000,1000);
  c_eff->Divide(3,3);
  TH1D *h_frame = new TH1D("h_frame","h_frame",100,0.0,1.0);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10);
    h_frame->SetBinError(i_bin+1,1);
  }
  string legEnergy = Form("AuAu %s 20%%-60%%",vmsa::mBeamEnergy[energy].c_str());
  //h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetTitle("|cos(#theta*)|");
  h_frame->GetXaxis()->CenterTitle();
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetNdivisions(505);

  h_frame->GetYaxis()->SetTitle("efficiency x acceptance (RC/MC)");
  h_frame->GetYaxis()->SetTitleSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->SetNdivisions(505);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetRangeUser(1e-3,0.6);
  if(pid == 2) h_frame->GetYaxis()->SetRangeUser(0.0,1.0);

  double intVal[10][20][20];
  double intErr[10][20][20];
  double slopeVal[10][20][20];
  double slopeErr[10][20][20];

  string outputname = Form("figures/AnalysisNote/Efficiency/c_effAuAu%s_com_%s_centrality_etamode%d_order%d.pdf",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode,order);
  string output_start = Form("%s[",outputname.c_str());
  c_eff->Print(output_start.c_str());

  for(int i_cent = 0; i_cent < 9; i_cent++)
  {
    c_eff->cd(i_cent+1);     
    //c_eff->cd(i_cent+1)->SetTopMargin(0.15);
    c_eff->cd(i_cent+1)->SetLeftMargin(0.15);
    c_eff->cd(i_cent+1)->SetBottomMargin(0.15);
    c_eff->cd(i_cent+1)->SetGrid(0,0);
    c_eff->cd(i_cent+1)->SetTicks(1,1);
    //c_eff->cd(i_y+1-2)->SetLogy();
    h_frame->SetTitle(Form("Centrality %d-%d \%",centval[i_cent+1],centval[i_cent]));
    h_frame->GetXaxis()->SetTitle(Form("|cos(#theta*)|"));
    h_frame->DrawCopy("pE");
    for(int i_pt = 0/*vmsa::pt_rebin_first_y[energy]*/; i_pt <=2/* vmsa::pt_rebin_last_y[energy]*/; ++i_pt)
    {
      h_mEff[i_cent][i_pt]->SetMarkerStyle(vmsa::Style[i_pt]);
      h_mEff[i_cent][i_pt]->SetMarkerColor(vmsa::Color[i_pt]);
      h_mEff[i_cent][i_pt]->SetMarkerSize(1.4);
      //h_mEff[i_cent][i_pt][i_y]->SetTitle(Form("Cent %d Rapidity %d",i_cent,i_y));
      h_mEff[i_cent][i_pt]->Draw("pEX0 same");
      //TF1 *f_poly = new TF1("f_poly","pol1",0.0,1.0);
      //h_mEff[i_cent][i_pt][i_y]->Fit(f_poly,"QN");

      //intVal[i_cent][i_pt][i_y] = f_poly->GetParameter(0);
      //intErr[i_cent][i_pt][i_y] = f_poly->GetParError(0);
      //slopeVal[i_cent][i_pt][i_y] = f_poly->GetParameter(1);
      //slopeErr[i_cent][i_pt][i_y] = f_poly->GetParError(1);

      //f_poly->SetLineColor(vmsa::Color[i_pt]);
      //f_poly->SetLineStyle(2);
      //f_poly->SetLineWidth(2);
      //f_poly->Draw("l same");
      string pt_range = Form("p_{T} = %1.1f-%1.1f GeV/c",vmsa::pt_low_cent[energy][i_pt],vmsa::pt_up_cent[energy][i_pt]);
      string cent_range = Form("Centrality %d-%d",centval[i_cent+1],centval[i_cent]);
      Draw_TGAE_Point_new_Symbol(pos_x[i_pt],pos_y[i_pt]-0.05,0.0,0.0,0.0,0.0,vmsa::Style[i_pt],vmsa::Color[i_pt],1.2);
      plotTopLegend((char*)pt_range.c_str(),pos_x[i_pt]+0.03,pos_y[i_pt]-0.055,0.04,1,0.0,42,0,1);
      plotTopLegend((char*)cent_range.c_str(),0.3,0.05,0.04,1,0.0,42,0,1);
    }
    //plotTopLegend((char*)legEnergy.c_str(),0.17,0.5,0.05,1,0.0,42,0,1);

    //for(int i_pt = 0/*(vmsa::pt_rebin_first_y[energy]*/; i_pt <= 1/*vmsa::pt_rebin_last_y[energy]*/; ++i_pt)
    //{
    //  cout << "cent = " << i_cent << ",  pT bin = : " << i_pt << ",   y = " << i_y << endl;
    //  cout << "y-int = " << intVal[i_cent][i_pt][i_y] << " +/- " << intErr[i_cent][i_pt][i_y] << endl;
    //  cout << "slope = " << slopeVal[i_cent][i_pt][i_y] << " +/- " << slopeErr[i_cent][i_pt][i_y] << endl;
    //}
  }

  c_eff->Update();
  c_eff->Print(outputname.c_str());

  cout << "Finished Loop" << endl;
  string output_stop = Form("%s]",outputname.c_str());
  c_eff->Print(output_stop.c_str()); // close pdf file
  cout << "Closed the File" << endl;
  //string FigName = Form("figures/AnalysisNote/Efficiency/c_effAuAu%s_com_%s_cent.pdf",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  //c_eff->SaveAs(FigName.c_str());
  // c_eff->SaveAs("../figures/effPt.png");
}
