#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include <TRandom3.h>
#include <map>
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/draw.h"
//#include "../../Utility/type.h"

typedef std::map<std::string,TH1D*> TH1DMap;

using namespace std;

void comTpcTrackEfficiency_Mother(int energy = 4, int mpid = 0, int pid = 0, int year = 0, int EP = 0)
{
  gStyle->SetOptDate(0);
  gRandom->SetSeed();

  std::string eventplane = "";
  if(EP == 1) eventplane = "_EP";
 
  int NCentralityMax = 9;
  int NEtaMax = 10;
  int NPhiMax = 12;

  string innput_0080 = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/TPC/Eff_%s_%s_Rc5Rc2.root",vmsa::mPID[mpid].c_str(),vmsa::mPID[mpid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  //string innput_0080 = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_pr_2060_root5.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str());
  // string innput_0080 = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_pr.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str());
  TFile *File_InPut_0080 = TFile::Open(innput_0080.c_str());

  TH1DMap h_mEfficiency_0080;
  File_InPut_0080->cd();
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    string HistNameEff;

    HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEfficiency_0080[HistNameEff] = (TH1D*)File_InPut_0080->Get(HistNameEff.c_str());
    cout << "read in => " << HistNameEff.c_str() << endl;

    HistNameEff = Form("h_mEffEta_Cent_%d",i_cent);
    h_mEfficiency_0080[HistNameEff] = (TH1D*)File_InPut_0080->Get(HistNameEff.c_str());
    cout << "read in => " << HistNameEff.c_str() << endl;

    HistNameEff = Form("h_mEffPhi_Cent_%d",i_cent);
    h_mEfficiency_0080[HistNameEff] = (TH1D*)File_InPut_0080->Get(HistNameEff.c_str());
    cout << "read in => " << HistNameEff.c_str() << endl;
    for(int i_eta = 0; i_eta < 2.0*vmsa::BinEta; ++i_eta)
    {
      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
	HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mEfficiency_0080[HistNameEff] = (TH1D*)File_InPut_0080->Get(HistNameEff.c_str());
	cout << "read in => " << HistNameEff.c_str() << endl;
      }
    }
  }

  //string innput_2060 = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_pr_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str());
  //TFile *File_InPut_2060 = TFile::Open(innput_2060.c_str());

  //TH1DMap h_mEfficiency_2060;
  //File_InPut_2060->cd();
  //for(int i_cent = 0; i_cent < 10; ++i_cent)
  //{
  //  string HistNameEff;

  //  HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
  //  h_mEfficiency_2060[HistNameEff] = (TH1D*)File_InPut_2060->Get(HistNameEff.c_str());
  //  cout << "read in => " << HistNameEff.c_str() << endl;

  //  HistNameEff = Form("h_mEffEta_Cent_%d",i_cent);
  //  h_mEfficiency_2060[HistNameEff] = (TH1D*)File_InPut_2060->Get(HistNameEff.c_str());
  //  cout << "read in => " << HistNameEff.c_str() << endl;

  //  HistNameEff = Form("h_mEffPhi_Cent_%d",i_cent);
  //  h_mEfficiency_2060[HistNameEff] = (TH1D*)File_InPut_2060->Get(HistNameEff.c_str());
  //  cout << "read in => " << HistNameEff.c_str() << endl;
  //  for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
  //  {
  //    for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
  //    {
  //      HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
  //      h_mEfficiency_2060[HistNameEff] = (TH1D*)File_InPut_2060->Get(HistNameEff.c_str());
  //      cout << "read in => " << HistNameEff.c_str() << endl;
  //    }
  //  }
  //}

  // int const cent_bin = floor(NCentralityMax * gRandom->Rndm());
  //int const cent_bin = 9;
  //int const eta_bin = floor(NEtaMax * gRandom->Rndm());
  //cout << "cent_bin = " << cent_bin << ", eta_bin = " << eta_bin << endl;

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,1600,1200);
  c_play->Divide(4,3);
 
  double etavalue[21] = {-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0};
  int centvalue[10] = {80,70,60,50,40,30,20,10,5,0};

  string outputname = Form("figures/%s/%s/TPC_Efficiency_%s_%s_Rc5Rc2.pdf",vmsa::mPID[mpid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[mpid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string output_start = Form("%s[",outputname.c_str());

  c_play->Print(output_start.c_str()); 

  for(int icent = 0; icent < 10; icent++)
  {
    for(int ieta = 0; ieta < 2.0*vmsa::BinEta; ieta++) 
    {
      for(int i_pad = 0; i_pad < 12; ++i_pad)
      {
        c_play->cd(i_pad+1)->SetLeftMargin(0.15);
        c_play->cd(i_pad+1)->SetBottomMargin(0.15);
        c_play->cd(i_pad+1)->SetGrid(0,0);
        c_play->cd(i_pad+1)->SetTicks(1,1);
        string HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",icent,ieta,i_pad);

        //h_mEfficiency_2060[HistNameEff]->SetMarkerStyle(24);
        //h_mEfficiency_2060[HistNameEff]->SetMarkerSize(1.2);
        //h_mEfficiency_2060[HistNameEff]->SetMarkerColor(2);
        //h_mEfficiency_2060[HistNameEff]->SetLineColor(2);
        //h_mEfficiency_2060[HistNameEff]->Draw("pE");
        cout << Form("Cent(%d-%d),%.1f<#eta<%.1f,%d#pi/6<#phi<%d#pi/6",centvalue[icent+1],centvalue[icent],etavalue[ieta],etavalue[ieta+1],i_pad-6,i_pad-5) << endl;
        if(icent < 9)  h_mEfficiency_0080[HistNameEff]->SetTitle(Form("Cent(%d-%d), %.1f<#eta<%.1f, %d#pi/6<#phi<%d#pi/6",centvalue[icent+1],centvalue[icent],etavalue[ieta],etavalue[ieta+1],i_pad,i_pad+1));
        if(icent == 9) h_mEfficiency_0080[HistNameEff]->SetTitle(Form("Cent(%d-%d), %.1f<#eta<%.1f, %d#pi/6<#phi<%d#pi/6",int(20),int(60),                    etavalue[ieta],etavalue[ieta+1],i_pad,i_pad+1));
        h_mEfficiency_0080[HistNameEff]->SetMarkerStyle(20);
        h_mEfficiency_0080[HistNameEff]->SetMarkerSize(1.0);
        h_mEfficiency_0080[HistNameEff]->SetMarkerColor(1);
        h_mEfficiency_0080[HistNameEff]->SetLineColor(1);
        h_mEfficiency_0080[HistNameEff]->Draw("pE");
      } 
      c_play->Update();
      c_play->Print(outputname.c_str());
    }
  }
  string output_stop = Form("%s]",outputname.c_str());
  c_play->Print(output_stop.c_str()); // close pdf file  
}
