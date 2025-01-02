#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include <TH3D.h>
// #include "../../Utility/functions.h"
#include "../../Utility/draw.h"
#include "../../Utility/type.h"
#include "../../Utility/StSpinAlignmentCons.h"

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

void calTpcTrackEfficiency(int energy = 6, int pid = 0, int year = 0)
{
  // string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_pr.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str());
  string inputfile = Form("/star/data01/pwg/sunxuhit//AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_pr.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  
  TH1D *h_FrameEta = (TH1D*)File_InPut->Get("h_FrameEta");
  TH1D *h_FramePhi = (TH1D*)File_InPut->Get("h_FramePhi");

  TH3D *h_mMcTracks[10];
  TH3D *h_mRcTracks[10];
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    std::string HistName;

    HistName = Form("h_mMcTracks_%d",i_cent);
    h_mMcTracks[i_cent] = (TH3D*)File_InPut->Get(HistName.c_str())->Clone();
    h_mMcTracks[i_cent]->Sumw2();

    HistName = Form("h_mRcTracks_%d",i_cent);
    h_mRcTracks[i_cent] = (TH3D*)File_InPut->Get(HistName.c_str())->Clone();
    h_mRcTracks[i_cent]->Sumw2();
  }

  for(int i_cent = 2; i_cent <= 5; ++i_cent) // make centrality 9 to be 20-60%
  {
    if(i_cent == 2) 
    {
      h_mMcTracks[9] = (TH3D*)h_mMcTracks[i_cent]->Clone("h_mMcTracks_9");
      h_mMcTracks[9]->Sumw2();
      h_mRcTracks[9] = (TH3D*)h_mRcTracks[i_cent]->Clone("h_mRcTracks_9");
      h_mRcTracks[9]->Sumw2();
    }
    else
    {
      h_mMcTracks[9]->Add(h_mMcTracks[i_cent]);
      h_mRcTracks[9]->Add(h_mRcTracks[i_cent]);
    }
  }

  TH1D *h_mMcEffPt[10]; // pt distritbution as a function of centrality
  TH1D *h_mRcEffPt[10];
  TH1D *h_mEffPt[10];

  TH1D *h_mMcEffEta[10]; // eta distritbution as a function of centrality
  TH1D *h_mRcEffEta[10];
  TH1D *h_mEffEta[10];

  TH1D *h_mMcEffPhi[10]; // phi distritbution as a function of centrality
  TH1D *h_mRcEffPhi[10];
  TH1D *h_mEffPhi[10];

  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    std::string HistNameMc, HistNameRc, HistNameEff;

    HistNameMc = Form("h_mMcEffPt_Cent_%d",i_cent);
    h_mMcEffPt[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("x")->Clone(HistNameMc.c_str());
    HistNameRc = Form("h_mRcEffPt_Cent_%d",i_cent);
    h_mRcEffPt[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("x")->Clone(HistNameRc.c_str());

    HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEffPt[i_cent] = (TH1D*)h_mMcEffPt[i_cent]->Clone(HistNameEff.c_str());
    h_mEffPt[i_cent]->SetTitle(HistNameEff.c_str());
    h_mEffPt[i_cent]->Reset();
    h_mEffPt[i_cent]->Divide(h_mRcEffPt[i_cent],h_mMcEffPt[i_cent],1,1,"B");
    cout << "Calculating => " << HistNameEff.c_str() << endl;

    HistNameMc = Form("h_mMcEffEta_Cent_%d",i_cent);
    h_mMcEffEta[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("y")->Clone(HistNameMc.c_str());
    HistNameRc = Form("h_mRcEffEta_Cent_%d",i_cent);
    h_mRcEffEta[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("y")->Clone(HistNameRc.c_str());

    HistNameEff = Form("h_mEffEta_Cent_%d",i_cent);
    h_mEffEta[i_cent] = (TH1D*)h_mMcEffEta[i_cent]->Clone(HistNameEff.c_str());
    h_mEffEta[i_cent]->SetTitle(HistNameEff.c_str());
    h_mEffEta[i_cent]->Reset();
    h_mEffEta[i_cent]->Divide(h_mRcEffEta[i_cent],h_mMcEffEta[i_cent],1,1,"B");
    cout << "Calculating => " << HistNameEff.c_str() << endl;

    HistNameMc = Form("h_mMcEffPhi_Cent_%d",i_cent);
    h_mMcEffPhi[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("z")->Clone(HistNameMc.c_str());
    HistNameRc = Form("h_mRcEffPhi_Cent_%d",i_cent);
    h_mRcEffPhi[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("z")->Clone(HistNameRc.c_str());

    HistNameEff = Form("h_mEffPhi_Cent_%d",i_cent);
    h_mEffPhi[i_cent] = (TH1D*)h_mMcEffPhi[i_cent]->Clone(HistNameEff.c_str());
    h_mEffPhi[i_cent]->SetTitle(HistNameEff.c_str());
    h_mEffPhi[i_cent]->Reset();
    h_mEffPhi[i_cent]->Divide(h_mRcEffPhi[i_cent],h_mMcEffPhi[i_cent],1,1,"B");
    cout << "Calculating => " << HistNameEff.c_str() << endl;
  }

  TH1DMap h_mMcEffPEP;  
  TH1DMap h_mRcEffPEP;
  TH1DMap h_mEfficiency; // efficiency as a fucntion of centrality, pt, eta and phi
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
    {
      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
	std::string HistNameMc = Form("h_mMcEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mMcEffPEP[HistNameMc] = (TH1D*)h_mMcTracks[i_cent]->ProjectionX(HistNameMc.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1);

	std::string HistNameRc = Form("h_mRcEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mRcEffPEP[HistNameRc] = (TH1D*)h_mRcTracks[i_cent]->ProjectionX(HistNameRc.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1);

	std::string HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mEfficiency[HistNameEff] = (TH1D*)h_mMcEffPEP[HistNameMc]->Clone(HistNameEff.c_str());
        h_mEfficiency[HistNameEff]->SetTitle(HistNameEff.c_str());
        h_mEfficiency[HistNameEff]->Reset();
        h_mEfficiency[HistNameEff]->Divide(h_mRcEffPEP[HistNameRc],h_mMcEffPEP[HistNameMc],1,1,"B");
	cout << "Calculating => " << HistNameEff.c_str() << endl;
      }
    }
  }


  // string outputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_pr_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str());
  string outputfile = Form("/star/data01/pwg/sunxuhit//AuAu%s/SpinAlignment/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_pr_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  // write histogram
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    h_mMcTracks[i_cent]->Write();
    h_mRcTracks[i_cent]->Write();
    h_mEffPt[i_cent]->Write();
    h_mEffEta[i_cent]->Write();
    h_mEffPhi[i_cent]->Write();

    for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
    {
      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
	std::string HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mEfficiency[HistNameEff]->Write();
      }
    }
  }
  h_FrameEta->Write();
  h_FramePhi->Write();
  File_OutPut->Close();

  File_InPut->cd();
  File_InPut->Close();
}
