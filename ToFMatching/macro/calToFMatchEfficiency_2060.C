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
#include "../../Utility/functions.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
//#include "../../Utility/type.h"
#include "../StRoot/StToFMatchMaker/StToFMatchCons.h"


typedef std::map<std::string,TH1D*> TH1DMap;

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

void calToFMatchEfficiency_2060(int energy = 4, int pid = 0, int EP = 5)
{

  string eventplane = "";
  if(EP == 1) eventplane = "_EP"; 
  if(EP == 2) eventplane = "_FinerPhiBins"; 
  if(EP == 3) eventplane = "_DCA2_nsigma0p5"; 
  if(EP == 4) eventplane = "_FinerEta_shiftedphi"; 
  if(EP == 5) eventplane = "_FinerEta_240bins"; 

  string inputfile = Form("../../data/file_%s_ToFMatching%s.root",vmsa::mBeamEnergy[energy].c_str(),eventplane.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH1D *h_FrameEta_ToF = (TH1D*)File_InPut->Get("h_FrameEta_ToF");
  TH1D *h_FramePhi_ToF = (TH1D*)File_InPut->Get("h_FramePhi_ToF");

  cout << "eta bins = " << h_FrameEta_ToF->GetNbinsX() << endl;;
  cout << "phi bins = " << h_FramePhi_ToF->GetNbinsX() << endl;;

  const int neta = h_FrameEta_ToF->GetNbinsX();
  const int nphi = h_FramePhi_ToF->GetNbinsX();

  TH3D *h_mTracks_TPC[2][3][10]; // pt, eta, phi distribution as a function of charge | pid | centrality
  TH3D *h_mTracks_ToF[2][3][10];
  TH3D *h_mTracks_TPC_both[3][10]; // pt, eta, phi distribution as a function of charge | pid | centrality
  TH3D *h_mTracks_ToF_both[3][10];

  int pid_start, pid_stop;
  if(pid == 0) { pid_start = 0; pid_stop = 1; }
  if(pid == 1) { pid_start = 1; pid_stop = 2; }
  if(pid == 2) { pid_start = 0; pid_stop = 2; }

  for(int i_pid = pid_start; i_pid < pid_stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 9; ++i_cent)
      {
	string HistName_Trakcs_TPC = Form("h_mTracks_TPC_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mTracks_TPC[i_charge][i_pid][i_cent] = (TH3D*)File_InPut->Get(HistName_Trakcs_TPC.c_str());

	string HistName_Trakcs_ToF = Form("h_mTracks_ToF_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mTracks_ToF[i_charge][i_pid][i_cent] = (TH3D*)File_InPut->Get(HistName_Trakcs_ToF.c_str());
      }
      for(int i_cent = 2; i_cent <= 5; ++i_cent) // make centrality 9 to be 20-60%
      {
	if(i_cent == 2)
	{
	  string HistName_Trakcs_TPC = Form("h_mTracks_TPC_%s%s_Cent_9",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str());
	  h_mTracks_TPC[i_charge][i_pid][9] = (TH3D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->Clone(HistName_Trakcs_TPC.c_str());
	  h_mTracks_TPC[i_charge][i_pid][9]->Sumw2();

	  string HistName_Trakcs_ToF = Form("h_mTracks_ToF_%s%s_Cent_9",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str());
	  h_mTracks_ToF[i_charge][i_pid][9] = (TH3D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->Clone(HistName_Trakcs_ToF.c_str());
	  h_mTracks_ToF[i_charge][i_pid][9]->Sumw2();
	}
	else
	{
	  h_mTracks_TPC[i_charge][i_pid][9]->Add(h_mTracks_TPC[i_charge][i_pid][i_cent]);
	  h_mTracks_ToF[i_charge][i_pid][9]->Add(h_mTracks_ToF[i_charge][i_pid][i_cent]);
	}
      }
    }
  }
  for(int i_pid = pid_start; i_pid < pid_stop; ++i_pid)
  {
    for(int i_cent = 0; i_cent < 9; ++i_cent)
    {
      string histname = Form("TPC_counts_cent_%d",i_cent);
      h_mTracks_TPC_both[i_pid][i_cent] = (TH3D*) h_mTracks_TPC[0][i_pid][i_cent]->Clone(histname.c_str());
      h_mTracks_TPC_both[i_pid][i_cent]->Sumw2();
      h_mTracks_TPC_both[i_pid][i_cent]->Add(h_mTracks_TPC[1][i_pid][i_cent]);

      histname = Form("ToF_counts_cent_%d",i_cent);
      h_mTracks_ToF_both[i_pid][i_cent] = (TH3D*) h_mTracks_ToF[0][i_pid][i_cent]->Clone(histname.c_str());
      h_mTracks_ToF_both[i_pid][i_cent]->Sumw2();
      h_mTracks_ToF_both[i_pid][i_cent]->Add(h_mTracks_ToF[1][i_pid][i_cent]);

    }
    for(int i_cent = 2; i_cent <= 5; ++i_cent) // make centrality 9 to be 20-60%
    {
      if(i_cent == 2)
      {
        string histname = Form("TPC_counts_cent_9");
        h_mTracks_TPC_both[i_pid][9] = (TH3D*)h_mTracks_TPC_both[i_pid][i_cent]->Clone(histname.c_str());
        h_mTracks_TPC_both[i_pid][9]->Sumw2();

        histname = Form("ToF_counts_cent_9");
        h_mTracks_ToF_both[i_pid][9] = (TH3D*)h_mTracks_ToF_both[i_pid][i_cent]->Clone(histname.c_str());
        h_mTracks_ToF_both[i_pid][9]->Sumw2();
      }
      else
      {
        h_mTracks_TPC_both[i_pid][9]->Add(h_mTracks_TPC_both[i_pid][i_cent]);
        h_mTracks_ToF_both[i_pid][9]->Add(h_mTracks_ToF_both[i_pid][i_cent]);
      }
    } 
  }

  cout << "Right before plotting" << endl;
  TH1D *h_mEffPhi[10];

  TCanvas *c1 = new TCanvas("phieff","phieff",10,10,400,400);
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetBottomMargin(0.15);
  c1->cd(1)->SetGrid(0,0);
  c1->cd(1)->SetTicks(1,1);

  string HistNameEff = Form("h_mEffPhi_Cent_%d",9);
  TH1D* tempTPCphi = (TH1D*) h_mTracks_TPC_both[0][9]->ProjectionZ("tempTPCphi",0,-1,0,-1);
  TH1D* tempTOFphi = (TH1D*) h_mTracks_ToF_both[0][9]->ProjectionZ("tempTOFphi",0,-1,0,-1);
  h_mEffPhi[9] = (TH1D*) tempTOFphi->Clone(HistNameEff.c_str());
  h_mEffPhi[9]->Sumw2();
  h_mEffPhi[9]->Reset();
  h_mEffPhi[9]->Divide(tempTOFphi,tempTPCphi,1,1,"B");   
 
  h_mEffPhi[9]->SetStats(0);
  if(EP == 0 || EP == 2 || EP == 3 || EP == 4) h_mEffPhi[9]->GetXaxis()->SetTitle("#phi");
  if(EP == 1) h_mEffPhi[9]->GetXaxis()->SetTitle("#phi-#Psi_{2}");
  h_mEffPhi[9]->GetYaxis()->SetTitle("RC/MC");
  h_mEffPhi[9]->Draw("pE");
  c1->SaveAs(Form("phiToFMatchingeff%s_cent9.pdf",eventplane.c_str()));


  // calculate efficiency vs centrality
  TH1DMap h_mCounts_TPC_cent; // efficiency vs. cent
  TH1DMap h_mCounts_ToF_cent;
  TH1DMap h_mEfficiency_cent;
  for(int i_pid = pid_start; i_pid < pid_stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_TPC_cent[HistName_TPC] = (TH1D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->ProjectionX(HistName_TPC.c_str())->Clone();

	string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_ToF_cent[HistName_ToF] = (TH1D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->ProjectionX(HistName_ToF.c_str())->Clone();

	string HistName = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_cent[HistName] = (TH1D*)h_mCounts_ToF_cent[HistName_ToF]->Clone(HistName.c_str());
	h_mEfficiency_cent[HistName]->SetTitle(HistName.c_str());
	h_mEfficiency_cent[HistName]->Reset();
	h_mEfficiency_cent[HistName]->Divide(h_mCounts_ToF_cent[HistName_ToF],h_mCounts_TPC_cent[HistName_TPC],1,1,"B");
      }
    }
  }

  // calculate efficiencya vs. eta
  TH1DMap h_mCounts_TPC_eta; // efficiency vs. cent & eta
  TH1DMap h_mCounts_ToF_eta;
  TH1DMap h_mEfficiency_eta;
  for(int i_pid = pid_start; i_pid < pid_stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	for(int i_eta = 0; i_eta < neta; ++i_eta)
	{
	  string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta);
	  h_mCounts_TPC_eta[HistName_TPC] = (TH1D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->ProjectionX(HistName_TPC.c_str(),i_eta+1,i_eta+1)->Clone();

	  string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta);
	  h_mCounts_ToF_eta[HistName_ToF] = (TH1D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->ProjectionX(HistName_ToF.c_str(),i_eta+1,i_eta+1)->Clone();

	  string HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta);
	  h_mEfficiency_eta[HistName] = (TH1D*)h_mCounts_ToF_eta[HistName_ToF]->Clone(HistName.c_str());
	  h_mEfficiency_eta[HistName]->SetTitle(HistName.c_str());
	  h_mEfficiency_eta[HistName]->Reset();
	  h_mEfficiency_eta[HistName]->Divide(h_mCounts_ToF_eta[HistName_ToF],h_mCounts_TPC_eta[HistName_TPC],1,1,"B");
	}
      }
    }
  }

  TH1DMap h_mCounts_TPC; // efficiency vs. cent & eta & phi
  TH1DMap h_mCounts_ToF;
  TH1DMap h_mEfficiency;
  for(int i_pid = pid_start; i_pid < pid_stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	for(int i_eta = 0; i_eta < neta; ++i_eta)
	{
	  for(int i_phi = 0; i_phi < nphi; ++i_phi)
	  {
	    string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mCounts_TPC[HistName_TPC] = (TH1D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->ProjectionX(HistName_TPC.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1)->Clone();
            //h_mCounts_TPC[HistName_TPC]->Print();
	    string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mCounts_ToF[HistName_ToF] = (TH1D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->ProjectionX(HistName_ToF.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1)->Clone();
            //h_mCounts_ToF[HistName_ToF]->Print();
	    string HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mEfficiency[HistName] = (TH1D*)h_mCounts_ToF[HistName_ToF]->Clone(HistName.c_str());
	    h_mEfficiency[HistName]->SetTitle(HistName.c_str());
	    h_mEfficiency[HistName]->Reset();
	    h_mEfficiency[HistName]->Divide(h_mCounts_ToF[HistName_ToF],h_mCounts_TPC[HistName_TPC],1,1,"B");
	    cout << "Calculating => " << HistName.c_str() << endl;
	  }
	}
      }
    }
  }

  // calculate efficiencya vs. eta
  TH2D *h_mCounts_TPC_pteta[10]; // efficiency vs. cent & eta
  TH2D *h_mCounts_ToF_pteta[10];
  TH2D *h_mEfficiency_pteta[10];
  TH2D *h_mCounts_TPC_ptphi[10]; // efficiency vs. cent & eta
  TH2D *h_mCounts_ToF_ptphi[10];
  TH2D *h_mEfficiency_ptphi[10];
  TH2D *h_mCounts_TPC_etaphi[10]; // efficiency vs. cent & eta
  TH2D *h_mCounts_ToF_etaphi[10];
  TH2D *h_mEfficiency_etaphi[10];

  TH1D *h_mCounts_TPC_pt1D[10]; // efficiency vs. cent & eta
  TH1D *h_mCounts_ToF_pt1D[10];
  TH1D *h_mEfficiency_pt1D[10];
  TH1D *h_mCounts_TPC_eta1D[10]; // efficiency vs. cent & eta
  TH1D *h_mCounts_ToF_eta1D[10];
  TH1D *h_mEfficiency_eta1D[10];
  TH1D *h_mCounts_TPC_phi1D[10]; // efficiency vs. cent & eta
  TH1D *h_mCounts_ToF_phi1D[10];
  TH1D *h_mEfficiency_phi1D[10];

  for(int i_pid = pid_start; i_pid < pid_stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
        //PT
	string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_pt",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_TPC_pt1D[i_cent] = (TH1D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->ProjectionX(HistName_TPC.c_str(),0,-1,0,-1)->Clone();

	string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_pt",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_ToF_pt1D[i_cent] = (TH1D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->ProjectionX(HistName_ToF.c_str(),0,-1,0,-1)->Clone();

	string HistName = Form("h_mEfficiency_%s%s_Cent_%d_pt",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_pt1D[i_cent] = (TH1D*)h_mCounts_ToF_pt1D[i_cent]->Clone(HistName.c_str());
	h_mEfficiency_pt1D[i_cent]->SetTitle(HistName.c_str());
	h_mEfficiency_pt1D[i_cent]->Reset();
	h_mEfficiency_pt1D[i_cent]->Divide(h_mCounts_ToF_pt1D[i_cent],h_mCounts_TPC_pt1D[i_cent],1,1,"B");

        //Eta
	HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_eta",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_TPC_eta1D[i_cent] = (TH1D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->ProjectionY(HistName_TPC.c_str(),0,-1,0,-1)->Clone();

	HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_eta",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_ToF_eta1D[i_cent] = (TH1D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->ProjectionY(HistName_ToF.c_str(),0,-1,0,-1)->Clone();

	HistName = Form("h_mEfficiency_%s%s_Cent_%d_eta",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_eta1D[i_cent] = (TH1D*)h_mCounts_ToF_eta1D[i_cent]->Clone(HistName.c_str());
	h_mEfficiency_eta1D[i_cent]->SetTitle(HistName.c_str());
	h_mEfficiency_eta1D[i_cent]->Reset();
	h_mEfficiency_eta1D[i_cent]->Divide(h_mCounts_ToF_eta1D[i_cent],h_mCounts_TPC_eta1D[i_cent],1,1,"B");
      
        //Phi
	HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_phi",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_TPC_phi1D[i_cent] = (TH1D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->ProjectionZ(HistName_TPC.c_str(),0,-1,0,-1)->Clone();

	HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_phi",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_ToF_phi1D[i_cent] = (TH1D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->ProjectionZ(HistName_ToF.c_str(),0,-1,0,-1)->Clone();

	HistName = Form("h_mEfficiency_%s%s_Cent_%d_phi",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_phi1D[i_cent] = (TH1D*)h_mCounts_ToF_phi1D[i_cent]->Clone(HistName.c_str());
	h_mEfficiency_phi1D[i_cent]->SetTitle(HistName.c_str());
	h_mEfficiency_phi1D[i_cent]->Reset();
	h_mEfficiency_phi1D[i_cent]->Divide(h_mCounts_ToF_phi1D[i_cent],h_mCounts_TPC_phi1D[i_cent],1,1,"B");

        //PT Eta
	HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_pteta",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_TPC_pteta[i_cent] = (TH2D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->Project3D("yx")->Clone();

	HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_pteta",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_ToF_pteta[i_cent] = (TH2D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->Project3D("yx")->Clone();

	HistName = Form("h_mEfficiency_%s%s_Cent_%d_pteta",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_pteta[i_cent] = (TH2D*)h_mCounts_ToF_pteta[i_cent]->Clone(HistName.c_str());
	h_mEfficiency_pteta[i_cent]->SetTitle(HistName.c_str());
	h_mEfficiency_pteta[i_cent]->Reset();
	h_mEfficiency_pteta[i_cent]->Divide(h_mCounts_ToF_pteta[i_cent],h_mCounts_TPC_pteta[i_cent],1,1,"B");
	
        //PT Phi
	HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_ptphi",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_TPC_ptphi[i_cent] = (TH2D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->Project3D("zx")->Clone();

	HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_ptphi",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_ToF_ptphi[i_cent] = (TH2D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->Project3D("zx")->Clone();

	HistName = Form("h_mEfficiency_%s%s_Cent_%d_pthpi",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_ptphi[i_cent] = (TH2D*)h_mCounts_ToF_ptphi[i_cent]->Clone(HistName.c_str());
	h_mEfficiency_ptphi[i_cent]->SetTitle(HistName.c_str());
	h_mEfficiency_ptphi[i_cent]->Reset();
	h_mEfficiency_ptphi[i_cent]->Divide(h_mCounts_ToF_ptphi[i_cent],h_mCounts_TPC_ptphi[i_cent],1,1,"B");

        //Eta Phi
	HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_etaphi",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_TPC_etaphi[i_cent] = (TH2D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->Project3D("zy")->Clone();

	HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_etaphi",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_ToF_etaphi[i_cent] = (TH2D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->Project3D("zy")->Clone();

	HistName = Form("h_mEfficiency_%s%s_Cent_%d_pthpi",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_etaphi[i_cent] = (TH2D*)h_mCounts_ToF_etaphi[i_cent]->Clone(HistName.c_str());
	h_mEfficiency_etaphi[i_cent]->SetTitle(HistName.c_str());
	h_mEfficiency_etaphi[i_cent]->Reset();
	h_mEfficiency_etaphi[i_cent]->Divide(h_mCounts_ToF_etaphi[i_cent],h_mCounts_TPC_etaphi[i_cent],1,1,"B");
	
      }
    }
  }


  //---------------output------------------------
  // string outputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string outputfile = Form("../../output/%s/Eff_%s_ToFMatch%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),eventplane.c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_pid = pid_start; i_pid < pid_stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	string HistName = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_cent[HistName]->Write();
	for(int i_eta = 0; i_eta < neta; ++i_eta)
	{
	  HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta);
	  h_mEfficiency_eta[HistName]->Write();
	  for(int i_phi = 0; i_phi < nphi; ++i_phi)
	  {
	    HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mEfficiency[HistName]->Write();
	  }
	}
      }
    }
  }
  h_FrameEta_ToF->Write();
  h_FramePhi_ToF->Write();
  File_OutPut->Close();
}
