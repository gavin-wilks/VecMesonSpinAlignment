#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TProfile2D.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
//#ifdef MAKECINT
//#pragma link C++ class std::map<std::string,TH1F*>+;
//#endif

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

#ifndef _SaveQA_
#define _SaveQA_  0
#endif

using namespace std;

void subBackGround_ptyspectra_TPCOnly(int energy = 4, int pid = 0, int year = 0, string date = "20241219", bool random3D = false, int order = 2, string etamode = "eta1_eta1")
{
  //old spectra has 2D at the end adn is 20240507
   
  std::string EP[2] = {"","2nd"};
  TGaxis::SetMaxDigits(4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  //string InPutFile_SE = Form("../data/Yields_Phi_SE_19GeV_20220527.root"); //original eta < 1.0
  //string InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s_2D.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  //string InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_eta1p5_eta1p5_ptyspectra_TPCOnly.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  //string InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_ptyspectra_TPCOnly_70vz_25phipsi_eta1_eta1_y1p5.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  //string InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_ptyspectra_TPCOnly_eta1_eta1_dipangle.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  string InPutFile_SE = Form("../data/Yields_Phi_SE_19GeV_20241211_ptyspectra_TPCOnly_70vz_25phipsi_eta1_eta1_y1p5.root");
  //if(energy == 3) InPutFile_SE = Form("../data/Yields_Phi_SE_14GeV_%s.root",etamode.c_str());
  if(order == 1) InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(random3D) InPutFile_SE = Form("../data/3DRandom/Yields_Phi_SE_%s_%s_3DRandom.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_19GeV_20220408.root"); //original eta < 1.0
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_eta1p5_eta1p5_ptyspectra_TPCOnly.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_ptyspectra_TPCOnly_70vz_25phipsi_eta1_eta1_y1p5.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  //string InPutFile_ME = Form("../data/Yields_Phi_SE_%s_20241218_ptyspectra_TPCOnly_likesign_eta1_eta1_y1p5.root",vmsa::mBeamEnergy[energy].c_str());
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_ptyspectra_TPCOnly_eta1_eta1_dipangle.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_ptyspectra_TPCOnly_eta1_eta1_dipangle.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  string InPutFile_ME = Form("../data/Yields_Phi_ME_19GeV_20241211_ptyspectra_TPCOnly_70vz_25phipsi_eta1_eta1_y1p5.root");
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s_2D.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  //if(energy == 3) InPutFile_ME = Form("../data/Yields_Phi_ME_14GeV_%s.root",etamode.c_str());
  if(order == 1) InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(random3D) InPutFile_ME = Form("../data/3DRandom/Yields_Phi_ME_%s_%s_3DRandom.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH1FMap h_mInPut_SE, h_mInPut_ME;
  TH1FMap h_mMass_SE, h_mMass_ME, h_mMass_SM;

  float scalingratio[vmsa::rebinpttotal][vmsa::rebinytotal] = {0.0};

  for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin 
  {
    for(Int_t i_cent = 9; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_eta = 0; i_eta < 16; i_eta++) // phi-psi bin
      {
	for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
	{
	  for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
	  {
            if( i_dca != 0 && i_sig != 0 ) continue;
	    string KEY_InPutSE = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_SE",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
	    string KEY_InPutSE_temp = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_SE_temp",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
            h_mInPut_SE[KEY_InPutSE] = (TH1F*)File_SE->Get(KEY_InPutSE.c_str())->Clone(KEY_InPutSE_temp.c_str()); 
            //h_mInPut_SE[KEY_InPutSE]->Rebin(2);
            cout << "SE Bin number = " << h_mInPut_SE[KEY_InPutSE]->GetNbinsX() << endl;

	    string KEY_InPutME = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_ME",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
	    //string KEY_InPutME = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_ME",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
	    string KEY_InPutME_temp = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_ME_temp",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str());
	    h_mInPut_ME[KEY_InPutME] = (TH1F*)File_ME->Get(KEY_InPutME.c_str())->Clone(KEY_InPutME_temp.c_str()); 
            //h_mInPut_ME[KEY_InPutME]->Rebin(2);          
            cout << "ME Bin number = " << h_mInPut_ME[KEY_InPutME]->GetNbinsX() << endl;
	    for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	    {
	      string KEY_SE = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	      h_mMass_SE[KEY_SE] = (TH1F*)h_mInPut_SE[KEY_InPutSE]->Clone(KEY_SE.c_str());
	      string KEY_ME = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	      h_mMass_ME[KEY_ME] = (TH1F*)h_mInPut_ME[KEY_InPutME]->Clone(KEY_ME.c_str());
            }
	  }
	}
      }
    }
  }
  for(Int_t i_cent = 9; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
      {
        for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
        {
          for(int ptbin = 0; ptbin < vmsa::rebinpttotal; ptbin++)
          {
            for(Int_t i_pt = vmsa::rebinpt[ptbin]; i_pt < vmsa::rebinpt[ptbin+1]; i_pt++) // pt bin 
            {
              for(int etabin = 0; etabin < vmsa::rebinytotal; etabin++)
              {
                for(Int_t i_eta = vmsa::rebiny[etabin]; i_eta < vmsa::rebiny[etabin+1]; i_eta++) // phi-psi bin
                {
                  string KEY_SE = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                  //cout << KEY_SE << endl;
                  string KEY_SE_rebin = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE_Rebin",ptbin,etabin,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                  //cout << KEY_SE_rebin << endl;
                  string KEY_ME = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                  //cout << KEY_ME << endl;
                  string KEY_ME_rebin = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME_Rebin",ptbin,etabin,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
                  //cout << KEY_ME_rebin << endl;
                  if(i_pt == vmsa::rebinpt[ptbin] && i_eta == vmsa::rebiny[etabin])
                  { 
                    cout << "Create: ptbin = " << ptbin << ", etabin = " << etabin << "     using pt = " << i_pt << ", eta = " << i_eta << endl;
                    h_mMass_SE[KEY_SE_rebin] = (TH1F*)h_mMass_SE[KEY_SE]->Clone(KEY_SE_rebin.c_str());
                    h_mMass_ME[KEY_ME_rebin] = (TH1F*)h_mMass_ME[KEY_ME]->Clone(KEY_ME_rebin.c_str());
                  }
                  else
                  {
                    cout << "Add: ptbin = " << ptbin << ", etabin = " << etabin << "     using pt = " << i_pt << ", eta = " << i_eta << endl;
                    h_mMass_SE[KEY_SE_rebin]->Add(h_mMass_SE[KEY_SE]);
                    h_mMass_ME[KEY_ME_rebin]->Add(h_mMass_ME[KEY_ME]);
                  }
                }
              }
            }
          }
        }
        for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
        {
          for(int ptbin = 0; ptbin < vmsa::rebinpttotal; ptbin++)
          {
            for(int etabin = 0; etabin < vmsa::rebinytotal; etabin++)
            {
              string KEY_SM = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM_Rebin",ptbin,etabin,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
              //cout << KEY_SM << endl;
              string KEY_SE = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE_Rebin",ptbin,etabin,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
              //cout << KEY_SE << endl;
              string KEY_ME = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME_Rebin",ptbin,etabin,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
              //cout << KEY_ME << endl;


	      h_mMass_SM[KEY_SM] = (TH1F*)h_mMass_SE[KEY_SE]->Clone();
	      if(i_norm < 2)
	      {
		int Norm_bin_start = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Start[pid][i_norm]);
		int Norm_bin_stop  = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Stop[pid][i_norm]);

		float Inte_SE = h_mMass_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);
		float Inte_ME = h_mMass_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);

		if(Inte_ME != 0.0) h_mMass_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
		h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);

                scalingratio[ptbin][etabin] = Inte_SE/Inte_ME;
	      }
	      else
	      {
		float Inte_SE = 0.0;
		float Inte_ME = 0.0;

		for(int i_inte = 0; i_inte < 2; ++i_inte)
		{
		  int Norm_bin_start = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Start[pid][i_inte]);
		  int Norm_bin_stop  = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Stop[pid][i_inte]);
		  Inte_SE += h_mMass_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);
		  Inte_ME += h_mMass_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);
		}

		if(Inte_ME != 0.0) h_mMass_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
		h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);
	      }
	    }
          }
        }
      }
    }
  }
  cout << endl << "const float scalingME[vmsa::rebinpttotal][vmsa::rebinytotal] = {";
  for(int ptbin = 0; ptbin < vmsa::rebinpttotal; ptbin++)
  {
    cout << "{";
    for(int etabin = 0; etabin < vmsa::rebinytotal; etabin++)
    {
      if(etabin < vmsa::rebinytotal-1) cout << scalingratio[ptbin][etabin] << ",";
      else                             cout << scalingratio[ptbin][etabin];
    }
    if(ptbin < vmsa::rebinpttotal -1 ) cout << "}," << endl;
    else                                cout << "}";
  }
  cout << "};" << endl << endl;
  // QA Plots for pT bins
  for(int i_cent = 9; i_cent < 10; i_cent++)
  {
    TCanvas *c_pT = new TCanvas("c_pT","c_pT",10,10,1400,1400);
    c_pT->Divide(5,5);

    string outputname = Form("figures/%s/%s/pTstudy/TPCOnly_ptyspectra_InvMassDistributions_%s_Order%d_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,i_cent);
    string output_start = Form("%s[",outputname.c_str());
    
    c_pT->Print(output_start.c_str());
    for(int ieta = 0; ieta < vmsa::rebinytotal; ieta++)
    {
      for(int i_pt = 0; i_pt < vmsa::rebinpttotal; i_pt++) // pt loop
      {
        c_pT->cd(i_pt+1);
        c_pT->cd(i_pt+1)->SetLeftMargin(0.15);
        c_pT->cd(i_pt+1)->SetBottomMargin(0.15);
        c_pT->cd(i_pt+1)->SetTicks(1,1);
        c_pT->cd(i_pt+1)->SetGrid(0,0);
        string KEY_SE_QA = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SE_Rebin",i_pt,ieta,i_cent,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
        h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
        h_mMass_SE[KEY_SE_QA]->SetTitle(Form("%.2f<y<%.2f",vmsa::rebinyval[ieta],vmsa::rebinyval[ieta+1]));
        h_mMass_SE[KEY_SE_QA]->DrawCopy();

        string KEY_ME_QA = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_ME_Rebin",i_pt,ieta,i_cent,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
        h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
        h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
        h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
        h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

        string KEY_SM_QA = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM_Rebin",i_pt,ieta,i_cent,EP[order-1].c_str(),vmsa::Dca_start,vmsa::nSigKaon_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
        h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
        h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
        h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
        h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");

        string pT_range = Form("[%.2f,%.2f]",vmsa::rebinptval[i_pt],vmsa::rebinptval[i_pt+1]);
        plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);

        if(vmsa::Norm_QA == 0 || vmsa::Norm_QA == 2)
        {
          PlotLine(vmsa::Norm_Start[pid][0],vmsa::Norm_Start[pid][0],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
          PlotLine(vmsa::Norm_Stop[pid][0],vmsa::Norm_Stop[pid][0],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
        }
        if(vmsa::Norm_QA == 1 || vmsa::Norm_QA == 2)
        {
          PlotLine(vmsa::Norm_Start[pid][1],vmsa::Norm_Start[pid][1],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
          PlotLine(vmsa::Norm_Stop[pid][1],vmsa::Norm_Stop[pid][1],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
        }
      }
      c_pT->Update();
      c_pT->Print(outputname.c_str());
    }
    string output_stop = Form("%s]",outputname.c_str());
    c_pT->Print(output_stop.c_str()); // close pdf file
  }
 
  // write background subtracted histograms to output file
  string outputfile = Form("../output/AuAu%s/%s/TPCOnly_ptyspectra_InvMassSubBg_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  //string outputfile = Form("../output/AuAu%s/%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  if(order == 1) outputfile = Form("../output/AuAu%s/%s/InvMassSubBg_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if(random3D) outputfile = Form("../output/AuAu%s/%s/3DRandom/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_pt = 0; i_pt < vmsa::rebinpttotal; i_pt++) // pt loop
  {
    for(int i_cent = 9; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
      {
	for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
	{
          if( i_dca != 0 && i_sig != 0 ) continue;
	  for(int i_eta = 0; i_eta < vmsa::rebinytotal; i_eta++) // cos(theta*) loop
	  {
	    for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	    {
	      string KEY_SM = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_SM_Rebin",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	      string KEY = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d",i_pt,i_eta,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm);
	      h_mMass_SM[KEY_SM]->SetName(KEY.c_str());;
	      h_mMass_SM[KEY_SM]->Write();
	    }
	  }
	}
      }
    }
  }
  File_OutPut->Close();
}
