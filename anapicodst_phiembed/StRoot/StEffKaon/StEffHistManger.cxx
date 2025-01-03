#include "StEffHistManger.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include <iostream>
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StRoot/Utility/phi_data_constants_19GeV.h"
#include "TString.h"
#include <string>
#include "TFile.h" 

using namespace std;


ClassImp(StEffHistManger)
//
StEffHistManger::StEffHistManger()
{
  /* */
}

StEffHistManger::~StEffHistManger()
{
  /* */
}

void StEffHistManger::InitHistQA()
{
  string charge[2] = {"plus","minus"};
  {
    for(int i_cent = 0; i_cent < 10; ++i_cent)
    {
      string HistName = Form("h_mVertex_Cent_%d",i_cent);
      h_mVertex[i_cent] = new TH3F(HistName.c_str(),HistName.c_str(),50,-2.5,2.5,50,-2.5,2.5,280,-70.0,70.0);

      HistName = Form("h_mNToFMatch_Cent_%d",i_cent);
      h_mNToFMatch[i_cent] = new TH1F(HistName.c_str(),HistName.c_str(),501,-0.5,500.5);

      HistName = Form("h_mRefMult_Cent_%d",i_cent);
      h_mRefMult[i_cent] = new TH1F(HistName.c_str(),HistName.c_str(),501,-0.5,500.5);

      for(int i_step = 0; i_step < 10; ++i_step)
      {
        for(int i_charge = 0; i_charge < 2; ++i_charge)
        {
          HistName = Form("h_mDca_K%s_Cent_%d_%d",charge[i_charge].c_str(),i_cent,i_step);
          h_mDca[i_step][i_charge][i_cent] = new TH1F(HistName.c_str(),HistName.c_str(),60,0.0,3.0);
          
          HistName = Form("h_mNHits_K%s_Cent_%d_%d",charge[i_charge].c_str(),i_cent,i_step);
          h_mNHits[i_step][i_charge][i_cent] = new TH1F(HistName.c_str(),HistName.c_str(),101,-0.5,100.5);

          HistName = Form("h_mNHitsRatio_K%s_Cent_%d_%d",charge[i_charge].c_str(),i_cent,i_step);
          h_mNHitsRatio[i_step][i_charge][i_cent] = new TH1F(HistName.c_str(),HistName.c_str(),100,0.0,1.0);

          HistName = Form("h_mDEdx_K%s_Cent_%d_%d",charge[i_charge].c_str(),i_cent,i_step);
          h_mDEdx_Kaon[i_step][i_charge][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),500,0.0,10.0,1000,0,40);
        }
      }
    }
  }
}

void StEffHistManger::FillEventHistQA(float reweight, int cent, float vx, float vy, float vz, float ntof, float refmult)
{
  h_mVertex[cent]->Fill(vx,vy,vz,reweight);
  h_mNToFMatch[cent]->Fill(ntof,reweight);
  h_mRefMult[cent]->Fill(refmult,reweight);
  if(cent >= 2 && cent <= 5)
  {
    h_mVertex[9]->Fill(vx,vy,vz,reweight);
    h_mNToFMatch[9]->Fill(ntof,reweight);
    h_mRefMult[9]->Fill(refmult,reweight);
  }
} 

void StEffHistManger::FillTrackHistQA(float reweight, int cent, int charge, float dca, float nhits, float nhitsratio, float p, float dedx, int step)
{
  h_mDca[step][charge][cent]->Fill(dca,reweight);
  h_mNHits[step][charge][cent]->Fill(nhits,reweight);
  h_mNHitsRatio[step][charge][cent]->Fill(nhitsratio,reweight);
  h_mDEdx_Kaon[step][charge][cent]->Fill(p,dedx,reweight);
  if(cent >= 2 && cent <= 5)
  {
    h_mDca[step][charge][9]->Fill(dca,reweight);
    h_mNHits[step][charge][9]->Fill(nhits,reweight);
    h_mNHitsRatio[step][charge][9]->Fill(nhitsratio,reweight);
    h_mDEdx_Kaon[step][charge][9]->Fill(p,dedx,reweight);
  }
} 

void StEffHistManger::WriteHistQA()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    h_mVertex[i_cent]->Write();
    h_mNToFMatch[i_cent]->Write();
    h_mRefMult[i_cent]->Write();
    for(int i_step = 0; i_step < 10; ++i_step)
    {
      for(int i_charge = 0; i_charge < 2; ++i_charge)
      {
        h_mDca[i_step][i_charge][i_cent]->Write();
        h_mNHits[i_step][i_charge][i_cent]->Write();
        h_mNHitsRatio[i_step][i_charge][i_cent]->Write();     
        h_mDEdx_Kaon[i_step][i_charge][i_cent]->Write();
      }
    }
  }
}

void StEffHistManger::InitHist()
{
  //TFile *inputweights = TFile::Open(Form("StRoot/Utility/pt_phi_weights_v2_function_75.root"));
  //TFile *inputweights = TFile::Open(Form("StRoot/Utility/pt_phi_weights_v2_FunctionEmbed.root"));
  TFile *inputweightsptyphi = TFile::Open(Form("StRoot/Utility/pt_y_phi_weights_v2_FunctionEmbed.root"));

  const int nybins = 16;
  double ybins[nybins+1] = {0.0};
  for(int iy = 0; iy < nybins+1; iy++)
  {
    ybins[iy] = double(iy-8.)/8.;
    cout << "iy = " << iy << ", y = " << ybins[iy] << endl;
  } 
  
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    if(i_cent >= 2 && i_cent <= 5)
    {
      std::string HistName = Form("h_ratio_cent%d",i_cent);
      //std::string CloneName = Form("h_ratio_cent%d_v2%d",i_cent);
      //h_mpTv2Weights[i_cent] = (TH2F*) inputweights->Get(HistName.c_str());
      h_mpTyv2Weights[i_cent] = (TH3F*) inputweightsptyphi->Get(HistName.c_str());
      //h_mpTv2Weights[i_cent] = ((TH2F*) inputweights->Get(HistName.c_str()))->Clone(;
    }
  }

  //TFile *inputweights_v2var[10];// = TFile::Open(Form("StRoot/Utility/pt_phi_weights_v2_function_75.root"));
  //std::string v2values[10] = {"0.0000","0.0333","0.0667","0.1000","0.1333","0.1667","0.2000","0.2333","0.2667","0.3000"};
  //for(int iv2 = 0; iv2 < 10; iv2++)
  //{
    //inputweights_v2var = TFile::Open(Form("StRoot/Utility/pt_phi_weights_v2_function_75_v2%s.root",v2values.c_str()));
    //for(int i_cent = 0; i_cent < 10; ++i_cent)
    //{
    //  if(i_cent >= 2 && i_cent <= 5)
    //  {
    //    std::string HistName = Form("h_ratio_cent%d",i_cent);
    //    std::string CloneName = Form("h_ratio_cent%d_v2%d",i_cent);
    //    //h_mpTv2Weights[i_cent] = (TH2F*) inputweights->Get(HistName.c_str());
    //    h_mpTv2Weights_Var[v2][i_cent] = (TH2F*)((TH2F*) inputweights->Get(HistName.c_str()))->Clone(CloneName.c_str());
    //  }
    //}
  //  TString KEY_MC = Form("mc_kplus_cent%d",i_cent);
  //  h_mKaon_MC[0][i_cent] = new TH2F(KEY_MC.Data(),KEY_MC.Data(),data_constants::kaon_pt_bins,data_constants::kaon_pt,int(2.0*data_constants::kaon_rapidity_bins),-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
  //  h_mKaon_MC[0][i_cent]->Sumw2();                                                                                                                                                             

  //  KEY_MC = Form("mc_kminus_cent%d",i_cent);
  //  h_mKaon_MC[1][i_cent] = new TH2F(KEY_MC.Data(),KEY_MC.Data(),data_constants::kaon_pt_bins,data_constants::kaon_pt,int(2.0*data_constants::kaon_rapidity_bins),-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
  //  h_mKaon_MC[1][i_cent]->Sumw2();                                                                                                                                                             

  //  string HistName = Form("h_mMcTracks_%d",i_cent);
  //  h_mMcTracks[i_cent] = new TH3F(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,2.0*vmsa::BinEta,-2.0*vmsa::mEtaMax,2.0*vmsa::mEtaMax,vmsa::BinPhi,0.0,2.0*TMath::Pi());
  //  h_mMcTracks[i_cent]->Sumw2();

  //  HistName = Form("h_mMcTracks_%d_Kplus",i_cent);
  //  h_mMcTracksKplus[i_cent] = new TH3F(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,vmsa::BinEta,vmsa::mEtaMax,vmsa::mEtaMax,vmsa::BinPhi,0.0,2.0*TMath::Pi());
  //  h_mMcTracksKplus[i_cent]->Sumw2();

  //  HistName = Form("h_mMcTracks_%d_Kminus",i_cent);
  //  h_mMcTracksKminus[i_cent] = new TH3F(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,vmsa::BinEta,vmsa::mEtaMax,vmsa::mEtaMax,vmsa::BinPhi,0.0,2.0*TMath::Pi());
  //  h_mMcTracksKminus[i_cent]->Sumw2();

  //  HistName = Form("h_mMcEffPtY_Cent_%d",i_cent);
  //  h_mMcEffPtY[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::pt_total,vmsa::pt_bin,nybins,ybins);
  //  h_mMcEffPtY[i_cent]->Sumw2();
  //  
  //  HistName = Form("h_mMcEffPtPhiPsi_Cent_%d",i_cent);
  //  h_mMcEffPtPhiPsi[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::pt_total,vmsa::pt_bin,40,0.0,2.0*TMath::Pi());
  //  h_mMcEffPtPhiPsi[i_cent]->Sumw2();
  //
  //  for(int i_pt = 2; i_pt < 6; i_pt++)
  //  { 
  //    HistName = Form("h_mMcEffPhiS_Cent_%d_Pt_%d",i_cent,i_pt);
  //    h_mMcEffPhiS[i_cent][i_pt] = new TH1F(HistName.c_str(),HistName.c_str(),10,-1.0,1.0);
  //    h_mMcEffPhiS[i_cent][i_pt]->Sumw2();

  //    HistName = Form("h_mMcEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
  //    h_mMcEffCos[i_cent][i_pt] = new TH1F(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2.,-1.0,1.0);
  //    h_mMcEffCos[i_cent][i_pt]->Sumw2();

  //    HistName = Form("h_mMcEffCosH_Cent_%d_Pt_%d",i_cent,i_pt);
  //    h_mMcEffCosH[i_cent][i_pt] = new TH1F(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2.,-1.0,1.0);
  //    h_mMcEffCosH[i_cent][i_pt]->Sumw2();

  //    HistName = Form("h_mMcEffCosCosH_Cent_%d_Pt_%d",i_cent,i_pt);
  //    h_mMcEffCosCosH[i_cent][i_pt] = new TH2F(HistName.c_str(),HistName.c_str(),9/*vmsa::BinCos*2.*/,-1.0,1.0,9/*vmsa::BinCos*2.*/,-1.0,1.0);
  //    h_mMcEffCosCosH[i_cent][i_pt]->Sumw2();

  //    HistName = Form("h_mMcEffCosPhiPrime_Cent_%d_Pt_%d",i_cent,i_pt);
  //    h_mMcEffCosPhiPrime[i_cent][i_pt] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2,-1.0,1.0,vmsa::BinCos*2,0.0,2.0*TMath::Pi());
  //    h_mMcEffCosPhiPrime[i_cent][i_pt]->Sumw2();

  //    HistName = Form("h_mMcEffCosHPhiPrime_Cent_%d_Pt_%d",i_cent,i_pt);
  //    h_mMcEffCosHPhiPrime[i_cent][i_pt] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2,-1.0,1.0,vmsa::BinCos*2,0.0,2.0*TMath::Pi());
  //    h_mMcEffCosHPhiPrime[i_cent][i_pt]->Sumw2();
  //  }
  //}
  double nsigmaarray[101];
  for(int i = 0; i < 101; i++)
  {
    nsigmaarray[i] = -5.0 + 0.1*i;
    cout << nsigmaarray[i] << endl;
  }  
  
  double m2array[51];
  for(int i = 0; i < 51; i++)
  {
    m2array[i] = 0.01 + 0.01*i;
    cout << m2array[i] << endl;
  }  
  int kaonybins = data_constants::kaon_rapidity_bins*2;
  double yarray[kaonybins+1];
  for(int i = 0; i < kaonybins+1; i++)
  {
    yarray[i] = -data_constants::kaon_rapidity_max + i*data_constants::kaon_rapidity_max/double(data_constants::kaon_rapidity_bins);
    cout << yarray[i] << endl;
  }
  //double phipsiarray[51];
  //for(int i = 0; i < 51; i++)
  //{
  //  phipsiarray[i] = TMath::Pi()*2.0/50.0*(float(i));
  //  cout << "phi-psi = " << phipsiarray[i] << endl;
  //}
  //double phiyarray[25];
  //for(int i = 0; i < 25; i++)
  //{
  //  phiyarray[i] = 1.5*(float(i-12))/12.;
  //  cout << "y = " << phiyarray[i] << endl;
  //}


  for(int i_step = 0; i_step < 10; ++i_step)
  {
    for(int i_cent = 0; i_cent < 10; ++i_cent)
    {
      //std::string HistName = Form("h_mRc%dEffPtPhiPsi_Cent_%d",i_step,i_cent);
      //h_mMcEffPtPhiPsi[i_step][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::pt_total,vmsa::pt_bin,50,0.0,2.0*TMath::Pi());
      //h_mMcEffPtPhiPsi[i_step][i_cent]->Sumw2();

      std::string HistName = Form("h_mRc%dEffPtYPhiPsi_Cent_%d",i_step,i_cent);
      //h_mMcEffPtYPhiPsi[i_step][i_cent] = new TH3F(HistName.c_str(),HistName.c_str(),vmsa::pt_total,vmsa::pt_bin,24,phiyarray,50,phipsiarray);
      h_mMcEffPtYPhiPsi[i_step][i_cent] = new TH3F(HistName.c_str(),HistName.c_str(),100,0.0,10.0,50,-1.5,1.5,120,0,2.0*TMath::Pi());
      h_mMcEffPtYPhiPsi[i_step][i_cent]->Sumw2();

      TString KEY_RC = Form("rc%d_kplus_cent%d",i_step,i_cent);                  
      h_mKaon_RC[i_step][0][i_cent] = new TH2F(KEY_RC.Data(),KEY_RC.Data(),data_constants::kaon_pt_bins,data_constants::kaon_pt,int(2.0*data_constants::kaon_rapidity_bins),-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
      h_mKaon_RC[i_step][0][i_cent]->Sumw2();

      KEY_RC = Form("rc%d_kminus_cent%d",i_step,i_cent);                  
      h_mKaon_RC[i_step][1][i_cent] = new TH2F(KEY_RC.Data(),KEY_RC.Data(),data_constants::kaon_pt_bins,data_constants::kaon_pt,int(2.0*data_constants::kaon_rapidity_bins),-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
      h_mKaon_RC[i_step][1][i_cent]->Sumw2();

      //KEY_RC = Form("rc%d_kplus_m2_cent%d",i_step,i_cent);                  
      //h_mKaon_RC_m2[i_step][0][i_cent] = new TH3F(KEY_RC.Data(),KEY_RC.Data(),data_constants::kaon_pt_bins,data_constants::kaon_pt,50,m2array,100,nsigmaarray);
      //h_mKaon_RC_m2[i_step][0][i_cent]->Sumw2();

      //KEY_RC = Form("rc%d_kminus_m2_cent%d",i_step,i_cent);                  
      //h_mKaon_RC_m2[i_step][1][i_cent] = new TH3F(KEY_RC.Data(),KEY_RC.Data(),data_constants::kaon_pt_bins,data_constants::kaon_pt,50,m2array,100,nsigmaarray);
      //h_mKaon_RC_m2[i_step][1][i_cent]->Sumw2();

      //KEY_RC = Form("rc%d_kplus_m2_cent%d",i_step,i_cent);                  
      //h_mKaon_RC_m2[i_step][0][i_cent] = new TH3F(KEY_RC.Data(),KEY_RC.Data(),data_constants::kaon_pt_bins,data_constants::kaon_pt,int(2.0*data_constants::kaon_rapidity_bins),yarray,50,m2array,100,nsigmaarray);
      //h_mKaon_RC_m2[i_step][0][i_cent]->Sumw2();

      //KEY_RC = Form("rc%d_kminus_m2_cent%d",i_step,i_cent);                  
      //h_mKaon_RC_m2[i_step][1][i_cent] = new TH3F(KEY_RC.Data(),KEY_RC.Data(),data_constants::kaon_pt_bins,data_constants::kaon_pt,int(2.0*data_constants::kaon_rapidity_bins),yarray,50,m2array,100,nsigmaarray);
      //h_mKaon_RC_m2[i_step][1][i_cent]->Sumw2();

      //string HistName = Form("h_mRc%dTracks_%d",i_step,i_cent);
      //h_mRcTracks[i_step][i_cent] = new TH3F(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,2.0*vmsa::BinEta,-2.0*vmsa::mEtaMax,2.0*vmsa::mEtaMax,vmsa::BinPhi,0.0,2.0*TMath::Pi());
      //h_mRcTracks[i_step][i_cent]->Sumw2();

      //HistName = Form("h_mRc%dTracks_%d_Kplus",i_step,i_cent);
      //h_mRcTracksKplus[i_step][i_cent] = new TH3F(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,vmsa::BinEta,vmsa::mEtaMax,vmsa::mEtaMax,vmsa::BinPhi,0.0,2.0*TMath::Pi());
      //h_mRcTracksKplus[i_step][i_cent]->Sumw2();

      //HistName = Form("h_mRc%dTracks_%d_Kminus",i_step,i_cent);
      //h_mRcTracksKminus[i_step][i_cent] = new TH3F(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,vmsa::BinEta,vmsa::mEtaMax,vmsa::mEtaMax,vmsa::BinPhi,0.0,2.0*TMath::Pi());
      //h_mRcTracksKminus[i_step][i_cent]->Sumw2();

      HistName = Form("h_mRc%dEffPtY_Cent_%d",i_step,i_cent);
      h_mRcEffPtY[i_step][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::pt_total,vmsa::pt_bin,nybins,ybins);
      h_mRcEffPtY[i_step][i_cent]->Sumw2();
      
      //for(int i_pt = 2; i_pt < 6; i_pt++)
      //{ 
      //  HistName = Form("h_mRc%dEffPhiS_Cent_%d_Pt_%d",i_step,i_cent,i_pt);
      //  h_mRcEffPhiS[i_step][i_cent][i_pt] = new TH1F(HistName.c_str(),HistName.c_str(),10,-1.0,1.0);
      //  h_mRcEffPhiS[i_step][i_cent][i_pt]->Sumw2();

      //  //HistName = Form("h_mRc%dEffCos_Cent_%d_Pt_%d",i_step,i_cent,i_pt);
      //  //h_mRcEffCos[i_step][i_cent][i_pt] = new TH1F(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2.,-1.0,1.0);
      //  //h_mRcEffCos[i_step][i_cent][i_pt]->Sumw2();

      //  //HistName = Form("h_mRc%dEffCosH_Cent_%d_Pt_%d",i_step,i_cent,i_pt);
      //  //h_mRcEffCosH[i_step][i_cent][i_pt] = new TH1F(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2.,-1.0,1.0);
      //  //h_mRcEffCosH[i_step][i_cent][i_pt]->Sumw2();

      //  //HistName = Form("h_mRc%dEffCosCosH_Cent_%d_Pt_%d",i_step,i_cent,i_pt);
      //  //h_mRcEffCosCosH[i_step][i_cent][i_pt] = new TH2F(HistName.c_str(),HistName.c_str(),9/*vmsa::BinCos*2.*/,-1.0,1.0,9/*vmsa::BinCos*2.*/,-1.0,1.0);
      //  //h_mRcEffCosCosH[i_step][i_cent][i_pt]->Sumw2();

      //  HistName = Form("h_mRc%dEffCosPhiPrime_Cent_%d_Pt_%d",i_step,i_cent,i_pt);
      //  h_mRcEffCosPhiPrime[i_step][i_cent][i_pt] = new TH2F(HistName.c_str(),HistName.c_str(),10,-1.0,1.0,10,0.0,2.0*TMath::Pi());
      //  h_mRcEffCosPhiPrime[i_step][i_cent][i_pt]->Sumw2();

      //  HistName = Form("h_mRc%dEffCosHPhiPrime_Cent_%d_Pt_%d",i_step,i_cent,i_pt);
      //  h_mRcEffCosHPhiPrime[i_step][i_cent][i_pt] = new TH2F(HistName.c_str(),HistName.c_str(),9,-1.0,1.0,12,0.0,2.0*TMath::Pi());
      //  h_mRcEffCosHPhiPrime[i_step][i_cent][i_pt]->Sumw2();
      //}
      //for(int i_y = 0; i_y < 10; i_y++)
      //{ 
      //  HistName = Form("h_mRc%dEffPhiS_Cent_%d_Y_%d",i_step,i_cent,i_y);
      //  h_mRcEffPhiSY[i_step][i_cent][i_y] = new TH1F(HistName.c_str(),HistName.c_str(),10,-1.0,1.0);
      //  h_mRcEffPhiSY[i_step][i_cent][i_y]->Sumw2();

      //  //HistName = Form("h_mRc%dEffCos_Cent_%d_Y_%d",i_step,i_cent,i_y);
      //  //h_mRcEffCosY[i_step][i_cent][i_y] = new TH1F(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2.,-1.0,1.0);
      //  //h_mRcEffCosY[i_step][i_cent][i_y]->Sumw2();

      //  //HistName = Form("h_mRc%dEffCosH_Cent_%d_Y_%d",i_step,i_cent,i_y);
      //  //h_mRcEffCosHY[i_step][i_cent][i_y] = new TH1F(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2.,-1.0,1.0);
      //  //h_mRcEffCosHY[i_step][i_cent][i_y]->Sumw2();

      //  //HistName = Form("h_mRc%dEffCosCosH_Cent_%d_Y_%d",i_step,i_cent,i_y);
      //  //h_mRcEffCosCosHY[i_step][i_cent][i_y] = new TH2F(HistName.c_str(),HistName.c_str(),9/*vmsa::BinCos*2.*/,-1.0,1.0,9/*vmsa::BinCos*2.*/,-1.0,1.0);
      //  //h_mRcEffCosCosHY[i_step][i_cent][i_y]->Sumw2();

      //  HistName = Form("h_mRc%dEffCosPhiPrime_Cent_%d_Y_%d",i_step,i_cent,i_y);
      //  h_mRcEffCosPhiPrimeY[i_step][i_cent][i_y] = new TH2F(HistName.c_str(),HistName.c_str(),10,-1.0,1.0,10,0.0,2.0*TMath::Pi());
      //  h_mRcEffCosPhiPrimeY[i_step][i_cent][i_y]->Sumw2();

      //  HistName = Form("h_mRc%dEffCosHPhiPrime_Cent_%d_Y_%d",i_step,i_cent,i_y);
      //  h_mRcEffCosHPhiPrimeY[i_step][i_cent][i_y] = new TH2F(HistName.c_str(),HistName.c_str(),9,-1.0,1.0,12,0.0,2.0*TMath::Pi());
      //  h_mRcEffCosHPhiPrimeY[i_step][i_cent][i_y]->Sumw2();
      //}
    }
  }
  h_FrameEta = new TH1F("h_FrameEta","h_FrameEta",2.0*vmsa::BinEta,-2.0*vmsa::mEtaMax,2.0*vmsa::mEtaMax);
  h_FramePhi = new TH1F("h_FramePhi","h_FramePhi",vmsa::BinPhi,0.0,2.0*TMath::Pi());
  h_FrameEtaKaon = new TH1F("h_FrameEtaKaon","h_FrameEtaKaon",vmsa::BinEta,vmsa::mEtaMax,vmsa::mEtaMax);
  h_FramePhiKaon = new TH1F("h_FramePhiKaon","h_FramePhiKaon",vmsa::BinPhi,0.0,2.0*TMath::Pi());
  flag_eff = 0;
  flag_eff_PtEtaPhi = 0;
}

float StEffHistManger::AngleShift(float phi)
{
  double const Psi2_low[3] = {-3.0*TMath::Pi()/2.0,-1.0*TMath::Pi()/2.0,1.0*TMath::Pi()/2.0};
  double const Psi2_up[3]  = {-1.0*TMath::Pi()/2.0, 1.0*TMath::Pi()/2.0,3.0*TMath::Pi()/2.0};

  float phi_shift = -999.0;
  for(int psi_bin = 0; psi_bin < 3; ++psi_bin)
  {
    if(phi >= Psi2_low[psi_bin] && phi < Psi2_up[psi_bin])
    {
      phi_shift = phi - (psi_bin-1)*2.0*TMath::Pi()/2.0;
    }
  }

  return phi_shift;
}

//void StEffHistManger::FillHistMc(double reweight, int cent, float pt, float eta, float y, float phi, float psi, float kppt, float kpeta, float kpy, float kpphi, float kmpt, float kmeta, float kmy, float kmphi, float phistar, float cos, float phiprime, float cosH, float phihelicity)
//{
//  reweight = 1.0;
//  if(phi < 0.0) phi += 2.0*TMath::Pi();
//  if(kpphi < 0.0) kpphi += 2.0*TMath::Pi();
//  if(kmphi < 0.0) kmphi += 2.0*TMath::Pi();
//
//  double phidiff = TMath::Cos(2.0*(phistar-phi));
//  if(pt < 1.2 || pt > 4.2) return;
//  double phipsi = phi-psi;
//  while(phipsi < 0.0) phipsi += 2.0*TMath::Pi();
//  while(phipsi > 2.0*TMath::Pi()) phipsi -= 2.0*TMath::Pi();
//  if(cent >= 2 && cent <= 5) 
//  {
//    double weight = h_mpTv2Weights[cent]->GetBinContent(h_mpTv2Weights[cent]->FindBin(pt,phipsi));
//    //h_mMcTracks[cent]->Fill(pt,eta,phi,weight);
//    //h_mMcTracks[9]->Fill(pt,eta,phi,weight);
//    //h_mMcTracksKplus[cent]->Fill(kppt,kpeta,kpphi,weight);
//    //h_mMcTracksKplus[9]->Fill(kppt,kpeta,kpphi,weight);
//    //h_mMcTracksKminus[cent]->Fill(kmpt,kmeta,kmphi,weight);
//    //h_mMcTracksKminus[9]->Fill(kmpt,kmeta,kmphi,weight);
//    h_mMcEffPtY[cent]->Fill(pt,y,weight);
//    h_mMcEffPtY[9]->Fill(pt,y,weight);
//    h_mKaon_MC[0][cent]->Fill(kppt,kpy,weight);
//    h_mKaon_MC[1][cent]->Fill(kmpt,kmy,weight);
//    h_mKaon_MC[0][9]->Fill(kppt,kpy,weight);
//    h_mKaon_MC[1][9]->Fill(kmpt,kmy,weight);
//
//    h_mMcEffPtPhiPsi[cent]->Fill(pt,phipsi,weight);
//    h_mMcEffPtPhiPsi[9]->Fill(pt,phipsi,weight);
//    for(int i_pt = 2; i_pt < 6; i_pt++)
//    {
//      if(pt >= vmsa::pt_low[4][i_pt] && pt < vmsa::pt_up[4][i_pt])
//      {
//        h_mMcEffPhiS[cent][i_pt]->Fill(phidiff,weight);
//        h_mMcEffPhiS[9][i_pt]->Fill(phidiff,weight);
//        h_mMcEffCos[cent][i_pt]->Fill(cos,weight);
//        h_mMcEffCos[9][i_pt]->Fill(cos,weight);
//        h_mMcEffCosH[cent][i_pt]->Fill(cosH,weight);
//        h_mMcEffCosH[9][i_pt]->Fill(cosH,weight);
//        h_mMcEffCosCosH[cent][i_pt]->Fill(cos,cosH,weight);
//        h_mMcEffCosCosH[9][i_pt]->Fill(cos,cosH,weight);
//        h_mMcEffCosPhiPrime[cent][i_pt]->Fill(cos,phiprime,weight);
//        h_mMcEffCosPhiPrime[9][i_pt]->Fill(cos,phiprime,weight);
//        h_mMcEffCosHPhiPrime[cent][i_pt]->Fill(cosH,phihelicity,weight);
//        h_mMcEffCosHPhiPrime[9][i_pt]->Fill(cosH,phihelicity,weight);
//      }
//    }
//  }
//  //else
//  //{
//  //  h_mMcTracks[cent]->Fill(pt,eta,phi);
//  //  h_mMcTracksKplus[cent]->Fill(kppt,kpeta,kpphi);
//  //  h_mMcTracksKminus[cent]->Fill(kmpt,kmeta,kmphi);
//  //  h_mMcEffPtY[cent]->Fill(pt,y);
//  //  h_mKaon_MC[0][cent]->Fill(kppt,kpy);
//  //  h_mKaon_MC[1][cent]->Fill(kmpt,kmy);
//
//  //  h_mMcEffPtPhiPsi[cent]->Fill(pt,phipsi);
//  //}
//  //if(eta > -1.0*vmsa::mEtaMax && eta < 0.0)
//  //
//  //{ // east eta cut
//  //  float phi_shift = AngleShift(phi-psiWest);
//  //  h_mMcTracks[cent]->Fill(pt,eta,phi_shift);
//  //  h_mMcTracks[9]->Fill(pt,eta,phi_shift);
//  //}
//  //if(eta >= 0.0 && eta < vmsa::mEtaMax)
//  //{ // west eta cut
//  //  float phi_shift = AngleShift(phi-psiEast);
//  //  h_mMcTracks[cent]->Fill(pt,eta,phi_shift);
//  //  h_mMcTracks[9]->Fill(pt,eta,phi_shift);
//  //}
//}
//
//void StEffHistManger::FillHistRc(double reweight, int cent, float pt, float eta, float y, float phi, float psi, float kppt, float kpeta, float kpy, float kpphi, float kmpt, float kmeta, float kmy, float kmphi, float mcpt, float mcphi, float phistar, float cos, float phiprime, float cosH, float phihelicity, float kp_m2, float km_m2, float kp_ns, float km_ns, int step)
//{
//  reweight = 1.0;
//
//  if(phi < 0.0) phi += 2.0*TMath::Pi();
//  if(mcphi < 0.0) mcphi += 2.0*TMath::Pi();
//  if(kpphi < 0.0) kpphi += 2.0*TMath::Pi();
//  if(kmphi < 0.0) kmphi += 2.0*TMath::Pi();
//
//  double phidiff = TMath::Cos(2.0*(phistar-phi));
//  if(pt < 1.2 || pt > 4.2) return;
//  if(mcpt < 1.2 || mcpt > 4.2) return;
//  double phipsi = mcphi-psi;
//  while(phipsi < 0.0) phipsi += 2.0*TMath::Pi();
//  while(phipsi > 2.0*TMath::Pi()) phipsi -= 2.0*TMath::Pi();
//  if(cent >= 2 && cent <= 5) 
//  {
//    double weight = h_mpTv2Weights[cent]->GetBinContent(h_mpTv2Weights[cent]->FindBin(mcpt,phipsi));
//    //h_mRcTracks[step][cent]->Fill(pt,eta,phi,weight);
//    //h_mRcTracks[step][9]->Fill(pt,eta,phi,weight);
//    //h_mRcTracksKplus[step][cent]->Fill(kppt,kpeta,kpphi,weight);
//    //h_mRcTracksKplus[step][9]->Fill(kppt,kpeta,kpphi,weight);
//    //h_mRcTracksKminus[step][cent]->Fill(kmpt,kmeta,kmphi,weight);
//    //h_mRcTracksKminus[step][9]->Fill(kmpt,kmeta,kmphi,weight);
//    h_mRcEffPtY[step][cent]->Fill(pt,y,weight*reweight);
//    h_mRcEffPtY[step][9]->Fill(pt,y,weight*reweight);
//    h_mKaon_RC[step][0][cent]->Fill(kppt,kpy,weight*reweight);
//    h_mKaon_RC[step][1][cent]->Fill(kmpt,kmy,weight*reweight);
//    h_mKaon_RC[step][0][9]->Fill(kppt,kpy,weight*reweight);
//    h_mKaon_RC[step][1][9]->Fill(kmpt,kmy,weight*reweight);
//
//    h_mKaon_RC_m2[step][0][cent]->Fill(kppt,kp_m2,kp_ns,weight*reweight);
//    h_mKaon_RC_m2[step][1][cent]->Fill(kmpt,km_m2,km_ns,weight*reweight);
//    h_mKaon_RC_m2[step][0][9]->Fill(kppt,kp_m2,kp_ns,weight*reweight);
//    h_mKaon_RC_m2[step][1][9]->Fill(kmpt,km_m2,km_ns,weight*reweight);
//
//    for(int i_pt = 2; i_pt < 6; i_pt++)
//    {
//      if(pt >= vmsa::pt_low[4][i_pt] && pt < vmsa::pt_up[4][i_pt])
//      {
//        h_mRcEffPhiS[step][cent][i_pt]->Fill(phidiff,weight*reweight);
//        h_mRcEffPhiS[step][9][i_pt]->Fill(phidiff,weight*reweight);
//        h_mRcEffCos[step][cent][i_pt]->Fill(cos,weight*reweight);
//        h_mRcEffCos[step][9][i_pt]->Fill(cos,weight*reweight);
//        h_mRcEffCosH[step][cent][i_pt]->Fill(cosH,weight*reweight);
//        h_mRcEffCosH[step][9][i_pt]->Fill(cosH,weight*reweight);
//        h_mRcEffCosCosH[step][cent][i_pt]->Fill(cos,cosH,weight*reweight);
//        h_mRcEffCosCosH[step][9][i_pt]->Fill(cos,cosH,weight*reweight);
//        h_mRcEffCosPhiPrime[step][cent][i_pt]->Fill(cos,phiprime,weight*reweight);
//        h_mRcEffCosPhiPrime[step][9][i_pt]->Fill(cos,phiprime,weight*reweight);
//        h_mRcEffCosHPhiPrime[step][cent][i_pt]->Fill(cosH,phihelicity,weight*reweight);
//        h_mRcEffCosHPhiPrime[step][9][i_pt]->Fill(cosH,phihelicity,weight*reweight);
//      }
//    }
//  }
//  //else
//  //{
//  //  h_mRcTracks[step][cent]->Fill(pt,eta,phi);
//  //  h_mRcTracksKplus[step][cent]->Fill(kppt,kpeta,kpphi);
//  //  h_mRcTracksKminus[step][cent]->Fill(kmpt,kmeta,kmphi);
//  //  h_mRcEffPtY[step][cent]->Fill(pt,y);
//  //  h_mKaon_RC[step][0][cent]->Fill(kppt,kpy);
//  //  h_mKaon_RC[step][1][cent]->Fill(kmpt,kmy);
//  //}
//  //if(eta > -1.0*vmsa::mEtaMax && eta < 0.0)
//  //{ // east eta cut
//  //  float phi_shift = AngleShift(phi-psiWest);
//  //  h_mRcTracks[cent]->Fill(pt,eta,phi_shift);
//  //  h_mRcTracks[9]->Fill(pt,eta,phi_shift);
//  //}
//  //if(eta >= 0.0 && eta < vmsa::mEtaMax)
//  //{ // west eta cut
//  //  float phi_shift = AngleShift(phi-psiEast);
//  //  h_mRcTracks[cent]->Fill(pt,eta,phi_shift);
//  //  h_mRcTracks[9]->Fill(pt,eta,phi_shift);
//  //}
//}

void StEffHistManger::FillHistMc(double reweight, int cent, float pt, float eta, float y, float phi, float psi, float kppt, float kpeta, float kpy, float kpphi, float kmpt, float kmeta, float kmy, float kmphi, float phistar, float cos, float phiprime, float cosH, float phihelicity)
{
  reweight = 1.0;
  if(phi < 0.0) phi += 2.0*TMath::Pi();
  if(kpphi < 0.0) kpphi += 2.0*TMath::Pi();
  if(kmphi < 0.0) kmphi += 2.0*TMath::Pi();

  double phidiff = TMath::Cos(2.0*(phistar-phi));
  if(pt < 1.2 || pt > 4.2) return;
  double phipsi = phi-psi;
  while(phipsi < 0.0) phipsi += 2.0*TMath::Pi();
  while(phipsi > 2.0*TMath::Pi()) phipsi -= 2.0*TMath::Pi();

  if(cent >= 2 && cent <= 5) 
  {
      double weight  = 1.0;//h_mpTv2Weights[cent]->GetBinContent(h_mpTv2Weights[cent]->FindBin(pt,phipsi));
      //h_mMcTracks[cent]->Fill(pt,eta,phi,weight);
      //h_mMcTracks[9]->Fill(pt,eta,phi,weight);
      //h_mMcTracksKplus[cent]->Fill(kppt,kpeta,kpphi,weight);
      //h_mMcTracksKplus[9]->Fill(kppt,kpeta,kpphi,weight);
      //h_mMcTracksKminus[cent]->Fill(kmpt,kmeta,kmphi,weight);
      //h_mMcTracksKminus[9]->Fill(kmpt,kmeta,kmphi,weight);
      h_mMcEffPtY[cent]->Fill(pt,y,weight);
      h_mMcEffPtY[9]->Fill(pt,y,weight);
      h_mKaon_MC[0][cent]->Fill(kppt,kpy,weight);
      h_mKaon_MC[1][cent]->Fill(kmpt,kmy,weight);
      h_mKaon_MC[0][9]->Fill(kppt,kpy,weight);
      h_mKaon_MC[1][9]->Fill(kmpt,kmy,weight);

      //h_mMcEffPtPhiPsi[cent]->Fill(pt,phipsi,weight);
      //h_mMcEffPtPhiPsi[9]->Fill(pt,phipsi,weight);
      for(int i_pt = 2; i_pt < 6; i_pt++)
      {
        if(pt >= vmsa::pt_low[4][i_pt] && pt < vmsa::pt_up[4][i_pt])
        {
          h_mMcEffPhiS[cent][i_pt]->Fill(phidiff,weight);
          h_mMcEffPhiS[9][i_pt]->Fill(phidiff,weight);
          h_mMcEffCos[cent][i_pt]->Fill(cos,weight);
          h_mMcEffCos[9][i_pt]->Fill(cos,weight);
          h_mMcEffCosH[cent][i_pt]->Fill(cosH,weight);
          h_mMcEffCosH[9][i_pt]->Fill(cosH,weight);
          h_mMcEffCosCosH[cent][i_pt]->Fill(cos,cosH,weight);
          h_mMcEffCosCosH[9][i_pt]->Fill(cos,cosH,weight);
          h_mMcEffCosPhiPrime[cent][i_pt]->Fill(cos,phiprime,weight);
          h_mMcEffCosPhiPrime[9][i_pt]->Fill(cos,phiprime,weight);
          h_mMcEffCosHPhiPrime[cent][i_pt]->Fill(cosH,phihelicity,weight);
          h_mMcEffCosHPhiPrime[9][i_pt]->Fill(cosH,phihelicity,weight);
        }
      }
  }
  //else
  //{
  //  h_mMcTracks[cent]->Fill(pt,eta,phi);
  //  h_mMcTracksKplus[cent]->Fill(kppt,kpeta,kpphi);
  //  h_mMcTracksKminus[cent]->Fill(kmpt,kmeta,kmphi);
  //  h_mMcEffPtY[cent]->Fill(pt,y);
  //  h_mKaon_MC[0][cent]->Fill(kppt,kpy);
  //  h_mKaon_MC[1][cent]->Fill(kmpt,kmy);

  //  h_mMcEffPtPhiPsi[cent]->Fill(pt,phipsi);
  //}
  //if(eta > -1.0*vmsa::mEtaMax && eta < 0.0)
  //
  //{ // east eta cut
  //  float phi_shift = AngleShift(phi-psiWest);
  //  h_mMcTracks[cent]->Fill(pt,eta,phi_shift);
  //  h_mMcTracks[9]->Fill(pt,eta,phi_shift);
  //}
  //if(eta >= 0.0 && eta < vmsa::mEtaMax)
  //{ // west eta cut
  //  float phi_shift = AngleShift(phi-psiEast);
  //  h_mMcTracks[cent]->Fill(pt,eta,phi_shift);
  //  h_mMcTracks[9]->Fill(pt,eta,phi_shift);
  //}
}

void StEffHistManger::FillHistRc(double reweight, int cent, float pt, float eta, float y, float phi, float psi, float kppt, float kpeta, float kpy, float kpphi, float kmpt, float kmeta, float kmy, float kmphi, float mcpt, float mcphi, float mcy, float phistar, float cos, float phiprime, float cosH, float phihelicity, float kp_m2, float km_m2, float kp_ns, float km_ns, int step)
{
  //reweight = 1.0;

  if(phi < 0.0) phi += 2.0*TMath::Pi();
  if(mcphi < 0.0) mcphi += 2.0*TMath::Pi();
  if(kpphi < 0.0) kpphi += 2.0*TMath::Pi();
  if(kmphi < 0.0) kmphi += 2.0*TMath::Pi();

  double phidiff = TMath::Cos(2.0*(phistar-phi));
  //if(pt < 1.2 || pt > 4.2) return;
  //if(mcpt < 1.2 || mcpt > 4.2) return;
  double phipsi = mcphi-psi;
  while(phipsi < 0.0) phipsi += 2.0*TMath::Pi();
  while(phipsi > 2.0*TMath::Pi()) phipsi -= 2.0*TMath::Pi();
  //if(cent >= 2 && cent <= 5) 
  //{
      double weight = 1.0;//h_mpTyv2Weights[cent]->GetBinContent(h_mpTyv2Weights[cent]->FindBin(mcpt,mcy,phipsi));
      //h_mRcTracks[step][cent]->Fill(pt,eta,phi,weight);
      //h_mRcTracks[step][9]->Fill(pt,eta,phi,weight);
      //h_mRcTracksKplus[step][cent]->Fill(kppt,kpeta,kpphi,weight);
      //h_mRcTracksKplus[step][9]->Fill(kppt,kpeta,kpphi,weight);
      //h_mRcTracksKminus[step][cent]->Fill(kmpt,kmeta,kmphi,weight);
      //h_mRcTracksKminus[step][9]->Fill(kmpt,kmeta,kmphi,weight);
      h_mRcEffPtY[step][cent]->Fill(pt,y,weight*reweight);
      h_mRcEffPtY[step][9]->Fill(pt,y,weight*reweight);
      h_mKaon_RC[step][0][cent]->Fill(kppt,kpy,weight*reweight);
      h_mKaon_RC[step][1][cent]->Fill(kmpt,kmy,weight*reweight);
      h_mKaon_RC[step][0][9]->Fill(kppt,kpy,weight*reweight);
      h_mKaon_RC[step][1][9]->Fill(kmpt,kmy,weight*reweight);

      //h_mKaon_RC_m2[step][0][cent]->Fill(kppt,kp_m2,kp_ns,weight*reweight);
      //h_mKaon_RC_m2[step][1][cent]->Fill(kmpt,km_m2,km_ns,weight*reweight);
      //h_mKaon_RC_m2[step][0][9]->Fill(kppt,kp_m2,kp_ns,weight*reweight);
      //h_mKaon_RC_m2[step][1][9]->Fill(kmpt,km_m2,km_ns,weight*reweight);

      //h_mMcEffPtPhiPsi[step][cent]->Fill(pt,phipsi,reweight);
      //h_mMcEffPtPhiPsi[step][9]->Fill(pt,phipsi,reweight);
      h_mMcEffPtYPhiPsi[step][cent]->Fill(pt,y,phipsi,reweight*weight);
      h_mMcEffPtYPhiPsi[step][9]->Fill(pt,y,phipsi,reweight*weight);
      //for(int i_pt = 2; i_pt < 6; i_pt++)
      //{
      //  if(pt >= vmsa::pt_low[4][i_pt] && pt < vmsa::pt_up[4][i_pt])
      //  {
      //    h_mRcEffPhiS[step][cent][i_pt]->Fill(phidiff,weight*reweight);
      //    h_mRcEffPhiS[step][9][i_pt]->Fill(phidiff,weight*reweight);
      //    //h_mRcEffCos[step][cent][i_pt]->Fill(cos,weight*reweight);
      //    //h_mRcEffCos[step][9][i_pt]->Fill(cos,weight*reweight);
      //    //h_mRcEffCosH[step][cent][i_pt]->Fill(cosH,weight*reweight);
      //    //h_mRcEffCosH[step][9][i_pt]->Fill(cosH,weight*reweight);
      //    //h_mRcEffCosCosH[step][cent][i_pt]->Fill(cos,cosH,weight*reweight);
      //    //h_mRcEffCosCosH[step][9][i_pt]->Fill(cos,cosH,weight*reweight);
      //    h_mRcEffCosPhiPrime[step][cent][i_pt]->Fill(cos,phiprime,weight*reweight);
      //    h_mRcEffCosPhiPrime[step][9][i_pt]->Fill(cos,phiprime,weight*reweight);
      //    h_mRcEffCosHPhiPrime[step][cent][i_pt]->Fill(cosH,phihelicity,weight*reweight);
      //    h_mRcEffCosHPhiPrime[step][9][i_pt]->Fill(cosH,phihelicity,weight*reweight);
      //  }
      //}
      //for(int i_y = 0; i_y < 10; i_y++)
      //{
      //  if(y >= (float(i_y)-5.)/5. && y < (float(i_y)-4.)/5.)
      //  {
      //    h_mRcEffPhiSY[step][cent][i_y]->Fill(phidiff,weight*reweight);
      //    h_mRcEffPhiSY[step][9][i_y]->Fill(phidiff,weight*reweight);
      //    //h_mRcEffCosY[step][cent][i_y]->Fill(cos,weight*reweight);
      //    //h_mRcEffCosY[step][9][i_y]->Fill(cos,weight*reweight);
      //    //h_mRcEffCosHY[step][cent][i_y]->Fill(cosH,weight*reweight);
      //    //h_mRcEffCosHY[step][9][i_y]->Fill(cosH,weight*reweight);
      //    //h_mRcEffCosCosHY[step][cent][i_y]->Fill(cos,cosH,weight*reweight);
      //    //h_mRcEffCosCosHY[step][9][i_y]->Fill(cos,cosH,weight*reweight);
      //    h_mRcEffCosPhiPrimeY[step][cent][i_y]->Fill(cos,phiprime,weight*reweight);
      //    h_mRcEffCosPhiPrimeY[step][9][i_y]->Fill(cos,phiprime,weight*reweight);
      //    h_mRcEffCosHPhiPrimeY[step][cent][i_y]->Fill(cosH,phihelicity,weight*reweight);
      //    h_mRcEffCosHPhiPrimeY[step][9][i_y]->Fill(cosH,phihelicity,weight*reweight);
      //  }
      //}
  //}
  //else
  //{
  //  h_mRcTracks[step][cent]->Fill(pt,eta,phi);
  //  h_mRcTracksKplus[step][cent]->Fill(kppt,kpeta,kpphi);
  //  h_mRcTracksKminus[step][cent]->Fill(kmpt,kmeta,kmphi);
  //  h_mRcEffPtY[step][cent]->Fill(pt,y);
  //  h_mKaon_RC[step][0][cent]->Fill(kppt,kpy);
  //  h_mKaon_RC[step][1][cent]->Fill(kmpt,kmy);
  //}
  //if(eta > -1.0*vmsa::mEtaMax && eta < 0.0)
  //{ // east eta cut
  //  float phi_shift = AngleShift(phi-psiWest);
  //  h_mRcTracks[cent]->Fill(pt,eta,phi_shift);
  //  h_mRcTracks[9]->Fill(pt,eta,phi_shift);
  //}
  //if(eta >= 0.0 && eta < vmsa::mEtaMax)
  //{ // west eta cut
  //  float phi_shift = AngleShift(phi-psiEast);
  //  h_mRcTracks[cent]->Fill(pt,eta,phi_shift);
  //  h_mRcTracks[9]->Fill(pt,eta,phi_shift);
  //}
}

void StEffHistManger::FillHistPt(int cent, float McPt, float gRcPt, float pRcPt)
{
  h_mPtGl[cent]->Fill(McPt,gRcPt-McPt);
  if(cent >= 2 && cent <= 5) h_mPtGl[9]->Fill(McPt,gRcPt-McPt);
  h_mPtPr[cent]->Fill(McPt,pRcPt-McPt);
  if(cent >= 2 && cent <= 5) h_mPtPr[9]->Fill(McPt,pRcPt-McPt);
}

TH1F* StEffHistManger::CalEffError(TH1F *h_Mc, TH1F *h_Rc, std::string HistName)
{
  TH1F* h_ratio = (TH1F*)h_Rc->Clone();
  //h_ratio->Divide(h_Mc);
  h_ratio->Divide(h_Rc,h_Mc,1,1,"B");
  //for(int i_bin = 1; i_bin < h_ratio->GetNbinsX()+1; ++i_bin)
  //{
  //  double n = h_Mc->GetBinContent(i_bin);
  //  double k = h_Rc->GetBinContent(i_bin);
  //  double variance = (k+1.0)*(k+2.0)/((n+2.0)*(n+3.0))-(k+1.0)*(k+1.0)/((n+2.0)*(n+2.0));
  //  double sigma = TMath::Sqrt(variance);
  //  if(n > 0.0 && k > 0.0) h_ratio->SetBinError(i_bin,sigma);
  //}
  h_ratio->SetName(HistName.c_str());

  return h_ratio;
}

//void StEffHistManger::CalEfficiency()
//{
//  for(int i_cent = 0; i_cent < 10; ++i_cent)
//  {
//    std::string HistName;
//
//    HistName = Form("h_mMcEffPt_Cent_%d",i_cent);
//    h_mMcEffPt[i_cent] = (TH1F*)h_mMcTracks[i_cent]->Project3D("x")->Clone(HistName.c_str());
//    HistName = Form("h_mRcEffPt_Cent_%d",i_cent);
//    h_mRcEffPt[i_cent] = (TH1F*)h_mRcTracks[i_cent]->Project3D("x")->Clone(HistName.c_str());
//    HistName = Form("h_mEffPt_Cent_%d",i_cent);
//    h_mEffPt[i_cent] = CalEffError(h_mMcEffPt[i_cent],h_mRcEffPt[i_cent],HistName.c_str());
//
//    HistName = Form("h_mMcEffEta_Cent_%d",i_cent);
//    h_mMcEffEta[i_cent] = (TH1F*)h_mMcTracks[i_cent]->Project3D("y")->Clone(HistName.c_str());
//    HistName = Form("h_mRcEffEta_Cent_%d",i_cent);
//    h_mRcEffEta[i_cent] = (TH1F*)h_mRcTracks[i_cent]->Project3D("y")->Clone(HistName.c_str());
//    HistName = Form("h_mEffEta_Cent_%d",i_cent);
//    h_mEffEta[i_cent] = CalEffError(h_mMcEffEta[i_cent],h_mRcEffEta[i_cent],HistName.c_str());
//
//    HistName = Form("h_mMcEffPhi_Cent_%d",i_cent);
//    h_mMcEffPhi[i_cent] = (TH1F*)h_mMcTracks[i_cent]->Project3D("z")->Clone(HistName.c_str());
//    HistName = Form("h_mRcEffPhi_Cent_%d",i_cent);
//    h_mRcEffPhi[i_cent] = (TH1F*)h_mRcTracks[i_cent]->Project3D("z")->Clone(HistName.c_str());
//    HistName = Form("h_mEffPhi_Cent_%d",i_cent);
//    h_mEffPhi[i_cent] = CalEffError(h_mMcEffPhi[i_cent],h_mRcEffPhi[i_cent],HistName.c_str());
//  }
//  flag_eff = 1;
//}
//
//void StEffHistManger::CalEffPtEtaPhi()
//{
//  for(int i_cent = 0; i_cent < 10; ++i_cent)
//  {
//    for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
//    {
//      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
//      {
//	std::string HistNameMc = Form("h_mMcEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
//	h_mMcEffPEP[HistNameMc] = (TH1F*)h_mMcTracks[i_cent]->ProjectionX(HistNameMc.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1);
//
//	std::string HistNameRc = Form("h_mRcEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
//	h_mRcEffPEP[HistNameRc] = (TH1F*)h_mRcTracks[i_cent]->ProjectionX(HistNameRc.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1);
//
//	std::string HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
//	h_mEffPEP[HistNameEff] = CalEffError(h_mMcEffPEP[HistNameMc],h_mRcEffPEP[HistNameRc],HistNameEff.c_str());
//      }
//    }
//  }
//  flag_eff_PtEtaPhi = 1;
//}

void StEffHistManger::WriteHist()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    //h_mMcTracks[i_cent]->Write();
    //h_mMcTracksKplus[i_cent]->Write();
    //h_mMcTracksKminus[i_cent]->Write();
    //h_mMcEffPtY[i_cent]->Write();

    //h_mMcEffPtPhiPsi[i_cent]->Write();
    //
    //h_mKaon_MC[0][i_cent]->Write();
    //h_mKaon_MC[1][i_cent]->Write();

    //if(i_cent >= 2 && i_cent<=5)
    //{
    //  for(int i_pt = 2; i_pt < 6; i_pt++)
    //  {
    //    h_mMcEffPhiS[i_cent][i_pt]->Write();
    //    h_mMcEffCos[i_cent][i_pt]->Write();
    //    h_mMcEffCosH[i_cent][i_pt]->Write();
    //    h_mMcEffCosCosH[i_cent][i_pt]->Write();
    //    h_mMcEffCosPhiPrime[i_cent][i_pt]->Write();
    //    h_mMcEffCosHPhiPrime[i_cent][i_pt]->Write();
    //    if(i_cent == 2)
    //    {
    //      h_mMcEffPhiS[9][i_pt]->Write();
    //      h_mMcEffCos[9][i_pt]->Write();
    //      h_mMcEffCosH[9][i_pt]->Write();
    //      h_mMcEffCosCosH[9][i_pt]->Write();
    //      h_mMcEffCosPhiPrime[9][i_pt]->Write();
    //      h_mMcEffCosHPhiPrime[9][i_pt]->Write();
    //    }
    //  }
    //}
  }
  for(int i_step = 0; i_step < 10; ++i_step)
  {
    for(int i_cent = 0; i_cent < 10; ++i_cent)
    {
      //h_mRcTracks[i_step][i_cent]->Write();
      //h_mRcTracksKplus[i_step][i_cent]->Write();
      //h_mRcTracksKminus[i_step][i_cent]->Write();
      h_mRcEffPtY[i_step][i_cent]->Write();

      h_mKaon_RC[i_step][0][i_cent]->Write();
      h_mKaon_RC[i_step][1][i_cent]->Write();
      //h_mKaon_RC_m2[i_step][0][i_cent]->Write();
      //h_mKaon_RC_m2[i_step][1][i_cent]->Write();

      //h_mMcEffPtPhiPsi[i_step][i_cent]->Write();
      h_mMcEffPtYPhiPsi[i_step][i_cent]->Write();
      //if(i_cent >= 2 && i_cent<=5)
      //{
        //for(int i_pt = 2; i_pt < 6; i_pt++)
        //{
        //  //std::string HistName = Form("h_m%dEffPhiS_Cent_%d_Pt_%d",i_step,i_cent,i_pt);
        //  //h_mEffPhiS[i_step][i_cent][i_pt] = CalEffError(h_mMcEffPhiS[i_cent][i_pt],h_mRcEffPhiS[i_step][i_cent][i_pt],HistName.c_str());
        //  h_mRcEffPhiS[i_step][i_cent][i_pt]->Write();
        //  //h_mRcEffCos[i_step][i_cent][i_pt]->Write();
        //  //h_mRcEffCosH[i_step][i_cent][i_pt]->Write();
        //  //h_mRcEffCosCosH[i_step][i_cent][i_pt]->Write();
        //  h_mRcEffCosPhiPrime[i_step][i_cent][i_pt]->Write();
        //  h_mRcEffCosHPhiPrime[i_step][i_cent][i_pt]->Write();
        //  //if(i_cent == 2)
        //  //{
        //  //  //std::string HistName = Form("h_m%dEffPhiS_Cent_%d_Pt_%d",i_step,9,i_pt);
        //  //  //h_mEffPhiS[i_step][9][i_pt] = CalEffError(h_mMcEffPhiS[9][i_pt],h_mRcEffPhiS[i_step][9][i_pt],HistName.c_str());
        //  //  h_mRcEffPhiS[i_step][9][i_pt]->Write();
        //  //  //h_mRcEffCos[i_step][9][i_pt]->Write();
        //  //  //h_mRcEffCosH[i_step][9][i_pt]->Write();
        //  //  //h_mRcEffCosCosH[i_step][9][i_pt]->Write();
        //  //  h_mRcEffCosPhiPrime[i_step][9][i_pt]->Write();
        //  //  h_mRcEffCosHPhiPrime[i_step][9][i_pt]->Write();
        //  //}
        //}
        //for(int i_y = 0; i_y < 10; i_y++)
        //{
        //  //std::string HistName = Form("h_m%dEffPhiS_Cent_%d_Pt_%d",i_step,i_cent,i_y);
        //  //h_mEffPhiS[i_step][i_cent][i_y] = CalEffError(h_mMcEffPhiS[i_cent][i_y],h_mRcEffPhiS[i_step][i_cent][i_y],HistName.c_str());
        //  h_mRcEffPhiSY[i_step][i_cent][i_y]->Write();
        //  //h_mRcEffCosY[i_step][i_cent][i_y]->Write();
        //  //h_mRcEffCosHY[i_step][i_cent][i_y]->Write();
        //  //h_mRcEffCosCosHY[i_step][i_cent][i_y]->Write();
        //  h_mRcEffCosPhiPrimeY[i_step][i_cent][i_y]->Write();
        //  h_mRcEffCosHPhiPrimeY[i_step][i_cent][i_y]->Write();
        //  ///if(i_cent == 2)
        //  ///{
        //  ///  //std::string HistName = Form("h_m%dEffPhiS_Cent_%d_Pt_%d",i_step,9,i_y);
        //  ///  //h_mEffPhiS[i_step][9][i_y] = CalEffError(h_mMcEffPhiS[9][i_y],h_mRcEffPhiS[i_step][9][i_y],HistName.c_str());
        //  ///  h_mRcEffPhiSY[i_step][9][i_y]->Write();
        //  ///  h_mRcEffCosY[i_step][9][i_y]->Write();
        //  ///  h_mRcEffCosHY[i_step][9][i_y]->Write();
        //  ///  h_mRcEffCosCosHY[i_step][9][i_y]->Write();
        //  ///  h_mRcEffCosPhiPrimeY[i_step][9][i_y]->Write();
        //  ///  h_mRcEffCosHPhiPrimeY[i_step][9][i_y]->Write();
        //  ///}
        //}
      //}
    }
  }
  h_FrameEta->Write();
  h_FramePhi->Write();
  h_FrameEtaKaon->Write();
  h_FramePhiKaon->Write();
}
