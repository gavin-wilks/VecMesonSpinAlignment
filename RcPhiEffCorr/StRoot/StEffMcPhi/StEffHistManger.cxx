#include "StEffHistManger.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3D.h"
#include "TMath.h"
#include <iostream>
#include "TLorentzVector.h"
#include "TVector3.h"
#include <iostream>
#include "TFile.h"

using namespace std;

ClassImp(StEffHistManger)
//
StEffHistManger::StEffHistManger(int energy, int pid, int mode, int startpt, int stoppt)
{
  mEnergy = energy;
  mMode = mode;
  mStartPt = startpt;
  mStopPt = stoppt;
  //mOrder = order;
  
  //string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/ptyspectra_datarcratio_%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  //TFile *File_Spectra = TFile::Open(inputfile.c_str());
  //h_mSpectraRatio = (TH2D*) ((TH2D*) File_Spectra->Get("pty_datarcratio"))->Clone();
 
  if( pid == 0 && (mode == 0 || mode == 3))
  {
    mpt_first = mStartPt;//vmsa::pt_rebin_first[mEnergy];
    mpt_last  = mStopPt;//vmsa::pt_rebin_last[mEnergy];

    for(int i = mpt_first; i <= mpt_last; i++)
    {
      mpt_low[i] = vmsa::pt_low[mEnergy][i];
      mpt_up[i]  = vmsa::pt_up[mEnergy][i];
    }
  }

  if( pid == 0 && mode == 1 )
  {
    mpt_first = mStartPt;
    mpt_last  = mStopPt;

    for(int i = mpt_first; i <= mpt_last; i++)
    {
      mpt_low[i] = vmsa::pt_low_cent[mEnergy][i];
      mpt_up[i]  = vmsa::pt_up_cent[mEnergy][i];
    }
  }

  if( pid == 0 && mode == 2 )
  {
    mpt_first = mStartPt;
    mpt_last  = mStopPt;

    for(int i = mpt_first; i <= mpt_last; i++)
    {
      mpt_low[i] = vmsa::pt_low_y[mEnergy][i];
      mpt_up[i]  = vmsa::pt_up_y[mEnergy][i];
    }
  }

  if( pid == 2 )
  {
    mpt_first = vmsa::pt_rebin_firstKS[mEnergy];
    mpt_last  = vmsa::pt_rebin_lastKS[mEnergy];

    for(int i = mpt_first; i < mpt_last; i++)
    {
      mpt_low[i] = vmsa::pt_lowKS[mEnergy][i];
      mpt_up[i]  = vmsa::pt_upKS[mEnergy][i];
    }
  }
  
  std::cout << "mpt_first = " << mpt_first << std::endl;
  std::cout << "mpt_last = " << mpt_last << std::endl;
}

StEffHistManger::~StEffHistManger()
{
  /* */
}

void StEffHistManger::InitPhiHist()
{
  for(Int_t i_pt_phi = mpt_first-2; i_pt_phi <= mpt_last-2; ++i_pt_phi) // use rebinned pt
  {
    for(Int_t i_pt = data_constants::phi_pt_start[i_pt_phi]; i_pt <= data_constants::phi_pt_stop[i_pt_phi]; i_pt++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++)
      {
        TString KEY_MC = Form("mcphi_phipt_%d_cos2phistarphi_%d",i_pt,i_phi);
        h_mPhi_MC[KEY_MC] = new TH2F(KEY_MC.Data(),KEY_MC.Data(),data_constants::kaon_azimuthal_bins*6,-TMath::Pi(),TMath::Pi(),data_constants::kaon_rapidity_bins*3,-1,1);
        h_mPhi_MC[KEY_MC]->Sumw2();
        //std::cout << "INITIALIZED: " << KEY_MC << std::endl;
        TString KEY_RC = Form("rcphi_phipt_%d_cos2phistarphi_%d",i_pt,i_phi);
        h_mPhi_RC[KEY_RC] = new TH2F(KEY_RC.Data(),KEY_RC.Data(),data_constants::kaon_azimuthal_bins*6,-TMath::Pi(),TMath::Pi(),data_constants::kaon_rapidity_bins*3,-1,1);
        h_mPhi_RC[KEY_RC]->Sumw2();
        //std::cout << "INITIALIZED: " << KEY_RC << std::endl;
      }

      TString KEY_MC = Form("mcphi_phistarphi_phipt_%d",i_pt);
      h_mPhi_MC[KEY_MC] = new TH2F(KEY_MC.Data(),KEY_MC.Data(),data_constants::kaon_azimuthal_bins*6,-TMath::Pi(),TMath::Pi(),data_constants::kaon_azimuthal_bins*6,-TMath::Pi(),TMath::Pi());
      h_mPhi_MC[KEY_MC]->Sumw2();
      //std::cout << "INITIALIZED: " << KEY_MC << std::endl;
      TString KEY_RC = Form("rcphi_phistarphi_phipt_%d",i_pt);
      h_mPhi_RC[KEY_RC] = new TH2F(KEY_RC.Data(),KEY_RC.Data(),data_constants::kaon_azimuthal_bins*6,-TMath::Pi(),TMath::Pi(),data_constants::kaon_azimuthal_bins*6,-TMath::Pi(),TMath::Pi());
      h_mPhi_RC[KEY_RC]->Sumw2();
      //std::cout << "INITIALIZED: " << KEY_RC << std::endl;

      KEY_MC = Form("mcphi_phistarmphi_phipt_%d",i_pt);
      h_mPhi1D_MC[KEY_MC] = new TH1F(KEY_MC.Data(),KEY_MC.Data(),data_constants::kaon_azimuthal_bins*6,-2*TMath::Pi(),2*TMath::Pi());
      h_mPhi1D_MC[KEY_MC]->Sumw2();
      //std::cout << "INITIALIZED: " << KEY_MC << std::endl;
      KEY_RC = Form("rcphi_phistarmphi_phipt_%d",i_pt);
      h_mPhi1D_RC[KEY_RC] = new TH1F(KEY_RC.Data(),KEY_RC.Data(),data_constants::kaon_azimuthal_bins*6,-2*TMath::Pi(),2*TMath::Pi());
      h_mPhi1D_RC[KEY_RC]->Sumw2();
      //std::cout << "INITIALIZED: " << KEY_RC << std::endl;
    }
  }
}

void StEffHistManger::InitKaonHist()
{
  for(Int_t i_pt_phi = mpt_first-2; i_pt_phi <= mpt_last-2; ++i_pt_phi) // use rebinned pt
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++)
    {
      for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
      {
        TString KEY_MC = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        h_mKaon_MC[KEY_MC] = new TH2F(KEY_MC.Data(),KEY_MC.Data(),data_constants::kaon_azimuthal_bins*10,-TMath::Pi(),TMath::Pi(),data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
        h_mKaon_MC[KEY_MC]->Sumw2();                                                                                                                                                                    
        //std::cout << "INITIALIZED: " << KEY_MC << std::endl;                                                                                                                                            
        TString KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);                                                                                                  
        h_mKaon_RC[KEY_RC] = new TH2F(KEY_RC.Data(),KEY_RC.Data(),data_constants::kaon_azimuthal_bins*10,-TMath::Pi(),TMath::Pi(),data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
        h_mKaon_RC[KEY_RC]->Sumw2();
        //std::cout << "INITIALIZED: " << KEY_RC << std::endl;

        KEY_MC = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);
        h_mKaon_MC[KEY_MC] = new TH2F(KEY_MC.Data(),KEY_MC.Data(),data_constants::kaon_azimuthal_bins*10,-TMath::Pi(),TMath::Pi(),data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
        h_mKaon_MC[KEY_MC]->Sumw2();                                                                                                                                                                    
        //std::cout << "INITIALIZED: " << KEY_MC << std::endl;                                                                                                                                            
        KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);                                                                                                  
        h_mKaon_RC[KEY_RC] = new TH2F(KEY_RC.Data(),KEY_RC.Data(),data_constants::kaon_azimuthal_bins*10,-TMath::Pi(),TMath::Pi(),data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
        h_mKaon_RC[KEY_RC]->Sumw2();
        //std::cout << "INITIALIZED: " << KEY_RC << std::endl;
      }
      TString KEY_MC1 = Form("mckaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_MC[KEY_MC1] = new TH2F(KEY_MC1.Data(),KEY_MC1.Data(),data_constants::kaon_azimuthal_bins*6,-0.6,0.6,data_constants::kaon_rapidity_bins*3,-0.6,0.6);
      h_mKaon_MC[KEY_MC1]->Sumw2();
      //std::cout << "INITIALIZED: " << KEY_MC1 << std::endl;
      TString KEY_RC1 = Form("rckaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_RC[KEY_RC1] = new TH2F(KEY_RC1.Data(),KEY_RC1.Data(),data_constants::kaon_azimuthal_bins*6,-0.6,0.6,data_constants::kaon_rapidity_bins*3,-0.6,0.6);
      h_mKaon_RC[KEY_RC1]->Sumw2();
      //std::cout << "INITIALIZED: " << KEY_RC1 << std::endl;
      KEY_MC1 = Form("mckaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_MC[KEY_MC1] = new TH2F(KEY_MC1.Data(),KEY_MC1.Data(),720,0.0,2.0*TMath::Pi(),200,0.0,0.9);
      h_mKaon_MC[KEY_MC1]->Sumw2();
      //std::cout << "INITIALIZED: " << KEY_MC1 << std::endl;
      KEY_RC1 = Form("rckaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_RC[KEY_RC1] = new TH2F(KEY_RC1.Data(),KEY_RC1.Data(),720,0.0,2.0*TMath::Pi(),200,0.0,0.9);
      h_mKaon_RC[KEY_RC1]->Sumw2();
      //std::cout << "INITIALIZED: " << KEY_RC1 << std::endl;
      KEY_MC1 = Form("mckaonplusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_MC[KEY_MC1] = new TH2F(KEY_MC1.Data(),KEY_MC1.Data(),data_constants::kaon_azimuthal_bins*6,-0.5,0.5,data_constants::kaon_rapidity_bins*3,-0.4,0.4);
      h_mKaon_MC[KEY_MC1]->Sumw2();
      //std::cout << "INITIALIZED: " << KEY_MC1 << std::endl;
      KEY_RC1 = Form("rckaonplusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_RC[KEY_RC1] = new TH2F(KEY_RC1.Data(),KEY_RC1.Data(),data_constants::kaon_azimuthal_bins*6,-0.5,0.5,data_constants::kaon_rapidity_bins*3,-0.4,0.4);
      h_mKaon_RC[KEY_RC1]->Sumw2();
      //std::cout << "INITIALIZED: " << KEY_RC1 << std::endl;
      KEY_MC1 = Form("mckaonminusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_MC[KEY_MC1] = new TH2F(KEY_MC1.Data(),KEY_MC1.Data(),data_constants::kaon_azimuthal_bins*6,-0.5,0.5,data_constants::kaon_rapidity_bins*3,-0.4,0.4);
      h_mKaon_MC[KEY_MC1]->Sumw2();
      //std::cout << "INITIALIZED: " << KEY_MC1 << std::endl;
      KEY_RC1 = Form("rckaonminusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_RC[KEY_RC1] = new TH2F(KEY_RC1.Data(),KEY_RC1.Data(),data_constants::kaon_azimuthal_bins*6,-0.5,0.5,data_constants::kaon_rapidity_bins*3,-0.4,0.4);
      h_mKaon_RC[KEY_RC1]->Sumw2();
      //std::cout << "INITIALIZED: " << KEY_RC1 << std::endl;
    }
  }
}

void StEffHistManger::InitHist()
{
  //for(int i_pt = mpt_first; i_pt < mpt_last; ++i_pt) // use rebinned pt
  //{
  //  for(int i_eta = 0; i_eta < vmsa::eta_total; ++i_eta)
  //  {
  //    std::string HistName = Form("h_mMcEffCos_Cent_10_Pt_%d_Eta_%d",i_pt,i_eta);
  //    h_mMcEffEtaCos[i_pt][i_eta] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
  //    h_mMcEffEtaCos[i_pt][i_eta]->Sumw2();
  //    HistName = Form("h_mRcEffCos_Cent_10_Pt_%d_Eta_%d",i_pt,i_eta);
  //    h_mRcEffEtaCos[i_pt][i_eta] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2.0,-1.0,1.0);
  //    h_mRcEffEtaCos[i_pt][i_eta]->Sumw2();
  //  }
  //}

  const int nybins = 16;
  double ybins[nybins+1] = {0.0};
  for(int iy = 0; iy < nybins+1; iy++)
  {
    ybins[iy] = double(iy-8.)/8.;
    cout << "iy = " << iy << ", y = " << ybins[iy] << endl;
  } 

  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    std::string HistName = Form("h_mMcTracks_%d",i_cent);
    h_mMcTracks[i_cent] = new TH3D(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,vmsa::BinEta,-1.0*vmsa::mEtaMax,vmsa::mEtaMax,vmsa::BinPhi,-1.0*TMath::Pi(),TMath::Pi());
    h_mMcTracks[i_cent]->Sumw2();
    HistName = Form("h_mRcTracks_%d",i_cent);
    h_mRcTracks[i_cent] = new TH3D(HistName.c_str(),HistName.c_str(),vmsa::BinPt,0.0,vmsa::ptEffMax,vmsa::BinEta,-1.0*vmsa::mEtaMax,vmsa::mEtaMax,vmsa::BinPhi,-1.0*TMath::Pi(),TMath::Pi());
    h_mRcTracks[i_cent]->Sumw2();

    HistName = Form("h_mMcEffPtY_Cent_%d",i_cent);
    h_mMcEffPtY[i_cent] = new TH2D(HistName.c_str(),HistName.c_str(),vmsa::pt_total,vmsa::pt_bin,nybins,ybins);
    h_mMcEffPtY[i_cent]->Sumw2();
    HistName = Form("h_mRcEffPtY_Cent_%d",i_cent);
    h_mRcEffPtY[i_cent] = new TH2D(HistName.c_str(),HistName.c_str(),vmsa::pt_total,vmsa::pt_bin,nybins,ybins);
    h_mRcEffPtY[i_cent]->Sumw2();
    
    for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
    {
      HistName = Form("h_mMcEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcEffCos[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      h_mMcEffCos[i_cent][i_pt]->Sumw2();
      HistName = Form("h_mRcEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcEffCos[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
      h_mRcEffCos[i_cent][i_pt]->Sumw2();
      HistName = Form("h_mMcEffPhiS_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcEffPhiS[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),10,-1.0,1.0);
      h_mMcEffPhiS[i_cent][i_pt]->Sumw2();
      HistName = Form("h_mRcEffPhiS_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcEffPhiS[i_cent][i_pt] = new TH1D(HistName.c_str(),HistName.c_str(),10,-1.0,1.0);
      h_mRcEffPhiS[i_cent][i_pt]->Sumw2();
      HistName = Form("h_mMcEffCosEP_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcEffCosEP[i_cent][i_pt] = new TH2D(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2.0,-1.0,1.0,7,-TMath::Pi()/2.0,TMath::Pi()/2.0);
      h_mMcEffCosEP[i_cent][i_pt]->Sumw2();
      HistName = Form("h_mRcEffCosEP_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcEffCosEP[i_cent][i_pt] = new TH2D(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2.0,-1.0,1.0,7,-TMath::Pi()/2.0,TMath::Pi()/2.0);
      h_mRcEffCosEP[i_cent][i_pt]->Sumw2();
      //////// phiprime vs cos(theta*)
      HistName = Form("h_mMcEffCosPhiPrime_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcEffCosPhiPrime[i_cent][i_pt] = new TH2D(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2,-1.0,1.0,vmsa::BinCos*2,0.0,2.0*TMath::Pi());
      h_mMcEffCosPhiPrime[i_cent][i_pt]->Sumw2();
      HistName = Form("h_mRcEffCosPhiPrime_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcEffCosPhiPrime[i_cent][i_pt] = new TH2D(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2,-1.0,1.0,vmsa::BinCos*2,0.0,2.0*TMath::Pi());
      h_mRcEffCosPhiPrime[i_cent][i_pt]->Sumw2();

      for(int i_y = 0; i_y < vmsa::y_total; i_y++)
      {
        //cos(theta)
        HistName = Form("h_mMcEffCos_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        h_mMcEffCosY[i_cent][i_pt][i_y] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        h_mMcEffCosY[i_cent][i_pt][i_y]->Sumw2();
        HistName = Form("h_mRcEffCos_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        h_mRcEffCosY[i_cent][i_pt][i_y] = new TH1D(HistName.c_str(),HistName.c_str(),vmsa::BinCos,0.0,1.0);
        h_mRcEffCosY[i_cent][i_pt][i_y]->Sumw2();
        HistName = Form("h_mMcEffPhiS_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        h_mMcEffPhiSY[i_cent][i_pt][i_y] = new TH1D(HistName.c_str(),HistName.c_str(),10,-1.0,1.0);
        h_mMcEffPhiSY[i_cent][i_pt][i_y]->Sumw2();
        HistName = Form("h_mRcEffPhiS_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        h_mRcEffPhiSY[i_cent][i_pt][i_y] = new TH1D(HistName.c_str(),HistName.c_str(),10,-1.0,1.0);
        h_mRcEffPhiSY[i_cent][i_pt][i_y]->Sumw2();
        HistName = Form("h_mMcEffCosEP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        h_mMcEffCosEPY[i_cent][i_pt][i_y] = new TH2D(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2.0,-1.0,1.0,7,-TMath::Pi()/2.0,TMath::Pi()/2.0);
        h_mMcEffCosEPY[i_cent][i_pt][i_y]->Sumw2();
        HistName = Form("h_mRcEffCosEP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
        h_mRcEffCosEPY[i_cent][i_pt][i_y] = new TH2D(HistName.c_str(),HistName.c_str(),vmsa::BinCos*2.0,-1.0,1.0,7,-TMath::Pi()/2.0,TMath::Pi()/2.0);
        h_mRcEffCosEPY[i_cent][i_pt][i_y]->Sumw2();
      }
    }
  }
  flag_eff = 0;
  flag_eff_PtEtaPhi = 0;
  flag_eff_Cos = 0;
}

void StEffHistManger::FillPhiHistMc(int cent, float pt, float phi, float y, float phistar)
{
  if(cent < vmsa::cent_low[0] || cent > vmsa::cent_up[0]) return;

  double Cos2PhiStarPhi = TMath::Cos(2.0*(phistar-phi));

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    if(pt >= mpt_low[i_pt] && pt < mpt_up[i_pt])
    {
      int phi_pt_bin = i_pt - 2;
      for(int i_pt_hist = data_constants::phi_pt_start[phi_pt_bin]; i_pt_hist <= data_constants::phi_pt_stop[phi_pt_bin]; i_pt_hist++)
      {
        if(pt >= data_constants::phi_pt_low[i_pt_hist] && pt < data_constants::phi_pt_high[i_pt_hist])
        { 
          for(Int_t i_phi = 0; i_phi < 10; i_phi++)
          {
            //std::cout << "Cos2PhiStarPhi = " << Cos2PhiStarPhi << ", (float(i_phi)-5.)/5. = " << (float(i_phi)-5.)/5. << ", (float(i_phi)-4.)/5. = " << (float(i_phi)-4.)/5 << std::endl;
            if(Cos2PhiStarPhi >= (float(i_phi)-5.)/5. && Cos2PhiStarPhi < (float(i_phi)-4.)/5.)
            {
              TString KEY = Form("mcphi_phipt_%d_cos2phistarphi_%d",i_pt_hist,i_phi);
              h_mPhi_MC[KEY]->Fill(phi,y);
              KEY = Form("mcphi_phistarphi_phipt_%d",i_pt_hist);
              h_mPhi_MC[KEY]->Fill(phi,phistar);
              KEY = Form("mcphi_phistarmphi_phipt_%d",i_pt_hist);
              h_mPhi1D_MC[KEY]->Fill(phistar-phi);
              break; 
            }
          }
          break;
        }
      }
      break;
    }
  }
}

void StEffHistManger::FillPhiHistRc(int cent, float pt, float phi, float y, float phistar)
{
  if(cent < vmsa::cent_low[0] || cent > vmsa::cent_up[0]) return;

  double Cos2PhiStarPhi = TMath::Cos(2.0*(phistar-phi));

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    if(pt >= mpt_low[i_pt] && pt < mpt_up[i_pt])
    {
      int phi_pt_bin = i_pt - 2;
      for(int i_pt_hist = data_constants::phi_pt_start[phi_pt_bin]; i_pt_hist <= data_constants::phi_pt_stop[phi_pt_bin]; i_pt_hist++)
      {
        if(pt >= data_constants::phi_pt_low[i_pt_hist] && pt < data_constants::phi_pt_high[i_pt_hist])
        { 
          for(Int_t i_phi = 0; i_phi < 10; i_phi++)
          {
            //std::cout << "Cos2PhiStarPhi = " << Cos2PhiStarPhi << ", (float(i_phi)-5.)/5. = " << (float(i_phi)-5.)/5. << ", (float(i_phi)-4.)/5. = " << (float(i_phi)-4.)/5 << std::endl;
            if(Cos2PhiStarPhi >= (float(i_phi)-5.)/5. && Cos2PhiStarPhi < (float(i_phi)-4)/5.)
            {
              TString KEY = Form("rcphi_phipt_%d_cos2phistarphi_%d",i_pt_hist,i_phi);
              h_mPhi_RC[KEY]->Fill(phi,y);
              KEY = Form("rcphi_phistarphi_phipt_%d",i_pt_hist);
              h_mPhi_RC[KEY]->Fill(phi,phistar);
              KEY = Form("rcphi_phistarmphi_phipt_%d",i_pt_hist);
              h_mPhi1D_RC[KEY]->Fill(phistar-phi);
              break; 
            }
          }
          break;
        }
      }
      break;
    }
  }
}

void StEffHistManger::FillKaonHistMc(int cent, float phi, float pt, float kpt, float ky, float keta, float kphi, float phistar)
{
  if(cent < vmsa::cent_low[0] || cent > vmsa::cent_up[0]) return;

  double Cos2PhiStarPhi = TMath::Cos(2.0*(phistar-phi));

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    //std::cout << "pt = " << pt << ", i_pt = " << i_pt << ", mpt_low[i_pt] = " << mpt_low[i_pt] << ", mpt_up[i_pt] = " << mpt_up[i_pt] << std::endl;
    if(pt >= mpt_low[i_pt] && pt < mpt_up[i_pt])
    {
      int phi_pt_bin = i_pt - 2;
      for(Int_t i_phi = 0; i_phi < 10; i_phi++)
      {
        //std::cout << "Cos2PhiStarPhi = " << Cos2PhiStarPhi << ", (float(i_phi)-5.)/5. = " << (float(i_phi)-5.)/5. << ", (float(i_phi)-4.)/5. = " << (float(i_phi)-4.)/5 << std::endl;
        if(Cos2PhiStarPhi >= (float(i_phi)-5.)/5. && Cos2PhiStarPhi < (float(i_phi)-4)/5.)
        {
          //std::cout << "Passed Cos2PhiStarPhi Cut" << std::endl;
          for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
          {
            if(kpt >= data_constants::kaon_pt_low[i_pt_kaon] && kpt < data_constants::kaon_pt_high[i_pt_kaon])   
            {
              //std::cout << "Passed kpt cut" << std::endl;
              TString KEY = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",phi_pt_bin,i_phi,i_pt_kaon);
              //std::cout << "Filling with phi = " << kphi << std::endl;
              h_mKaon_MC[KEY]->Fill(kphi,ky);
              KEY = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",phi_pt_bin,i_phi,i_pt_kaon);
              //std::cout << "Filling with phi = " << kphi << std::endl;
              h_mKaon_MC[KEY]->Fill(kphi,keta);
              break;
            }
          }
          break;
        }
      }
      break;
    }
  }
}

void StEffHistManger::FillKaonDeltaHistMc(int cent, float phi, float eta, float pt, float kppt, float kpeta, float kpphi, float kmpt, float kmeta, float kmphi, float phistar)
{
  if(cent < vmsa::cent_low[0] || cent > vmsa::cent_up[0]) return;

  double Cos2PhiStarPhi = TMath::Cos(2.0*(phistar-phi));
  float phiphidiffp = phi-kpphi;
  float phietadiffp = eta-kpeta;
  float phiphidiffm = phi-kmphi;
  float phietadiffm = eta-kmeta;
  float phidiff = kpphi-kmphi;
  float etadiff = kpeta-kmeta;
  while(phiphidiffp < -TMath::Pi()) phiphidiffp += 2*TMath::Pi();
  while(phiphidiffp > TMath::Pi() ) phiphidiffp -= 2*TMath::Pi();
  while(phiphidiffm < -TMath::Pi()) phiphidiffm += 2*TMath::Pi();
  while(phiphidiffm > TMath::Pi() ) phiphidiffm -= 2*TMath::Pi();
  while(phidiff < -TMath::Pi()) phidiff += 2*TMath::Pi();
  while(phidiff > TMath::Pi() ) phidiff -= 2*TMath::Pi();

  float deltar = TMath::Sqrt(etadiff*etadiff+phidiff*phidiff);
  float phi_r = TMath::ATan2(etadiff,phidiff);
  if(phi_r < 0.0) phi_r += 2.0*TMath::Pi();
  //if(phidiff < 0.0) phi_r += TMath::Pi(); 
  //if(phi_r < 0.0) phi_r += 2.0*TMath::Pi();
  if(deltar > 0.9) 
  {
    cout << "Outside of r Range" << endl;
    cout << "r = " << deltar << endl;
    cout << "delta(eta) = " << etadiff << ", delta(phi) = " << phidiff << endl << endl;
  }
  if(phi_r > 2.0*TMath::Pi() || phi_r < 0.0)
  {
    cout << "Outside of phi Range" << endl;
    cout << "phi = " << phi_r << endl;
    cout << "delta(eta) = " << etadiff << ", delta(phi) = " << phidiff << endl << endl;
  }

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    //std::cout << "pt = " << pt << ", i_pt = " << i_pt << ", mpt_low[i_pt] = " << mpt_low[i_pt] << ", mpt_up[i_pt] = " << mpt_up[i_pt] << std::endl;
    if(pt >= mpt_low[i_pt] && pt < mpt_up[i_pt])
    {
      int phi_pt_bin = i_pt - 2;
      for(Int_t i_phi = 0; i_phi < 10; i_phi++)
      {
        //std::cout << "Cos2PhiStarPhi = " << Cos2PhiStarPhi << ", (float(i_phi)-5.)/5. = " << (float(i_phi)-5.)/5. << ", (float(i_phi)-4.)/5. = " << (float(i_phi)-4.)/5 << std::endl;
        if(Cos2PhiStarPhi >= (float(i_phi)-5.)/5. && Cos2PhiStarPhi < (float(i_phi)-4)/5.)
        {
          TString KEY = Form("mckaon_delta_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_MC[KEY]->Fill(phidiff,etadiff);
          KEY = Form("mckaon_deltarphi_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_MC[KEY]->Fill(phi_r,deltar);
          KEY = Form("mckaonplusphi_delta_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_MC[KEY]->Fill(phiphidiffp,phietadiffp);
          KEY = Form("mckaonminusphi_delta_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_MC[KEY]->Fill(phiphidiffm,phietadiffm);
          break;
        }
      }
      break;
    }
  }
}

void StEffHistManger::FillKaonHistRc(int cent, float phi, float pt, float kpt, float ky, float keta, float kphi, float phistar)
{
  if(cent < vmsa::cent_low[0] || cent > vmsa::cent_up[0]) return;

  double Cos2PhiStarPhi = TMath::Cos(2.0*(phistar-phi));

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    if(pt >= mpt_low[i_pt] && pt < mpt_up[i_pt])
    {
      int phi_pt_bin = i_pt - 2;
      for(Int_t i_phi = 0; i_phi < 10; i_phi++)
      {
        if(Cos2PhiStarPhi >= (float(i_phi)-5.)/5. && Cos2PhiStarPhi < (float(i_phi)-4)/5.)
        {
          for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
          {
            if(kpt >= data_constants::kaon_pt_low[i_pt_kaon] && kpt < data_constants::kaon_pt_high[i_pt_kaon])   
            {
              TString KEY = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",phi_pt_bin,i_phi,i_pt_kaon);
              h_mKaon_RC[KEY]->Fill(kphi,ky);
              KEY = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",phi_pt_bin,i_phi,i_pt_kaon);
              h_mKaon_RC[KEY]->Fill(kphi,keta);
              break;
            }
          }
          break;
        }
      }
      break;
    }
  }
}
void StEffHistManger::FillKaonDeltaHistRc(int cent, float phi, float eta, float pt, float kppt, float kpeta, float kpphi, float kmpt, float kmeta, float kmphi, float phistar)
{
  if(cent < vmsa::cent_low[0] || cent > vmsa::cent_up[0]) return;

  double Cos2PhiStarPhi = TMath::Cos(2.0*(phistar-phi));
  float phiphidiffp = phi-kpphi;
  float phietadiffp = eta-kpeta;
  float phiphidiffm = phi-kmphi;
  float phietadiffm = eta-kmeta;
  float phidiff = kpphi-kmphi;
  float etadiff = kpeta-kmeta;
  while(phiphidiffp < -TMath::Pi()) phiphidiffp += 2*TMath::Pi();
  while(phiphidiffp > TMath::Pi() ) phiphidiffp -= 2*TMath::Pi();
  while(phiphidiffm < -TMath::Pi()) phiphidiffm += 2*TMath::Pi();
  while(phiphidiffm > TMath::Pi() ) phiphidiffm -= 2*TMath::Pi();
  while(phidiff < -TMath::Pi()) phidiff += 2*TMath::Pi();
  while(phidiff > TMath::Pi() ) phidiff -= 2*TMath::Pi();
 
  float deltar = TMath::Sqrt(etadiff*etadiff+phidiff*phidiff);
  float phi_r = TMath::ATan2(etadiff,phidiff);
  if(phi_r < 0.0) phi_r += 2.0*TMath::Pi();
  //if(phidiff < 0.0) phi_r += TMath::Pi(); 
  //if(phi_r < 0.0) phi_r += 2.0*TMath::Pi();
  if(deltar > 0.9) 
  {
    cout << "Outside of r Range" << endl;
    cout << "r = " << deltar << endl;
    cout << "delta(eta) = " << etadiff << ", delta(phi) = " << phidiff << endl << endl;
  }
  if(phi_r > 2.0*TMath::Pi() || phi_r < 0.0)
  {
    cout << "Outside of phi Range" << endl;
    cout << "phi = " << phi_r << endl;
    cout << "delta(eta) = " << etadiff << ", delta(phi) = " << phidiff << endl << endl;
  }


  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    //std::cout << "pt = " << pt << ", i_pt = " << i_pt << ", mpt_low[i_pt] = " << mpt_low[i_pt] << ", mpt_up[i_pt] = " << mpt_up[i_pt] << std::endl;
    if(pt >= mpt_low[i_pt] && pt < mpt_up[i_pt])
    {
      int phi_pt_bin = i_pt - 2;
      for(Int_t i_phi = 0; i_phi < 10; i_phi++)
      {
        //std::cout << "Cos2PhiStarPhi = " << Cos2PhiStarPhi << ", (float(i_phi)-5.)/5. = " << (float(i_phi)-5.)/5. << ", (float(i_phi)-4.)/5. = " << (float(i_phi)-4.)/5 << std::endl;
        if(Cos2PhiStarPhi >= (float(i_phi)-5.)/5. && Cos2PhiStarPhi < (float(i_phi)-4)/5.)
        {
          TString KEY = Form("rckaon_delta_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_RC[KEY]->Fill(phidiff,etadiff);
          KEY = Form("rckaon_deltarphi_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_RC[KEY]->Fill(phi_r,deltar);
          KEY = Form("rckaonplusphi_delta_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_RC[KEY]->Fill(phiphidiffp,phietadiffp);
          KEY = Form("rckaonminusphi_delta_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_RC[KEY]->Fill(phiphidiffm,phietadiffm);
          break;
        }
      }
      break;
    }
  }
}

void StEffHistManger::WritePhiHist()
{
  for(Int_t i_pt_phi = mpt_first-2; i_pt_phi <= mpt_last-2; ++i_pt_phi) // use rebinned pt
  {
    for(Int_t i_pt = data_constants::phi_pt_start[i_pt_phi]; i_pt <= data_constants::phi_pt_stop[i_pt_phi]; i_pt++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++)
      {
        TString KEY_MC = Form("mcphi_phipt_%d_cos2phistarphi_%d",i_pt,i_phi);
        h_mPhi_MC[KEY_MC]->Write();
        std::cout << "WROTE: " << KEY_MC << std::endl;
        TString KEY_RC = Form("rcphi_phipt_%d_cos2phistarphi_%d",i_pt,i_phi);
        h_mPhi_RC[KEY_RC]->Write();
        std::cout << "WROTE: " << KEY_RC << std::endl;
      }

      TString KEY_MC = Form("mcphi_phistarphi_phipt_%d",i_pt);
      h_mPhi_MC[KEY_MC]->Write();
      std::cout << "WROTE: " << KEY_MC << std::endl;
      TString KEY_RC = Form("rcphi_phistarphi_phipt_%d",i_pt);
      h_mPhi_RC[KEY_RC]->Write();
      std::cout << "WROTE: " << KEY_RC << std::endl;

      KEY_MC = Form("mcphi_phistarmphi_phipt_%d",i_pt);
      h_mPhi1D_MC[KEY_MC]->Write();
      std::cout << "WROTE: " << KEY_MC << std::endl;
      KEY_RC = Form("rcphi_phistarmphi_phipt_%d",i_pt);
      h_mPhi1D_RC[KEY_RC]->Write();
      std::cout << "WROTE: " << KEY_RC << std::endl;
    }
  }
}
void StEffHistManger::WriteKaonHist()
{
  for(Int_t i_pt_phi = mpt_first-2; i_pt_phi <= mpt_last-2; ++i_pt_phi) // use rebinned pt
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++)
    {
      for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
      {
        TString KEY_MC = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        h_mKaon_MC[KEY_MC]->Write();
        std::cout << "WROTE: " << KEY_MC << std::endl;
        TString KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        h_mKaon_RC[KEY_RC]->Write();
        std::cout << "WROTE: " << KEY_RC << std::endl;

        KEY_MC = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);
        h_mKaon_MC[KEY_MC]->Write();
        std::cout << "WROTE: " << KEY_MC << std::endl;
        KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);
        h_mKaon_RC[KEY_RC]->Write();
        std::cout << "WROTE: " << KEY_RC << std::endl;
      }
      TString KEY_MC1 = Form("mckaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_MC[KEY_MC1]->Write();
      TString KEY_RC1 = Form("rckaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_RC[KEY_RC1]->Write();
      KEY_MC1 = Form("mckaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_MC[KEY_MC1]->Write();
      KEY_RC1 = Form("rckaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_RC[KEY_RC1]->Write();
      KEY_MC1 = Form("mckaonplusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_MC[KEY_MC1]->Write();
      KEY_RC1 = Form("rckaonplusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_RC[KEY_RC1]->Write();
      KEY_MC1 = Form("mckaonminusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_MC[KEY_MC1]->Write();
      KEY_RC1 = Form("rckaonminusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mKaon_RC[KEY_RC1]->Write();
    }
  }
}

void StEffHistManger::FillHistMc(int cent, float pt, float y, float phi, float cos, float Psi, float phistar, float phiprime)
{
  double ratio_weight = 1.0;//h_mSpectraRatio->GetBinContent(h_mSpectraRatio->FindBin(pt,y));
  double phidiff = TMath::Cos(2.0*(phistar-phi));
  double phiPsi = AngleShift(phi-Psi);
  h_mMcTracks[cent]->Fill(pt,y,phi);
  if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
  {//20-60%
    h_mMcTracks[9]->Fill(pt,y,phi,vmsa::weight[cent]*ratio_weight);
  }

  h_mMcEffPtY[cent]->Fill(pt,y);
  if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
  {//20-60%
    h_mMcEffPtY[9]->Fill(pt,y,ratio_weight);
  }

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    if(pt > mpt_low[i_pt] && pt <= mpt_up[i_pt])
    {
      if(mMode == 1) 
      {
        h_mMcEffCos[cent][i_pt]->Fill(cos);
        h_mMcEffPhiS[cent][i_pt]->Fill(phidiff);
        h_mMcEffCosEP[cent][i_pt]->Fill(cos,phiPsi);
      }
      if(mMode == 2)
      {
        for(int i_y = 0; i_y < vmsa::y_total; ++i_y)
        {
          if(y > vmsa::y_bin[i_y] && y < vmsa::y_bin[i_y+1])
          {
            //std::cout << "y = " << y << "     y_bin[i_y] = " << vmsa::y_bin[i_y] << "     y_bin[i_y+1] = " << vmsa::y_bin[i_y+1] << std::endl;
            h_mMcEffCosY[cent][i_pt][i_y]->Fill(cos);
            h_mMcEffPhiSY[cent][i_pt][i_y]->Fill(phidiff);
            h_mMcEffCosEPY[cent][i_pt][i_y]->Fill(cos,phiPsi);
          }
        }
      }
    }
  }
  
  if(mMode == 0)
  {
    for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
    {
      if(pt > mpt_low[i_pt] && pt <= mpt_up[i_pt])
      {
        h_mMcEffCos[cent][i_pt]->Fill(cos);
        h_mMcEffPhiS[cent][i_pt]->Fill(phidiff);
        h_mMcEffCosEP[cent][i_pt]->Fill(cos,phiPsi);
        h_mMcEffCosPhiPrime[cent][i_pt]->Fill(cos,phiprime);
        if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
        {//20-60%
          h_mMcEffCos[9][i_pt]->Fill(cos,vmsa::weight[cent]*ratio_weight);
          h_mMcEffPhiS[9][i_pt]->Fill(phidiff,vmsa::weight[cent]*ratio_weight);
          h_mMcEffCosEP[9][i_pt]->Fill(cos,phiPsi,vmsa::weight[cent]*ratio_weight);
          h_mMcEffCosPhiPrime[9][i_pt]->Fill(cos,phiprime,ratio_weight);
//  	h_mMcCosEP[9][i_pt]->Fill(cos,phi_shift,vmsa::weight[cent]);
          //cout << "Filling MC cent9" << endl;
        }
      }
    }
  }
  if(mMode == 3)
  {
    for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
    {
      if(pt > mpt_low[i_pt] && pt <= mpt_up[i_pt])
      {
        for(int i_y = 0; i_y < vmsa::y_rebin; ++i_y)
        {
          if(y > vmsa::y_pt_rebin[i_y] && y < vmsa::y_pt_rebin[i_y+1])
          {
            h_mMcEffCosY[cent][i_pt][i_y]->Fill(cos);
            h_mMcEffPhiSY[cent][i_pt][i_y]->Fill(phidiff);
            h_mMcEffCosEPY[cent][i_pt][i_y]->Fill(cos,phiPsi);
      //  h_mMcCosEP[cent][i_pt]->Fill(cos,phi_shift);
            if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
            {//20-60%
              h_mMcEffCosY[9][i_pt][i_y]->Fill(cos,vmsa::weight[cent]);
              h_mMcEffPhiSY[9][i_pt][i_y]->Fill(phidiff,vmsa::weight[cent]);
              h_mMcEffCosEPY[9][i_pt][i_y]->Fill(cos,phiPsi,vmsa::weight[cent]);
//  	    h_mMcCosEP[9][i_pt]->Fill(cos,phi_shift,vmsa::weight[cent]);
              //std::cout << "Filled MC" << std::endl;
            }
          }
        }
      }
    }
  }
}

void StEffHistManger::FillHistRc(int cent, float pt, float y, float phi, float cos, float Psi, float phistar, float phiprime)
{
  double ratio_weight = 1.0;//h_mSpectraRatio->GetBinContent(h_mSpectraRatio->FindBin(pt,y));
  double phidiff = TMath::Cos(2.0*(phistar-phi));
  double phiPsi = AngleShift(phi-Psi);
  h_mRcTracks[cent]->Fill(pt,y,phi);
  if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
  {//20-60%
    h_mRcTracks[9]->Fill(pt,y,phi,vmsa::weight[cent]*ratio_weight);
  }

  h_mRcEffPtY[cent]->Fill(pt,y);
  if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
  {//20-60%
    h_mRcEffPtY[9]->Fill(pt,y,ratio_weight);
  }

  //float phi_shift = AngleShift(phi-Psi2);
  

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    //std::cout << "pt = " << pt << "    mpt_low[i_pt] = " << mpt_low[i_pt] << "       mpt_up[i_pt] = " << mpt_up[i_pt] << std::endl;
    if(pt > mpt_low[i_pt] && pt <= mpt_up[i_pt])
    {
      if(mMode == 1) 
      {
        h_mRcEffCos[cent][i_pt]->Fill(cos);
        h_mRcEffPhiS[cent][i_pt]->Fill(phidiff);
        h_mRcEffCosEP[cent][i_pt]->Fill(cos,phiPsi);
      }
      if(mMode == 2)
      {
        for(int i_y = 0; i_y < vmsa::y_total; ++i_y)
        {
          if(y > vmsa::y_bin[i_y] && y < vmsa::y_bin[i_y+1])
          {
            //std::cout << "y = " << y << "     y_bin[i_y] = " << vmsa::y_bin[i_y] << "     y_bin[i_y+1] = " << vmsa::y_bin[i_y+1] << std::endl;
            h_mRcEffCosY[cent][i_pt][i_y]->Fill(cos);
            h_mRcEffPhiSY[cent][i_pt][i_y]->Fill(phidiff);
            h_mRcEffCosEPY[cent][i_pt][i_y]->Fill(cos,phiPsi);
          }
        }
      }
    }
  }
  
  if(mMode == 0)
  {
    for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
    {
      if(pt > mpt_low[i_pt]&& pt <= mpt_up[i_pt])
      {
        h_mRcEffCos[cent][i_pt]->Fill(cos);
        h_mRcEffPhiS[cent][i_pt]->Fill(phidiff);
        h_mRcEffCosEP[cent][i_pt]->Fill(cos,phiPsi);
        h_mRcEffCosPhiPrime[cent][i_pt]->Fill(cos,phiprime);
      //  h_mRcCosEP[cent][i_pt]->Fill(cos,phi_shift);
        if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
        {//20-60%
          h_mRcEffCos[9][i_pt]->Fill(cos,vmsa::weight[cent]*ratio_weight);
          h_mRcEffPhiS[9][i_pt]->Fill(phidiff,vmsa::weight[cent]*ratio_weight);
          h_mRcEffCosEP[9][i_pt]->Fill(cos,phiPsi,vmsa::weight[cent]*ratio_weight);
          h_mRcEffCosPhiPrime[9][i_pt]->Fill(cos,phiprime,ratio_weight);
//  	h_mRcCosEP[9][i_pt]->Fill(cos,phi_shift,vmsa::weight[cent]);
          //cout << "Filling RC cent9" << endl;
        }
      }
    }
  }
  if(mMode == 3)
  {
    for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
    {
      if(pt > mpt_low[i_pt] && pt <= mpt_up[i_pt])
      {
        for(int i_y = 0; i_y < vmsa::y_rebin; ++i_y)
        {
          if(y > vmsa::y_pt_rebin[i_y] && y < vmsa::y_pt_rebin[i_y+1])
          {
            h_mRcEffCosY[cent][i_pt][i_y]->Fill(cos);
            h_mRcEffPhiSY[cent][i_pt][i_y]->Fill(phidiff);
            h_mRcEffCosEPY[cent][i_pt][i_y]->Fill(cos,phiPsi);
      //  h_mMcCosEP[cent][i_pt]->Fill(cos,phi_shift);
            if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
            {//20-60%
              h_mRcEffCosY[9][i_pt][i_y]->Fill(cos,vmsa::weight[cent]);
              h_mRcEffPhiSY[9][i_pt][i_y]->Fill(phidiff,vmsa::weight[cent]);
              h_mRcEffCosEPY[9][i_pt][i_y]->Fill(phiPsi,vmsa::weight[cent]);
//  	    h_mMcCosEP[9][i_pt]->Fill(cos,phi_shift,vmsa::weight[cent]);
              //std::cout << "Filled RC" << std::endl;
            }
          }
        }
      }
    }
  }
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

TH1D* StEffHistManger::CalEffError(TH1D *h_Mc, TH1D *h_Rc, std::string HistName)
{
  TH1D* h_ratio = (TH1D*)h_Rc->Clone();
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

TH2D* StEffHistManger::CalEffError(TH2D *h_Mc, TH2D *h_Rc, std::string HistName)
{
  TH2D* h_ratio = (TH2D*)h_Rc->Clone();
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

void StEffHistManger::CalEfficiency()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    std::string HistName;
    HistName = Form("h_mMcEffPt_Cent_%d",i_cent);
    h_mMcEffPt[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("x")->Clone(HistName.c_str());
    HistName = Form("h_mRcEffPt_Cent_%d",i_cent);
    h_mRcEffPt[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("x")->Clone(HistName.c_str());
    HistName = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEffPt[i_cent] = CalEffError(h_mMcEffPt[i_cent],h_mRcEffPt[i_cent],HistName.c_str());

    HistName = Form("h_mMcEffEta_Cent_%d",i_cent);
    h_mMcEffEta[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("y")->Clone(HistName.c_str());
    HistName = Form("h_mRcEffEta_Cent_%d",i_cent);
    h_mRcEffEta[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("y")->Clone(HistName.c_str());
    HistName = Form("h_mEffEta_Cent_%d",i_cent);
    h_mEffEta[i_cent] = CalEffError(h_mMcEffEta[i_cent],h_mRcEffEta[i_cent],HistName.c_str());

    HistName = Form("h_mMcEffPhi_Cent_%d",i_cent);
    h_mMcEffPhi[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("z")->Clone(HistName.c_str());
    HistName = Form("h_mRcEffPhi_Cent_%d",i_cent);
    h_mRcEffPhi[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("z")->Clone(HistName.c_str());
    HistName = Form("h_mEffPhi_Cent_%d",i_cent);
    h_mEffPhi[i_cent] = CalEffError(h_mMcEffPhi[i_cent],h_mRcEffPhi[i_cent],HistName.c_str());
  }
  flag_eff = 1;
}

void StEffHistManger::CalEffPtEtaPhi()
{
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
	h_mEffPEP[HistNameEff] = CalEffError(h_mMcEffPEP[HistNameMc],h_mRcEffPEP[HistNameRc],HistNameEff.c_str());
      }
    }
  }
  flag_eff_PtEtaPhi = 1;
}

void StEffHistManger::CalEffCosThetaStar()
{
  if(mMode == 0)
  {
    for(int i_cent = 0; i_cent < 10; ++i_cent)
    {
    
      std::string HistName = Form("h_mEffPtY_Cent_%d",i_cent);
      h_mEffPtY[i_cent] = CalEffError(h_mMcEffPtY[i_cent],h_mRcEffPtY[i_cent],HistName.c_str());
      for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
      {
        HistName = Form("h_mEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCos[i_cent][i_pt] = CalEffError(h_mMcEffCos[i_cent][i_pt],h_mRcEffCos[i_cent][i_pt],HistName.c_str());
        HistName = Form("h_mEffPhiS_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffPhiS[i_cent][i_pt] = CalEffError(h_mMcEffPhiS[i_cent][i_pt],h_mRcEffPhiS[i_cent][i_pt],HistName.c_str());
        HistName = Form("h_mEffCosEP_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCosEP[i_cent][i_pt] = CalEffError(h_mMcEffCosEP[i_cent][i_pt],h_mRcEffCosEP[i_cent][i_pt],HistName.c_str());
        HistName = Form("h_mEffCosPhiPrime_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCosPhiPrime[i_cent][i_pt] = CalEffError(h_mMcEffCosPhiPrime[i_cent][i_pt],h_mRcEffCosPhiPrime[i_cent][i_pt],HistName.c_str());
      }
    }
  }
  if(mMode == 3)
  {
    for(int i_cent = 0; i_cent < 10; ++i_cent)
    {
      for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
      {
        for(int i_y = 0; i_y < vmsa::y_rebin; ++i_y)
        {      
          std::string HistName = Form("h_mEffCos_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
          h_mEffCosY[i_cent][i_pt][i_y] = CalEffError(h_mMcEffCosY[i_cent][i_pt][i_y],h_mRcEffCosY[i_cent][i_pt][i_y],HistName.c_str());
          HistName = Form("h_mEffPhiS_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
          h_mEffPhiSY[i_cent][i_pt][i_y] = CalEffError(h_mMcEffPhiSY[i_cent][i_pt][i_y],h_mRcEffPhiSY[i_cent][i_pt][i_y],HistName.c_str());
          HistName = Form("h_mEffCosEP_Cent_%d_Pt_%d_Y_%d",i_cent,i_pt,i_y);
          h_mEffCosEPY[i_cent][i_pt][i_y] = CalEffError(h_mMcEffCosEPY[i_cent][i_pt][i_y],h_mRcEffCosEPY[i_cent][i_pt][i_y],HistName.c_str());
        }
      }
    }
  }
 
  int cent_bin = 0;
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
    {
      if(mMode == 1)
      {
        std::string HistName = Form("h_mEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCos[i_cent][i_pt] = CalEffError(h_mMcEffCos[i_cent][i_pt],h_mRcEffCos[i_cent][i_pt],HistName.c_str());
        HistName = Form("h_mEffPhiS_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffPhiS[i_cent][i_pt] = CalEffError(h_mMcEffPhiS[i_cent][i_pt],h_mRcEffPhiS[i_cent][i_pt],HistName.c_str());
        HistName = Form("h_mEffCosEP_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCosEP[i_cent][i_pt] = CalEffError(h_mMcEffCosEP[i_cent][i_pt],h_mRcEffCosEP[i_cent][i_pt],HistName.c_str());
      }
      if(mMode == 2)
      {
        if(i_cent == 9) continue;
        for(int i_y = 0; i_y < vmsa::y_total; ++i_y)
        {      
          std::string HistName = Form("h_mEffCos_Cent_%d_Pt_%d_Y_%d",cent_bin,i_pt,i_y);

          std::cout << "Cent = " << i_cent << std::endl;
          std::cout << "Cent rebin = " << cent_bin << std::endl;
          if(i_cent > vmsa::cent_rebin[cent_bin] && i_cent < vmsa::cent_rebin[cent_bin+1])
          {
            std::cout << "add cent = " << i_cent << std::endl;
            h_mMcEffCosY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mMcEffCosY[i_cent][i_pt][i_y],1.0);
            h_mRcEffCosY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mRcEffCosY[i_cent][i_pt][i_y],1.0);
          }
          if(i_cent == vmsa::cent_rebin[cent_bin + 1] - 1) 
          {
            std::cout << "cal eff for cent = " << cent_bin << std::endl;
            h_mEffCosY[cent_bin][i_pt][i_y] = CalEffError(h_mMcEffCosY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],h_mRcEffCosY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],HistName.c_str());
          }

          HistName = Form("h_mEffPhiS_Cent_%d_Pt_%d_Y_%d",cent_bin,i_pt,i_y);

          if(i_cent > vmsa::cent_rebin[cent_bin] && i_cent < vmsa::cent_rebin[cent_bin+1])
          {
            std::cout << "add cent = " << i_cent << std::endl;
            h_mMcEffPhiSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mMcEffPhiSY[i_cent][i_pt][i_y],1.0);
            h_mRcEffPhiSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mRcEffPhiSY[i_cent][i_pt][i_y],1.0);
          }
          if(i_cent == vmsa::cent_rebin[cent_bin + 1] - 1) 
          {
            std::cout << "cal eff for cent = " << cent_bin << std::endl;
            h_mEffPhiSY[cent_bin][i_pt][i_y] = CalEffError(h_mMcEffPhiSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],h_mRcEffPhiSY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],HistName.c_str());
          }

          HistName = Form("h_mEffCosEP_Cent_%d_Pt_%d_Y_%d",cent_bin,i_pt,i_y);

          if(i_cent > vmsa::cent_rebin[cent_bin] && i_cent < vmsa::cent_rebin[cent_bin+1])
          {
            std::cout << "add cent = " << i_cent << std::endl;
            h_mMcEffCosEPY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mMcEffCosEPY[i_cent][i_pt][i_y],1.0);
            h_mRcEffCosEPY[vmsa::cent_rebin[cent_bin]][i_pt][i_y]->Add(h_mRcEffCosEPY[i_cent][i_pt][i_y],1.0);
          }
          if(i_cent == vmsa::cent_rebin[cent_bin + 1] - 1) 
          {
            std::cout << "cal eff for cent = " << cent_bin << std::endl;
            h_mEffCosEPY[cent_bin][i_pt][i_y] = CalEffError(h_mMcEffCosEPY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],h_mRcEffCosEPY[vmsa::cent_rebin[cent_bin]][i_pt][i_y],HistName.c_str());
          }
        }
      }
    }
    if(i_cent == vmsa::cent_rebin[cent_bin + 1] - 1) cent_bin++;
  }


  flag_eff_Cos = 1;
}

void StEffHistManger::WriteHist()
{
  if(mMode == 1)
  {
    if(flag_eff_Cos > 0.5)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
        for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
        {
          h_mMcEffCos[i_cent][i_pt]->Write();
          h_mRcEffCos[i_cent][i_pt]->Write();
          h_mEffCos[i_cent][i_pt]->Write();
          h_mMcEffPhiS[i_cent][i_pt]->Write();
          h_mRcEffPhiS[i_cent][i_pt]->Write();
          h_mEffPhiS[i_cent][i_pt]->Write();
          h_mMcEffCosEP[i_cent][i_pt]->Write();
          h_mRcEffCosEP[i_cent][i_pt]->Write();
          h_mEffCosEP[i_cent][i_pt]->Write();
        }
      }
    }
  }
 
  if(mMode == 2)
  {
    if(flag_eff_Cos > 0.5)
    {
      for(int i_cent = 0; i_cent < vmsa::cent_rebin_total; ++i_cent)
      {
        for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
        {
          for(int i_y = 0; i_y < vmsa::y_total; ++i_y)
          {
            h_mMcEffCosY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->Write();
            h_mRcEffCosY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->Write();
              h_mEffCosY[i_cent][i_pt][i_y]->Write();
            h_mMcEffPhiSY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->Write();
            h_mRcEffPhiSY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->Write();
              h_mEffPhiSY[i_cent][i_pt][i_y]->Write();
            h_mMcEffCosEPY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->Write();
            h_mRcEffCosEPY[vmsa::cent_rebin[i_cent]][i_pt][i_y]->Write();
              h_mEffCosEPY[i_cent][i_pt][i_y]->Write();
          }
        }
      }
    }
  }

  if(mMode == 3)
  {
    if(flag_eff_Cos > 0.5)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
        for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
        {
          for(int i_y = 0; i_y < vmsa::y_rebin; ++i_y)
          {
            h_mMcEffCosY[i_cent][i_pt][i_y]->Write();
            h_mRcEffCosY[i_cent][i_pt][i_y]->Write();
              h_mEffCosY[i_cent][i_pt][i_y]->Write();
            h_mMcEffPhiSY[i_cent][i_pt][i_y]->Write();
            h_mRcEffPhiSY[i_cent][i_pt][i_y]->Write();
              h_mEffPhiSY[i_cent][i_pt][i_y]->Write();
            h_mMcEffCosEPY[i_cent][i_pt][i_y]->Write();
            h_mRcEffCosEPY[i_cent][i_pt][i_y]->Write();
              h_mEffCosEPY[i_cent][i_pt][i_y]->Write();
          }
        }
      }
    }
  }
  if(mMode == 0)
  {
    for(int i_cent = 0; i_cent < 10; ++i_cent)
    {
      h_mMcTracks[i_cent]->Write();
      h_mRcTracks[i_cent]->Write();

      if(flag_eff > 0.5)
      {
        h_mEffPt[i_cent]->Write();
        h_mEffEta[i_cent]->Write();
        h_mEffPhi[i_cent]->Write();
      }

      if(flag_eff_PtEtaPhi > 0.5)
      {
        for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
        {
          for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
          {
            std::string HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
            h_mEffPEP[HistNameEff]->Write();
          }
        }
      }

      if(flag_eff_Cos > 0.5)
      {
        h_mMcEffPtY[i_cent]->Write();
        h_mRcEffPtY[i_cent]->Write();
        h_mEffPtY[i_cent]->Write();
        for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt)
        {
          h_mMcEffCos[i_cent][i_pt]->Write();
          h_mRcEffCos[i_cent][i_pt]->Write();
            h_mEffCos[i_cent][i_pt]->Write();
          h_mMcEffPhiS[i_cent][i_pt]->Write();
          h_mRcEffPhiS[i_cent][i_pt]->Write();
            h_mEffPhiS[i_cent][i_pt]->Write();
          h_mMcEffCosEP[i_cent][i_pt]->Write();
          h_mRcEffCosEP[i_cent][i_pt]->Write();
            h_mEffCosEP[i_cent][i_pt]->Write();
          h_mMcEffCosPhiPrime[i_cent][i_pt]->Write();
          h_mRcEffCosPhiPrime[i_cent][i_pt]->Write();
            h_mEffCosPhiPrime[i_cent][i_pt]->Write();
          //cout << "Write histograms" << endl;

        }
      }
    }
  }
}
