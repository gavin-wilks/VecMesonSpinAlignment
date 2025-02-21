#include "StEffHistMangerHelicityGlobal.h"
#include "StRoot/Utility/phi_data_constants_19GeV.h"
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

ClassImp(StEffHistMangerHelicityGlobal)
//
StEffHistMangerHelicityGlobal::StEffHistMangerHelicityGlobal(int energy, int pid, int mode, int startpt, int stoppt, int ptfixed, int yfixed)
{
  mEnergy = energy;
  mMode = mode;
  mStartPt = startpt;
  mStopPt = stoppt;
  //mOrder = order;
  
  string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/ptyspectra_datarcratio_%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Spectra = TFile::Open(inputfile.c_str());
  h_mSpectraRatio = (TH2D*) ((TH2D*) File_Spectra->Get("pty_datarcratio"))->Clone();
  cout << "Loaded the phi ptyspectra" << endl;

  //string inputfilekaon = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/kaon_ptyratio_%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  //TFile *File_Kaon_Spectra = TFile::Open(inputfilekaon.c_str());
  //for(int i = 0; i < 4; i++)
  //{
  //  h_mKaonSpectraRatio[i] = (TH2D*) ((TH2D*) File_Kaon_Spectra->Get(Form("ratio_pty_phipt_%d",i)))->Clone();
  //}
  //cout << "Loaded the kaon ptyspectra" << endl;

  string inputfilecos = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/CosCosHDataRcRatio_%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Cos = TFile::Open(inputfilecos.c_str());
  for(int i = 0; i < 4; i++)
  {
    h_mCosCosHRatio[i] = (TH2D*) ((TH2D*) File_Cos->Get(Form("CosCosHDataRcRatio_pt_%d",i+2)))->Clone();
  }
  cout << "Loaded the 2D cos cos " << endl;

  //inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/v2weighting_%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  //TFile *File_v2 = TFile::Open(inputfile.c_str());
  //for(int i = 0; i < 11; i++)
  //{
  //  h_mV2[i] = (TH2D*) ((TH2D*) File_v2->Get(Form("Weight_v2_%d_pt_%d_y_%d",i,ptfixed,yfixed)))->Clone();
  //}
 
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

StEffHistMangerHelicityGlobal::~StEffHistMangerHelicityGlobal()
{
  /* */
}

void StEffHistMangerHelicityGlobal::InitPhiHist()
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

void StEffHistMangerHelicityGlobal::InitKaonHist()
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

void StEffHistMangerHelicityGlobal::InitHist(int bincos = 20, int binphi = 20)
{
  for(int iv2 = 0; iv2 < 1; iv2++)
  {
    for(int irhoinput = 0; irhoinput < 5; irhoinput++)
    {
      for(int irho = 0; irho < 5; irho++)
      {
        for(int ieta = 0; ieta < 11; ieta+=10)
        {
          for(int ipt = 0; ipt < 5; ipt+=4)
          {
//            string HistName = Form("h_mMcEffCosPhiPrime_v2_%d_rhoinput_%d_rho_%d_eta_%d_pt%d",iv2,irhoinput,irho,ieta,ipt);
//            h_mMcEffCosPhiPrime[iv2][irhoinput][irho][ieta][ipt] = new TH2D(HistName.c_str(),HistName.c_str(),20,-1.0,1.0,20,0.0,2.0*TMath::Pi());
//            h_mMcEffCosPhiPrime[iv2][irhoinput][irho][ieta][ipt]->Sumw2();
//            HistName = Form("h_mMcEffCosPhiPrimeH_v2_%d_rhoinput_%d_rho_%d_eta_%d_pt%d",iv2,irhoinput,irho,ieta,ipt);
//            h_mMcEffCosPhiPrimeH[iv2][irhoinput][irho][ieta][ipt] = new TH2D(HistName.c_str(),HistName.c_str(),20,-1.0,1.0,20,0.0,2.0*TMath::Pi());
//            h_mMcEffCosPhiPrimeH[iv2][irhoinput][irho][ieta][ipt]->Sumw2();
//
            for(int icos = 0; icos < 5; icos++)
            {
              for(int iphi = 0; iphi < 5; iphi++) 
              {
            
                string HistName = Form("h3_mMcEffCosPhiPrime_v2_%d_rhoinput_%d_rho_%d_eta_%d_pt%d_cos%d_phi%d",iv2,irhoinput,irho,ieta,ipt,icos+8,iphi+8);
                h3_mMcEffCosPhiPrime[irhoinput][irho][ieta][ipt][icos*5+iphi] = new TH3D(HistName.c_str(),HistName.c_str(),1,0.0,5.0,icos+8,-1.0,1.0,iphi+8,0.0,2.0*TMath::Pi());
                h3_mMcEffCosPhiPrime[irhoinput][irho][ieta][ipt][icos*5+iphi]->Sumw2();
              //HistName = Form("h3_mMcEffCosPhiPrimeH_v2_%d_rhoinput_%d_rho_%d_eta_%d_pt%d",iv2,irhoinput,irho,ieta,ipt);
              //h3_mMcEffCosPhiPrimeH[iv2][irhoinput][irho][ieta][ipt] = new TH3D(HistName.c_str(),HistName.c_str(),1,0.0,5.0,bincos,-1.0,1.0,binphi,0.0,2.0*TMath::Pi());
              //h3_mMcEffCosPhiPrimeH[iv2][irhoinput][irho][ieta][ipt]->Sumw2();
              }
            }
//            string HistName = Form("h3_mMcEffCosPhiPrime_v2_%d_rhoinput_%d_rho_%d_eta_%d_pt%d",iv2,irhoinput,irho,ieta,ipt);
//            h3_mMcEffCosPhiPrime[iv2][irhoinput][irho][ieta][ipt] = new TH3D(HistName.c_str(),HistName.c_str(),1,0.0,5.0,bincos,-1.0,1.0,binphi,0.0,2.0*TMath::Pi());
//            h3_mMcEffCosPhiPrime[iv2][irhoinput][irho][ieta][ipt]->Sumw2();
//            HistName = Form("h3_mMcEffCosPhiPrimeH_v2_%d_rhoinput_%d_rho_%d_eta_%d_pt%d",iv2,irhoinput,irho,ieta,ipt);
//            h3_mMcEffCosPhiPrimeH[iv2][irhoinput][irho][ieta][ipt] = new TH3D(HistName.c_str(),HistName.c_str(),1,0.0,5.0,bincos,-1.0,1.0,binphi,0.0,2.0*TMath::Pi());
//            h3_mMcEffCosPhiPrimeH[iv2][irhoinput][irho][ieta][ipt]->Sumw2();
//
//            HistName = Form("h3_mMcEffCosPhiPrimeY_v2_%d_rhoinput_%d_rho_%d_eta_%d_pt%d",iv2,irhoinput,irho,ieta,ipt);
//            h3_mMcEffCosPhiPrimeY[iv2][irhoinput][irho][ieta][ipt] = new TH3D(HistName.c_str(),HistName.c_str(),10,-1.0,1.0,7,0.0,1.0,binphi,0.0,2.0*TMath::Pi());
//            h3_mMcEffCosPhiPrimeY[iv2][irhoinput][irho][ieta][ipt]->Sumw2();
//            HistName = Form("h3_mMcEffCosPhiPrimeHY_v2_%d_rhoinput_%d_rho_%d_eta_%d_pt%d",iv2,irhoinput,irho,ieta,ipt);
//            h3_mMcEffCosPhiPrimeHY[iv2][irhoinput][irho][ieta][ipt] = new TH3D(HistName.c_str(),HistName.c_str(),10,-1.0,1.0,7,0.0,1.0,binphi,0.0,2.0*TMath::Pi());
//            h3_mMcEffCosPhiPrimeHY[iv2][irhoinput][irho][ieta][ipt]->Sumw2();
       
        //for(int iphipsi = 0; iphipsi < 20; iphipsi++)
        //{
        //  HistName = Form("h_mMcEffCosPhiPrime_v2_%d_rhoinput_%d_rho_%d_phipsi_%d",iv2,irhoinput,irho,iphipsi);
        //  h_mMcEffCosPhiPrimePsi[iv2][irhoinput][irho][iphipsi] = new TH2D(HistName.c_str(),HistName.c_str(),20,-1.0,1.0,20,0.0,2.0*TMath::Pi());
        //  h_mMcEffCosPhiPrimePsi[iv2][irhoinput][irho][iphipsi]->Sumw2();
        //  HistName = Form("h_mMcEffCosPhiPrimeH_v2_%d_rhoinput_%d_rho_%d_phipsi_%d",iv2,irhoinput,irho,iphipsi);
        //  h_mMcEffCosPhiPrimeHPsi[iv2][irhoinput][irho][iphipsi] = new TH2D(HistName.c_str(),HistName.c_str(),20,-1.0,1.0,20,0.0,2.0*TMath::Pi());
        //  h_mMcEffCosPhiPrimeHPsi[iv2][irhoinput][irho][iphipsi]->Sumw2();
        //}
          }
        }
      }
    }
  }
}

void StEffHistMangerHelicityGlobal::FillPhiHistMc(int cent, float pt, float phi, float y, float phistar, float cos, float cosH, float kpt, float ky)
{
  if(cent < vmsa::cent_low[0] || cent > vmsa::cent_up[0]) return;
  //cout << "pt = " << pt << ", y = " << y << endl;
  double  ratio_weight = 1.0;//h_mSpectraRatio->GetBinContent(h_mSpectraRatio->FindBin(pt,y));

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
          double  ratio_cos = 1.0;//h_mCosCosHRatio[phi_pt_bin]->GetBinContent(h_mCosCosHRatio[phi_pt_bin]->FindBin(cos,cosH));
          double  ratio_weight_kaon = 1.0;//h_mKaonSpectraRatio[phi_pt_bin]->GetBinContent(h_mKaonSpectraRatio[phi_pt_bin]->FindBin(kpt,ky));
          for(Int_t i_phi = 0; i_phi < 10; i_phi++)
          {
            //std::cout << "Cos2PhiStarPhi = " << Cos2PhiStarPhi << ", (float(i_phi)-5.)/5. = " << (float(i_phi)-5.)/5. << ", (float(i_phi)-4.)/5. = " << (float(i_phi)-4.)/5 << std::endl;
            if(Cos2PhiStarPhi >= (float(i_phi)-5.)/5. && Cos2PhiStarPhi < (float(i_phi)-4.)/5.)
            {
              TString KEY = Form("mcphi_phipt_%d_cos2phistarphi_%d",i_pt_hist,i_phi);
              h_mPhi_MC[KEY]->Fill(phi,y,ratio_weight*ratio_cos*ratio_weight_kaon);
              KEY = Form("mcphi_phistarphi_phipt_%d",i_pt_hist);
              h_mPhi_MC[KEY]->Fill(phi,phistar,ratio_weight*ratio_cos*ratio_weight_kaon);
              KEY = Form("mcphi_phistarmphi_phipt_%d",i_pt_hist);
              h_mPhi1D_MC[KEY]->Fill(phistar-phi,ratio_weight*ratio_cos*ratio_weight_kaon);
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

void StEffHistMangerHelicityGlobal::FillPhiHistRc(int cent, float pt, float phi, float y, float phistar, float cos, float cosH, float kpt, float ky)
{
  if(cent < vmsa::cent_low[0] || cent > vmsa::cent_up[0]) return;
  double  ratio_weight = 1.0; //h_mSpectraRatio->GetBinContent(h_mSpectraRatio->FindBin(pt,y));

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
          double  ratio_cos = 1.0;//h_mCosCosHRatio[phi_pt_bin]->GetBinContent(h_mCosCosHRatio[phi_pt_bin]->FindBin(cos,cosH));
          double  ratio_weight_kaon = 1.0;//h_mKaonSpectraRatio[phi_pt_bin]->GetBinContent(h_mKaonSpectraRatio[phi_pt_bin]->FindBin(kpt,ky));
          for(Int_t i_phi = 0; i_phi < 10; i_phi++)
          {
            //std::cout << "Cos2PhiStarPhi = " << Cos2PhiStarPhi << ", (float(i_phi)-5.)/5. = " << (float(i_phi)-5.)/5. << ", (float(i_phi)-4.)/5. = " << (float(i_phi)-4.)/5 << std::endl;
            if(Cos2PhiStarPhi >= (float(i_phi)-5.)/5. && Cos2PhiStarPhi < (float(i_phi)-4)/5.)
            {
              TString KEY = Form("rcphi_phipt_%d_cos2phistarphi_%d",i_pt_hist,i_phi);
              h_mPhi_RC[KEY]->Fill(phi,y,ratio_weight*ratio_cos*ratio_weight_kaon);
              KEY = Form("rcphi_phistarphi_phipt_%d",i_pt_hist);
              h_mPhi_RC[KEY]->Fill(phi,phistar,ratio_weight*ratio_cos*ratio_weight_kaon);
              KEY = Form("rcphi_phistarmphi_phipt_%d",i_pt_hist);
              h_mPhi1D_RC[KEY]->Fill(phistar-phi,ratio_weight*ratio_cos*ratio_weight_kaon);
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

void StEffHistMangerHelicityGlobal::FillKaonHistMc(int cent, float phi, float pt, float y, float kpt, float ky, float keta, float kphi, float phistar, float cos, float cosH)
{
  if(cent < vmsa::cent_low[0] || cent > vmsa::cent_up[0]) return;
  double  ratio_weight = 1.0;//h_mSpectraRatio->GetBinContent(h_mSpectraRatio->FindBin(pt,y));

  double Cos2PhiStarPhi = TMath::Cos(2.0*(phistar-phi));

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    //std::cout << "pt = " << pt << ", i_pt = " << i_pt << ", mpt_low[i_pt] = " << mpt_low[i_pt] << ", mpt_up[i_pt] = " << mpt_up[i_pt] << std::endl;
    if(pt >= mpt_low[i_pt] && pt < mpt_up[i_pt])
    {
      int phi_pt_bin = i_pt - 2;
      double  ratio_cos = 1.0;//h_mCosCosHRatio[phi_pt_bin]->GetBinContent(h_mCosCosHRatio[phi_pt_bin]->FindBin(cos,cosH));
      double  ratio_weight_kaon = 1.0;//h_mKaonSpectraRatio[phi_pt_bin]->GetBinContent(h_mKaonSpectraRatio[phi_pt_bin]->FindBin(kpt,ky));
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
              h_mKaon_MC[KEY]->Fill(kphi,ky,ratio_weight*ratio_cos*ratio_weight_kaon);
              KEY = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",phi_pt_bin,i_phi,i_pt_kaon);
              //std::cout << "Filling with phi = " << kphi << std::endl;
              h_mKaon_MC[KEY]->Fill(kphi,keta,ratio_weight*ratio_cos*ratio_weight_kaon);
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

void StEffHistMangerHelicityGlobal::FillKaonDeltaHistMc(int cent, float phi, float eta, float pt, float y, float kppt, float kpeta, float kpphi, float kmpt, float kmeta, float kmphi, float phistar)
{
  if(cent < vmsa::cent_low[0] || cent > vmsa::cent_up[0]) return;
  double  ratio_weight = 1.0;//h_mSpectraRatio->GetBinContent(h_mSpectraRatio->FindBin(pt,y));

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
          h_mKaon_MC[KEY]->Fill(phidiff,etadiff,ratio_weight);
          KEY = Form("mckaon_deltarphi_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_MC[KEY]->Fill(phi_r,deltar,ratio_weight);
          KEY = Form("mckaonplusphi_delta_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_MC[KEY]->Fill(phiphidiffp,phietadiffp,ratio_weight);
          KEY = Form("mckaonminusphi_delta_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_MC[KEY]->Fill(phiphidiffm,phietadiffm,ratio_weight);
          break;
        }
      }
      break;
    }
  }
}

void StEffHistMangerHelicityGlobal::FillKaonHistRc(int cent, float phi, float pt, float y, float kpt, float ky, float keta, float kphi, float phistar, float cos, float cosH)
{
  if(cent < vmsa::cent_low[0] || cent > vmsa::cent_up[0]) return;
  double  ratio_weight = 1.0; //h_mSpectraRatio->GetBinContent(h_mSpectraRatio->FindBin(pt,y));

  double Cos2PhiStarPhi = TMath::Cos(2.0*(phistar-phi));

  for(int i_pt = mpt_first; i_pt <= mpt_last; ++i_pt) // use rebinned pt
  {
    if(pt >= mpt_low[i_pt] && pt < mpt_up[i_pt])
    {
      int phi_pt_bin = i_pt - 2;
      double  ratio_cos = 1.0;//h_mCosCosHRatio[phi_pt_bin]->GetBinContent(h_mCosCosHRatio[phi_pt_bin]->FindBin(cos,cosH));
      double  ratio_weight_kaon = 1.0;//h_mKaonSpectraRatio[phi_pt_bin]->GetBinContent(h_mKaonSpectraRatio[phi_pt_bin]->FindBin(kpt,ky));
      for(Int_t i_phi = 0; i_phi < 10; i_phi++)
      {
        if(Cos2PhiStarPhi >= (float(i_phi)-5.)/5. && Cos2PhiStarPhi < (float(i_phi)-4)/5.)
        {
          for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
          {
            if(kpt >= data_constants::kaon_pt_low[i_pt_kaon] && kpt < data_constants::kaon_pt_high[i_pt_kaon])   
            {
              TString KEY = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",phi_pt_bin,i_phi,i_pt_kaon);
              h_mKaon_RC[KEY]->Fill(kphi,ky,ratio_weight*ratio_cos*ratio_weight_kaon);
              KEY = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",phi_pt_bin,i_phi,i_pt_kaon);
              h_mKaon_RC[KEY]->Fill(kphi,keta,ratio_weight*ratio_cos*ratio_weight_kaon);
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
void StEffHistMangerHelicityGlobal::FillKaonDeltaHistRc(int cent, float phi, float eta, float pt, float y, float kppt, float kpeta, float kpphi, float kmpt, float kmeta, float kmphi, float phistar)
{
  if(cent < vmsa::cent_low[0] || cent > vmsa::cent_up[0]) return;
  double  ratio_weight =1.0;// h_mSpectraRatio->GetBinContent(h_mSpectraRatio->FindBin(pt,y));

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
          h_mKaon_RC[KEY]->Fill(phidiff,etadiff,ratio_weight);
          KEY = Form("rckaon_deltarphi_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_RC[KEY]->Fill(phi_r,deltar,ratio_weight);
          KEY = Form("rckaonplusphi_delta_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_RC[KEY]->Fill(phiphidiffp,phietadiffp,ratio_weight);
          KEY = Form("rckaonminusphi_delta_phipt_%d_cos2phistarphi_%d",phi_pt_bin,i_phi);
          h_mKaon_RC[KEY]->Fill(phiphidiffm,phietadiffm,ratio_weight);
          break;
        }
      }
      break;
    }
  }
}

void StEffHistMangerHelicityGlobal::WritePhiHist()
{
  for(Int_t i_pt_phi = mpt_first-2; i_pt_phi <= mpt_last-2; ++i_pt_phi) // use rebinned pt
  {
    for(Int_t i_pt = data_constants::phi_pt_start[i_pt_phi]; i_pt <= data_constants::phi_pt_stop[i_pt_phi]; i_pt++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++)
      {
        TString KEY_MC = Form("mcphi_phipt_%d_cos2phistarphi_%d",i_pt,i_phi);
        h_mPhi_MC[KEY_MC]->Write();
        //std::cout << "WROTE: " << KEY_MC << std::endl;
        TString KEY_RC = Form("rcphi_phipt_%d_cos2phistarphi_%d",i_pt,i_phi);
        h_mPhi_RC[KEY_RC]->Write();
        //std::cout << "WROTE: " << KEY_RC << std::endl;
      }

      TString KEY_MC = Form("mcphi_phistarphi_phipt_%d",i_pt);
      h_mPhi_MC[KEY_MC]->Write();
      //std::cout << "WROTE: " << KEY_MC << std::endl;
      TString KEY_RC = Form("rcphi_phistarphi_phipt_%d",i_pt);
      h_mPhi_RC[KEY_RC]->Write();
      //std::cout << "WROTE: " << KEY_RC << std::endl;

      KEY_MC = Form("mcphi_phistarmphi_phipt_%d",i_pt);
      h_mPhi1D_MC[KEY_MC]->Write();
      //std::cout << "WROTE: " << KEY_MC << std::endl;
      KEY_RC = Form("rcphi_phistarmphi_phipt_%d",i_pt);
      h_mPhi1D_RC[KEY_RC]->Write();
      //std::cout << "WROTE: " << KEY_RC << std::endl;
    }
  }
}
void StEffHistMangerHelicityGlobal::WriteKaonHist()
{
  for(Int_t i_pt_phi = mpt_first-2; i_pt_phi <= mpt_last-2; ++i_pt_phi) // use rebinned pt
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++)
    {
      for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
      {
        TString KEY_MC = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        h_mKaon_MC[KEY_MC]->Write();
        //std::cout << "WROTE: " << KEY_MC << std::endl;
        TString KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        h_mKaon_RC[KEY_RC]->Write();
        //std::cout << "WROTE: " << KEY_RC << std::endl;

        KEY_MC = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);
        h_mKaon_MC[KEY_MC]->Write();
        //std::cout << "WROTE: " << KEY_MC << std::endl;
        KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);
        h_mKaon_RC[KEY_RC]->Write();
        //std::cout << "WROTE: " << KEY_RC << std::endl;
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

void StEffHistMangerHelicityGlobal::FillHistMc(int cent, float pt, float y, float phi, float cos, float Psi, float phistar, float phiprime, float kpt, float ky, float cosH, float phiprimeH, float weight, int iv2, int irhoinput, int irho, int ieta, int ipt)
{
  
  double weightv2 = 1.0;//h_mV2[iv2]->GetBinContent(cosH,phiprimeH);
  //h_mMcEffCosPhiPrime[iv2][irhoinput][irho][ieta][ipt]->Fill(cos,phiprime,weight*weightv2);
  //h_mMcEffCosPhiPrimeH[iv2][irhoinput][irho][ieta][ipt]->Fill(cosH,phiprimeH,weight*weightv2);
  for(int icos = 0; icos < 5; icos++)
  {
    for(int iphi = 0; iphi < 5; iphi++) 
    {
      h3_mMcEffCosPhiPrime[irhoinput][irho][ieta][ipt][icos*5+iphi]->Fill(pt,cos,phiprime,weight*weightv2);
      //h3_mMcEffCosPhiPrimeH[iv2][irhoinput][irho][ieta][ipt]->Fill(pt,cosH,phiprimeH,weight*weightv2);
    }
  }
//  h3_mMcEffCosPhiPrime[iv2][irhoinput][irho][ieta][ipt]->Fill(pt,cos,phiprime,weight*weightv2);
//  h3_mMcEffCosPhiPrimeH[iv2][irhoinput][irho][ieta][ipt]->Fill(pt,cosH,phiprimeH,weight*weightv2);
//
//  h3_mMcEffCosPhiPrimeY[iv2][irhoinput][irho][ieta][ipt]->Fill(y,TMath::Abs(cos),phiprime,weight*weightv2);
//  h3_mMcEffCosPhiPrimeHY[iv2][irhoinput][irho][ieta][ipt]->Fill(y,TMath::Abs(cosH),phiprimeH,weight*weightv2);

  //for(int iphipsi = 0; iphipsi < 20; iphipsi++)
  //{
  //  if( phi > double(iphipsi)*TMath::Pi()/20. && phi < double(iphipsi+1)*TMath::Pi()/20.)
  //  {
  //    h_mMcEffCosPhiPrimePsi[iv2][irhoinput][irho][iphipsi]->Fill(cos,phiprime,weight);
  //    h_mMcEffCosPhiPrimeHPsi[iv2][irhoinput][irho][iphipsi]->Fill(cosH,phiprimeH,weight);
  //    break;
  //  }
  //}
}

void StEffHistMangerHelicityGlobal::FillHistRc(int cent, float pt, float y, float phi, float cos, float Psi, float phistar, float phiprime, float kpt, float ky, float cosH)
{
  double  ratio_weight = 1.0;//h_mSpectraRatio->GetBinContent(h_mSpectraRatio->FindBin(pt,y));
  //cout << "pt = " << pt << ", y = " << y << " weight = " << ratio_weight << endl;
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
        double  ratio_cos = 1.0; //h_mCosCosHRatio[i_pt-2]->GetBinContent(h_mCosCosHRatio[i_pt-2]->FindBin(cos,cosH));
        double  ratio_weight_kaon = 1.0;//h_mKaonSpectraRatio[i_pt-2]->GetBinContent(h_mKaonSpectraRatio[i_pt-2]->FindBin(kpt,ky));
        //cout << "kpt = " << kpt << ", ky = " << ky << " weight = " << ratio_weight_kaon << endl << endl;
        h_mRcEffCos[cent][i_pt]->Fill(cos);
        h_mRcEffCosH[cent][i_pt]->Fill(cosH);
        h_mRcEffCosCosH[cent][i_pt]->Fill(cos,cosH);
        h_mRcEffPhiS[cent][i_pt]->Fill(phidiff);
        h_mRcEffCosEP[cent][i_pt]->Fill(cos,phiPsi);
        h_mRcEffCosPhiPrime[cent][i_pt]->Fill(cos,phiprime);
      //  h_mRcCosEP[cent][i_pt]->Fill(cos,phi_shift);
        if(cent >= vmsa::cent_low[0] && cent <= vmsa::cent_up[0])
        {//20-60%
          h_mRcEffCos[9][i_pt]->Fill(cos,vmsa::weight[cent]*ratio_weight*ratio_weight_kaon*ratio_cos);
          h_mRcEffCosH[9][i_pt]->Fill(cosH,vmsa::weight[cent]*ratio_weight*ratio_weight_kaon*ratio_cos);
          h_mRcEffCosCosH[9][i_pt]->Fill(cos,cosH,vmsa::weight[cent]*ratio_weight*ratio_weight_kaon*ratio_cos);
          h_mRcEffPhiS[9][i_pt]->Fill(phidiff,vmsa::weight[cent]*ratio_weight*ratio_weight_kaon*ratio_cos);
          h_mRcEffCosEP[9][i_pt]->Fill(cos,phiPsi,vmsa::weight[cent]*ratio_weight*ratio_weight_kaon*ratio_cos);
          h_mRcEffCosPhiPrime[9][i_pt]->Fill(cos,phiprime,ratio_weight*ratio_weight_kaon*ratio_cos);
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

float StEffHistMangerHelicityGlobal::AngleShift(float phi)
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

TH1D* StEffHistMangerHelicityGlobal::CalEffError(TH1D *h_Mc, TH1D *h_Rc, std::string HistName)
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

TH2D* StEffHistMangerHelicityGlobal::CalEffError(TH2D *h_Mc, TH2D *h_Rc, std::string HistName)
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

void StEffHistMangerHelicityGlobal::CalEfficiency()
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

void StEffHistMangerHelicityGlobal::CalEffPtEtaPhi()
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

void StEffHistMangerHelicityGlobal::CalEffCosThetaStar()
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
        HistName = Form("h_mEffCosH_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCosH[i_cent][i_pt] = CalEffError(h_mMcEffCosH[i_cent][i_pt],h_mRcEffCosH[i_cent][i_pt],HistName.c_str());
        HistName = Form("h_mEffCosCosH_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCosCosH[i_cent][i_pt] = CalEffError(h_mMcEffCosCosH[i_cent][i_pt],h_mRcEffCosCosH[i_cent][i_pt],HistName.c_str());
        HistName = Form("h_mEffPhiS_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffPhiS[i_cent][i_pt] = CalEffError(h_mMcEffPhiS[i_cent][i_pt],h_mRcEffPhiS[i_cent][i_pt],HistName.c_str());
        HistName = Form("h_mEffCosEP_Cent_%d_Pt_%d",i_cent,i_pt);
        h_mEffCosEP[i_cent][i_pt] = CalEffError(h_mMcEffCosEP[i_cent][i_pt],h_mRcEffCosEP[i_cent][i_pt],HistName.c_str());
        //HistName = Form("h_mEffCosPhiPrime_Cent_%d_Pt_%d",i_cent,i_pt);
        //h_mEffCosPhiPrime[i_cent][i_pt] = CalEffError(h_mMcEffCosPhiPrime[i_cent][i_pt],h_mRcEffCosPhiPrime[i_cent][i_pt],HistName.c_str());
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

void StEffHistMangerHelicityGlobal::WriteHist()
{
  for(int iv2 = 0; iv2 < 1; iv2++)
  {
    for(int irhoinput = 0; irhoinput < 5; irhoinput++)
    {
      for(int irho = 0; irho < 5; irho++)
      {
        for(int ieta = 0; ieta < 11; ieta+=10)
        {
          for(int ipt = 0; ipt < 5; ipt+=4)
          {
            //h_mMcEffCosPhiPrime[iv2][irhoinput][irho][ieta][ipt]->Write();
            //h_mMcEffCosPhiPrimeH[iv2][irhoinput][irho][ieta][ipt]->Write();
            for(int icos = 0; icos < 5; icos++)
            {
              for(int iphi = 0; iphi < 5; iphi++) 
              {
                h3_mMcEffCosPhiPrime[irhoinput][irho][ieta][ipt][icos*5+iphi]->Write();
              }
            }
//            h3_mMcEffCosPhiPrime[iv2][irhoinput][irho][ieta][ipt]->Write();
//            h3_mMcEffCosPhiPrimeH[iv2][irhoinput][irho][ieta][ipt]->Write();
//            h3_mMcEffCosPhiPrimeY[iv2][irhoinput][irho][ieta][ipt]->Write();
//            h3_mMcEffCosPhiPrimeHY[iv2][irhoinput][irho][ieta][ipt]->Write();
        //for(int iphipsi = 0; iphipsi < 20; iphipsi++)
        //{
        //  h_mMcEffCosPhiPrimePsi[iv2][irhoinput][irho][iphipsi]->Write();
        //  h_mMcEffCosPhiPrimeHPsi[iv2][irhoinput][irho][iphipsi]->Write();
        //}
          }
        }
      }
    }
  }
}
