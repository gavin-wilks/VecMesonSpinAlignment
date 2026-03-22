#include "StRoot/StVecMesonAna/StVecMesonHistoManger.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "TProfile.h"

ClassImp(StVecMesonHistoManger)

//-------------------------------------------------------------
StVecMesonHistoManger::StVecMesonHistoManger(Int_t EP_mode, Int_t study)
{
  if(EP_mode == 1) EP_string = "1st";
  if(EP_mode == 2) EP_string = "2nd";

  mStudy = study;
}

StVecMesonHistoManger::~StVecMesonHistoManger()
{
}
//-------------------------------------------------------------

void StVecMesonHistoManger::InitSys(Int_t X_flag, Int_t mode) // 0 for Same Event, 1 for Mixed Event
{
  // spin alignment analysis
  if(mStudy != 2)
  {
    for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin 
    {
      for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
      {
        for(Int_t i_phi = 0; i_phi < 12/*2*vmsa::CTS_total*/; i_phi++) // phi-psi bin
        {
          for(Int_t i_thetah = 0; i_thetah < 9/*2*vmsa::CTS_total*/; i_thetah++) // phi-psi bin
          {
            for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
            {
              for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
              {
                if( i_dca != 0 && i_sig != 0) continue;
                TString Mode[2] = {"SE","ME"};
                TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_cent,i_thetah,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                if(mStudy == 0)
                {
                  h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),100,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
                  h_mMass2[KEY_Mass2]->Sumw2();
                } 
              }
            }
          }
        }
      }
    }
    for(Int_t i_pt = 0; i_pt < vmsa::pty_total; i_pt++) // pt bin 
    {
      for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
      {
        for(Int_t i_phi = 0; i_phi < 12/*2*vmsa::CTS_total*/; i_phi++) // phi-psi bin
        {
          for(Int_t i_thetah = 0; i_thetah < 9/*2*vmsa::CTS_total*/; i_thetah++) // phi-psi bin
          {
            for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
            {
              for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
              {
                if( i_dca != 0 && i_sig != 0) continue;
                TString Mode[2] = {"SE","ME"};
                TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_cent,i_thetah,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                if(mStudy == 1)
                {
                  for(int i_eta = 0; i_eta < vmsa::eta_total; i_eta++)
                  {
                    KEY_Mass2 = Form("pt_%d_eta_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_eta,i_cent,i_thetah,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                    h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),50,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
                    h_mMass2[KEY_Mass2]->Sumw2();
                  }
                }
                if(mStudy == 3) 
                {
                  for(int i_phipsi = 0; i_phipsi < 12; i_phipsi++)
                  {
                    KEY_Mass2 = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_phipsi,i_cent,i_thetah,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                    //cout << KEY_Mass2 << endl;
                    h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),50,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
                    h_mMass2[KEY_Mass2]->Sumw2();
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if(mStudy == 2)
  {
    for(Int_t i_p = 0; i_p < vmsa::p_total; i_p++) // pt bin 
    {
      for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
      {
        for(Int_t i_phi = 0; i_phi < 12/*2*vmsa::CTS_total*/; i_phi++) // phi-psi bin
        {
          for(Int_t i_thetah = 0; i_thetah < 9/*2*vmsa::CTS_total*/; i_thetah++) // phi-psi bin
          {
            for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
            {
              for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
              {
                if( i_dca != 0 && i_sig != 0) continue;
                TString Mode[2] = {"SE","ME"};
                TString KEY_Mass2 = Form("p_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_p,i_cent,i_thetah,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),100,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
                h_mMass2[KEY_Mass2]->Sumw2();
              }
            }
          }
        }
      }
    }
  }
  // spin alignment analysis
  //for(Int_t i_eta = 0; i_eta < 10; i_eta++) // pt bin 
  //{
  //  for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
  //  {
  //    for(Int_t i_phi = 0; i_phi < 12/*2*vmsa::CTS_total*/; i_phi++) // phi-psi bin
  //    {
  //      for(Int_t i_thetah = 0; i_thetah < 9/*2*vmsa::CTS_total*/; i_thetah++) // phi-psi bin
  //      {
  //        for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  //        {
  //          for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
  //          {
  //            if( i_dca != 0 && i_sig != 0) continue;
  //            TString Mode[2] = {"SE","ME"};
  //            TString KEY_Mass2 = Form("eta_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_%s",i_eta,i_cent,i_thetah,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
  //            h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
  //            h_mMass2[KEY_Mass2]->Sumw2();
  //          }
  //        }
  //      }
  //    }
  //  }
  //}
}
//-------------------------------------------------------------
void StVecMesonHistoManger::FillSys(Float_t pt, Float_t eta, Float_t p, Int_t cent9, Float_t CosThetaStar, Float_t CosThetaStarH, Float_t phihelicity, Float_t PhiPsi, Float_t Cos2PhiStarPhi, Int_t dcaSys, Int_t nSigSys, Float_t Res2, Float_t InvMass, Double_t reweight, Int_t X_flag, Int_t mode)
{
  TString Mode[2] = {"SE","ME"};

  if(mStudy != 2)
  {
    if(mStudy == 0)
    {
      for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt_bin
      {
        if(pt >= vmsa::ptRawStart[i_pt] && pt < vmsa::ptRawStop[i_pt])
        {
          for(Int_t i_phi = 0; i_phi < 12/*2*vmsa::CTS_total*/; i_phi++) // phi-psi2 bin
          {
            if(phihelicity >= TMath::Pi()*float(i_phi)/6.0 && phihelicity < TMath::Pi()*float(i_phi+1)/6.0)
            {
              for(Int_t i_thetah = 0; i_thetah < 9; i_thetah++)
              {
                if(CosThetaStarH >= float(i_thetah-4.5)/4.5 && CosThetaStarH < float(i_thetah-3.5)/4.5)
                {
                  TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,cent9,i_thetah,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                  h_mMass2[KEY_Mass2]->Fill(InvMass,reweight);
                  if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
                  {
                    TString KEY_Mass2Sys = Form("pt_%d_Centrality_9_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_thetah,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                    h_mMass2[KEY_Mass2Sys]->Fill(InvMass,reweight);
                  }
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
    else
    {
      for(Int_t i_pt = 0; i_pt < vmsa::pty_total; i_pt++) // pt_bin
      {
        if(pt >= vmsa::ptRawStart[i_pt] && pt < vmsa::ptRawStop[i_pt])
        {
          for(Int_t i_phi = 0; i_phi < 12/*2*vmsa::CTS_total*/; i_phi++) // phi-psi2 bin
          {
            if(phihelicity >= TMath::Pi()*float(i_phi)/6.0 && phihelicity < TMath::Pi()*float(i_phi+1)/6.0)
            {
              for(Int_t i_thetah = 0; i_thetah < 9; i_thetah++)
              {
                if(CosThetaStarH >= float(i_thetah-4.5)/4.5 && CosThetaStarH < float(i_thetah-3.5)/4.5)
                {
                  if(mStudy == 1)
                  {
                    for(int i_eta = 0; i_eta < vmsa::eta_total; i_eta++)
                    {
                      if(eta >= vmsa::eta_raw[i_eta] && eta < vmsa::eta_raw[i_eta+1])
                      {
                        TString KEY_Mass2 = Form("pt_%d_eta_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_eta,cent9,i_thetah,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                        h_mMass2[KEY_Mass2]->Fill(InvMass,reweight);
                        if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
                        {
                          TString KEY_Mass2Sys = Form("pt_%d_eta_%d_Centrality_9_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_eta,i_thetah,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                          h_mMass2[KEY_Mass2Sys]->Fill(InvMass,reweight);
                        }
                        break;
                      }
                    }
                  }
                  if(mStudy == 3) 
                  {
                    for(Int_t i_phipsi = 0; i_phipsi < 12; i_phipsi++)
                    {
                      if(PhiPsi >= float(i_phipsi)*TMath::Pi()/6. && PhiPsi < float(i_phipsi+1)*TMath::Pi()/6.)
                      {
                        TString KEY_Mass2 = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_phipsi,cent9,i_thetah,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                        h_mMass2[KEY_Mass2]->Fill(InvMass,reweight);
                        if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
                        {
                          TString KEY_Mass2Sys = Form("pt_%d_phipsi_%d_Centrality_9_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_phipsi,i_thetah,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                          h_mMass2[KEY_Mass2Sys]->Fill(InvMass,reweight);
                        }
                        break;
                      }
                    }
                  }
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
  }
   
  if(mStudy == 2)
  {
    for(Int_t i_p = 0; i_p < vmsa::p_total; i_p++) // pt_bin
    {
      if(p >= vmsa::pRawStart[i_p] && p < vmsa::pRawStop[i_p])
      {
        for(Int_t i_phi = 0; i_phi < 12/*2*vmsa::CTS_total*/; i_phi++) // phi-psi2 bin
        {
          if(phihelicity >= TMath::Pi()*float(i_phi)/6.0 && phihelicity < TMath::Pi()*float(i_phi+1)/6.0)
          {
            for(Int_t i_thetah = 0; i_thetah < 9; i_thetah++)
            {
              if(CosThetaStarH >= float(i_thetah-4.5)/4.5 && CosThetaStarH < float(i_thetah-3.5)/4.5)
              {
                TString KEY_Mass2 = Form("p_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_p,cent9,i_thetah,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                h_mMass2[KEY_Mass2]->Fill(InvMass,reweight);
                if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
                {
                  TString KEY_Mass2Sys = Form("p_%d_Centrality_9_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_p,i_thetah,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                  h_mMass2[KEY_Mass2Sys]->Fill(InvMass,reweight);
                }
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
  //if(pt < 1.2 || pt >= 4.2) return; 

  //for(Int_t i_eta = 0; i_eta < 10; i_eta++) // pt_bin
  //{
  //  if(eta >= float(i_eta-5)/5. && eta < float(i_eta-4)/5.)
  //  {
  //    for(Int_t i_phi = 0; i_phi < 12/*2*vmsa::CTS_total*/; i_phi++) // phi-psi2 bin
  //    {
  //      if(phihelicity >= TMath::Pi()*float(i_phi)/6.0 && phihelicity < TMath::Pi()*float(i_phi+1)/6.0)
  //      {
  //        for(Int_t i_thetah = 0; i_thetah < 9; i_thetah++)
  //        {
  //          if(CosThetaStarH >= float(float(i_thetah)-4.5)/4.5 && CosThetaStarH < float(float(i_thetah)-3.5)/4.5)
  //          {
  //            TString KEY_Mass2 = Form("eta_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_%s",i_eta,cent9,i_thetah,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
  //            h_mMass2[KEY_Mass2]->Fill(InvMass,reweight);
  //            if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
  //            {
  //              TString KEY_Mass2Sys = Form("eta_%d_Centrality_9_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_%s",i_eta,i_thetah,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
  //              h_mMass2[KEY_Mass2Sys]->Fill(InvMass,reweight);
  //            }
  //            break;
  //          }
  //        }
  //        break;
  //      }
  //    }
  //  }	  
  //}
}

void StVecMesonHistoManger::FillDcaSys(Float_t dcaA, Float_t dcaB, Int_t dcaSys)
{
  TString KEY_DcaA = Form("Tracks_DcaA_%d",dcaSys);
  h_mDca[KEY_DcaA]->Fill(dcaA);
  TString KEY_DcaB = Form("Tracks_DcaB_%d",dcaSys);
  h_mDca[KEY_DcaB]->Fill(dcaB);
}

void StVecMesonHistoManger::FillSigSys(Float_t nsA, Float_t nsB, Int_t nSigSys)
{
  TString KEY_nSigKaonA = Form("TracksSigA_%d",nSigSys);
  h_mSigKaon[KEY_nSigKaonA]->Fill(nsA); 
  TString KEY_nSigKaonB = Form("TracksSigB_%d",nSigSys);
  h_mSigKaon[KEY_nSigKaonB]->Fill(nsB); 
}

//-------------------------------------------------------------
void StVecMesonHistoManger::WriteSys(Int_t X_flag, Int_t mode)
{
  TString Mode[2] = {"SE","ME"};
  cout << "About to Write Histograms" << endl;
 
  if(mStudy != 2)
  {
    if(mStudy == 0)
    {
      for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin
      {
        for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
        {
          for(Int_t i_phi = 0; i_phi < 12/*2*vmsa::CTS_total*/; i_phi ++) // cos(theta*) bin
          {
            for(Int_t i_thetah = 0; i_thetah < 9/*2*vmsa::CTS_total*/; i_thetah ++) // cos(theta*) bin
            {
              for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
              {
                for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
                {
                  if( i_dca != 0 && i_sig != 0) continue;
                  TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_cent,i_thetah,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                  h_mMass2[KEY_Mass2]->Write();
                }
              }
            }
          }
        }
      }
    }
    else 
    {
      for(Int_t i_pt = 0; i_pt < vmsa::pty_total; i_pt++) // pt bin
      {
        for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
        {
          for(Int_t i_phi = 0; i_phi < 12/*2*vmsa::CTS_total*/; i_phi ++) // cos(theta*) bin
          {
            for(Int_t i_thetah = 0; i_thetah < 9/*2*vmsa::CTS_total*/; i_thetah ++) // cos(theta*) bin
            {
              for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
              {
                for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
                {
                  if( i_dca != 0 && i_sig != 0) continue;
                  if(mStudy == 1)
                  {
                    for(int i_eta = 0; i_eta < vmsa::eta_total; i_eta++)
                    {
                      TString KEY_Mass2 = Form("pt_%d_eta_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_eta,i_cent,i_thetah,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                      h_mMass2[KEY_Mass2]->Write();
                    }
                  }
                  if(mStudy == 3)
                  {
                    for(Int_t i_phipsi = 0; i_phipsi < 12; i_phipsi++)
                    {
                      TString KEY_Mass2 = Form("pt_%d_phipsi_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_phipsi,i_cent,i_thetah,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                      h_mMass2[KEY_Mass2]->Write();
                    }
                  }
                  //cout << "Wrote: " << KEY_Mass2 << endl;
                }
              }
            }
          }
        }
      }
    }
  }
 
  if(mStudy == 2)
  {
    for(Int_t i_p = 0; i_p < vmsa::p_total; i_p++) // pt bin
    {
      for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
      {
        for(Int_t i_phi = 0; i_phi < 12/*2*vmsa::CTS_total*/; i_phi ++) // cos(theta*) bin
        {
          for(Int_t i_thetah = 0; i_thetah < 9/*2*vmsa::CTS_total*/; i_thetah ++) // cos(theta*) bin
          {
            for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
            {
              for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
              {
                if( i_dca != 0 && i_sig != 0) continue;
                TString KEY_Mass2 = Form("p_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_p,i_cent,i_thetah,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                h_mMass2[KEY_Mass2]->Write();
                //cout << "Wrote: " << KEY_Mass2 << endl;
              }
            }
          }
        }
      }
    }
  }
  //for(Int_t i_eta = 0; i_eta < 10; i_eta++) // pt bin
  //{
  //  for(Int_t i_cent = 9/*vmsa::Cent_start*/; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
  //  {
  //    for(Int_t i_phi = 0; i_phi < 12/*2*vmsa::CTS_total*/; i_phi ++) // cos(theta*) bin
  //    {
  //      for(Int_t i_thetah = 0; i_thetah < 9/*2*vmsa::CTS_total*/; i_thetah ++) // cos(theta*) bin
  //      {
  //        for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  //        {
  //          for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
  //          {
  //            if( i_dca != 0 && i_sig != 0) continue;
  //            TString KEY_Mass2 = Form("eta_%d_Centrality_%d_CosThetaStar_%d_Phi_%d_%s_Dca_%d_Sig_%d_%s_%s",i_eta,i_cent,i_thetah,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
  //            h_mMass2[KEY_Mass2]->Write();
  //            cout << "Wrote: " << KEY_Mass2 << endl;
  //          }
  //        }
  //      }
  //    }
  //  }
  //}
  //cout << "Write invmass" << endl;
}

//-----------------------EP Direction----------------------------

void StVecMesonHistoManger::InitSys_EP(Int_t X_flag, Int_t mode) // 0 for Same Event, 1 for Mixed Event
{
  // spin alignment analysis
  for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin 
  {
    for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_theta = 0; i_theta < vmsa::CTS_total; i_theta++) // phi-psi bin
      {
	for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
	{
	  for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	  {
	    TString Mode[2] = {"SE","ME"};
	    TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_%s_EP",i_pt,i_cent,i_theta,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	    h_mMass2_EP[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
	    h_mMass2_EP[KEY_Mass2]->Sumw2();
	  }
	}
      }
    }
  }
}
//-------------------------------------------------------------
void StVecMesonHistoManger::FillSys_EP(Float_t pt, Int_t cent9, Float_t CosThetaStar, Int_t dcaSys, Int_t nSigSys, Float_t Res2, Float_t InvMass, Double_t reweight, Int_t X_flag, Int_t mode)
{
  TString Mode[2] = {"SE","ME"};
  if(Res2 > 0.0)
  {
    for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt_bin
    {
      if(pt >= vmsa::ptRawStart[i_pt] && pt < vmsa::ptRawStop[i_pt])
      {
	for(Int_t i_theta = 0; i_theta < vmsa::CTS_total; i_theta++) // phi-psi2 bin
	{
	  if(TMath::Abs(CosThetaStar) >= vmsa::CTS_low[i_theta] && TMath::Abs(CosThetaStar) < vmsa::CTS_up[i_theta])
	  {
	    // spin alignment
	    TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_%s_EP",i_pt,cent9,i_theta,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	    h_mMass2_EP[KEY_Mass2]->Fill(InvMass,reweight);
	    if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
	    {
	      TString KEY_Mass2Sys = Form("pt_%d_Centrality_9_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_%s_EP",i_pt,i_theta,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	      h_mMass2_EP[KEY_Mass2Sys]->Fill(InvMass,reweight);
	    }
	  }
	}
      }
    }
  }
}

//-------------------------------------------------------------
void StVecMesonHistoManger::WriteSys_EP(Int_t X_flag, Int_t mode)
{
  TString Mode[2] = {"SE","ME"};
  // flow
  for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin
  {
    for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_theta = 0; i_theta < vmsa::CTS_total; i_theta ++) // cos(theta*) bin
      {
	for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
	{
	  for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	  {
	    TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_%s_Dca_%d_Sig_%d_%s_%s_EP",i_pt,i_cent,i_theta,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	    h_mMass2_EP[KEY_Mass2]->Write();
	  }
	}
      }
    }
  }
}
