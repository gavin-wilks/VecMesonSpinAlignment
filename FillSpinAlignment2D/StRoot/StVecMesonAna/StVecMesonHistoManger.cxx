#include "StRoot/StVecMesonAna/StVecMesonHistoManger.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "TProfile.h"

ClassImp(StVecMesonHistoManger)

//-------------------------------------------------------------
StVecMesonHistoManger::StVecMesonHistoManger(Int_t EP_mode)
{
  if(EP_mode == 1) EP_string = "1st";
  if(EP_mode == 2) EP_string = "2nd";
}

StVecMesonHistoManger::~StVecMesonHistoManger()
{
}
//-------------------------------------------------------------

void StVecMesonHistoManger::InitSys(Int_t X_flag, Int_t mode) // 0 for Same Event, 1 for Mixed Event
{
  // for spectra 
  for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin 
  {
    for(Int_t i_eta = 0; i_eta < 16; i_eta++)
    {
      //for(Int_t i_phi = 0; i_phi < 24; i_phi++)
      //{
        for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
        {
          for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
          {
            for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
            {
              if( i_dca != 0 && i_sig != 0) continue;
              TString Mode[2] = {"SE","ME"};
              
              //TString KEY_Mass2 = Form("pt_%d_eta_%d_phi_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_eta,i_phi,i_cent,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
              //h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),100,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
              //h_mMass2[KEY_Mass2]->Sumw2();
              TString KEY_Mass2 = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_eta,i_cent,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
              h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),100,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
              h_mMass2[KEY_Mass2]->Sumw2();
            }
          }
        }
      //}
    }
  }

  // spin alignment analysis
  //for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin 
  //{
  //  for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
  //  {
  //    for(Int_t i_theta = 0; i_theta < 10; i_theta++) // cos(theta*) bins
  //    {
  //      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // cos(theta*) bins
  //      {
  //        for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  //        {
  //          for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
  //          {
  //            if( i_dca != 0 && i_sig != 0) continue;
  //            TString Mode[2] = {"SE","ME"};
  //            
  //            TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_cent,i_theta,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
  //            h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),100,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
  //            h_mMass2[KEY_Mass2]->Sumw2();
  //          }
  //        }
  //      }
  //    }
  //  }
  //}
}
//-------------------------------------------------------------
void StVecMesonHistoManger::FillSys(Float_t pt, Float_t eta, Float_t phi, Int_t cent9, Float_t CosThetaStar, Float_t PhiPsi, Float_t Cos2PhiStarPhi, Float_t phiprime, Int_t dcaSys, Int_t nSigSys, Float_t Res2, Float_t InvMass, Double_t reweight, Int_t X_flag, Int_t mode)
{
  TString Mode[2] = {"SE","ME"};
  if(Res2 > 0.0)
  {
    for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt_bin
    {
      if(pt >= vmsa::ptRawStart[i_pt] && pt < vmsa::ptRawStop[i_pt])
      {
        for(Int_t i_eta = 0; i_eta < 16; i_eta++)
        {
          if(eta >= float(i_eta-8.)/8. && eta < float(i_eta-7.)/8.)
          {
            //for(Int_t i_phi = 0; i_phi < 24; i_phi++)
            //{
            //  if(phi >= (float(i_phi)-12.)*TMath::Pi()/12. && phi < (float(i_phi)-11.)*TMath::Pi()/12.)
            //  {
                // spectra 
                //TString KEY_Mass2 = Form("pt_%d_eta_%d_phi_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_eta,i_phi,cent9,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                //h_mMass2[KEY_Mass2]->Fill(InvMass,reweight);
                //if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
                //{
                //  TString KEY_Mass2Sys = Form("pt_%d_eta_%d_phi_%d_Centrality_9_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_eta,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                //  h_mMass2[KEY_Mass2Sys]->Fill(InvMass,reweight);
                //}
                TString KEY_Mass2 = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_eta,cent9,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                h_mMass2[KEY_Mass2]->Fill(InvMass,reweight);
                if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
                {
                  TString KEY_Mass2Sys = Form("pt_%d_eta_%d_Centrality_9_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_eta,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                  h_mMass2[KEY_Mass2Sys]->Fill(InvMass,reweight);
                }
                break;
            //  }
            //}
            //break;
          }
        }
        break;
      }	  
    }
  }

  //if(Res2 > 0.0)
  //{
  //  for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt_bin
  //  {
  //    if(pt >= vmsa::ptRawStart[i_pt] && pt < vmsa::ptRawStop[i_pt])
  //    {
  //      for(Int_t i_theta = 0; i_theta < 10; i_theta++) // phi-psi2 bin
  //      {
  //        if(CosThetaStar >= float(i_theta-5.)/5. && CosThetaStar < float(i_theta-4.)/5.)
  //        {
  //          for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi2 bin
  //          {
  //            if(phiprime >= TMath::Pi()*2.*float(i_phi)/10. && phiprime < TMath::Pi()*2.*float(i_phi+1.)/10.)
  //            {
  //              // spin alignment
  //              TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,cent9,i_theta,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
  //              h_mMass2[KEY_Mass2]->Fill(InvMass,reweight);
  //              if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
  //              {
  //                TString KEY_Mass2Sys = Form("pt_%d_Centrality_9_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_theta,i_phi,EP_string.c_str(),dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
  //                h_mMass2[KEY_Mass2Sys]->Fill(InvMass,reweight);
  //              }
  //              break;
  //            }
  //          }
  //          break;
  //        }
  //      }
  //      break;
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
  //for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin
  //{
  //  for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
  //  {
  //    for(Int_t i_theta = 0; i_theta < 10; i_theta++) // cos(theta*) bin
  //    {
  //      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // cos(theta*) bin
  //      {
  //        for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  //        {
  //          for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
  //          {
  //            if( i_dca != 0 && i_sig != 0) continue;
  //            TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_PhiPrime_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_cent,i_theta,i_phi,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
  //            h_mMass2[KEY_Mass2]->Write();
  //            cout << "Wrote: " << KEY_Mass2 << endl;     
  //          }
  //        }
  //      }
  //    }
  //  }
  //}
  for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin
  {
    for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_eta = 0; i_eta < 16; i_eta++) // cos(theta*) bin
      {
        //for(Int_t i_phi = 0; i_phi < 24; i_phi++) // cos(theta*) bin
        //{
          for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
          {
            for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
            {
              if( i_dca != 0 && i_sig != 0) continue;
              //TString KEY_Mass2 = Form("pt_%d_eta_%d_phi_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_eta,i_phi,i_cent,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
              //h_mMass2[KEY_Mass2]->Write();
              TString KEY_Mass2 = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_%s",i_pt,i_eta,i_cent,EP_string.c_str(),i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
              h_mMass2[KEY_Mass2]->Write();
              //cout << "Wrote: " << KEY_Mass2 << endl;     
            }
          }
        //}
      }
    }
  }
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
