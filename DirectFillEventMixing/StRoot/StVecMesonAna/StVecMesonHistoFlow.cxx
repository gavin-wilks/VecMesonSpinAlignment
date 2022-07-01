#include "StRoot/StVecMesonAna/StVecMesonHistoFlow.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"

ClassImp(StVecMesonHistoFlow)

//-------------------------------------------------------------
StVecMesonHistoFlow::StVecMesonHistoFlow()
{
}

StVecMesonHistoFlow::~StVecMesonHistoFlow()
{
}
//-------------------------------------------------------------

void StVecMesonHistoFlow::InitSys(Int_t X_flag, Int_t mode) // 0 for Same Event, 1 for Mixed Event
{
  // flow analysis
  for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin 
  {
    for(Int_t i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // centrality bin
    {
      for(Int_t i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // phi-psi bin
      {
	for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
	{
	  for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	  {
	    if(i_dca != 0 && i_sig != 0) continue;
            for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
            {
              if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
	      TString Mode[2] = {"SE","ME"};
	      TString KEY_Mass2 = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	      h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
	      h_mMass2[KEY_Mass2]->Sumw2();
            }
	  }
	}
      }
    }
  }
}
//-------------------------------------------------------------
void StVecMesonHistoFlow::FillSys(Float_t pt, Int_t cent9, Float_t PhiPsi, Int_t dcaSys, Int_t nSigSys, Int_t nHitSys, Float_t Res2, Float_t InvMass, Double_t reweight, Int_t X_flag, Int_t mode)
{
  TString Mode[2] = {"SE","ME"};
  if(Res2 > 0.0)
  {
    for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt_bin
    {
      if(pt >= vmsa::ptRawStart[i_pt] && pt < vmsa::ptRawStop[i_pt])
      {
	for(Int_t i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta++) // phi-psi2 bin
	{
	  if(PhiPsi >= vmsa::PhiPsi_low[i_theta] && PhiPsi < vmsa::PhiPsi_up[i_theta])
	  {
	    // spin alignment
	    TString KEY_Mass2 = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,cent9,i_theta,dcaSys,nSigSys,nHitSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	    h_mMass2[KEY_Mass2]->Fill(InvMass,(reweight/Res2));
          
	    TString KEY_Mass2Sys = Form("pt_%d_Centrality_9_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,i_theta,dcaSys,nSigSys,nHitSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	    h_mMass2[KEY_Mass2Sys]->Fill(InvMass,(reweight/Res2));
	    if(cent9 == 8 || cent9 == 7)
            {
	      TString KEY_Mass2Sys = Form("pt_%d_Centrality_10_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,i_theta,dcaSys,nSigSys,nHitSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	      h_mMass2[KEY_Mass2Sys]->Fill(InvMass,(reweight/Res2));
            }
	    if(cent9 <= 6 && cent9 >= 4 )
            {
	      TString KEY_Mass2Sys = Form("pt_%d_Centrality_11_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,i_theta,dcaSys,nSigSys,nHitSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	      h_mMass2[KEY_Mass2Sys]->Fill(InvMass,(reweight/Res2));
            }
	    if(cent9 <= 3 && cent9 >= 0)
            {
	      TString KEY_Mass2Sys = Form("pt_%d_Centrality_12_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,i_theta,dcaSys,nSigSys,nHitSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	      h_mMass2[KEY_Mass2Sys]->Fill(InvMass,(reweight/Res2));
            }     
	  }
	}
      }
    }
  }
}

//-------------------------------------------------------------
void StVecMesonHistoFlow::WriteSys(Int_t X_flag, Int_t mode)
{
  TString Mode[2] = {"SE","ME"};
  // flow
  for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin
  {
    for(Int_t i_cent = vmsa::Cent_start; i_cent < 13; i_cent++) // centrality bin
    {
      for(Int_t i_theta = 0; i_theta < vmsa::PhiPsi_total; i_theta ++) // cos(theta*) bin
      {
	for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
	{
	  for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
	  {
            if(i_dca != 0 && i_sig != 0) continue;
            for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
            {
              if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
	      TString KEY_Mass2 = Form("pt_%d_Centrality_%d_PhiPsi_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	      h_mMass2[KEY_Mass2]->Write();
            }
	  }
	}
      }
    }
  }
}
