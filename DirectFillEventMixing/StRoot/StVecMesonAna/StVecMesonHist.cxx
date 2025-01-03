#include "StRoot/StVecMesonAna/StVecMesonHist.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"

ClassImp(StVecMesonHist)

//-------------------------------------------------------------
StVecMesonHist::StVecMesonHist()
{
}

StVecMesonHist::~StVecMesonHist()
{
}
//-------------------------------------------------------------

void StVecMesonHist::InitSys(Int_t X_flag, Int_t mode) // 0 for Same Event, 1 for Mixed Event
{
  // spin alignment analysis
  if(mode == 0) // phi
  {
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
              if(i_dca != 0 && i_sig != 0) continue;
              TString Mode[2] = {"SE","ME"};
              TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_%s",i_pt,i_cent,i_theta,i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
              h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
              h_mMass2[KEY_Mass2]->Sumw2();
            }
          }
        }
      }
    }
  }
  if(mode == 2) // K*0
  {
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
              if(i_dca != 0 && i_sig != 0) continue;
              for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
              {
                if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
                TString Mode[2] = {"SE","ME"};
                TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
                h_mMass2[KEY_Mass2]->Sumw2();
              }
            }
          }
        }
      }
    }
  }

  // raw pt spectra
  if(mode == 0)
  { 
    for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin
    {
      for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
      {
        for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
        {
          for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
          {
            if(i_dca != 0 && i_sig != 0) continue;
            for(Int_t i_pT = 0; i_pT < 2; i_pT++)
            {
              TString Mode[2] = {"SE","ME"};
              TString pT[2] = {"low","high"};
              TString KEY_Mass2_Spec = Form("Spec_pt_%d_%s_Centrality_%d_Dca_%d_Sig_%d_%s_%s",i_pt,pT[i_pT].Data(),i_cent,i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
              h_mMass_Spec[KEY_Mass2_Spec] = new TH1F(KEY_Mass2_Spec.Data(),KEY_Mass2_Spec.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
              h_mMass_Spec[KEY_Mass2_Spec]->Sumw2();
            } 
          }
        }
      }
    }
  }

  if(mode == 2)
  { 
    for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin
    {
      for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
      {
        for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
        {
          for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
          {
            if(i_dca != 0 && i_sig != 0) continue;
            for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
            {
              if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
              for(Int_t i_pT = 0; i_pT < 2; i_pT++)
              {
                TString Mode[2] = {"SE","ME"};
                TString pT[2] = {"low","high"};
                TString KEY_Mass2_Spec = Form("Spec_pt_%d_%s_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,pT[i_pT].Data(),i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                h_mMass_Spec[KEY_Mass2_Spec] = new TH1F(KEY_Mass2_Spec.Data(),KEY_Mass2_Spec.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
                h_mMass_Spec[KEY_Mass2_Spec]->Sumw2();
              }
            }
          }
        }
      }
    }
  }

  
  // Yields
  if(mode == 0)
  {
    for(Int_t i_cent = 0; i_cent < 9; i_cent++) // centrality bin
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
      {
        for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
        {
          if(i_dca != 0 && i_sig != 0) continue;
          TString Mode[2] = {"SE","ME"};
          TString KEY_Mass2_Yields = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_%s_%s",i_cent,i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
          h_mMass_Yields[KEY_Mass2_Yields] = new TH1F(KEY_Mass2_Yields.Data(),KEY_Mass2_Yields.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
          h_mMass_Yields[KEY_Mass2_Yields]->Sumw2();
        }
      }
    }
  }

  if(mode == 2)
  {
    for(Int_t i_cent = 0; i_cent < 9; i_cent++) // centrality bin
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
            TString KEY_Mass2_Yields = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
            h_mMass_Yields[KEY_Mass2_Yields] = new TH1F(KEY_Mass2_Yields.Data(),KEY_Mass2_Yields.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
            h_mMass_Yields[KEY_Mass2_Yields]->Sumw2();
          }
        }
      }
    }
  }

  // Dca QA
  for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  {
    for(Int_t i_charge = vmsa::Charge_start; i_charge < vmsa::Charge_stop; i_charge++)
    {
      TString Charge[2] = {"A","B"};
      TString KEY_Dca = Form("Tracks_Dca%s_%d",Charge[i_charge].Data(),i_dca);
      h_mDca[KEY_Dca] = new TH1F(KEY_Dca.Data(),KEY_Dca.Data(),100,-5.0,5.0);
      h_mDca[KEY_Dca]->Sumw2();
    }
  }

  // nSigKaon QA
  for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
  {
    for(Int_t i_charge = vmsa::Charge_start; i_charge < vmsa::Charge_stop; i_charge++)
    {
      TString Charge[2] = {"A","B"};
      TString KEY_nSigKaon = Form("Tracks_Sig%s_%d",Charge[i_charge].Data(),i_sig);
      h_mSigKaon[KEY_nSigKaon] = new TH1F(KEY_nSigKaon.Data(),KEY_nSigKaon.Data(),100,-5.0,5.0);
      h_mSigKaon[KEY_nSigKaon]->Sumw2();
    }
  }
  
  // nHit QA
  if(mode == 2)
  {
    for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
    {
      for(Int_t i_charge = vmsa::Charge_start; i_charge < vmsa::Charge_stop; i_charge++)
      {
        TString Charge[2] = {"A","B"};
        TString KEY_NHit = Form("Tracks_nHit%s_%d",Charge[i_charge].Data(),i_nhit);
        h_mNHit[KEY_NHit] = new TH1F(KEY_NHit.Data(),KEY_NHit.Data(),100,10,80);
        h_mNHit[KEY_NHit]->Sumw2();
      }
    }
  }
}
//-------------------------------------------------------------
void StVecMesonHist::FillSys(Float_t pt, Int_t cent9, Float_t CosThetaStar, Int_t dcaSys, Int_t nSigSys, Int_t nHitSys, Float_t Res2, Float_t InvMass, Double_t reweight, Int_t X_flag, Int_t mode)
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
	    if(mode == 0)
            {
	      TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_%s",i_pt,cent9,i_theta,dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	      h_mMass2[KEY_Mass2]->Fill(InvMass,reweight);
	      if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
	      {
	        TString KEY_Mass2Sys = Form("pt_%d_Centrality_9_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_%s",i_pt,i_theta,dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	        h_mMass2[KEY_Mass2Sys]->Fill(InvMass,reweight);
	      }
	      // raw pt spectra
	      if(pt < 0.5*(vmsa::ptRawStart[i_pt]+vmsa::ptRawStop[i_pt])) 
	      {
	        TString KEY_Mass2_Spec = Form("Spec_pt_%d_low_Centrality_%d_Dca_%d_Sig_%d_%s_%s",i_pt,cent9,dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	        h_mMass_Spec[KEY_Mass2_Spec]->Fill(InvMass,reweight);
	        if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
	        {
	          TString KEY_Mass2_SpecSys = Form("Spec_pt_%d_low_Centrality_9_Dca_%d_Sig_%d_%s_%s",i_pt,dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	          h_mMass_Spec[KEY_Mass2_SpecSys]->Fill(InvMass,reweight);
	        }
	      }
	      else
	      {
	        TString KEY_Mass2_Spec = Form("Spec_pt_%d_high_Centrality_%d_Dca_%d_Sig_%d_%s_%s",i_pt,cent9,dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	        h_mMass_Spec[KEY_Mass2_Spec]->Fill(InvMass,reweight);
	        if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
	        {
	          TString KEY_Mass2_SpecSys = Form("Spec_pt_%d_high_Centrality_9_Dca_%d_Sig_%d_%s_%s",i_pt,dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	          h_mMass_Spec[KEY_Mass2_SpecSys]->Fill(InvMass,reweight);
	        }
	      }
            }
	    if(mode == 2)
            {
	      TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,cent9,i_theta,dcaSys,nSigSys,nHitSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	      h_mMass2[KEY_Mass2]->Fill(InvMass,reweight);
	      if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
	      {
	        TString KEY_Mass2Sys = Form("pt_%d_Centrality_9_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,i_theta,dcaSys,nSigSys,nHitSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	        h_mMass2[KEY_Mass2Sys]->Fill(InvMass,reweight);
	      }
	      // raw pt spectra
	      if(pt < 0.5*(vmsa::ptRawStart[i_pt]+vmsa::ptRawStop[i_pt])) 
	      {
	        TString KEY_Mass2_Spec = Form("Spec_pt_%d_low_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,cent9,dcaSys,nSigSys,nHitSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	        h_mMass_Spec[KEY_Mass2_Spec]->Fill(InvMass,reweight);
	        if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
	        {
	          TString KEY_Mass2_SpecSys = Form("Spec_pt_%d_low_Centrality_9_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,dcaSys,nSigSys,nHitSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	          h_mMass_Spec[KEY_Mass2_SpecSys]->Fill(InvMass,reweight);
	        }
	      }
	      else
	      {
	        TString KEY_Mass2_Spec = Form("Spec_pt_%d_high_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,cent9,dcaSys,nSigSys,nHitSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	        h_mMass_Spec[KEY_Mass2_Spec]->Fill(InvMass,reweight);
	        if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
	        {
	          TString KEY_Mass2_SpecSys = Form("Spec_pt_%d_high_Centrality_9_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,dcaSys,nSigSys,nHitSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	          h_mMass_Spec[KEY_Mass2_SpecSys]->Fill(InvMass,reweight);
	        }
	      }
            }
	  }
	}
      }
    }
  }
  if(mode == 0)
  {
    TString KEY_Mass2_Yields = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_%s_%s",cent9,dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
    h_mMass_Yields[KEY_Mass2_Yields]->Fill(InvMass,reweight);
  }
  if(mode == 2)
  {
    TString KEY_Mass2_Yields = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_%s",cent9,dcaSys,nSigSys,nHitSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
    h_mMass_Yields[KEY_Mass2_Yields]->Fill(InvMass,reweight);
  }
}

void StVecMesonHist::FillDcaSys(Float_t dcaA, Float_t dcaB, Int_t dcaSys)
{
  TString KEY_DcaA = Form("Tracks_DcaA_%d",dcaSys);
  h_mDca[KEY_DcaA]->Fill(dcaA);
  TString KEY_DcaB = Form("Tracks_DcaB_%d",dcaSys);
  h_mDca[KEY_DcaB]->Fill(dcaB);
}

void StVecMesonHist::FillSigSys(Float_t nsA, Float_t nsB, Int_t nSigSys)
{
  TString KEY_nSigKaonA = Form("Tracks_SigA_%d",nSigSys);
  h_mSigKaon[KEY_nSigKaonA]->Fill(nsA); 
  TString KEY_nSigKaonB = Form("Tracks_SigB_%d",nSigSys);
  h_mSigKaon[KEY_nSigKaonB]->Fill(nsB); 
}

void StVecMesonHist::FillNHitSys(Int_t nsA, Int_t nsB, Int_t nNHitSys)
{
  TString KEY_nNHitA = Form("Tracks_nHitA_%d",nNHitSys);
  h_mNHit[KEY_nNHitA]->Fill(nsA); 
  TString KEY_nNHitB = Form("Tracks_nHitB_%d",nNHitSys);
  h_mNHit[KEY_nNHitB]->Fill(nsB); 
}

//-------------------------------------------------------------
void StVecMesonHist::WriteSys(Int_t X_flag, Int_t mode)
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
            if(i_dca != 0 && i_sig != 0) continue;
	    if(mode == 0)
            {
              TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_%s",i_pt,i_cent,i_theta,i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	      h_mMass2[KEY_Mass2]->Write();
            }
            if(mode == 2)
            {
              for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
              {
                if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
	        TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,i_cent,i_theta,i_dca,i_sig,i_nhit,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	        h_mMass2[KEY_Mass2]->Write();
              }
            }
	  }
	}
      }
    }
  }

  // Yields
  for(Int_t i_cent = 0; i_cent < 9; i_cent++) // centrality bin
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if(i_dca != 0 && i_sig != 0) continue;
        if(mode == 0)
        {
	  TString KEY_Mass2_Yields = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_%s_%s",i_cent,i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	  h_mMass_Yields[KEY_Mass2_Yields]->Write();
        }
        if(mode == 2)
        {
          for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
          {
            if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
	    TString KEY_Mass2_Yields = Form("Yields_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	    h_mMass_Yields[KEY_Mass2_Yields]->Write();
          }
        }
      }
    }
  }

  // raw pt spectra
  for(Int_t i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin
  {
    for(Int_t i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
      {
        for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
        {
          if(i_dca != 0 && i_sig != 0) continue;
          if(mode == 0)
          {
            for(Int_t i_pT = 0; i_pT < 2; i_pT++)
            {
              TString pT[2] = {"low","high"};
              TString KEY_Mass2_Spec = Form("Spec_pt_%d_%s_Centrality_%d_Dca_%d_Sig_%d_%s_%s",i_pt,pT[i_pT].Data(),i_cent,i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
              h_mMass_Spec[KEY_Mass2_Spec]->Write();
            }
          }
          if(mode == 2)
          {
            for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
            {
              if((i_sig != 0 && i_nhit != 0) || (i_dca != 0 && i_nhit != 0)) continue;
              for(Int_t i_pT = 0; i_pT < 2; i_pT++)
              {
                TString pT[2] = {"low","high"};
                TString KEY_Mass2_Spec = Form("Spec_pt_%d_%s_Centrality_%d_Dca_%d_Sig_%d_NHit_%d_%s_%s",i_pt,pT[i_pT].Data(),i_cent,i_dca,i_sig,i_nhit,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
                h_mMass_Spec[KEY_Mass2_Spec]->Write();
              }
            }
          }
        }
      }
    }
  }

  TString Charge[2] = {"A","B"};
  // Dca QA
  for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  {
    for(Int_t i_charge = vmsa::Charge_start; i_charge < vmsa::Charge_stop; i_charge++)
    {
      TString KEY_Dca = Form("Tracks_Dca%s_%d",Charge[i_charge].Data(),i_dca);
      h_mDca[KEY_Dca]->Write();
    }
  }

  // nSigKaon QA
  for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
  {
    for(Int_t i_charge = vmsa::Charge_start; i_charge < vmsa::Charge_stop; i_charge++)
    {
      TString KEY_nSigKaon = Form("Tracks_Sig%s_%d",Charge[i_charge].Data(),i_sig);
      h_mSigKaon[KEY_nSigKaon]->Write();
    }
  }

  if(mode == 2)
  {
    for(Int_t i_nhit = vmsa::mNHit_start; i_nhit < vmsa::mNHit_stop; i_nhit++)
    {
      for(Int_t i_charge = vmsa::Charge_start; i_charge < vmsa::Charge_stop; i_charge++)
      {
        TString KEY_NHit = Form("Tracks_nHit%s_%d",Charge[i_charge].Data(),i_nhit);
        h_mNHit[KEY_NHit]->Write();
      }
    }
  }
}
//-----------------------EP Direction----------------------------

void StVecMesonHist::InitSys_EP(Int_t X_flag, Int_t mode) // 0 for Same Event, 1 for Mixed Event
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
	    TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_%s_EP",i_pt,i_cent,i_theta,i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	    h_mMass2_EP[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
	    h_mMass2_EP[KEY_Mass2]->Sumw2();
	  }
	}
      }
    }
  }
}
//-------------------------------------------------------------
void StVecMesonHist::FillSys_EP(Float_t pt, Int_t cent9, Float_t CosThetaStar, Int_t dcaSys, Int_t nSigSys, Float_t Res2, Float_t InvMass, Double_t reweight, Int_t X_flag, Int_t mode)
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
	    TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_%s_EP",i_pt,cent9,i_theta,dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	    h_mMass2_EP[KEY_Mass2]->Fill(InvMass,reweight);
	    if(cent9 >= vmsa::cent_low[0] && cent9 <= vmsa::cent_up[0]) // 20%-60%
	    {
	      TString KEY_Mass2Sys = Form("pt_%d_Centrality_9_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_%s_EP",i_pt,i_theta,dcaSys,nSigSys,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	      h_mMass2_EP[KEY_Mass2Sys]->Fill(InvMass,reweight);
	    }
	  }
	}
      }
    }
  }
}

//-------------------------------------------------------------
void StVecMesonHist::WriteSys_EP(Int_t X_flag, Int_t mode)
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
	    TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_2nd_Dca_%d_Sig_%d_%s_%s_EP",i_pt,i_cent,i_theta,i_dca,i_sig,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	    h_mMass2_EP[KEY_Mass2]->Write();
	  }
	}
      }
    }
  }
}
