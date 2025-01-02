#include "StRoot/StVecMesonAna/StVecMesonHistoManger.h"
#include "../Utility/StSpinAlignmentCons.h"
//#include "../Utility/phi_data_constants_19GeV.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
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
void StVecMesonHistoManger::InitPhiSys(Int_t X_flag, Int_t mode) // 0 for Same Event, 1 for Mixed Event
{
  TString Mode[2] = {"SE","ME"};
  for(Int_t i_pt_phi = 0; i_pt_phi < 12; i_pt_phi++)
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++)
    {
      TString KEY1 = Form("phi_phipt_%d_cos2phistarphi_%d_%s",i_pt_phi,i_phi,Mode[X_flag].Data());
      h_mPhi[KEY1] = new TH2F(KEY1.Data(),KEY1.Data(),data_constants::kaon_azimuthal_bins*6,-TMath::Pi(),TMath::Pi(),data_constants::kaon_rapidity_bins*3,-1,1);
      h_mPhi[KEY1]->Sumw2();
    }
    TString KEY_psp = Form("phi_phistarphi_phipt_%d_%s",i_pt_phi,Mode[X_flag].Data());
    TString KEY_psmp = Form("phi_phistarmphi_phipt_%d_%s",i_pt_phi,Mode[X_flag].Data());
    h_mPhi[KEY_psp] = new TH2F(KEY_psp.Data(),KEY_psp.Data(),data_constants::kaon_azimuthal_bins*6,-TMath::Pi(),TMath::Pi(),data_constants::kaon_azimuthal_bins*6,-TMath::Pi(),TMath::Pi());
    h_mPhi[KEY_psp]->Sumw2();
    h_mPhi1D[KEY_psmp] = new TH1F(KEY_psmp.Data(),KEY_psmp.Data(),data_constants::kaon_azimuthal_bins*6,-2*TMath::Pi(),2*TMath::Pi());
    h_mPhi1D[KEY_psmp]->Sumw2();
  }
}


void StVecMesonHistoManger::InitHistQA(Int_t X_flag, Int_t mode) // 0 for Same Event, 1 for Mixed Event
{
  TString Mode[2] = {"SE","ME"};
  string charge[2] = {"plus","minus"};

  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    string HistName = Form("h_mVertex_Cent_%d",i_cent);
    h_mVertex[i_cent] = new TH3F(HistName.c_str(),HistName.c_str(),50,-2.5,2.5,50,-2.5,2.5,280,-70.0,70.0);

    HistName = Form("h_mNToFMatch_Cent_%d",i_cent);
    h_mNToFMatch[i_cent] = new TH1F(HistName.c_str(),HistName.c_str(),501,-0.5,500.5);

    HistName = Form("h_mRefMult_Cent_%d",i_cent);
    h_mRefMult[i_cent] = new TH1F(HistName.c_str(),HistName.c_str(),501,-0.5,500.5);

    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      HistName = Form("h_mDca_K%s_Cent_%d_%s",charge[i_charge].c_str(),i_cent,Mode[X_flag].Data());
      h_mDcaTrack[i_charge][i_cent] = new TH1F(HistName.c_str(),HistName.c_str(),60,0.0,3.0);
      
      HistName = Form("h_mNHits_K%s_Cent_%d_%s",charge[i_charge].c_str(),i_cent,Mode[X_flag].Data());
      h_mNHits[i_charge][i_cent] = new TH1F(HistName.c_str(),HistName.c_str(),101,-0.5,100.5);

      HistName = Form("h_mNHitsRatio_K%s_Cent_%d_%s",charge[i_charge].c_str(),i_cent,Mode[X_flag].Data());
      h_mNHitsRatio[i_charge][i_cent] = new TH1F(HistName.c_str(),HistName.c_str(),100,0.0,1.0);

      HistName = Form("h_mDEdx_K%s_Cent_%d_%s",charge[i_charge].c_str(),i_cent,Mode[X_flag].Data());
      h_mDEdx_Kaon[i_charge][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),500,0.0,10.0,1000,0,40);
    }
  }
}

void StVecMesonHistoManger::FillEventHistQA(float reweight, int cent, float vx, float vy, float vz, float ntof, float refmult)
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

void StVecMesonHistoManger::FillTrackHistQA(float phi_pt, float phi_y, float InvMass, float reweight, int cent, int charge, float dca, float nhits, float nhitsratio, float p, float dedx, int X_flag)
{
  for(Int_t i_pt_phi = 0; i_pt_phi < data_constants::rebinpttotal; i_pt_phi++)
  {
    if(phi_pt >= data_constants::rebinptval[i_pt_phi] && phi_pt < data_constants::rebinptval[i_pt_phi+1])
    {
      //final_phi_pt_bin = data_constants::rebinptfinal[i_pt_phi];
      for(Int_t i_y_phi = 0; i_y_phi < data_constants::rebinytotal; i_y_phi++)
      {
        if(phi_y >= data_constants::rebinyval[i_y_phi] && phi_y < data_constants::rebinyval[i_y_phi+1])
        {
          //float scale = (1.0 - data_constants::fpoly[i_pt_phi][i_y_phi]);
          //if(X_flag == 1) scale *= data_constants::scalingME[i_pt_phi][i_y_phi];    
          //if(   InvMass < data_constants::InvMassLow[i_pt_phi][i_y_phi]   
          //   || InvMass > data_constants::InvMassHigh[i_pt_phi][i_y_phi]) return;
          //h_mDcaTrack[charge][cent]->Fill(dca,reweight*scale);
          //h_mNHits[charge][cent]->Fill(nhits,reweight*scale);
          //h_mNHitsRatio[charge][cent]->Fill(nhitsratio,reweight*scale);
          //h_mDEdx_Kaon[charge][cent]->Fill(p,dedx,reweight*scale);
          //if(cent >= 2 && cent <= 5)
          //{
          //  h_mDcaTrack[charge][9]->Fill(dca,reweight*scale);
          //  h_mNHits[charge][9]->Fill(nhits,reweight*scale);
          //  h_mNHitsRatio[charge][9]->Fill(nhitsratio,reweight*scale);
          //  h_mDEdx_Kaon[charge][9]->Fill(p,dedx,reweight*scale);
          //}
          //break;
        }
      }
      break;
    }
  }
} 

void StVecMesonHistoManger::WriteHistQA()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    h_mVertex[i_cent]->Write();
    h_mNToFMatch[i_cent]->Write();
    h_mRefMult[i_cent]->Write();
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      h_mDcaTrack[i_charge][i_cent]->Write();
      h_mNHits[i_charge][i_cent]->Write();
      h_mNHitsRatio[i_charge][i_cent]->Write();     
      h_mDEdx_Kaon[i_charge][i_cent]->Write();
    }
  }
}


void StVecMesonHistoManger::InitSys(Int_t X_flag, Int_t mode) // 0 for Same Event, 1 for Mixed Event
{
  TString Mode[2] = {"SE","ME"};
  for(Int_t i_pt_phi = 0; i_pt_phi < 4; i_pt_phi++)
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++)
    {
      for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
      {
        TString KEY = Form("kaonminus_phipt_%d_cos2phistarphi_%d_kaonpt_%d_%s",i_pt_phi,i_phi,i_pt_kaon,Mode[X_flag].Data());
        h_mKaon[KEY] = new TH2F(KEY.Data(),KEY.Data(),data_constants::kaon_azimuthal_bins*10,-TMath::Pi(),TMath::Pi(),data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
        h_mKaon[KEY]->Sumw2();
        cout << "INITIALIZED: " << KEY << endl;
        //KEY = Form("kaonminus_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta_%s",i_pt_phi,i_phi,i_pt_kaon,Mode[X_flag].Data());
        //h_mKaon[KEY] = new TH2F(KEY.Data(),KEY.Data(),data_constants::kaon_azimuthal_bins*10,-TMath::Pi(),TMath::Pi(),data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
        //h_mKaon[KEY]->Sumw2();
        //cout << "INITIALIZED: " << KEY << endl;
        KEY = Form("kaonplus_phipt_%d_cos2phistarphi_%d_kaonpt_%d_%s",i_pt_phi,i_phi,i_pt_kaon,Mode[X_flag].Data());
        h_mKaon[KEY] = new TH2F(KEY.Data(),KEY.Data(),data_constants::kaon_azimuthal_bins*10,-TMath::Pi(),TMath::Pi(),data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
        h_mKaon[KEY]->Sumw2();
        cout << "INITIALIZED: " << KEY << endl;
        //KEY = Form("kaonplus_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta_%s",i_pt_phi,i_phi,i_pt_kaon,Mode[X_flag].Data());
        //h_mKaon[KEY] = new TH2F(KEY.Data(),KEY.Data(),data_constants::kaon_azimuthal_bins*10,-TMath::Pi(),TMath::Pi(),data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
        //h_mKaon[KEY]->Sumw2();
        //cout << "INITIALIZED: " << KEY << endl;
        //h_mKaon[i_pt_phi][i_phi][i_pt_kaon] = new TH2F(KEY.Data(),KEY.Data(),data_constants::kaon_azimuthal_bins,-TMath::Pi(),TMath::Pi(),data_constants::kaon_rapidity_bins,0.0,data_constants::kaon_rapidity_max);
        //h_mKaon[i_pt_phi][i_phi][i_pt_kaon]->Sumw2();
      }
      //TString KEY1 = Form("kaon_delta_phipt_%d_cos2phistarphi_%d_%s",i_pt_phi,i_phi,Mode[X_flag].Data());
      //h_mKaon[KEY1] = new TH2F(KEY1.Data(),KEY1.Data(),data_constants::kaon_azimuthal_bins*6,-0.6,0.6,data_constants::kaon_rapidity_bins*3,-0.6,0.6);
      //h_mKaon[KEY1]->Sumw2();
      //KEY1 = Form("kaon_deltarphi_phipt_%d_cos2phistarphi_%d_%s",i_pt_phi,i_phi,Mode[X_flag].Data());
      //h_mKaon[KEY1] = new TH2F(KEY1.Data(),KEY1.Data(),720,0.0,2.0*TMath::Pi(),200,0,0.9);
      //h_mKaon[KEY1]->Sumw2();
      //KEY1 = Form("kaonplusphi_delta_phipt_%d_cos2phistarphi_%d_%s",i_pt_phi,i_phi,Mode[X_flag].Data());
      //h_mKaon[KEY1] = new TH2F(KEY1.Data(),KEY1.Data(),data_constants::kaon_azimuthal_bins*6,-0.5,0.5,data_constants::kaon_rapidity_bins*3,-0.4,0.4);
      //h_mKaon[KEY1]->Sumw2();
      //KEY1 = Form("kaonminusphi_delta_phipt_%d_cos2phistarphi_%d_%s",i_pt_phi,i_phi,Mode[X_flag].Data());
      //h_mKaon[KEY1] = new TH2F(KEY1.Data(),KEY1.Data(),data_constants::kaon_azimuthal_bins*6,-0.5,0.5,data_constants::kaon_rapidity_bins*3,-0.4,0.4);
      //h_mKaon[KEY1]->Sumw2();
    }
  }
}
//-------------------------------------------------------------
void StVecMesonHistoManger::FillPhiSys(Int_t cent9, Float_t pt, Float_t phi, Float_t y, Float_t phi_star, Float_t Cos2PhiStarPhi, Float_t InvMass, Double_t reweight, Int_t X_flag)
{
  TString Mode[2] = {"SE","ME"};
  if(cent9 < vmsa::cent_low[0] || cent9 > vmsa::cent_up[0]) return;// 20%-60%
  int final_phi_pt_bin = -1;
  for(Int_t i_pt_phi = 0; i_pt_phi < 12; i_pt_phi++)
  {
    if(pt >= data_constants::phi_pt_low[i_pt_phi] && pt < data_constants::phi_pt_high[i_pt_phi])
    {
      final_phi_pt_bin = data_constants::phi_final_pt_bins[i_pt_phi];
      //if(   InvMass < (data_constants::phi_center[final_phi_pt_bin] - 2.0*data_constants::phi_width[final_phi_pt_bin]) 
      //   || InvMass > (data_constants::phi_center[final_phi_pt_bin] + 2.0*data_constants::phi_width[final_phi_pt_bin])) return;
      //for(Int_t i_phi = 0; i_phi < 10; i_phi++)
      //{
      //  if(Cos2PhiStarPhi >= (float(i_phi)-5.)/5. && Cos2PhiStarPhi < (float(i_phi)-4.)/5.)
      //  {
      //    TString KEY = Form("phi_phipt_%d_cos2phistarphi_%d_%s",i_pt_phi,i_phi,Mode[X_flag].Data());
      //    TString KEY_psmp = Form("phi_phistarmphi_phipt_%d_%s",i_pt_phi,Mode[X_flag].Data());
      //    TString KEY_psp = Form("phi_phistarphi_phipt_%d_%s",i_pt_phi,Mode[X_flag].Data());
      //    if(X_flag == 0) {
      //      h_mPhi[KEY]->Fill(phi,y,reweight);
      //      h_mPhi[KEY_psp]->Fill(phi,phi_star,reweight);
      //      h_mPhi1D[KEY_psmp]->Fill(phi_star-phi,reweight);
      //    }
      //    if(X_flag == 1) {
      //      h_mPhi[KEY]->Fill(phi,y,reweight*data_constants::phi_norm[i_pt_phi][i_phi]);
      //      h_mPhi[KEY_psp]->Fill(phi,phi_star,reweight*data_constants::phi_norm[i_pt_phi][i_phi]);
      //      h_mPhi1D[KEY_psmp]->Fill(phi_star-phi,reweight*data_constants::phi_norm[i_pt_phi][i_phi]);
      //    }
      //    break;
      //  }
      //}
      //break;
    }
  }
}

void StVecMesonHistoManger::FillSys(Int_t cent9, Float_t phi_pt, Float_t phi_y, Float_t Cos2PhiStarPhi, Float_t kaon_pt, Float_t kaon_phi, Float_t kaon_y, Float_t kaon_eta, Int_t kcharge, Float_t InvMass, Double_t reweight, Int_t X_flag)
{
  TString Mode[2] = {"SE","ME"}; 
  string charge[2] = {"plus","minus"};
  if(cent9 < vmsa::cent_low[0] || cent9 > vmsa::cent_up[0]) return;// 20%-60%
  int final_phi_pt_bin = -1;
  for(Int_t i_pt_phi = 0; i_pt_phi < data_constants::rebinpttotal; i_pt_phi++)
  {
    if(phi_pt >= data_constants::rebinptval[i_pt_phi] && phi_pt < data_constants::rebinptval[i_pt_phi+1])
    {
      final_phi_pt_bin = data_constants::rebinptfinal[i_pt_phi];
      for(Int_t i_y_phi = 0; i_y_phi < data_constants::rebinytotal; i_y_phi++)
      {
        if(phi_y >= data_constants::rebinyval[i_y_phi] && phi_y < data_constants::rebinyval[i_y_phi+1])
        {
          
          //if(   InvMass < data_constants::InvMassLow[i_pt_phi][i_y_phi]   
          //   || InvMass > data_constants::InvMassHigh[i_pt_phi][i_y_phi]) return;
          //for(Int_t i_phi = 0; i_phi < 10; i_phi++)
          //{
          //  if(Cos2PhiStarPhi >= (float(i_phi)-5.)/5. && Cos2PhiStarPhi < (float(i_phi)-4.)/5.)
          //  {
          //    for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
          //    {
          //      if(kaon_pt >= data_constants::kaon_pt_low[i_pt_kaon] && kaon_pt < data_constants::kaon_pt_high[i_pt_kaon])   
          //      {
          //        TString KEY = Form("kaon%s_phipt_%d_cos2phistarphi_%d_kaonpt_%d_%s",charge[kcharge].c_str(),final_phi_pt_bin,i_phi,i_pt_kaon,Mode[X_flag].Data());
          //        if(X_flag == 0) {
          //          h_mKaon[KEY]->Fill(kaon_phi,kaon_y,reweight*(1.0 - data_constants::fpoly[i_pt_phi][i_y_phi]));
          //        }
          //        if(X_flag == 1) {
          //          h_mKaon[KEY]->Fill(kaon_phi,kaon_y,reweight*data_constants::scalingME[i_pt_phi][i_y_phi]*(1.0 - data_constants::fpoly[i_pt_phi][i_y_phi]));
          //        }
          //        break;
          //      }
          //    }
          //    break;
          //  }
          //}
          //break;
        }
      }
      break;
    }
  }
}

void StVecMesonHistoManger::FillDeltaSys(Int_t cent9, Float_t phi_pt, Float_t Cos2PhiStarPhi, Float_t phi, Float_t eta, Float_t kaonp_pt, Float_t kaonp_phi, Float_t kaonp_eta, Float_t kaonm_pt, Float_t kaonm_phi, Float_t kaonm_eta, Float_t InvMass, Double_t reweight, Int_t X_flag)
{
  TString Mode[2] = {"SE","ME"};
  
  if(cent9 < vmsa::cent_low[0] || cent9 > vmsa::cent_up[0]) return;// 20%-60%
  //while(kaonm_phi < 0.0) kaonm_phi += 2*TMath::Pi();
  //while(kaonp_phi < 0.0) kaonp_phi += 2*TMath::Pi();
  float phiphidiffp = phi-kaonp_phi;
  float phietadiffp = eta-kaonp_eta;
  float phiphidiffm = phi-kaonm_phi;
  float phietadiffm = eta-kaonm_eta;
  float phidiff = kaonp_phi-kaonm_phi;
  float etadiff = kaonp_eta-kaonm_eta;


  while(phiphidiffp < -TMath::Pi()) phiphidiffp += 2*TMath::Pi();
  while(phiphidiffp > TMath::Pi() ) phiphidiffp -= 2*TMath::Pi();
  while(phiphidiffm < -TMath::Pi()) phiphidiffm += 2*TMath::Pi();
  while(phiphidiffm > TMath::Pi() ) phiphidiffm -= 2*TMath::Pi();
  while(phidiff < -TMath::Pi()) phidiff += 2*TMath::Pi();
  while(phidiff > TMath::Pi() ) phidiff -= 2*TMath::Pi();

  if(phidiff < -TMath::Pi() || phidiff > TMath::Pi()) 
  {
    cout << "etadiff " << etadiff << ", phidiff = " << phidiff << " WARNING!" << endl;
    cout << "K+ phi = " << kaonp_phi << endl;
    cout << "K- phi = " << kaonm_phi << endl << endl;
  } 

  float deltar = TMath::Sqrt(etadiff*etadiff+phidiff*phidiff);
  float phi_r = TMath::ATan2(etadiff,phidiff);
  if(phi_r < 0.0) phi_r += 2.0*TMath::Pi();
  //if(phidiff < 0.0) phi_r += TMath::Pi(); 
  //if(phi_r < 0.0) phi_r += 2.0*TMath::Pi();

  int final_phi_pt_bin = -1;
  for(Int_t i_pt_phi = 0; i_pt_phi < 12; i_pt_phi++)
  {
    if(phi_pt >= data_constants::phi_pt_low[i_pt_phi] && phi_pt < data_constants::phi_pt_high[i_pt_phi])
    {
      final_phi_pt_bin = data_constants::phi_final_pt_bins[i_pt_phi];
      if(   InvMass < (data_constants::phi_center[final_phi_pt_bin] - 2.0*data_constants::phi_width[final_phi_pt_bin]) 
         || InvMass > (data_constants::phi_center[final_phi_pt_bin] + 2.0*data_constants::phi_width[final_phi_pt_bin])) return;
      for(Int_t i_phi = 0; i_phi < 10; i_phi++)
      {
        if(Cos2PhiStarPhi >= (float(i_phi)-5.)/5. && Cos2PhiStarPhi < (float(i_phi)-4.)/5.)
        {
          //cout << "phi_pt = " << phi_pt << ", i_pt_phi = " << i_pt_phi << ", final_phi_pt_bin = " << final_phi_pt_bin << ", Cos2PhiStarPhi = " << Cos2PhiStarPhi << ", i_phi = " << i_phi << ", kaon_pt = " << kaon_pt << ", i_pt_kaon = " << i_pt_kaon << endl;
          TString KEY = Form("kaon_delta_phipt_%d_cos2phistarphi_%d_%s",final_phi_pt_bin,i_phi,Mode[X_flag].Data());
          if(X_flag == 0) h_mKaon[KEY]->Fill(phidiff,etadiff,reweight);
          if(X_flag == 1) h_mKaon[KEY]->Fill(phidiff,etadiff,reweight*data_constants::phi_norm[i_pt_phi][i_phi]);
          KEY = Form("kaon_deltarphi_phipt_%d_cos2phistarphi_%d_%s",final_phi_pt_bin,i_phi,Mode[X_flag].Data());
          if(X_flag == 0) h_mKaon[KEY]->Fill(phi_r,deltar,reweight);
          if(X_flag == 1) h_mKaon[KEY]->Fill(phi_r,deltar,reweight*data_constants::phi_norm[i_pt_phi][i_phi]);
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
          KEY = Form("kaonplusphi_delta_phipt_%d_cos2phistarphi_%d_%s",final_phi_pt_bin,i_phi,Mode[X_flag].Data());
          if(X_flag == 0) h_mKaon[KEY]->Fill(phiphidiffp,phietadiffp,reweight);
          if(X_flag == 1) h_mKaon[KEY]->Fill(phiphidiffp,phietadiffp,reweight*data_constants::phi_norm[i_pt_phi][i_phi]);
          KEY = Form("kaonminusphi_delta_phipt_%d_cos2phistarphi_%d_%s",final_phi_pt_bin,i_phi,Mode[X_flag].Data());
          if(X_flag == 0) h_mKaon[KEY]->Fill(phiphidiffm,phietadiffm,reweight);
          if(X_flag == 1) h_mKaon[KEY]->Fill(phiphidiffm,phietadiffm,reweight*data_constants::phi_norm[i_pt_phi][i_phi]);
          //if(X_flag == 0) h_mKaon[final_phi_pt_bin][i_phi][i_pt_kaon]->Fill(kaon_phi,TMath::Abs(kaon_y),reweight);
          //if(X_flag == 1) h_mKaon[final_phi_pt_bin][i_phi][i_pt_kaon]->Fill(kaon_phi,TMath::Abs(kaon_y),reweight*data_constants::phi_norm[i_pt_phi][i_phi]);
          break;
        }
      }
      break;
    }
  }
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
  for(Int_t i_pt_phi = 0; i_pt_phi < 4; i_pt_phi++)
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++)
    {
      for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
      {
        TString KEY = Form("kaonplus_phipt_%d_cos2phistarphi_%d_kaonpt_%d_%s",i_pt_phi,i_phi,i_pt_kaon,Mode[X_flag].Data());
        //h_mKaon[i_pt_phi][i_phi][i_pt_kaon]->Write(); 
        h_mKaon[KEY]->Write(); 
        //TString KEY_eta = Form("kaonplus_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta_%s",i_pt_phi,i_phi,i_pt_kaon,Mode[X_flag].Data());
        //h_mKaon[i_pt_phi][i_phi][i_pt_kaon]->Write(); 
        //h_mKaon[KEY_eta]->Write(); 
        KEY = Form("kaonminus_phipt_%d_cos2phistarphi_%d_kaonpt_%d_%s",i_pt_phi,i_phi,i_pt_kaon,Mode[X_flag].Data());
        //h_mKaon[i_pt_phi][i_phi][i_pt_kaon]->Write(); 
        h_mKaon[KEY]->Write(); 
        //KEY_eta = Form("kaonminus_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta_%s",i_pt_phi,i_phi,i_pt_kaon,Mode[X_flag].Data());
        //h_mKaon[i_pt_phi][i_phi][i_pt_kaon]->Write(); 
        //h_mKaon[KEY_eta]->Write(); 
      } 
      //TString KEY1 = Form("kaon_delta_phipt_%d_cos2phistarphi_%d_%s",i_pt_phi,i_phi,Mode[X_flag].Data());
      //h_mKaon[KEY1]->Write(); 
      //KEY1 = Form("kaon_deltarphi_phipt_%d_cos2phistarphi_%d_%s",i_pt_phi,i_phi,Mode[X_flag].Data());
      //h_mKaon[KEY1]->Write(); 
      //KEY1 = Form("kaonplusphi_delta_phipt_%d_cos2phistarphi_%d_%s",i_pt_phi,i_phi,Mode[X_flag].Data());
      //h_mKaon[KEY1]->Write(); 
      //KEY1 = Form("kaonminusphi_delta_phipt_%d_cos2phistarphi_%d_%s",i_pt_phi,i_phi,Mode[X_flag].Data());
      //h_mKaon[KEY1]->Write(); 
    }
  }
}

void StVecMesonHistoManger::WritePhiSys(Int_t X_flag, Int_t mode)
{
  TString Mode[2] = {"SE","ME"};
  cout << "About to Write Histograms" << endl;
  for(Int_t i_pt_phi = 0; i_pt_phi < 12; i_pt_phi++)
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++)
    {
      TString KEY1 = Form("phi_phipt_%d_cos2phistarphi_%d_%s",i_pt_phi,i_phi,Mode[X_flag].Data());
      h_mPhi[KEY1]->Write(); 
    }
    TString KEY1 = Form("phi_phistarphi_phipt_%d_%s",i_pt_phi,Mode[X_flag].Data());
    h_mPhi[KEY1]->Write(); 
    KEY1 = Form("phi_phistarmphi_phipt_%d_%s",i_pt_phi,Mode[X_flag].Data());
    h_mPhi1D[KEY1]->Write(); 
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
