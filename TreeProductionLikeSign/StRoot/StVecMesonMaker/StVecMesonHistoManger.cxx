#include "StRoot/StVecMesonMaker/StVecMesonHistoManger.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TString.h"

ClassImp(StVecMesonHistoManger)

//-------------------------------------------------------------------------------------------

StVecMesonHistoManger::StVecMesonHistoManger()
{
}

//-------------------------------------------------------------------------------------------

StVecMesonHistoManger::~StVecMesonHistoManger()
{
  /* */
}

//-------------------------------------------------------------------------------------------
void StVecMesonHistoManger::InitQA()
{
  h_mDEdx = new TH2F("h_mDEdx","h_mDEdx",1000,0,4.0,1000,0,40);
  h_mMass2 = new TH2F("h_mMass2","h_mMass2",1000,0,4.0,1000,-0.3,1.7);
  h_mFullRaw = new TH1F("h_mFullRaw","h_mFullRaw",360,-TMath::Pi(),TMath::Pi());
  h_mEastRaw = new TH1F("h_mEastRaw","h_mEastRaw",360,-TMath::Pi(),TMath::Pi());
  h_mWestRaw = new TH1F("h_mWestRaw","h_mWestRaw",360,-TMath::Pi(),TMath::Pi());
  h_mVz   = new TH1F("h_mVz","h_mVz",201,-100.5,100.5);
  h_mRefMult = new TH1F("h_mRefMult","h_mRefMult",1000,-0.5,999.5);
}

void StVecMesonHistoManger::InitPID()
{
  //for(int icent = 0; icent < 10; icent++)
  //{
    //for(int ic = 0; ic < 2; ic++)
    //{
    //  for(int ipt = 0; ipt < 20; ipt++)
    //  {
    //    for(int ieta = 0; ieta < 30; ieta++)
    //    {
    //      for(int iphi = 0; iphi < 1; iphi++)
    //      {
    //        string HistName = Form("h_mPID_cent%d_charge%d_pt%d_eta%d_phi%d",9,ic,ipt,ieta,iphi);
    //        h_mPID[HistName] = new TH2F(HistName.c_str(),HistName.c_str(),150,-7.5,7.5,300,-0.5,1.5);
    //        //h_mPID[ic][ipt][ieta][iphi] = new TH2F(HistName.c_str(),HistName.c_str(),150,-7.5,7.5,300,-0.5,1.5);
    //        //h_mPID[icent][ic][ipt][ieta] = new TH2F(HistName.c_str(),HistName.c_str(),200,-10,10,500,-0.5,1.5);
    //      }
    //    }
    //  }
    //}
  //}
    for(int ic = 0; ic < 2; ic++)
    {
      for(int ipt = 0; ipt < 40; ipt++)
      {
        for(int ieta = 0; ieta < 24; ieta++)
        {
          for(int iphi = 0; iphi < 12; iphi++)
          {
            string HistName = Form("h_mPID_cent%d_charge%d_pt%d_eta%d_phi%d",9,ic,ipt,ieta,iphi);
            h_mPID_1D[HistName] = new TH1F(HistName.c_str(),HistName.c_str(),200,-10.0,10.0);
            HistName = Form("h_mPID_cent%d_charge%d_pt%d_eta%d_phi%d_m2",9,ic,ipt,ieta,iphi);
            h_mPID_1D[HistName] = new TH1F(HistName.c_str(),HistName.c_str(),200,-10.0,10.0);
            //h_mPID[HistName] = new TH2F(HistName.c_str(),HistName.c_str(),150,-7.5,7.5,300,-0.5,1.5);
            //h_mPID[ic][ipt][ieta][iphi] = new TH2F(HistName.c_str(),HistName.c_str(),150,-7.5,7.5,300,-0.5,1.5);
            //h_mPID[icent][ic][ipt][ieta] = new TH2F(HistName.c_str(),HistName.c_str(),200,-10,10,500,-0.5,1.5);
          }
        }
      }
    }
}

void StVecMesonHistoManger::FillPID(int cent, short charge, double pt, double eta, double phi, double nsig, double mass2)
{

  int chargeidx = -1; 
  if(charge < 0) chargeidx = 0;  
  else if(charge > 0) chargeidx = 1;  

  for(int ipt = 0; ipt < 40; ipt++)
  {
    if(pt >= 0.1 + double(ipt)*(4.9/40.) && pt < 0.1 + double(ipt+1)*(4.9/40.))
    {
      for(int ieta = 0; ieta < 24; ieta++)
      {
        if(eta >= (double(ieta)-12.)/12. && eta < (double(ieta)-11.)/12. )
        {
          for(int iphi = 0; iphi < 12; iphi++)
          {
            if(phi >= TMath::Pi()*(double(iphi)-6.)/6. && phi < TMath::Pi()*(double(iphi)-5.)/6. )
            {
              //string HistName = Form("h_mPID_cent%d_charge%d_pt%d_eta%d",cent,chargeidx,ipt,ieta);
              //h_mPID[HistName]->Fill(nsig,mass2);
              //h_mPID[cent][chargeidx][ipt][ieta]->Fill(nsig,mass2);
              if(cent >= 2 && cent <= 5)
              {
                //string HistName = Form("h_mPID_cent9_charge%d_pt%d_eta%d",chargeidx,ipt,ieta);
                string HistName = Form("h_mPID_cent%d_charge%d_pt%d_eta%d_phi%d",9,chargeidx,ipt,ieta,iphi);
                h_mPID_1D[HistName]->Fill(nsig);
                if(mass2 > -900.0)
                {
                  HistName = Form("h_mPID_cent%d_charge%d_pt%d_eta%d_phi%d_m2",9,chargeidx,ipt,ieta,iphi);
                  h_mPID_1D[HistName]->Fill(nsig);
                }
                //h_mPID[HistName]->Fill(nsig,mass2);
                //h_mPID[chargeidx][ipt][ieta][iphi]->Fill(nsig,mass2);
                //h_mPID[HistName]->Fill(nsig,mass2);
              }
              //cout << "Filled for cent = " << cent << endl;
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

void StVecMesonHistoManger::WritePID()
{
  //for(int icent = 0; icent < 10; icent++)
  //{
    for(int ic = 0; ic < 2; ic++)
    {
      for(int ipt = 0; ipt < 40; ipt++)
      {
        for(int ieta = 0; ieta < 24; ieta++)
        {
          for(int iphi = 0; iphi < 12; iphi++)
          {
            //string HistName = Form("h_mPID_cent%d_charge%d_pt%d_eta%d",icent,ic,ipt,ieta);
            string HistName = Form("h_mPID_cent%d_charge%d_pt%d_eta%d_phi%d",9,ic,ipt,ieta,iphi);
            //h_mPID[HistName]->Write();
            h_mPID_1D[HistName]->Write();
            HistName = Form("h_mPID_cent%d_charge%d_pt%d_eta%d_phi%d_m2",9,ic,ipt,ieta,iphi);
            //h_mPID[HistName]->Write();
            h_mPID_1D[HistName]->Write();
            //h_mPID[ic][ipt][ieta][iphi]->Write();
          }
        }
      }
    }
  //}
}

void StVecMesonHistoManger::FillQA_Detector(Float_t dEdx, Float_t Mass2, Float_t p)
{
  h_mDEdx->Fill(p,dEdx);
  h_mMass2->Fill(p,Mass2);
}

void StVecMesonHistoManger::FillQA_Event(Float_t vz, Float_t refMult)
{
  h_mVz->Fill(vz);
  h_mRefMult->Fill(refMult);
}

void StVecMesonHistoManger::FillEP_Eta(Float_t Psi2_East, Float_t Psi2_West)
{
  h_mEastRaw->Fill(Psi2_East);
  h_mWestRaw->Fill(Psi2_West);
}

void StVecMesonHistoManger::FillEP_Full(Float_t Psi2_Full)
{
  h_mFullRaw->Fill(Psi2_Full);
}

void StVecMesonHistoManger::WriteQA()
{
  h_mDEdx->Write();
  h_mMass2->Write();
  h_mVz->Write();
  h_mRefMult->Write();
  h_mFullRaw->Write();
  h_mEastRaw->Write();
  h_mWestRaw->Write();
}
//-------------------------------------------------------------------------------------------
void StVecMesonHistoManger::InitEP()
{
  h_mEastReCenter = new TH1F("h_mEastReCenter","h_mEastReCenter",360,-TMath::Pi(),TMath::Pi());
  h_mWestReCenter = new TH1F("h_mWestReCenter","h_mWestReCenter",360,-TMath::Pi(),TMath::Pi());
  h_mRanAReCenter = new TH1F("h_mRanAReCenter","h_mRanAReCenter",360,-TMath::Pi(),TMath::Pi());
  h_mRanBReCenter = new TH1F("h_mRanBReCenter","h_mRanBReCenter",360,-TMath::Pi(),TMath::Pi());
  h_mFullReCenter = new TH1F("h_mFullReCenter","h_mFullReCenter",360,-TMath::Pi(),TMath::Pi());

  h_mEastShift = new TH1F("h_mEastShift","h_mEastShift",360,-TMath::Pi(),TMath::Pi());
  h_mWestShift = new TH1F("h_mWestShift","h_mWestShift",360,-TMath::Pi(),TMath::Pi());
  h_mRanAShift = new TH1F("h_mRanAShift","h_mRanAShift",360,-TMath::Pi(),TMath::Pi());
  h_mRanBShift = new TH1F("h_mRanBShift","h_mRanBShift",360,-TMath::Pi(),TMath::Pi());
  h_mFullShift = new TH1F("h_mFullShift","h_mFullShift",360,-TMath::Pi(),TMath::Pi());
}

void StVecMesonHistoManger::FillEP_Sub(Float_t Psi2East_ReCenter, Float_t Psi2East_Shift, Float_t Psi2West_ReCenter, Float_t Psi2West_Shift)
{
  h_mEastReCenter->Fill(Psi2East_ReCenter);
  h_mEastShift->Fill(Psi2East_Shift);

  h_mWestReCenter->Fill(Psi2West_ReCenter);
  h_mWestShift->Fill(Psi2West_Shift);
}

void StVecMesonHistoManger::FillEP_Ran(Float_t Psi2RanA_ReCenter, Float_t Psi2RanA_Shift, Float_t Psi2RanB_ReCenter, Float_t Psi2RanB_Shift, Float_t Psi2Full_ReCenter, Float_t Psi2Full_Shift)
{
  h_mRanAReCenter->Fill(Psi2RanA_ReCenter);
  h_mRanAShift->Fill(Psi2RanA_Shift);

  h_mRanBReCenter->Fill(Psi2RanB_ReCenter);
  h_mRanBShift->Fill(Psi2RanB_Shift);

  h_mFullReCenter->Fill(Psi2Full_ReCenter);
  h_mFullShift->Fill(Psi2Full_Shift);
}

void StVecMesonHistoManger::WriteEP()
{
  h_mEastReCenter->Write();
  h_mWestReCenter->Write();
  h_mRanAReCenter->Write();
  h_mRanBReCenter->Write();
  h_mFullReCenter->Write();

  h_mEastShift->Write();
  h_mWestShift->Write();
  h_mRanAShift->Write();
  h_mRanBShift->Write();
  h_mFullShift->Write();
}
