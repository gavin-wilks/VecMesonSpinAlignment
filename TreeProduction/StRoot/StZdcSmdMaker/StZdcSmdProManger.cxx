#include "StRoot/StZdcSmdMaker/StZdcSmdProManger.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TMath.h"

ClassImp(StZdcSmdProManger)

string StZdcSmdProManger::mVStr[2] = {"pos","neg"};

//---------------------------------------------------------------------------------

StZdcSmdProManger::StZdcSmdProManger()
{
}

//---------------------------------------------------------------------------------

StZdcSmdProManger::~StZdcSmdProManger()
{
  /* */
}

//---------------------------------------------------------------------------------

void StZdcSmdProManger::InitReCenter()
{
  string ProName;

  ProName = "p_mQEastVertical";
  p_mQEastVertical = new TProfile2D(ProName.c_str(),ProName.c_str(),6000, 1212000.5,1218000.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
  ProName = "p_mQEastHorizontal";
  p_mQEastHorizontal = new TProfile2D(ProName.c_str(),ProName.c_str(),6000, 1212000.5,1218000.5,9,-0.5,8.5);

  ProName = "p_mQWestVertical";
  p_mQWestVertical = new TProfile2D(ProName.c_str(),ProName.c_str(),6000, 1212000.5,1218000.5,9,-0.5,8.5);
  ProName = "p_mQWestHorizontal";
  p_mQWestHorizontal = new TProfile2D(ProName.c_str(),ProName.c_str(),6000, 1212000.5,1218000.5,9,-0.5,8.5);
}

void StZdcSmdProManger::FillReCenterEast(TVector2 qVector, int Cent9, int RunIndex) // vz_sign = vertex pos/neg 
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  // Event Plane method
  p_mQEastVertical->Fill((double)RunIndex,(double)Cent9,(double)qx);
  p_mQEastHorizontal->Fill((double)RunIndex,(double)Cent9,(double)qy);
}

void StZdcSmdProManger::FillReCenterWest(TVector2 qVector, int Cent9, int RunIndex) // vz_sign = vertex pos/neg 
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  p_mQWestVertical->Fill((double)RunIndex,(double)Cent9,(double)qx);
  p_mQWestHorizontal->Fill((double)RunIndex,(double)Cent9,(double)qy);
}

void StZdcSmdProManger::WriteReCenter()
{
  p_mQEastVertical->Write();
  p_mQEastHorizontal->Write();
  p_mQWestVertical->Write();
  p_mQWestHorizontal->Write();
}

//---------------------------------------------------------------------------------

void StZdcSmdProManger::InitShift()
{
  for(int i_shift = 0; i_shift < 20; ++i_shift) // Shift Order
  {
    string ProName;

    ProName = Form("p_mQEastCos_%d",i_shift);
    p_mQEastCos[i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),6000, 1212000.5,1218000.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
    ProName = Form("p_mQEastSin_%d",i_shift);
    p_mQEastSin[i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),6000, 1212000.5,1218000.5,9,-0.5,8.5);

    ProName = Form("p_mQWestCos_%d",i_shift);
    p_mQWestCos[i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),6000, 1212000.5,1218000.5,9,-0.5,8.5);
    ProName = Form("p_mQWestSin_%d",i_shift);
    p_mQWestSin[i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),6000, 1212000.5,1218000.5,9,-0.5,8.5);
  }
}

void StZdcSmdProManger::FillShiftEast(TVector2 qVector, int Cent9, int RunIndex)
{
  float Psi = TMath::ATan2(qVector.Y(),qVector.X());
  for(int i_shift = 0; i_shift < 20; ++i_shift)
  {
    p_mQEastCos[i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos((i_shift+1)*Psi));
    p_mQEastSin[i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin((i_shift+1)*Psi));
  }
}

void StZdcSmdProManger::FillShiftWest(TVector2 qVector, int Cent9, int RunIndex)
{
  float Psi = TMath::ATan2(qVector.Y(),qVector.X());
  for(int i_shift = 0; i_shift < 20; ++i_shift)
  {
    p_mQWestCos[i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos((i_shift+1)*Psi));
    p_mQWestSin[i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin((i_shift+1)*Psi));
  }
}

void StZdcSmdProManger::WriteShift()
{
  for(int i_shift = 0; i_shift < 20; ++i_shift) // Shift Order
  {
    p_mQEastCos[i_shift]->Write();
    p_mQEastSin[i_shift]->Write();
    p_mQWestCos[i_shift]->Write();
    p_mQWestSin[i_shift]->Write();
  }
}

//---------------------------------------------------------------------------------

void StZdcSmdProManger::InitShiftFull()
{
  for(int i_shift = 0; i_shift < 20; ++i_shift) // Shift Order
  {
    string ProName;

    ProName = Form("p_mQFullCos_%d",i_shift);
    p_mQFullCos[i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),6000, 1212000.5,1218000.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
    ProName = Form("p_mQFullSin_%d",i_shift);
    p_mQFullSin[i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),6000, 1212000.5,1218000.5,9,-0.5,8.5);
  }
}

void StZdcSmdProManger::FillShiftFull(TVector2 qVector, int Cent9, int RunIndex)
{
  float Psi = TMath::ATan2(qVector.Y(),qVector.X());
  for(int i_shift = 0; i_shift < 20; ++i_shift)
  {
    p_mQFullCos[i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos((i_shift+1)*Psi));
    p_mQFullSin[i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin((i_shift+1)*Psi));
  }
}

void StZdcSmdProManger::WriteShiftFull()
{
  for(int i_shift = 0; i_shift < 20; ++i_shift) // Shift Order
  {
    p_mQFullCos[i_shift]->Write();
    p_mQFullSin[i_shift]->Write();
  }
}
//---------------------------------------------------------------------------------

void StZdcSmdProManger::InitResolution()
{
  p_mResolution = new TProfile("p_mResolution","p_mResolution",9,-0.5,8.5);
}

void StZdcSmdProManger::FillResolution(TVector2 QEast, TVector2 QWest, int Cent9)
{
  float Psi_East = TMath::ATan2(QEast.Y(),QEast.X());
  float Psi_West = TMath::ATan2(QWest.Y(),QWest.X());
  float resolution = TMath::Cos(Psi_West-Psi_East+TMath::Pi());
  p_mResolution->Fill((double)Cent9,resolution);
}

void StZdcSmdProManger::WriteResolution()
{
  p_mResolution->Write();
}

//---------------------------------------------------------------------------------

void StZdcSmdProManger::InitDirectedFlow()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string ProName = Form("p_mDirectedFlow_%d",i_cent);
    p_mDirectedFlow[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),10,-1.0,1.0);
  }
  p_mDirectedFlowCom = new TProfile("p_mDirectedFlowCom","p_mDirectedFlowCom",10,-1.0,1.0);
}

void StZdcSmdProManger::FillDirectedFlow(int Cent9, float eta, float pt, float v1, float resolution, float reweight)
{
  if(pt > 0.15 && pt < 2.0)
  {
    p_mDirectedFlow[Cent9]->Fill(eta,v1/resolution,reweight);
    if(Cent9 >= 2 && Cent9 <= 4) p_mDirectedFlowCom->Fill(eta,v1/resolution,reweight); // 30-60%
  }
}

void StZdcSmdProManger::WriteDirectedFlow()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent) p_mDirectedFlow[i_cent]->Write();
  p_mDirectedFlowCom->Write();
}
//---------------------------------------------------------------------------------
