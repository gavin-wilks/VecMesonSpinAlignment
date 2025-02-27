#include "StRoot/StVecMesonEpdMaker/StVecMesonEpdProManger.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TMath.h"

ClassImp(StVecMesonEpdProManger)

TString StVecMesonEpdProManger::mVStr[4] = {"neg70","neg30","pos30","pos70"};

//---------------------------------------------------------------------------------

StVecMesonEpdProManger::StVecMesonEpdProManger()
{
}

//---------------------------------------------------------------------------------

StVecMesonEpdProManger::~StVecMesonEpdProManger()
{
  /* */
}

//---------------------------------------------------------------------------------

void StVecMesonEpdProManger::InitReCenter()
{
  for(Int_t i = 0; i < 4; i++) // vertex pos/neg
  {
    TString ProName;
    ProName = Form("qx_2nd_Vertex_%s_East",mVStr[i].Data());
    p_mq2x_East_EP[i] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
    ProName = Form("qy_2nd_Vertex_%s_East",mVStr[i].Data());
    p_mq2y_East_EP[i] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
    ProName = Form("qx_2nd_Vertex_%s_West",mVStr[i].Data());
    p_mq2x_West_EP[i] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality
    ProName = Form("qy_2nd_Vertex_%s_West",mVStr[i].Data());
    p_mq2y_West_EP[i] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality

    ProName = Form("qx_2nd_Vertex_%s_Full",mVStr[i].Data());
    p_mq2x_Full_EP[i] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
    ProName = Form("qy_2nd_Vertex_%s_Full",mVStr[i].Data());
    p_mq2y_Full_EP[i] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
  }
}

void StVecMesonEpdProManger::FillTrackEast(TVector2 q2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Float_t pt) // i = vertex pos/neg 
{
  const Float_t q2x = q2Vector.X();
  const Float_t q2y = q2Vector.Y();

  Float_t w;
  if(pt <= vmsa::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > vmsa::mPrimPtWeight)
  {
    w = vmsa::mPrimPtWeight;
  }

  // Event Plane method
  p_mq2x_East_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x,(Double_t)w);
  p_mq2y_East_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y,(Double_t)w);
}

void StVecMesonEpdProManger::FillTrackWest(TVector2 q2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Float_t pt) // i = vertex pos/neg, j = eta_gap
{
  const Float_t q2x = q2Vector.X();
  const Float_t q2y = q2Vector.Y();

  Float_t w;
  if(pt <= vmsa::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > vmsa::mPrimPtWeight)
  {
    w = vmsa::mPrimPtWeight;
  }

  // Event Plane method
  p_mq2x_West_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x,(Double_t)w);
  p_mq2y_West_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y,(Double_t)w);
}

void StVecMesonEpdProManger::FillTrackFull(TVector2 q2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Float_t pt) // i = vertex pos/neg
{
  const Float_t q2x = q2Vector.X();
  const Float_t q2y = q2Vector.Y();

  Float_t w;
  if(pt <= vmsa::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > vmsa::mPrimPtWeight)
  {
    w = vmsa::mPrimPtWeight;
  }

  // Event Plane method
  p_mq2x_Full_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x,(Double_t)w);
  p_mq2y_Full_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y,(Double_t)w);
}

void StVecMesonEpdProManger::WriteReCenter()
{
  for(Int_t i = 0; i < 4; i++) // vertex pos/neg
  {
    p_mq2x_East_EP[i]->Write();
    p_mq2y_East_EP[i]->Write();
    p_mq2x_West_EP[i]->Write();
    p_mq2y_West_EP[i]->Write();

    p_mq2x_Full_EP[i]->Write();
    p_mq2y_Full_EP[i]->Write();
  }
}
//----------------------------------------------------------------------------

void StVecMesonEpdProManger::InitShift()
{
  for(Int_t i = 0; i < 4; i++) // vertex pos/neg
  {
    for(Int_t k = 0; k < 5; k++) // Shift Order
    {
      TString ProName;
      ProName = Form("CosPsi2_Vertex_%s_Order_%d_East",mVStr[i].Data(),k);
      p_mcos2_East_EP[i][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
      ProName = Form("SinPsi2_Vertex_%s_Order_%d_East",mVStr[i].Data(),k);
      p_msin2_East_EP[i][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
      ProName = Form("CosPsi2_Vertex_%s_Order_%d_West",mVStr[i].Data(),k);
      p_mcos2_West_EP[i][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
      ProName = Form("SinPsi2_Vertex_%s_Order_%d_West",mVStr[i].Data(),k);
      p_msin2_West_EP[i][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

      ProName = Form("CosPsi2_Vertex_%s_Order_%d_Full",mVStr[i].Data(),k);
      p_mcos2_Full_EP[i][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
      ProName = Form("SinPsi2_Vertex_%s_Order_%d_Full",mVStr[i].Data(),k);
      p_msin2_Full_EP[i][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
    }
  }
}

void StVecMesonEpdProManger::FillEventEast_EP(TVector2 Psi2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
  const Float_t cos2 = Psi2Vector.X();
  const Float_t sin2 = Psi2Vector.Y();

  p_mcos2_East_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
  p_msin2_East_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);
}

void StVecMesonEpdProManger::FillEventWest_EP(TVector2 Psi2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
  const Float_t cos2 = Psi2Vector.X();
  const Float_t sin2 = Psi2Vector.Y();

  p_mcos2_West_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
  p_msin2_West_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);
}

void StVecMesonEpdProManger::FillEventFull_EP(TVector2 Psi2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k) // i = vertex pos/neg, k = ShiftOrder
{
  const Float_t cos2 = Psi2Vector.X();
  const Float_t sin2 = Psi2Vector.Y();

  p_mcos2_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
  p_msin2_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);
}

void StVecMesonEpdProManger::WriteShift()
{
  for(Int_t i = 0; i < 4; i++) // vertex pos/neg
  {
    for(Int_t k = 0; k < 5; k++) // Shift Order
    {
      p_mcos2_East_EP[i][k]->Write();
      p_msin2_East_EP[i][k]->Write();
      p_mcos2_West_EP[i][k]->Write();
      p_msin2_West_EP[i][k]->Write();

      p_mcos2_Full_EP[i][k]->Write();
      p_msin2_Full_EP[i][k]->Write();
    }
  }
}

//----------------------------------------------------------------------------
void StVecMesonEpdProManger::InitResolution()
{
  p_mRes2_Sub = new TProfile("p_mRes2_Sub","p_mRes2_Sub",9,-0.5,8.5);
  p_mRes2_Ran = new TProfile("p_mRes2_Ran","p_mRes2_Ran",9,-0.5,8.5);
}

void StVecMesonEpdProManger::FillRes_Sub(Int_t Centrality, Float_t Psi2_East, Float_t Psi2_West)
{
  Float_t Res2_Sub = TMath::Cos(2.0*(Psi2_East-Psi2_West));
  p_mRes2_Sub->Fill(Centrality,Res2_Sub);
}

void StVecMesonEpdProManger::FillRes_Ran(Int_t Centrality, Float_t Psi2_RanA, Float_t Psi2_RanB)
{
  Float_t Res2_Ran = TMath::Cos(2.0*(Psi2_RanA-Psi2_RanB));
  p_mRes2_Ran->Fill(Centrality,Res2_Ran);
}

void StVecMesonEpdProManger::WriteResolution()
{
  p_mRes2_Sub->Write();
  p_mRes2_Ran->Write();
}
//----------------------------------------------------------------------------

void StVecMesonEpdProManger::InitD12()
{
  p_mD12 = new TProfile("p_mD12","p_mD12",9,-0.5,8.5);
}

void StVecMesonEpdProManger::FillD12(int cent9, double Psi1, double Psi2)
{
  double D12 = TMath::Cos(2*(Psi1-Psi2));
  p_mD12->Fill(cent9,D12);
}

void StVecMesonEpdProManger::WriteD12()
{
  p_mD12->Write();
}
