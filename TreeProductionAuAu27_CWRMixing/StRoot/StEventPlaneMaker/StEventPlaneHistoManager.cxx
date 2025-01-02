#include <TH2F.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TMath.h>
#include <TString.h>
#include <TProfile2D.h>
#include "StRoot/StEventPlaneMaker/StEventPlaneHistoManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StEpdUtil/StEpdEpInfo.h"

ClassImp(StEventPlaneHistoManager)

//-------------------------------------------------------------------------------------------

StEventPlaneHistoManager::StEventPlaneHistoManager()
{
}

//-------------------------------------------------------------------------------------------

StEventPlaneHistoManager::~StEventPlaneHistoManager()
{
  /* */
}

//-------------------------------------------------------------------------------------------
// ZDC EP
void StEventPlaneHistoManager::initZdcGainCorr()
{
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	string HistName = Form("h_mZdcGainCorr%s%s_%d",recoEP::mEastWest[i_eastwest].c_str(),recoEP::mVertHori[i_verthori].c_str(),i_slat);
	h_mZdcGainCorr[i_eastwest][i_verthori][i_slat] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,5000,-4.5,4995.5);
      }
    }
  }
}

void StEventPlaneHistoManager::fillZdcGainCorr(int i_eastwest, int i_verthori, int i_slat, int runIndex, float zdcsmd)
{
  h_mZdcGainCorr[i_eastwest][i_verthori][i_slat]->Fill((double)runIndex,zdcsmd);
}

void StEventPlaneHistoManager::writeZdcGainCorr()
{
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	h_mZdcGainCorr[i_eastwest][i_verthori][i_slat]->Write();
      }
    }
  }
}

// raw ZDC-SMD EP
void StEventPlaneHistoManager::initZdcRawEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mZdcRawEpEast_%d",i_cent);
    h_mZdcRawEpEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcRawEpWest_%d",i_cent);
    h_mZdcRawEpWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcRawEpFull_%d",i_cent);
    h_mZdcRawEpFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StEventPlaneHistoManager::fillZdcRawSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex)
{
  double PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mZdcRawEpEast[Cent9]->Fill(runIndex,PsiEast);
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcRawEpWest[Cent9]->Fill(runIndex,PsiWest);
}

void StEventPlaneHistoManager::fillZdcRawFullEP(TVector2 QFull, int Cent9, int runIndex)
{
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcRawEpFull[Cent9]->Fill(runIndex,PsiFull);
}

void StEventPlaneHistoManager::writeZdcRawEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mZdcRawEpEast[i_cent]->Write();
    h_mZdcRawEpWest[i_cent]->Write();
    h_mZdcRawEpFull[i_cent]->Write();
  }
}

// recenter ZDC-SMD EP
void StEventPlaneHistoManager::initZdcReCenterEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mZdcReCenterEpEast_%d",i_cent);
    h_mZdcReCenterEpEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcReCenterEpWest_%d",i_cent);
    h_mZdcReCenterEpWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcReCenterEpFull_%d",i_cent);
    h_mZdcReCenterEpFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StEventPlaneHistoManager::fillZdcReCenterSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex)
{
  double PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mZdcReCenterEpEast[Cent9]->Fill(runIndex,PsiEast);
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcReCenterEpWest[Cent9]->Fill(runIndex,PsiWest);
}

void StEventPlaneHistoManager::fillZdcReCenterFullEP(TVector2 QFull, int Cent9, int runIndex)
{
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcReCenterEpFull[Cent9]->Fill(runIndex,PsiFull);
}

void StEventPlaneHistoManager::writeZdcReCenterEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mZdcReCenterEpEast[i_cent]->Write();
    h_mZdcReCenterEpWest[i_cent]->Write();
    h_mZdcReCenterEpFull[i_cent]->Write();
  }
}

// shift ZDC-SMD EP
void StEventPlaneHistoManager::initZdcShiftEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mZdcShiftEpEast_%d",i_cent);
    h_mZdcShiftEpEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcShiftEpWest_%d",i_cent);
    h_mZdcShiftEpWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcShiftEpDiff_%d",i_cent);
    h_mZdcShiftEpDiff[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcShiftEpFull_%d",i_cent);
    h_mZdcShiftEpFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StEventPlaneHistoManager::fillZdcShiftSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex)
{
  double PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mZdcShiftEpEast[Cent9]->Fill(runIndex,PsiEast);
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcShiftEpWest[Cent9]->Fill(runIndex,PsiWest);
}

void StEventPlaneHistoManager::fillZdcShiftFullEP(TVector2 QDiff, TVector2 QFull, int Cent9, int runIndex)
{
  double PsiDiff = TMath::ATan2(QDiff.Y(),QDiff.X()); h_mZdcShiftEpDiff[Cent9]->Fill(runIndex,PsiDiff);
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcShiftEpFull[Cent9]->Fill(runIndex,PsiFull);
}

void StEventPlaneHistoManager::writeZdcShiftEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mZdcShiftEpEast[i_cent]->Write();
    h_mZdcShiftEpWest[i_cent]->Write();
    h_mZdcShiftEpDiff[i_cent]->Write();
    h_mZdcShiftEpFull[i_cent]->Write();
  }
}
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
// TPC EP
// raw TPC EP
//
void StEventPlaneHistoManager::initTpcBRQA()
{
  string ProName = Form("p_mTpcQ2xEast");
  p_mTpcQ2xEast = new TProfile2D(ProName.c_str(),ProName.c_str(),vmsa::mNumOfRunIndex,-0.5,(float)vmsa::mNumOfRunIndex-0.5,9,-0.5,8.5);
  ProName = Form("p_mTpcQ2yEast");
  p_mTpcQ2yEast= new TProfile2D(ProName.c_str(),ProName.c_str(),vmsa::mNumOfRunIndex,-0.5,(float)vmsa::mNumOfRunIndex-0.5,9,-0.5,8.5);
  ProName = Form("p_mTpcQ2xWest");
  p_mTpcQ2xWest = new TProfile2D(ProName.c_str(),ProName.c_str(),vmsa::mNumOfRunIndex,-0.5,(float)vmsa::mNumOfRunIndex-0.5,9,-0.5,8.5);
  ProName = Form("p_mTpcQ2yWest");
  p_mTpcQ2yWest= new TProfile2D(ProName.c_str(),ProName.c_str(),vmsa::mNumOfRunIndex,-0.5,(float)vmsa::mNumOfRunIndex-0.5,9,-0.5,8.5);
}

void StEventPlaneHistoManager::fillTpcBRQAEast(int runIndex, int cent9, TVector2 Q2)
{
  p_mTpcQ2xEast->Fill(runIndex,cent9,Q2.X());
  p_mTpcQ2yEast->Fill(runIndex,cent9,Q2.Y());
}

void StEventPlaneHistoManager::fillTpcBRQAWest(int runIndex, int cent9, TVector2 Q2)
{
  p_mTpcQ2xWest->Fill(runIndex,cent9,Q2.X());
  p_mTpcQ2yWest->Fill(runIndex,cent9,Q2.Y());
}

void StEventPlaneHistoManager::writeTpcBRQA()
{
  p_mTpcQ2xEast->Write();
  p_mTpcQ2yEast->Write();
  p_mTpcQ2xWest->Write();
  p_mTpcQ2yWest->Write();
}

void StEventPlaneHistoManager::initTpcRawEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    for(int order = 1; order <= 3; ++order)
    {
      string HistName = Form("h_mTpcRawEpEast_%d_%d",order,i_cent);
      h_mTpcRawEpEast[order-1][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::mNumOfRunIndex,-0.5,(double)vmsa::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
      HistName = Form("h_mTpcRawEpWest_%d_%d",order,i_cent);
      h_mTpcRawEpWest[order-1][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::mNumOfRunIndex,-0.5,(double)vmsa::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
      HistName = Form("h_mTpcRawEpFull_%d_%d",order,i_cent);
      h_mTpcRawEpFull[order-1][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::mNumOfRunIndex,-0.5,(double)vmsa::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    }
  }
}

void StEventPlaneHistoManager::fillTpcRawSubEP(int order, TVector2 QEast, TVector2 QWest, int Cent9, int runIndex)
{
  double PsiEast = (1.0/double(order))*TMath::ATan2(QEast.Y(),QEast.X()); h_mTpcRawEpEast[order-1][Cent9]->Fill(runIndex,PsiEast);
  double PsiWest = (1.0/double(order))*TMath::ATan2(QWest.Y(),QWest.X()); h_mTpcRawEpWest[order-1][Cent9]->Fill(runIndex,PsiWest);
}

void StEventPlaneHistoManager::fillTpcRawFullEP(int order, TVector2 QFull, int Cent9, int runIndex)
{
  double PsiFull = (1.0/double(order))*TMath::ATan2(QFull.Y(),QFull.X()); h_mTpcRawEpFull[order-1][Cent9]->Fill(runIndex,PsiFull);
}

void StEventPlaneHistoManager::writeTpcRawEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    for(int order = 1; order <= 3; ++order)
    {
      h_mTpcRawEpEast[order-1][i_cent]->Write();
      h_mTpcRawEpWest[order-1][i_cent]->Write();
      h_mTpcRawEpFull[order-1][i_cent]->Write();
    }
  }
}

// recenter TPC EP
void StEventPlaneHistoManager::initTpcReCenterEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    for(int order = 1; order <=3; ++order)
    {
      string HistName = Form("h_mTpcReCenterEpEast_%d_%d",order,i_cent);
      h_mTpcReCenterEpEast[order-1][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::mNumOfRunIndex,-0.5,(double)vmsa::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
      HistName = Form("h_mTpcReCenterEpWest_%d_%d",order,i_cent);
      h_mTpcReCenterEpWest[order-1][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::mNumOfRunIndex,-0.5,(double)vmsa::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
      HistName = Form("h_mTpcReCenterEpFull_%d_%d",order,i_cent);
      h_mTpcReCenterEpFull[order-1][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::mNumOfRunIndex,-0.5,(double)vmsa::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    }
  }
}

void StEventPlaneHistoManager::fillTpcReCenterSubEP(int order, TVector2 QEast, TVector2 QWest, int Cent9, int runIndex)
{
  double PsiEast = (1.0/double(order))*TMath::ATan2(QEast.Y(),QEast.X()); h_mTpcReCenterEpEast[order-1][Cent9]->Fill(runIndex,PsiEast);
  double PsiWest = (1.0/double(order))*TMath::ATan2(QWest.Y(),QWest.X()); h_mTpcReCenterEpWest[order-1][Cent9]->Fill(runIndex,PsiWest);
}

void StEventPlaneHistoManager::fillTpcReCenterFullEP(int order, TVector2 QFull, int Cent9, int runIndex)
{
  double PsiFull = (1.0/double(order))*TMath::ATan2(QFull.Y(),QFull.X()); h_mTpcReCenterEpFull[order-1][Cent9]->Fill(runIndex,PsiFull);
}

void StEventPlaneHistoManager::writeTpcReCenterEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    for(int order = 1; order <= 3; ++order)
    {
      h_mTpcReCenterEpEast[order-1][i_cent]->Write();
      h_mTpcReCenterEpWest[order-1][i_cent]->Write();
      h_mTpcReCenterEpFull[order-1][i_cent]->Write();
    }
  }
}

// shift TPC EP
void StEventPlaneHistoManager::initTpcShiftEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    for(int order = 1; order <= 3; ++order)
    {
      string HistName = Form("h_mTpcShiftEpEast_%d_%d",order,i_cent);
      h_mTpcShiftEpEast[order-1][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::mNumOfRunIndex,-0.5,(double)vmsa::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
      HistName = Form("h_mTpcShiftEpWest_%d_%d",order,i_cent);
      h_mTpcShiftEpWest[order-1][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::mNumOfRunIndex,-0.5,(double)vmsa::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
      HistName = Form("h_mTpcShiftEpRanA_%d_%d",order,i_cent);
      h_mTpcShiftEpRanA[order-1][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::mNumOfRunIndex,-0.5,(double)vmsa::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
      HistName = Form("h_mTpcShiftEpRanB_%d_%d",order,i_cent);
      h_mTpcShiftEpRanB[order-1][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::mNumOfRunIndex,-0.5,(double)vmsa::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
      HistName = Form("h_mTpcShiftEpFull_%d_%d",order,i_cent);
      h_mTpcShiftEpFull[order-1][i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),vmsa::mNumOfRunIndex,-0.5,(double)vmsa::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    }
  }
}

void StEventPlaneHistoManager::fillTpcShiftSubEP(int order, double PsiEast, double PsiWest, int Cent9, int runIndex)
{
  h_mTpcShiftEpEast[order-1][Cent9]->Fill(runIndex,PsiEast);
  h_mTpcShiftEpWest[order-1][Cent9]->Fill(runIndex,PsiWest);
}

void StEventPlaneHistoManager::fillTpcShiftRanEP(int order, double PsiRanA, double PsiRanB, int Cent9, int runIndex)
{
  h_mTpcShiftEpRanA[order-1][Cent9]->Fill(runIndex,PsiRanA);
  h_mTpcShiftEpRanB[order-1][Cent9]->Fill(runIndex,PsiRanB);
}

void StEventPlaneHistoManager::fillTpcShiftFullEP(int order, double PsiFull, int Cent9, int runIndex)
{
  h_mTpcShiftEpFull[order-1][Cent9]->Fill(runIndex,PsiFull);
}

void StEventPlaneHistoManager::writeTpcShiftEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    for(int order = 1; order <= 3; ++order)
    {
      h_mTpcShiftEpEast[order-1][i_cent]->Write();
      h_mTpcShiftEpWest[order-1][i_cent]->Write();
      h_mTpcShiftEpRanA[order-1][i_cent]->Write();
      h_mTpcShiftEpRanB[order-1][i_cent]->Write();
      h_mTpcShiftEpFull[order-1][i_cent]->Write();
    }
  }
}

void StEventPlaneHistoManager::initEpdEp()
{
  // histograms and profiles of event planes
  // EPD "psi" (first index for Q histograms): 0 = full psi | 1 = east psi | 2 = west psi
  // EPD "correction level" (last index):  0 = raw | 1 = weighted | 2 = weighted and shifted
  std::string corrName[3] = {"Raw", "Pw", "PwS"};
  std::string corrTitle[3] = {"Raw", "#phi-weighted", "#phi-weighted & Shifted"}; 

  for (int  icor=0; icor<3; icor++)
  {
    for (int order=1; order<=vmsa::mEpdEpOrder; order++)
    {
      double PhiMax = 2.0*TMath::Pi()/((double)order);
      h_mEpdEwPsi[order-1][icor]  = new TH2F(Form("EpdEwPsi%d%s",order,corrName[icor].c_str()),  Form("EPD %s #Psi_{%d,W} vs #Psi_{%d,E}",corrTitle[icor].c_str(),order,order),  100,0.0,PhiMax,100,0.0,PhiMax);
      
      h_mEpdFullPsi[order-1][icor]  = new TH1F(Form("EpdFullPsi%d%s",order,corrName[icor].c_str()), Form("EPD %s #Psi_{%d}",corrTitle[icor].c_str(),order), 100,0.0,PhiMax);
      
      h_mEpdEwPsi[order-1][icor]->GetXaxis()->SetTitle(Form("#Psi_{%d,E}",order));                 
      h_mEpdEwPsi[order-1][icor]->GetYaxis()->SetTitle(Form("#Psi_{%d,W}",order));
      h_mEpdFullPsi[order-1][icor]->GetXaxis()->SetTitle(Form("#Psi_{%d}",order));
    } 
    //h_mEpdAveCosD12[icor] = new TProfile(Form("EpdD12Cos%s",corrName[icor].c_str()),Form("EPD D12 %s #LT cos(2#Psi_{1^{st}}-2#Psi_{2^{nd}}) #GT",corrTitle[icor].c_str()),9,-0.5,8.5); 

  }
}
 
void StEventPlaneHistoManager::fillEpdEp(StEpdEpInfo result, int CentId)
{
  for (int order=1; order<=vmsa::mEpdEpOrder; order++)
  {
    h_mEpdEwPsi[order-1][0]->Fill(result.EastRawPsi(order),result.WestRawPsi(order));
    h_mEpdEwPsi[order-1][1]->Fill(result.EastPhiWeightedPsi(order),result.WestPhiWeightedPsi(order));
    h_mEpdEwPsi[order-1][2]->Fill(result.EastPhiWeightedAndShiftedPsi(order),result.WestPhiWeightedAndShiftedPsi(order));
    h_mEpdFullPsi[order-1][0]->Fill(result.FullRawPsi(order));
    h_mEpdFullPsi[order-1][1]->Fill(result.FullPhiWeightedPsi(order));
    h_mEpdFullPsi[order-1][2]->Fill(result.FullPhiWeightedAndShiftedPsi(order));
  }

  //h_mEpdAveCosD12[0]->Fill(double(CentId),cos((double)2.0*(result.FullRawPsi(1)-result.FullRawPsi(2))));
  //h_mEpdAveCosD12[1]->Fill(double(CentId),cos((double)2.0*(result.FullPhiWeightedPsi(1)-result.FullPhiWeightedPsi(2))));
  //h_mEpdAveCosD12[2]->Fill(double(CentId),cos((double)2.0*(result.FullPhiWeightedAndShiftedPsi(1)-result.FullPhiWeightedAndShiftedPsi(2))));

}

void StEventPlaneHistoManager::writeEpdEp()
{
  for (int icor=0; icor<3; icor++)
  {
    for (int order=1; order<=vmsa::mEpdEpOrder; order++)
    { 
      h_mEpdEwPsi[order-1][icor]->Write();
      h_mEpdFullPsi[order-1][icor]->Write(); 
    }
    //h_mEpdAveCosD12[icor]->Write();
  } 
}
//-------------------------------------------------------------------------------------------
void StEventPlaneHistoManager::initEpdQ()
{
  // histograms and profiles of event planes
  // EPD "psi" (first index for Q histograms): 0 = east psi | 1 = west psi
  // EPD "correction level" (last index):  0 = raw | 1 = phi-weighted 
  std::string corrName[2] = {"Raw", "Pw"};
  std::string corrTitle[2] = {"Raw", "#phi-weighted"}; 
  std::string psiName[2] = {"East","West"};  

  for (int  icor=0; icor<2; icor++)
  {
    for (int order=1; order<=vmsa::mEpdEpOrder; order++)
    {
      for (int psi = 0; psi < 2; psi++)
      { 
        h_mEpdQyQx[psi][order-1][icor] = new TH2F(Form("EpdQyQx%s%d%s",psiName[psi].c_str(),order,corrName[icor].c_str()),  Form("EPD %s %s%d #Q_{y} vs #Q_{x}",corrTitle[icor].c_str(),psiName[psi].c_str(),order),  100,-1.5,1.5,100,-1.5,1.5);
        h_mEpdQyQx[psi][order-1][icor]->GetXaxis()->SetTitle("Q_{x}");
        h_mEpdQyQx[psi][order-1][icor]->GetYaxis()->SetTitle("Q_{y}");
      }
    } 
  }
}
 
void StEventPlaneHistoManager::fillEpdQ(StEpdEpInfo result, int CentId)
{
  for (int order=1; order<=vmsa::mEpdEpOrder; order++)
  {   
    h_mEpdQyQx[0][order-1][0]->Fill(result.EastRawQ(order).X(),result.EastRawQ(order).Y());
    h_mEpdQyQx[0][order-1][1]->Fill(result.EastPhiWeightedQ(order).X(),result.EastPhiWeightedQ(order).Y()); 
    h_mEpdQyQx[1][order-1][0]->Fill(result.WestRawQ(order).X(),result.WestRawQ(order).Y());
    h_mEpdQyQx[1][order-1][1]->Fill(result.WestPhiWeightedQ(order).X(),result.WestPhiWeightedQ(order).Y()); 
  }
}

void StEventPlaneHistoManager::writeEpdQ()
{
  for (int icor=0; icor<2; icor++)
  {
    for (int order=1; order<=vmsa::mEpdEpOrder; order++)
    { 
      for (int psi = 0; psi < 2; psi++) // psi = 0 (east) | psi = 1 (west)
      {
        h_mEpdQyQx[psi][order-1][icor]->Write();
      }
    }
  } 
}


