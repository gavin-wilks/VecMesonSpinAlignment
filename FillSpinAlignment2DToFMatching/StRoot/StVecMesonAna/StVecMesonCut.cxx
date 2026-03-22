#include "StRoot/StVecMesonAna/StVecMesonCut.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StMessMgr.h"

ClassImp(StVecMesonCut)

//---------------------------------------------------------------------------------

StVecMesonCut::StVecMesonCut(Int_t energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StVecMesonCut::~StVecMesonCut()
{
}

//---------------------------------------------------------------------------------
double StVecMesonCut::getRefMultReweight(double vz, int refMult)
{
  if(mEnergy == 4)
  { 
    if(vz > -5.0 && vz < 5.0) return refMult*vmsa::vz_corr[14];

    for(int ivz = 0; ivz < 14; ivz++)
    {
      if(vz <  5.0+(ivz+1)*10.0 && vz >=  5.0+ivz*10.0) return refMult*vmsa::vz_corr[15+ivz];
      if(vz > -5.0-(ivz+1)*10.0 && vz <= -5.0-ivz*10.0) return refMult*vmsa::vz_corr[13-ivz];
    }
  }
}

double StVecMesonCut::getEventWeight(int cent9, double refMult)
{
  if(mEnergy == 4)
  {
    if(cent9 >= 6) return 1.0;

    if(cent9 >= 0 && cent9 < 6) 
    {
      return 1.04433e+00+-1.13957e-01/(4.54889e-01*refMult + -3.43209e-01) + -1.55692e-03*(4.54889e-01*refMult + -3.43209e-01) + 4.00872e+00/TMath::Power(4.54889e-01*refMult+-3.43209e-01 ,2) + 1.03440e-05*TMath::Power(4.54889e-01*refMult+-3.43209e-01 ,2);
    }
  }
  return 1.0;
}


bool StVecMesonCut::passTrackEP(TLorentzVector lTrack, Float_t dca)
{
  // only used for primary track
  // dca cut for event plane reconstruction
  if(fabs(dca) > vmsa::mDcaEPMax[mEnergy])
  {
    return kFALSE;
  }

  // pt cut 0.2 - 2.0 GeV/c
  Float_t pt = lTrack.Perp();
  //Float_t p = lTrack.P();
  if(!(pt > vmsa::mPrimPtMin[mEnergy] && pt < vmsa::mPrimPtMax))
  {
    return kFALSE;
  }

  // eta cut -> only important for particles reconstructed by global tracks
  Float_t eta = lTrack.Eta();
  if(fabs(eta) > vmsa::mEtaMax)
  {
    return kFALSE;
  }

  return kTRUE;
}
//---------------------------------------------------------------------------------
bool StVecMesonCut::passTrackEtaEast(TLorentzVector lTrack) // neg
{
  Float_t eta = lTrack.Eta();
  
  // eta cut
  // eta_gap between two sub event plane is 2*mEta_Gap
  if(!(eta > -1.0*vmsa::mEtaMax && eta < -1.0*vmsa::mEta_Gap))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StVecMesonCut::passTrackEtaWest(TLorentzVector lTrack) // pos
{
  Float_t eta = lTrack.Eta();

  // eta cut
  // eta_gap between two sub event plane is 2*mEta_Gap
  if(!(eta > vmsa::mEta_Gap && eta < vmsa::mEtaMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StVecMesonCut::passEtaEast(TLorentzVector lTrack) // neg
{
  Float_t eta = lTrack.Eta();
  
  // eta cut
  // eta_gap between mother particle and sub event plane is mEta_Gap
  if(!(eta > -1.0*vmsa::mEtaMax && eta < 0.0))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StVecMesonCut::passEtaWest(TLorentzVector lTrack) // pos
{
  Float_t eta = lTrack.Eta();

  // eta cut
  // eta_gap between mother particle and sub event plane is mEta_Gap
  if(!(eta > 0.0 && eta < vmsa::mEtaMax))
  {
    return kFALSE;
  }

  return kTRUE;
}
//---------------------------------------------------------------------------------
bool StVecMesonCut::passTrackDcaSys(Float_t dcaA, Float_t dcaB, Int_t dcaSys, Int_t mMode) 
{
  if( !(fabs(dcaA) <= vmsa::mDcaSys[dcaSys] && fabs(dcaB) <= vmsa::mDcaSys[dcaSys]) ) return kFALSE;

  return kTRUE;
}

bool StVecMesonCut::passTrackSigSys(Float_t nsA, Float_t nsB, Int_t sigSys, Int_t mMode) 
{
  if( mMode == 0 && !(fabs(nsA) <= vmsa::mNSigmaKaonSys[sigSys] && fabs(nsB) <= vmsa::mNSigmaKaonSys[sigSys]) ) return kFALSE; // phi-meson  
  if( mMode == 1 && !(fabs(nsA) <= vmsa::mNSigmaPionSys[sigSys] && fabs(nsB) <= vmsa::mNSigmaPionSys[sigSys]) ) return kFALSE; // rho-meson
  if( mMode == 2 && !(fabs(nsA) <= vmsa::mNSigmaKaonSysKStar[sigSys] && fabs(nsB) <= vmsa::mNSigmaPionSys[sigSys]) ) return kFALSE; // Kstar-meson
  
  return kTRUE;
}

bool StVecMesonCut::passDipAngleCut(TVector3 kp, TVector3 km)
{
  double costheta = kp.Dot(km)/kp.Mag()/km.Mag();
  //costheta = TMath::Clamp(costheta,-1.0,1.0);
  double theta = TMath::ACos(costheta);
  if(theta < 0.04) return kFALSE;
  
  return kTRUE;
}
//---------------------------------------------------------------------------------
