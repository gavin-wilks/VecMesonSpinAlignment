#include "StRoot/StVecMesonAna/StVecMesonCut.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
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
  Float_t p = lTrack.P();
  if(!(pt > vmsa::mPrimPtMin[mEnergy] && pt < vmsa::mPrimPtMax && p < vmsa::mPrimMomMax))
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
  if( mMode == 2 && !(fabs(nsA) <= vmsa::mNSigmaKaonSys[sigSys] && fabs(nsB) <= vmsa::mNSigmaPionSys[sigSys]) ) return kFALSE; // Kstar-meson
  
  return kTRUE;
}
//---------------------------------------------------------------------------------
