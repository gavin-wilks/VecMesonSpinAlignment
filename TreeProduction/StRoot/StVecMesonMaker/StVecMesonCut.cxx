#include "StRoot/StVecMesonMaker/StVecMesonCut.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StMessMgr.h"

ClassImp(StVecMesonCut)

StRefMultCorr* StVecMesonCut::mRefMultCorr = NULL;
//---------------------------------------------------------------------------------

StVecMesonCut::StVecMesonCut(const int energy)
{
  mEnergy = energy;
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
}

//---------------------------------------------------------------------------------

StVecMesonCut::~StVecMesonCut()
{
  /* */
}

//---------------------------------------------------------------------------------

bool StVecMesonCut::isMinBias(StPicoEvent *picoEvent)
{
  // std::cout << "year: " << picoEvent->year() << std::endl;
  // std::cout << "day: " << picoEvent->day() << std::endl;
  // std::cout << "triggerIds: " << picoEvent->triggerIds()[0] << std::endl;  
  if(mEnergy == 0 && vmsa::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(640001) || picoEvent->isTrigger(640011) || picoEvent->isTrigger(640021) || picoEvent->isTrigger(640031) || picoEvent->isTrigger(640041) || picoEvent->isTrigger(640051) )) return false; // 7.7GeV_2019
  if(mEnergy == 1 && vmsa::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(640001) || picoEvent->isTrigger(640011) || picoEvent->isTrigger(640021) || picoEvent->isTrigger(640031) || picoEvent->isTrigger(640041) || picoEvent->isTrigger(640051) )) return false; // 9.1GeV_2019
  if(mEnergy == 2 && vmsa::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(640001) || picoEvent->isTrigger(640011) || picoEvent->isTrigger(640021) || picoEvent->isTrigger(640031) || picoEvent->isTrigger(640041) || picoEvent->isTrigger(640051) )) return false; // 11.5GeV_2019
  if(mEnergy == 3 && vmsa::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(650001) || picoEvent->isTrigger(650011) || picoEvent->isTrigger(650021) || picoEvent->isTrigger(650031) || picoEvent->isTrigger(650041) || picoEvent->isTrigger(650051) )) return false; // 14.6GeV_2019
  if(mEnergy == 4 && vmsa::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(640001) || picoEvent->isTrigger(640011) || picoEvent->isTrigger(640021) || picoEvent->isTrigger(640031) || picoEvent->isTrigger(640041) || picoEvent->isTrigger(640051) )) return false; // 19.6GeV_2019
  return true;
}

bool StVecMesonCut::isPileUpEvent(int refMult, int numOfBTofMatch)
{
  if(mEnergy == 0 || mEnergy == 1 || mEnergy == 2 || mEnergy == 3)
  {        
    if(numOfBTofMatch > vmsa::mMatchedToFMin) return kFALSE;         
  }

  if(mEnergy == 4)
  {
    if(numOfBTofMatch <= vmsa::mMatchedToFMin) return kTRUE; 
        
    // 5th order polynomial coefficients for lower cut, l
    double p0l = -13.11    ;   
    double p1l = 0.8207    ;   
    double p2l = -4.241e-3 ;
    double p3l = 2.81e-5   ;   
    double p4l = -6.434e-8 ;
    double p5l = 4.833e-11 ;
    // 5th order polynomial coefficients for higher cut, h
    double p0h = 10.07     ;   
    double p1h = 1.417     ;   
    double p2h = 1.979e-4  ;
    double p3h = -4.87e-6  ;
    double p4h = 1.109e-8  ;
    double p5h = -1.111e-11;

    double refLow  = p0l + p1l*numOfBTofMatch + p2l*pow(numOfBTofMatch,2) + p3l*pow(numOfBTofMatch,3) + p4l*pow(numOfBTofMatch,4) + p5l*pow(numOfBTofMatch,5);
    double refHigh = p0h + p1h*numOfBTofMatch + p2h*pow(numOfBTofMatch,2) + p3h*pow(numOfBTofMatch,3) + p4h*pow(numOfBTofMatch,4) + p5h*pow(numOfBTofMatch,5); 
   
    if(refMult > refLow && refMult < refHigh) return kFALSE;
  }     

  return kTRUE;
}

bool StVecMesonCut::passEventCut(StPicoDst *pico)
{
  // initialize mMatchedToF
  mMatchedToF = 0;
  mN_prim = 0;
  mN_non_prim = 0;

  StPicoEvent *event = pico->event();
  if(!event)
  {
    return kFALSE;
  }

  // initialize StRefMultCorr
  const Int_t runId = event->runId();
  const Int_t refMult = event->refMult();
  const unsigned short numOfBTofMatch = event->nBTOFMatch();
  const Float_t vx = event->primaryVertex().x();
  const Float_t vy = event->primaryVertex().y();
  const Float_t vz = event->primaryVertex().z();
  //const Float_t zdcX = event->ZDCx();
  const Float_t vzVpd = event->vzVpd();
  //const Bool_t isBES = (event->energy() < 200.);
  mRefMultCorr->init(runId);

  // StRefMultCorr bad run cut
  if(mRefMultCorr->isBadRun(runId))
  {
    return kFALSE;
  }

  // minBias event cut
  if(isMinBias(event))
  {
    return kFALSE;
  }

  // event vertex cut
  // vz cut
  if(fabs(vz) > vmsa::mVzMaxMap[mEnergy])
  {
    return kFALSE;
  }
  // vr cut
  if(sqrt(vx*vx+vy*vy) > vmsa::mVrMax)
  {
    return kFALSE;
  }
  // vz-vzVpd cut for 200 GeV
  //if(!isBES)
  //{
  if(fabs(vz-vzVpd) > vmsa::mVzVpdDiffMax)
  {
    return kFALSE;
  }
  //}

  // refMult (0-80%) cut
  //if(!isBES) mRefMultCorr->initEvent(refMult,vz,zdcX); // 200GeV
  //if(isBES) 
  mRefMultCorr->initEvent(refMult,vz); // BES
  //if(!mRefMultCorr->isRefMultOk())
  //{
  //  return kFALSE;
  //}

  // ToF matched points cut
  /*Int_t nMatchedToF = 0;
  Int_t nN_prim = 0;
  Int_t nN_non_prim = 0;
  const Int_t nTracks = pico->numberOfTracks();
  for(Int_t i = 0; i < nTracks; i++)
  {
    StPicoTrack *track = (StPicoTrack*)pico->track(i);
    if(!track)
    {
      continue;
    }
    // stop loop if already have enough TOF matched points
//    if(nMatchedToF >= vmsa::mMatchedToFMin)
//    {
//      return kTRUE;
//    }
    if(track->gDCA(vx,vy,vz) > 3) // global track
    {
      nN_non_prim++;
    }
    else
    {
      nN_prim++;
      if(track->btofMatchFlag() > 0 && track->btof() != 0 && track->btofBeta() != 0)
      {
	nMatchedToF++;
      }
    }
  }

  mMatchedToF = nMatchedToF;
  mN_prim = nN_prim;
  mN_non_prim = nN_non_prim;


  if(nMatchedToF < vmsa::mMatchedToFMin)
  {
    return kFALSE;
  }
  */
  // pileUpEvent cut
  if(isPileUpEvent(refMult,numOfBTofMatch))
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

Int_t StVecMesonCut::getMatchedToF()
{
  return mMatchedToF;
}

Int_t StVecMesonCut::getNpirm()
{
  return mN_prim;
}

Int_t StVecMesonCut::getNnonprim()
{
  return mN_non_prim;
}
//---------------------------------------------------------------------------------

Float_t StVecMesonCut::getBeta(StPicoTrack *picoTrack, StPicoDst *picoDst)
{
  double beta = -999.9;
  // StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    beta = tofTrack->btofBeta();
  }

  return beta;
}

Float_t StVecMesonCut::getPrimaryMass2(StPicoTrack *picoTrack, StPicoDst *picoDst)
{
  double mass2 = -999.9;
  // StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    const double beta = tofTrack->btofBeta();
    // const double Momentum = picoTrack->pMom().mag(); // primary momentum for 54GeV_2017
    // const double Momentum = picoTrack->pMom().Mag(); // primary momentum for 27GeV_2018
    TVector3 primMom; // temp fix for StThreeVectorF & TVector3
    const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
    const double primPy    = picoTrack->pMom().y();
    const double primPz    = picoTrack->pMom().z();
    primMom.SetXYZ(primPx,primPy,primPz);
    const double Momentum = primMom.Mag(); // primary momentum

    if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && beta != 0)
    {
      mass2 = Momentum*Momentum*(1.0/(beta*beta) - 1.0);
    }
  }

  return mass2;
}

Float_t StVecMesonCut::getV0Mass2(StPicoTrack *picoTrack, StPicoDst *picoDst)
{
  double mass2 = -999.9;
  // StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    const double beta = tofTrack->btofBeta();
    // const double Momentum = picoTrack->pMom().mag(); // primary momentum for 54GeV_2017
    // const double Momentum = picoTrack->pMom().Mag(); // primary momentum for 27GeV_2018
    TVector3 primMom; // temp fix for StThreeVectorF & TVector3
    const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
    const double primPy    = picoTrack->pMom().y();
    const double primPz    = picoTrack->pMom().z();
    primMom.SetXYZ(primPx,primPy,primPz);
    const double Momentum = primMom.Mag(); // primary momentum

    if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && beta != 0)
    {
      mass2 = Momentum*Momentum*(1.0/(beta*beta) - 1.0);
    }
  }

  return mass2;
}

bool StVecMesonCut::passSigPionCut(StPicoTrack* track, Float_t scale_nSigma_factor)
{
  Float_t nSigmaPion = track->nSigmaPion();
  if(fabs(nSigmaPion*scale_nSigma_factor) > vmsa::mNSigmaPionMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

bool StVecMesonCut::passSigKaonCut(StPicoTrack* track, Float_t scale_nSigma_factor)
{
  Float_t nSigmaKaon = track->nSigmaKaon();
  if(fabs(nSigmaKaon*scale_nSigma_factor) > vmsa::mNSigmaKaonMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

bool StVecMesonCut::passSigProntonCut(StPicoTrack* track, Float_t scale_nSigma_factor)
{
  Float_t nSigmaProton = track->nSigmaProton();
  if(fabs(nSigmaProton*scale_nSigma_factor) > vmsa::mNSigmaProtonMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StVecMesonCut::passTrackBasic(StPicoTrack *track)
{
  // nHitsFit cut
  if(track->nHitsFit() < vmsa::mHitsFitTPCMin)
  {
    return kFALSE;
  }

  // nHitsRatio cut
  if(track->nHitsMax() <= vmsa::mHitsMaxTPCMin)
  {
    return kFALSE;
  }
  if((Float_t)track->nHitsFit()/(Float_t)track->nHitsMax() < vmsa::mHitsRatioTPCMin)
  {
    return kFALSE;
  }

  // eta cut
  Float_t eta = track->pMom().PseudoRapidity();
  if(fabs(eta) > vmsa::mEtaMax)
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StVecMesonCut::passTrackEP(StPicoTrack *track, StPicoEvent *picoEvent)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track)) return kFALSE;

  const double vx = picoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const double vy = picoEvent->primaryVertex().y();
  const double vz = picoEvent->primaryVertex().z();

  // dca cut for event plane reconstruction: 200GeV = 3.0, BES = 1.0
  if(track->gDCA(vx,vy,vz) > vmsa::mDcaEPMax[mEnergy])
  {
    return kFALSE;
  }  

  // pt cut 0.2 - 2.0 GeV/c
  Float_t pt = track->pMom().Perp();
  Float_t p  = track->pMom().Mag();
  if(!(pt > vmsa::mPrimPtMin[mEnergy] && pt < vmsa::mPrimPtMax && p < vmsa::mPrimMomMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------
bool StVecMesonCut::passTrackPhi(StPicoTrack *track, StPicoEvent *picoEvent)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track)) return kFALSE;

  const double vx = picoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const double vy = picoEvent->primaryVertex().y();
  const double vz = picoEvent->primaryVertex().z();

  // dca cut for event plane reconstruction: 200GeV = 3.0, BES = 1.0
  if(track->gDCA(vx,vy,vz) > vmsa::mDcaEPMax[mEnergy])
  {
    return kFALSE;
  }

  // primary pt and momentum cut: PtMin = 0.1
  if(!(track->pMom().Perp() > vmsa::mGlobPtMin && track->pMom().Mag() < vmsa::mPrimMomMax))
  {
    return kFALSE;
  }

  return kTRUE;
}
