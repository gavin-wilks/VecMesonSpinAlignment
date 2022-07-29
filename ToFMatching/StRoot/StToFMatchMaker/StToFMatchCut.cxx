#include "StRoot/StToFMatchMaker/StToFMatchCut.h"
#include "StRoot/StToFMatchMaker/StToFMatchCons.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoETofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StMessMgr.h"

ClassImp(StToFMatchCut)

StRefMultCorr* StToFMatchCut::mRefMultCorr = NULL;
//---------------------------------------------------------------------------------

StToFMatchCut::StToFMatchCut(Int_t energy)
{
  mEnergy = energy;
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
}

//---------------------------------------------------------------------------------

StToFMatchCut::~StToFMatchCut()
{
  /* */
}

//---------------------------------------------------------------------------------

bool StToFMatchCut::isMinBias(StPicoEvent *picoEvent)
{
  //std::cout << "year: " << picoEvent->year() << std::endl;
  //std::cout << "day: " << picoEvent->day() << std::endl;
  //std::cout << "triggerIds: " << picoEvent->triggerIds()[0] << std::endl;  
  if(mEnergy == 0 && vmsa::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(640001) || picoEvent->isTrigger(640011) || picoEvent->isTrigger(640021) || picoEvent->isTrigger(640031) || picoEvent->isTrigger(640041) || picoEvent->isTrigger(640051) )) return false; // 7.7GeV_2019
  if(mEnergy == 1 && vmsa::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(640001) || picoEvent->isTrigger(640011) || picoEvent->isTrigger(640021) || picoEvent->isTrigger(640031) || picoEvent->isTrigger(640041) || picoEvent->isTrigger(640051) )) return false; // 9.1GeV_2019
  if(mEnergy == 2 && vmsa::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(640001) || picoEvent->isTrigger(640011) || picoEvent->isTrigger(640021) || picoEvent->isTrigger(640031) || picoEvent->isTrigger(640041) || picoEvent->isTrigger(640051) )) return false; // 11.5GeV_2019
  if(mEnergy == 3 && vmsa::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(650001) || picoEvent->isTrigger(650011) || picoEvent->isTrigger(650021) || picoEvent->isTrigger(650031) || picoEvent->isTrigger(650041) || picoEvent->isTrigger(650051) )) return false; // 14.6GeV_2019
  if(mEnergy == 4 && vmsa::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(640001) || picoEvent->isTrigger(640011) || picoEvent->isTrigger(640021) || picoEvent->isTrigger(640031) || picoEvent->isTrigger(640041) || picoEvent->isTrigger(640051) )) return false; // 19.6GeV_2019
  return true;
}

bool StToFMatchCut::isPileUpEvent(int refMult, int numOfBTofMatch, double vz)
{
  if(mEnergy == 0 || mEnergy == 1 || mEnergy == 2 || mEnergy == 3)
  {
    if(numOfBTofMatch > vmsa::mMatchedToFMin) return kFALSE;
  }

  if(mEnergy == 4)
  {
    if(numOfBTofMatch <= vmsa::mMatchedToFMin) return kTRUE;

    int vzbin = -1;
    if(      vz >  -145 && vz < -87 ) vzbin = 0;
    else if( vz >= -87  && vz < -29 ) vzbin = 1;
    else if( vz >= -29  && vz <= 29 ) vzbin = 2;
    else if( vz >   29  && vz <= 87 ) vzbin = 3;
    else if( vz >   87  && vz <  145) vzbin = 4;

    double refLow  = vmsa::pl19[vzbin][0] + vmsa::pl19[vzbin][1]*numOfBTofMatch + vmsa::pl19[vzbin][2]*pow(numOfBTofMatch,2) + vmsa::pl19[vzbin][3]*pow(numOfBTofMatch,3) + vmsa::pl19[vzbin][4]*pow(numOfBTofMatch,4);
    double refHigh = vmsa::ph19[vzbin][0] + vmsa::ph19[vzbin][1]*numOfBTofMatch + vmsa::ph19[vzbin][2]*pow(numOfBTofMatch,2) + vmsa::ph19[vzbin][3]*pow(numOfBTofMatch,3) + vmsa::ph19[vzbin][4]*pow(numOfBTofMatch,4);

    if(refMult > refLow && refMult < refHigh) return kFALSE;
  }

  return kTRUE;
}

double StToFMatchCut::getRefMultReweight(double vz, int refMult)
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

int StToFMatchCut::getCentrality(double refMultCorr) // 9 Centrality bins 
{
  if(mEnergy == 4)
  {
    if(refMultCorr < 500 && refMultCorr >= 296) return 8; // 0-5%
    if(refMultCorr < 296 && refMultCorr >= 243) return 7; // 5-10%
    if(refMultCorr < 243 && refMultCorr >= 165) return 6; // 10-20%
    if(refMultCorr < 165 && refMultCorr >= 110) return 5; // 20-30%
    if(refMultCorr < 110 && refMultCorr >= 70 ) return 4; // 30-40%
    if(refMultCorr < 70  && refMultCorr >= 43 ) return 3; // 40-50%
    if(refMultCorr < 43  && refMultCorr >= 24 ) return 2; // 50-60%
    if(refMultCorr < 24  && refMultCorr >= 13 ) return 1; // 60-70%
    if(refMultCorr < 13  && refMultCorr >= 6  ) return 0; // 70-80%
    return -1;
  }
  return -1;
}

double StToFMatchCut::getEventWeight(int cent9, double refMult)
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


bool StToFMatchCut::passEventCut(StPicoDst *pico)
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
  //const Int_t runId = event->runId();
  //const Int_t refMult = event->refMult();
  const unsigned short numOfBTofMatch = event->nBTOFMatch();
  const Float_t vx = event->primaryVertex().x();
  const Float_t vy = event->primaryVertex().y();
  const Float_t vz = event->primaryVertex().z();
  //const Float_t zdcX = event->ZDCx();
  //const Float_t vzVpd = event->vzVpd();
  //const Bool_t isBES = (event->energy() < 200.);
  //mRefMultCorr->init(runId);

  // StRefMultCorr bad run cut
  //if(mRefMultCorr->isBadRun(runId))
  //{
  //  return kFALSE;
  //}

  // minBias event cut
  if(!isMinBias(event))
  {
    //cout << "Not minbias: SKIP!" << endl;
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
  //if(fabs(vz-vzVpd) > vmsa::mVzVpdDiffMax)
  //{
  //  return kFALSE;
  //}
  //}

  // refMult (0-80%) cut
  //if(!isBES) mRefMultCorr->initEvent(refMult,vz,zdcX); // 200GeV
  //if(isBES) 
  //mRefMultCorr->initEvent(refMult,vz); // BES
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

  //if(isPileUpEvent(refMult,numOfBTofMatch,vz))
  //{
  //  return kFALSE;
  //}
  if(numOfBTofMatch <= vmsa::mMatchedToFMin) return kFALSE;

  return kTRUE;
}

//---------------------------------------------------------------------------------

Int_t StToFMatchCut::getMatchedToF()
{
  return mMatchedToF;
}

Int_t StToFMatchCut::getNpirm()
{
  return mN_prim;
}

Int_t StToFMatchCut::getNnonprim()
{
  return mN_non_prim;
}
//---------------------------------------------------------------------------------

Float_t StToFMatchCut::getMass2(StPicoTrack *picoTrack, StPicoDst *picoDst)
{
  double mass2 = -999.9;
  // StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  int etofIndex = picoTrack->eTofPidTraitsIndex(); // return ToF PID traits
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
  if(etofIndex >= 0)
  {
    StPicoETofPidTraits *etofTrack = picoDst->etofPidTraits(etofIndex);
    const double beta = etofTrack->beta();
    // const double Momentum = picoTrack->pMom().mag(); // primary momentum for 54GeV_2017
    // const double Momentum = picoTrack->pMom().Mag(); // primary momentum for 27GeV_2018
    TVector3 primMom; // temp fix for StThreeVectorF & TVector3
    const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
    const double primPy    = picoTrack->pMom().y();
    const double primPz    = picoTrack->pMom().z();
    primMom.SetXYZ(primPx,primPy,primPz);
    const double Momentum = primMom.Mag(); // primary momentum

    if(etofTrack->matchFlag() > 0 && etofTrack->tof() != 0 && beta != 0)
    {
      mass2 = Momentum*Momentum*(1.0/(beta*beta) - 1.0);
    }
  }

  return mass2;
}

Float_t StToFMatchCut::getV0Mass2(StPicoTrack *picoTrack, StPicoDst *picoDst)
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

//---------------------------------------------------------------------------------

bool StToFMatchCut::passTrackBasic(StPicoTrack *track, StPicoEvent *picoEvent)
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

  // primary pt and momentum cut: PtMin = 0.1
  if(!(track->pMom().Perp() > vmsa::mGlobPtMin && track->pMom().Perp() < vmsa::mGlobPtMax))
  {
    return kFALSE;
  }

  const double vx = picoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const double vy = picoEvent->primaryVertex().y();
  const double vz = picoEvent->primaryVertex().z();

  if(track->gDCA(vx,vy,vz) > vmsa::mDcaTrMax_phi)
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StToFMatchCut::passSigPionCut(StPicoTrack* track, StPicoEvent* picoEvent, Float_t scale_nSigma_factor)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track,picoEvent)) return kFALSE;

  Float_t nSigmaPion = track->nSigmaPion();
  if(fabs(nSigmaPion*scale_nSigma_factor) > vmsa::mNSigmaToF)
  {
    return kFALSE;
  }
  return kTRUE;
}

bool StToFMatchCut::passSigProntonCut(StPicoTrack* track, StPicoEvent* picoEvent, Float_t scale_nSigma_factor)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track,picoEvent)) return kFALSE;

  Float_t nSigmaProton = track->nSigmaProton();
  if(fabs(nSigmaProton*scale_nSigma_factor) > vmsa::mNSigmaToF)
  {
    return kFALSE;
  }
  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StToFMatchCut::passSigKaonCut(StPicoTrack* track, StPicoEvent* picoEvent, Float_t scale_nSigma_factor)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track,picoEvent)) return kFALSE;

  Float_t nSigmaKaon = track->nSigmaKaon();
  if(fabs(nSigmaKaon*scale_nSigma_factor) > vmsa::mNSigmaToF)
  {
    return kFALSE;
  }
  return kTRUE;
}

bool StToFMatchCut::passToFMatchCut(StPicoTrack* picoTrack, StPicoEvent* picoEvent, StPicoDst *picoDst)
{
  if(!picoTrack) return kFALSE;

  if(!passTrackBasic(picoTrack,picoEvent)) return kFALSE;

  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  int etofIndex = picoTrack->eTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex < 0 && etofIndex < 0) return kFALSE;
  
  StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
  StPicoETofPidTraits *etofTrack = picoDst->etofPidTraits(etofIndex);
  const double betaBToF = tofTrack->btofBeta();
  const double betaEToF = etofTrack->beta();

  if(! (tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && betaBToF != 0) && ! (etofTrack->matchFlag() > 0 && etofTrack->tof() != 0 && betaEToF != 0)) return kFALSE; 

  return kTRUE;
}

//---------------------------------------------------------------------------------
