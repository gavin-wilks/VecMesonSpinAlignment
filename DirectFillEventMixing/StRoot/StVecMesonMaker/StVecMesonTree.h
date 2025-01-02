#ifndef StVecMesonTree_h
#define StVecMesonTree_h

#include "StMessMgr.h"
#include "StPhysicalHelixD.hh"
#include "StThreeVectorF.hh"
#include "StVecMesonMEKey.h"
#include <vector>
#include <map>
#include "TVector2.h"
#include "TLorentzVector.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"

class StPicoDst;
class StMesonEvent;
class StMesonTrack;
class StRefMultCorr;
//class StV0TofCorrection;
//class StV0Event;
//class StV0Track;
class StVecMesonCut;
class StVecMesonFillCorr;
class StUtility;
class TH1F;
class TH2F;
class TTree;
class TVector2;
class StVecMesonHist;
class StVecMesonHistoFlow;
class TLorentzVector;
class MEKey;

typedef std::map<MEKey, std::vector<Int_t> > vectorIntMap;
typedef std::map<MEKey, std::vector<Float_t> > vectorFloatMap;
typedef std::map<MEKey, std::vector<StPhysicalHelixD> > vectorHelixMap;
typedef std::map<MEKey, std::vector<TLorentzVector> > vectorLorentzMap;
typedef std::map<MEKey, std::vector<StThreeVectorF> > vectorStThreeFMap;

class StVecMesonTree
{
  public:
    StVecMesonTree(Int_t energy, Int_t flag);
    virtual ~StVecMesonTree();

    //KStar-meson specific functions
    void InitPhi();
    void doPhi(Int_t,Int_t,Int_t,Int_t);
    void MixEvent_Phi(Int_t,StPicoDst*,Int_t,Float_t,Float_t,Int_t,Float_t,Int_t);
    void FillPhi(StMesonTrack*,Int_t,Int_t,Int_t);
    void clear_phi(Int_t,Int_t,Int_t,Int_t);
    void WriteMass2Phi();
    /////////////////////////////////
    //KStar-meson specific functions
    void InitKStar();
    void doKStar(Int_t,Int_t,Int_t,Int_t);
    void MixEvent_KStar(Int_t,StPicoDst*,Int_t,Float_t,Float_t,Int_t,Float_t,Int_t);
    void FillKStar(StMesonTrack*,Int_t,Int_t,Int_t);
    void clear_kstar(Int_t,Int_t,Int_t,Int_t);
    void WriteMass2KStar();
    /////////////////////////////////

    void clearEvent();
    void passEvent(Int_t,Int_t,Int_t); // N_prim,N_non_prim,N_Tof_match
    void passEventPlane(TVector2,TVector2,TVector2); // qVector ater re-center: east, west, full
    void passNumTrack(Int_t,Int_t,Int_t,Int_t,Int_t); // Number of Tracks: east, west, full, full_east, full_west
    //void passExtraEvent(Int_t,Int_t,Int_t,Double_t);    

  private:
    StUtility *mUtility;
    StVecMesonFillCorr *mVecMesonCorr;
    StVecMesonCut *mVecMesonCut;
    StVecMesonHist *mVecMesonHistoManger;
    StVecMesonHistoFlow *mVecMesonHistoFlow;
    StRefMultCorr *mRefMultCorr;
    //StV0TofCorrection *mTofCorr;
    TH2F *h_Mass2;
    TH2F *h_KdEdxRig;
    TH2F *h_KM2Rig;
    TH2F *h_KInvBetaRig;
    TH2F *h_PidEdxRig;
    TH2F *h_PiM2Rig;
    TH2F *h_PiInvBetaRig;

    // 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin, 3 = mixed event bin, 4 = charge bin(0 for pos, 1 for neg) || push_back->track
    TTree *mTree;
    StMesonEvent *mMesonEvent;
    StMesonEvent *mMesonEventRK;
    StMesonEvent *mMesonEventRP;
    StMesonTrack *mMesonTrack;

    vectorFloatMap mMass2;
    vectorFloatMap mDca;
    vectorFloatMap mNHitsFit;
    vectorFloatMap mNSigma;
    vectorIntMap mCharge;
    vectorHelixMap mHelix;
    vectorFloatMap mMomentum;


    Int_t mEventCounter2[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi]; // 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin
    // event information | 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin || push_back->event
    std::vector<StThreeVectorF> mPrimaryvertex[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Int_t> mRefMult[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Int_t> mCentrality[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    std::vector<Int_t> mRunIdx[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    std::vector<Float_t> mReweight[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    std::vector<Int_t> mVzEPBin[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Int_t> mEventId[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Int_t> mN_prim[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Int_t> mN_non_prim[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Int_t> mN_Tof_match[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Float_t> mZDCx[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Float_t> mBBCx[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Float_t> mVzVpd[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Float_t> mField[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<UShort_t> mNumTracks[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    std::vector<TVector2> mQ2East[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    std::vector<TVector2> mQ2West[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<TVector2> mQ2Full[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    std::vector<Int_t> mNumTrackEast[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    std::vector<Int_t> mNumTrackWest[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Int_t> mNumTrackFull[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Int_t> mNumTrackFullEast[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];
    //std::vector<Int_t> mNumTrackFullWest[9][vmsa::Bin_VertexZ][vmsa::Bin_Psi];

    // passing variable
    Int_t mNumber_prim, mNumber_non_prim, mNumber_Tof_match;
    TVector2 mQVector2East, mQVector2West, mQVector2Full;
    Int_t mTrackEtaEast, mTrackEtaWest, mTrackFull, mTrackFullEast, mTrackFullWest;
    //Int_t mRunIndex, mCent, mVzSign;
    //Double_t mReweight;
    Int_t mEnergy;
    Int_t mX_flag;

  ClassDef(StVecMesonTree,1)
};
#endif
