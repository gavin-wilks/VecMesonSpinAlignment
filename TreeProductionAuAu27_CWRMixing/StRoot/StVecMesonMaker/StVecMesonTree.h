#ifndef StVecMesonTree_h
#define StVecMesonTree_h

#include "StMessMgr.h"
#include "StPhysicalHelixD.hh"
#include "StThreeVectorF.hh"
#include "StVecMesonMEKey.h"
#include <vector>
#include "TVector2.h"
#include "../Utility/StSpinAlignmentCons.h"

class StPicoDst;
class StMesonEvent;
class StMesonTrack;
//class StV0TofCorrection;
//class StV0Event;
//class StV0Track;
class StVecMesonCut;
class TH1F;
class TH2F;
class TTree;
class TVector2;

class StVecMesonTree
{
  public:
    StVecMesonTree(Int_t energy);
    virtual ~StVecMesonTree();

    //phi-meson specific functions
    void InitPhi();

    void doPhi(Int_t,Int_t,Int_t,Int_t,Int_t&);
    void MixEvent_Phi(Int_t,StPicoDst*,Int_t,Float_t,Float_t);
    void FillTree_Phi();
    void clear_phi(Int_t,Int_t,Int_t);
    void size_phi(Int_t,Int_t,Int_t);
    ////////////////////////////////

    void WriteMass2Phi();

    void clearEvent();
    void passEvent(Int_t,Int_t,Int_t); // N_prim,N_non_prim,N_Tof_match
    void passEventPlane(TVector2,TVector2,TVector2); // qVector ater re-center: east, west, full
    void passNumTrack(Int_t,Int_t,Int_t,Int_t,Int_t); // Number of Tracks: east, west, full, full_east, full_west

  private:
    StVecMesonCut *mVecMesonCut;
    //StV0TofCorrection *mTofCorr;
    TH2F *h_Mass2;
    TH2F *h_KdEdxRig;
    TH2F *h_KM2Rig;
    TH2F *h_KInvBetaRig;
    TH2F *h_PidEdxRig;
    TH2F *h_PiM2Rig;
    TH2F *h_PiInvBetaRig;
    Int_t mEventCounter2[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi]; // 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin

    // 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin, 3 = mixed event bin, 4 = charge bin(0 for pos, 1 for neg) || push_back->track
    vectorHelixMap mHelix;
    vectorFloatMap mMomentum;
    vectorFloatMap mMass2;
    vectorFloatMap mDca;
    vectorFloatMap mNHitsFit;
    vectorFloatMap mNSigma;
    vectorIntMap mCharge;
    vectorIntMap mEventBin;

    TTree *mTree;
    StMesonEvent *mMesonEvent;
    StMesonEvent *mMesonEvents[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi][2];
    StMesonTrack *mMesonTrack;

    // event information | 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin || push_back->event
    Int_t mMixed[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi][2];
    std::vector<StThreeVectorF> mPrimaryvertex[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Int_t> mRefMult[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Int_t> mCentrality[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Int_t> mRunId[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Int_t> mEventId[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Int_t> mN_prim[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Int_t> mN_non_prim[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Int_t> mN_Tof_match[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Float_t> mZDCx[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Float_t> mBBCx[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Float_t> mVzVpd[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Float_t> mField[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<UShort_t> mNumTracks[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<TVector2> mQ2East[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<TVector2> mQ2West[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<TVector2> mQ2Full[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Int_t> mNumTrackEast[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Int_t> mNumTrackWest[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Int_t> mNumTrackFull[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Int_t> mNumTrackFullEast[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];
    std::vector<Int_t> mNumTrackFullWest[vmsa::Bin_Centrality][vmsa::Bin_VertexZ][vmsa::Bin_Phi_Psi];

    // passing variable
    Int_t mNumber_prim, mNumber_non_prim, mNumber_Tof_match;
    TVector2 mQVector2East, mQVector2West, mQVector2Full;
    Int_t mTrackEtaEast, mTrackEtaWest, mTrackFull, mTrackFullEast, mTrackFullWest;
    Int_t mEnergy;

  ClassDef(StVecMesonTree,1)
};
#endif
