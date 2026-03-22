#ifndef __STMESONEVENT_H__
#define __STMESONEVENT_H__

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "StThreeVectorF.hh"
#include "TVector2.h"
#include "TLorentzVector.h"

// Adapted from StAlexPhiMesonEvent (A. Schmah 19.12.2011)
// G. Wilks 21.01.2022

class StMesonTrack : public TObject
{
  private:
    // Track properties
    // TrackA: for phi => K_plus  | for rho => pi_plus  | for KStar => K+/-
    // TrackB: for phi => K_minus | for rho => pi_minus | for KStar => pi+/-
    Float_t mMass2A; // mass2 of track A
    Float_t mMass2B;
    Float_t mNSigA; // nsigma dE/dx of particle A
    Float_t mNSigB;
    Float_t mDcaA; // distance of closest approach of particle A * charge
    Float_t mDcaB;
    Float_t mNHitsFitA; // distance of closest approach of particle A * charge
    Float_t mNHitsFitB;
    Float_t mNHitsMaxA; // distance of closest approach of particle A * charge
    Float_t mNHitsMaxB;
    Float_t mDEdxA; // distance of closest approach of particle A * charge
    Float_t mDEdxB;
    Int_t mChargeA; // charge, 0 => + | 1 => -
    Int_t mChargeB; 
    TLorentzVector mTrackA; // Lorentz Vector for track A (px, py, pz, mass2)
    TLorentzVector mTrackB;
    Int_t mFlagA; // Flag for Event A: 0 for Same Event, others for Mixed Event
    Int_t mFlagB;
    Bool_t mWhichTag;

  public:
    StMesonTrack() :
      mMass2A(-999.9),mMass2B(-999.9),mNSigA(-999.9),mNSigB(-999.9),mDcaA(-999.9),mDcaB(-999.9),mNHitsFitA(-999.9),mNHitsFitB(-999.9),mNHitsMaxA(-999.9),mNHitsMaxB(-999.9),mDEdxA(-999.9),mDEdxB(-999.9),mTrackA(0,0,0,0),mTrackB(0,0,0,0),mFlagA(-1),mFlagB(-1),mWhichTag(false)
  {
  }
    ~StMesonTrack() {}

    // setters
    void setMass2A(Float_t f)                  { mMass2A = f;       }
    void setMass2B(Float_t f)                  { mMass2B = f;       }
    void setNSigA(Float_t f)                   { mNSigA = f;    }
    void setNSigB(Float_t f)                   { mNSigB = f;    }
    void setDcaA(Float_t f)                    { mDcaA = f;         }
    void setDcaB(Float_t f)                    { mDcaB = f;         }
    void setNHitsFitA(Float_t f)               { mNHitsFitA = f;         }
    void setNHitsFitB(Float_t f)               { mNHitsFitB = f;         }
    void setNHitsMaxA(Float_t f)               { mNHitsMaxA = f;         }
    void setNHitsMaxB(Float_t f)               { mNHitsMaxB = f;         }
    void setDEdxA(Float_t f)                   { mDEdxA = f;         }
    void setDEdxB(Float_t f)                   { mDEdxB = f;         }
    void setChargeA(Int_t f)                   { mChargeA = f;      }
    void setChargeB(Int_t f)                   { mChargeB = f;      }
    void setTrackA(TLorentzVector f)           { mTrackA = f;       }
    void setTrackB(TLorentzVector f)           { mTrackB = f;       }
    void setFlagA(Int_t f)                     { mFlagA = f;        }
    void setFlagB(Int_t f)                     { mFlagB = f;        }
    void setWhichTagged(Bool_t f)              { mWhichTag = f;        }

    // getters
    Float_t getMass2A() const                  { return mMass2A;    }
    Float_t getMass2B() const                  { return mMass2B;    }
    Float_t getNSigA() const                   { return mNSigA; }
    Float_t getNSigB() const                   { return mNSigB; }
    Float_t getDcaA() const                    { return mDcaA;      }
    Float_t getDcaB() const                    { return mDcaB;      }
    Float_t getNHitsFitA() const               { return mNHitsFitA;      }
    Float_t getNHitsFitB() const               { return mNHitsFitB;      }
    Float_t getNHitsMaxA() const               { return mNHitsMaxA;      }
    Float_t getNHitsMaxB() const               { return mNHitsMaxB;      }
    Float_t getDEdxA() const                   { return mDEdxA;      }
    Float_t getDEdxB() const                   { return mDEdxB;      }
    Int_t getChargeA() const                   { return mChargeA;   }
    Int_t getChargeB() const                   { return mChargeB;   }
    TLorentzVector getTrackA() const           { return mTrackA;    }
    TLorentzVector getTrackB() const           { return mTrackB;    }
    Int_t getFlagA() const                     { return mFlagA;     }
    Int_t getFlagB() const                     { return mFlagB;     }
    Bool_t getWhichTagged() const              { return mWhichTag;     }

    ClassDef(StMesonTrack,1)  // A simple track of a particle
};

class StMesonEvent : public TObject
{
  private:
    StThreeVectorF mPrimaryvertex;
    Int_t   mRunId;
    Int_t   mEventId;
    Int_t   mRefMult;
    Int_t   mCentrality;
    Int_t   mN_prim;
    Int_t   mN_non_prim;
    Int_t   mN_Tof_match;

    Float_t mZDCx;
    Float_t mBBCx;
    Float_t mVzVpd;

    UShort_t      fNumTracks;

    TVector2 mQ2East;
    TVector2 mQ2West;
    TVector2 mQ2Full;
    double mPsiEast;
    double mPsiWest;
    double mPsiFull;
    Int_t   mNumTrackEast;
    Int_t   mNumTrackWest;
    Int_t   mNumTrackFull;
    Int_t   mNumTrackFullEast;
    Int_t   mNumTrackFullWest;

    TClonesArray* fTracks;      //->

  public:
    StMesonEvent() :
      mPrimaryvertex(-1.0,-1.0,-1.0),mRunId(-1),mEventId(-1),mRefMult(-1),mCentrality(-1),mN_prim(-1),mN_non_prim(-1),mN_Tof_match(-1),mZDCx(-1),mBBCx(-1),mVzVpd(-1),fNumTracks(0)
  {
    mQ2East.Set(-999.9,-999.9); // QVector2 East
    mQ2West.Set(-999.9,-999.9); // QVector2 West
    mQ2Full.Set(-999.9,-999.9); // QVector2 West
    mPsiEast = 0.0;
    mPsiWest = 0.0;
    mPsiFull = 0.0;
     
    mNumTrackEast = 0;
    mNumTrackWest = 0;
    mNumTrackFull = 0;
    mNumTrackFullEast = 0;
    mNumTrackFullWest = 0;

    fTracks      = new TClonesArray( "StMesonTrack", 10 );
  }
    ~StMesonEvent()
    {
      delete fTracks;
      fTracks = NULL;
    }

    void       setPrimaryVertex(StThreeVectorF r)      { mPrimaryvertex = r;       }
    StThreeVectorF    getPrimaryVertex() const         { return mPrimaryvertex;    }

    void       setRunId(Int_t  r)                      { mRunId = r;               }
    Int_t      getRunId() const                        { return mRunId;            }

    void       setEventId(Int_t  r)                    { mEventId = r;             }
    Int_t      getEventId() const                      { return mEventId;          }

    void       setRefMult(Int_t r)                     { mRefMult = r;             }
    Int_t      getRefMult() const                      { return mRefMult;          }

    void       setCentrality(Int_t r)                  { mCentrality = r;          }
    Int_t      getCentrality() const                   { return mCentrality;       }

    void       setN_prim(Int_t r)                      { mN_prim = r;              }
    Int_t      getN_prim() const                       { return mN_prim;           }

    void       setN_non_prim(Int_t r)                  { mN_non_prim = r;          }
    Int_t      getN_non_prim() const                   { return mN_non_prim;       }

    void       setN_Tof_match(Int_t r)                 { mN_Tof_match = r;         }
    Int_t      getN_Tof_match() const                  { return mN_Tof_match;      }


    void       setZDCx(Float_t r)                      { mZDCx = r;                }
    Float_t    getZDCx() const                         { return mZDCx;             }

    void       setBBCx(Float_t r)                      { mBBCx = r;                }
    Float_t    getBBCx() const                         { return mBBCx;             }

    void       setVzVpd(Float_t r)                     { mVzVpd = r;               }
    Float_t    getVzVpd() const                        { return mVzVpd;            }

    // ---------------------------------------QVector---------------------------------------------
    // QVector2 East
    void       setQ2East(TVector2 r)                   { mQ2East = r;              }
    TVector2   getQ2East() const                       { return mQ2East;           }
    // QVector2 West
    void       setQ2West(TVector2 r)                   { mQ2West = r;              }
    TVector2   getQ2West() const                       { return mQ2West;           }
    // QVector2 Full 
    void       setQ2Full(TVector2 r)                   { mQ2Full = r;              }
    TVector2   getQ2Full() const                       { return mQ2Full;           }
    // ---------------------------------------QVector---------------------------------------------
    void       setPsiEast(double r)                    { mPsiEast = r;            }
    double     getPsiEast() const                      { return mPsiEast;         }

    void       setPsiWest(double r)                    { mPsiWest = r;            }
    double     getPsiWest() const                      { return mPsiWest;         }

    void       setPsiFull(double r)                    { mPsiFull = r;            }
    double     getPsiFull() const                      { return mPsiFull;         }
    // ---------------------------------------QVector---------------------------------------------

    // -----------------------------------Number of Tracks----------------------------------------
    // East
    void       setNumTrackEast(Int_t r)                { mNumTrackEast = r;        }
    Int_t      getNumTrackEast() const                 { return mNumTrackEast;     }
    // West
    void       setNumTrackWest(Int_t r)                { mNumTrackWest = r;        }
    Int_t      getNumTrackWest() const                 { return mNumTrackWest;     }
    // Full 
    void       setNumTrackFull(Int_t r)                { mNumTrackFull = r;        }
    Int_t      getNumTrackFull() const                 { return mNumTrackFull;     }
    // Full East
    void       setNumTrackFullEast(Int_t r)            { mNumTrackFullEast = r;    }
    Int_t      getNumTrackFullEast() const             { return mNumTrackFullEast; }
    // Full West
    void       setNumTrackFullWest(Int_t r)            { mNumTrackFullWest = r;    }
    Int_t      getNumTrackFullWest() const             { return mNumTrackFullWest; }
    // -----------------------------------Number of Tracks----------------------------------------
    StMesonTrack* createTrack()
    {
      if (fNumTracks == fTracks->GetSize())
	fTracks->Expand( fNumTracks + 10 );
      if (fNumTracks >= 50000)
      {
	Fatal( "StMesonEvent::createTrack()", "ERROR: Too many tracks (>50000)!" );
	exit( 2 );
      }

      new((*fTracks)[fNumTracks++]) StMesonTrack;
      return (StMesonTrack*)((*fTracks)[fNumTracks - 1]);
    }
    void clearTrackList()
    {
      fNumTracks   = 0;
      fTracks      ->Clear();
    }
    UShort_t getNumTracks() const
    {
      return fNumTracks;
    }
    StMesonTrack* getTrack(UShort_t i) const
    {
      return i < fNumTracks ? (StMesonTrack*)((*fTracks)[i]) : NULL;
    }

    ClassDef(StMesonEvent,1)  // A simple event compiled of tracks
};


#endif // __STMESONEVENT_H__
