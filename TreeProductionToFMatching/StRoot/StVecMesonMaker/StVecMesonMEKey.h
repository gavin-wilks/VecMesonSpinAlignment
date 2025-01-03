#ifndef StVecMesonMEKey_h
#define StVecMesonMEKey_h
#include "StPhysicalHelixD.hh"
#include "Rtypes.h"
#include <map>
#include <vector>
#include "TLorentzVector.h"
#include "StThreeVectorF.hh"

class MEKey
{
  public:
    Int_t Centrality;
    Int_t vertex;
    Int_t EventPlane;
    Int_t Event_depth;
    Int_t charge;
    Int_t tag;

    MEKey() {}
    MEKey(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m, Int_t n) : Centrality(i), vertex(j), EventPlane(k), Event_depth(l), charge(m), tag(n) {}
    virtual ~MEKey() {/* */};

    bool operator<(const MEKey &r) const
    { // "r == right"
      if(Centrality == r.Centrality)
      {
        if(vertex == r.vertex)
	{
	  if(EventPlane == r.EventPlane)
	  {
	    if(Event_depth == r.Event_depth)
	    {
              if(charge == r.charge)
              {
	        return tag < r.tag;
              }
              else  
              {
                return charge < r.charge;
              }
	    }
 	    else
	    {
	      return Event_depth < r.Event_depth;
	    }
	  }
	  else
	  {
	    return EventPlane < r.EventPlane;
	  }
	}
	else
	{
	  return vertex < r.vertex;
	}
      }
      else
      {
        return Centrality < r.Centrality;
      }
    }

    void print()
    {
      cout << "Centrality = " << MEKey::Centrality << endl;
      cout << "vertex = " << MEKey::vertex << endl;
      cout << "EventPlane = " << MEKey::EventPlane << endl;
      cout << "Event_depth = " << MEKey::Event_depth << endl;
      cout << "charge = " << MEKey::charge << endl;
      cout << "tag = " << MEKey::tag << endl;
    }

    ClassDef(MEKey, 1)
};

typedef std::map<MEKey, std::vector<Bool_t> > vectorBoolMap;
typedef std::map<MEKey, std::vector<Int_t> > vectorIntMap;
typedef std::map<MEKey, std::vector<Float_t> > vectorFloatMap;
typedef std::map<MEKey, std::vector<StPhysicalHelixD> > vectorHelixMap;
typedef std::map<MEKey, std::vector<TLorentzVector> > vectorLorentzMap;
typedef std::map<MEKey, std::vector<StThreeVectorF> > vectorStThreeFMap;

#endif
