#ifndef StUtility_h
#define StUtility_h

#include "StMessMgr.h"
#include <map>

class StUtility
{
  public:
    StUtility(int energy);
    virtual ~StUtility();

    void initRunIndex();
    void initEventPlane();

    bool read_in_runIndex();
    int findRunIndex(int runId);

    bool read_in_EventPlane();
    float findEPwest(int eventId);
    float findEPeast(int eventId);
    float findEPfull(int eventId);
    
    bool read_in_badRunList();
    bool isBadRun(int runId);

  private:
    int mEnergy;
    std::map<int,int> map_runIndex;
    std::map<int,float> map_ep_west;
    std::map<int,float> map_ep_east;
    std::map<int,float> map_ep_full;
    std::vector<int> vec_badRunId;

    ClassDef(StUtility,1)
};

#endif
