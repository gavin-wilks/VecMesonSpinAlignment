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

    bool read_in_runIndex();
    int findRunIndex(int runId);
    
    bool read_in_badRunList();
    bool isBadRun(int runId);

  private:
    int mEnergy;
    std::map<int,int> map_runIndex;
    std::vector<int> vec_badRunId;

    ClassDef(StUtility,1)
};

#endif
