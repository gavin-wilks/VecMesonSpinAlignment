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

  private:
    int mEnergy;
    std::map<int,int> map_runIndex;

    ClassDef(StUtility,1)
};

#endif
