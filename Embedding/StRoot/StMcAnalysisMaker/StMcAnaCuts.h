#ifndef StMcAnaCuts_H
#define StMcAnaCuts_H

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */
#include <vector>

#include "Rtypes.h"
#include "StEvent/StEnumerations.h"

namespace McAnaCuts
{
  std::vector<unsigned int> getAllTriggers()
  {
    std::vector<unsigned int> t;
    t.push_back(640001); // miniBias trigger @ 19.6 GeV
    t.push_back(640011);  
    t.push_back(640021); 
    t.push_back(640031); 
    t.push_back(640041);
    t.push_back(640051);

    return t;
  }

  std::vector<unsigned int> const interesting_triggers = getAllTriggers();

  float const mcTrackStartVtxR = 1.0; // maximum
  int const geantId = 11; // K+
  // int const geantId = 12; // K-

  StDedxMethod dedxMethod = kLikelihoodFitId;

  int const maxNumberOfTriggers = 6;
}
#endif
