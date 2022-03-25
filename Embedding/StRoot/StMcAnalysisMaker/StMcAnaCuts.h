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
  const double pl19[5][5] = {{-22.706,0.863809,-0.00036812,5.97123e-06,-1.27439e-08},   //-145 < vz < -87  cm
                       {-15.5924,0.604207,0.00131806,-2.04754e-06,5.73182e-10},   //-87  < vz < -29  cm 
                       {-13.1158,0.504708,0.00187998,-4.7317e-06,4.82342e-09},    //-29  < vz <  29  cm
                       {-15.9411,0.615061,0.00118242,-1.48902e-06,-2.29371e-10},  // 29  < vz <  87  cm 
                       {-22.468,0.839062,-0.000106213,4.93946e-06,-1.14503e-08}}; // 87  < vz <  145 cm

  const double ph19[5][5] = {{33.7733,1.75938,-0.00285868,8.50344e-06,-1.19215e-08},  //-145 < vz < -87  cm
                     {20.5875,1.67371,-0.00307534,7.93756e-06,-8.46257e-09},  //-87  < vz < -29  cm 
                     {15.1015,1.53929,-0.00269203,6.48876e-06,-6.06074e-09},  //-29  < vz <  29  cm
                     {20.7719,1.67316,-0.00315093,8.35824e-06,-9.14822e-09},  // 27  < vz <  87  cm 
                     {33.4926,1.79373,-0.00319461,9.56613e-06,-1.31049e-08}}; // 87  < vz <  145 cm


  const double vz_corr[29] = {
  0.998171,  // (-145,-135)cm
  0.99364,   // (-135,-125)cm
  0.991578,  // (-125,-115)cm
  0.990252,  // (-115,-105)cm
  0.990494,  // (-105,-95)cm
  0.990065,  // (-95,-85)cm
  0.990332,  // (-85,-75)cm
  0.996478,  // (-75,-65)cm
  0.999687,  // (-65,-55)cm
  0.998645,  // (-55,-45)cm
  0.993835, // (-45,-35)cm
  0.996273, // (-35,-25)cm
  0.998307, // (-25,-15)cm
  0.999295, // (-15,-5)cm
  1,        // (-5,5)cm
  1.00056,  // (5,15)cm
  1.00019,  // (15,25)cm
  0.999894, // (25,35)cm
  0.998907, // (35,45)cm
  1.00468,  // (45,55)cm
  1.0055,   // (55,65)cm
  1.00142,  // (65,75)cm
  0.996247, // (75,85)cm
  0.995789, // (85,95)cm
  0.996513, // (95,105)cm
  0.996928, // (105,115)cm
  0.998196, // (115,125)cm
  1.00097,  // (125,135)cm
  1.00699}; // (135,145)cm 
 
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

  double getRefMultReweight(double vz, int refMult, int mEnergy)
  {
    if(mEnergy == 4)
    { 
      if(vz > -5.0 && vz < 5.0) return refMult*vz_corr[14];

      for(int ivz = 0; ivz < 14; ivz++)
      {
        if(vz <  5.0+(ivz+1)*10.0 && vz >=  5.0+ivz*10.0) return refMult*vz_corr[15+ivz];
        if(vz > -5.0-(ivz+1)*10.0 && vz <= -5.0-ivz*10.0) return refMult*vz_corr[13-ivz];
      }
    }
    return refMult;
  }

  int getCentrality(double refMultCorr, int mEnergy) // 9 Centrality bins 
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
  
  double getEventWeight(int cent9, double refMult, int mEnergy)
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

  std::vector<unsigned int> const interesting_triggers = getAllTriggers();

  float const mcTrackStartVtxR = 1.0; // maximum
  int const geantId = 11; // K+
  // int const geantId = 12; // K-

  StDedxMethod dedxMethod = kLikelihoodFitId;

  int const maxNumberOfTriggers = 6;
}
#endif
