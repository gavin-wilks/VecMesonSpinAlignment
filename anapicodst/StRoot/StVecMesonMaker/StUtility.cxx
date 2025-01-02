#include "StRoot/StVecMesonMaker/StUtility.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TMath.h"

ClassImp(StUtility)

//---------------------------------------------------------------------------------

StUtility::StUtility(int energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StUtility::~StUtility()
{
  /* */
}

//---------------------------------------------------------------------------------


void StUtility::initRunIndex()
{ 
  bool isOpen_runIndex = read_in_runIndex();
  if(isOpen_runIndex) std::cout << "Run Index read in!" << std::endl;

  bool isOpen_badRunList = read_in_badRunList(); // read in Bad Run List
  if(isOpen_badRunList) std::cout << "Bad Run List read in!" << std::endl;
}

void StUtility::initEventPlane()
{ 
  bool isOpen_ep = read_in_EventPlane();
  if(isOpen_ep) std::cout << "Event Planes read in!" << std::endl;
}

int StUtility::findRunIndex(int runId)
{
  // print map_runIndex content:
  /*
     for (std::map<int,int>::iterator it=map_runIndex.begin(); it!=map_runIndex.end(); ++it)
     {
     std::cout << it->first << " => " << it->second << '\n';
     }
     */

  std::map<int,int>::iterator it_runId = map_runIndex.find(runId);
  if(it_runId == map_runIndex.end())
  {
    // std::cout << "StUtility -> could not find in full run list! & send signal to kill the run!" << std::endl;
    return -999;
  }
  else
  {
    // std::cout << "StUtility -> runId: " << it_runId->first << " => runIndex: " << it_runId->second << std::endl;
    return it_runId->second;
  }

  return -999;
}

bool StUtility::read_in_runIndex()
{
  map_runIndex.clear();
  // std::string inputfile = Form("StRoot/StEventPlaneMaker/runIndex_%s.txt",vmsa::mBeamEnergy[mEnergy].c_str());
  std::string inputfile = Form("StRoot/Utility/RunIndex/runIndex_%s_%d.txt",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mBeamYear[mEnergy]);
  std::cout << "inputfile = " << inputfile.c_str() << std::endl;
  std::ifstream file_runIndex ( inputfile.c_str() );
  if ( !file_runIndex.is_open() )
  {
    std::cout << "Abort. Fail to read in run Index file: " << inputfile << std::endl;
    return false;
  }

  int temp_runId = 0, temp_runIndex = 0;
  std::cout << "reading run Index: " << std::endl;
  while (file_runIndex >> temp_runId >> temp_runIndex)
  {
    // std::cout << "runId = " << temp_runId << ", runIndex = " << temp_runIndex << std::endl;
    map_runIndex[temp_runId] = temp_runIndex;
  }
  file_runIndex.close();

  return true;
}

bool StUtility::read_in_badRunList()
{
  vec_badRunId.clear();

  // std::string inputfile = Form("StRoot/StRunQAMaker/runIndex_%s.txt",runQA::mBeamEnergy[mEnergy].c_str());
  std::string inputfile = Form("StRoot/Utility/RunIndex/badRunList_%s_%d.txt",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mBeamYear[mEnergy]);
  std::cout << "inputfile = " << inputfile.c_str() << std::endl;
  std::ifstream file_badRunList ( inputfile.c_str() );
  if ( !file_badRunList.is_open() )
  {
    std::cout << "Abort. Fail to read in bad Run List file: " << inputfile << std::endl;
    return false;
  }

  int runId = 0;
  std::cout << "reading bad Bun List: " << inputfile.c_str() << std::endl;
  while (file_badRunList >> runId)
  {
    vec_badRunId.push_back(runId);
  }
  file_badRunList.close();

  // print vec_badRunId content:
  // for (std::vector<int>::iterator it=vec_badRunId.begin(); it!=vec_badRunId.end(); ++it)
  // {
  //   std::cout << (*it) << std::endl;;
  // }

  return true;
}

bool StUtility::isBadRun(int runId)
{
  // print vec_badRunId content:
  // for (std::vector<int>::iterator it=vec_badRunId.begin(); it!=vec_badRunId.end(); ++it)
  // {
  //   std::cout << (*it) << std::endl;;
  // }

  std::vector<int>::iterator it_runId = std::find(vec_badRunId.begin(), vec_badRunId.end(), runId);

  return ( it_runId != vec_badRunId.end() ); // true if can be found in bad run list
}

bool StUtility::read_in_EventPlane()
{
  map_ep_west.clear();
  map_ep_east.clear();
  map_ep_full.clear();
  // std::string inputfile = Form("StRoot/StEventPlaneMaker/runIndex_%s.txt",vmsa::mBeamEnergy[mEnergy].c_str());
  std::string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Phi/EventPlane/EP_%s_%d.txt",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mBeamYear[mEnergy]);
  std::cout << "inputfile = " << inputfile.c_str() << std::endl;
  std::ifstream file_ep ( inputfile.c_str() );
  if ( !file_ep.is_open() )
  {
    std::cout << "Abort. Fail to read in run Index file: " << inputfile << std::endl;
    return false;
  }

  int temp_runIndex = 0, temp_eventId = 0;
  float ep_west = 0.0, ep_east = 0.0, ep_full = 0.0;
  std::cout << "reading event planes: " << std::endl;
  while (file_ep >> temp_runIndex >> temp_eventId >> ep_west >> ep_east >> ep_full)
  {
    // std::cout << "runId = " << temp_runId << ", runIndex = " << temp_runIndex << std::endl;
    map_ep_west[temp_eventId] = ep_west;
    map_ep_east[temp_eventId] = ep_east;
    map_ep_full[temp_eventId] = ep_full;
  }
  file_ep.close();

  return true;
}

float StUtility::findEPwest(int eventId)
{
  std::map<int,float>::iterator it_eventId = map_ep_west.find(eventId);
  if(it_eventId == map_ep_west.end())
  {
    //std::cout << "StUtility -> could not find in full event list! & send signal to kill the run!" << std::endl;
    return -999;
  }
  else
  {
    //std::cout << "StUtility -> eventId: " << it_eventId->first << " => ep_west: " << it_eventId->second << std::endl;
    return it_eventId->second;
  }

  return -999;
}

float StUtility::findEPeast(int eventId)
{
  std::map<int,float>::iterator it_eventId = map_ep_east.find(eventId);
  if(it_eventId == map_ep_east.end())
  {
    //std::cout << "StUtility -> could not find in full event list! & send signal to kill the run!" << std::endl;
    return -999;
  }
  else
  {
    //std::cout << "StUtility -> eventId: " << it_eventId->first << " => ep_east: " << it_eventId->second << std::endl;
    return it_eventId->second;
  }

  return -999;
}

float StUtility::findEPfull(int eventId)
{
  std::map<int,float>::iterator it_eventId = map_ep_full.find(eventId);
  if(it_eventId == map_ep_full.end())
  {
    //std::cout << "StUtility -> could not find in full event list! & send signal to kill the run!" << std::endl;
    return -999;
  }
  else
  {
    //std::cout << "StUtility -> eventId: " << it_eventId->first << " => ep_full: " << it_eventId->second << std::endl;
    return it_eventId->second;
  }

  return -999;
}

