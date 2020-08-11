#include "StRoot/StRunIdMatching/StRunIdMatching.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StAlexPhiMesonEvent/StAlexPhiMesonEvent.h"
#include "StRoot/StVecMesonAna/StVecMesonCorr.h"
#include "StRoot/StRunIdEventsDb/StRunIdEventsDb.h"
#include "StThreeVectorF.hh"
#include "StMessMgr.h"
#include "TFile.h"
#include "TChain.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <fstream>

ClassImp(StRunIdMatching)

StRefMultCorr* StRunIdMatching::mRefMultCorr = NULL;
Int_t StRunIdMatching::mInPut_flag = 1;
char* StRunIdMatching::VM_EVENT_TREE = NULL;
char* StRunIdMatching::VM_EVENT_BRANCH = NULL;

//----------------------------------------------------
StRunIdMatching::StRunIdMatching(Int_t energy, Int_t X_flag, Int_t List, Long64_t start_event, Long64_t stop_event, Int_t mode)
{
  mEnergy = energy;
  mX_flag = X_flag;
  mList = List;
  mStart_Event = start_event;
  mStop_Event = stop_event;
  mMode = mode;
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  mVecMesonCorr = new StVecMesonCorr(mEnergy);
}

StRunIdMatching::~StRunIdMatching()
{
}
//----------------------------------------------------
// set Input/Output
void StRunIdMatching::setInputDir(const TString inputdir)
{
  mInputdir = inputdir.Copy();
  cout << "Input directory was set to: " << mInputdir.Data() << endl;
}
void StRunIdMatching::setOutputfile(const TString outputfile)
{
  mOutputfile = outputfile.Copy();
  cout << "Output file was set to: " << mOutputfile.Data() << endl;
}
void StRunIdMatching::setInPutList(const TString iInPutList)
{
  mInPutList = iInPutList.Copy();
  TString InFo_InPutList = Form("InPut %s list was set to: %s",vmsa::MixEvent[mX_flag].Data(),mInPutList.Data());
  cout << InFo_InPutList.Data() << endl;
}
void StRunIdMatching::setStopEvent(const Long64_t StopEvent)
{
    mStopEvent = StopEvent;
    cout << "nStopEvent = " << mStopEvent << endl;
}
void StRunIdMatching::setStartEvent(const Long64_t StartEvent)
{
    mStartEvent = StartEvent;
    cout << "nStartEvent = " << mStartEvent << endl;
}
//----------------------------------------------------
// initial functions
void StRunIdMatching::Init()
{
  mVecMesonCorr->InitReCenterCorrection();
  mVecMesonCorr->InitShiftCorrection();

  TString inputdir = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/Forest/",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str());
  setInputDir(inputdir);

  const Int_t list_start = vmsa::mList_Delta*mList + 1; // start list
  const Int_t list_stop  = vmsa::mList_Delta*(mList+1); // stop list

  TString InPutList = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/List/Split_%s_%s_%d_%d.list",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data(),vmsa::mBeamEnergy[mEnergy].c_str(),list_start,list_stop);
  setInPutList(InPutList);

  TString outputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/RunId/RunId_%s_%s_%d.txt",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data(),vmsa::mBeamEnergy[mEnergy].c_str(),mList);
  setOutputfile(outputfile);

  setStartEvent(Long64_t(mStart_Event));
  setStopEvent(Long64_t(mStop_Event));
  //----------------------------------------------------------------------------------------------------

  TString Notification = Form("Initializing parameters and input/output for %s %s",vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data());
  cout << Notification.Data() << endl;

  VM_EVENT_TREE       = (char*)vmsa::vm_tree[mMode].Data();
  VM_EVENT_BRANCH     = (char*)vmsa::vm_branch[mMode].Data();

  //----------------------------------------------------------------------------------------------------
  // input
  if (!mInPutList.IsNull())   // if input file is ok
  {
    TString InFo_List = Form("Open %s file list ",vmsa::MixEvent[mX_flag].Data());
    cout << InFo_List.Data() << endl;
    ifstream in(mInPutList);  // input stream
    if(in)
    {
      cout << "input file list is ok" << endl;
      mInPut = new TChain( VM_EVENT_TREE, VM_EVENT_TREE );
      char str[255];       // char array for each file name
      Long64_t entries_save = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  TString addfile;
	  addfile = str;
	  addfile = mInputdir+addfile;
	  mInPut->AddFile(addfile.Data(),-1, VM_EVENT_TREE );
	  Long64_t file_entries = mInPut->GetEntries();
	  cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      TString InFo_Warning = Form("WARNING: %s file input is problemtic",vmsa::MixEvent[mX_flag].Data());
      cout << InFo_Warning.Data() << endl;
      mInPut_flag = 0;
    }
  }

  // Set the input tree
  if (mInPut_flag == 1 && !mInPut->GetBranch( VM_EVENT_BRANCH ))
  {
    cerr << "ERROR: Could not find branch '"
      << VM_EVENT_BRANCH << "'in tree!" << endl;
  }

  if(mMode == 0) mPhiMeson_event = new StAlexPhiMesonEvent();

  if(mInPut_flag == 1)
  {
    if(mMode == 0) mInPut->SetBranchAddress( VM_EVENT_BRANCH, &mPhiMeson_event );

    Int_t num_events = mInPut->GetEntriesFast();
    cout << "Number of events in file(s) = " << num_events << endl;
    if(mStartEvent > num_events) mStartEvent = num_events;
    if(mStopEvent > num_events) mStopEvent   = num_events;

    cout << "New nStartEvent = " << mStartEvent << ", new nStopEvent = " << mStopEvent << endl;
  }

  InitMap();
}

void StRunIdMatching::InitMap()
{
  TString InPutFile = "/star/u/sunxuhit/AuAu200GeV/SpinAlignment/Phi/RunId/evtId.txt"; // embedded runId & evtId
  cout << "Input mapping was set to: " << InPutFile.Data() << endl;
  FILE *fp = fopen(InPutFile.Data(),"r");
  if(fp == NULL)
  {
    perror("Error opening mapping file");
  }
  else
  {
    Int_t runId, eventId;
    char line[80];

    Int_t line_counter = 0;
    while(fgets(line,80,fp))
    {
      sscanf(&line[0],"%d %d", &runId, &eventId);
      line_counter++;
      string KEY = Form("%d_%d",runId,eventId);
      mEmbdRunId[KEY] = runId;
      mEmbdEvtId[KEY] = eventId;
      // cout << "runId = " << runId << ", eventId = " << eventId << ", KEY = " << KEY.c_str() << ", mEmbdRunId = " << mEmbdRunId[KEY] << ", mEmbdEvtId = " << mEmbdEvtId[KEY] << endl;
    }
  }
}

void StRunIdMatching::Make()
{
  if(mMode == 0) MakePhi();
}


// loop phi meson Same Event
void StRunIdMatching::MakePhi()
{
  Long64_t start_event_use;
  Long64_t stop_event_use;

  start_event_use = mStartEvent;
  stop_event_use  = mStopEvent;
  mInPut->SetBranchAddress( VM_EVENT_BRANCH, &mPhiMeson_event);
  mInPut->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

  ofstream file_EventPlane;
  file_EventPlane.open(mOutputfile.Data());

  // Initialise Event Head
  StThreeVectorF PrimaryVertex(0.0,0.0,0.0);
  Int_t          RunId = 0;
  Int_t          EventId = 0;
  Int_t          RefMult = 0;
  Int_t          Centrality = 0;
  Int_t          N_prim = 0;
  Int_t          N_non_prim = 0;
  Int_t          N_Tof_match = 0;
  Float_t        ZDCx = 0.0; 
  Float_t        BBCx = 0.0; 
  Float_t        VzVpd = 0.0;
  Int_t          NumTrackUsed = 0;
  // ---------------------------------------QVector---------------------------------------------
  TVector2 Q2East(0.0,0.0);
  TVector2 Q2West(0.0,0.0);
  TVector2 Q2Full(0.0,0.0);
  // -----------------------------------Number of Tracks----------------------------------------
  Int_t   NumTrackEast = 0;
  Int_t   NumTrackWest = 0;
  Int_t   NumTrackFull = 0;
  Int_t   NumTrackFullEast = 0;
  Int_t   NumTrackFullWest = 0;

  for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
  {
    if (!mInPut->GetEntry( counter )) // take the event -> information is stored in event
      break;  // end of data chunk

    if (counter != 0  &&  counter % 1000 == 0)
      cout << "." << flush;
    if (counter != 0  &&  counter % 10000 == 0)
    {
      if((stop_event_use-start_event_use) > 0)
      {
	Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
	cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (VecMesonSpinAlignment) " << flush;
      }
    }

    // get Event Header
    PrimaryVertex    = mPhiMeson_event->getPrimaryVertex();
    RunId            = mPhiMeson_event->getRunId();
    EventId          = mPhiMeson_event->getEventId();
    RefMult          = mPhiMeson_event->getRefMult();
    Centrality       = mPhiMeson_event->getCentrality();
    N_prim           = mPhiMeson_event->getN_prim();
    N_non_prim       = mPhiMeson_event->getN_non_prim();
    N_Tof_match      = mPhiMeson_event->getN_Tof_match();
    ZDCx             = mPhiMeson_event->getZDCx(); 
    BBCx             = mPhiMeson_event->getBBCx(); 
    VzVpd            = mPhiMeson_event->getVzVpd();
    NumTrackUsed     = mPhiMeson_event->getNumTracks();
    Q2East           = mPhiMeson_event->getQ2East();
    Q2West           = mPhiMeson_event->getQ2West();
    Q2Full           = mPhiMeson_event->getQ2Full();
    NumTrackEast     = mPhiMeson_event->getNumTrackEast();
    NumTrackWest     = mPhiMeson_event->getNumTrackWest();
    NumTrackFull     = mPhiMeson_event->getNumTrackFull();
    NumTrackFullEast = mPhiMeson_event->getNumTrackFullEast();
    NumTrackFullWest = mPhiMeson_event->getNumTrackFullWest();

    // vz sign
    Int_t vz_sign;
    if(PrimaryVertex.z() > 0.0)
    {
      vz_sign = 0;
    }
    else
    {
      vz_sign = 1;
    }

    // Centrality
    mRefMultCorr->init(RunId);
    if(mEnergy == 6) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx); // 200 GeV
    if(mEnergy != 6) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z());
    const Int_t cent9 = Centrality;
    const Double_t reweight = mRefMultCorr->getWeight();

    // runIndex
    mRunIdEventsDb = StRunIdEventsDb::Instance(vmsa::mEnergyValue[mEnergy],vmsa::mBeamYear[mEnergy]);
    const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(RunId); // expensive
    // cout << runIndex << endl;

    // find embedding event
    string KEY = Form("%d_%d",RunId,EventId);
    if( !(mEmbdRunId[KEY] && mEmbdEvtId[KEY]) ) continue;

    // get Event Plane Info
    Float_t Psi2_west = mVecMesonCorr->calShiftAngle2West_EP(Q2West,runIndex,cent9,vz_sign);
    Float_t Psi2_east = mVecMesonCorr->calShiftAngle2East_EP(Q2East,runIndex,cent9,vz_sign);
    // cout << "runId = " << RunId << ", runIndex = " << runIndex << ", eventId = " << EventId << ", Psi2_east = " << Psi2_east << ", Psi2_west = " << Psi2_west << endl;
    file_EventPlane << RunId << "    " << EventId << "    " << Psi2_east << "    " << Psi2_west << endl;
  }

  cout << "." << flush;
  cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;

  file_EventPlane.close();
}

//-------------------------------------------------------------------
void StRunIdMatching::Finish()
{
  cout << "StRunIdMatching::Finish!!" << endl;
}
