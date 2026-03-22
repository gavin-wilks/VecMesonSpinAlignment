#include "StRoot/StVecMesonAna/StVecMesonAna.h"
#include "StRoot/StVecMesonAna/StVecMesonCut.h"
#include "StRoot/StVecMesonAna/StVecMesonCorr.h"
#include "StRoot/StVecMesonAna/StVecMesonHistoManger.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StMesonEvent/StMesonEvent.h"
//#include "StRoot/StRunIdEventsDb/StRunIdEventsDb.h"
#include "StRoot/StVecMesonAna/StUtility.h"
#include "StarGenerator/UTIL/StarRandom.h"
#include "StThreeVectorF.hh"
#include "StMessMgr.h"
#include "TFile.h"
#include "TChain.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <fstream>
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TTree.h"

ClassImp(StVecMesonAna)

StRefMultCorr* StVecMesonAna::mRefMultCorr = NULL;
Int_t StVecMesonAna::mInPut_flag = 1;
char* StVecMesonAna::VM_EVENT_TREE_SE = NULL;
char* StVecMesonAna::VM_EVENT_TREE_ME = NULL;
char* StVecMesonAna::VM_EVENT_BRANCH = NULL;

//----------------------------------------------------
StVecMesonAna::StVecMesonAna(const Char_t *list, const Char_t *jobId, Int_t energy, Int_t X_flag, Int_t mode, Int_t etamode)
{
  mEnergy = energy;
  mEtaMode = etamode;
  if(mEtaMode == 0) mEtaCut = 1.0;
  if(mEtaMode == 3) mEtaCut = 0.4;
  if(mEtaMode == 4) mEtaCut = 0.6;
  if(mEtaMode == 5) mEtaCut = 0.8;
  mX_flag = X_flag;
  mList = list;
  mJobId = jobId;
  //mStart_Event = start_event;
  //mStop_Event = stop_event;
  mMode = mode;
  //if(!mRefMultCorr)
  //{
  //  mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  //}
  mVecMesonCorr = new StVecMesonCorr(mEnergy);
  mVecMesonCut = new StVecMesonCut(mEnergy);
  mVecMesonHistoManger = new StVecMesonHistoManger(2); // input 1 for 1st order EP, input 2 for 2nd order EP
}

StVecMesonAna::~StVecMesonAna()
{
}
//----------------------------------------------------
// set Input/Output
void StVecMesonAna::setInputDir(const TString inputdir)
{
  mInputdir = inputdir.Copy();
  cout << "Input directory was set to: " << mInputdir.Data() << endl;
}
void StVecMesonAna::setOutputfile(const TString outputfile)
{
  mOutputfile = outputfile.Copy();
  cout << "Output file was set to: " << mOutputfile.Data() << endl;
}
void StVecMesonAna::setInPutList(const TString iInPutList)
{
  mInPutList = iInPutList.Copy();
  TString InFo_InPutList = Form("InPut %s list was set to: %s",vmsa::MixEvent[mX_flag].Data(),mInPutList.Data());
  cout << InFo_InPutList.Data() << endl;
}
//----------------------------------------------------
// initial functions
void StVecMesonAna::Init()
{
  if(gRandom) delete gRandom;
  gRandom = new TRandom3();
  gRandom->SetSeed();

  //mVecMesonHistoManger->InitSys(0,mMode);
  //mVecMesonHistoManger->InitSys(1,mMode);

  mKaonTreeSE = new TTree("kaontreeSE","Tree for kaon info");
  mKaonTreeSE->Branch("cent", &mCentSE, "cent/I");
  mKaonTreeSE->Branch("weight", &mWeightSE, "weight/F");
  mKaonTreeSE->Branch("charge", &mChargeSE, "charge/I");
  mKaonTreeSE->Branch("pt", &mPtSE, "pt/F");
  mKaonTreeSE->Branch("rapidity", &mRapiditySE, "rapidity/F");
  mKaonTreeSE->Branch("eta", &mEtaSE, "eta/F");
  mKaonTreeSE->Branch("phi", &mPhiSE, "phi/F");
  mKaonTreeSE->Branch("hastofinfo",&mHasTofInfoSE,"hastofinfo/O");

  mKaonTreeME = new TTree("kaontreeME","Tree for kaon info");
  mKaonTreeME->Branch("cent", &mCentME, "cent/I");
  mKaonTreeME->Branch("weight", &mWeightME, "weight/F");
  mKaonTreeME->Branch("charge", &mChargeME, "charge/I");
  mKaonTreeME->Branch("pt", &mPtME, "pt/F");
  mKaonTreeME->Branch("rapidity", &mRapidityME, "rapidity/F");
  mKaonTreeME->Branch("eta", &mEtaME, "eta/F");
  mKaonTreeME->Branch("phi", &mPhiME, "phi/F");
  mKaonTreeME->Branch("hastofinfo",&mHasTofInfoME,"hastofinfo/O");

  cout << "Initialized the Branches" << endl;

  mUtility = new StUtility(mEnergy);
  mUtility->initRunIndex(); // initialize std::map for run index
  mRefMultCorr = new StRefMultCorr("refmult");
 

 // mVecMesonCorr->InitReCenterCorrection();
 // mVecMesonCorr->InitShiftCorrection();
 // mVecMesonCorr->InitResolutionCorr();
  cout << "mVecMesonHistoManger->InitSys(mX_flag,mMode);" << endl;
  //mVecMesonHistoManger->InitSys(mX_flag,mMode);
  //mVecMesonHistoManger->InitHistQA(mX_flag,mMode);
  //mVecMesonHistoManger->InitPhiSys(mX_flag,mMode);
  //cout << "mVecMesonHistoManger->InitSys(mX_flag,mMode);" << endl;
 // mVecMesonHistoManger->InitSys_EP(mX_flag,mMode);

  // TString inputdir = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Forest/",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str());
  //TString inputdir = Form("/gpfs01/star/pwg/gwilks3/AuAu%s/SpinAlignment/%s/Forest/",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str());
  //setInputDir(inputdir);

  //const Int_t list_start = vmsa::mList_Delta*mList + 1; // start list
  //const Int_t list_stop  = vmsa::mList_Delta*(mList+1); // stop list

  // TString InPutList = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/List/Split_%s_%s_%d_%d.list",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data(),vmsa::mBeamEnergy[mEnergy].c_str(),list_start,list_stop);
  //TString InPutList = Form("/gpfs01/star/pwg/gwilks3/AuAu%s/SpinAlignment/%s/List/Split_%s_%s_%d_%d.list",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data(),vmsa::mBeamEnergy[mEnergy].c_str(),list_start,list_stop);
  //setInPutList(InPutList);

  // TString outputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Yields/Yields_%s_%s_%d.root",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data(),vmsa::mBeamEnergy[mEnergy].c_str(),mList);
  TString outputfile = Form("Yields_%s_%s_%s_%s.root",vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data(),vmsa::mBeamEnergy[mEnergy].c_str(),mJobId);
  //setOutputfile(outputfile);

  //setStartEvent(Long64_t(mStart_Event));
  //setStopEvent(Long64_t(mStop_Event));
  //----------------------------------------------------------------------------------------------------

  TString Notification = Form("Initializing parameters and input/output for %s %s",vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data());
  cout << Notification.Data() << endl;
  mFile_OutPut = new TFile(outputfile.Data(),"RECREATE");

  VM_EVENT_TREE_SE    = (char*)Form("%s_SE",vmsa::vm_tree[mMode].Data());
  VM_EVENT_TREE_ME    = (char*)Form("%s_ME",vmsa::vm_tree[mMode].Data());
  VM_EVENT_BRANCH     = (char*)vmsa::vm_branch[mMode].Data();

  //----------------------------------------------------------------------------------------------------
  // input
  if (mList != NULL)   // if input file is ok
  {
    TString InFo_List = Form("Open %s file list ",vmsa::MixEvent[mX_flag].Data());
    cout << InFo_List.Data() << endl;
    ifstream in(mList);  // input stream
    if(in)
    {
      cout << "input file list is ok" << endl;
      mInPutSE = new TChain( VM_EVENT_TREE_SE, VM_EVENT_TREE_SE );
      char str[255];       // char array for each file name
      Long64_t entries_save = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  TString addfile;
	  addfile = str;
	  //addfile = mInputdir+addfile;
	  mInPutSE->AddFile(addfile.Data(),-1, VM_EVENT_TREE_SE );
	  Long64_t file_entries = mInPutSE->GetEntries();
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

  if (mList != NULL)   // if input file is ok
  {
    TString InFo_List = Form("Open %s file list ",vmsa::MixEvent[mX_flag].Data());
    cout << InFo_List.Data() << endl;
    ifstream in(mList);  // input stream
    if(in)
    {
      cout << "input file list is ok" << endl;
      mInPutME = new TChain( VM_EVENT_TREE_ME, VM_EVENT_TREE_ME );
      char str[255];       // char array for each file name
      Long64_t entries_save = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  TString addfile;
	  addfile = str;
	  //addfile = mInputdir+addfile;
	  mInPutME->AddFile(addfile.Data(),-1, VM_EVENT_TREE_ME );
	  Long64_t file_entries = mInPutME->GetEntries();
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
  if (mInPut_flag == 1 && !mInPutSE->GetBranch( VM_EVENT_BRANCH ))
  {
    cerr << "ERROR: Could not find branch '"
      << VM_EVENT_BRANCH << "'in tree!" << endl;
  }
  if (mInPut_flag == 1 && !mInPutME->GetBranch( VM_EVENT_BRANCH ))
  {
    cerr << "ERROR: Could not find branch '"
      << VM_EVENT_BRANCH << "'in tree!" << endl;
  }

  mMeson_eventSE = new StMesonEvent();
  mMeson_eventME = new StMesonEvent();

  if(mInPut_flag == 1)
  {
    mInPutSE->SetBranchAddress( VM_EVENT_BRANCH, &mMeson_eventSE );
  }
  if(mInPut_flag == 1)
  {
    mInPutME->SetBranchAddress( VM_EVENT_BRANCH, &mMeson_eventME );
  }
}

void StVecMesonAna::Make()
{
  if(mMode == 0) MakePhi();
}


// loop phi meson Same Event
void StVecMesonAna::MakePhi()
{
  //Long64_t start_event_use;
  //Long64_t stop_event_use;

  //start_event_use = mStartEvent;
  //stop_event_use  = mStopEvent;
  mInPutSE->SetBranchAddress( VM_EVENT_BRANCH, &mMeson_eventSE);
  mInPutME->SetBranchAddress( VM_EVENT_BRANCH, &mMeson_eventME);
  long numOfEventsSE = (long)mInPutSE->GetEntries();
  long numOfEventsME = (long)mInPutME->GetEntries();
  cout << "numOfEventSE" << numOfEventsSE << endl;
  mInPutSE->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

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

  //for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
  //{
  //  if (!mInPut->GetEntry( counter )) // take the event -> information is stored in event
  //   break;  // end of data chunk
  for(int i_event = 0; i_event < numOfEventsSE; ++i_event)
  {
    if(i_event%1000==0) cout << "processing events:  " << i_event << "/" << numOfEventsSE << endl;
    mInPutSE->GetEntry(i_event+1); 
    // get Event Header
    PrimaryVertex    = mMeson_eventSE->getPrimaryVertex();
    RunId            = mMeson_eventSE->getRunId();
    EventId          = mMeson_eventSE->getEventId();
    RefMult          = mMeson_eventSE->getRefMult();
    Centrality       = mMeson_eventSE->getCentrality();
    N_prim           = mMeson_eventSE->getN_prim();
    N_non_prim       = mMeson_eventSE->getN_non_prim();
    N_Tof_match      = mMeson_eventSE->getN_Tof_match();
    ZDCx             = mMeson_eventSE->getZDCx(); 
    BBCx             = mMeson_eventSE->getBBCx(); 
    VzVpd            = mMeson_eventSE->getVzVpd();
    NumTrackUsed     = mMeson_eventSE->getNumTracks();
    Q2East           = mMeson_eventSE->getQ2East();
    Q2West           = mMeson_eventSE->getQ2West();
    Q2Full           = mMeson_eventSE->getQ2Full();
    NumTrackEast     = mMeson_eventSE->getNumTrackEast();
    NumTrackWest     = mMeson_eventSE->getNumTrackWest();
    NumTrackFull     = mMeson_eventSE->getNumTrackFull();
    NumTrackFullEast = mMeson_eventSE->getNumTrackFullEast();
    NumTrackFullWest = mMeson_eventSE->getNumTrackFullWest();
   
    // Initialise Track 
    Float_t m2A = -999.9;
    Float_t m2B = -999.9;
    Float_t nsA = -999.9;
    Float_t nsB = -999.9;
    Float_t dcaA = -999.9;
    Float_t dcaB = -999.9;
    Float_t nhitsfitA = -999.9;
    Float_t nhitsfitB = -999.9;
    Float_t nhitsmaxA = -999.9;
    Float_t nhitsmaxB = -999.9;
    Float_t dEdxA = -999.9;
    Float_t dEdxB = -999.9;
    TLorentzVector lTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lTrackB(0.0,0.0,0.0,0.0);
    Int_t flagA = -1;
    Int_t flagB = -1;
    Bool_t tag = false;

    mRefMultCorr->init(RunId);

    if(mRefMultCorr->isBadRun( RunId ))
    {
      LOG_ERROR << "Bad Run! Skip!" << endm;
      continue;
    }

    //bool isPileUpEvent = false;
    // IMPORTANT: vertex position is needed for Au+Au 19.6 GeV 2019
    //if (mRefMultCorr->isPileUpEvent( RefMult, N_Tof_match, PrimaryVertex.z() ) ) isPileUpEvent = true;
    mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx);
 

    // vz sign 
    int vz_sign = 0;
    if(PrimaryVertex.z() > -70.0 && PrimaryVertex.z() <= -30.0) vz_sign = 0;
    if(PrimaryVertex.z() > -30.0 && PrimaryVertex.z() <= 0.0  ) vz_sign = 1;
    if(PrimaryVertex.z() > 0.0   && PrimaryVertex.z() <= +30.0) vz_sign = 2;
    if(PrimaryVertex.z() < +70.0 && PrimaryVertex.z() >  +30.0) vz_sign = 3;
    // Centrality
    const Int_t cent9 = Centrality;
    //const Double_t refMultCorr = mVecMesonCut->getRefMultReweight(PrimaryVertex.z(), RefMult);
    const Double_t reweight = mRefMultCorr->getWeight();

    const int runIndex = mUtility->findRunIndex(RunId); // find run index for a specific run
    //cout << "runIndex = " << runIndex << endl;
    //cout << "Centrality = " << Centrality << endl;
    // get Track Information
    if(mVecMesonCorr->passTrackEtaNumCut(NumTrackEast,NumTrackWest))
    {
      //PUT EVENT HISTOGRAMS HERE
      //mVecMesonHistoManger->FillEventHistQA(reweight,cent9,PrimaryVertex.x(),PrimaryVertex.y(),PrimaryVertex.z(),N_prim,RefMult);

      //cout << "After eta num cut" << endl;
      //cout << "NumTrackUsed = " << NumTrackUsed << endl;
      for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
      {
	mMeson_track = mMeson_eventSE->getTrack(nTracks);
	m2A = mMeson_track->getMass2A();
	m2B = mMeson_track->getMass2B();
	nsA = mMeson_track->getNSigA();
	nsB = mMeson_track->getNSigB();
	dcaA = mMeson_track->getDcaA();
	dcaB = mMeson_track->getDcaB();
	//nhitsfitA = mMeson_track->getNHitsFitA();
	//nhitsfitB = mMeson_track->getNHitsFitB();
	//nhitsmaxA = mMeson_track->getNHitsMaxA();
	//nhitsmaxB = mMeson_track->getNHitsMaxB();
	//dEdxA = mMeson_track->getDEdxA();
	//dEdxB = mMeson_track->getDEdxB();
	lTrackA = mMeson_track->getTrackA();
	lTrackB = mMeson_track->getTrackB();
	flagA = mMeson_track->getFlagA();
	flagB = mMeson_track->getFlagB();
        tag = mMeson_track->getWhichTagged();

	Float_t pA = lTrackA.P();
	Float_t pB = lTrackB.P();
	TLorentzVector lTrack = lTrackA + lTrackB; // phi-meson
	Float_t pt_lTrack = lTrack.Perp();
        Float_t eta_TrackA = lTrackA.PseudoRapidity();
        Float_t eta_TrackB = lTrackB.PseudoRapidity();
        Float_t phi_lTrack = lTrack.Phi();
        Float_t phi_lTrackA = lTrackA.Phi();
        Float_t phi_lTrackB = lTrackB.Phi();
        Float_t y_lTrackA = lTrackA.Rapidity();
        Float_t y_lTrackB = lTrackB.Rapidity();
        Float_t pt_lTrackA = lTrackA.Pt();
        Float_t pt_lTrackB = lTrackB.Pt();
	//Float_t eta_lTrack = lTrack.Rapidity();
        //cout << "m2A = " << m2A << "   m2B = " << m2B << endl;
	//if(
	//    ((m2A > 0.16 && m2A < 0.36) && (m2B > 0.16 && m2B < 0.36)) //||
        //  )
	{
          if( fabs(eta_TrackA) > mEtaCut || fabs(eta_TrackB) > mEtaCut ) continue;
          //cout << "pass mass cut" << endl;
	  // Float_t eta_lTrack = lTrack.Eta();
	  // if(TMath::Abs(eta_lTrack) > 1.0) continue;
	  Float_t rapidity_lTrack = lTrack.Rapidity();
	  Float_t eta_lTrack = lTrack.PseudoRapidity();
	  if(TMath::Abs(eta_TrackA) > vmsa::mEtaMax) continue;
	  if(TMath::Abs(eta_TrackB) > vmsa::mEtaMax) continue;
	  if(TMath::Abs(rapidity_lTrack) > 1.0) continue;
          //cout << "pass rapidity cut" << endl;
	  Float_t InvMass_lTrack = lTrack.M();
          TVector3 kp_vect = lTrackA.Vect();
          TVector3 km_vect = lTrackB.Vect();

	  for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++) // systematic loop for dca
	  {
	    if( !(mVecMesonCut->passTrackDcaSys(dcaA,dcaB,i_dca,mMode)) ) continue;

	    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1 /*vmsa::nSigKaon_stop*/; i_sig++) // systematic loop for nSigmaKaon
	    {
              if( i_dca != 0 && i_sig != 0) continue;
	      if( !(mVecMesonCut->passTrackSigSys(nsA,nsB,i_sig,mMode)) ) continue;
              if(!mVecMesonCut->passDipAngleCut(kp_vect,km_vect)) 
              {
                cout << "Failed dip angle cut" << endl;
                continue;
              }

              float scalingweight = 1.0;
              int ptbin = -1;
              int ybin = -1;
              for(Int_t i_pt_phi = 0; i_pt_phi < data_constants::rebinpttotal; i_pt_phi++)
              {
                if(pt_lTrack >= data_constants::rebinptval[i_pt_phi] && pt_lTrack < data_constants::rebinptval[i_pt_phi+1])
                {
                  for(Int_t i_y_phi = 0; i_y_phi < data_constants::rebinytotal; i_y_phi++)
                  {
                    if(rapidity_lTrack >= data_constants::rebinyval[i_y_phi] && rapidity_lTrack < data_constants::rebinyval[i_y_phi+1])
                    {    
                       ptbin = i_pt_phi;
                       ybin  = i_y_phi;
                       break;
                    }
                  }
                  break;
                }
              }

              if(ptbin < 0 || ybin < 0) continue;
     
              //cout << "ptbin = " << ptbin << ", ybin = " << ybin << endl;
             
              if(   InvMass_lTrack <  data_constants::InvMassLow[tag][ptbin][ybin]   
                 || InvMass_lTrack > data_constants::InvMassHigh[tag][ptbin][ybin]) continue;


              scalingweight = reweight * (1. - data_constants::fpoly[tag][ptbin][ybin]);

              mCentSE = cent9; 
              mWeightSE = scalingweight; 
              if(!tag) // tag = 0, K+ is tagged, so K- is probe
              {
                mChargeSE = 1;
                mPtSE = pt_lTrackB;
                mRapiditySE = y_lTrackB;
                mEtaSE = eta_TrackB;
                mPhiSE = phi_lTrackB;
                mHasTofInfoSE = (m2B > -900.0);
              }
              else // tag = 1, K- is tagged, so K+ is probe
              {
                mChargeSE = 0;
                mPtSE = pt_lTrackA;
                mRapiditySE = y_lTrackA;
                mEtaSE = eta_TrackA;
                mPhiSE = phi_lTrackA;
                mHasTofInfoSE = (m2A > -900.0);
              }

              mKaonTreeSE->Fill();
	    }
	  }
	}
      }
    }
  }

  for(int i_event = 0; i_event < numOfEventsME; ++i_event)
  {
    if(i_event%1000==0) cout << "processing events:  " << i_event << "/" << numOfEventsME << endl;
    mInPutME->GetEntry(i_event+1); 
    // get Event Header
    PrimaryVertex    = mMeson_eventME->getPrimaryVertex();
    RunId            = mMeson_eventME->getRunId();
    EventId          = mMeson_eventME->getEventId();
    RefMult          = mMeson_eventME->getRefMult();
    Centrality       = mMeson_eventME->getCentrality();
    N_prim           = mMeson_eventME->getN_prim();
    N_non_prim       = mMeson_eventME->getN_non_prim();
    N_Tof_match      = mMeson_eventME->getN_Tof_match();
    ZDCx             = mMeson_eventME->getZDCx(); 
    BBCx             = mMeson_eventME->getBBCx(); 
    VzVpd            = mMeson_eventME->getVzVpd();
    NumTrackUsed     = mMeson_eventME->getNumTracks();
    Q2East           = mMeson_eventME->getQ2East();
    Q2West           = mMeson_eventME->getQ2West();
    Q2Full           = mMeson_eventME->getQ2Full();
    NumTrackEast     = mMeson_eventME->getNumTrackEast();
    NumTrackWest     = mMeson_eventME->getNumTrackWest();
    NumTrackFull     = mMeson_eventME->getNumTrackFull();
    NumTrackFullEast = mMeson_eventME->getNumTrackFullEast();
    NumTrackFullWest = mMeson_eventME->getNumTrackFullWest();
   
    // Initialise Track 
    Float_t m2A = -999.9;
    Float_t m2B = -999.9;
    Float_t nsA = -999.9;
    Float_t nsB = -999.9;
    Float_t dcaA = -999.9;
    Float_t dcaB = -999.9;
    Float_t nhitsfitA = -999.9;
    Float_t nhitsfitB = -999.9;
    Float_t nhitsmaxA = -999.9;
    Float_t nhitsmaxB = -999.9;
    Float_t dEdxA = -999.9;
    Float_t dEdxB = -999.9;
    TLorentzVector lTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lTrackB(0.0,0.0,0.0,0.0);
    Int_t flagA = -1;
    Int_t flagB = -1;
    Bool_t tag = false;

    mRefMultCorr->init(RunId);

    if(mRefMultCorr->isBadRun( RunId ))
    {
      LOG_ERROR << "Bad Run! Skip!" << endm;
      continue;
    }

    //bool isPileUpEvent = false;
    // IMPORTANT: vertex position is needed for Au+Au 19.6 GeV 2019
    //if (mRefMultCorr->isPileUpEvent( RefMult, N_Tof_match, PrimaryVertex.z() ) ) isPileUpEvent = true;
    mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx);
 

    // vz sign 
    int vz_sign = 0;
    if(PrimaryVertex.z() > -70.0 && PrimaryVertex.z() <= -30.0) vz_sign = 0;
    if(PrimaryVertex.z() > -30.0 && PrimaryVertex.z() <= 0.0  ) vz_sign = 1;
    if(PrimaryVertex.z() > 0.0   && PrimaryVertex.z() <= +30.0) vz_sign = 2;
    if(PrimaryVertex.z() < +70.0 && PrimaryVertex.z() >  +30.0) vz_sign = 3;
    // Centrality
    const Int_t cent9 = Centrality;
    //const Double_t refMultCorr = mVecMesonCut->getRefMultReweight(PrimaryVertex.z(), RefMult);
    const Double_t reweight = mRefMultCorr->getWeight();

    const int runIndex = mUtility->findRunIndex(RunId); // find run index for a specific run
    //cout << "runIndex = " << runIndex << endl;
    //cout << "Centrality = " << Centrality << endl;
    // get Track Information
    if(mVecMesonCorr->passTrackEtaNumCut(NumTrackEast,NumTrackWest))
    {
      //PUT EVENT HISTOGRAMS HERE
      //mVecMesonHistoManger->FillEventHistQA(reweight,cent9,PrimaryVertex.x(),PrimaryVertex.y(),PrimaryVertex.z(),N_prim,RefMult);

      //cout << "After eta num cut" << endl;
      //cout << "NumTrackUsed = " << NumTrackUsed << endl;
      for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
      {
	mMeson_track = mMeson_eventME->getTrack(nTracks);
	m2A = mMeson_track->getMass2A();
	m2B = mMeson_track->getMass2B();
	nsA = mMeson_track->getNSigA();
	nsB = mMeson_track->getNSigB();
	dcaA = mMeson_track->getDcaA();
	dcaB = mMeson_track->getDcaB();
	//nhitsfitA = mMeson_track->getNHitsFitA();
	//nhitsfitB = mMeson_track->getNHitsFitB();
	//nhitsmaxA = mMeson_track->getNHitsMaxA();
	//nhitsmaxB = mMeson_track->getNHitsMaxB();
	//dEdxA = mMeson_track->getDEdxA();
	//dEdxB = mMeson_track->getDEdxB();
	lTrackA = mMeson_track->getTrackA();
	lTrackB = mMeson_track->getTrackB();
	flagA = mMeson_track->getFlagA();
	flagB = mMeson_track->getFlagB();
        tag = mMeson_track->getWhichTagged();

	Float_t pA = lTrackA.P();
	Float_t pB = lTrackB.P();
	TLorentzVector lTrack = lTrackA + lTrackB; // phi-meson
	Float_t pt_lTrack = lTrack.Perp();
        Float_t eta_TrackA = lTrackA.PseudoRapidity();
        Float_t eta_TrackB = lTrackB.PseudoRapidity();
        Float_t phi_lTrack = lTrack.Phi();
        Float_t phi_lTrackA = lTrackA.Phi();
        Float_t phi_lTrackB = lTrackB.Phi();
        Float_t y_lTrackA = lTrackA.Rapidity();
        Float_t y_lTrackB = lTrackB.Rapidity();
        Float_t pt_lTrackA = lTrackA.Pt();
        Float_t pt_lTrackB = lTrackB.Pt();
	//Float_t eta_lTrack = lTrack.Rapidity();
        //cout << "m2A = " << m2A << "   m2B = " << m2B << endl;
	//if(
	//    ((m2A > 0.16 && m2A < 0.36) && (m2B > 0.16 && m2B < 0.36)) //||
        //  )
	{
          if( fabs(eta_TrackA) > mEtaCut || fabs(eta_TrackB) > mEtaCut ) continue;
          //cout << "pass mass cut" << endl;
	  // Float_t eta_lTrack = lTrack.Eta();
	  // if(TMath::Abs(eta_lTrack) > 1.0) continue;
	  Float_t rapidity_lTrack = lTrack.Rapidity();
	  Float_t eta_lTrack = lTrack.PseudoRapidity();
	  if(TMath::Abs(eta_TrackA) > vmsa::mEtaMax) continue;
	  if(TMath::Abs(eta_TrackB) > vmsa::mEtaMax) continue;
	  if(TMath::Abs(rapidity_lTrack) > 1.0) continue;
          //cout << "pass rapidity cut" << endl;
	  Float_t InvMass_lTrack = lTrack.M();
          TVector3 kp_vect = lTrackA.Vect();
          TVector3 km_vect = lTrackB.Vect();

	  for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++) // systematic loop for dca
	  {
	    if( !(mVecMesonCut->passTrackDcaSys(dcaA,dcaB,i_dca,mMode)) ) continue;

	    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1 /*vmsa::nSigKaon_stop*/; i_sig++) // systematic loop for nSigmaKaon
	    {
              if( i_dca != 0 && i_sig != 0) continue;
	      if( !(mVecMesonCut->passTrackSigSys(nsA,nsB,i_sig,mMode)) ) continue;
              if(!mVecMesonCut->passDipAngleCut(kp_vect,km_vect)) 
              {
                cout << "Failed dip angle cut" << endl;
                continue;
              }

              float scalingweight = 1.0;
              int ptbin = -1;
              int ybin = -1;
              for(Int_t i_pt_phi = 0; i_pt_phi < data_constants::rebinpttotal; i_pt_phi++)
              {
                if(pt_lTrack >= data_constants::rebinptval[i_pt_phi] && pt_lTrack < data_constants::rebinptval[i_pt_phi+1])
                {
                  for(Int_t i_y_phi = 0; i_y_phi < data_constants::rebinytotal; i_y_phi++)
                  {
                    if(rapidity_lTrack >= data_constants::rebinyval[i_y_phi] && rapidity_lTrack < data_constants::rebinyval[i_y_phi+1])
                    {    
                       ptbin = i_pt_phi;
                       ybin  = i_y_phi;
                       break;
                    }
                  }
                  break;
                }
              }

              if(ptbin < 0 || ybin < 0) continue;
     
              //cout << "ptbin = " << ptbin << ", ybin = " << ybin << endl;
             
              if(   InvMass_lTrack <  data_constants::InvMassLow[tag][ptbin][ybin]   
                 || InvMass_lTrack > data_constants::InvMassHigh[tag][ptbin][ybin]) continue;


              scalingweight = reweight * (1. - data_constants::fpoly[tag][ptbin][ybin]) * data_constants::scalingME[tag][ptbin][ybin];

 
              mCentME = cent9; 
              mWeightME = scalingweight; 
              if(!tag) // tag = 0, K+ is tagged, so K- is probe
              {
                mChargeME = 1;
                mPtME = pt_lTrackB;
                mRapidityME = y_lTrackB;
                mEtaME = eta_TrackB;
                mPhiME = phi_lTrackB;
                mHasTofInfoME = (m2B > -900.0);
              }
              else // tag = 1, K- is tagged, so K+ is probe
              {
                mChargeME = 0;
                mPtME = pt_lTrackA;
                mRapidityME = y_lTrackA;
                mEtaME = eta_TrackA;
                mPhiME = phi_lTrackA;
                mHasTofInfoME = (m2A > -900.0);
              }

              mKaonTreeME->Fill();
	    }
	  }
	}
      }
    }
  }

  cout << endl;
}

//-------------------------------------------------------------------
void StVecMesonAna::Finish()
{
  mFile_OutPut->cd();
  mKaonTreeSE->Write();
  mKaonTreeME->Write();
  //mVecMesonHistoManger->WritePhiSys(mX_flag,mMode);
  //mVecMesonHistoManger->WriteSys(mX_flag,mMode);
  //mVecMesonHistoManger->WriteHistQA();
  //mVecMesonHistoManger->WriteSys_EP(mX_flag,mMode);
  mFile_OutPut->Close();
}
