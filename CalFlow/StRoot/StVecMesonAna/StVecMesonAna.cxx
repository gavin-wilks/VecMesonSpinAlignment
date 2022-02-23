#include "StRoot/StVecMesonAna/StVecMesonAna.h"
#include "StRoot/StVecMesonAna/StVecMesonCut.h"
#include "StRoot/StVecMesonAna/StVecMesonCorr.h"
#include "StRoot/StVecMesonAna/StVecMesonHistoManger.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StMesonEvent/StMesonEvent.h"
#include "StRoot/StVecMesonAna/StUtility.h"
#include "StThreeVectorF.hh"
#include "StMessMgr.h"
#include "TFile.h"
#include "TChain.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <fstream>

ClassImp(StVecMesonAna)

StRefMultCorr* StVecMesonAna::mRefMultCorr = NULL;
Int_t StVecMesonAna::mInPut_flag = 1;
char* StVecMesonAna::VM_EVENT_TREE = NULL;
char* StVecMesonAna::VM_EVENT_BRANCH = NULL;

//----------------------------------------------------
StVecMesonAna::StVecMesonAna(const Char_t *list, const Char_t *jobId, Int_t energy, Int_t X_flag, Int_t mode)
{
  mEnergy = energy;
  mX_flag = X_flag;
  mList = list;
  mJobId = jobId;
  mMode = mode;
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  mVecMesonCorr = new StVecMesonCorr(mEnergy);
  mVecMesonCut = new StVecMesonCut(mEnergy);
  mVecMesonHistoManger = new StVecMesonHistoManger();
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
  mUtility = new StUtility(mEnergy);
  mUtility->initRunIndex(); // initialize std::map for run index
 
  mVecMesonCorr->InitReCenterCorrection();
  mVecMesonCorr->InitShiftCorrection();
  mVecMesonCorr->InitResolutionCorr();
  mVecMesonHistoManger->InitSys(mX_flag,mMode);
  mVecMesonHistoManger->InitSys_EP(mX_flag,mMode);

  TString outputfile = Form("Yields_%s_%s_%s_%s.root",vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data(),vmsa::mBeamEnergy[mEnergy].c_str(),mJobId);

  TString Notification = Form("Initializing parameters and input/output for %s %s",vmsa::mPID[mMode].c_str(),vmsa::MixEvent[mX_flag].Data());
  cout << Notification.Data() << endl;
  mFile_OutPut = new TFile(outputfile.Data(),"RECREATE");

  VM_EVENT_TREE       = (char*)vmsa::vm_tree[mMode].Data();
  VM_EVENT_BRANCH     = (char*)vmsa::vm_branch[mMode].Data();

  // input
  if (mList != NULL)   // if input file is ok
  {
    TString InFo_List = Form("Open %s file list ",vmsa::MixEvent[mX_flag].Data());
    cout << InFo_List.Data() << endl;
    ifstream in(mList);  // input stream
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

  if(mMode == 0) mMeson_event = new StMesonEvent();

  if(mInPut_flag == 1)
  {
    if(mMode == 0) mInPut->SetBranchAddress( VM_EVENT_BRANCH, &mMeson_event );
  }
}

void StVecMesonAna::Make()
{
  if(mMode == 0) MakePhi();
  if(mMode == 1) MakeRho();
  if(mMode == 2) MakeKStar();
}


// loop phi meson Same Event
void StVecMesonAna::MakePhi()
{
  mInPut->SetBranchAddress( VM_EVENT_BRANCH, &mMeson_event);
  long numOfEvents = (long)mInPut->GetEntries();
  cout << "numOfEvent" << numOfEvents << endl;
  mInPut->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

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

  for(int i_event = 0; i_event < numOfEvents; ++i_event)
  {
    if(i_event%100==0) cout << "processing events:  " << i_event << "/" << numOfEvents << endl;
    mInPut->GetEntry(i_event+1); 
    // get Event Header
    PrimaryVertex    = mMeson_event->getPrimaryVertex();
    RunId            = mMeson_event->getRunId();
    EventId          = mMeson_event->getEventId();
    RefMult          = mMeson_event->getRefMult();
    Centrality       = mMeson_event->getCentrality();
    N_prim           = mMeson_event->getN_prim();
    N_non_prim       = mMeson_event->getN_non_prim();
    N_Tof_match      = mMeson_event->getN_Tof_match();
    ZDCx             = mMeson_event->getZDCx(); 
    BBCx             = mMeson_event->getBBCx(); 
    VzVpd            = mMeson_event->getVzVpd();
    NumTrackUsed     = mMeson_event->getNumTracks();
    Q2East           = mMeson_event->getQ2East();
    Q2West           = mMeson_event->getQ2West();
    Q2Full           = mMeson_event->getQ2Full();
    NumTrackEast     = mMeson_event->getNumTrackEast();
    NumTrackWest     = mMeson_event->getNumTrackWest();
    NumTrackFull     = mMeson_event->getNumTrackFull();
    NumTrackFullEast = mMeson_event->getNumTrackFullEast();
    NumTrackFullWest = mMeson_event->getNumTrackFullWest();
   
    // Initialise Track 
    Float_t m2A = -999.9;
    Float_t m2B = -999.9;
    Float_t nsA = -999.9;
    Float_t nsB = -999.9;
    Float_t dcaA = -999.9;
    Float_t dcaB = -999.9;
    TLorentzVector lTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lTrackB(0.0,0.0,0.0,0.0);
    Int_t flagA = -1;
    Int_t flagB = -1;

    // vz sign 
    int vz_sign = 0;
    if(PrimaryVertex.z() > -70.0 && PrimaryVertex.z() <= -30.0) vz_sign = 0;
    if(PrimaryVertex.z() > -30.0 && PrimaryVertex.z() <= 0.0  ) vz_sign = 1;
    if(PrimaryVertex.z() > 0.0   && PrimaryVertex.z() <= +30.0) vz_sign = 2;
    if(PrimaryVertex.z() > +30.0 && PrimaryVertex.z() <  +70.0) vz_sign = 3;
    // Centrality
    //mRefMultCorr->init(RunId);
    //mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx); 
    const Int_t cent9 = Centrality;
    const Double_t refMultCorr = mVecMesonCut->getRefMultReweight(PrimaryVertex.z(), RefMult);
    const Double_t reweight = mVecMesonCut->getEventWeight(cent9, refMultCorr);
    // runIndex

    const int runIndex = mUtility->findRunIndex(RunId); // find run index for a specific run
    
    // get Track Information
    if(mVecMesonCorr->passTrackEtaNumCut(NumTrackEast,NumTrackWest))
    {
      for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
      {
	mMeson_track = mMeson_event->getTrack(nTracks);
	m2A = mMeson_track->getMass2A();
	m2B = mMeson_track->getMass2B();
	nsA = mMeson_track->getNSigA();
	nsB = mMeson_track->getNSigB();
	dcaA = mMeson_track->getDcaA();
	dcaB = mMeson_track->getDcaB();
	lTrackA = mMeson_track->getTrackA();
	lTrackB = mMeson_track->getTrackB();
	flagA = mMeson_track->getFlagA();
	flagB = mMeson_track->getFlagB();

	Float_t pA = lTrackA.P();
	Float_t pB = lTrackB.P();
	TLorentzVector lTrack = lTrackA + lTrackB; // phi-meson
	Float_t pt_lTrack = lTrack.Perp();
	if(
	    (m2A > 0.16 && m2A < 0.36) && (m2B > 0.16 && m2B < 0.36)
	  )
	{
	  Float_t rapidity_lTrack = lTrack.Rapidity();
	  if(TMath::Abs(rapidity_lTrack) > vmsa::mEtaMax) continue;
	  Float_t InvMass_lTrack = lTrack.M();
	  TVector3 vBetaPhi = -1.0*lTrack.BoostVector(); // get phi beta
	  TLorentzVector lKpRest = lTrackA;

	  //for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++) // systematic loop for dca
          Int_t i_dca = 0;
	  {
	    if( !(mVecMesonCut->passTrackDcaSys(dcaA,dcaB,i_dca,mMode)) ) continue;
	    //for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++) // systematic loop for nSigmaKaon
            Int_t i_sig = 0;
	    {
	      if( !(mVecMesonCut->passTrackSigSys(nsA,nsB,i_sig,mMode)) ) continue;
	      if(mVecMesonCut->passEtaEast(lTrackA)) // K+ neg eta(east)
	      { // Below is West Only
		TVector2 Q2Vector = Q2West;
		// subtract auto-correlation from pos eta(west) event plane
		if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaWest(lTrackB)) // trackB
		{
		  Float_t  w = mVecMesonCorr->getWeight(lTrackB);
		  TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
		  TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_West(cent9,runIndex,vz_sign);
		  Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
		}
		Float_t Res2 = mVecMesonCorr->getResolution2_EP(cent9);
		Float_t Psi2_west = mVecMesonCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign);
		Float_t PhiMinusPsi = lTrack.Phi() - Psi2_west;
                if(PhiMinusPsi < -TMath::Pi()) PhiMinusPsi += (2.0*TMath::Pi());
                if(PhiMinusPsi > +TMath::Pi()) PhiMinusPsi -= (2.0*TMath::Pi());
                mVecMesonHistoManger->FillSys(pt_lTrack,cent9,PhiMinusPsi,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
	      }

	      if(mVecMesonCut->passEtaWest(lTrackA)) // K+ pos eta (west)
	      { // Below is East Only
		TVector2 Q2Vector = Q2East;
		// subtract auto-correlation from pos eta(west) event plane
		if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaEast(lTrackB)) // trackB
		{
		  Float_t  w = mVecMesonCorr->getWeight(lTrackB);
		  TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
		  TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_East(cent9,runIndex,vz_sign);
		  Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
		}
		Float_t Res2 = mVecMesonCorr->getResolution2_EP(cent9);
		Float_t Psi2_east = mVecMesonCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign);
		Float_t PhiMinusPsi = lTrack.Phi() - Psi2_east;
                if(PhiMinusPsi < -TMath::Pi()) PhiMinusPsi += (2.0*TMath::Pi());
                if(PhiMinusPsi > +TMath::Pi()) PhiMinusPsi -= (2.0*TMath::Pi());
                mVecMesonHistoManger->FillSys(pt_lTrack,cent9,PhiMinusPsi,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
	      }
	    }
	  }
	}
      }
    }
  }

  cout << endl;
}

// loop rho meson Same Event
void StVecMesonAna::MakeRho()
{
  //Long64_t start_event_use;
  //Long64_t stop_event_use;

  //start_event_use = mStartEvent;
  //stop_event_use  = mStopEvent;
  mInPut->SetBranchAddress( VM_EVENT_BRANCH, &mMeson_event);
  long numOfEvents = (long)mInPut->GetEntries();
  mInPut->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

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
  for(int i_event = 0; i_event < numOfEvents; ++i_event)
  {
    if(i_event%1000==0) cout << "processing events:  " << i_event << "/" << numOfEvents << endl;
 
    // get Event Header
    PrimaryVertex    = mMeson_event->getPrimaryVertex();
    RunId            = mMeson_event->getRunId();
    EventId          = mMeson_event->getEventId();
    RefMult          = mMeson_event->getRefMult();
    Centrality       = mMeson_event->getCentrality();
    N_prim           = mMeson_event->getN_prim();
    N_non_prim       = mMeson_event->getN_non_prim();
    N_Tof_match      = mMeson_event->getN_Tof_match();
    ZDCx             = mMeson_event->getZDCx(); 
    BBCx             = mMeson_event->getBBCx(); 
    VzVpd            = mMeson_event->getVzVpd();
    NumTrackUsed     = mMeson_event->getNumTracks();
    Q2East           = mMeson_event->getQ2East();
    Q2West           = mMeson_event->getQ2West();
    Q2Full           = mMeson_event->getQ2Full();
    NumTrackEast     = mMeson_event->getNumTrackEast();
    NumTrackWest     = mMeson_event->getNumTrackWest();
    NumTrackFull     = mMeson_event->getNumTrackFull();
    NumTrackFullEast = mMeson_event->getNumTrackFullEast();
    NumTrackFullWest = mMeson_event->getNumTrackFullWest();

    // Initialise Track 
    Float_t m2A = -999.9;
    Float_t m2B = -999.9;
    Float_t nsA = -999.9;
    Float_t nsB = -999.9;
    Float_t dcaA = -999.9;
    Float_t dcaB = -999.9;
    TLorentzVector lTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lTrackB(0.0,0.0,0.0,0.0);
    Int_t flagA = -1;
    Int_t flagB = -1;

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
    //mRunIdEventsDb = StRunIdEventsDb::Instance(vmsa::mEnergyValue[mEnergy],vmsa::mBeamYear[mEnergy]);
    //const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(RunId); // expensive
    // cout << runIndex << endl;

    const int runIndex = mUtility->findRunIndex(RunId); // find run index for a specific run
 

    //if (counter != 0  &&  counter % 1000 == 0)
    //  cout << "." << flush;
    //if (counter != 0  &&  counter % 10000 == 0)
    //{
    //  if((stop_event_use-start_event_use) > 0)
    //  {
    //	Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
    //	cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (VecMesonSpinAlignment) " << flush;
    //  }
    //}

    // get Track Information
    if(mVecMesonCorr->passTrackEtaNumCut(NumTrackEast,NumTrackWest))
    {
      for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
      {
	mMeson_track = mMeson_event->getTrack(nTracks);
	m2A = mMeson_track->getMass2A();
	m2B = mMeson_track->getMass2B();
	nsA = mMeson_track->getNSigA();
	nsB = mMeson_track->getNSigB();
	dcaA = mMeson_track->getDcaA();
	dcaB = mMeson_track->getDcaB();
	lTrackA = mMeson_track->getTrackA();
	lTrackB = mMeson_track->getTrackB();
	flagA = mMeson_track->getFlagA();
	flagB = mMeson_track->getFlagB();

	Float_t pA = lTrackA.P();
	Float_t pB = lTrackB.P();
	TLorentzVector lTrack = lTrackA + lTrackB; // phi-meson
	Float_t pt_lTrack = lTrack.Perp();

	if(
	    // ((fabs(pA) <= 0.65 && m2A < -10) || (m2A > 0 && ((fabs(pA) < 1.5 && m2A > 0.16 && m2A < 0.36) || (fabs(pA) >= 1.5 && m2A > 0.125 && m2A < 0.36)) )) &&
	    // ((fabs(pB) <= 0.65 && m2B < -10) || (m2B > 0 && ((fabs(pB) < 1.5 && m2B > 0.16 && m2B < 0.36) || (fabs(pB) >= 1.5 && m2B > 0.125 && m2B < 0.36)) )) &&
	    // // (pt_lTrack < 0.8 || (pt_lTrack >= 0.8 && ( (m2A > 0.16 && m2A < 0.36) || (m2B > 0.16 && m2B < 0.36)))) &&
	    // (
	    //  ((m2A < -10 && nsA < 3.0 && nsA > -1.5) || (m2A > 0.16 && m2A < 0.36)) &&
	    //  ((m2B < -10 && nsB < 3.0 && nsB > -1.5) || (m2B > 0.16 && m2B < 0.36))
	    // )
	    (m2A > -0.2 && m2A < 0.15) && (m2B > -0.2 && m2B < 0.15)
	  )
	{
	  // Float_t eta_lTrack = lTrack.Eta();
	  // if(TMath::Abs(eta_lTrack) > 1.0) continue;
	  Float_t rapidity_lTrack = lTrack.Rapidity();
	  if(TMath::Abs(rapidity_lTrack) > 1.0) continue;

	  Float_t InvMass_lTrack = lTrack.M();
	  TVector3 vBetaRho = -1.0*lTrack.BoostVector(); // get phi beta
	  TLorentzVector lPipRest = lTrackA;
	  lPipRest.Boost(vBetaRho); // boost pi+ back to phi rest frame
	  TVector3 vPipRest = lPipRest.Vect().Unit(); // pi+ momentum direction in phi rest frame

	  for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++) // systematic loop for dca
	  {
	    if( !(mVecMesonCut->passTrackDcaSys(dcaA,dcaB,i_dca,mMode)) ) continue;
	    mVecMesonHistoManger->FillDcaSys(dcaA,dcaB,i_dca); // fill QA for dcaA and dcaB

	    for(Int_t i_sig = vmsa::nSigPion_start; i_sig < vmsa::nSigPion_stop; i_sig++) // systematic loop for nSigmaPion
	    {
	      if( !(mVecMesonCut->passTrackSigSys(nsA,nsB,i_sig,mMode)) ) continue;
	      mVecMesonHistoManger->FillSigSys(nsA,nsB,i_sig); // fill QA for nsA and nsB 

	      if(mVecMesonCut->passEtaEast(lTrackA)) // Pi+ neg eta(east)
	      { // Below is West Only
		TVector2 Q2Vector = Q2West;
		// subtract auto-correlation from pos eta(west) event plane
		if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaWest(lTrackB)) // trackB
		{
		  Float_t  w = mVecMesonCorr->getWeight(lTrackB);
		  TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
		  TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_West(cent9,runIndex,vz_sign);
		  Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
		}
		Float_t Res2 = mVecMesonCorr->getResolution2_EP(cent9);
		Float_t Psi2_west = mVecMesonCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign);

		TVector3 nQ_West(TMath::Sin(Psi2_west),-1.0*TMath::Cos(Psi2_west),0.0); // normal vector of 2nd Event Plane
		TVector3 nQ = nQ_West.Unit();
		Double_t CosThetaStar = vPipRest.Dot(nQ);
		mVecMesonHistoManger->FillSys(pt_lTrack,cent9,CosThetaStar,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);

		TVector3 nQ_West_EP(TMath::Cos(Psi2_west),TMath::Sin(Psi2_west),0.0); // tangent vector of 2nd Event Plane
		TVector3 nQ_EP = nQ_West_EP.Unit();
		Double_t CosThetaStar_EP = vPipRest.Dot(nQ_EP);
		mVecMesonHistoManger->FillSys_EP(pt_lTrack,cent9,CosThetaStar_EP,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
	      }

	      if(mVecMesonCut->passEtaWest(lTrackA)) // Pi+ pos eta (west)
	      { // Below is East Only
		TVector2 Q2Vector = Q2East;
		// subtract auto-correlation from pos eta(west) event plane
		if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaEast(lTrackB)) // trackB
		{
		  Float_t  w = mVecMesonCorr->getWeight(lTrackB);
		  TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
		  TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_East(cent9,runIndex,vz_sign);
		  Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
		}
		Float_t Res2 = mVecMesonCorr->getResolution2_EP(cent9);
		Float_t Psi2_east = mVecMesonCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign);

		TVector3 nQ_East(TMath::Sin(Psi2_east),-1.0*TMath::Cos(Psi2_east),0.0); // normal vector of 2nd Event Plane
		TVector3 nQ = nQ_East.Unit();
		Double_t CosThetaStar = vPipRest.Dot(nQ);
		mVecMesonHistoManger->FillSys(pt_lTrack,cent9,CosThetaStar,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);

		TVector3 nQ_East_EP(TMath::Cos(Psi2_east),TMath::Sin(Psi2_east),0.0); // tangent vector of 2nd Event Plane
		TVector3 nQ_EP = nQ_East_EP.Unit();
		Double_t CosThetaStar_EP = vPipRest.Dot(nQ_EP);
		mVecMesonHistoManger->FillSys_EP(pt_lTrack,cent9,CosThetaStar_EP,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
	      }
	    }
	  }
	}
      }
    }
  }

  //cout << "." << flush;
  //cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
}

// loop rho meson Same Event
void StVecMesonAna::MakeKStar()
{
  //Long64_t start_event_use;
  //Long64_t stop_event_use;

  //start_event_use = mStartEvent;
  //stop_event_use  = mStopEvent;
  mInPut->SetBranchAddress( VM_EVENT_BRANCH, &mMeson_event);
  long numOfEvents = (long)mInPut->GetEntries();
  mInPut->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

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
  for(int i_event = 0; i_event < numOfEvents; ++i_event)
  {
    if(i_event%1000==0) cout << "processing events:  " << i_event << "/" << numOfEvents << endl;
 
    // get Event Header
    PrimaryVertex    = mMeson_event->getPrimaryVertex();
    RunId            = mMeson_event->getRunId();
    EventId          = mMeson_event->getEventId();
    RefMult          = mMeson_event->getRefMult();
    Centrality       = mMeson_event->getCentrality();
    N_prim           = mMeson_event->getN_prim();
    N_non_prim       = mMeson_event->getN_non_prim();
    N_Tof_match      = mMeson_event->getN_Tof_match();
    ZDCx             = mMeson_event->getZDCx(); 
    BBCx             = mMeson_event->getBBCx(); 
    VzVpd            = mMeson_event->getVzVpd();
    NumTrackUsed     = mMeson_event->getNumTracks();
    Q2East           = mMeson_event->getQ2East();
    Q2West           = mMeson_event->getQ2West();
    Q2Full           = mMeson_event->getQ2Full();
    NumTrackEast     = mMeson_event->getNumTrackEast();
    NumTrackWest     = mMeson_event->getNumTrackWest();
    NumTrackFull     = mMeson_event->getNumTrackFull();
    NumTrackFullEast = mMeson_event->getNumTrackFullEast();
    NumTrackFullWest = mMeson_event->getNumTrackFullWest();

    // Initialise Track 
    Float_t m2A = -999.9;
    Float_t m2B = -999.9;
    Float_t nsA = -999.9;
    Float_t nsB = -999.9;
    Float_t dcaA = -999.9;
    Float_t dcaB = -999.9;
    TLorentzVector lTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lTrackB(0.0,0.0,0.0,0.0);
    Int_t flagA = -1;
    Int_t flagB = -1;

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
    //mRunIdEventsDb = StRunIdEventsDb::Instance(vmsa::mEnergyValue[mEnergy],vmsa::mBeamYear[mEnergy]);
    //const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(RunId); // expensive
    // cout << runIndex << endl;

    const int runIndex = mUtility->findRunIndex(RunId); // find run index for a specific run
 

    //if (counter != 0  &&  counter % 1000 == 0)
    //  cout << "." << flush;
    //if (counter != 0  &&  counter % 10000 == 0)
    //{
    //  if((stop_event_use-start_event_use) > 0)
    //  {
    //	Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
    //	cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (VecMesonSpinAlignment) " << flush;
    //  }
    //}

    // get Track Information
    if(mVecMesonCorr->passTrackEtaNumCut(NumTrackEast,NumTrackWest))
    {
      for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
      {
	mMeson_track = mMeson_event->getTrack(nTracks);
	m2A = mMeson_track->getMass2A();
	m2B = mMeson_track->getMass2B();
	nsA = mMeson_track->getNSigA();
	nsB = mMeson_track->getNSigB();
	dcaA = mMeson_track->getDcaA();
	dcaB = mMeson_track->getDcaB();
	lTrackA = mMeson_track->getTrackA();
	lTrackB = mMeson_track->getTrackB();
	flagA = mMeson_track->getFlagA();
	flagB = mMeson_track->getFlagB();

	Float_t pA = lTrackA.P();
	Float_t pB = lTrackB.P();
	TLorentzVector lTrack = lTrackA + lTrackB; // KStar-meson
	Float_t pt_lTrack = lTrack.Perp();

	if(
	    // ((fabs(pA) <= 0.65 && m2A < -10) || (m2A > 0 && ((fabs(pA) < 1.5 && m2A > 0.16 && m2A < 0.36) || (fabs(pA) >= 1.5 && m2A > 0.125 && m2A < 0.36)) )) &&
	    // ((fabs(pB) <= 0.65 && m2B < -10) || (m2B > 0 && ((fabs(pB) < 1.5 && m2B > 0.16 && m2B < 0.36) || (fabs(pB) >= 1.5 && m2B > 0.125 && m2B < 0.36)) )) &&
	    // // (pt_lTrack < 0.8 || (pt_lTrack >= 0.8 && ( (m2A > 0.16 && m2A < 0.36) || (m2B > 0.16 && m2B < 0.36)))) &&
	    // (
	    //  ((m2A < -10 && nsA < 3.0 && nsA > -1.5) || (m2A > 0.16 && m2A < 0.36)) &&
	    //  ((m2B < -10 && nsB < 3.0 && nsB > -1.5) || (m2B > 0.16 && m2B < 0.36))
	    // )
	    (m2A > 0.16 && m2A < 0.36) && (m2B > -0.2 && m2B < 0.15)
	  )
	{
	  // Float_t eta_lTrack = lTrack.Eta();
	  // if(TMath::Abs(eta_lTrack) > 1.0) continue;
	  Float_t rapidity_lTrack = lTrack.Rapidity();
	  if(TMath::Abs(rapidity_lTrack) > 1.0) continue;

	  Float_t InvMass_lTrack = lTrack.M();
	  TVector3 vBetaKStar = -1.0*lTrack.BoostVector(); // get phi beta
	  TLorentzVector lKRest = lTrackA;
	  lKRest.Boost(vBetaKStar); // boost pi+ back to phi rest frame
	  TVector3 vKRest = lKRest.Vect().Unit(); // pi+ momentum direction in phi rest frame

	  for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++) // systematic loop for dca
	  {
	    if( !(mVecMesonCut->passTrackDcaSys(dcaA,dcaB,i_dca,mMode)) ) continue;
	    mVecMesonHistoManger->FillDcaSys(dcaA,dcaB,i_dca); // fill QA for dcaA and dcaB

	    for(Int_t i_sig = vmsa::nSigPion_start; i_sig < vmsa::nSigPion_stop; i_sig++) // systematic loop for nSigmaPion
	    {
	      if( !(mVecMesonCut->passTrackSigSys(nsA,nsB,i_sig,mMode)) ) continue;
	      mVecMesonHistoManger->FillSigSys(nsA,nsB,i_sig); // fill QA for nsA and nsB 

	      if(mVecMesonCut->passEtaEast(lTrackA)) // Pi+ neg eta(east)
	      { // Below is West Only
		TVector2 Q2Vector = Q2West;
		// subtract auto-correlation from pos eta(west) event plane
		if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaWest(lTrackB)) // trackB
		{
		  Float_t  w = mVecMesonCorr->getWeight(lTrackB);
		  TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
		  TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_West(cent9,runIndex,vz_sign);
		  Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
		}
		Float_t Res2 = mVecMesonCorr->getResolution2_EP(cent9);
		Float_t Psi2_west = mVecMesonCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign);

		TVector3 nQ_West(TMath::Sin(Psi2_west),-1.0*TMath::Cos(Psi2_west),0.0); // normal vector of 2nd Event Plane
		TVector3 nQ = nQ_West.Unit();
		Double_t CosThetaStar = vKRest.Dot(nQ);
		mVecMesonHistoManger->FillSys(pt_lTrack,cent9,CosThetaStar,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);

		TVector3 nQ_West_EP(TMath::Cos(Psi2_west),TMath::Sin(Psi2_west),0.0); // tangent vector of 2nd Event Plane
		TVector3 nQ_EP = nQ_West_EP.Unit();
		Double_t CosThetaStar_EP = vKRest.Dot(nQ_EP);
		mVecMesonHistoManger->FillSys_EP(pt_lTrack,cent9,CosThetaStar_EP,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
	      }

	      if(mVecMesonCut->passEtaWest(lTrackA)) // Pi+ pos eta (west)
	      { // Below is East Only
		TVector2 Q2Vector = Q2East;
		// subtract auto-correlation from pos eta(west) event plane
		if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaEast(lTrackB)) // trackB
		{
		  Float_t  w = mVecMesonCorr->getWeight(lTrackB);
		  TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
		  TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_East(cent9,runIndex,vz_sign);
		  Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
		}
		Float_t Res2 = mVecMesonCorr->getResolution2_EP(cent9);
		Float_t Psi2_east = mVecMesonCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign);

		TVector3 nQ_East(TMath::Sin(Psi2_east),-1.0*TMath::Cos(Psi2_east),0.0); // normal vector of 2nd Event Plane
		TVector3 nQ = nQ_East.Unit();
		Double_t CosThetaStar = vKRest.Dot(nQ);
		mVecMesonHistoManger->FillSys(pt_lTrack,cent9,CosThetaStar,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);

		TVector3 nQ_East_EP(TMath::Cos(Psi2_east),TMath::Sin(Psi2_east),0.0); // tangent vector of 2nd Event Plane
		TVector3 nQ_EP = nQ_East_EP.Unit();
		Double_t CosThetaStar_EP = vKRest.Dot(nQ_EP);
		mVecMesonHistoManger->FillSys_EP(pt_lTrack,cent9,CosThetaStar_EP,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
	      }
	    }
	  }
	}
      }
    }
  }

  //cout << "." << flush;
  //cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
}
//-------------------------------------------------------------------
void StVecMesonAna::Finish()
{
  mFile_OutPut->cd();
  mVecMesonHistoManger->WriteSys(mX_flag,mMode);
  mVecMesonHistoManger->WriteSys_EP(mX_flag,mMode);
  mFile_OutPut->Close();
}
