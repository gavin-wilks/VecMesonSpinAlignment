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


  mUtility = new StUtility(mEnergy);
  mUtility->initRunIndex(); // initialize std::map for run index
  mRefMultCorr = new StRefMultCorr("refmult");
 
  mVecMesonCorr->InitReCenterCorrection();
  mVecMesonCorr->InitShiftCorrection();
  mVecMesonCorr->InitResolutionCorr();
  mVecMesonHistoManger->InitSys(0,mMode);
  mVecMesonHistoManger->InitSys(1,mMode);
  //mVecMesonHistoManger->InitSys_EP(mX_flag,mMode);

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
	lTrackA = mMeson_track->getTrackA();
	lTrackB = mMeson_track->getTrackB();
	flagA = mMeson_track->getFlagA();
	flagB = mMeson_track->getFlagB();
        tag = mMeson_track->getWhichTagged();

	TLorentzVector lTrack = lTrackA + lTrackB; // phi-meson
	//TLorentzVector lTrackOriginal = lTrackA + lTrackB; // phi-meson
        //TVector3 lTrackMomPhi(lTrackOriginal.Px(),lTrackOriginal.Py(),lTrackOriginal.Pz());
        //TVector3 lTrackMomA(lTrackA.Px(),lTrackA.Py(),lTrackA.Pz());
        //TVector3 lTrackMomB(lTrackB.Px(),lTrackB.Py(),lTrackB.Pz());
        //TVector3 lTrackMomBRotate(-lTrackB.Px(),-lTrackB.Py(),lTrackB.Pz());
        //TVector3 lTrackMomPhiRotate = lTrackMomA + lTrackMomBRotate;
        //double energyA = lTrackA.Energy();
        //double massA = lTrackA.M();
        //double energyB = lTrackB.Energy();
        //double massB = lTrackB.M();
        //double energyPhi = lTrackOriginal.Energy();
 
        //TLorentzVector lTrackBRotate;
        //lTrackBRotate.SetPxPyPzE(-lTrackB.Px(),-lTrackB.Py(),lTrackB.Pz(),energyB);
        double InvMass_lTrack = lTrack.M();
        TVector3 kp_vect = lTrackA.Vect();
        TVector3 km_vect = lTrackB.Vect();
        //double InvMass_lTrack = TMath::Sqrt(massA*massA + massB*massB + 2.0*energyA*energyB - 2.0*lTrackMomA.Dot(lTrackMomBRotate));
        //TLorentzVector lTrackRotate;
        //lTrackRotate.SetPxPyPzE(lTrackMomPhiRotate.X(),lTrackMomPhiRotate.Y(),lTrackMomPhiRotate.Z(),energyPhi);

	//Float_t pA = lTrackA.P();
	//Float_t pB = lTrackB.P();

        Float_t eta_TrackA = lTrackA.PseudoRapidity();
        Float_t eta_TrackB = lTrackB.PseudoRapidity();
        Float_t phi_lTrackA = lTrackA.Phi();
        Float_t phi_lTrackB = lTrackB.Phi();
	Float_t rapidity_lTrack = lTrack.Rapidity();
       	Float_t pt_lTrack = lTrack.Pt();
        Float_t phi_lTrack = lTrack.Phi();
	//Float_t rapidity_lTrack = lTrackRotate.Rapidity();
       	//Float_t pt_lTrack = lTrackMomPhiRotate.Perp();
        //Float_t phi_lTrack = lTrackMomPhiRotate.Phi();
       
        //cout << "rapidity = " << rapidity_lTrack << endl;
        //cout << "pt = " << pt_lTrack << endl;
        //cout << "M = " << InvMass_lTrack << endl << endl;

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
	  //Float_t rapidity_lTrackOriginal = lTrackOriginal.Rapidity();
	  if(TMath::Abs(eta_TrackA) > vmsa::mEtaMax) continue;
	  if(TMath::Abs(eta_TrackB) > vmsa::mEtaMax) continue;
	  //if(TMath::Abs(rapidity_lTrackOriginal) > 1.0) continue;
	  if(TMath::Abs(rapidity_lTrack) > 1.0) continue;
          //cout << "pass rapidity cut" << endl;

	  //TVector3 vBetaPhi = -1.0*lTrackOriginal.BoostVector(); // get phi beta
	  //TLorentzVector lKmRest = lTrackB;
	  //lKmRest.Boost(vBetaPhi); // boost K+ back to phi rest frame
          //if(mX_flag == 1) lKmRest.RotateZ(TMath::Pi()/4.);
          //lKmRest.Boost(-vBetaPhi);
          //TLorentzVector lTrack = lKmRest + lTrackA;

          while(phi_lTrack < -TMath::Pi()) phi_lTrack += 2.0*TMath::Pi();
          while(phi_lTrack > TMath::Pi())  phi_lTrack -= 2.0*TMath::Pi();
	  //Float_t InvMass_lTrack = lTrack.M();
 
          //cout << "Track " << nTracks << endl;
          //cout << "Original k-: " << endl;
          //std::cout << "px: " << lTrackB.Px() << std::endl;
          //std::cout << "py: " << lTrackB.Py() << std::endl;
          //std::cout << "pz: " << lTrackB.Pz() << std::endl;
          //std::cout << "Energy (E): " << lTrackB.E() << std::endl;
          //std::cout << "Mass: " << lTrackB.M() << std::endl;   // Mass of the particle
          //std::cout << "pT: " << lTrackB.Pt() << std::endl;    // Transverse momentum
          //std::cout << "Rapidity: " << lTrackB.Rapidity() << std::endl;
          //std::cout << "Phi: " << lTrackB.Phi() << std::endl;

          //cout << "Rotated k-: " << endl;
          //std::cout << "px: " << lKmRest.Px() << std::endl;
          //std::cout << "py: " << lKmRest.Py() << std::endl;
          //std::cout << "pz: " << lKmRest.Pz() << std::endl;
          //std::cout << "Energy (E): " << lKmRest.E() << std::endl;
          //std::cout << "Mass: " << lKmRest.M() << std::endl;   // Mass of the particle
          //std::cout << "pT: " << lKmRest.Pt() << std::endl;    // Transverse momentum
          //std::cout << "Rapidity: " << lKmRest.Rapidity() << std::endl;
          //std::cout << "Phi: " << lKmRest.Phi() << std::endl << std::endl;


	  //TVector3 vKpRest = lKpRest.Vect().Unit(); // K+ momentum direction in phi rest frame

	  for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++) // systematic loop for dca
	  {
	    if( !(mVecMesonCut->passTrackDcaSys(dcaA,dcaB,i_dca,mMode)) ) continue;
	    //mVecMesonHistoManger->FillDcaSys(dcaA,dcaB,i_dca); // fill QA for dcaA and dcaB

	    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++) // systematic loop for nSigmaKaon
	    {
              if( i_dca != 0 && i_sig != 0) continue;
	      if( !(mVecMesonCut->passTrackSigSys(nsA,nsB,i_sig,mMode)) ) continue;
	      //mVecMesonHistoManger->FillSigSys(nsA,nsB,i_sig); // fill QA for nsA and nsB 

              //double randtheta = gRandom->Uniform(0.,TMath::Pi());
              //double randPsi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());

              //cout << "Rand Theta: " << randtheta << "    rand PSi = " << randPsi << endl;

	      //if(mVecMesonCut->passEtaEast(lTrackA)) // K+ neg eta(east)
	      //{ // Below is West Only
	      //  TVector2 Q2Vector = Q2West;
	      //  // subtract auto-correlation from pos eta(west) event plane
	      //  if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaWest(lTrackB)) // trackB
	      //  {
	      //    Float_t  w = mVecMesonCorr->getWeight(lTrackB);
	      //    TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
	      //    TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_West(cent9,runIndex,vz_sign);
	      //    Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
	      //  }
	      //  Float_t Res2 = mVecMesonCorr->getResolution2_EP(cent9);
	      //  Float_t Psi2_west = mVecMesonCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign);

	      //  TVector3 nQ_West(TMath::Sin(Psi2_west),-1.0*TMath::Cos(Psi2_west),0.0); // normal vector of 2nd Event Plane
	      //  //TVector3 nQ_West(TMath::Sin(randPsi)*TMath::Cos(randtheta),-1.0*TMath::Cos(randPsi)*TMath::Cos(randtheta),TMath::Sin(randtheta)); // normal vector of 2nd Event Plane
	      //  TVector3 nQ = nQ_West.Unit();
	      //  Double_t CosThetaStar = vKpRest.Dot(nQ);
              //  Float_t PhiPsi = mVecMesonCorr->AngleShift(phi_lTrack-Psi2_west);
              //  Float_t Cos2PhiStarPhi = TMath::Cos(2.*(vKpRest.Phi()-phi_lTrack));
              //  TVector3 xprime(TMath::Cos(Psi2_west),TMath::Sin(Psi2_west),0.0);
              //  TVector3 zprime(0.0,0.0,1.0);
              //  Double_t proj_xprime = vKpRest.Dot(xprime);
              //  Double_t proj_zprime = vKpRest.Dot(zprime);
              //  Float_t zxangle = TMath::ATan2(TMath::Abs(proj_zprime),TMath::Abs(proj_xprime));
              //  Float_t phiprime = 0.0;
              //  if(proj_zprime > 0.0)
              //  {
              //    if(proj_xprime > 0.0)  phiprime = 2.0*TMath::Pi()-zxangle; 
              //    if(proj_xprime < 0.0)  phiprime = zxangle;
              //    if(proj_xprime == 0.0) phiprime = 0.0; 
              //  }               
              //  if(proj_zprime < 0.0)
              //  {
              //    if(proj_xprime > 0.0)  phiprime = TMath::Pi()+zxangle;
              //    if(proj_xprime < 0.0)  phiprime = TMath::Pi()-zxangle;
              //    if(proj_xprime == 0.0) phiprime = TMath::Pi();
              //  }
              //  if(proj_zprime == 0.0)
              //  {
              //    if(proj_xprime > 0.0)  phiprime = 3.0*TMath::Pi()/2.0;
              //    if(proj_xprime < 0.0)  phiprime = TMath::Pi()/2.0;
              //    if(proj_xprime == 0.0) phiprime = 0.0;
              //  }
                //Float_t Cos2PhiStarPhi = TMath::Cos(2.*(phistar-phi_lTrack));
		//mVecMesonHistoManger->FillSys(pt_lTrack,rapidity_lTrack,phi_lTrack,cent9,CosThetaStar,PhiPsi,Cos2PhiStarPhi,phiprime,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
                if(!mVecMesonCut->passDipAngleCut(kp_vect,km_vect)) 
                {
                  cout << "Failed dip angle cut" << endl;
                  continue;
                }
		mVecMesonHistoManger->FillSys(pt_lTrack,rapidity_lTrack,phi_lTrack,cent9,0.0,0.0,0.0,0.0,i_dca,i_sig,1.0,InvMass_lTrack,reweight,0,mMode,tag);

		//TVector3 nQ_West_EP(TMath::Cos(Psi2_west),TMath::Sin(Psi2_west),0.0); // tangent vector of 2nd Event Plane
		//TVector3 nQ_EP = nQ_West_EP.Unit();
		//Double_t CosThetaStar_EP = vKpRest.Dot(nQ_EP);
		//mVecMesonHistoManger->FillSys_EP(pt_lTrack,cent9,CosThetaStar_EP,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
	      //}

	      //if(mVecMesonCut->passEtaWest(lTrackA)) // K+ pos eta (west)
	      //{ // Below is East Only
	      //  TVector2 Q2Vector = Q2East;
	      //  // subtract auto-correlation from pos eta(west) event plane
	      //  if(flagB == 0 && mVecMesonCut->passTrackEP(lTrackB,dcaB) && mVecMesonCut->passTrackEtaEast(lTrackB)) // trackB
	      //  {
	      //    Float_t  w = mVecMesonCorr->getWeight(lTrackB);
	      //    TVector2 q2VectorB = mVecMesonCorr->calq2Vector(lTrackB);
	      //    TVector2 q2CorrB   = mVecMesonCorr->getReCenterPar_East(cent9,runIndex,vz_sign);
	      //    Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);
	      //  }
	      //  Float_t Res2 = mVecMesonCorr->getResolution2_EP(cent9);
	      //  Float_t Psi2_east = mVecMesonCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign);

	      //  TVector3 nQ_East(TMath::Sin(Psi2_east),-1.0*TMath::Cos(Psi2_east),0.0); // normal vector of 2nd Event Plane*/
	      //  //TVector3 nQ_East(TMath::Sin(randPsi)*TMath::Cos(randtheta),-1.0*TMath::Cos(randPsi)*TMath::Cos(randtheta),TMath::Sin(randtheta)); // normal vector of 2nd Event Plane
	      //  TVector3 nQ = nQ_East.Unit();
	      //  Double_t CosThetaStar = vKpRest.Dot(nQ);
              //  Float_t PhiPsi = mVecMesonCorr->AngleShift(phi_lTrack-Psi2_east);
              //  Float_t Cos2PhiStarPhi = TMath::Cos(2.*(vKpRest.Phi()-phi_lTrack));
              //  TVector3 xprime(TMath::Cos(Psi2_east),TMath::Sin(Psi2_east),0.0);
              //  TVector3 zprime(0.0,0.0,1.0);
              //  Double_t proj_xprime = vKpRest.Dot(xprime);
              //  Double_t proj_zprime = vKpRest.Dot(zprime);
              //  Float_t zxangle = TMath::ATan2(TMath::Abs(proj_zprime),TMath::Abs(proj_xprime));
              //  Float_t phiprime = 0.0;
              //  if(proj_zprime > 0.0)
              //  {
              //    if(proj_xprime > 0.0)  phiprime = 2.0*TMath::Pi()-zxangle; 
              //    if(proj_xprime < 0.0)  phiprime = zxangle;
              //    if(proj_xprime == 0.0) phiprime = 0.0; 
              //  }               
              //  if(proj_zprime < 0.0)
              //  {
              //    if(proj_xprime > 0.0)  phiprime = TMath::Pi()+zxangle;
              //    if(proj_xprime < 0.0)  phiprime = TMath::Pi()-zxangle;
              //    if(proj_xprime == 0.0) phiprime = TMath::Pi();
              //  }
              //  if(proj_zprime == 0.0)
              //  {
              //    if(proj_xprime > 0.0)  phiprime = 3.0*TMath::Pi()/2.0;
              //    if(proj_xprime < 0.0)  phiprime = TMath::Pi()/2.0;
              //    if(proj_xprime == 0.0) phiprime = 0.0;
              //  }
              //  //TVector3 xprime(TMath::Cos(Psi2_east),TMath::Sin(Psi2_east),0.0);
              //  //TVector3 zprime(0.0,0.0,1.0);
              //  //Double_t proj_xprime = vKpRest.Dot(xprime);
              //  //Double_t proj_zprime = vKpRest.Dot(zprime);
              //  //Float_t zxangle = TMath::ATan2(TMath::Abs(proj_zprime),TMath::Abs(proj_xprime));
              //  //Float_t phistar = 0.0;
              //  //if(proj_zprime > 0.0)
              //  //{
              //  //  if(proj_xprime > 0.0)  phistar = zxangle; 
              //  //  if(proj_xprime < 0.0)  phistar = 2.0*TMath::Pi()-zxangle;
              //  //  if(proj_xprime == 0.0) phistar = 0.0; 
              //  //}               
              //  //if(proj_zprime < 0.0)
              //  //{
              //  //  if(proj_xprime > 0.0)  phistar = TMath::Pi()-zxangle;
              //  //  if(proj_xprime < 0.0)  phistar = TMath::Pi()+zxangle;
              //  //  if(proj_xprime == 0.0) phistar = TMath::Pi();
              //  //}
              //  //if(proj_zprime == 0.0)
              //  //{
              //  //  if(proj_xprime > 0.0)  phistar = TMath::Pi()/2.0;
              //  //  if(proj_xprime < 0.0)  phistar = 3.0*TMath::Pi()/2.0;
              //  //  if(proj_xprime == 0.0) phistar = 0.0;
              //  //}
              //  //Float_t Cos2PhiStarPhi = TMath::Cos(2.*(phistar-phi_lTrack));
		//mVecMesonHistoManger->FillSys(pt_lTrack,rapidity_lTrack,phi_lTrack,cent9,CosThetaStar,PhiPsi,Cos2PhiStarPhi,phiprime,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);

		//TVector3 nQ_East_EP(TMath::Cos(Psi2_east),TMath::Sin(Psi2_east),0.0); // tangent vector of 2nd Event Plane
		//TVector3 nQ_EP = nQ_East_EP.Unit();
		//Double_t CosThetaStar_EP = vKpRest.Dot(nQ_EP);
		//mVecMesonHistoManger->FillSys_EP(pt_lTrack,cent9,CosThetaStar_EP,i_dca,i_sig,Res2,InvMass_lTrack,reweight,mX_flag,mMode);
	      //}
	    }
	  }
	}
      }
    }
  }

  cout << "numOfEventME" << numOfEventsME << endl;
  mInPutME->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry
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

    mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx);
 

    // vz sign 
    int vz_sign = 0;
    if(PrimaryVertex.z() > -70.0 && PrimaryVertex.z() <= -30.0) vz_sign = 0;
    if(PrimaryVertex.z() > -30.0 && PrimaryVertex.z() <= 0.0  ) vz_sign = 1;
    if(PrimaryVertex.z() > 0.0   && PrimaryVertex.z() <= +30.0) vz_sign = 2;
    if(PrimaryVertex.z() < +70.0 && PrimaryVertex.z() >  +30.0) vz_sign = 3;
    const Int_t cent9 = Centrality;
    const Double_t reweight = mRefMultCorr->getWeight();

    const int runIndex = mUtility->findRunIndex(RunId); // find run index for a specific run
    // get Track Information
    if(mVecMesonCorr->passTrackEtaNumCut(NumTrackEast,NumTrackWest))
    {
      for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
      {
	mMeson_track = mMeson_eventME->getTrack(nTracks);
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
        tag = mMeson_track->getWhichTagged();

	TLorentzVector lTrack = lTrackA + lTrackB; // phi-meson
        double InvMass_lTrack = lTrack.M();
        TVector3 kp_vect = lTrackA.Vect();
        TVector3 km_vect = lTrackB.Vect();

        Float_t eta_TrackA = lTrackA.PseudoRapidity();
        Float_t eta_TrackB = lTrackB.PseudoRapidity();
        Float_t phi_lTrackA = lTrackA.Phi();
        Float_t phi_lTrackB = lTrackB.Phi();
	Float_t rapidity_lTrack = lTrack.Rapidity();
       	Float_t pt_lTrack = lTrack.Pt();
        Float_t phi_lTrack = lTrack.Phi();
	{
          if( fabs(eta_TrackA) > mEtaCut || fabs(eta_TrackB) > mEtaCut ) continue;
	  if(TMath::Abs(eta_TrackA) > vmsa::mEtaMax) continue;
	  if(TMath::Abs(eta_TrackB) > vmsa::mEtaMax) continue;
	  if(TMath::Abs(rapidity_lTrack) > 1.0) continue;

          while(phi_lTrack < -TMath::Pi()) phi_lTrack += 2.0*TMath::Pi();
          while(phi_lTrack > TMath::Pi())  phi_lTrack -= 2.0*TMath::Pi();

	  for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++) // systematic loop for dca
	  {
	    if( !(mVecMesonCut->passTrackDcaSys(dcaA,dcaB,i_dca,mMode)) ) continue;
	    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++) // systematic loop for nSigmaKaon
	    {
              if( i_dca != 0 && i_sig != 0) continue;
	      if( !(mVecMesonCut->passTrackSigSys(nsA,nsB,i_sig,mMode)) ) continue;
              if(!mVecMesonCut->passDipAngleCut(kp_vect,km_vect)) 
              {
                cout << "Failed dip angle cut" << endl;
                continue;
              }
	      mVecMesonHistoManger->FillSys(pt_lTrack,rapidity_lTrack,phi_lTrack,cent9,0.0,0.0,0.0,0.0,i_dca,i_sig,1.0,InvMass_lTrack,reweight,1,mMode,tag);
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
  mVecMesonHistoManger->WriteSys(0,mMode);
  mVecMesonHistoManger->WriteSys(1,mMode);
  //mVecMesonHistoManger->WriteSys_EP(mX_flag,mMode);
  mFile_OutPut->Close();
}
