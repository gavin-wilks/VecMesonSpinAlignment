#include "StRoot/StEffMcPhi/StEffMcPhi.h"
//#include "StRoot/StEffMcPhi/StEffHistManger.h"
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include "TNtuple.h"
#include "TBranch.h"
#include "TMath.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StMessMgr.h"
#include "TFile.h"
#include "TChain.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TString.h"

using namespace std;

ClassImp(StEffMcPhi)

int StEffMcPhi::mInPut_flag = 1;

StEffMcPhi::StEffMcPhi(const char* input, int Energy, long StartEvent, long StopEvent, int PID, int mode, int etamode, int inputpt, int startpt, int stoppt)
{
  energy = Energy;
  pid = PID;
  mInputPt = inputpt;
  mStartPt = startpt;
  mStopPt = stoppt;
  mEtaMode = etamode;
  mMode = mode;
  //cout << "mode = " << mode << endl; 
  //cout << "mMode = " << mMode << endl; 
  mList = input;


  // mode = 0 ==> centrality dependence bins, mode = 1 ==> rapidity dependence bins
  
  // string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Efficiency/Eff_%s_SingleKaon_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str(),vmsa::mCuts[cut].c_str());
  std::string InPutFile = input;
  //Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/Lists/%s_%s_Tuples%s.list",vmsa::mPID[pid].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),setting);

  SetInPutFile(InPutFile); // set input list

  SetStartEvent(StartEvent); // set start event
  SetStopEvent(StopEvent); // set stop event

  // string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Efficiency/Eff_%s_SingleKaon.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  std::string OutPutFile = Form("Acceptance_%s_%s_Mode%d_EtaMode%d_ptBins%d_%d.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),mode,etamode,startpt,stoppt);
  SetOutPutFile(OutPutFile); // set output file

//mEffCut = new StEffCut();
  mEffHistManger = new StEffHistManger(energy, pid, mode, mStartPt, mStopPt);
}

StEffMcPhi::~StEffMcPhi()
{
}

//------------------------------------------------------------
void StEffMcPhi::SetInPutFile(const std::string inputfile)
{
  mInPutFile = inputfile;
  cout << "Input file was set to: " << mInPutFile.c_str() << endl;
}

void StEffMcPhi::SetOutPutFile(const std::string outputfile)
{
  mOutPutFile = outputfile;
  cout << "Output file was set to: " << mOutPutFile.c_str() << endl;
}

void StEffMcPhi::SetStartEvent(const long StartEvent)
{
  mStartEvent = StartEvent;
  cout << "nStartEvent = " << mStartEvent << endl;
}

void StEffMcPhi::SetStopEvent(const long StopEvent)
{
  mStopEvent = StopEvent;
  cout << "nStopEvent = " << mStopEvent << endl;
}
//------------------------------------------------------------

void StEffMcPhi::Init()
{
  if(gRandom) delete gRandom;
  gRandom = new TRandom3();
  gRandom->SetSeed();

  mEffHistManger->InitHist();

  // initialize the TNtuple
  //mFile_InPut = TFile::Open(mInPutFile.c_str());
  //cout << "OPEN InPut File: " << mInPutFile.c_str() << endl;
  
  //if(pid == 0) mNtuple = (TNtuple*)mFile_InPut->Get("McPhiMeson");
  //if(pid == 2) mNtuple = (TNtuple*)mFile_InPut->Get("McKStarMeson");


  TString Notification = Form("Initializing parameters and input/output");
  cout << Notification.Data() << endl;
  mFile_OutPut = new TFile(mOutPutFile.c_str(),"RECREATE");

  //VM_EVENT_TREE       = (char*)vmsa::vm_tree[mMode].Data();
  //VM_EVENT_BRANCH     = (char*)vmsa::vm_branch[mMode].Data();

  //----------------------------------------------------------------------------------------------------
  // input
  if (mList != NULL)   // if input file is ok
  {
    mInPut_flag = 1;
    TString InFo_List = "Open file list";
    cout << InFo_List.Data() << endl;
    ifstream in(mList);  // input stream
    if(in)
    {
      cout << "input file list is ok" << endl;
      mInPut = new TChain( "McPhiMeson", "McPhiMeson" );
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
	  mInPut->AddFile(addfile.Data(),-1, "McPhiMeson" );
	  Long64_t file_entries = mInPut->GetEntries();
	  cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      TString InFo_Warning = "WARNING: file input is problemtic";
      cout << InFo_Warning.Data() << endl;
      mInPut_flag = 0;
    }
  }

  // Set the input tree
  if (mInPut_flag == 1)
  {
    //cerr << "ERROR: Could not find branch '"
    //  << VM_EVENT_BRANCH << "'in tree!" << endl;
    mInPut->SetBranchAddress("Centrality",&mCentrality);
    mInPut->SetBranchAddress("McPt",&mMcPt);
    mInPut->SetBranchAddress("McP",&mMcP);
    mInPut->SetBranchAddress("McEta",&mMcEta);
    mInPut->SetBranchAddress("McY",&mMcY);
    mInPut->SetBranchAddress("McPhi",&mMcPhi);
    mInPut->SetBranchAddress("McCos",&mMcCos);
    mInPut->SetBranchAddress("McCosTheta",&mMcCosTheta);
    mInPut->SetBranchAddress("McKpEta",&mMcKpEta);
    mInPut->SetBranchAddress("McKmEta",&mMcKmEta);
    mInPut->SetBranchAddress("McCosRP",&mMcCosRP);


    int num_tracks = mInPut->GetEntriesFast();
    cout << "Number of tracks in McPhiMeson = " << num_tracks<< endl;

    if(mStartEvent > num_tracks) mStartEvent = num_tracks;
    if(mStopEvent  > num_tracks) mStopEvent  = num_tracks;
    cout << "New nStartEvent = " << mStartEvent << ", new nStopEvent = " << mStopEvent << endl;

    //mFile_OutPut = new TFile(mOutPutFile.c_str(),"RECREATE");
  }
  //cout << "BEFORE INITIALIZING NTUPLES" << endl;
  // initialize Ntuple
  cout << "FINISHED INIT" << endl;
}

void StEffMcPhi::Make()
{
  long start_event_use = mStartEvent;
  long stop_event_use  = mStopEvent;

  mInPut->SetBranchAddress("Centrality",&mCentrality);
  mInPut->SetBranchAddress("McPt",&mMcPt);
  mInPut->SetBranchAddress("McP",&mMcP);
  mInPut->SetBranchAddress("McEta",&mMcEta);
  mInPut->SetBranchAddress("McY",&mMcY);
  mInPut->SetBranchAddress("McPhi",&mMcPhi);
  mInPut->SetBranchAddress("McCos",&mMcCos);
  mInPut->SetBranchAddress("McCosTheta",&mMcCosTheta);
  mInPut->SetBranchAddress("McKpEta",&mMcKpEta);
  mInPut->SetBranchAddress("McKmEta",&mMcKmEta);
  mInPut->SetBranchAddress("McCosRP",&mMcCosRP);
  mInPut->SetBranchAddress("McBeta",&mMcBeta);
  mInPut->SetBranchAddress("McBetaP",&mMcBetaP);

  //gRandom->SetSeed();
  mInPut->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

  for(long i_track = start_event_use; i_track < stop_event_use; ++i_track)
  {
    if (!mInPut->GetEntry(i_track)) // take track information
      break;  // end of data chunk
    if (floor(10.0*i_track/ static_cast<float>(stop_event_use)) > floor(10.0*(i_track-1)/ static_cast<float>(stop_event_use)))
      cout << "=> processing data: " << 100.0*i_track/ static_cast<float>(stop_event_use) << "%" << endl;
    //cout << "mMcPt = " << mMcPt << endl;
    //cout << mMode << endl;
    //if( mMode == 2 ) cout << "???????????????????mMcPt = " << mMcPt << "    is in range [" << vmsa::pt_low[energy][mStartPt] << "," << vmsa::pt_up[energy][mStopPt] << "]" << endl;
    if(mMode == 0 && (mMcPt <= vmsa::pt_low_cent[energy][mStartPt] || mMcPt > vmsa::pt_up_cent[energy][mStopPt])) continue;
    //if( mMode == 0) cout << "mMcPt = " << mMcPt << "    is in range [" << vmsa::pt_low_cent[energy][mStartPt] << "," << vmsa::pt_up_cent[energy][mStopPt] << "]" << endl;
    if(mMode == 1 && (mMcPt <= vmsa::pt_low_y[energy][mStartPt]    || mMcPt > vmsa::pt_up_y[energy][mStopPt])   ) continue;
    if(mMode == 2 && (mMcPt <= vmsa::pt_low[energy][mStartPt]      || mMcPt > vmsa::pt_up[energy][mStopPt])   ) continue;
    //if( mMode == 2 ) cout << "mMcPt = " << mMcPt << "    is in range [" << vmsa::pt_low[energy][mStartPt] << "," << vmsa::pt_up[energy][mStopPt] << "]" << endl;
    //cout << "Centrality = " << mCentrality << endl;
    //cout << "pT = " << mMcPt << endl;
    //cout << "y = " << mMcY << endl;
    //cout << "BEFORE FILLING HIST MC" << endl;
    
    if(mEtaMode == 0 && TMath::Abs(mMcY) > 1.0) continue; 
    if(mEtaMode == 3 && TMath::Abs(mMcY) > 0.4) continue; 
    if(mEtaMode == 4 && TMath::Abs(mMcY) > 0.6) continue; 
    if(mEtaMode == 5 && TMath::Abs(mMcY) > 0.8) continue; 
    //cout << "AFTER CUTTING ON Y" << endl;
    
    if(mMode != 2) mEffHistManger->FillHistMc(int(mCentrality),mMcPt,mMcY,TMath::Abs(mMcCos),TMath::Abs(mMcCosTheta),TMath::Abs(mMcCosRP),mMcBetaP,mMcBeta);
    if(mMode == 2) mEffHistManger->FillHistMc(9,mMcPt,mMcY,TMath::Abs(mMcCos),TMath::Abs(mMcCosTheta),TMath::Abs(mMcCosRP),mMcBetaP,mMcBeta);
    //cout << "BEFORE CUTTING ON ETA" << endl;

    if(mEtaMode == 0)
    {
      if( fabs(mMcKpEta) > 1.0 || fabs(mMcKmEta) > 1.0 ) continue; // eta cuts for McPhi 
    }
    else if(mEtaMode == 1)
    {
      if(    ! (fabs(mMcKpEta) <= 1.0 && (fabs(mMcKmEta) < 1.5 && fabs(mMcKmEta) > 1.0) ) 
          && ! (fabs(mMcKmEta) <= 1.0 && (fabs(mMcKpEta) < 1.5 && fabs(mMcKpEta) > 1.0) ) ) continue; // eta cuts for McPhi 
    }
    else if(mEtaMode == 2)
    {
      if( ! (fabs(mMcKpEta) < 1.5 && fabs(mMcKpEta) > 1.0 && fabs(mMcKmEta) < 1.5 && fabs(mMcKmEta) > 1.0 ) ) continue; // eta cuts for McPhi 
    }
    else if(mEtaMode == 3)
    {
      if( fabs(mMcKpEta) > 0.4 || fabs(mMcKmEta) > 0.4 ) continue; // eta cuts for McPhi 
    }
    else if(mEtaMode == 4)
    {
      if( fabs(mMcKpEta) > 0.6 || fabs(mMcKmEta) > 0.6 ) continue; // eta cuts for McPhi 
    }
    else if(mEtaMode == 5)
    {
      if( fabs(mMcKpEta) > 0.8 || fabs(mMcKmEta) > 0.8 ) continue; // eta cuts for McPhi 
    }
    //cout << "BEFORE FILLING HIST RC" << endl;
    if(mMode != 2) mEffHistManger->FillHistRc(int(mCentrality),mMcPt,mMcY,TMath::Abs(mMcCos),TMath::Abs(mMcCosTheta),TMath::Abs(mMcCosRP),mMcBetaP,mMcBeta);
    if(mMode == 2) mEffHistManger->FillHistRc(9,mMcPt,mMcY,TMath::Abs(mMcCos),TMath::Abs(mMcCosTheta),TMath::Abs(mMcCosRP),mMcBetaP,mMcBeta);
    //cout << "AFTER EVERYTHING" << endl;
  }

  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  mEffHistManger->CalEffCosThetaStar();

}

void StEffMcPhi::Finish()
{
  mFile_OutPut->cd();
  mEffHistManger->WriteHist();
  mFile_OutPut->Close();
}
