#include <string>
#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "../../Utility/StSpinAlignmentCons.h"

void plotQA_InvMass(int energy = 4, int flag_ME = 0, int flag_PID = 0, int flag_dPID = 0)
{
  TString dPID[2]={"K","pi"};
  string inputlist = Form("../../FileLists/%s_%d/%s_%s_Forest.list",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamYear[energy],vmsa::mPID[flag_PID].c_str(),vmsa::MixEvent[flag_ME].Data());
  //string inputdir = Form("/global/homes/x/xusun/AuAu%s/SpinAlignment/Phi/Forest/",vmsa::mBeamEnergy[energy].c_str());
  //string inputlist = Form("/global/homes/x/xusun/AuAu%s/SpinAlignment/Phi/List/Phi_%s_tree.list",vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[flag_ME].Data());
  TH1F *h_mInvMass_Tot;
  TH2F *h_mDEdxRig_Tot;
  TH2F *h_mM2Rig_Tot;
  TH2F *h_mInvBetaRig_Tot;

  TCanvas *c_InvMass = new TCanvas("c_InvMass","c_InvMass",10,10,1600,400);
  c_InvMass->Divide(4,1);
  for(int i_pad = 0; i_pad < 4; ++i_pad)
  {
    c_InvMass->cd(i_pad+1);
    c_InvMass->cd(i_pad+1)->SetLeftMargin(0.1);
    c_InvMass->cd(i_pad+1)->SetBottomMargin(0.1);
    c_InvMass->cd(i_pad+1)->SetGrid(0,0);
    c_InvMass->cd(i_pad+1)->SetTicks(1,1);
  }
  //string output_start = Form("../../figures/TreeProduction_AuAu%s_%s_QA.pdf[",vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[flag_ME].Data());
  //c_InvMass->Print(output_start.c_str());
  //int i_pad = 0;
  //int i_page = 1;
  int i_file = 0;
  //string outputname = Form("../../figures/TreeProduction_AuAu%s_%s_QA.pdf",vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[flag_ME].Data());
  if (!inputlist.empty())   // if input file is ok
  {
    cout << "Open input file list: " << inputlist.c_str() << endl;
    ifstream in(inputlist.c_str());  // input stream
    if(in)
    {
      cout << "input database file list is ok" << endl;
      char str[255];       // char array for each file name
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  string inputfile;
	  inputfile = str;
	  cout << "open file: " << inputfile.c_str() << endl;
	  TFile *File_InPut = TFile::Open(inputfile.c_str());
	  //cout << "opened file" << endl;
          TH2F *h_mMass2_pt = (TH2F*)File_InPut->Get("Mass2_pt");//->Clone();
	  TH1F *h_mInvMass = (TH1F*)h_mMass2_pt->ProjectionY();//->Clone("h_mInvMass");
          //h_mInvMass->Rebin(2);
          //cout << "Grab Mass2_pt" << endl;
          TH2F *h_mDEdxRig    = (TH2F*)File_InPut->Get(dPID[flag_dPID]+"_dEdx_Rig");//->Clone(); cout << "Grab dEdx_Rig" << endl;
          TH2F *h_mM2Rig      = (TH2F*)File_InPut->Get(dPID[flag_dPID]+"_Mass2_Rig");//->Clone(); cout << "Grab Mass2_Rig" << endl;
          TH2F *h_mInvBetaRig = (TH2F*)File_InPut->Get(dPID[flag_dPID]+"_InvBeta_Rig");//->Clone(); cout << "InvBeta_Rig" << endl;
          if(i_file == 0) 
          { 
            h_mInvMass_Tot = (TH1F*)h_mInvMass->Clone("h_mInvMass_Tot");
            h_mDEdxRig_Tot = (TH2F*)h_mDEdxRig->Clone("h_mDEdxRig_Tot");
            h_mM2Rig_Tot = (TH2F*)h_mM2Rig->Clone("h_mM2Rig_Tot");
            h_mInvBetaRig_Tot = (TH2F*)h_mInvBetaRig->Clone("h_mInvBetaRig_Tot"); //cout << "Loaded first set of files" << endl;
            h_mInvMass_Tot->SetDirectory(0);            
            h_mDEdxRig_Tot->SetDirectory(0); 
            h_mM2Rig_Tot->SetDirectory(0); 
            h_mInvBetaRig_Tot->SetDirectory(0); 
          }
          if(i_file != 0)
          {
            h_mInvMass_Tot->Add(h_mInvMass);
            h_mDEdxRig_Tot->Add(h_mDEdxRig);
            h_mM2Rig_Tot->Add(h_mM2Rig);
            h_mInvBetaRig_Tot->Add(h_mInvBetaRig);
          }
          i_file++;
	  //c_InvMass->cd(i_pad+1);
	  //h_mInvMass->DrawCopy("hE");

	  //i_pad++;
	  //int NumOfTracks = h_mInvMass->GetEntries();
	  //cout << "In page " << i_page << " pad " << i_pad << " with " << NumOfTracks << " tracks!" << endl;
	  //cout << endl;
	  h_mInvMass_Tot->Print();  
          h_mDEdxRig_Tot->Print();          
          h_mM2Rig_Tot->Print();
          h_mInvBetaRig_Tot->Print();

	  h_mMass2_pt->Reset();
	  h_mInvMass->Reset();
          h_mDEdxRig->Reset();   
	  h_mM2Rig->Reset();    
	  h_mInvBetaRig->Reset();
          File_InPut->Close(); //cout << "Closed the file" << endl;
	}
      }
      //cout << "Exited the loop" << endl;
      c_InvMass->cd(1); cout << //"first pad" << endl;
      h_mInvMass_Tot->Draw("hE"); //cout << "Draw first plot" << endl;     
      c_InvMass->cd(2);
      h_mDEdxRig_Tot->Draw("colz"); //cout << "Draw second plot" << endl;
      c_InvMass->cd(3);
      h_mM2Rig_Tot->Draw("colz"); //cout << "Draw third plot" << endl;
      c_InvMass->cd(4); //cout << "fourth pad" << endl;
      h_mInvBetaRig_Tot->Draw("colz"); //cout << "Draw fourth plot" << endl;
      //c_InvMass->Update(); cout << "Is this updated" << endl;
      c_InvMass->SaveAs(Form("../../figures/%s_%s_TreeProduction_AuAu%s_%s_QA.pdf",vmsa::mPID[flag_PID].c_str(),dPID[flag_dPID].Data(),vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[flag_ME].Data())); //cout << "print the stuff" << endl;
    }
    else
    {
      cout << "WARNING: input database file input is problemtic" << endl;
    }
  }
  //cout << "Prepare to print the file" << endl;
  //string output_stop = Form("../../figures/TreeProduction_AuAu%s_%s_QA.pdf]",vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[flag_ME].Data());
  //c_InvMass->Print(output_stop.c_str());
}
