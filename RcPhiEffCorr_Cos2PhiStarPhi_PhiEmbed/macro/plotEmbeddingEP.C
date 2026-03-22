#include "TString.h"
#include "TH1F.h"
#include "TMath.h"
#include <fstream>

void plotEmbeddingEP()
{
  TH1F *h_mEast = new TH1F("h_mEast","h_mEast",360,-1.0*TMath::Pi(),TMath::Pi());
  TH1F *h_mWest = new TH1F("h_mWest","h_mWest",360,-1.0*TMath::Pi(),TMath::Pi());

  TString InPutFile = "/star/u/sunxuhit/AuAu200GeV/SpinAlignment/Phi/RunId/embedding_map_pico.txt"; // embedding map pico
  cout << "Input mapping was set to: " << InPutFile.Data() << endl;
  FILE *fp = fopen(InPutFile.Data(),"r");
  if(fp == NULL)
  {
    perror("Error opening mapping file");
  }
  else
  {
    Int_t runId, eventId;
    Float_t psiEast, psiWest;
    char line[80];

    Int_t line_counter = 0;
    while(fgets(line,80,fp))
    {
      sscanf(&line[0],"%d %d %f %f", &runId, &eventId, &psiEast, &psiWest);
      line_counter++;
      h_mEast->Fill(psiEast);
      h_mWest->Fill(psiWest);
      // cout << "eventId  = " << mEventId[eventId] << ", psiEast = " << mPsiEast[eventId] << ", psiWest = " << mPsiWest[eventId] << endl;
    }
  }
  TString OutPutFile = "/star/u/sunxuhit/AuAu200GeV/SpinAlignment/Phi/RunId/EP_Embedding.root";
  TFile *File_OutPut = new TFile(OutPutFile.Data(),"RECREATE");
  File_OutPut->cd();
  h_mEast->Write();
  h_mWest->Write();
  File_OutPut->Close();
}
