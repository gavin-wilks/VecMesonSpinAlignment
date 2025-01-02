#include <fstream>
#include <string>
#include <iostream>
#include <THn.h>
#include <TChain.h>
#include <TH3F.h>
#include <TH3F.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TMath.h>

void processKaonTTrees(int energy = 4, bool datarcweight = false, const Char_t *SM, const Char_t *fileinput, const Char_t *jobId) 
{


    TH3F* ptyphibasic[2];

    TH3F* ptynhf[2];
    TH3F* ptphinhf[2];
    TH3F* yphinhf[2];
 
    TH3F* ptynhm[2];
    TH3F* ptphinhm[2];
    TH3F* yphinhm[2];
 
    TH3F* ptynhr[2];
    TH3F* ptphinhr[2];
    TH3F* yphinhr[2];
 
    TH3F* ptyde[2];
    TH3F* ptphide[2];
    TH3F* yphide[2];
 
    TH3F* ptydca[2];
    TH3F* ptphidca[2];
    TH3F* yphidca[2];
 
    const int npt = 35;
    const float ptmin = 0.0, ptmax = 3.5; 
    const int ny = 48;
    const float ymin = -1.2, ymax = 1.2;
    const int neta = 30;
    const float etamin = -1.5, etamax = 1.5;
    const int nphi = 48;
    const float phimin = -TMath::Pi(), phimax = TMath::Pi();

    const int nhf = 50;
    const float nhfmin = -0.5, nhfmax = 100.5; 
    const int nhm = 50;
    const float nhmmin = -0.5, nhmmax = 100.5; 
    const int nhr = 20;
    const float nhrmin = 0.5, nhrmax = 1.1; 
    const int nde = 50;
    const float demin = 0.0, demax = 10.0; 
    const int ndca = 50;
    const float dcamin = 0.0, dcamax = 2.1; 


    Int_t bins[4] = {npt, ny, neta, nphi};
    Double_t xmin[4] = {ptmin, ymin, etamin, phimin};
    Double_t xmax[4] = {ptmax, ymax, etamax, phimax};

    string charge[2] = {"plus","minus"};
    for(int ic = 0; ic < 2; ic++)
    {
      string hist = Form("ptyphibasic_k%s",charge[ic].c_str());
      ptyphibasic[ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, ny, ymin, ymax, nphi, phimin, phimax);
      ptyphibasic[ic]->Sumw2();

      hist = Form("ptynhf_k%s",charge[ic].c_str());
      ptynhf[ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, ny, ymin, ymax, nhf, nhfmin, nhfmax);
      ptynhf[ic]->Sumw2();
      hist = Form("ptphinhf_k%s",charge[ic].c_str());
      ptphinhf[ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nphi, phimin, phimax, nhf, nhfmin, nhfmax);
      ptphinhf[ic]->Sumw2();
      hist = Form("yphinhf_k%s",charge[ic].c_str());
      yphinhf[ic] = new TH3F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nphi, phimin, phimax, nhf, nhfmin, nhfmax);
      yphinhf[ic]->Sumw2();
 
      hist = Form("ptynhm_k%s",charge[ic].c_str());
      ptynhm[ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, ny, ymin, ymax, nhm, nhmmin, nhmmax);
      ptynhm[ic]->Sumw2();
      hist = Form("ptphinhm_k%s",charge[ic].c_str());
      ptphinhm[ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nphi, phimin, phimax, nhm, nhmmin, nhmmax);
      ptphinhm[ic]->Sumw2();
      hist = Form("yphinhm_k%s",charge[ic].c_str());
      yphinhm[ic] = new TH3F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nphi, phimin, phimax, nhm, nhmmin, nhmmax);
      yphinhm[ic]->Sumw2();
 
      hist = Form("ptynhr_k%s",charge[ic].c_str());
      ptynhr[ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, ny, ymin, ymax, nhr, nhrmin, nhrmax);
      ptynhr[ic]->Sumw2();
      hist = Form("ptphinhr_k%s",charge[ic].c_str());
      ptphinhr[ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nphi, phimin, phimax, nhr, nhrmin, nhrmax);
      ptphinhr[ic]->Sumw2();
      hist = Form("yphinhr_k%s",charge[ic].c_str());
      yphinhr[ic] = new TH3F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nphi, phimin, phimax, nhr, nhrmin, nhrmax);
      yphinhr[ic]->Sumw2();

      hist = Form("ptyde_k%s",charge[ic].c_str());
      ptyde[ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, ny, ymin, ymax, nde, demin, demax);
      ptyde[ic]->Sumw2();
      hist = Form("ptphide_k%s",charge[ic].c_str());
      ptphide[ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nphi, phimin, phimax, nde, demin, demax);
      ptphide[ic]->Sumw2();
      hist = Form("yphide_k%s",charge[ic].c_str());
      yphide[ic] = new TH3F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nphi, phimin, phimax, nde, demin, demax);
      yphide[ic]->Sumw2();

      hist = Form("ptydca_k%s",charge[ic].c_str());
      ptydca[ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, ny, ymin, ymax, ndca, dcamin, dcamax);
      ptydca[ic]->Sumw2();
      hist = Form("ptphidca_k%s",charge[ic].c_str());
      ptphidca[ic] = new TH3F(hist.c_str(), hist.c_str(), npt, ptmin, ptmax, nphi, phimin, phimax, ndca, dcamin, dcamax);
      ptphidca[ic]->Sumw2();
      hist = Form("yphidca_k%s",charge[ic].c_str());
      yphidca[ic] = new TH3F(hist.c_str(), hist.c_str(), ny, ymin, ymax, nphi, phimin, phimax, ndca, dcamin, dcamax);
      yphidca[ic]->Sumw2();

    }

    TChain *mKaonTree = new TChain("kaontree");
    // INITIALIZE TTREES 
    std::ifstream filelist(fileinput);
    std::string filename;
    while (std::getline(filelist,filename)) {
      if(!filename.empty()) {
        mKaonTree->Add(filename.c_str());
        std::cout << "Added file: " << filename << std::endl;
      }
    }
    filelist.close();
 
    int   mCent;
    float mWeight;
    int   mCharge;
    float mPt;
    float mRapidity;
    float mEta;
    float mPhi;
    int   mNHitsFit;
    int   mNHitsMax;
    float mDEdx;
    float mDca;
    mKaonTree->SetBranchAddress("cent", &mCent);
    mKaonTree->SetBranchAddress("weight", &mWeight);
    mKaonTree->SetBranchAddress("charge", &mCharge);
    mKaonTree->SetBranchAddress("pt", &mPt);
    mKaonTree->SetBranchAddress("rapidity", &mRapidity);
    mKaonTree->SetBranchAddress("eta", &mEta);
    mKaonTree->SetBranchAddress("phi", &mPhi);
    mKaonTree->SetBranchAddress("nhitsfit", &mNHitsFit);
    mKaonTree->SetBranchAddress("nhitsmax", &mNHitsMax);
    mKaonTree->SetBranchAddress("dedx", &mDEdx);
    mKaonTree->SetBranchAddress("dca", &mDca);

    Long64_t nentries = mKaonTree->GetEntries();
    //Long64_t nentries = 10000;//mKaonTree->GetEntries();
    cout << "nentries = " << nentries << endl;
    for (Long64_t i = 0; i < nentries; i++)
    {
      mKaonTree->GetEntry(i);
      if (i%1000000 == 0)
      {
         cout << "Proccessed " << i << endl;
      }

      //Double_t values[4] = {mPt,mRapidity,mEta,mPhi};
      //ptyetaphibasic[mCharge]->Fill(values,mWeight);
      ptyphibasic[mCharge]->Fill(mPt,mRapidity,mPhi,mWeight);

      ptynhf[mCharge]->Fill(mPt,mRapidity,mNHitsFit,mWeight);
      ptphinhf[mCharge]->Fill(mPt,mPhi,mNHitsFit,mWeight);
      yphinhf[mCharge]->Fill(mRapidity,mPhi,mNHitsFit,mWeight);
 
      ptynhm[mCharge]->Fill(mPt,mRapidity,mNHitsMax,mWeight);
      ptphinhm[mCharge]->Fill(mPt,mPhi,mNHitsMax,mWeight);
      yphinhm[mCharge]->Fill(mRapidity,mPhi,mNHitsMax,mWeight);
 
      ptynhr[mCharge]->Fill(mPt,mRapidity,float(mNHitsFit)/float(mNHitsMax),mWeight);
      ptphinhr[mCharge]->Fill(mPt,mPhi,float(mNHitsFit)/float(mNHitsMax),mWeight);
      yphinhr[mCharge]->Fill(mRapidity,mPhi,float(mNHitsFit)/float(mNHitsMax),mWeight);
 
      ptyde[mCharge]->Fill(mPt,mRapidity,mDEdx,mWeight);
      ptphide[mCharge]->Fill(mPt,mPhi,mDEdx,mWeight);
      yphide[mCharge]->Fill(mRapidity,mPhi,mDEdx,mWeight);
 
      ptydca[mCharge]->Fill(mPt,mRapidity,mDca,mWeight);
      ptphidca[mCharge]->Fill(mPt,mPhi,mDca,mWeight);
      yphidca[mCharge]->Fill(mRapidity,mPhi,mDca,mWeight);
    }
    
    delete mKaonTree;

    TFile *outputFile = new TFile(Form("output_%s_%s.root",SM,jobId), "RECREATE"); // "RECREATE" will overwrite an existing file
    for(int ic = 0; ic < 2; ic++)
    {
      ptyphibasic[ic]->Write();
 
      ptynhf[ic]->Write();
      ptphinhf[ic]->Write();
      yphinhf[ic]->Write();
 
      ptynhm[ic]->Write();
      ptphinhm[ic]->Write();
      yphinhm[ic]->Write();
 
      ptynhr[ic]->Write();
      ptphinhr[ic]->Write();
      yphinhr[ic]->Write();
 
      ptyde[ic]->Write();
      ptphide[ic]->Write();
      yphide[ic]->Write();
 
      ptydca[ic]->Write();
      ptphidca[ic]->Write();
      yphidca[ic]->Write();
     
    }

    outputFile->Close(); 
 
}

