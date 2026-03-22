#include "TPythia6.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TMCParticle.h"
#include "TH1F.h"
#include "TH3F.h"

void scan_phi_mesons(int nevents = 1000, char* jobid = "1") {
    // Load the Pythia 6 library
    gSystem->Load("libEGPythia6.so");

    // Create an instance of TPythia6
    TPythia6* pythia = TPythia6::Instance();

    // Set random seed using clock
    //pythia->SetMRPY(1, 0);

    // Set up the p+p collision at 19.6 GeV center-of-mass energy
    pythia->SetMSEL(1);  // Select hard QCD processes
    pythia->SetPARP(171, 19.6);  // Set the center-of-mass energy (sqrt(s)) in GeV

    // Initialize Pythia
    pythia->Initialize("CMS", "p", "p", 19.6);  // p + p, sqrt(s) = 19.6 GeV

    //pythia->SetPrintEvery(1);  

    //if (!pythia->Initialize("CMS", "p", "p", 19.6)) {
    //  std::cerr << "Pythia initialization failed!" << std::endl;
    //  return; // Handle the error appropriately
    //}

    //pythia->ReadString("HardQCD:all = on"); // Ensure hard QCD processes are enabled
    //pythia->ReadString("Charm:all = on");   // Ensure charm quark production is enabled if you're looking for kaons
    //pythia->ReadString("Strange:all = on");  // Ensure strange quark production is enabled
    //pythia->ReadString("Beams:eCM = 19.6"); // Set the collision energy
    //pythia->ReadString("HardQCD:all = on"); // Enable hard QCD
    //pythia->Initialize(2212, 2212); // Initialize for p+p collisions

    // Create a file to store histograms or data
    TFile* file = new TFile(Form("phi_mesons_%s.root",jobid), "RECREATE");

    // Histograms for \phi-meson analysis
    TH3F* h_phi_ptyphi = new TH3F("h_phi_ptyphi", "Phi-meson pT, y, phi; pT (GeV/c); y; #phi; Counts", 100, 0, 10, 100, -5, 5, 100, 0, 2.0*TMath::Pi());
    TH1F* h_phi_mass = new TH1F("h_phi_mass", "Phi-meson mass; mass (GeV/c^2); Counts", 100, 0.99, 1.05);

    // Number of phi-mesons
    int nphi = 0;

    // Event loop to generate and scan 1000 events
    for (int i = 0; i < nevents; ++i) {
        pythia->GenerateEvent();  // Generate an event
        //if(!pythia->Next()) continue;  // Generate an event

        // Loop over the particles in the event
        for (int j = 0; j < pythia->GetN(); ++j) {
            TParticle* particle = (TParticle*)pythia->GetParticle(j);

            // Check if the particle is a phi-meson (PDG code 333)
            //if (particle->GetPdgCode() == 333) {
                // Get transverse momentum (pT), pseudorapidity (eta), and mass
                cout << "Found " << ++nphi << " phi-mesons" << endl;
                  
                particle->Print();
      
                Double_t pt = particle->Pt();
                Double_t y = particle->Y();
                Double_t phi = particle->Phi();
                Double_t mass = particle->GetMass();

                //cout << "pT = " << pt << ", y = " << y << ", phi = " << phi << endl;
                cout << "px = " << particle->Px() << ", py = " << particle->Py() << ", pz = " << particle->Pz() << endl;

                // Fill histograms
                h_phi_ptyphi->Fill(pt,y,phi);
                h_phi_mass->Fill(mass);
            //}
        }
    }

    // Save histograms to the file and close the file
    file->Write();
    file->Close();

    delete pythia;
}

