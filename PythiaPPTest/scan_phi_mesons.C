#include "TPythia8.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TH1F.h"
#include "TH3F.h"

void scan_phi_mesons(int nevents = 10000000, char* jobid = "1") {
    // Initialize TPythia8 object
    TPythia8 *pythia = new TPythia8();

    // Enable and set the random seed
    pythia->ReadString("Random:setSeed = on");
    pythia->ReadString(Form("Random:seed = %s",jobid));

    // Set up p+p collision at sqrt(s) = 19.6 GeV
    pythia->ReadString("Beams:idA = 2212");  // Proton A
    pythia->ReadString("Beams:idB = 2212");  // Proton B
    pythia->ReadString("Beams:eCM = 19.6");  // Center-of-mass energy

    // Enable soft and hard QCD processes
    pythia->ReadString("SoftQCD:all = on"); 

    // Set random seed using clock
    //pythia->SetMRPY(1, 0);

    // Set up the p+p collision at 19.6 GeV center-of-mass energy
    //pythia->SetMSEL(1);  // Select hard QCD processes
    //pythia->SetPARP(171, 19.6);  // Set the center-of-mass energy (sqrt(s)) in GeV

    cout << "Initialize?" << endl;
    // Initialize Pythia
    pythia->Initialize(2212, 2212, 19.6);  // p + p, sqrt(s) = 19.6 GeV
    //pythia->Initialize(2112, 2112, 19.6);  // n + n, sqrt(s) = 19.6 GeV
    cout << "YES" << endl;

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
    cout << "Create file?" << endl;
    TFile* file = new TFile(Form("phi_mesons_%s.root",jobid), "RECREATE");
    cout << "YES" << endl;

    // Histograms for \phi-meson analysis
    TH3F* h_phi_ptyphi = new TH3F("h_phi_ptyphi", "Phi-meson pT, y, phi; pT (GeV/c); y; #phi; Counts", 100, 0, 10, 100, -5, 5, 100, 0, 2.0*TMath::Pi());
    TH1F* h_phi_mass = new TH1F("h_phi_mass", "Phi-meson mass; mass (GeV/c^2); Counts", 100, 0.99, 1.05);

    // Number of phi-mesons
    int nphi = 0;

    // Event loop to generate and scan 1000 events
    for (int i = 0; i < nevents; ++i) {
        pythia->GenerateEvent();  // Generate an event

        //cout << "Generated the Event " << endl;
        // Loop over the particles in the event
        for (int j = 0; j < pythia->GetN(); ++j) {
            //cout << "Look for particle " << j << endl;
            TParticle* particle = pythia->GetParticle(j);
            if (!particle) {
              //std::cerr << "Null particle pointer at index " << j << std::endl;
              continue;
            }

            //cout << "Found a particle" << endl;
            // Check if the particle is a phi-meson (PDG code 333)
            //particle->Print();
            //cout << "The PDG code is " << particle->GetPdgCode() << endl;
            if (particle->GetPdgCode() == 333) {
                //cout << "Found a phi-meson" << endl;
                // Get transverse momentum (pT), pseudorapidity (eta), and mass
                if(++nphi % 100 == 0) cout << "Found " << nphi << " phi-mesons" << endl;


                Double_t pt = particle->Pt();
                Double_t y = particle->Y();
                Double_t phi = particle->Phi();
                Double_t mass = particle->GetMass();

                //cout << "pT = " << pt << ", y = " << y << ", phi = " << phi << endl;
                //cout << "px = " << particle->Px() << ", py = " << particle->Py() << ", pz = " << particle->Pz() << endl;

                // Fill histograms
                h_phi_ptyphi->Fill(pt,y,phi);
                h_phi_mass->Fill(mass);
            }
        }
    }

    // Save histograms to the file and close the file
    file->Write();
    file->Close();

    delete pythia;
}

