#include <TH3F.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMath.h>
#include "../Utility/functions.h"
#include "../Utility/StSpinAlignmentCons.h"

int tableNum[5][9] = {{6,6,5,5,4,3,2,1,1},
                      {0,0,0,0,0,0,0,0,0},
                      {12,12,11,11,10,9,8,7,7},
                      {0,0,0,0,0,0,0,0,0},
                      {18,18,17,17,16,15,14,13,13}};

int tableNumV2[5][9] = {{0,0,0,0,0,0,0,0,0},
                        {0,0,0,0,0,0,0,0,0},
                        {0,0,0,0,0,0,0,0,0},
                        {0,0,0,0,0,0,0,0,0},
                        {239,239,239,239,141,141,141,43,43}};

TF1* readspec(int energy, int centrality)
{
 
  double ptlow[5][11] = { {0.4,0.5,0.6,0.7,0.9,1.0,1.3,0.,0.,0.,0.}, 
                          {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                          {0.4,0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,0.},
                          {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                          {0.4,0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.0} };
 
  double pthigh[5][11] = { {0.5,0.6,0.7,0.9,1.0,1.3,1.7,0.,0.,0.,0.}, 
                           {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                           {0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.5,0.},
                           {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                           {0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.0,4.0} };

  double ptvalues[5][9][11] = { { {0.4541,0.5510,0.6500,0.7973,0.9490,1.1399,1.4818,0.,0.,0.,0.},
                                  {0.4541,0.5510,0.6500,0.7973,0.9490,1.1399,1.4818,0.,0.,0.,0.},
                                  {0.4498,0.5562,0.6518,0.8012,0.9497,1.1433,1.4818,0.,0.,0.,0.},
                                  {0.4498,0.5562,0.6518,0.8012,0.9497,1.1433,1.4818,0.,0.,0.,0.},
                                  {0.4437,0.5493,0.6523,0.8020,0.9498,1.1446,1.4839,0.,0.,0.,0.},
                                  {0.4430,0.5493,0.6521,0.8018,0.9498,1.1442,1.4834,0.,0.,0.,0.},                 
                                  {0.4449,0.5494,0.6531,0.8027,0.9500,1.1456,1.4857,0.,0.,0.,0.},
                                  {0.4461,0.5496,0.6543,0.8037,0.9502,1.1467,1.4876,0.,0.,0.,0.},  
                                  {0.4461,0.5496,0.6543,0.8037,0.9502,1.1467,1.4876,0.,0.,0.,0.} },
                                { {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},                 
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},  
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.} },
                                { {0.4500,0.5556,0.6518,0.8013,0.9497,1.1438,1.4835,1.8392,2.2188,2.8795,0.},
                                  {0.4500,0.5556,0.6518,0.8013,0.9497,1.1438,1.4835,1.8392,2.2188,2.8795,0.},
                                  {0.4500,0.5551,0.6516,0.8008,0.9496,1.1427,1.4813,1.8377,2.2137,2.8589,0.},
                                  {0.4500,0.5551,0.6516,0.8008,0.9496,1.1427,1.4813,1.8377,2.2137,2.8589,0.},
                                  {0.4437,0.5494,0.6523,0.8021,0.9498,1.1446,1.4840,1.8390,2.2165,2.8633,0.},
                                  {0.4441,0.5494,0.6525,0.8022,0.9499,1.1449,1.4845,1.8393,2.2173,2.8660,0.},                 
                                  {0.4459,0.5494,0.6539,0.8035,0.9501,1.1464,1.4871,1.8408,2.2214,2.8798,0.},
                                  {0.4457,0.5494,0.6535,0.8033,0.9501,1.1463,1.4868,1.8406,2.2210,2.8785,0.},  
                                  {0.4457,0.5494,0.6535,0.8033,0.9501,1.1463,1.4868,1.8406,2.2210,2.8785,0.} },
                                { {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},                 
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},  
                                  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.} },
                                { {0.4500,0.5547,0.6515,0.8006,0.9496,1.1427,1.4814,1.8380,2.2150,2.7142,3.3666},
                                  {0.4500,0.5547,0.6515,0.8006,0.9496,1.1427,1.4814,1.8380,2.2150,2.7142,3.3666},
                                  {0.4434,0.5495,0.6523,0.8020,0.9498,1.1446,1.4845,1.8396,2.2191,2.7178,3.3767},
                                  {0.4434,0.5495,0.6523,0.8020,0.9498,1.1446,1.4845,1.8396,2.2191,2.7178,3.3767},
                                  {0.4455,0.5494,0.6535,0.8031,0.9500,1.1460,1.4864,1.8404,2.2203,2.7178,3.3718},
                                  {0.4453,0.5494,0.6531,0.8029,0.9500,1.1458,1.4860,1.8402,2.2198,2.7173,3.3699},                 
                                  {0.4463,0.5494,0.6543,0.8039,0.9502,1.1469,1.4879,1.8412,2.2227,2.7202,3.3798},
                                  {0.4457,0.5494,0.6535,0.8033,0.9501,1.1463,1.4869,1.8407,2.2211,2.7186,3.3745},  
                                  {0.4457,0.5494,0.6535,0.8033,0.9501,1.1463,1.4869,1.8407,2.2211,2.7186,3.3745} } };

  string centlabel = "6080";
  if(centrality >= 2 && centrality <= 3) centlabel = "4060";
  if(centrality >= 4 && centrality <= 4) centlabel = "3040";
  if(centrality >= 5 && centrality <= 5) centlabel = "2030";
  if(centrality >= 6 && centrality <= 6) centlabel = "1020";
  if(centrality >= 7 && centrality <= 8) centlabel = "0010";

  TF1 *f_spec = new TF1("f_spec",pTLevy,0.1,5.0,3);
  TGraphAsymmErrors *g_spec;
  if(energy != 3)
  {
    string InPutSpec;
    if((energy == 4 || energy == 2 || energy == 0)) InPutSpec = "pTspectra/HEPData-ins1378002-v1-root.root";
    TFile *File_Spec = TFile::Open(InPutSpec.c_str());
    cout << "Input spectra" << InPutSpec << endl;
 
    g_spec = (TGraphAsymmErrors*)File_Spec->Get("g_spec");
    if((energy == 4 || energy == 2 || energy == 0)) 
    {
      TDirectory *dir = (TDirectory*) File_Spec->Get(Form("Table %d",tableNum[energy][centrality]));
      dir->cd(); 
      g_spec = (TGraphAsymmErrors*)dir->Get("Graph1D_y1");
      for(int i = 0; i < g_spec->GetN(); i++)
      {
        double x,y;
        g_spec->GetPoint(i,x,y);
        g_spec->SetPoint(i,ptvalues[energy][centrality][i],y);
        g_spec->SetPointError(i,0.0,0.0,g_spec->GetErrorYlow(i),g_spec->GetErrorYhigh(i));
      }
    }

    TF1 *f_Levy = new TF1("f_Levy",Levy,0.1,5.0,3);
    f_Levy->SetParameter(0,1);
    if(energy == 4 && centrality >= 0 && centrality <= 1) f_Levy->SetParameter(0,0.01);
    if(energy == 0 ) f_Levy->SetParameter(0,0.01);
    f_Levy->SetParameter(1,10);
    f_Levy->SetParameter(2,0.1);
    g_spec->Print();
    cout << "Fitting full pT distribution" << endl;
    g_spec->Fit(f_Levy,"N");


    f_spec->SetParameter(0,f_Levy->GetParameter(0));
    f_spec->SetParameter(1,f_Levy->GetParameter(1));
    f_spec->SetParameter(2,f_Levy->GetParameter(2));
  }

  //double dNdy14[9] = {0.009126,0.009126,0.042299,0.042299,0.095731,0.152881,0.223975,0.329114,0.329114};    
  //double Texp14[9] = {0.211483,0.211483,0.228253,0.228253,0.241851,0.245252,0.266930,0.268719,0.268719};

  //TF1 *f_Expo = new TF1("f_Expo",specExp,0.0,5.0,2);
  //if(energy == 3)
  //{
  //  f_spec = new TF1("f_specExp",pTspecExp,0.0,5.0,2);
  //  f_spec->SetParameter(0,dNdy14[centrality]);
  //  f_spec->SetParameter(1,Texp14[centrality]);
  //  f_Expo->SetParameter(0,dNdy14[centrality]);
  //  f_Expo->SetParameter(1,Texp14[centrality]);
  //}
 
  return f_spec;
}

TF1* readv2(int energy, int centrality){
  
  string centlabel = "4080";
  if(centrality >= 4 && centrality <= 6) centlabel = "1040";
  if(centrality >= 7 && centrality <= 8) centlabel = "0010";

  TGraphAsymmErrors *g_v2;
  if(energy != 3)
  {
    string InPutV2 = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Phi_v2_1040.root",vmsa::mBeamEnergy[energy].c_str());
    if((energy == 2 || energy == 0)) InPutV2 = "v2files/HEPData-ins1395151-v2-root.root";
    if(energy == 4) InPutV2 = Form("v2files/OutPhi_v2_Cent%s.root",centlabel.c_str());
    TFile *File_v2 = TFile::Open(InPutV2.c_str());
    std::cout << "v2 file: " << InPutV2 << endl;


    g_v2 = (TGraphAsymmErrors*)File_v2->Get("g_v2");

    if((energy == 2 || energy == 0)) 
    {
      TDirectory *dir = (TDirectory*) File_v2->Get(Form("Table %d",tableNumV2[energy][centrality]));
      dir->cd(); 
      g_v2 = (TGraphAsymmErrors*)dir->Get("Graph1D_y1");
    }
    if(energy == 4) 
    {
      g_v2 = (TGraphAsymmErrors*) File_v2->Get("Graph");
      g_v2->Print();
    }
  }
  
  //if( energy == 3)
  //{
  //  int centidx = 0;
  //  if(centrality >= 0 && centrality <= 3) centidx = 3;
  //  if(centrality >= 4 && centrality <= 6) centidx = 2;
  //  if(centrality >= 7 && centrality <= 8) centidx = 1;
  //  g_v2 = new TGraphAsymmErrors();
  //  for(int ipt = 0; ipt < phiv2_14::ptbins; ipt++)
  //  {
  //    g_v2->SetPoint(ipt, phiv2_14::pt[centidx][ipt], phiv2_14::v2[centidx][ipt]);
  //    double stat = phiv2_14::stat[centidx][ipt];
  //    double sys  = phiv2_14::sys[centidx][ipt];
  //    double totalerr = TMath::Sqrt( stat*stat + sys*sys );
  //    g_v2->SetPointError(ipt, 0.0, 0.0, totalerr, totalerr);  
  //  }

  //}


  TF1 *f_v2 = new TF1("f_v2",v2_pT_FitFunc,0.1,5.0,5);
 // if(centrality <= 3) f_v2->SetRange(0.1,2.8); 
  f_v2->FixParameter(0,2);
  f_v2->SetParameter(1,0.1);
  f_v2->SetParameter(2,0.1);
  f_v2->SetParameter(3,0.1);
  f_v2->SetParameter(4,0.1);
  f_v2->SetLineColor(2);
  f_v2->SetLineWidth(2);
  f_v2->SetLineStyle(2);
  cout << "Fitting v2" << endl;
  
  g_v2->Fit(f_v2,"N");

  TCanvas *c_v2 = new TCanvas("c_v2","c_v2",10,10,800,800);
  c_v2->cd()->SetLeftMargin(0.15);
  c_v2->cd()->SetBottomMargin(0.15);
  c_v2->cd()->SetTicks(1,1);
  c_v2->cd()->SetGrid(0,0);
  TH1F *h_v2 = new TH1F("h_v2","h_v2",100,0.0,10.0);
  for(int i_bin = 1; i_bin < 101; ++i_bin)
  {
    h_v2->SetBinContent(i_bin,-10.0);
    h_v2->SetBinError(i_bin,1.0);
  }
  h_v2->SetTitle("");
  h_v2->SetStats(0);
  h_v2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_v2->GetXaxis()->CenterTitle();
  h_v2->GetYaxis()->SetTitle("v_{2}");
  h_v2->GetYaxis()->CenterTitle();
  h_v2->GetYaxis()->SetRangeUser(0.0,0.50);
  h_v2->Draw("pE");
  g_v2->Draw("pE same");
  f_v2->Draw("l same");
  c_v2->SaveAs(Form("v2_%s_cent%s.pdf",vmsa::mBeamEnergy[energy].c_str(),centlabel.c_str()));

  return f_v2;
}

void calculate_2D_weights_FromEmbed_addy_allcent(int energy = 4, double inputv2 = 500) 
{

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);    


    // Create histograms for flat and desired pT-phi distributions
    const int nBinsPt = 50;
    const int nBinsY  = 30;
    const int nBinsPhi = 50;
    double ptMin = 0.0, ptMax = 5.0;
    double yMin = -1.5, yMax = 1.5;
    double phiMin = 0, phiMax = 2 * TMath::Pi();
    
    double resolution[9]    = {0.221173,    0.296191,    0.411608,    0.535956,    0.630288,    0.675484,    0.648544,    0.539959,    0.377757  };  
    double resolutionerr[9] = {0.000461111, 0.000319516, 0.000202537, 0.000141295, 0.000102773, 8.85656e-05, 9.71308e-05, 0.000193934, 0.00031656};


    //string spectra = "PhiEmbedding_noweights_nomccuts_rapidity_allcent";
    //string filename = "PhiEmbedding_noweights_nomccuts_rapidity_allcent.root";
    //string spectra  = "PhiEmbedding_20240830_All";
    //string filename = "PhiEmbedding_20240830_All.root";
    string spectra  = "PhiEmbedding_20240830_All_NoReweight";
    string filename = "PhiEmbedding_20240830_All_NoReweight.root";

    string inputfile = Form("effaccfiles/%s/%s/%s/%s",vmsa::mPID[0].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),filename.c_str());
    TFile *File_Input = TFile::Open(inputfile.c_str());

    TH3F *h_flat_pt_phi[10];
    TH3F *h_desired_pt_phi[10];
    TH1F *h_flat_pt[10];
    TH1F *h_flat_y[10];
    TH1F *h_flat_phi[10];
    TH1F *h_desired_pt[10];
    TH1F *h_desired_y[10];
    TH1F *h_desired_phi[10];

    double sum = 0.0;
    for(int icent = 0; icent < 10; icent++)
    {
      string flat = Form("h_flat_pt_phi_cent%d",icent);
      h_flat_pt_phi[icent] = (TH3F*) ((TH3F*)File_Input->Get(Form("h_mRc0EffPtYPhiPsi_Cent_%d",icent)))->Clone(flat.c_str());
      cout << "Centrality = " << icent << ", NEntries = " << h_flat_pt_phi[icent]->GetEntries() << endl;
      cout << "NBinsX = " << h_flat_pt_phi[icent]->GetNbinsX() << endl;
      cout << "NBinsY = " << h_flat_pt_phi[icent]->GetNbinsY() << endl;
      cout << "NBinsZ = " << h_flat_pt_phi[icent]->GetNbinsZ() << endl;

      if(icent < 9) sum += h_flat_pt_phi[icent]->GetEntries();
    }
    cout << "Centrality SUM, NEntries = " << sum << endl;
  

    // Create a histogram for the ratio
    TH2F *h_flat_pty[10];   
    TH2F *h_flat_ptphi[10];   
    TH2F *h_flat_yphi[10];   
    for(int icent = 0; icent < 10; icent++)
    {
      string flat = Form("h_flat_pt_cent%d",icent);
      h_flat_pt[icent] = (TH1F*) h_flat_pt_phi[icent]->ProjectionX(flat.c_str(),0,-1,0,-1,"e");    
      flat = Form("h_flat_y_cent%d",icent);
      h_flat_y[icent] = (TH1F*) h_flat_pt_phi[icent]->ProjectionY(flat.c_str(),0,-1,0,-1,"e");    
      flat = Form("h_flat_phi_cent%d",icent);
      h_flat_phi[icent] = (TH1F*) h_flat_pt_phi[icent]->ProjectionZ(flat.c_str(),0,-1,0,-1,"e");    
    
      h_flat_pty[icent] = (TH2F*) h_flat_pt_phi[icent]->Project3D("yx");
      h_flat_ptphi[icent] = (TH2F*) h_flat_pt_phi[icent]->Project3D("zx");
      h_flat_yphi[icent] = (TH2F*) h_flat_pt_phi[icent]->Project3D("zy");
    }

    // Optionally, draw the histograms
    TCanvas *c1 = new TCanvas("c1", "pT-phi Weights", 1200, 800);
    c1->Divide(3,2);
    for(int i = 0; i < 6; i++)
    {
      c1->cd(i+1);
      c1->cd(i+1)->SetLeftMargin(0.15);
      c1->cd(i+1)->SetRightMargin(0.15);
      c1->cd(i+1)->SetBottomMargin(0.15);
      c1->cd(i+1)->SetTicks(1,1);
      c1->cd(i+1)->SetGrid(0,0); 
    }
     
    string outputname = Form("figures/EmbeddingWeights/allcent_pt_y_phi_%s.pdf",spectra.c_str());
    string outputstart = Form("%s[",outputname.c_str());
    string outputstop = Form("%s]",outputname.c_str());

    int centup[10] = {80,70,60,50,40,30,20,10,5,80};
    int centlow[10] = {70,60,50,40,30,20,10,5,0,0};

    c1->Print(outputstart.c_str());
    for(int icent = 0; icent < 10; icent++)
    {
      c1->cd(1);
      h_flat_pt[icent]->SetStats(0);
      h_flat_pt[icent]->SetTitle(Form("Centrality %d-%d%%",centlow[icent],centup[icent]));
      h_flat_pt[icent]->GetXaxis()->SetTitle("p_{T} GeV/c");
      h_flat_pt[icent]->GetYaxis()->SetTitle("Counts");
      h_flat_pt[icent]->Draw("pE");

      c1->cd(2);
      h_flat_y[icent]->SetStats(0);
      h_flat_y[icent]->SetTitle(Form("Centrality %d-%d%%",centlow[icent],centup[icent]));
      h_flat_y[icent]->GetXaxis()->SetTitle("y");
      h_flat_y[icent]->GetYaxis()->SetTitle("Counts");
      h_flat_y[icent]->Draw("pE");

      c1->cd(3);
      h_flat_phi[icent]->SetStats(0);
      h_flat_phi[icent]->SetTitle(Form("Centrality %d-%d%%",centlow[icent],centup[icent]));
      h_flat_phi[icent]->GetXaxis()->SetTitle("#phi-#Psi_{2}");
      h_flat_phi[icent]->GetYaxis()->SetTitle("Counts");
      h_flat_phi[icent]->Draw("pE");

      c1->cd(4);
      h_flat_pty[icent]->SetStats(0);
      h_flat_pty[icent]->SetTitle(Form("Centrality %d-%d%%",centlow[icent],centup[icent]));
      h_flat_pty[icent]->GetXaxis()->SetTitle("p_{T} GeV/c");
      h_flat_pty[icent]->GetYaxis()->SetTitle("y");
      h_flat_pty[icent]->Draw("colz");

      c1->cd(5);
      h_flat_ptphi[icent]->SetStats(0);
      h_flat_ptphi[icent]->SetTitle(Form("Centrality %d-%d%%",centlow[icent],centup[icent]));
      h_flat_ptphi[icent]->GetXaxis()->SetTitle("p_{T} GeV/c");
      h_flat_ptphi[icent]->GetYaxis()->SetTitle("#phi-#Psi_{2}");
      h_flat_ptphi[icent]->Draw("colz");

      c1->cd(6);
      h_flat_yphi[icent]->SetStats(0);
      h_flat_yphi[icent]->SetTitle(Form("Centrality %d-%d%%",centlow[icent],centup[icent]));
      h_flat_yphi[icent]->GetXaxis()->SetTitle("y");
      h_flat_yphi[icent]->GetYaxis()->SetTitle("#phi-#Psi_{2}");
      h_flat_yphi[icent]->Draw("colz");


      c1->Update();
      c1->Print(outputname.c_str());
    }
    c1->Print(outputstop.c_str());   

}

