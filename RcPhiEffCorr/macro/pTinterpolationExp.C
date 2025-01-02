#include <iomanip>
#include <iostream>
#include <string> 
#include <map>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TGraphAsymmErrors.h"
#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "StRoot/Utility/functions.h"

using namespace std;

// 7.7, 11.5, 19.6, 27, 39 GeV

void pTinterpolationExp()
{

  double snn[5] = {7.7,11.5,19.6,27,39};

  int tableNum[5][9] = {{6,6,5,5,4,3,2,1,1},
                        {12,12,11,11,10,9,8,7,7},
                        {18,18,17,17,16,15,14,13,13},
                        {24,24,23,23,22,21,20,19,19},
                        {30,30,29,29,28,27,26,25,25}};

 
  double ptlow[5][12] = { {0.4,0.5,0.6,0.7,0.9,1.0,1.3,0. ,0. ,0. ,0. ,0. }, 
                          {0.4,0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,0. ,0. },
                          {0.4,0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.0,0. },
                          {0.4,0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.0,4.0},
                          {0.4,0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.0,4.0} };
 
  double pthigh[5][12] = { {0.5,0.6,0.7,0.9,1.0,1.3,1.7,0. ,0. ,0. ,0. ,0. }, 
                           {0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.5,0. ,0. },
                           {0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.5,4.0,0. },
                           {0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.5,4.0,5.0},
                           {0.5,0.6,0.7,0.9,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0} };

  double ptvalues[5][9][12] = { { {0.4541,0.5510,0.6500,0.7973,0.9490,1.1399,1.4818,0.    ,0.    ,0.    ,0.    ,0.    },
                                  {0.4541,0.5510,0.6500,0.7973,0.9490,1.1399,1.4818,0.    ,0.    ,0.    ,0.    ,0.    },
                                  {0.4498,0.5562,0.6518,0.8012,0.9497,1.1433,1.4818,0.    ,0.    ,0.    ,0.    ,0.    },
                                  {0.4498,0.5562,0.6518,0.8012,0.9497,1.1433,1.4818,0.    ,0.    ,0.    ,0.    ,0.    },
                                  {0.4437,0.5493,0.6523,0.8020,0.9498,1.1446,1.4839,0.    ,0.    ,0.    ,0.    ,0.    },
                                  {0.4430,0.5493,0.6521,0.8018,0.9498,1.1442,1.4834,0.    ,0.    ,0.    ,0.    ,0.    },                 
                                  {0.4449,0.5494,0.6531,0.8027,0.9500,1.1456,1.4857,0.    ,0.    ,0.    ,0.    ,0.    },
                                  {0.4461,0.5496,0.6543,0.8037,0.9502,1.1467,1.4876,0.    ,0.    ,0.    ,0.    ,0.    },  
                                  {0.4461,0.5496,0.6543,0.8037,0.9502,1.1467,1.4876,0.    ,0.    ,0.    ,0.    ,0.    } },
                                { {0.4500,0.5556,0.6518,0.8013,0.9497,1.1438,1.4835,1.8392,2.2188,2.8795,0.    ,0.    },
                                  {0.4500,0.5556,0.6518,0.8013,0.9497,1.1438,1.4835,1.8392,2.2188,2.8795,0.    ,0.    },
                                  {0.4500,0.5551,0.6516,0.8008,0.9496,1.1427,1.4813,1.8377,2.2137,2.8589,0.    ,0.    },
                                  {0.4500,0.5551,0.6516,0.8008,0.9496,1.1427,1.4813,1.8377,2.2137,2.8589,0.    ,0.    },
                                  {0.4437,0.5494,0.6523,0.8021,0.9498,1.1446,1.4840,1.8390,2.2165,2.8633,0.    ,0.    },
                                  {0.4441,0.5494,0.6525,0.8022,0.9499,1.1449,1.4845,1.8393,2.2173,2.8660,0.    ,0.    },                 
                                  {0.4459,0.5494,0.6539,0.8035,0.9501,1.1464,1.4871,1.8408,2.2214,2.8798,0.    ,0.    },
                                  {0.4457,0.5494,0.6535,0.8033,0.9501,1.1463,1.4868,1.8406,2.2210,2.8785,0.    ,0.    },  
                                  {0.4457,0.5494,0.6535,0.8033,0.9501,1.1463,1.4868,1.8406,2.2210,2.8785,0.    ,0.    } },
                                { {0.4500,0.5547,0.6515,0.8006,0.9496,1.1427,1.4814,1.8380,2.2150,2.7142,3.3666,0.    },
                                  {0.4500,0.5547,0.6515,0.8006,0.9496,1.1427,1.4814,1.8380,2.2150,2.7142,3.3666,0.    },
                                  {0.4434,0.5495,0.6523,0.8020,0.9498,1.1446,1.4845,1.8396,2.2191,2.7178,3.3767,0.    },
                                  {0.4434,0.5495,0.6523,0.8020,0.9498,1.1446,1.4845,1.8396,2.2191,2.7178,3.3767,0.    },
                                  {0.4455,0.5494,0.6535,0.8031,0.9500,1.1460,1.4864,1.8404,2.2203,2.7178,3.3718,0.    },
                                  {0.4453,0.5494,0.6531,0.8029,0.9500,1.1458,1.4860,1.8402,2.2198,2.7173,3.3699,0.    },                 
                                  {0.4463,0.5494,0.6543,0.8039,0.9502,1.1469,1.4879,1.8412,2.2227,2.7202,3.3798,0.    },
                                  {0.4457,0.5494,0.6535,0.8033,0.9501,1.1463,1.4869,1.8407,2.2211,2.7186,3.3745,0.    },  
                                  {0.4457,0.5494,0.6535,0.8033,0.9501,1.1463,1.4869,1.8407,2.2211,2.7186,3.3745,0.    } }, 
                                { {0.4497,0.5555,0.6516,0.8009,0.9496,1.1434,1.4832,1.8393,2.2195,2.7197,3.3887,4.3973},
                                  {0.4497,0.5555,0.6516,0.8009,0.9496,1.1434,1.4832,1.8393,2.2195,2.7197,3.3887,4.3973},
                                  {0.4430,0.5493,0.6521,0.8018,0.9498,1.1444,1.4848,1.8400,2.2210,2.7206,3.3889,4.3962},
                                  {0.4430,0.5493,0.6521,0.8018,0.9498,1.1444,1.4848,1.8400,2.2210,2.7206,3.3889,4.3962},
                                  {0.4457,0.5495,0.6535,0.8035,0.9501,1.1463,1.4870,1.8407,2.2213,2.7188,3.3749,4.3698},
                                  {0.4461,0.5496,0.6543,0.8039,0.9502,1.1468,1.4877,1.8411,2.2224,2.7199,3.3789,4.3738}, 
                                  {0.4465,0.5496,0.6551,0.8043,0.9502,1.1472,1.4885,1.8416,2.2235,2.7210,3.3828,4.3776},
                                  {0.4465,0.5496,0.6555,0.8045,0.9502,1.1474,1.4887,1.8417,2.2238,2.7213,3.3839,4.3787}, 
                                  {0.4465,0.5496,0.6555,0.8045,0.9502,1.1474,1.4887,1.8417,2.2238,2.7213,3.3839,4.3787} },
                                { {0.4426,0.5493,0.6520,0.8016,0.9498,1.1444,1.4848,1.8401,2.2216,2.7214,3.3937,4.4011},
                                  {0.4426,0.5493,0.6520,0.8016,0.9498,1.1444,1.4848,1.8401,2.2216,2.7214,3.3937,4.4011},
                                  {0.4449,0.5494,0.6527,0.8027,0.9500,1.1457,1.4865,1.8409,2.2231,2.7223,3.3944,4.3989},
                                  {0.4449,0.5494,0.6527,0.8027,0.9500,1.1457,1.4865,1.8409,2.2231,2.7223,3.3944,4.3989},
                                  {0.4469,0.5497,0.6563,0.8051,0.9503,1.1478,1.4894,1.8421,2.2250,2.7225,3.3879,4.3828},
                                  {0.4449,0.5497,0.6563,0.8051,0.9503,1.1479,1.4895,1.8421,2.2250,2.7225,3.3881,4.3830}, 
                                  {0.4469,0.5497,0.6563,0.8051,0.9503,1.1479,1.4895,1.8421,2.2250,2.7226,3.3882,4.3831},
                                  {0.4467,0.5497,0.6559,0.8049,0.9502,1.1477,1.4892,1.8419,2.2245,2.7221,3.3865,4.3813}, 
                                  {0.4467,0.5497,0.6559,0.8049,0.9502,1.1477,1.4892,1.8419,2.2245,2.7221,3.3865,4.3813} } };
   
  double par[5][9][3] = {0.0}; // 5 energies, 9 centrality bins, 3 parameters each

  std::string InPutSpec = "pTdata/HEPData-ins1378002-v1-root.root";
  TFile *File_Spec = TFile::Open(InPutSpec.c_str());
  cout << "Input spectra" << InPutSpec << endl;
 
  TGraphAsymmErrors *g_spec[5][9]; // 5 energies, 9 centrality bins
  TGraphAsymmErrors *g_specerror[5][9]; // 5 energies, 9 centrality bins
  TGraphAsymmErrors *g_par[9][3]; // 9 centralities, 3 parameters, energy is the x variable

  TCanvas *cl = new TCanvas("c_l","c_l",10,10,1200,1200);
  cl->SetFillColor(0);
  cl->SetGrid(0,0);
  cl->SetTitle(0);
  cl->SetBottomMargin(0.15);
  cl->SetLeftMargin(0.15);
  cl->SetLogy();

  for(int ien = 0; ien < 5; ien++)
  {
    std::string outputname = Form("figures/ptspecExpFits_en%d.pdf",ien);
    std::string output_start = Form("%s[",outputname.c_str());
    cl->Print(output_start.c_str());
    for(int icent = 0; icent < 9; icent++)
    { 
      TDirectory *dir = (TDirectory*) File_Spec->Get(Form("Table %d",tableNum[ien][icent]));
      dir->cd(); 
      g_spec[ien][icent] = (TGraphAsymmErrors*)dir->Get("Graph1D_y1");
      g_specerror[ien][icent] = (TGraphAsymmErrors*)dir->Get("Graph1D_y2");
      for(int i = 0; i < g_spec[ien][icent]->GetN(); i++)
      {
        double x,y;
        g_spec[ien][icent]->GetPoint(i,x,y);
        //cout << "i = " << i << ", x = " << x << ", y = " << y << endl;
        double xerr,yerr;
        g_specerror[ien][icent]->GetPoint(i,xerr,yerr);
        //cout << "i = " << i << ", xerr = " << xerr << ", yerr = " << yerr << endl;
        g_spec[ien][icent]->SetPoint(i,ptvalues[ien][icent][i],y);
        //g_spec[ien][icent]->SetPointError(i,fabs(ptvalues[ien][icent][i]-ptlow[ien][i]),fabs(ptvalues[ien][icent][i]-pthigh[ien][i]),g_spec[ien][icent]->GetErrorYlow(i),g_spec[ien][icent]->GetErrorYhigh(i));
        //g_spec[ien][icent]->SetPointError(i,0.0,0.0,g_spec[ien][icent]->GetErrorYlow(i),g_spec[ien][icent]->GetErrorYhigh(i));
        //g_spec[ien][icent]->SetPointError(i,fabs(ptvalues[ien][icent][i]-ptlow[ien][i]),fabs(ptvalues[ien][icent][i]-pthigh[ien][i]),yerr,yerr);
        g_spec[ien][icent]->SetPointError(i,0.0,0.0,yerr,yerr);
      }
      TF1 *f_specExp = new TF1("f_specExp",specExp,0.1,5.0,3);
      f_specExp->SetParameter(0,0);
      f_specExp->SetParameter(1,10);
      f_specExp->SetParameter(2,0.2);
      f_specExp->SetLineStyle(2);
      f_specExp->SetLineColor(4);
      f_specExp->SetLineWidth(2);
      g_spec[ien][icent]->Fit(f_specExp,"N");

      g_spec[ien][icent]->SetMarkerStyle(20);
      g_spec[ien][icent]->GetHistogram()->SetTitle(Form("Centrality %d",icent));
      g_spec[ien][icent]->GetHistogram()->GetXaxis()->SetTitle(Form("p_{T}"));
      g_spec[ien][icent]->Draw("APE");
      f_specExp->Draw("same");
       
      cl->Update();
      cl->Print(outputname.c_str());

      for(int ipar = 0; ipar < 2; ipar++)
      {
        if(ien == 0) g_par[icent][ipar] = new TGraphAsymmErrors();
        g_par[icent][ipar]->SetPoint(ien, snn[ien], f_specExp->GetParameter(ipar));
        g_par[icent][ipar]->SetPointError(ien, 0.0, 0.0, f_specExp->GetParError(ipar), f_specExp->GetParError(ipar));
      }
      cout << "Finished Energy " << ien << " Centrality " << icent << endl;
    }
    string output_stop = Form("%s]",outputname.c_str());
    cl->Print(output_stop.c_str()); // close pdf file
  }
  

   
  TCanvas *c1 = new TCanvas("c_pt","c_pt",10,10,1200,1200);
  c1->SetFillColor(0);
  c1->SetGrid(0,0);
  c1->SetTitle(0);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.15);
  //c1->SetLogy();

  double Texp[9] = {0.0};
  double dNdy[9] = {0.0};

  TF1 *f_spec[9]; 
  for(int ipar = 0; ipar < 2; ipar++)
  {
    std::string outputname = Form("figures/ptExp_par%d.pdf",ipar);
    std::string output_start = Form("%s[",outputname.c_str());
    c1->Print(output_start.c_str());
    for(int icent = 0; icent < 9; icent++)
    {
      TF1 *f_parfit;
      f_parfit = new TF1("f_parfit","[0]+[1]*x+[2]*x*x",5,40);
      g_par[icent][ipar]->SetMarkerStyle(20);
      g_par[icent][ipar]->GetHistogram()->SetTitle(Form("Centrality %d",icent));
      g_par[icent][ipar]->GetHistogram()->GetXaxis()->SetTitle(Form("#sqrt{s_{NN}}"));
      g_par[icent][ipar]->Draw("APE");
      g_par[icent][ipar]->Fit(f_parfit,"NQ");

      if(ipar == 0) f_spec[icent] = new TF1("spec",pTspecExp,0.1,5.0,2); 
      f_spec[icent]->SetParameter(ipar, f_parfit->Eval(14.6));

      if(ipar == 0) dNdy[icent] = f_parfit->Eval(14.6);
      if(ipar == 1) Texp[icent] = f_parfit->Eval(14.6);

      f_parfit->SetLineStyle(2);
      f_parfit->SetLineColor(4);
      f_parfit->SetLineWidth(2);
      f_parfit->Draw("same");
      c1->Update();
      c1->Print(outputname.c_str());
    }
    string output_stop = Form("%s]",outputname.c_str());
    c1->Print(output_stop.c_str()); // close pdf file
  }

  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-0.05,9.95);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,1.0);
  }
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,5.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.05);
  h_frame->GetXaxis()->SetTitleOffset(1.2);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("dN/dy");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();

  TCanvas *cspec = new TCanvas("c_spec","c_spec",10,10,1200,1200);
  cspec->SetFillColor(0);
  cspec->SetGrid(0,0);
  cspec->SetTitle(0);
  cspec->SetBottomMargin(0.15);
  cspec->SetLeftMargin(0.15);
  cspec->SetLogy();

  std::string outputname = Form("figures/spectraFor14p6.pdf");
  std::string output_start = Form("%s[",outputname.c_str());
  cspec->Print(output_start.c_str());
  
  cout << "Centrality    dNdy         Texp" << endl;
  for(int icent = 0; icent < 9; icent++) 
  {
    cout << icent << "             " << std::fixed << std::showpoint << std::setprecision(6) << dNdy[icent] << "     " << Texp[icent] << endl;
  }

  for(int icent = 0; icent < 9; icent++)
  {
    h_frame->SetTitle(Form("Centrality %d",icent));
    h_frame->GetYaxis()->SetRangeUser(f_spec[icent]->GetMinimum()*0.5,f_spec[icent]->GetMaximum()*1.1);
    h_frame->DrawCopy("pE");
    f_spec[icent]->SetLineStyle(2);
    f_spec[icent]->SetLineColor(4);
    f_spec[icent]->SetLineWidth(2);
    f_spec[icent]->Draw("same");

    cspec->Update();
    cspec->Print(outputname.c_str());
  }
  string output_stop = Form("%s]",outputname.c_str());
  cspec->Print(output_stop.c_str()); // close pdf file


  //TF1 *f_spec = new TF1("f_spec",pTLevy,vmsa::ptMin,vmsa::ptMax,3);
  //TF1 *f_spec = new TF1("f_spec",pTLevy,pt_low[mPtBin], pt_high[mPtBin], 3);
  //f_spec->SetParameter(0,f_Levy->GetParameter(0));
  //f_spec->SetParameter(1,f_Levy->GetParameter(1));
  //f_spec->SetParameter(2,f_Levy->GetParameter(2));
  //f_spec->SetLineStyle(2);
  //f_spec->SetLineColor(2);
  //f_spec->SetLineWidth(2);


  //c1->SetLogy();
  //g_spec->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  //g_spec->GetXaxis()->SetLabelSize(0.05);
  //g_spec->GetXaxis()->SetTitleSize(0.05);
  //g_spec->GetXaxis()->SetTitleOffset(1.0);
  //g_spec->GetYaxis()->SetTitle("d^{2}N/2#pip_{T}dp_{T}dy");
  //g_spec->GetYaxis()->SetLabelSize(0.05);
  //g_spec->GetYaxis()->SetTitleSize(0.05);
  //g_spec->GetYaxis()->SetTitleOffset(1.1);
  //g_spec->Draw("ap");
  //f_Levy->Draw("same");
  //c1->SaveAs("pt.pdf");

  //TCanvas *c10 = new TCanvas("c_ptdist","c_ptdist",10,10,800,800);
  //c10->cd();
  //c10->SetFillColor(0);
  //c10->SetGrid(0,0);
  //c10->SetTitle(0);
  //c10->SetBottomMargin(0.15);
  //c10->SetLeftMargin(0.15);
  //f_spec->Draw();
  //f_spec->SaveAs("indivitualptspec.pdf");
  
}

