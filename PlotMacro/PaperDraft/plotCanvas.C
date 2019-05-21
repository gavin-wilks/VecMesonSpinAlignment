#include <string>
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "../../Utility/draw.h"

using namespace std;

static TString LegendA[3] = {"spec","v2","spec+v2"}; 
static TString LegendB[4] = {"0%-80%","0%-10%","10%-40%","40%-80%"};
static TString LegendC[2] = {"particles/anti-particles","particles+anti-pariticles"};
static TString LegendD[4] = {"all","w/o pions","pions+kions+protons","strange"};
static TString LegendE[7] = {"#sqrt{s_{NN}} = 7.7GeV","#sqrt{s_{NN}} = 11.5GeV","#sqrt{s_{NN}} = 39GeV","#sqrt{s_{NN}} = 62.4GeV","#sqrt{s_{NN}} = 19.6GeV","#sqrt{s_{NN}} = 27GeV","#sqrt{s_{NN}} = 200GeV"};
static Int_t DrawMap[7] = {0,1,4,5,2,3,6};

void plotCanvas()
{
  gStyle->SetOptDate(0);

  TCanvas* c_rho00_single = new TCanvas("c_rho00_single","c_rho00_single",1000,20,450,960);
  c_rho00_single->SetFillColor(10);
  c_rho00_single->SetTopMargin(0.2);
  c_rho00_single->SetBottomMargin(0.22);
  c_rho00_single->SetRightMargin(0.05);
  c_rho00_single->SetLeftMargin(0.15);
  c_rho00_single->Divide(1,3,0.0,0.0,10);
  TH1F* h_test_sim_three[3];
  for(Int_t i = 0; i < 3; i ++)
  {
    TString HistName = "h_test_sim_three_";
    HistName += i;
    h_test_sim_three[i] = new TH1F(HistName.Data(),HistName.Data(),1000,-0.5,6);
    for(Int_t j = 0; j < h_test_sim_three[i]->GetNbinsX(); j++)
    {
      h_test_sim_three[i]->SetBinContent(j,-10);
    }
  }

  for(Int_t i = 0; i < 3; i++)
  {
    Float_t scaling_factor = 1.0;
    if(i == 2) scaling_factor = 0.78;
    c_rho00_single->cd(i+1)->SetFillColor(10);
    c_rho00_single->cd(i+1)->SetRightMargin(0.02);
    c_rho00_single->cd(i+1)->SetLeftMargin(0.18);

    c_rho00_single->cd(i+1);
    c_rho00_single->cd(i+1)->SetTicks(1,1);
    c_rho00_single->cd(i+1)->SetGrid(0,0);
    c_rho00_single->cd(i+1)->SetFillColor(10);

    h_test_sim_three[i]->SetStats(0);
    h_test_sim_three[i]->SetTitle("");
    h_test_sim_three[i]->GetXaxis()->SetTitleOffset(0.85/scaling_factor);
    h_test_sim_three[i]->GetYaxis()->SetTitleOffset(0.78/scaling_factor);
    h_test_sim_three[i]->GetYaxis()->SetLabelOffset(0.01/scaling_factor);
    h_test_sim_three[i]->GetXaxis()->SetLabelSize(0.12*scaling_factor);
    h_test_sim_three[i]->GetYaxis()->SetLabelSize(0.12*scaling_factor);
    h_test_sim_three[i]->GetXaxis()->SetTitleSize(0.12*scaling_factor);
    h_test_sim_three[i]->GetYaxis()->SetTitleSize(0.12*scaling_factor);
    h_test_sim_three[i]->GetXaxis()->SetNdivisions(505,'N');
    h_test_sim_three[i]->GetYaxis()->SetNdivisions(505,'N');
    h_test_sim_three[i]->GetXaxis()->CenterTitle();
    h_test_sim_three[i]->GetYaxis()->CenterTitle();
    h_test_sim_three[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_test_sim_three[i]->GetYaxis()->SetTitle("ratio");
    h_test_sim_three[i]->GetYaxis()->SetRangeUser(0.2,0.5);
    h_test_sim_three[i]->DrawCopy("h");
  }


  {
    TString HistName;
    //-----------------------------Set the Canvas---------------------------------
    const Int_t N_x_pads = 2; //number of x-pads
    const Int_t N_y_pads = 3; //number of y-pads
    const Int_t N_total_pads = N_x_pads*N_y_pads; //number of total pads

    TCanvas *c_rho00_double = new TCanvas("c_rho00_double","c_rho00_double",1400,10,800,1000);
    c_rho00_double->SetFillColor(10);
    c_rho00_double->Divide(N_x_pads,N_y_pads,0.0,0.0,10);

    //-----------------------------Set the coordinate-------------------------------  

    TH1D* h_play[N_total_pads];
    for(Int_t i = 0; i < N_total_pads; i++)
    {
      HistName = "h_play_";
      HistName += i;
      h_play[i] = new TH1D(HistName.Data(),HistName.Data(),1000,-0.5,5.5);
      for(Int_t bin_x = 1; bin_x < h_play[i]->GetNbinsX(); bin_x++)
      {
	h_play[i]->SetBinContent(bin_x,-10.0);
      }
      h_play[i]->GetYaxis() -> SetRangeUser(0.2,0.5);
      h_play[i]->GetXaxis() -> SetRangeUser(0.5,5.5);
      h_play[i]->SetStats(0);
      h_play[i]->SetTitle("");
      h_play[i]->GetXaxis()->SetNdivisions(505,'N');
      h_play[i]->GetYaxis()->SetNdivisions(505,'N');
      h_play[i]->GetYaxis()->SetLabelOffset(0.0);
      h_play[i]->GetXaxis()->SetLabelOffset(0.0);
      h_play[i]->GetYaxis()->SetLabelSize(0.0);
      h_play[i]->GetXaxis()->SetLabelSize(0.0);
    }

    //-----------------------------Set the coordinate--end--------------------------  

    Double_t Margin_t = 0.01; // top offset in percent of total size
    Double_t Margin_b = 0.1; // bottom offset in percent of total size
    Double_t Margin_l = 0.19; // left offset in percent of total size
    Double_t Margin_r = 0.01; // right offset in percent of total size
    Double_t delta_y; // fraction in percent per pad in y direction
    Double_t delta_x; // fraction in percent per pad in x direction
    Double_t yUp_array[N_y_pads];
    Double_t yDown_array[N_y_pads];
    Double_t xLeft_array[N_x_pads];
    Double_t xRight_array[N_x_pads];

    delta_y = (1.0-Margin_b-Margin_t)/(Double_t)N_y_pads;
    delta_x = (1.0-Margin_l-Margin_r)/(Double_t)N_x_pads;

    for(Int_t y_pads = 0; y_pads < N_y_pads; y_pads++)//caculate pad size in y direction
    {
      yUp_array[y_pads] = 1.0 - Margin_t - y_pads*delta_y; 
      yDown_array[y_pads] = 1.0 - Margin_t - (y_pads+1.0)*delta_y; 
    }
    for(Int_t x_pads = 0; x_pads < N_x_pads; x_pads++)
    {
      xLeft_array[x_pads] = Margin_l + x_pads*delta_x;
      xRight_array[x_pads] = Margin_l + (x_pads+1)*delta_x;
    }

    xLeft_array[0] = xLeft_array[0] - Margin_l;
    yDown_array[N_y_pads-1] = yDown_array[N_y_pads-1] - Margin_b;
    xRight_array[N_x_pads-1] = xRight_array[N_x_pads-1] + Margin_r;
    yUp_array[0] = yUp_array[0] + Margin_t;


    Double_t scaling_pre = (yUp_array[N_y_pads-1]-yDown_array[N_y_pads-1])*(xRight_array[0]-xLeft_array[0]);
    //  Double_t scaling_pre = (yUp_array[0]-yDown_array[0]);
    //  Double_t scaling_pre = (xRight_array[0]-xLeft_array[0]);
    //    Double_t scaling_pre = 1.0;

    for(Int_t y_pads = 0; y_pads < N_y_pads; y_pads++)
    {
      for(Int_t x_pads = 0; x_pads < N_x_pads; x_pads++)
      {
	Int_t total_pad = (x_pads+1) + y_pads*N_x_pads;
	Double_t scaling_factor = scaling_pre/((yUp_array[y_pads]-yDown_array[y_pads])*(xRight_array[x_pads]-xLeft_array[x_pads]));
	//      Double_t scaling_factor = scaling_pre/(yUp_array[y_pads]-yDown_array[y_pads]);
	//      Double_t scaling_factor = scaling_pre/(xRight_array[x_pads]-xLeft_array[x_pads]);
	//      cout << "scaling_factor = " << scaling_factor << endl;

	c_rho00_double->cd(total_pad);
	c_rho00_double->cd(total_pad)->SetTicks(1,1);
	c_rho00_double->cd(total_pad)->SetGrid(0,0);
	c_rho00_double->cd(total_pad)->SetFillColor(10);;
	c_rho00_double->cd(total_pad)->SetPad(xLeft_array[x_pads],yDown_array[y_pads],xRight_array[x_pads],yUp_array[y_pads]);


	if(x_pads == 0)
	{
	  c_rho00_double->cd(total_pad)->SetLeftMargin(Margin_l/(xRight_array[x_pads]-xLeft_array[x_pads]));
	  h_play[total_pad-1]->GetYaxis()->SetLabelSize(0.08*scaling_factor);
	  h_play[total_pad-1]->GetYaxis()->SetLabelOffset(0.02*scaling_factor);
	  h_play[total_pad-1]->SetTickLength(0.02*scaling_factor,"X");
	  h_play[total_pad-1]->SetTickLength(0.02*scaling_factor,"Y");
	}
	if(x_pads == 0 && y_pads == N_y_pads-1)
	{
	  h_play[total_pad-1]->SetTickLength(0.04*scaling_factor,"Y");
	}
	if(y_pads == N_y_pads-1)
	{
	  c_rho00_double->cd(total_pad)->SetBottomMargin(Margin_b/(yUp_array[y_pads]-yDown_array[y_pads]));
	  h_play[total_pad-1]->GetXaxis()->SetLabelSize(0.08*scaling_factor);
	  //	h_play[total_pad-1]->GetYaxis()->SetLabelSize(0.085*scaling_factor);
	  h_play[total_pad-1]->GetXaxis()->SetLabelOffset(0.015*scaling_factor);
	  h_play[total_pad-1]->GetYaxis()->SetLabelOffset(0.03*scaling_factor);
	}
	if(x_pads == N_x_pads-1)
	{
	  c_rho00_double->cd(total_pad)->SetRightMargin(Margin_r/(xRight_array[x_pads]-xLeft_array[x_pads]));
	  h_play[total_pad-1]->GetXaxis()->SetLabelSize(0.065*scaling_factor);
	  h_play[total_pad-1]->GetXaxis()->SetLabelOffset(0.001*scaling_factor);
	  h_play[total_pad-1]->SetTickLength(0.01*scaling_factor,"X");
	  h_play[total_pad-1]->SetTickLength(0.02*scaling_factor,"Y");
	}
	if(x_pads == N_x_pads-1 && y_pads == N_y_pads-1)
	{
	  h_play[total_pad-1]->SetTickLength(0.04*scaling_factor,"Y");
	}
	if(y_pads == 0)
	{
	  c_rho00_double->cd(total_pad)->SetTopMargin(Margin_t/(yUp_array[y_pads]-yDown_array[y_pads]));
	}

	h_play[total_pad-1]->DrawCopy("PE");
	if(x_pads == 0 && y_pads != N_y_pads-1)
	{
	  plotTopLegend((char*)LegendE[DrawMap[total_pad-1]].Data(),0.0,0.120,0.06*scaling_factor,1,0.0,42,0);
	}
	if(x_pads == 0 && y_pads == N_y_pads-1)
	{
	  plotTopLegend((char*)LegendE[DrawMap[total_pad-1]].Data(),0.0,0.120,0.06*scaling_factor,1,0.0,42,0);
	  h_play[total_pad-1]->SetTickLength(0.03);
	}
	if(x_pads == N_x_pads-1 && y_pads != 0)
	{
	  if(y_pads == 1) plotTopLegend((char*)LegendE[DrawMap[total_pad-1]].Data(),0.0,0.120,0.042*scaling_factor,1,0.0,42,0);
	  if(y_pads == 2) plotTopLegend((char*)LegendE[DrawMap[total_pad-1]].Data(),0.0,0.120,0.050*scaling_factor,1,0.0,42,0);
	}
	if(x_pads == N_x_pads-1 && y_pads == 0)
	{
	  plotTopLegend((char*)LegendE[DrawMap[total_pad-1]].Data(),0.0,0.120,0.042*scaling_factor,1,0.0,42,0);
	}
	if(x_pads == 0 && y_pads == 0)
	{
	  plotTopLegend("Au+Au",0.0,0.100,0.06*scaling_factor,1,0.0,42,0);
	  plotTopLegend((char*)LegendB[0].Data(),0.0,0.080,0.06*scaling_factor,1,0.0,42,0);
	}
      }
    }
  }

}
