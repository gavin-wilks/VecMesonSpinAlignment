#include <string>

#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>

#include "../../../../draw.h"
#include "../../Utility/StSpinAlignmentCons.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void plotExtFig6_rho00PtPhi_Poly1()
{
  float pt[] = {0.0,
                0.1,
                0.2,
                0.3,
                0.4,
                0.5,
                0.6,
                0.7,
                0.8,
                0.9,
                1.0, 
                1.1,
                1.2,
                1.3,
                1.4,
                1.5,
                1.6,
                1.7,
                1.8,
                1.9,
                2.0, 
                2.1,
                2.2,
                2.3,
                2.4,
                2.5,
                2.6,
                2.7,
                2.8,
                2.9,
                3.0, 
                3.1,
                3.2,
                3.3,
                3.4,
                3.5,
                3.6,
                3.7,
                3.8,
                3.9,
                4.0, 
                4.1,
                4.2,
                4.3,
                4.4,
                4.5,
                4.6,
                4.7,
                4.8,
                4.9,
                5.0};

  float rhoTheory[] = {0.35483902825973845,
                   0.35480898330111477,
                   0.3547203350632933 ,
                   0.3545751040088587 ,
                   0.35437507372018945,
                   0.35412185325215145,
                   0.35381695577641453,
                   0.3536118445102056 ,
                   0.35365800907394396,
                   0.3538780251564116 ,
                   0.3543081857987083 ,
                   0.3548306995618629 ,
                   0.35558406406724724,
                   0.35660012666216534,
                   0.3577065383985127 ,
                   0.3584284710677728 ,
                   0.3592627246826893 ,
                   0.3602177048784971 ,
                   0.36134442288995405,
                   0.3627819146926556 ,
                   0.364411492835298  ,
                   0.36624523374207285,
                   0.3682952404864106 ,
                   0.3699815763545366 ,
                   0.3708154902401001 ,
                   0.3716215986242615 ,
                   0.37239497363280993,
                   0.3731306993072521 ,
                   0.37382386874519447,
                   0.37446958018946847,
                   0.3750629439691616 ,
                   0.37559907395838843,
                   0.3760730893093502 ,
                   0.3764801154778445 ,
                   0.3768152788224172 ,
                   0.37707371083733177,
                   0.377250546214963  ,
                   0.3773409224364039 ,
                   0.3773399795773913 ,
                   0.37724285999699436,
                   0.37704470803609647,
                   0.3767406697889577 ,
                   0.37632589294889723,
                   0.37579552676724015,
                   0.37514472176471714,
                   0.3743686294637884 ,
                   0.3734624026216854 ,
                   0.3724211949956433 ,
                   0.3712401617432428 ,
                   0.3699144572211928 ,
                   0.3684392380198456 };

  float rhoTheoryErr[] = {0.025977431244681073,
                      0.026025364720095805,
                      0.026168573817707597,
                      0.026405153353903893,
                      0.026732039352791663,
                      0.027145328832677097,
                      0.02764062686068166,
                      0.02822395015017426,
                      0.028899674243032213,
                      0.029656243195078518,
                      0.030491126440395016,
                      0.03139189685487565,
                      0.03236477542618756,
                      0.03340980600442218,
                      0.03451231581633563,
                      0.03563176671714653,
                      0.03680416420540372,
                      0.03802849077538553,
                      0.039307731551918516,
                      0.040653442494363985,
                      0.0420559371539935,
                      0.04351677261581958,
                      0.045037986230800305,
                      0.04655999493549137,
                      0.048027262510742866,
                      0.04952353248668059,
                      0.05104772135124614,
                      0.0525988449123866,
                      0.05417601185739452,
                      0.05577841788071937,
                      0.05740534606016474,
                      0.059056160467815345,
                      0.06073030545756969,
                      0.06242730518246895,
                      0.06414676036599441,
                      0.06588835027391031,
                      0.06765183192124657,
                      0.06943704033015752,
                      0.07124388905786722,
                      0.0730723708370545,
                      0.07492255846014956,
                      0.07679460572513451,
                      0.07868874846492556,
                      0.08060530565824563,
                      0.0825446805565428,
                      0.08450736181741333,
                      0.08649392462777764,
                      0.08850503171820856,
                      0.09054143488986514,
                      0.09260397422910079,
                      0.09469358132884248 };


  TGraphErrors *g_rhoTheory = new TGraphErrors();
  for(int i = 0; i < 51; i++)
  {
    g_rhoTheory->SetPoint(i,pt[i],rhoTheory[i]);
    g_rhoTheory->SetPointError(i,0.0,rhoTheoryErr[i]);
  }


  gStyle->SetOptDate(0);
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;

  const float size_marker = 1.4;
  
  float rho00_low[6] = {0.24,0.24,0.24,0.24,0.24,0.24}; // 11.5 -- 200 GeV
  float rho00_high[6] = {0.52,0.52,0.54,0.52,0.52,0.52};
  // float rho00_low[6] = {0.20,0.20,0.20,0.24,0.24,0.29}; // 11.5 -- 200 GeV
  // float rho00_high[6] = {0.50,0.52,0.64,0.40,0.40,0.35};
  float pt_low = 0.6;
  float pt_high = 4.8;
  float font_size = 0.07;
  float leg_size = 0.0475;

  string mBeanEnergy[6] = {"14.6 GeV","19.6 GeV","27 GeV","39 GeV","62.4 GeV","200 GeV"};
  int mEnergy[6] = {14,19,27,39,62,200};

  TGraphAsymmErrors *g_rho_1st_stat_temp[6];
  TGraphAsymmErrors *g_rho_1st_sys_temp[6];
  TGraphAsymmErrors *g_rho_1st_stat[6];
  TGraphAsymmErrors *g_rho_1st_sys[6];
  TGraphAsymmErrors *g_rho_2nd_stat[6];
  TGraphAsymmErrors *g_rho_2nd_sys[6];

  gStyle->SetImageScaling(3.0);


  TFile *File_Input = TFile::Open("/Users/gavinwilks/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/output/RapidityDependentPaper/rho00pt_all_Poly.root");
  for(int i_energy = 0; i_energy < 2; ++i_energy)
  {
    string GrapName_1st_stat = Form("g_rho00_order%d_%s_%s_StatError",1,vmsa::mBeamEnergy[i_energy+3].c_str(),vmsa::mPID[0].c_str());
    g_rho_1st_stat_temp[i_energy] = (TGraphAsymmErrors*)File_Input->Get(GrapName_1st_stat.c_str());
    g_rho_1st_stat[i_energy] = new TGraphAsymmErrors();
    for(int i_pt = 0; i_pt < g_rho_1st_stat_temp[i_energy]->GetN(); ++i_pt) // shift pT of 1st EP results by 0.03
    {
      double pt, rho;
      g_rho_1st_stat_temp[i_energy]->GetPoint(i_pt,pt,rho);
      double err = g_rho_1st_stat_temp[i_energy]->GetErrorYhigh(i_pt);
      g_rho_1st_stat[i_energy]->SetPoint(i_pt,pt+0.15,rho);
      g_rho_1st_stat[i_energy]->SetPointError(i_pt,0.0,0.0,err,err);
    }

    string GrapName_1st_sys = Form("g_rho00_order%d_%s_%s_SysError",1,vmsa::mBeamEnergy[i_energy+3].c_str(),vmsa::mPID[0].c_str());
    g_rho_1st_sys_temp[i_energy] = (TGraphAsymmErrors*)File_Input->Get(GrapName_1st_sys.c_str());
    g_rho_1st_sys[i_energy] = new TGraphAsymmErrors();
    for(int i_pt = 0; i_pt < g_rho_1st_sys_temp[i_energy]->GetN(); ++i_pt) // shift pT of 1st EP results by 0.03
    {
      double pt, rho;
      g_rho_1st_sys_temp[i_energy]->GetPoint(i_pt,pt,rho);
      double err = g_rho_1st_sys_temp[i_energy]->GetErrorYhigh(i_pt);
      g_rho_1st_sys[i_energy]->SetPoint(i_pt,pt+0.15,rho);
      g_rho_1st_sys[i_energy]->SetPointError(i_pt,0.0,0.0,err,err);
    }

    string GrapName_2nd_stat = Form("g_rho00_order%d_%s_%s_StatError",2,vmsa::mBeamEnergy[i_energy+3].c_str(),vmsa::mPID[0].c_str());
    g_rho_2nd_stat[i_energy] = (TGraphAsymmErrors*)File_Input->Get(GrapName_2nd_stat.c_str());
    for(int i_point = 0; i_point < g_rho_2nd_stat[i_energy]->GetN(); ++i_point)
    {
      g_rho_2nd_stat[i_energy]->SetPointEXlow(i_point,0.0);
      g_rho_2nd_stat[i_energy]->SetPointEXhigh(i_point,0.0);
    }

    string GrapName_2nd_sys = Form("g_rho00_order%d_%s_%s_SysError",2,vmsa::mBeamEnergy[i_energy+3].c_str(),vmsa::mPID[0].c_str());
    g_rho_2nd_sys[i_energy] = (TGraphAsymmErrors*)File_Input->Get(GrapName_2nd_sys.c_str());
  }

  //-----------------------------Set the Canvas---------------------------------
  const int N_x_pads = 1; //number of x-pads
  const int N_y_pads = 2; //number of y-pads
  const int N_total_pads = N_x_pads*N_y_pads; //number of total pads

  TCanvas *c_rho00_double = new TCanvas("c_rho00_double","c_rho00_double",10,10,400,667);
  c_rho00_double->SetFillColor(10);
  c_rho00_double->Divide(N_x_pads,N_y_pads,0.0,0.0,10);

  //-----------------------------Set the coordinate-------------------------------  

  TH1F* h_frame[N_total_pads];
  for(int i_pad = 0; i_pad < N_total_pads; i_pad++)
  {
    string HistName = Form("h_frame_%d",i_pad);
    h_frame[i_pad] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.1,6.0);
    for(int bin_x = 1; bin_x < h_frame[i_pad]->GetNbinsX(); bin_x++)
    {
      h_frame[i_pad]->SetBinContent(bin_x,-10.0);
    }
    h_frame[i_pad]->GetXaxis()->SetRangeUser(0.6,4.8);
    h_frame[i_pad]->GetYaxis()->SetRangeUser(0.23,0.435);
    h_frame[i_pad]->SetStats(0);
    h_frame[i_pad]->SetTitle("");
    h_frame[i_pad]->GetXaxis()->SetNdivisions(505,'N');
    h_frame[i_pad]->GetYaxis()->SetNdivisions(503,'N');
    h_frame[i_pad]->GetYaxis()->SetLabelOffset(0.0);
    h_frame[i_pad]->GetXaxis()->SetLabelOffset(0.0);
    h_frame[i_pad]->GetYaxis()->SetLabelSize(0.0);
    h_frame[i_pad]->GetXaxis()->SetLabelSize(0.0);
  }

  //-----------------------------Set the coordinate--end--------------------------  

  double Margin_t = 0.01; // top offset in percent of total size
  double Margin_b = 0.12; // bottom offset in percent of total size
  //double Margin_l = 0.24; // left offset in percent of total size
  double Margin_l = 0.19; // left offset in percent of total size
  double Margin_r = 0.01; // right offset in percent of total size
  double delta_y; // fraction in percent per pad in y direction
  double delta_x; // fraction in percent per pad in x direction
  double yUp_array[N_y_pads];
  double yDown_array[N_y_pads];
  double xLeft_array[N_x_pads];
  double xRight_array[N_x_pads];

  delta_y = (1.0-Margin_b-Margin_t)/(double)N_y_pads;
  delta_x = (1.0-Margin_l-Margin_r)/(double)N_x_pads;

  for(int y_pads = 0; y_pads < N_y_pads; y_pads++)//caculate pad size in y direction
  {
    yUp_array[y_pads] = 1.0 - Margin_t - y_pads*delta_y; 
    yDown_array[y_pads] = 1.0 - Margin_t - (y_pads+1.0)*delta_y; 
  }
  for(int x_pads = 0; x_pads < N_x_pads; x_pads++)
  {
    xLeft_array[x_pads] = Margin_l + x_pads*delta_x;
    xRight_array[x_pads] = Margin_l + (x_pads+1)*delta_x;
  }

  xLeft_array[0] = xLeft_array[0] - Margin_l;
  yDown_array[N_y_pads-1] = yDown_array[N_y_pads-1] - Margin_b;
  xRight_array[N_x_pads-1] = xRight_array[N_x_pads-1] + Margin_r;
  yUp_array[0] = yUp_array[0] + Margin_t;


  double scaling_pre = (yUp_array[N_y_pads-1]-yDown_array[N_y_pads-1])*(xRight_array[0]-xLeft_array[0]);
  //  double scaling_pre = (yUp_array[0]-yDown_array[0]);
  //  double scaling_pre = (xRight_array[0]-xLeft_array[0]);
  //    double scaling_pre = 1.0;

  for(int y_pads = 0; y_pads < N_y_pads; y_pads++)
  {
    for(int x_pads = 0; x_pads < N_x_pads; x_pads++)
    {
      int total_pad = (x_pads+1) + y_pads*N_x_pads;
      double scaling_factor = scaling_pre/((yUp_array[y_pads]-yDown_array[y_pads])*(xRight_array[x_pads]-xLeft_array[x_pads]));
      //      double scaling_factor = scaling_pre/(yUp_array[y_pads]-yDown_array[y_pads]);
      //      double scaling_factor = scaling_pre/(xRight_array[x_pads]-xLeft_array[x_pads]);
      //      cout << "scaling_factor = " << scaling_factor << endl;

      c_rho00_double->cd(total_pad);
      c_rho00_double->cd(total_pad)->SetTicks(1,1);
      c_rho00_double->cd(total_pad)->SetGrid(0,0);
      c_rho00_double->cd(total_pad)->SetFillColor(10);;
      c_rho00_double->cd(total_pad)->SetPad(xLeft_array[x_pads],yDown_array[y_pads],xRight_array[x_pads],yUp_array[y_pads]);

      cout << "xleft  = " << xLeft_array[x_pads] << endl;
      cout << "xright = " << xRight_array[x_pads] << endl;
      cout << "yup    = " << yUp_array[y_pads] << endl;
      cout << "ydown  = " << yDown_array[y_pads] << endl;

      if(x_pads == 0)
      {
	c_rho00_double->cd(total_pad)->SetLeftMargin(Margin_l/(xRight_array[x_pads]-xLeft_array[x_pads]));
	h_frame[total_pad-1]->GetYaxis()->SetLabelSize(0.065*scaling_factor);
	h_frame[total_pad-1]->GetYaxis()->SetLabelOffset(0.01*scaling_factor);
	h_frame[total_pad-1]->SetTickLength(0.03*scaling_factor,"X");
	h_frame[total_pad-1]->SetTickLength(0.03*scaling_factor,"Y");
      }
      if(x_pads == 0 && y_pads == N_y_pads-1)
      {
        c_rho00_double->cd(total_pad)->SetBottomMargin(Margin_b/(yUp_array[y_pads]-yDown_array[y_pads]));
        h_frame[total_pad-1]->SetTickLength(0.03*scaling_factor,"Y");
        h_frame[total_pad-1]->GetXaxis()->SetLabelSize(0.065*scaling_factor);
        h_frame[total_pad-1]->GetXaxis()->SetLabelOffset(0.01*scaling_factor);
        h_frame[total_pad-1]->GetYaxis()->SetLabelOffset(0.015*scaling_factor);
      }
      if(x_pads == 0)
      {
        c_rho00_double->cd(total_pad)->SetRightMargin(Margin_r/(xRight_array[x_pads]-xLeft_array[x_pads]));
        h_frame[total_pad-1]->GetXaxis()->SetLabelSize(0.065*scaling_factor);
        h_frame[total_pad-1]->GetXaxis()->SetLabelOffset(0.001*scaling_factor);
        h_frame[total_pad-1]->SetTickLength(0.03*scaling_factor,"X");
        h_frame[total_pad-1]->SetTickLength(0.03*scaling_factor,"Y");
      }
      if(x_pads == 0 && y_pads == N_y_pads-1)
      {
        c_rho00_double->cd(total_pad)->SetBottomMargin(Margin_b/(yUp_array[y_pads]-yDown_array[y_pads]));
        h_frame[total_pad-1]->SetTickLength(0.03*scaling_factor,"Y");
        h_frame[total_pad-1]->GetXaxis()->SetLabelSize(0.065*scaling_factor);
        h_frame[total_pad-1]->GetXaxis()->SetLabelOffset(0.003*scaling_factor);
      }
      if(y_pads == 0)
      {
	c_rho00_double->cd(total_pad)->SetTopMargin(Margin_t/(yUp_array[y_pads]-yDown_array[y_pads]));
      }

      if(x_pads == 0 && y_pads == 1)
      {
        // h_frame[total_pad-1]->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
        //h_frame[total_pad-1]->GetYaxis()->SetTitle("#rho_{00}");
        //h_frame[total_pad-1]->GetYaxis()->CenterTitle();
        //h_frame[total_pad-1]->GetYaxis()->SetTitleSize(0.10*scaling_factor);
      }
      h_frame[total_pad-1]->DrawCopy("PE");
      PlotLine(pt_low,pt_high,1.0/3.0,1.0/3.0,1,3,2);

      if(x_pads == 0 && y_pads != N_y_pads-1)
      {
	plotTopLegend((char*)mBeanEnergy[total_pad-1].c_str(),0.8,0.41,font_size*scaling_factor,1,0.0,42,0);
      }
      if(x_pads == 0 && y_pads == N_y_pads-1)
      {
	plotTopLegend((char*)mBeanEnergy[total_pad-1].c_str(),0.8,0.41,font_size*scaling_factor,1,0.0,42,0);
	//h_frame[total_pad-1]->SetTickLength(0.03);
        g_rhoTheory->SetLineWidth(4.0);
        g_rhoTheory->SetLineColor(kBlue);
        g_rhoTheory->SetLineStyle(2);
        g_rhoTheory->SetFillColorAlpha(kBlue,0.40);
        g_rhoTheory->SetFillStyle(3002);
        //g_rhoTheory->Draw("same C3");
        //PlotLine(0.8,1.15,0.283,0.283,kBlue,4,2);
        //string leg_line2 = "Theory";// (20-60\% &";
        //string leg_line2 = "PRL #bf{131}, 402304";// (20-60\% &";
        //plotTopLegend((char*)leg_line2.c_str(),1.2,0.28,leg_size*scaling_factor,1,0.0,42,0);
        //string leg_line2b = "1.2 < p_{T} < 5.4 GeV/c)";
        //plotTopLegend((char*)leg_line2b.c_str(),1.2,0.268,leg_size*scaling_factor,1,0.0,42,0);
      }
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_1st_stat[total_pad-1],style_phi_1st,color_phi_1st,colorDiff_phi,size_marker-0.2);
      // plotSysErrors(g_rho_1st_sys[total_pad-1],color_phi_1st);
      plotSysErrorsBox(g_rho_1st_sys[total_pad-1],color_phi_1st);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_2nd_stat[total_pad-1],style_phi_2nd,color_phi_2nd,colorDiff_phi,size_marker+0.2);
      // plotSysErrors(g_rho_2nd_sys[total_pad-1],color_phi_2nd);
      plotSysErrorsBox(g_rho_2nd_sys[total_pad-1],color_phi_2nd);
      //if(x_pads == N_x_pads-1 && y_pads == 1)
      //{
      //  plotTopLegend((char*)mBeanEnergy[total_pad-1].c_str(),3.6,0.28,0.048*scaling_factor,1,0.0,42,0);
      //}
      //if(x_pads == N_x_pads-1 && y_pads == N_y_pads-1)
      //{
      //  plotTopLegend((char*)mBeanEnergy[total_pad-1].c_str(),3.6,0.28,0.048*scaling_factor,1,0.0,42,0);
      //}
      //if(x_pads == N_x_pads-1 && y_pads == 0)
      //{
      //  plotTopLegend((char*)mBeanEnergy[total_pad-1].c_str(),3.6,0.28,0.05*scaling_factor,1,0.0,42,0);
      //}
      if(x_pads == 0 && y_pads == 0)
      {
	plotTopLegend((char*)"Au+Au (20-60\% & |y| < 1.0)",0.9,0.295,leg_size*scaling_factor,1,0.0,42,0);
	//plotTopLegend((char*)"",0.035,0.46,leg_size*scaling_factor,1,0.0,42,0);
	//plotTopLegend((char*)"Au+Au (0-80\% & |y| < 1.0)",0.1,0.51,leg_size*scaling_factor,1,0.0,42,0);

	Draw_TGAE_Point_new_Symbol(1.0,0.28,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
	plotTopLegend((char*)"#phi (1^{st}-order EP)",1.2,0.277,leg_size*scaling_factor,1,0.0,42,0);

	Draw_TGAE_Point_new_Symbol(1.0,0.265,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2);
	plotTopLegend((char*)"#phi (2^{nd}-order EP)",1.2,0.262,leg_size*scaling_factor,1,0.0,42,0);
      }

      // if(x_pads == 0 && y_pads == N_y_pads-1) plotTopLegend("#font[12]{p}_{#font[132]{T}}",0.94,0.10,0.1,1,0.0,42,1);
      // if(x_pads == 1 && y_pads == N_y_pads-1) plotTopLegend("#font[132]{(GeV}/#font[12]{c}#font[132]{)}",0.0,0.09,0.09,1,0.0,42,1);
      if(x_pads == 0 && y_pads == 0) plotTopLegend((char*)"#rho_{00}",0.05,0.0,0.1,1,90,42,1);
      if(x_pads == 0 && y_pads == N_y_pads-1) plotTopLegend((char*)"p_{T} (GeV/c)",0.43,0.05,0.075,1,0.0,42,1);
      //if(x_pads == 1 && y_pads == N_y_pads-1) plotTopLegend((char*)"(GeV/c)",0.0,0.09,0.09,1,0.0,42,1);
    }
  }

  c_rho00_double->SaveAs("figures/c_rhoSys_pt_BESII.pdf");
  c_rho00_double->SaveAs("figures/c_rhoSys_pt_BESII.png");
}

//void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color)
//{
//  for(int i_pt = 0; i_pt < g_rho->GetN(); ++i_pt) // plot sys errors
//  {
//    double pt, rho;
//    g_rho->GetPoint(i_pt,pt,rho);
//    double err = g_rho->GetErrorYhigh(i_pt);
//
//    PlotLine(pt-0.02,pt+0.02,rho+err,rho+err,plot_color,2,1);
//    PlotLine(pt-0.02,pt-0.02,rho+err-0.005,rho+err,plot_color,2,1);
//    PlotLine(pt+0.02,pt+0.02,rho+err-0.005,rho+err,plot_color,2,1);
//    PlotLine(pt-0.02,pt+0.02,rho-err,rho-err,plot_color,2,1);
//    PlotLine(pt-0.02,pt-0.02,rho-err+0.005,rho-err,plot_color,2,1);
//    PlotLine(pt+0.02,pt+0.02,rho-err+0.005,rho-err,plot_color,2,1);
//  }
//}

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color)
{
  const int nPt = g_rho->GetN();
  TBox *bSys[nPt];
  for(int i_pt = 0; i_pt < g_rho->GetN(); ++i_pt) // plot sys errors
  {
    double pt, rho;
    g_rho->GetPoint(i_pt,pt,rho);
    double err = g_rho->GetErrorYhigh(i_pt);

    bSys[i_pt] = new TBox(pt-0.12,rho-err,pt+0.12,rho+err);
    bSys[i_pt]->SetFillColor(0);
    bSys[i_pt]->SetFillStyle(0);
    bSys[i_pt]->SetLineStyle(1);
    bSys[i_pt]->SetLineWidth(1);
    bSys[i_pt]->SetLineColor(plot_color);
    bSys[i_pt]->Draw("l Same");
  }
}
