#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../../Utility/draw.h"
// #include "../../Utility/StSpinAlignmentCons.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void plotExtFig8_rho00Cent()
{
  gStyle->SetOptDate(0);
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;

  const int style_Kstr = 20;
  const int color_Kstr = kAzure-9;

  const float size_marker = 1.2;

  // read in Phi Data for 27GeV (run18), 39GeV, 62.4GeV and 200GeV (run14)
  TFile *File_InPutPhi[4];
  std::string const mBeamEnergyPhi[4] = {"27GeV_Run18","39GeV","62GeV","200GeV_Run14"};
  // std::string const mLegEnergyPhi[4] = {"Au+Au 27 GeV (Run18)","Au+Au 39 GeV","Au+Au 62.4 GeV","Au+Au 200 GeV (Run14)"};
  std::string const mLegEnergyPhi[4] = {"Au+Au 27 GeV","Au+Au 39 GeV","Au+Au 62.4 GeV","Au+Au 200 GeV"};
  TGraphAsymmErrors *g_rhoPhiCent_1st_stat[4];
  TGraphAsymmErrors *g_rhoPhiCent_1st_sys[4];
  TGraphAsymmErrors *g_rhoPhiCent_2nd_stat[4];
  TGraphAsymmErrors *g_rhoPhiCent_2nd_sys[4];

  for(int i_energy = 0; i_energy < 4; ++i_energy)
  {
    string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rhoCent_%s_LXaxis.root",mBeamEnergyPhi[i_energy].c_str());
    File_InPutPhi[i_energy] = TFile::Open(inputfile.c_str());

    string GraName_1st_stat = Form("g_rhoCent_%s_1st_stat",mBeamEnergyPhi[i_energy].c_str());
    g_rhoPhiCent_1st_stat[i_energy] = (TGraphAsymmErrors*)File_InPutPhi[i_energy]->Get(GraName_1st_stat.c_str());

    string GraName_1st_sys = Form("g_rhoCent_%s_1st_sys",mBeamEnergyPhi[i_energy].c_str());
    g_rhoPhiCent_1st_sys[i_energy] = (TGraphAsymmErrors*)File_InPutPhi[i_energy]->Get(GraName_1st_sys.c_str());

    string GraName_2nd_stat = Form("g_rhoCent_%s_2nd_stat",mBeamEnergyPhi[i_energy].c_str());
    g_rhoPhiCent_2nd_stat[i_energy] = (TGraphAsymmErrors*)File_InPutPhi[i_energy]->Get(GraName_2nd_stat.c_str());

    string GraName_2nd_sys = Form("g_rhoCent_%s_2nd_sys",mBeamEnergyPhi[i_energy].c_str());
    g_rhoPhiCent_2nd_sys[i_energy] = (TGraphAsymmErrors*)File_InPutPhi[i_energy]->Get(GraName_2nd_sys.c_str());
  }

  TFile *File_InPutKstar;
  std::string const mBeamEnergyKstar[3] = {"39GeV","54GeV","200GeV"};
  std::string const mLegEnergyKstar[3] = {"Au+Au 39 GeV","Au+Au 54.4 GeV","Au+Au 200 GeV"};
  TGraphAsymmErrors *g_rhoKstarCent_2nd_stat[3];
  TGraphAsymmErrors *g_rhoKstarCent_2nd_sys[3];
  for(int i_energy = 0; i_energy < 3; ++i_energy)
  {
    string inputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_Cent_Nov22_2021.root";
    File_InPutKstar = TFile::Open(inputfile.c_str());

    string GraName_2nd_stat = Form("g_rhoCent_%s_2nd_stat",mBeamEnergyKstar[i_energy].c_str());
    g_rhoKstarCent_2nd_stat[i_energy] = (TGraphAsymmErrors*)File_InPutKstar->Get(GraName_2nd_stat.c_str());

    string GraName_2nd_sys = Form("g_rhoCent_%s_2nd_sys",mBeamEnergyKstar[i_energy].c_str());
    g_rhoKstarCent_2nd_sys[i_energy] = (TGraphAsymmErrors*)File_InPutKstar->Get(GraName_2nd_sys.c_str());
  }

  // Prepare Canvas
  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,1600,800);
  c_play->cd();

  const double width_left = 0.030;
  const double width_bottom = 0.04;
  const double delta_x = 0.26;
  const double delta_y = 0.50;
  const double length_x = delta_x + delta_x*0.9*3.0;
  const double length_y = delta_y + delta_y*0.9;

  const double cornerLL_x = width_left;
  const double cornerLL_y = width_bottom;
  const double cornerUR_x = width_left+length_x;
  const double cornerUR_y = width_bottom+length_y;

  // test size will use the shorter one in pad height or width
  // label size will use the pad size?
  const double size_font = 0.08;
  const double scale_top0 = (delta_y*0.9*800.0)/(delta_x*0.9*1600.0);
  const double scale_bottom0 = 0.9;
  const double size_offset = 0.02;

  // Margin Pad
  TPad *pad_LeftMargin   = new TPad("pad_LeftMargin","pad_LeftMargin",0.00,0.00,cornerLL_x,1.00,0);
  TPad *pad_RightMargin  = new TPad("pad_RightMargin","pad_RightMargin",cornerUR_x,0.00,1.0,1.0,0);
  TPad *pad_TopMargin    = new TPad("pad_TopMargin","pad_TopMargin",cornerLL_x,cornerUR_y,cornerUR_x,1.00,0);
  TPad *pad_BottomMargin = new TPad("pad_BottomMargin","pad_BottomMargin",cornerLL_x,0.00,cornerUR_x,cornerLL_y,0);
  pad_LeftMargin->Draw();
  pad_RightMargin->Draw();
  pad_TopMargin->Draw();
  pad_BottomMargin->Draw();

  TPad *pad_Bottom0 = new TPad("pad_Bottom0","pad_Bottom0",cornerLL_x+0.0*delta_x, cornerLL_y, cornerLL_x+1.0*delta_x, cornerLL_y+delta_y,0);
  TPad *pad_Bottom1 = new TPad("pad_Bottom1","pad_Bottom1",cornerLL_x+1.0*delta_x, cornerLL_y, cornerLL_x+1.9*delta_x, cornerLL_y+delta_y,0);
  TPad *pad_Bottom2 = new TPad("pad_Bottom2","pad_Bottom2",cornerLL_x+1.9*delta_x, cornerLL_y, cornerLL_x+2.8*delta_x, cornerLL_y+delta_y,0);
  TPad *pad_Bottom3 = new TPad("pad_Bottom3","pad_Bottom3",cornerLL_x+2.8*delta_x, cornerLL_y, cornerLL_x+3.7*delta_x, cornerLL_y+delta_y,0);
  pad_Bottom0->Draw();
  pad_Bottom1->Draw();
  pad_Bottom2->Draw();
  pad_Bottom3->Draw();

  TPad *pad_Top0 = new TPad("pad_Top0","pad_Top0",cornerLL_x+0.0*delta_x, cornerLL_y+delta_y, cornerLL_x+1.0*delta_x, cornerLL_y+1.9*delta_y,0);
  TPad *pad_Top1 = new TPad("pad_Top1","pad_Top1",cornerLL_x+1.0*delta_x, cornerLL_y+delta_y, cornerLL_x+1.9*delta_x, cornerLL_y+1.9*delta_y,0);
  TPad *pad_Top2 = new TPad("pad_Top2","pad_Top2",cornerLL_x+1.9*delta_x, cornerLL_y+delta_y, cornerLL_x+2.8*delta_x, cornerLL_y+1.9*delta_y,0);
  TPad *pad_Top3 = new TPad("pad_Top3","pad_Top3",cornerLL_x+2.8*delta_x, cornerLL_y+delta_y, cornerLL_x+3.7*delta_x, cornerLL_y+1.9*delta_y,0);
  pad_Top0->Draw();
  pad_Top1->Draw();
  pad_Top2->Draw();
  pad_Top3->Draw();

  const float cent_start = -3.5;
  const float cent_stop  = 89.5;
  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-5.5,94.5);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(cent_start,cent_stop);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(size_font);
  h_frame->GetXaxis()->SetLabelOffset(0.01);

  h_frame->GetYaxis()->SetRangeUser(0.22,0.44);
  h_frame->GetYaxis()->SetNdivisions(504,'N');
  h_frame->GetYaxis()->SetLabelSize(size_font);

  // cout << "Bottom0 Pad Width: " << pad_Bottom0->XtoPixel(pad_Bottom0->GetX2()) << endl;
  // cout << "Bottom0 Pad Hight : " << pad_Bottom0->YtoPixel(pad_Bottom0->GetY1()) << endl;
  // cout << "Bottom1 Pad Width: " << pad_Bottom1->XtoPixel(pad_Bottom1->GetX2()) << endl;
  // cout << "Bottom1 Pad Hight : " << pad_Bottom1->YtoPixel(pad_Bottom1->GetY1()) << endl;
  // cout << "Bottom2 Pad Width: " << pad_Bottom2->XtoPixel(pad_Bottom2->GetX2()) << endl;
  // cout << "Bottom2 Pad Hight : " << pad_Bottom2->YtoPixel(pad_Bottom2->GetY1()) << endl;
  // cout << "Bottom3 Pad Width: " << pad_Bottom3->XtoPixel(pad_Bottom3->GetX2()) << endl;
  // cout << "Bottom3 Pad Hight : " << pad_Bottom3->YtoPixel(pad_Bottom3->GetY1()) << endl;

  // cout << "Top0 Pad Width: " << pad_Top0->XtoPixel(pad_Top0->GetX2()) << endl;
  // cout << "Top0 Pad Hight : " << pad_Top0->YtoPixel(pad_Top0->GetY1()) << endl;
  // cout << "Top1 Pad Width: " << pad_Top1->XtoPixel(pad_Top1->GetX2()) << endl;
  // cout << "Top1 Pad Hight : " << pad_Top1->YtoPixel(pad_Top1->GetY1()) << endl;
  // cout << "Top2 Pad Width: " << pad_Top2->XtoPixel(pad_Top2->GetX2()) << endl;
  // cout << "Top2 Pad Hight : " << pad_Top2->YtoPixel(pad_Top2->GetY1()) << endl;
  // cout << "Top3 Pad Width: " << pad_Top3->XtoPixel(pad_Top3->GetX2()) << endl;
  // cout << "Top3 Pad Hight : " << pad_Top3->YtoPixel(pad_Top3->GetY1()) << endl;

  // Phi
  pad_Top0->cd();
  pad_Top0->cd()->SetLeftMargin(0.1);
  pad_Top0->cd()->SetRightMargin(0.0);
  pad_Top0->cd()->SetTopMargin(0.0);
  pad_Top0->cd()->SetBottomMargin(0.0);
  pad_Top0->cd()->SetTicks(1,1);
  pad_Top0->cd()->SetGrid(0,0);
  h_frame->GetYaxis()->SetLabelSize(size_font*scale_top0);
  h_frame->DrawCopy("pE");
  PlotLine(cent_start,cent_stop,1.0/3.0,1.0/3.0,1,3,2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhiCent_1st_stat[0],style_phi_1st,color_phi_1st,size_marker-0.2);
  // plotSysErrors(g_rhoPhiCent_1st_sys[0],color_phi_1st);
  plotSysErrorsBox(g_rhoPhiCent_1st_sys[0],color_phi_1st);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhiCent_2nd_stat[0],style_phi_2nd,color_phi_2nd,size_marker+0.2);
  // plotSysErrors(g_rhoPhiCent_2nd_sys[0],color_phi_2nd);
  plotSysErrorsBox(g_rhoPhiCent_2nd_sys[0],color_phi_2nd);
  plotTopLegend((char*)"a) #phi",5,0.415,size_font*scale_top0,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(25,0.4195,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
  plotTopLegend((char*)"1^{st}-order EP",28,0.415,size_font*scale_top0,1,0.0,42,0);
  plotTopLegend((char*)mLegEnergyPhi[0].c_str(),30,0.25,size_font*scale_top0*1.05,1,0.0,42,0);
  pad_Top0->Update();

  pad_Top1->cd();
  pad_Top1->cd()->SetLeftMargin(0.0);
  pad_Top1->cd()->SetRightMargin(0.0);
  pad_Top1->cd()->SetTopMargin(0.0);
  pad_Top1->cd()->SetBottomMargin(0.0);
  pad_Top1->cd()->SetTicks(1,1);
  pad_Top1->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  PlotLine(cent_start,cent_stop,1.0/3.0,1.0/3.0,1,3,2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhiCent_1st_stat[1],style_phi_1st,color_phi_1st,size_marker-0.2);
  // plotSysErrors(g_rhoPhiCent_1st_sys[1],color_phi_1st);
  plotSysErrorsBox(g_rhoPhiCent_1st_sys[1],color_phi_1st);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhiCent_2nd_stat[1],style_phi_2nd,color_phi_2nd,size_marker+0.2);
  // plotSysErrors(g_rhoPhiCent_2nd_sys[1],color_phi_2nd);
  plotSysErrorsBox(g_rhoPhiCent_2nd_sys[1],color_phi_2nd);
  plotTopLegend((char*)"b) #phi",5,0.415,size_font,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(25,0.4195,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  plotTopLegend((char*)"2^{nd}-order EP",28,0.415,size_font*scale_top0,1,0.0,42,0);
  plotTopLegend((char*)mLegEnergyPhi[1].c_str(),30,0.25,size_font*1.05,1,0.0,42,0);
  pad_Top1->Update();

  pad_Top2->cd();
  pad_Top2->cd()->SetLeftMargin(0.0);
  pad_Top2->cd()->SetRightMargin(0.0);
  pad_Top2->cd()->SetTopMargin(0.0);
  pad_Top2->cd()->SetBottomMargin(0.0);
  pad_Top2->cd()->SetTicks(1,1);
  pad_Top2->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  PlotLine(cent_start,cent_stop,1.0/3.0,1.0/3.0,1,3,2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhiCent_1st_stat[2],style_phi_1st,color_phi_1st,size_marker-0.2);
  // plotSysErrors(g_rhoPhiCent_1st_sys[2],color_phi_1st);
  plotSysErrorsBox(g_rhoPhiCent_1st_sys[2],color_phi_1st);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhiCent_2nd_stat[2],style_phi_2nd,color_phi_2nd,size_marker+0.2);
  // plotSysErrors(g_rhoPhiCent_2nd_sys[2],color_phi_2nd);
  plotSysErrorsBox(g_rhoPhiCent_2nd_sys[2],color_phi_2nd);
  plotTopLegend((char*)"c) #phi",5,0.415,size_font,1,0.0,42,0);
  plotTopLegend((char*)mLegEnergyPhi[2].c_str(),30,0.25,size_font*1.05,1,0.0,42,0);
  pad_Top2->Update();

  pad_Top3->cd();
  pad_Top3->cd()->SetLeftMargin(0.0);
  pad_Top3->cd()->SetRightMargin(0.0);
  pad_Top3->cd()->SetTopMargin(0.0);
  pad_Top3->cd()->SetBottomMargin(0.0);
  pad_Top3->cd()->SetTicks(1,1);
  pad_Top3->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  PlotLine(cent_start,cent_stop,1.0/3.0,1.0/3.0,1,3,2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhiCent_1st_stat[3],style_phi_1st,color_phi_1st,size_marker-0.2);
  // plotSysErrors(g_rhoPhiCent_1st_sys[3],color_phi_1st);
  plotSysErrorsBox(g_rhoPhiCent_1st_sys[3],color_phi_1st);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhiCent_2nd_stat[3],style_phi_2nd,color_phi_2nd,size_marker+0.2);
  // plotSysErrors(g_rhoPhiCent_2nd_sys[3],color_phi_2nd);
  plotSysErrorsBox(g_rhoPhiCent_2nd_sys[3],color_phi_2nd);
  plotTopLegend((char*)"d) #phi",5,0.415,size_font,1,0.0,42,0);
  plotTopLegend((char*)mLegEnergyPhi[3].c_str(),30,0.25,size_font*1.05,1,0.0,42,0);
  pad_Top3->Update();

  // K*
  pad_Bottom0->cd();
  pad_Bottom0->cd()->SetLeftMargin(0.1);
  pad_Bottom0->cd()->SetRightMargin(0.0);
  pad_Bottom0->cd()->SetTopMargin(0.0);
  pad_Bottom0->cd()->SetBottomMargin(0.10);
  pad_Bottom0->cd()->SetTicks(1,1);
  pad_Bottom0->cd()->SetGrid(0,0);
  h_frame->GetXaxis()->SetLabelSize(size_font*scale_bottom0);
  h_frame->GetYaxis()->SetLabelSize(size_font*scale_bottom0);
  h_frame->GetXaxis()->SetLabelOffset(size_offset/scale_bottom0/scale_bottom0);
  // h_frame->GetYaxis()->SetRangeUser(0.17,0.39);
  h_frame->DrawCopy("pE");
  PlotLine(cent_start,cent_stop,1.0/3.0,1.0/3.0,1,3,2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoKstarCent_2nd_stat[0],style_Kstr,color_Kstr,size_marker);
  // plotSysErrors(g_rhoKstarCent_2nd_sys[0],color_Kstr+2);
  plotSysErrorsBox(g_rhoKstarCent_2nd_sys[0],color_Kstr+2);
  plotTopLegend((char*)"e) K^{*0}",5,0.415,size_font*scale_bottom0,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(25,0.420,0.0,0.0,0.0,0.0,style_Kstr,color_Kstr,size_marker);
  plotTopLegend((char*)"2^{nd}-order EP",28,0.415,size_font*scale_top0,1,0.0,42,0);
  plotTopLegend((char*)mLegEnergyKstar[0].c_str(),30,0.25,size_font*scale_bottom0*1.05,1,0.0,42,0);
  pad_Bottom0->Update();

  pad_Bottom1->cd();
  pad_Bottom1->cd()->SetLeftMargin(0.0);
  pad_Bottom1->cd()->SetRightMargin(0.0);
  pad_Bottom1->cd()->SetTopMargin(0.0);
  pad_Bottom1->cd()->SetBottomMargin(0.10);
  pad_Bottom1->cd()->SetTicks(1,1);
  pad_Bottom1->cd()->SetGrid(0,0);
  h_frame->GetXaxis()->SetLabelSize(size_font);
  h_frame->GetXaxis()->SetLabelOffset(size_offset);
  h_frame->DrawCopy("pE");
  PlotLine(cent_start,cent_stop,1.0/3.0,1.0/3.0,1,3,2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoKstarCent_2nd_stat[1],style_Kstr,color_Kstr,size_marker);
  // plotSysErrors(g_rhoKstarCent_2nd_sys[1],color_Kstr+2);
  plotSysErrorsBox(g_rhoKstarCent_2nd_sys[1],color_Kstr+2);
  plotTopLegend((char*)"f) K^{*0}",5,0.415,size_font,1,0.0,42,0);
  plotTopLegend((char*)mLegEnergyKstar[1].c_str(),30,0.25,size_font*1.05,1,0.0,42,0);
  pad_Bottom1->Update();

  pad_Bottom2->cd();
  pad_Bottom2->cd()->SetLeftMargin(0.0);
  pad_Bottom2->cd()->SetRightMargin(0.0);
  pad_Bottom2->cd()->SetTopMargin(0.0);
  pad_Bottom2->cd()->SetBottomMargin(0.10);
  pad_Bottom2->cd()->SetTicks(1,1);
  pad_Bottom2->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  PlotLine(cent_start,cent_stop,1.0/3.0,1.0/3.0,1,3,2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoKstarCent_2nd_stat[2],style_Kstr,color_Kstr,size_marker);
  // plotSysErrors(g_rhoKstarCent_2nd_sys[2],color_Kstr+2);
  plotSysErrorsBox(g_rhoKstarCent_2nd_sys[2],color_Kstr+2);
  plotTopLegend((char*)"g) K^{*0}",5,0.415,size_font,1,0.0,42,0);
  plotTopLegend((char*)mLegEnergyKstar[2].c_str(),30,0.25,size_font*1.05,1,0.0,42,0);
  pad_Bottom2->Update();

  pad_Bottom3->cd();
  pad_Bottom3->cd()->SetLeftMargin(0.0);
  pad_Bottom3->cd()->SetRightMargin(0.0);
  pad_Bottom3->cd()->SetTopMargin(0.0);
  pad_Bottom3->cd()->SetBottomMargin(0.10);
  pad_Bottom3->cd()->SetTicks(1,1);
  pad_Bottom3->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  plotTopLegend((char*)"h)",5,0.415,size_font,1,0.0,42,0);

  // Phi Legend
  // Draw_TGAE_Point_new_Symbol(30,0.30,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
  // Draw_TGAE_Point_new_Symbol(35,0.30,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  plotTopLegend((char*)"#phi",10,0.38,size_font,1,0.0,42,0);
  plotTopLegend((char*)"|y| < 1.0",20,0.38,size_font,1,0.0,42,0);
  plotTopLegend((char*)"1.2 < p_{T}< 5.4 GeV/c",20,0.36,size_font,1,0.0,42,0);

  // plotTopLegend((char*)"#phi  (|y| < 1.0 & 1.2 < p_{T}< 5.4 GeV/c)",8,0.2775,size_font,1,0.0,42,0);
  // plotTopLegend((char*)"",3,0.26,size_font,1,0.0,42,0);

  // K* Legend
  // Draw_TGAE_Point_new_Symbol(5,0.26,0.0,0.0,0.0,0.0,style_Kstr,color_Kstr,size_marker);
  plotTopLegend((char*)"K^{*0}",10,0.33,size_font,1,0.0,42,0);
  plotTopLegend((char*)"|y| < 1.0",20,0.33,size_font,1,0.0,42,0);
  plotTopLegend((char*)"1.0 < p_{T}< 5.0 GeV/c",20,0.31,size_font,1,0.0,42,0);
  pad_Bottom3->Update();

  pad_LeftMargin->cd();
  plotTopLegend((char*)"#rho_{00}",0.5,0.51,0.8,1,90.0,42,1);
  pad_LeftMargin->Update();

  pad_BottomMargin->cd();
  plotTopLegend((char*)"Centrality (%)",0.46,0.4,0.8,1,0.0,42,1);
  pad_BottomMargin->Update();

  c_play->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/NatureSubmission/extFig8_rho00Centrality.eps");
  c_play->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/NatureSubmission/extFig8_rho00Centrality.png");
}

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color)
{
  for(int i_cent = 0; i_cent < g_rho->GetN(); ++i_cent) // plot sys errors
  {
    double centrality, rho;
    g_rho->GetPoint(i_cent,centrality,rho);
    double err = g_rho->GetErrorYhigh(i_cent);

    PlotLine(centrality-1.5,centrality+1.5,rho+err,rho+err,plot_color,1,1);
    PlotLine(centrality-1.5,centrality-1.5,rho+err-0.002,rho+err,plot_color,1,1);
    PlotLine(centrality+1.5,centrality+1.5,rho+err-0.002,rho+err,plot_color,1,1);
    PlotLine(centrality-1.5,centrality+1.5,rho-err,rho-err,plot_color,1,1);
    PlotLine(centrality-1.5,centrality-1.5,rho-err+0.002,rho-err,plot_color,1,1);
    PlotLine(centrality+1.5,centrality+1.5,rho-err+0.002,rho-err,plot_color,1,1);
  }
}

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color)
{
  const int nCent = g_rho->GetN();
  TBox *bSys[nCent];
  for(int i_cent = 0; i_cent < g_rho->GetN(); ++i_cent) // plot sys errors
  {
    double centrality, rho;
    g_rho->GetPoint(i_cent,centrality,rho);
    double err = g_rho->GetErrorYhigh(i_cent);

    bSys[i_cent] = new TBox(centrality-1.8,rho-err,centrality+1.8,rho+err);
    bSys[i_cent]->SetFillColor(0);
    bSys[i_cent]->SetFillStyle(0);
    bSys[i_cent]->SetLineStyle(1);
    bSys[i_cent]->SetLineWidth(1);
    bSys[i_cent]->SetLineColor(plot_color);
    bSys[i_cent]->Draw("l Same");
  }
}
