#include <TCanvas.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include "../../../Utility/draw.h"

using namespace std;

void plotCanvas()
{
  gStyle->SetOptDate(0);
  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,1600,800);
  c_play->cd();

  const double width_left = 0.020;
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
  const double size_font = 0.06;
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

  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-1.5,98.5);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(-1.0,89.5);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(size_font);
  h_frame->GetXaxis()->SetLabelOffset(0.01);

  h_frame->GetYaxis()->SetRangeUser(0.14,0.42);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetLabelSize(size_font);

  cout << "Bottom0 Pad Width: " << pad_Bottom0->XtoPixel(pad_Bottom0->GetX2()) << endl;
  cout << "Bottom0 Pad Hight : " << pad_Bottom0->YtoPixel(pad_Bottom0->GetY1()) << endl;
  cout << "Bottom1 Pad Width: " << pad_Bottom1->XtoPixel(pad_Bottom1->GetX2()) << endl;
  cout << "Bottom1 Pad Hight : " << pad_Bottom1->YtoPixel(pad_Bottom1->GetY1()) << endl;
  cout << "Bottom2 Pad Width: " << pad_Bottom2->XtoPixel(pad_Bottom2->GetX2()) << endl;
  cout << "Bottom2 Pad Hight : " << pad_Bottom2->YtoPixel(pad_Bottom2->GetY1()) << endl;
  cout << "Bottom3 Pad Width: " << pad_Bottom3->XtoPixel(pad_Bottom3->GetX2()) << endl;
  cout << "Bottom3 Pad Hight : " << pad_Bottom3->YtoPixel(pad_Bottom3->GetY1()) << endl;

  cout << "Top0 Pad Width: " << pad_Top0->XtoPixel(pad_Top0->GetX2()) << endl;
  cout << "Top0 Pad Hight : " << pad_Top0->YtoPixel(pad_Top0->GetY1()) << endl;
  cout << "Top1 Pad Width: " << pad_Top1->XtoPixel(pad_Top1->GetX2()) << endl;
  cout << "Top1 Pad Hight : " << pad_Top1->YtoPixel(pad_Top1->GetY1()) << endl;
  cout << "Top2 Pad Width: " << pad_Top2->XtoPixel(pad_Top2->GetX2()) << endl;
  cout << "Top2 Pad Hight : " << pad_Top2->YtoPixel(pad_Top2->GetY1()) << endl;
  cout << "Top3 Pad Width: " << pad_Top3->XtoPixel(pad_Top3->GetX2()) << endl;
  cout << "Top3 Pad Hight : " << pad_Top3->YtoPixel(pad_Top3->GetY1()) << endl;

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
  h_frame->DrawCopy("pE");
  PlotLine(0.0,100.0,1.0/3.0,1.0/3.0,1,3,2);
  plotTopLegend((char*)"This is a test",0.15,0.35,size_font*scale_bottom0,1,0.0,42,1);
  pad_Bottom0->Update();
  // pad_Bottom0->cd()->SetFillColor(2);
  // pad_Bottom0->cd()->SetBorderMode(1);
  // pad_Bottom0->cd()->SetFrameFillColor(10);

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
  PlotLine(0.0,100.0,1.0/3.0,1.0/3.0,1,3,2);
  plotTopLegend((char*)"This is a test",0.15,0.35,size_font,1,0.0,42,1);
  pad_Bottom1->Update();

  pad_Bottom2->cd();
  pad_Bottom2->cd()->SetLeftMargin(0.0);
  pad_Bottom2->cd()->SetRightMargin(0.0);
  pad_Bottom2->cd()->SetTopMargin(0.0);
  pad_Bottom2->cd()->SetBottomMargin(0.10);
  pad_Bottom2->cd()->SetTicks(1,1);
  pad_Bottom2->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  PlotLine(0.0,100.0,1.0/3.0,1.0/3.0,1,3,2);
  plotTopLegend((char*)"This is a test",0.15,0.35,size_font,1,0.0,42,1);
  pad_Bottom2->Update();

  pad_Bottom3->cd();
  pad_Bottom3->cd()->SetLeftMargin(0.0);
  pad_Bottom3->cd()->SetRightMargin(0.0);
  pad_Bottom3->cd()->SetTopMargin(0.0);
  pad_Bottom3->cd()->SetBottomMargin(0.10);
  pad_Bottom3->cd()->SetTicks(1,1);
  pad_Bottom3->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  PlotLine(0.0,100.0,1.0/3.0,1.0/3.0,1,3,2);
  plotTopLegend((char*)"This is a test",0.15,0.35,size_font,1,0.0,42,1);
  pad_Bottom3->Update();

  pad_Top0->cd();
  pad_Top0->cd()->SetLeftMargin(0.1);
  pad_Top0->cd()->SetRightMargin(0.0);
  pad_Top0->cd()->SetTopMargin(0.0);
  pad_Top0->cd()->SetBottomMargin(0.0);
  pad_Top0->cd()->SetTicks(1,1);
  pad_Top0->cd()->SetGrid(0,0);
  h_frame->GetYaxis()->SetLabelSize(size_font*scale_top0);
  h_frame->DrawCopy("pE");
  PlotLine(0.0,100.0,1.0/3.0,1.0/3.0,1,3,2);
  plotTopLegend((char*)"This is a test",0.15,0.35,size_font*scale_top0,1,0.0,42,1);
  pad_Top0->Update();

  pad_Top1->cd();
  pad_Top1->cd()->SetLeftMargin(0.0);
  pad_Top1->cd()->SetRightMargin(0.0);
  pad_Top1->cd()->SetTopMargin(0.0);
  pad_Top1->cd()->SetBottomMargin(0.0);
  pad_Top1->cd()->SetTicks(1,1);
  pad_Top1->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  PlotLine(0.0,100.0,1.0/3.0,1.0/3.0,1,3,2);
  plotTopLegend((char*)"This is a test",0.15,0.35,size_font,1,0.0,42,1);
  pad_Top1->Update();

  pad_Top2->cd();
  pad_Top2->cd()->SetLeftMargin(0.0);
  pad_Top2->cd()->SetRightMargin(0.0);
  pad_Top2->cd()->SetTopMargin(0.0);
  pad_Top2->cd()->SetBottomMargin(0.0);
  pad_Top2->cd()->SetTicks(1,1);
  pad_Top2->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  PlotLine(0.0,100.0,1.0/3.0,1.0/3.0,1,3,2);
  plotTopLegend((char*)"This is a test",0.15,0.35,size_font,1,0.0,42,1);
  pad_Top2->Update();

  pad_Top3->cd();
  pad_Top3->cd()->SetLeftMargin(0.0);
  pad_Top3->cd()->SetRightMargin(0.0);
  pad_Top3->cd()->SetTopMargin(0.0);
  pad_Top3->cd()->SetBottomMargin(0.0);
  pad_Top3->cd()->SetTicks(1,1);
  pad_Top3->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  PlotLine(0.0,100.0,1.0/3.0,1.0/3.0,1,3,2);
  plotTopLegend((char*)"This is a test",0.15,0.35,size_font,1,0.0,42,1);
  pad_Top3->Update();

  // c_play->SaveAs("./testCanvas.eps");
}
