

#include "TString.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGaxis.h"
//#include <iostream.h>
#include <fstream>
#include "TMath.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TExec.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TPolyLine.h"
#include "TMarker3DBox.h"
#include "TVector3.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "Math/MinimizerOptions.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLCamera.h"
#include "TGLPerspectiveCamera.h"
#include "TGFrame.h"
#include "TGLUtil.h"
#include "TGLLightSet.h"
#include "TGLCameraOverlay.h"
#include "TLorentzVector.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <vector>
#include "TF1.h"

#include "TEveManager.h"
#include "TEveGeoShape.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TEveTrans.h"

static TGeoVolume   *top;
static TGeoManager  *geom       = new TGeoManager("geom","My 3D Project");
static TGeoMaterial *geo_track  = new TGeoMaterial("geo_track",55.84,26.7,7.87);
static TGeoMedium   *med_track  = new TGeoMedium("med_track",1,geo_track);

//----------------------------------------------------------------------------------------
void Draw_Track()
{
    geo_track->SetTransparency(80); // higher value means more transparent, 100 is maximum

    TGeoVolume *track       = geom->MakeTube("track",med_track,5,6.0,600.0);  // r_min, r_max, dz (half of total length)
    track       ->SetLineColor(2);
    top->AddNodeOverlap(track,1,new TGeoTranslation(0,0,0));
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_STAR_3D()
{

    //TGeoManager *geom = new TGeoManager("geom","My 3D Project");
    //------------------Creat materials------------------------------
    TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
    vacuum->SetTransparency(0);
    TGeoMaterial *Fe = new TGeoMaterial("Fe",55.84,26.7,7.87);
    Fe->SetTransparency(90); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_outer_tube = new TGeoMaterial("M_outer_tube",55.84,26.7,7.87);
    //M_outer_tube->SetTransparency(93); // higher value means more transparent, 100 is maximum
    M_outer_tube->SetTransparency(90); // higher value means more transparent, 100 is maximum

    // TGeoMaterial *M_outer_cap = new TGeoMaterial("M_outer_cap",55.84,26.7,7.87);
    TGeoMaterial *M_outer_cap = new TGeoMaterial("M_outer_cap",55.84,26.7,7.87);
    M_outer_cap->SetTransparency(70); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_IDS = new TGeoMaterial("M_IDS",55.84,26.7,7.87);
    M_IDS       ->SetTransparency(80); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_beampipe = new TGeoMaterial("M_beampipe",55.84,26.7,7.87);
    M_beampipe       ->SetTransparency(70); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_Pixel_support = new TGeoMaterial("M_Pixel_support",55.84,26.7,7.87);
    M_Pixel_support    ->SetTransparency(70); // higher value means more transparent, 100 is maximum


    //------------------Create media---------------------------------
    TGeoMedium *Air = new TGeoMedium("Air",0,vacuum);
    TGeoMedium *Iron = new TGeoMedium("Iron",1,Fe);
    TGeoMedium *Me_outer_tube = new TGeoMedium("Me_outer_tube",1,M_outer_tube);
    TGeoMedium *Me_outer_cap = new TGeoMedium("Me_outer_cap",1,M_outer_cap);
    TGeoMedium *Me_IDS        = new TGeoMedium("Me_IDS",1,M_IDS);
    TGeoMedium *Me_beampipe   = new TGeoMedium("Me_beampipe",1,M_beampipe);
    TGeoMedium *Me_Pixel_support   = new TGeoMedium("Me_Pixel_support",1,M_Pixel_support);

    //------------------Create TOP volume----------------------------
    top = geom->MakeBox("top",Air,500,500,500);
    geom->SetTopVolume(top);
    geom->SetTopVisible(0);
    // If you want to see the boundary, please input the number, 1 instead of 0.
    // Like this, geom->SetTopVisible(1);


    TGeoVolume *inner_field_tube       = geom->MakeTube("inner_field_tube",Iron,49.5,50.0,200.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *outer_field_tube       = geom->MakeTube("outer_field_tube",Me_outer_tube,199.5,220.0,200.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *outer_field_cap        = geom->MakeTube("outer_field_cap",Me_outer_cap,199.5,220.0,2.0);  // r_min, r_max, dz (half of total length)
    // TGeoVolume *outer_field_tube       = geom->MakeTube("outer_field_tube",Me_outer_tube,219.5,220.0,200.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_central_part       = geom->MakeTube("IDS_central_part",Me_IDS,42.8/2.0,43.0/2.0,56.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_side_parts         = geom->MakeTube("IDS_side_parts",Me_IDS,79.3/2.0,79.5/2.0,(222.7-64.0)/2.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_connection_parts_R = geom->MakeCone("IDS_connection_parts_R",Me_IDS,(64.0-56.0)/2.0,42.8/2.0,43.0/2.0,79.3/2.0,79.5/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max
    TGeoVolume *IDS_connection_parts_L = geom->MakeCone("IDS_connection_parts_L",Me_IDS,(64.0-56.0)/2.0,79.3/2.0,79.5/2.0,42.8/2.0,43.0/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max

    TGeoVolume *beampipe_central_part       = geom->MakeTube("beampipe_central_part",Me_beampipe,4.05/2.0,4.15/2.0,141.5);  // r_min, r_max, dz (half of total length)
    TGeoVolume *beampipe_side_parts         = geom->MakeTube("beampipe_side_parts",Me_beampipe,9.52/2.0,9.62/2.0,100.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *beampipe_connection_parts_R = geom->MakeCone("beampipe_connection_parts_R",Me_beampipe,(191.5-141.5)/2.0,4.05/2.0,4.15/2.0,9.52/2.0,9.62/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max
    TGeoVolume *beampipe_connection_parts_L = geom->MakeCone("beampipe_connection_parts_L",Me_beampipe,(191.5-141.5)/2.0,9.52/2.0,9.62/2.0,4.05/2.0,4.15/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max

    TGeoVolume *Pixel_support       = geom->MakeTube("Pixel_support",Me_Pixel_support,21.8/2.0,22.0/2.0,56.0);  // r_min, r_max, dz (half of total length)

    inner_field_tube       ->SetLineColor(4);
    // outer_field_tube       ->SetLineColor(kRed-8);
    outer_field_tube       ->SetLineColor(4);
    outer_field_cap        ->SetLineColor(4);
    IDS_central_part       ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_side_parts         ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_connection_parts_R ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_connection_parts_L ->SetLineColor(2);  // Inner Detector Support (IDS)

    beampipe_central_part       ->SetLineColor(3);  // (beampipe)
    beampipe_side_parts         ->SetLineColor(3);  // (beampipe)
    beampipe_connection_parts_R ->SetLineColor(3);  // (beampipe)
    beampipe_connection_parts_L ->SetLineColor(3);  // (beampipe)

    Pixel_support ->SetLineColor(kYellow-3);  // (pixel support)

    top->AddNodeOverlap(inner_field_tube,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(outer_field_tube,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(outer_field_cap,1,new TGeoTranslation(0,0,201));
    //top->AddNodeOverlap(IDS_central_part,1,new TGeoTranslation(0,0,0));
    //top->AddNodeOverlap(IDS_side_parts,1,new TGeoTranslation(0,0,64.0+(222.7-64.0)/2.0));
    //top->AddNodeOverlap(IDS_side_parts,1,new TGeoTranslation(0,0,-(64.0+(222.7-64.0)/2.0)));
    //top->AddNodeOverlap(IDS_connection_parts_R,1,new TGeoTranslation(0,0,56.0+(64.0-56.0)/2.0));
    //top->AddNodeOverlap(IDS_connection_parts_L,1,new TGeoTranslation(0,0,-(56.0+(64.0-56.0)/2.0)));

    top->AddNodeOverlap(beampipe_central_part,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(beampipe_side_parts,1,new TGeoTranslation(0,0,191.5+100.0));
    top->AddNodeOverlap(beampipe_side_parts,1,new TGeoTranslation(0,0,-(191.5+100.0)));
    top->AddNodeOverlap(beampipe_connection_parts_R,1,new TGeoTranslation(0,0,141.4+(191.5-141.5)/2.0));
    top->AddNodeOverlap(beampipe_connection_parts_L,1,new TGeoTranslation(0,0,-(141.4+(191.5-141.5)/2.0)));

    //top->AddNodeOverlap(Pixel_support,1,new TGeoTranslation(0,0,0));

    top->DrawClone("ogl");

    Color_t TPC_color = kGray+2;

    const Int_t n_TPC_points = 50;
    TPolyLine3D   *TPC_endcaps[5];
    TPolyLine3D   *TPC_tube[4];
    TPolyLine3D   *TPC_tube_lines[n_TPC_points+1];

    Float_t radius_table[5] = {200,200,3.81,3.81,220};
    Float_t z_val_table[5]  = {200,-200,200,-200,200};

    Float_t radius_table_tube[4] = {50,50,50,50};
    Float_t z_val_table_tube[4]  = {200,-200,100,-100};

    for(Int_t r = 0; r < 5; r++)
    {
        if(r == 1) continue;
        TPC_endcaps[r] = new TPolyLine3D();
        Float_t radius   = radius_table[r];
        Float_t x_offset = 0.0;
        Float_t y_offset = 0.0;
        Float_t z_tpc_val   = z_val_table[r];
        for(Int_t t = 0; t < n_TPC_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_TPC_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            TPC_endcaps[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
        }
        TPC_endcaps[r]->SetLineStyle(0);
        // TPC_endcaps[r]->SetLineColor(TPC_color); // 28
        TPC_endcaps[r]->SetLineColor(4); // 28
        TPC_endcaps[r]->SetLineWidth(2);
        TPC_endcaps[r]->DrawClone("ogl");
    }

    for(Int_t r = 0; r < 1; r++)
    {
        TPC_tube[r] = new TPolyLine3D();
        Float_t radius   = radius_table_tube[r];
        Float_t x_offset = 0.0;
        Float_t y_offset = 0.0;
        Float_t z_tpc_val   = z_val_table_tube[r];
        for(Int_t t = 0; t < n_TPC_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_TPC_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            TPC_tube[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
            if(r == 0 && (t%4 == 0))
            {
                TPC_tube_lines[t] = new TPolyLine3D();
                TPC_tube_lines[t]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
                TPC_tube_lines[t]->SetNextPoint(x_tpc_val,y_tpc_val,z_val_table_tube[r+1]);
                TPC_tube_lines[t]->SetLineStyle(0);
                TPC_tube_lines[t]->SetLineColor(TPC_color); // 28
                TPC_tube_lines[t]->SetLineWidth(1);
                // TPC_tube_lines[t]->DrawClone("ogl");
            }
        }
        TPC_tube[r]->SetLineStyle(0);
        TPC_tube[r]->SetLineColor(TPC_color); // 28
        TPC_tube[r]->SetLineWidth(2);
        TPC_tube[r]->DrawClone("ogl");
    }

    TPolyLine3D   *BeamLine;
    BeamLine       = new TPolyLine3D(2);
    BeamLine   ->SetPoint(0,0,0,-750);
    BeamLine   ->SetPoint(1,0,0,750);
    BeamLine   ->SetLineStyle(0);
    BeamLine   ->SetLineColor(kOrange);
    BeamLine   ->SetLineWidth(6);
    BeamLine   ->DrawClone("ogl");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_Circle_Detector_3D(Float_t z, Float_t radius_in, Float_t radius_out,const Int_t n_radii,const Int_t n_delta_phi, Float_t color)
{
    const Int_t n_points = 50;
    TPolyLine3D   *tp_Circles[n_radii];
    TPolyLine3D   *tp_Radial[n_delta_phi];
    Float_t radius_table[n_radii];
    Float_t delta_radius;
    if(n_radii > 1) {delta_radius = (radius_out-radius_in)/((Float_t)(n_radii-1));}
    else{delta_radius = 0.0;}
    Float_t delta_phi    = 2.0*TMath::Pi()/((Float_t)n_delta_phi);
    Float_t z_tpc_val    = z;
    Float_t x_offset     = 0.0;
    Float_t y_offset     = 0.0;

    for(Int_t r = 0; r < n_radii; r++)
    {
        radius_table[r] = radius_in + r*delta_radius;
        tp_Circles[r] = new TPolyLine3D();
        Float_t radius   = radius_table[r];
        for(Int_t t = 0; t < n_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            tp_Circles[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
        }
        tp_Circles[r]->SetLineStyle(0);
        tp_Circles[r]->SetLineColor(color); // 28
        tp_Circles[r]->SetLineWidth(1);
        tp_Circles[r]->DrawClone("ogl");
    }

    for(Int_t r = 0; r < n_delta_phi; r++)
    {
        tp_Radial[r] = new TPolyLine3D();
        Float_t phi_val = r*delta_phi;
        for(Int_t t = 0; t < 2; t++)
        {
            Float_t radius;
            if(t == 0) {radius = radius_table[0];}
            else {radius = radius_table[n_radii-1];}
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            tp_Radial[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
        }
        tp_Radial[r]->SetLineStyle(0);
        tp_Radial[r]->SetLineColor(color); // 28
        tp_Radial[r]->SetLineWidth(1);
        tp_Radial[r]->DrawClone("ogl");
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_Pseudo_Rapidity_Lines(Float_t eta, Float_t z_min, Float_t z_max, const Int_t n_delta_phi, Int_t color)
{
    TPolyLine3D   *tp_eta_lines[n_delta_phi];
    TLorentzVector tl_vector;
    tl_vector.SetPtEtaPhiM(1.0,eta,0.0,0.1);
    Float_t theta = tl_vector.Theta();
    //cout << "theta = " << TMath::RadToDeg()*theta << endl;
    Float_t delta_phi    = 2.0*TMath::Pi()/((Float_t)n_delta_phi);
    Float_t z_val[2] = {z_min,z_max};
    for(Int_t r = 0; r < n_delta_phi; r++)
    {
        tp_eta_lines[r] = new TPolyLine3D();
        Float_t phi_val = r*delta_phi;
        for(Int_t t = 0; t < 2; t++)
        {
            Float_t radius = TMath::Tan(theta)*z_val[t];
            Float_t x_tpc_val = TMath::Cos(phi_val)*radius;
            Float_t y_tpc_val = TMath::Sin(phi_val)*radius;
            tp_eta_lines[r] ->SetNextPoint(x_tpc_val,y_tpc_val,z_val[t]);
            //cout << "r = " << r << ", t = " << t << ", x = " << x_tpc_val << ", y = " << y_tpc_val << ", z = " << z_val[t] << endl;
        }
        tp_eta_lines[r]->SetLineStyle(0);
        tp_eta_lines[r]->SetLineColor(color); // 28
        tp_eta_lines[r]->SetLineWidth(1);
        tp_eta_lines[r]->DrawClone("ogl");
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_Pseudo_Rapidity_Circles(Float_t eta, Float_t z, Int_t color)
{
    const Int_t n_TPC_points = 50;
    TPolyLine3D *tp_circle_line = new TPolyLine3D();
    TLorentzVector tl_vector;
    tl_vector.SetPtEtaPhiM(1.0,eta,0.0,0.1);
    Float_t theta = tl_vector.Theta();
    Float_t radius = TMath::Tan(theta)*z;
    Float_t x_offset     = 0.0;
    Float_t y_offset     = 0.0;

    for(Int_t t = 0; t < n_TPC_points+1; t++)
    {
        Float_t phi_val = ((Float_t)t/(Float_t)n_TPC_points)*(2.0*TMath::Pi());
        Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
        Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
        tp_circle_line->SetNextPoint(x_tpc_val,y_tpc_val,z);
    }

    tp_circle_line->SetLineStyle(0);
    tp_circle_line->SetLineColor(color); // 28
    tp_circle_line->SetLineWidth(1);
    tp_circle_line->DrawClone("ogl");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Float_t Get_Radius_from_eta_z(Float_t eta, Float_t z)
{
    TLorentzVector tl_vector;
    tl_vector.SetPtEtaPhiM(1.0,eta,0.0,0.1);
    Float_t theta = tl_vector.Theta();
    Float_t radius = fabs(TMath::Tan(theta)*z);

    return radius;
}
//----------------------------------------------------------------------------------------







//----------------------------------------------------------------------------------------
void plotFig2_STAREventDisplay(Int_t Event, Int_t eBeamTime)
{
    cout << "Plot_STAR_Event macro started" << endl;

    // BeamTime_num:
    // 0 = 7.7
    // 1 = 11.5
    // 2 = 39
    // 3 = 62.4
    // 4 = 19.6
    // 5 = 27
    // 6 = 200
    // 7 = 14.5
    // 8 = pp200
    // 9 = pp510
    // 10 = run 10 200 GeV Au+Au with BEMC

    Int_t flag_draw_towers  = 0; // 0 = do not draw BEMC towers
    Int_t flag_Draw_STAR_3D = 1;

    //**************************** Set graphic style ***************************************
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    //gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetLabelSize(0.07,"X");
    gStyle->SetLabelSize(0.07,"Y");
    gStyle->SetTitleSize(0.07,"X");
    gStyle->SetTitleSize(0.07,"Y");
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds,greens, blues, NCont);
    gStyle->SetNumberContours(NCont);
    //**************************************************************************************

    gStyle->SetEndErrorSize(3);
    TRandom3 r3b;
    r3b.SetSeed(0); // seed for random number generator changes every second
    gRandom->SetSeed(0);


    //----------------------------------------------------------------------------------------
    TString HistName;
    char NoP[50];
    //----------------------------------------------------------------------------------------


    //----------------------------------------------------------------------------------------
    cout << "Define input" << endl;
    //HistName = "./Data/Event_tracks_";
    //HistName += Event;
    //HistName += ".root";

    // if(eBeamTime == 0) HistName = "./Data/Event_tracks_event_12_7GeV.root";
    // if(eBeamTime == 1) HistName = "./Data/Event_tracks_event_19_11GeV.root";
    // if(eBeamTime == 2) HistName = "./Data/Event_tracks_event_20_39GeV.root";
    // if(eBeamTime == 3) HistName = "./Data/Event_tracks_event_6_62GeV.root";
    // if(eBeamTime == 4) HistName = "./Data/Event_tracks_event_3_19GeV.root";
    // if(eBeamTime == 5) HistName = "./Data/Event_tracks_event_6_27GeV.root";
    // //if(eBeamTime == 6) HistName = "./Data/Event_tracks_event_21_200GeV.root";
    // if(eBeamTime == 6) HistName = "./Data/Event_tracks_event_7_62GeV.root";
    // //if(eBeamTime == 7) HistName = "./Data/Event_tracks_event_1_15GeV.root";
    // if(eBeamTime == 7) HistName = "./Data/Event_tracks_event_5_15GeV.root";
    // if(eBeamTime == 8) HistName = "./Data/Event_tracks_event_88_pp200GeV.root";
    // //if(eBeamTime == 9) HistName = "./Data/Event_tracks_event_2_510GeV.root";
    // if(eBeamTime == 9) HistName = "./Data/Event_tracks_event_22_pp510GeV.root";
    // if(eBeamTime == 10) HistName = "./Data/Event_tracks_event_602_200GeV.root";

    if(eBeamTime == 0) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_12_7GeV.root";
    if(eBeamTime == 1) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_19_11GeV.root";
    if(eBeamTime == 2) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_20_39GeV.root";
    if(eBeamTime == 3) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_6_62GeV.root";
    if(eBeamTime == 4) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_3_19GeV.root";
    if(eBeamTime == 5) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_6_27GeV.root";
    //if(eBeamTime == 6) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_21_200GeV.root";
    if(eBeamTime == 6) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_7_62GeV.root";
    //if(eBeamTime == 7) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_1_15GeV.root";
    if(eBeamTime == 7) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_5_15GeV.root";
    if(eBeamTime == 8) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_88_pp200GeV.root";
    //if(eBeamTime == 9) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_2_510GeV.root";
    if(eBeamTime == 9) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_22_pp510GeV.root";
    if(eBeamTime == 10) HistName = "/Users/xusun/WorkSpace/STAR/Data/EventDisplay/Event_tracks_event_602_200GeV.root";

    cout << "Input file: " << HistName.Data() << endl;

    TFile* inputfile_tracks = TFile::Open(HistName.Data());  // open the file

    TH1F* h_params           = (TH1F*)inputfile_tracks->Get("h_params");
    TH1I* h_track_pid        = (TH1I*)inputfile_tracks->Get("h_track_pid");

    const Int_t n_tracks = h_params->GetBinContent(1);
    cout << "number of track in file: " << n_tracks << endl;

    TPolyLine3D   *pTrack[n_tracks];
    TPolyLine3D   *pTrack_extrapolate[n_tracks];
    TPolyLine3D   *pTrack_extrapolate_to_Tof[n_tracks];
    TMarker3DBox  *pTowerHit[n_tracks];
    TMarker3DBox  *pTofHit[n_tracks];

    for(Int_t i_track = 0; i_track < n_tracks; i_track++)
    {
        //printf("i_track: %d, out of %d \n",i_track,n_tracks);
        HistName  = "pTrack_";
        HistName  += i_track;
        pTrack[i_track] = (TPolyLine3D*)inputfile_tracks->Get(HistName.Data());
        HistName  = "pTrack_extrapolate_";
        HistName  += i_track;
        pTrack_extrapolate[i_track] = (TPolyLine3D*)inputfile_tracks->Get(HistName.Data());
        HistName  = "pTrack_extrapolate_to_Tof_";
        HistName  += i_track;
        pTrack_extrapolate_to_Tof[i_track] = (TPolyLine3D*)inputfile_tracks->Get(HistName.Data());
        HistName  = "pTowerHit_";
        HistName  += i_track;
        pTowerHit[i_track] = (TMarker3DBox*)inputfile_tracks->Get(HistName.Data());
        HistName  = "pTofHit_";
        HistName  += i_track;
        pTofHit[i_track] = (TMarker3DBox*)inputfile_tracks->Get(HistName.Data());
    }


    Int_t color_array_particle_pidA[16] = {16,16,kAzure-2,kMagenta+2,16,16,16,16,kCyan+2,kRed-6,16,4,9,16,kRed+1,kRed+4};
    Int_t color_array_particle_pidB[16] = {16,16,kGreen,kMagenta+0,16,16,16,16,kBlue-3,kRed+2,16,4,9,16,kOrange,kYellow+2};

    printf("Draw c_3D \n");
    TCanvas* c_3D    = new TCanvas("c_3D","c_3D",10,10,800,800);

    if(flag_Draw_STAR_3D)
    {
        Draw_STAR_3D();
        // Draw_Track();
        top->DrawClone("ogl");
    }



#if 0
    // TEve test with arrow
    //------------------------
    TEveManager::Create();

    TEveGeoShape* s = 0;

    s = new TEveGeoShape("Tube");
    s->SetShape(new TGeoTube(0, 10, 200));
    s->SetMainColor(kRed);
    // s->SetMainTransparency(50);
    // s->SetNSegments(12);
    gEve->AddElement(s);

    s = new TEveGeoShape("Cone");
    s->SetShape(new TGeoCone(20, 0, 40, 0, 0));
    s->SetMainColor(kRed);
    // s->SetMainTransparency(50);
    // s->SetNSegments(12);
    s->RefMainTrans().MoveLF(3, 200 + 10); // See docs of TEveTrans.
    gEve->AddElement(s);

    gEve->Redraw3D(kTRUE);
    //------------------------
#endif


    /*
    for(Int_t i_track = 0; i_track < n_tracks; i_track++)
    //for(Int_t i_track = 0; i_track < 20; i_track++)
    {
        Int_t i_pid = h_track_pid->GetBinContent(i_track);
        //if(i_pid == 2 || i_pid == 3)
        //if(i_pid == 8)
        {
            cout << "i_track = " << i_track << endl;
            pTrack[i_track]             ->SetLineWidth(2);
            pTrack[i_track]             ->SetLineColor(color_array_particle_pidA[i_pid]);
            pTrack[i_track]             ->DrawClone("ogl");

            pTrack[i_track]             ->SetLineWidth(3);
            pTrack[i_track]             ->SetLineColor(color_array_particle_pidB[i_pid]);
            pTrack[i_track]             ->DrawClone("ogl");

            pTrack[i_track]             ->SetLineWidth(2);
            pTrack[i_track]             ->SetLineColor(color_array_particle_pidA[i_pid]);
            pTrack[i_track]             ->DrawClone("ogl");


            //pTrack_extrapolate[i_track] ->DrawClone("ogl");
        }
    }
    */


    // PID:
    // 2  = e+
    // 3  = e-
    // 8  = pi+
    // 9  = pi-
    // 11 = K+
    // 12 = K-
    // 14 = p
    // 15 = anti-p

    Int_t color_array_particle_pidC[16] = {kAzure-1,kAzure-1,kRed+1,kRed+1,kAzure-1,kAzure-1,kAzure-1,kAzure-1,kAzure-1,kAzure-1,kAzure-1,kAzure-1,kAzure-1,kAzure-1,kCyan+3,kCyan+3};
    Int_t color_array_particle_pidC_inner[16] = {kAzure-3,kAzure-3,kRed+1,kRed+1,kAzure-3,kAzure-3,kAzure-3,kAzure-3,kAzure-3,kAzure-3,kAzure-3,kAzure-3,kAzure-3,kAzure-3,kCyan+3,kCyan+3};

    Int_t track_sel_A = 66;
    Int_t track_sel_B = 83;

    Int_t track_Kstar_A = 12;
    Int_t track_Kstar_B = 246;

    printf("Loop over tracks \n");
    for(Int_t i_track = 0; i_track < n_tracks; i_track++)
    // for(Int_t i_track = 0; i_track < 30; i_track++)
    {
        Int_t i_pid = h_track_pid->GetBinContent(i_track);
        {
            pTrack[i_track]             ->SetLineWidth(1);
            //pTrack[i_track]             ->SetLineColor(color_array_particle_pidC[i_pid]);
            // pTrack[i_track]             ->SetLineColorAlpha(kGray,0.2);
            pTrack[i_track]             ->SetLineColorAlpha(1,0.2);

	    // if(i_pid == 11 || i_pid == 12)
	    // {
	    //   cout << "i_track = " << i_track << ", i_pid: " << i_pid <<endl;
	    // }

	    // phi candidate
	    if(i_pid == 11 && i_track == track_sel_A)
	    {
	      pTrack[i_track] ->SetLineWidth(20);
	      pTrack[i_track] ->SetLineColorAlpha(kMagenta,1.0);
	    }
	    if(i_pid == 12 && i_track == track_sel_B)
	    {
	      pTrack[i_track] ->SetLineWidth(20);
	      pTrack[i_track] ->SetLineColorAlpha(kCyan,1.0);
	    }

	    // K* candidate
	    if(i_pid == 8 && i_track == track_Kstar_A)
	    {
	      pTrack[i_track] ->SetLineWidth(20);
	      pTrack[i_track] ->SetLineColorAlpha(kRed,1.0);
	    }
	    if(i_pid == 12 && i_track == track_Kstar_B)
	    {
	      pTrack[i_track] ->SetLineWidth(20);
	      pTrack[i_track] ->SetLineColorAlpha(kCyan,1.0);
	    }

            pTrack[i_track]->DrawClone("ogl");

	    // phi extrapolation
	    if(i_pid == 11 && i_track == track_sel_A)
	    {
	      pTrack_extrapolate[i_track] ->SetLineWidth(20);
	      pTrack_extrapolate[i_track] ->SetLineColorAlpha(kMagenta,1.0);
	      pTrack_extrapolate[i_track] ->DrawClone("ogl");
	    }
	    if(i_pid == 12 && i_track == track_sel_B)
	    {
	      pTrack_extrapolate[i_track] ->SetLineWidth(20);
	      pTrack_extrapolate[i_track] ->SetLineColorAlpha(kCyan,1.0);
	      pTrack_extrapolate[i_track] ->DrawClone("ogl");
	    }

	    // KStar extrapolation
	    if(i_pid == 8 && i_track == track_Kstar_A)
	    {
	      pTrack_extrapolate[i_track] ->SetLineWidth(20);
	      pTrack_extrapolate[i_track] ->SetLineColorAlpha(kRed,1.0);
	      pTrack_extrapolate[i_track] ->DrawClone("ogl");
	    }
	    if(i_pid == 12 && i_track == track_Kstar_B)
	    {
	      pTrack_extrapolate[i_track] ->SetLineWidth(20);
	      pTrack_extrapolate[i_track] ->SetLineColorAlpha(kCyan,1.0);
	      pTrack_extrapolate[i_track] ->DrawClone("ogl");
	    }

            Float_t x,y,z;
            pTowerHit[i_track]->GetPosition(x,y,z);
            Float_t dx,dy,dz;
            pTowerHit[i_track]->GetSize(dx,dy,dz);
            TVector3 tower_vec;
            //cout << "i_track = " << i_track << ", dx = " << dx << ", dy = " << dy << ", dz = " << dz << endl;
            tower_vec.SetXYZ(x,y,z);
            tower_vec.SetPerp(205.0 + 0.5+dx);
            pTowerHit[i_track]->SetSize(dx,3.0,3.0);
            pTowerHit[i_track]->SetPosition(tower_vec.X(),tower_vec.Y(),tower_vec.Z());
            //pTowerHit[i_track]->SetFillColorAlpha(kAzure-6,0.5);
            //pTowerHit[i_track]->SetLineColorAlpha(kAzure-6,0.5);
            pTowerHit[i_track]->SetFillColorAlpha(kOrange-2,0.5);
            pTowerHit[i_track]->SetLineColorAlpha(kOrange-2,0.5);
            if(flag_draw_towers) pTowerHit[i_track]->DrawClone("ogl");

            Float_t xTOF,yTOF,zTOF;
            pTofHit[i_track]->GetPosition(xTOF,yTOF,zTOF);
            Float_t dxTOF,dyTOF,dzTOF;
            pTofHit[i_track]->GetSize(dxTOF,dyTOF,dzTOF);
            TVector3 tof_vec;
            //cout << "i_track = " << i_track << ", dx = " << dx << ", dy = " << dy << ", dz = " << dz << endl;
            tof_vec.SetXYZ(xTOF,yTOF,zTOF);
            tof_vec.SetPerp(201.0 + 0.5+dxTOF);
            //pTofHit[i_track]->SetSize(dxTOF,4.0,4.0);
            pTofHit[i_track]->SetSize(2.0,2.0,2.0);
            pTofHit[i_track]->SetPosition(tof_vec.X(),tof_vec.Y(),tof_vec.Z());
            pTofHit[i_track]->SetFillColor(kAzure+3);
            pTofHit[i_track]->SetLineColor(kAzure+3);


            if(i_pid == 11 && i_track == track_sel_A)
            {
                pTofHit[i_track]->SetFillColor(kMagenta);
                pTofHit[i_track]->SetLineColor(kMagenta);
                pTofHit[i_track]->SetSize(4.0,4.0,4.0);
            }
            if(i_pid == 12 && i_track == track_sel_B)
            {
                pTofHit[i_track]->SetFillColor(kCyan);
                pTofHit[i_track]->SetLineColor(kCyan);
                pTofHit[i_track]->SetSize(4.0,4.0,4.0);
            }

            if(i_pid == 8 && i_track == track_Kstar_A)
            {
                pTofHit[i_track]->SetFillColor(kRed);
                pTofHit[i_track]->SetLineColor(kRed);
                pTofHit[i_track]->SetSize(4.0,4.0,4.0);
            }
            if(i_pid == 12 && i_track == track_Kstar_B)
            {
                pTofHit[i_track]->SetFillColor(kCyan);
                pTofHit[i_track]->SetLineColor(kCyan);
                pTofHit[i_track]->SetSize(4.0,4.0,4.0);
            }

            pTofHit[i_track]->DrawClone("ogl");
        }
    }
}



