#include "Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "Utility/type.h"
#include <string>

void fitcos_1st(const int energy = 4, const int pid = 0, bool doall = false, bool isBesI = false, bool random3D = false) {

  gROOT->Reset();
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(22);
  gStyle->SetLabelFont(22);
  gStyle->SetTitleFont(22);

  Double_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};
  Double_t pt_set[8] = {0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};
//  Double_t pt_set[8] = {3.0, 3.3, 3.6, 3.9, 4.3, 4.6, 4.9, 5.2};
//  Double_t pt_set[8] = {0.6, 1.4, 2.2, 3.0, 3.8, 4.6, 5.4, 7.2};

  TH1D *h_theta_star_before[6];
  TH1D *h_theta[6];
  TH1D *h_theta_star[6];
  TH1D *h_out1[6];
  TH1D *h_out2[6];

  TFile *MCFiles[6];
  double Fval[6] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Ferr[6] = {0.0};
  double Fpval[6]={0.0};
  double Fperr[6]={0.0};
  double FfromFp[6]={0.0};
  double FfromFperr[6]={0.0};
  double FvalBESI[6] = {0.0,-0.00989006, -0.0102287, -0.0102549, -0.00800524, -0.00767692};

  for(int ipt = 2; ipt < 6; ipt++)
  {
    MCFiles[ipt] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/50MEvents_OnlyEtaCut/McAcceptanceOutput_pt%d_energy%d_pid%d_cent9.root",vmsa::mPID[pid].c_str(),ipt+1,energy,pid),"READ");

    //if(ipt == 2) MCFiles[ipt] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/McAcceptanceOutput_pt%d_energy%d_pid%d_cent9_00004999.root",vmsa::mPID[pid].c_str(),ipt+1,energy,pid),"READ");

    h_theta_star_before[ipt] = (TH1D*) MCFiles[ipt]->Get("h_theta_star_before");
    h_theta[ipt] = (TH1D*) MCFiles[ipt]->Get("h_theta");
    h_theta_star[ipt] = (TH1D*) MCFiles[ipt]->Get("h_theta_star");
    h_out1[ipt] = (TH1D*) MCFiles[ipt]->Get("h_out1");
    h_out2[ipt] = (TH1D*) MCFiles[ipt]->Get("h_out2");
 
    TF1 *Func_rho = new TF1("Func_rho","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
    TF1 *Func_A = new TF1("Func_A","[0]*(1.+[1]*(x*x))",0,1);

    Func_rho->SetParameter(0,h_theta_star[ipt]->GetBinContent(1));
    Func_rho->SetParameter(1,1./3.);
    h_theta_star[ipt]->Fit(Func_rho,"ERQ");
    cout << "rho00 = " << Func_rho->GetParameter(1) << endl;;    


    Func_A->SetParameter(0,h_theta[ipt]->GetBinContent(1));
    Func_A->SetParameter(1,0);
    h_theta[ipt]->Fit(Func_A,"ER");
    Fpval[ipt] = Func_A->GetParameter(1);
    Fperr[ipt] = Func_A->GetParError(1);


    TH1D *h_theta_star_clone = (TH1D*)h_theta_star[ipt]->Clone("h_theta_star_clone");
    h_theta_star_clone->Sumw2();
    h_theta_star_before[ipt]->Sumw2();
    h_theta_star_clone->Divide(h_theta_star_before[ipt]); 
    Func_A->SetParameter(0,h_theta_star_clone->GetBinContent(1));
    Func_A->SetParameter(1,0);
    h_theta_star_clone->Fit(Func_A,"ER");

    Fval[ipt] = Func_A->GetParameter(1);
    Ferr[ipt] = Func_A->GetParError(1);

    FfromFp[ipt] = -Fpval[ipt]/(2.+Fpval[ipt]);
    FfromFperr[ipt] = TMath::Abs((Fpval[ipt]/(2.-Fpval[ipt])/(2.-Fpval[ipt])-1./(2.-Fpval[ipt]))*Fperr[ipt]);

    TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
    c_play->SetLeftMargin(0.15);
    c_play->SetBottomMargin(0.15);
    c_play->SetGrid(0,0);
    c_play->SetTicks(1,1);
    c_play->cd();
    h_theta_star_clone->Draw("pE");
    Func_A->Draw("same");
    c_play->SaveAs(Form("FValues/pt%d.pdf",ipt));

    delete c_play;
    delete h_theta_star_clone;
    delete Func_rho;
    delete Func_A;
  }
 
  if(isBesI)
  { 
    Fval[1] = -0.00989006;
    Fval[2] = -0.0102287;
    Fval[3] = -0.0102549;
    Fval[4] = -0.00800524;
    Fval[5] = -0.000652552;
  }
  //double FfromF[6] ={0.0};
  for(int ipt = 2; ipt < 6; ipt++)
  {
    //FfromF[ipt] = -Fpval[ipt]/(2.+Fpval[ipt]);
    cout << "pt bin: " <<  pt_set[ipt] << "-" << pt_set[ipt+1] << std::fixed << std::setprecision(7) << "GeV/c    F = " << Fval[ipt] << " +/- " << Ferr[ipt] << "    F* = " << Fpval[ipt] << " +/- " << Fperr[ipt] << "    F from F* = " << FfromFp[ipt] << " +/- " << FfromFperr[ipt] << endl;   
    //cout << "pt bin: " << std::fixed << std::setprecision(5) << pt_set[ipt] << "-" << pt_set[ipt+1] << "GeV/c    F = " << Fval[ipt] << " +/- " << Ferr[ipt] << "    F' = " << Fpval[ipt] << " +/- " << Fperr[ipt] << "    F from F' (swapping them)= " << -Fval[ipt]/(2.+Fval[ipt]) << " +/- " << TMath::Abs((Fval[ipt]/(2.-Fval[ipt])/(2.-Fval[ipt])-1./(2.-Fval[ipt]))*Ferr[ipt]) << endl;   
  //out<<"-D'/(2+D'): "<<-D_theta/(2.+D_theta)<<endl;
  }

  Char_t Char[9][10] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
  Double_t centCent[9] = {75.0,65.0,55.0,45.0,35.0,25.0,15.0,7.5,2.5};

  //double eta_D[6] = {0.0855382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};//19GeV


  TF1 *line = new TF1("line","1/3",-0.5,8.5);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
     
  TFile *fres1 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/Resolution_file_19GeV_EpdCorrections_4.root","READ");
  TFile *fres2 = new TFile("../TreeProduction/OldUtilityFilesEta1OfficialCent/file_19GeV_Resolution.root","READ");
  TFile  *fd12 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/file_19GeV_EpdFlow.root","READ");
  //TFile event("all_event_output.root");
  
  TH1D *DeltaPsi1  = (TH1D*) fres1->Get("AveCosDeltaPsi1");
  TH1D *DeltaPsi2  = (TH1D*) fres2->Get("p_mRes2_Sub");
  TH1D *DeltaPsi12 = (TH1D*)  fd12->Get("p_mD12");

  TH1D *resolution_1 = new TH1D("resolution_1","resolution_1", 9, -0.5, 8.5);
  TH1D *resolution_sub1 = new TH1D("resolution_sub1","resolution_sub1", 9, -0.5, 8.5);
  TH1D *resolution_sub2 = new TH1D("resolution_sub2","resolution_sub2", 9, -0.5, 8.5);
  TH1D *resolution_12 = new TH1D("resolution_12","resolution_12", 9, -0.5, 8.5);
  TH1D *resolution_2 = new TH1D("resolution_2","resolution_2", 9, -0.5, 8.5);

  TGraphAsymmErrors *g_res1  = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_ressub1  = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_ressub2  = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_res2  = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_res12 = new TGraphAsymmErrors();

  for(Int_t cent=0; cent<9; cent++) {

    Double_t CosMean = DeltaPsi1->GetBinContent(cent+1); // load from resolution file
    Double_t CosError = DeltaPsi1->GetBinError(cent+1);
    Double_t ZDCSMD_resSub = CosMean>0.? (TMath::Sqrt(CosMean)) : 0.;
    Double_t ZDCSMD_resSubErr = CosMean>0.? CosError/2./CosMean : 0.;
    Double_t ZDCSMD_chiSub = chi(ZDCSMD_resSub);
    Double_t ZDCSMD_chiSubDelta = chi(ZDCSMD_resSub+0.005);
    Double_t ResZ = resEventPlane(TMath::Sqrt(2.)*ZDCSMD_chiSub);
    Double_t ResZDelta = resEventPlane(TMath::Sqrt(2.)*ZDCSMD_chiSubDelta);
    Double_t ResErrZ = ZDCSMD_resSubErr * TMath::Abs((ResZ-ResZDelta)/0.005);

    resolution_sub1->SetBinContent(cent+1, ZDCSMD_resSub);
    resolution_sub1->SetBinError(cent+1, ZDCSMD_resSubErr);
    resolution_1->SetBinContent(cent+1, ResZ);
    resolution_1->SetBinError(cent+1, ResErrZ);

    Double_t Cos12 = DeltaPsi12->GetBinContent(cent+1)/TMath::Sqrt(2.); // load file for D12 
    Double_t CosError12 = DeltaPsi12->GetBinError(cent+1)/TMath::Sqrt(2.);
    Double_t mean = Cos12/ResZ;
    Double_t error = ResErrZ*ResErrZ/ResZ/ResZ + CosError12*CosError12/Cos12/Cos12;

    resolution_12->SetBinContent(cent+1,mean);
    resolution_12->SetBinError(cent+1, mean*TMath::Sqrt(error));

    Double_t CosMean_2 = DeltaPsi2->GetBinContent(cent+1);
    Double_t CosError_2 = DeltaPsi2->GetBinError(cent+1);
    Double_t ZDCSMD_resSub_2 = CosMean_2>0.? (TMath::Sqrt(CosMean_2)) : 0.;
    Double_t ZDCSMD_resSubErr_2 = CosMean_2>0.? CosError_2/2./CosMean_2 : 0.;
    Double_t ZDCSMD_chiSub_2 = chi(ZDCSMD_resSub_2);
    Double_t ZDCSMD_chiSubDelta_2 = chi(ZDCSMD_resSub_2+0.005);
    Double_t ResZ_2 = res1(TMath::Sqrt(2.)*ZDCSMD_chiSub_2);
    Double_t ResZDelta_2 = res1(TMath::Sqrt(2.)*ZDCSMD_chiSubDelta_2);
    Double_t ResErrZ_2 = ZDCSMD_resSubErr_2 * TMath::Abs((ResZ_2-ResZDelta_2)/0.005);

    resolution_2->SetBinContent(cent+1,ResZ_2);
    resolution_2->SetBinError(cent+1,ResErrZ_2);
    cout<<"First order Cent"<<cent_set[cent+1]<<"_"<<cent_set[cent]<<": "<<ZDCSMD_resSub<<" +/- " << ZDCSMD_resSubErr<<endl;
    cout<<"Cent"<<cent_set[cent+1]<<"_"<<cent_set[cent]<<": "<<ResZ<<" +/- "<<ResErrZ<<" "<<ResZ_2<<" +/- "<<ResErrZ_2<<" "<<mean<<" +/- "<<mean*TMath::Sqrt(error)<<endl;

    g_res1->SetPoint(cent,centCent[cent],ResZ);
    g_res1->SetPointError(cent,0.0,0.0,ResErrZ,ResErrZ);
    g_ressub1->SetPoint(cent,centCent[cent],ZDCSMD_resSub);
    g_ressub1->SetPointError(cent,0.0,0.0,ZDCSMD_resSubErr,ZDCSMD_resSubErr);
    g_ressub2->SetPoint(cent,centCent[cent],ZDCSMD_resSub_2);
    g_ressub2->SetPointError(cent,0.0,0.0,ZDCSMD_resSubErr_2,ZDCSMD_resSubErr_2);
    g_res2->SetPoint(cent,centCent[cent],ResZ_2);
    g_res2->SetPointError(cent,0.0,0.0,ResErrZ_2,ResErrZ_2);
    g_res12->SetPoint(cent,centCent[cent],mean);
    g_res12->SetPointError(cent,0.0,0.0,mean*TMath::Sqrt(error),mean*TMath::Sqrt(error));
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.15);
  c_play->SetGrid(0,0);
  c_play->SetTicks(1,1);
  c_play->cd();

  TH1F *h_play = new TH1F("h_play","h_play",100,0,100);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("Centrality (%)");
  h_play->GetYaxis()->SetTitle("Resolution");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(0,80);
  h_play->GetYaxis()->SetRangeUser(0.0,0.8);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");
  h_play->Draw("pE");

  //g_res1->SetMarkerStyle(24);
  //g_res1->SetMarkerColor(kGreen+2);
  //g_res1->SetMarkerSize(1.5);
  //g_res1->Draw("pE Same");

  //g_ressub1->SetMarkerStyle(20);
  //g_ressub1->SetMarkerColor(kBlue);
  //g_ressub1->SetMarkerSize(2.0);
  //g_ressub1->Draw("pE Same");

  g_ressub2->SetMarkerStyle(20);
  g_ressub2->SetMarkerColor(kBlue);
  //g_ressub2->SetMarkerColor(kOrange+7);
  g_ressub2->SetMarkerSize(2.0);
  g_ressub2->Draw("pE Same");

  //g_res2->SetMarkerStyle(20);
  //g_res2->SetMarkerColor(kGreen+2);
  //g_res2->SetMarkerSize(1.5);
  //g_res2->Draw("pE Same");

  //g_res12->SetMarkerStyle(20);
  //g_res12->SetMarkerColor(kRed);
  //g_res12->SetMarkerSize(1.5);
  //g_res12->Draw("pE Same");

  //p_mEpdSubRes->Draw("pE Same");
  //g_mTpcSubRes2->SetMarkerStyle(24);
  //g_mTpcSubRes2->SetMarkerColor(kAzure+2);
  //g_mTpcSubRes2->SetMarkerSize(1.5);
  //g_mTpcSubRes2->Draw("pE Same");

  //g_mTpcFullRes2->SetMarkerStyle(20);
  //g_mTpcFullRes2->SetMarkerColor(kAzure+2);
  //g_mTpcFullRes2->SetMarkerSize(1.5);
  //g_mTpcFullRes2->Draw("pE Same");

  //g_mTpcSubRes3->SetMarkerStyle(24);
  //g_mTpcSubRes3->SetMarkerColor(kGray+2);
  //g_mTpcSubRes3->SetMarkerSize(1.5);
  //g_mTpcSubRes3->Draw("pE Same");

  //g_mTpcFullRes3->SetMarkerStyle(20);
  //g_mTpcFullRes3->SetMarkerColor(kGray+2);
  //g_mTpcFullRes3->SetMarkerSize(1.5);
  //g_mTpcFullRes3->Draw("pE Same");

  TLegend *leg1 = new TLegend(0.60,0.70,0.85,0.85);
  leg1->SetFillColor(10);
  leg1->SetBorderSize(0);
  //leg1->AddEntry(g_res1,"R_{1} EPD","p");
  //leg1->AddEntry(g_ressub1,"R_{1}^{Sub} EPD","p");
  leg1->AddEntry(g_ressub2,"R_{2}^{Sub} TPC","p");
  //leg1->AddEntry(g_res2,"R_{2} TPC","p");
  //leg1->AddEntry(g_res12,"R_{21}^{Sub} = D_{12}/(R_{1}#sqrt{2})","p");
  //leg->AddEntry(g_mTpcFullRes2,"2^{nd} Full EPD EP","p");
  //leg->AddEntry(g_mTpcSubRes2,"2^{nd} #eta_{sub} EPD EP","p");
  //leg->AddEntry(g_mTpcFullRes3,"3^{rd} RanFull EPD","p");
  //leg->AddEntry(g_mTpcSubRes3,"3^{rd} #eta_{sub}#sqrt{2} EPD","p");//leg->AddEntry(g_mZdcFullRes1,"1^{st} ZDC-SMD Full EP","p");
  //leg->AddEntry(g_mZdcFullRes2,"2^{nd} ZDC-SMD Full EP","p");
  leg1->Draw("same");

  string FigureName = Form("./figures/Resolution/Resolutions_%s.pdf",vmsa::mBeamEnergy[energy].c_str());
  c_play->SaveAs(FigureName.c_str());

  
  double Res_12 = 0;
  double Res_12_weight = 0;
  if(!isBesI)
  {
  for(int cent=2; cent<=5; cent++) {
    Res_12 += resolution_sub1->GetBinContent(cent+1)/resolution_sub1->GetBinError(cent+1)/resolution_sub1->GetBinError(cent+1);
    Res_12_weight += 1./resolution_sub1->GetBinError(cent+1)/resolution_sub1->GetBinError(cent+1);
  }
  Res_12 = Res_12/Res_12_weight;
  cout << Res_12;
  }

  double R21_BESI[4] = {0.204623, 0.296309, 0.360465, 0.396092};
  double R21_BESIe[4] = {0.00703685, 0.00541934, 0.00497943, 0.00549284}; 

  if(isBesI)
  {
  Res_12 = 0.0;
  Res_12_weight = 0.0;
  for(int cent=0; cent<=3; cent++) {
    Res_12 += R21_BESI[cent]/R21_BESIe[cent]/R21_BESIe[cent];
    Res_12_weight += 1./R21_BESIe[cent]/R21_BESIe[cent];
  }
  Res_12 = Res_12/Res_12_weight;
  cout << Res_12;
  }


  if(doall)
  {
  double eff[10][7][7];
  double eff_error[10][7][7];
  //if(eff_mode==0)
  TFile *eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Efficiency/Cos/Eff_%s_SingleParticle_2060.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()),"READ");
  if(isBesI) eff_file = new TFile(Form("../../../FileTransfers/Eff_%s_SingleKaon_second.root",vmsa::mBeamEnergy[energy].c_str()),"READ");
  eff_file->Print();
  //else
  //  TFile *eff_file = new TFile("Eff_19GeV_SingleKaon_eta.root","READ");

  for(int i=0; i<10;i++)
  {
    for(int j=1; j<=6; j++)
    { 
     for(int k=1; k<=7; k++)
      {
        //if(j==7) {
        //  eff[i][j-1][k-1] = 1;
        //  eff_error[i][j-1][k-1] = 0;
        //  continue;
        //}

        Title = new TString("h_mEffCos_Cent_");
        *Title += 9;
        *Title += "_Pt_";
        *Title += (j-1);
        TH1D *eff_hist = (TH1D*)eff_file->Get(Title->Data());
        eff[i][j-1][k-1] = eff_hist->GetBinContent(k);
        eff_error[i][j-1][k-1] = eff_hist->GetBinError(k)/eff[i][j-1][k-1];
        delete Title;
        delete eff_hist;

      }
    }
  }
  eff_file->Close();

  TFile *input = new TFile(Form("rho00/%s/%s/Raw%sPtSys_Psi1.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()),"READ");
  if(random3D) *input = new TFile(Form("rho00/%s/%s/3DRandom/Raw%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()),"READ");
  if(isBesI) input = new TFile(Form("rho00/%s/%s/BESI/Raw%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str()),"READ");

  TF1 *Func_rho = new TF1("Func_rho",FuncAD,0,1,4);
  TF1 *Func_rdl = new TF1("Func_rdl","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);

  TH1D *Big_yield = new TH1D("Big_yield","Big_yield", 7, 0, 1);
  double big_cos[9][2] = {0};

  TString *Title;

  Double_t weight_rho00 = 0;
  Double_t weight_error_rho00 = 0;
  Double_t weight_all = 0;

  Double_t pt_rho00[6] = {0};
  Double_t pt_error_rho00[6] = {0};
  Double_t pt_all[6] = {0};

  TH1D *rho00_hist[6][6][6][6][6];
  TGraphAsymmErrors *g_rho00[6][6][6][6][6];
  //TH1DMap PtCos;
  TH1F *clonePt; 
  
  string outputname = Form("output/%s/%s/AccRes%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  if(isBesI) outputname = Form("output/%s/%s/BESI/AccRes%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *output = new TFile(outputname.c_str(),"RECREATE");
  output->cd();
  Double_t ptbin_rho00[7] = {0.4,0.8,1.2, 1.8, 2.4, 3.0, 4.2};

  for(int i_cent = 9/*vmsa::Cent_start*/; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
        for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
        {
          for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
          {
            for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
            {
            for(int i_F = 0; i_F < 2; i_F++)
            {
  
              Double_t weight_rho00 = 0;
              Double_t weight_error_rho00 = 0;
              Double_t weight_all = 0;

              Double_t pt_rho00[6] = {0};
              Double_t pt_error_rho00[6] = {0};
              Double_t pt_all[6] = {0};

              for(Int_t PtBin=1; PtBin<=6; PtBin++) 
              {
                string key = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",PtBin-1,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());

                TH1F *PtCos_raw = (TH1F*)input->Get(key.c_str())->Clone("PtCos_raw");
                TH1F *PtCos = new TH1F(key.c_str(),key.c_str(), 7, 0, 1);
                //delete Title;
                for(int i=1; i<=7; i++) {
                  float inte_mean = PtCos_raw->GetBinContent(i);
                  float inte_mean_error = PtCos_raw->GetBinError(i);
                  PtCos->SetBinContent(i, inte_mean/eff[i_cent][PtBin-1][i-1]);
                  PtCos->SetBinError(i, inte_mean/eff[i_cent][PtBin-1][i-1]*TMath::Sqrt(inte_mean_error*inte_mean_error/inte_mean/inte_mean+eff_error[i][PtBin-1][i-1]*eff_error[i][PtBin-1][i-1]));
                }

                PtCos->GetXaxis()->SetTitleOffset(1.2);
                PtCos->GetXaxis()->SetTitle("cos#theta*");
                PtCos->GetYaxis()->SetTitle("yield");
                PtCos->GetYaxis()->SetTitleOffset(1.0);
                PtCos->SetMarkerColor(2);
                PtCos->SetMarkerSize(1.8);
                PtCos->SetMarkerStyle(21);
                PtCos->Draw();
                Func_rho->SetParameter(0, PtCos->GetBinContent(5));
                Func_rho->SetParameter(1, 0.3);
                if(i_F == 0) Func_rho->FixParameter(2, Fval[PtBin-1]);
                if(i_F == 1) Func_rho->FixParameter(2, FvalBESI[PtBin-1]);
                Func_rho->FixParameter(3, Res_12);
                if(random3D) Func_rho->FixParameter(3, 0.);
                cout << "Resolution = " << Func_rho->GetParameter(3);
                //cout<<"Fit with real EP:"<<endl;
                PtCos->Fit(Func_rho, "NMI");
                Func_rho->Draw("same");

                Title = new TString(Form("fit/yield_pt_%d_cent_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s.pdf",PtBin-1,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str()));
           
                TCanvas *c1 = new TCanvas();
                c1->SetFillColor(0);
                c1->SetGrid(0,0);
                c1->SetTitle(0);
                c1->SetBottomMargin(0.15);
                c1->SetLeftMargin(0.15);
                PtCos->Draw();
                Func_rho->Draw("same");
                c1->SaveAs(Title->Data());
                delete c1;

                Float_t real_rho = Func_rho->GetParameter(1);
                Float_t real_rho_error = Func_rho->GetParError(1);
                Float_t weight = PtCos->Integral(1,7);

                if(PtBin>=3) {
                  weight_rho00 += real_rho*weight;
                  weight_error_rho00 += real_rho_error*real_rho_error*weight*weight;
                  weight_all += weight;
                }

                pt_rho00[PtBin-1] += real_rho*weight;
                pt_error_rho00[PtBin-1] += real_rho_error*real_rho_error*weight*weight;
                pt_all[PtBin-1] += weight;

              
                if(i_cent == 9 && PtBin == 3 && i_dca == 0 && i_sig == 0 && i_norm == 0 && i_sigma == 0 && i_method == 1)
                { 
                   string name = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d",PtBin-1,i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F);                             
                   cout << name.c_str() << endl; 
                   clonePt = new TH1F(name.c_str(),name.c_str(),7,0,1);
                   for(int i=1; i<=7; i++) {
                     float inte_mean = PtCos->GetBinContent(i);
                     float inte_mean_error = PtCos->GetBinError(i);
                     clonePt->SetBinContent(i, inte_mean);
                     clonePt->SetBinError(i, inte_mean_error);
                   }
                  
                }
                delete PtCos;
                delete PtCos_raw;
              }

              weight_rho00 = weight_rho00/weight_all;
              weight_error_rho00 = TMath::Sqrt(weight_error_rho00)/weight_all;

              for(int PtBin=1; PtBin<=6; PtBin++) {
                pt_rho00[PtBin-1] = pt_rho00[PtBin-1]/pt_all[PtBin-1];
                pt_error_rho00[PtBin-1] = TMath::Sqrt(pt_error_rho00[PtBin-1])/pt_all[PtBin-1];
              }

              cout<<"rho00"<<endl;
              cout<<weight_rho00<<" +/- "<<weight_error_rho00<<endl;

              for(int PtBin=1; PtBin<=6; PtBin++) {
                cout<<pt_rho00[PtBin-1];
                if(PtBin==6) cout<<" "<<weight_rho00<<endl;
                else cout<<" ";
              }

              for(int PtBin=1; PtBin<=6; PtBin++) {
                cout<<pt_error_rho00[PtBin-1];
                if(PtBin==6) cout<<" "<<weight_error_rho00<<endl;
                else cout<<" ";
              }
 
              Title = new TString(Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F));                             

              g_rho00[i_dca][i_sig][i_norm][i_sigma][i_method] = new TGraphAsymmErrors();
              g_rho00[i_dca][i_sig][i_norm][i_sigma][i_method]->SetName(Title->Data());
              for(Int_t PtBin=1; PtBin<=6; PtBin++) {
                double ptmean = (ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                g_rho00[i_dca][i_sig][i_norm][i_sigma][i_method]->SetPoint(PtBin-1,ptmean,pt_rho00[PtBin-1]);
                g_rho00[i_dca][i_sig][i_norm][i_sigma][i_method]->SetPointError(PtBin-1,0.0,0.0,pt_error_rho00[PtBin-1],pt_error_rho00[PtBin-1]);
              }
              delete Title;

              Title = new TString(Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F));                             

              rho00_hist[i_dca][i_sig][i_norm][i_sigma][i_method] = new TH1D(Title->Data(), Title->Data(), 7, 0.5, 7.5);
              rho00_hist[i_dca][i_sig][i_norm][i_sigma][i_method]->SetBinContent(7,weight_rho00);
              rho00_hist[i_dca][i_sig][i_norm][i_sigma][i_method]->SetBinError(7,weight_error_rho00);
              for(Int_t PtBin=1; PtBin<=6; PtBin++) {
                rho00_hist[i_dca][i_sig][i_norm][i_sigma][i_method]->SetBinContent(PtBin,pt_rho00[PtBin-1]);
                rho00_hist[i_dca][i_sig][i_norm][i_sigma][i_method]->SetBinError(PtBin,pt_error_rho00[PtBin-1]);
              }

              g_rho00[i_dca][i_sig][i_norm][i_sigma][i_method]->Write();
              rho00_hist[i_dca][i_sig][i_norm][i_sigma][i_method]->Write();
              cout << g_rho00[i_dca][i_sig][i_norm][i_sigma][i_method]->GetN() << endl;;
              delete Title;
              //delete rho00_hist;
              //delete g_rho00;
            }
            }
          }
        }
      }
    }
  }
  
  clonePt->Write();
  output->Close();
  }
}

void Correction(Double_t Res, Double_t Res_error, Double_t obv_rho, Double_t obv_rho_error, Double_t &real_rho, Double_t &real_rho_error) {

  real_rho = (Res-1+4*obv_rho)/(1+3*Res);

  Double_t real_rho_error_1 = obv_rho_error*4/(1+3*Res);
  Double_t real_rho_error_2 = Res_error*(4-12*obv_rho)/(1+3*Res)/(1+3*Res);
  real_rho_error = TMath::Sqrt(real_rho_error_1*real_rho_error_1+real_rho_error_2*real_rho_error_2);

}

Double_t chi(Double_t res) {

  Double_t chi = 2.0;
  Double_t delta = 1.0;

  for(Int_t i=0; i<15; i++) {
    chi = (res1(chi) < res) ? chi + delta : chi - delta;
    delta = delta/2.;
  }

  return chi;
}

Double_t res1(Double_t chi) {

  Double_t con = 0.626657;                   // TMath::Sqrt(pi/2)/2
  Double_t arg = chi * chi / 4.;

  Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

  return res;
}

Double_t resEventPlane(Double_t chi) {

  Double_t con = 0.626657;                   // TMath::Sqrt(pi/2)/2
  Double_t arg = chi * chi / 4.;
  Double_t halfpi = TMath::Pi()/2.;

  Double_t besselOneHalf = TMath::Sqrt(arg/halfpi) * TMath::SinH(arg)/arg;
  Double_t besselThreeHalfs = TMath::Sqrt(arg/halfpi) * (TMath::CosH(arg)/arg - TMath::SinH(arg)/(arg*arg));
  Double_t res = con * chi * exp(-arg) * (besselOneHalf + besselThreeHalfs);

  return res;
}

double FuncAD(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double D = par[2];
  double R = par[3];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double result = (1.+Bs*D/2.) + (As+D)*CosTheta*CosTheta + (As*D-Bs*D/2.)*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}

