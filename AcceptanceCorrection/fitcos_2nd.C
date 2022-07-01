void fitcos_2nd(int eff_mode=0, int dca_mode=0, int sigma_mode=1, int norm_mode=0, int inte_mode=-1) {

  gROOT->Reset();
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(22);
  gStyle->SetLabelFont(22);
  gStyle->SetTitleFont(22);

  Double_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};
  Double_t pt_set[8] = {0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4, 7.2};
//  Double_t pt_set[8] = {3.0, 3.3, 3.6, 3.9, 4.3, 4.6, 4.9, 5.2};
//  Double_t pt_set[8] = {0.6, 1.4, 2.2, 3.0, 3.8, 4.6, 5.4, 7.2};
  Char_t Char[9][10] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};

//  double eta_D[6] = {0.0262611, 0.0181768, 0.0107771, 0.00359328, 0.00353325, 0.00185303}; //200GeV
//  double eta_D[6] = {-0.0092467, -0.00955624, -0.0083202, -0.00649984, -0.0054873, -0.00480371};//11GeV
//  double eta_D[6] = {-0.00989006, -0.0102287, -0.0102549, -0.00800524, -0.00767692, -0.00384473};//19GeV
  double eta_D[6] = {-0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006, -0.00539932};//27GeV
//  double eta_D[6] = {-0.0092467, -0.00955624, -0.0083202, -0.00649984, -0.0054873, -0.00480371};//39GeV
//  double eta_D[6] = {0.00149495,0.00145323, 0.000128829, 0.000765007, -0.000433217, -0.000407656}; //62GeV

  TCanvas *c1 = new TCanvas();
  c1->SetFillColor(0);
  c1->SetGrid(0,0);
  c1->SetTitle(0);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.15);

  TF1 *line = new TF1("line","1/3",-0.5,8.5);
  line->SetLineStyle(2);
  line->SetLineWidth(1);

  TFile event("all_event_output.root");
  TH1D *resolution_1 = new TH1D("resolution_1","resolution_1", 9, -0.5, 8.5);
  TH1D *resolution_12 = new TH1D("resolution_12","resolution_12", 9, -0.5, 8.5);
  TH1D *resolution_2 = new TH1D("resolution_2","resolution_2", 9, -0.5, 8.5);

  for(Int_t cent=0; cent<9; cent++) {

    Double_t CosMean = DeltaPsi1->GetBinContent(cent+1);
    Double_t CosError = DeltaPsi1->GetBinError(cent+1);
    Double_t ZDCSMD_resSub = CosMean>0.? (TMath::Sqrt(CosMean)) : 0.;
    Double_t ZDCSMD_resSubErr = CosMean>0.? CosError/2./CosMean : 0.;
    Double_t ZDCSMD_chiSub = chi(ZDCSMD_resSub);
    Double_t ZDCSMD_chiSubDelta = chi(ZDCSMD_resSub+0.005);
    Double_t ResZ = resEventPlane(TMath::Sqrt(2.)*ZDCSMD_chiSub);
    Double_t ResZDelta = resEventPlane(TMath::Sqrt(2.)*ZDCSMD_chiSubDelta);
    Double_t ResErrZ = ZDCSMD_resSubErr * TMath::Abs((ResZ-ResZDelta)/0.005);

    resolution_1->SetBinContent(cent+1, ResZ);
    resolution_1->SetBinError(cent+1, ResErrZ);

    Double_t Cos12 = DeltaPsi12_2->GetBinContent(cent+1)/sqrt(2);
    Double_t CosError12 = DeltaPsi12_2->GetBinError(cent+1)/sqrt(2);
    Double_t mean = Cos12/ResZ;
    Double_t error = ResErrZ*ResErrZ/ResZ/ResZ + CosError12*CosError12/Cos12/Cos12;

    resolution_12->SetBinContent(cent+1,mean);
    resolution_12->SetBinError(cent+1, mean*sqrt(error));

    Double_t CosMean_2 = DeltaPsi2_2->GetBinContent(cent+1);
    Double_t CosError_2 = DeltaPsi2_2->GetBinError(cent+1);
    Double_t ZDCSMD_resSub_2 = CosMean_2>0.? (TMath::Sqrt(CosMean_2)) : 0.;
    Double_t ZDCSMD_resSubErr_2 = CosMean_2>0.? CosError_2/2./CosMean_2 : 0.;
    Double_t ZDCSMD_chiSub_2 = chi(ZDCSMD_resSub_2);
    Double_t ZDCSMD_chiSubDelta_2 = chi(ZDCSMD_resSub_2+0.005);
    Double_t ResZ_2 = res1(TMath::Sqrt(2.)*ZDCSMD_chiSub_2);
    Double_t ResZDelta_2 = res1(TMath::Sqrt(2.)*ZDCSMD_chiSubDelta_2);
    Double_t ResErrZ_2 = ZDCSMD_resSubErr_2 * TMath::Abs((ResZ_2-ResZDelta_2)/0.005);

    resolution_2->SetBinContent(cent+1,ResZ_2);
    resolution_2->SetBinError(cent+1,ResErrZ_2);

    cout<<"Cent"<<cent_set[cent+1]<<"_"<<cent_set[cent]<<": "<<ResZ<<" +/- "<<ResErrZ<<" "<<ResZ_2<<" +/- "<<ResErrZ_2<<" "<<mean<<" +/- "<<mean*sqrt(error)<<endl;

  }

  double Res_12 = 0;
  double Res_12_weight = 0;
  for(int cent=2; cent<=5; cent++) {
    Res_12 += resolution_12->GetBinContent(cent+1)/resolution_12->GetBinError(cent+1)/resolution_12->GetBinError(cent+1);
    Res_12_weight += 1./resolution_12->GetBinError(cent+1)/resolution_12->GetBinError(cent+1);
  }
  Res_12 = Res_12/Res_12_weight;

  double eff[9][7][7];
  double eff_error[9][7][7];
  if(eff_mode==0)
    TFile *eff_file = new TFile("Eff_27GeV_SingleKaon_second.root","READ");
  else
    TFile *eff_file = new TFile("Eff_27GeV_SingleKaon_eta.root","READ");

  for(int i=0; i<9;i++)
    for(int j=1; j<=5; j++)
      for(int k=1; k<=7; k++)
      {
        if(j==7) {
          eff[i][j-1][k-1] = 1;
          eff_error[i][j-1][k-1] = 0;
          continue;
        }

        Title = new TString("h_mEffCos_Cent_");
        *Title += 9;
        *Title += "_Pt_";
        *Title += j;
        TH1D *eff_hist = (TH1D*)eff_file->Get(Title->Data());
        eff[i][j-1][k-1] = eff_hist->GetBinContent(k);
        eff_error[i][j-1][k-1] = eff_hist->GetBinError(k)/eff[i][j-1][k-1];
        delete Title;
        delete eff_hist;

      }
  eff_file->Close();

  TFile *input = new TFile("RawRhoPtSys.root","READ");

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

    for(Int_t PtBin=1; PtBin<=5; PtBin++) {

      Title = new TString("pt_");
      *Title += PtBin;
      *Title += "_Centrality_9_2nd_Dca_";
      *Title += dca_mode;
      *Title += "_Sig_";
      *Title += sigma_mode;
      *Title += "_Phi_Norm_";
      *Title += norm_mode;
      *Title += "_Sigma_";
      if(inte_mode<0) *Title += "2_Inte";
      else {
        *Title += inte_mode;
      *Title += "_Count";
      }

      TH1F *PtCos_raw = (TH1F*)input->Get(Title->Data());
      delete Title;
      TH1F *PtCos = new TH1F("PtCos","PtCos", 7, 0, 1);
      for(int i=1; i<=7; i++) {
        float inte_mean = PtCos_raw->GetBinContent(i);
        float inte_mean_error = PtCos_raw->GetBinError(i);
        PtCos->SetBinContent(i, inte_mean/eff[i][PtBin-1][i-1]);
        PtCos->SetBinError(i, inte_mean/eff[i][PtBin-1][i-1]*sqrt(inte_mean_error*inte_mean_error/inte_mean/inte_mean+eff_error[i][PtBin-1][i-1]*eff_error[i][PtBin-1][i-1]));
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
      Func_rho->FixParameter(2, eta_D[PtBin-1]);
      Func_rho->FixParameter(3, Res_12);
//      Func_rho->FixParameter(2, 0);
//      Func_rho->FixParameter(3, 1);
      cout<<"Fit with real EP:"<<endl;
      PtCos->Fit(Func_rho, "NMI");
      Func_rho->Draw("same");

      Title = new TString("fit/yield_cent");
      *Title += cent;
      *Title += "_pt";
      *Title += PtBin;
      *Title += ".eps";
      c1->SaveAs(Title->Data());

      Float_t real_rho = Func_rho->GetParameter(1);
      Float_t real_rho_error = Func_rho->GetParError(1);
      Float_t weight = PtCos->Integral(1,7);

      if(PtBin>=2) {
        weight_rho00 += real_rho*weight;
        weight_error_rho00 += real_rho_error*real_rho_error*weight*weight;
        weight_all += weight;
      }

      pt_rho00[PtBin-1] += real_rho*weight;
      pt_error_rho00[PtBin-1] += real_rho_error*real_rho_error*weight*weight;
      pt_all[PtBin-1] += weight;

      delete PtCos;
      delete PtCos_raw;
    }

  weight_rho00 = weight_rho00/weight_all;
  weight_error_rho00 = sqrt(weight_error_rho00)/weight_all;

  for(int PtBin=1; PtBin<=5; PtBin++) {
    pt_rho00[PtBin-1] = pt_rho00[PtBin-1]/pt_all[PtBin-1];
    pt_error_rho00[PtBin-1] = sqrt(pt_error_rho00[PtBin-1])/pt_all[PtBin-1];
  }

  cout<<"rho00"<<endl;
  cout<<weight_rho00<<" +/- "<<weight_error_rho00<<endl;

  for(int PtBin=1; PtBin<=5; PtBin++) {
    cout<<pt_rho00[PtBin-1];
    if(PtBin==5) cout<<" "<<weight_rho00<<endl;
    else cout<<" ";
  }

  for(int PtBin=1; PtBin<=5; PtBin++) {
    cout<<pt_error_rho00[PtBin-1];
    if(PtBin==5) cout<<" "<<weight_error_rho00<<endl;
    else cout<<" ";
  }

  TFile *output = new TFile("rho00_27GeV.root","UPDATE");
  Double_t ptbin_rho00[7] = {0,1.2, 1.8, 2.4, 3.0, 4.2, 5.4};

  Title = new TString("EP_");
  *Title += "2";
  *Title += "_eff_";
  *Title += eff_mode;
  *Title += "_Dca_";
  *Title += dca_mode;
  *Title += "_Sig_";
  *Title += sigma_mode;
  *Title += "_Phi_Norm_";
  *Title += norm_mode;
  *Title += "_Sigma_";
  if(inte_mode<0) *Title += "2_Inte";
  else {
    *Title += inte_mode;
  *Title += "_Count";
  }

  TH1D *rho00_hist = new TH1D(Title->Data(), Title->Data(), 6, ptbin_rho00);
  rho00_hist->SetBinContent(1,weight_rho00);
  rho00_hist->SetBinError(1,weight_error_rho00);
  for(Int_t PtBin=2; PtBin<=5; PtBin++) {
    rho00_hist->SetBinContent(PtBin,pt_rho00[PtBin-1]);
    rho00_hist->SetBinError(PtBin,pt_error_rho00[PtBin-1]);
  }
//  output->Write();

}

void Correction(Double_t Res, Double_t Res_error, Double_t obv_rho, Double_t obv_rho_error, Double_t &real_rho, Double_t &real_rho_error) {

  real_rho = (Res-1+4*obv_rho)/(1+3*Res);

  Double_t real_rho_error_1 = obv_rho_error*4/(1+3*Res);
  Double_t real_rho_error_2 = Res_error*(4-12*obv_rho)/(1+3*Res)/(1+3*Res);
  real_rho_error = sqrt(real_rho_error_1*real_rho_error_1+real_rho_error_2*real_rho_error_2);

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

  Double_t con = 0.626657;                   // sqrt(pi/2)/2
  Double_t arg = chi * chi / 4.;

  Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

  return res;
}

Double_t resEventPlane(Double_t chi) {

  Double_t con = 0.626657;                   // sqrt(pi/2)/2
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

