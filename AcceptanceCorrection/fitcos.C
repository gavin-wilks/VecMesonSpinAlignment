void fitcos(int eff_mode=0, int dca_mode=0, int sigma_mode=1, int norm_mode=0, int inte_mode=-1,int plane_mode=1) {

  gROOT->Reset();
  gStyle->SetOptDate(0);
//  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(22);
  gStyle->SetLabelFont(22);
  gStyle->SetTitleFont(22);

  Double_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};
  Double_t pt_set[8] = {0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4, 7.2};
//  Double_t pt_set[8] = {3.0, 3.3, 3.6, 3.9, 4.3, 4.6, 4.9, 5.2};
//  Double_t pt_set[8] = {0.6, 1.4, 2.2, 3.0, 3.8, 4.6, 5.4, 7.2};
  Char_t Char[10][10] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%","20-60%"};

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
  TH1D *resolution_1 = new TH1D("resolution_1","resolution_1", 10, -0.5, 9.5);
  TH1D *resolution_12 = new TH1D("resolution_12","resolution_12", 10, -0.5, 9.5);
  TH1D *resolution_2 = new TH1D("resolution_2","resolution_2", 10, -0.5, 9.5);

  for(Int_t cent=0; cent<=9; cent++) {

    Double_t CosMean = DeltaPsi1->GetBinContent(cent+1);
    Double_t CosError = DeltaPsi1->GetBinError(cent+1);
    Double_t ZDCSMD_resSub = CosMean>0.? (TMath::Sqrt(CosMean)) : 0.;
    Double_t ZDCSMD_resSubErr = CosMean>0.? CosError/2./CosMean : 0.;
    Double_t ZDCSMD_chiSub = chi(ZDCSMD_resSub);
    Double_t ZDCSMD_chiSubDelta = chi(ZDCSMD_resSub+0.005);
    Double_t ResZ = resEventPlane(TMath::Sqrt(2.)*ZDCSMD_chiSub);
    Double_t ResZDelta = resEventPlane(TMath::Sqrt(2.)*ZDCSMD_chiSubDelta);
    Double_t ResErrZ = ZDCSMD_resSubErr * TMath::Abs((ResZ-ResZDelta)/0.005);

    if(cent<9) cout<<"Cent"<<cent_set[cent+1]<<"_"<<cent_set[cent]<<": "<<ResZ<<" +/- "<<ResErrZ<<endl;
    else cout<<"20-60%: "<<ResZ<<" +/- "<<ResErrZ<<endl;

    resolution_1->SetBinContent(cent+1, ResZ);
    resolution_1->SetBinError(cent+1, ResErrZ);

    Double_t Cos12 = DeltaPsi12_2->GetBinContent(cent+1)/TMath::Sqrt(2);
    Double_t CosError12 = DeltaPsi12_2->GetBinError(cent+1)/TMath::Sqrt(2);
    Double_t mean = Cos12/ResZ;
    Double_t error = ResErrZ*ResErrZ/ResZ/ResZ + CosError12*CosError12/Cos12/Cos12;

    resolution_12->SetBinContent(cent+1,mean);
    resolution_12->SetBinError(cent+1, mean*sqrt(error));

    if(cent<9) cout<<"Cent"<<cent_set[cent+1]<<"_"<<cent_set[cent]<<": "<<mean<<" +/- "<<error<<endl;
    else cout<<"20-60%: "<<mean<<" +/- "<<error<<endl;

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

  }

//  TFile *input = new TFile("massfile/yield012000.root","READ");
  TFile *input = new TFile("yield_pt.root","READ");

  TF1 *Func_rho = new TF1("Func_rho",FuncAD,0,1,4);
  TF1 *Func_rdl = new TF1("Func_rdl","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);

  TH1D *Big_yield = new TH1D("Big_yield","Big_yield", 7, 0, 1);
  double big_cos[9][2] = {0};

  TString *Title;

  Double_t weight_rho00[7] = {0};
  Double_t weight_error_rho00[7] = {0};
  Double_t weight_all[7] = {0};

  for(Int_t cent=5; cent<=5; cent++) {

    Title = new TString("yield_Cent");
    *Title += cent_set[cent+1];
    *Title += "_";
    *Title += cent_set[cent];
    cout<<"Read: "<<Title->Data()<<endl;
    TH2D *PtCos2D = (TH2D*)input->Get(Title->Data());
    delete Title;

    for(Int_t PtBin=1; PtBin<=6; PtBin++) {

      Double_t Res_1 = resolution_1->GetBinContent(10);
      Double_t Res_2 = resolution_2->GetBinContent(10);
      Double_t Res_12 = resolution_12->GetBinContent(10);

      TH1D *PtCos = PtCos2D->ProjectionY("PtCos",PtBin,PtBin);

      PtCos->GetXaxis()->SetTitleOffset(1.2);
      PtCos->GetXaxis()->SetTitle("cos#theta*");
      PtCos->GetXaxis()->SetLabelSize(0.05);
      PtCos->GetXaxis()->SetLabelFont();
      PtCos->GetXaxis()->SetTitleSize(0.05);
      PtCos->GetXaxis()->SetTitleOffset(1.2);
      PtCos->GetXaxis()->SetNdivisions(506);
      PtCos->GetYaxis()->SetRangeUser(PtCos->GetMinimum()/1.02, PtCos->GetMaximum()*1.02);
      PtCos->GetYaxis()->SetTitle("Yield");
      PtCos->GetYaxis()->SetLabelSize(0.05);
      PtCos->GetYaxis()->SetLabelFont();
      PtCos->GetYaxis()->SetTitleSize(0.05);
      PtCos->GetYaxis()->SetTitleOffset(1.2);
      PtCos->GetYaxis()->SetNdivisions(506);
      PtCos->SetMarkerColor(2);
      PtCos->SetMarkerSize(1.8);
      PtCos->SetMarkerStyle(21);
      PtCos->Draw();
      Func_rho->SetParameter(0, PtCos->GetBinContent(5));
      Func_rho->SetParameter(1, 0.33);
      if(plane_mode==0) Func_rho->FixParameter(2, 0);
      else Func_rho->FixParameter(2, eta_D[PtBin-1]);
      if(plane_mode==0) Func_rho->FixParameter(3, 1);
      if(plane_mode==1) Func_rho->FixParameter(3, Res_1);
      if(plane_mode==2) Func_rho->FixParameter(3, Res_12);
//      Func_rho->FixParameter(2, 0);
//      Func_rho->FixParameter(3, 1);
      cout<<"Fit with real EP:"<<endl;
      PtCos->Fit(Func_rho, "NQMI");
      Func_rho->Draw("same");

      if(PtBin>=2) {
cout<<PtBin<<endl;
cout<<Func_rho->GetChisquare()<<"/"<<Func_rho->GetNDF()<<endl;
}

      Title = new TString("fit/yield_cent");
      *Title += cent;
      *Title += "_pt";
      *Title += PtBin;
      *Title += ".eps";
      c1->SaveAs(Title->Data());

      if(PtBin==2) {
        TFile *fig2 = new TFile("fig2.root","RECREATE");
        PtCos->Write();
        Func_rho->Write();
        fig2->Close();
        c1->SaveAs("fig2.eps");
      }

      Double_t real_rho = Func_rho->GetParameter(1);
      Double_t real_rho_error = Func_rho->GetParError(1);
      Double_t weight = PtCos->Integral(1,7);

      weight_rho00[PtBin] += real_rho*weight;
      weight_error_rho00[PtBin] += real_rho_error*real_rho_error*weight*weight;
      weight_all[PtBin] += weight;

      if(PtBin>=2) {
        weight_rho00[0] += real_rho*weight;
        weight_error_rho00[0] += real_rho_error*real_rho_error*weight*weight;
        weight_all[0] += weight;

      if(plane_mode==1) TFile *note_fig = new TFile("fig_1stEP_200GeV.root","UPDATE");
      else TFile *note_fig = new TFile("fig_2ndEP_200GeV.root","UPDATE");
        Title = new TString("fit_cos_pt");
        *Title += PtBin;
        PtCos->SetName(Title->Data());
        PtCos->SetTitle(Title->Data());
        PtCos->Write();
        *Title += "_func";
        Func_rho->SetName(Title->Data());
        Func_rho->SetTitle(Title->Data());
        Func_rho->Write();
        note_fig->Close();
      }

      delete PtCos;
    }

    delete PtCos2D;

  }

  for(Int_t PtBin=0; PtBin<=6; PtBin++) {
    weight_rho00[PtBin] = weight_rho00[PtBin]/weight_all[PtBin];
    weight_error_rho00[PtBin] = sqrt(weight_error_rho00[PtBin])/weight_all[PtBin];
  }

  for(Int_t PtBin=1; PtBin<=6; PtBin++) {
    cout<<weight_rho00[PtBin]<<" ";
  }
  cout<<weight_rho00[0]<<endl;
  for(Int_t PtBin=1; PtBin<=6; PtBin++) {
    cout<<weight_error_rho00[PtBin]<<" ";
  }
  cout<<weight_error_rho00[0]<<endl;

  cout<<weight_rho00[0]<<" +/-"<<weight_error_rho00[0]<<endl;

  TFile *output = new TFile("rho00_27GeV.root","UPDATE");
  Double_t pt_rho00[7] = {0, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};

  Title = new TString("EP_");
  *Title += plane_mode;
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

  TH1D *rho00_hist = new TH1D(Title->Data(), Title->Data(), 6, pt_rho00);
  rho00_hist->SetBinContent(1,weight_rho00[0]);
  rho00_hist->SetBinError(1,weight_error_rho00[0]);
  for(Int_t PtBin=2; PtBin<=6; PtBin++) {
    rho00_hist->SetBinContent(PtBin,weight_rho00[PtBin]);
    rho00_hist->SetBinError(PtBin,weight_error_rho00[PtBin]);
  }
  output->Write();

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

