#include "Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"

void acceptanceQA_4th_Res(const int energy = 4, const int pid = 0, bool doall = true, std::string res = "0p1", double resval = 0.1, int cut = 1, int method = 1) {


  std::string cutoption = "";
  if(cut == 1) cutoption = "_EtaCut";

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

  TH1D *h_theta_star_before[5];
  TH1D *h_theta[5];
  TH1D *h_theta_before[5];
  TH1D *h_theta_star[5];
  TH1D *h_theta_star_RP[5];
  TH1D *h_theta_star_before_RP[5];
  TH1D *h_out1[5];
  TH1D *h_out2[5];

  TFile *MCFiles[5]; 
  double rho00in[5] = {0.2667,0.30,0.3333,0.3667,0.40};
  int rho00val[5] = {2667,3000,3333,3667,4000};
  TGraphAsymmErrors *g_rho00outBefore = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rho00outAfter  = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rho00outBeforeIn = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rho00outAfterIn  = new TGraphAsymmErrors();
  double Fval[5] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Ferr[5] = {0.0};
  double Gval[5] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Gerr[5] = {0.0};
  double Hval[5] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Herr[5] = {0.0};
  double Fpval[5] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Fperr[5] = {0.0};
  double Gpval[5] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Gperr[5] = {0.0};

  double FvalOrig[5] = {0.0};
  double FerrOrig[5] = {0.0};

  double rhoobs[5] = {0.0};
  double rhoreal[5] = {0.0};
  double rhorealErr[5] = {0.0};

  for(int irho = 0; irho < 5; irho++)
  {
    //MCFiles[irho] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/VaryRho/McAcceptanceOutput_pt3_energy%d_pid%d_cent9_rho%d.root",vmsa::mPID[pid].c_str(),energy,pid,irho),"READ");
    //MCFiles[irho] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/ResolutionTesting/yabs1_%s_EtaCut_noDelta_PsiWrapping_2D/McAcceptanceOutput_pt1_energy%d_pid%d_cent5_EtaMode_0_nrho%d.root",vmsa::mPID[pid].c_str(),res.c_str(),energy,pid,rho00val[irho]),"READ");
    //MCFiles[irho] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/ResolutionTesting/0p8y1_%s_EtaCut_noDelta/McAcceptanceOutput_pt1_energy%d_pid%d_cent5_EtaMode_0_nrho%d.root",vmsa::mPID[pid].c_str(),res.c_str(),energy,pid,rho00val[irho]),"READ");
    MCFiles[irho] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/ResolutionTesting/yabs1_0p1_EtaCut_nov2/McAcceptanceOutput_pt1_energy%d_pid%d_cent5_EtaMode_0_nrho%d.root",vmsa::mPID[pid].c_str(),energy,pid,rho00val[irho]),"READ");
    //MCFiles[irho] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/ResolutionTesting/yabs1_0p1_EtaCut_noDelta/McAcceptanceOutput_pt1_energy%d_pid%d_cent5_EtaMode_0_nrho%d.root",vmsa::mPID[pid].c_str(),energy,pid,rho00val[irho]),"READ");
    //if(cut == 1 && method == 2)  Func_AFG->FixParameter(5, 1.0);//Res_12);
    //MCFiles[irho] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/VaryRho_0p8y1/McAcceptanceOutput_pt1_energy%d_pid%d_cent5_rho%d.root",vmsa::mPID[pid].c_str(),energy,pid,irho),"READ");

    h_theta_star_before[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_theta_star_before"))->Clone(Form("h_theta_star_before_%d",irho));
    h_theta[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_theta"))->Clone(Form("h_theta_%d",irho));
    h_theta_before[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_theta_before"))->Clone(Form("h_theta_before_%d",irho));
    h_theta_star[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_theta_star"))->Clone(Form("h_theta_star_%d",irho));
    h_theta_star_RP[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_theta_star_RP"))->Clone(Form("h_theta_star_RP_%d",irho));
    h_theta_star_before_RP[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_theta_star_before_RP"))->Clone(Form("h_theta_star_before_RP_%d",irho));
    //h_out1[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_out1"))->Clone(Form("h_out1_%d",irho));
    //h_out2[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_out2"))->Clone(Form("h_out2_%d",irho));
 
    TF1 *Func_rho = new TF1("Func_rho","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
    TF1 *Func_A = new TF1("Func_A",Func4th,0,1,3);
    Func_A->SetLineColor(kBlue);
    TF1 *Func_Orig = new TF1("Func_Orig","[0]*(1 + [1]*x*x)",0,1);
    //Func_Orig->SetLineColor(kRed);
    TF1 *Func_Theta = new TF1("Func_Theta","[0]*(1 + [1]*x*x + [2]*x*x*x*x + [3]*x*x*x*x*x*x)",0,1);
    //TF1 *Func_AFG = new TF1("Func_AFG",FuncAFG,0,1,3);

    Func_rho->SetParameter(0,h_theta_star[irho]->GetBinContent(1));
    Func_rho->SetParameter(1,1./3.);
    h_theta_star[irho]->Fit(Func_rho,"NERQ");
    rhoobs[irho] = Func_rho->GetParameter(1);

    Func_rho->SetParameter(0,h_theta_star_RP[irho]->GetBinContent(1));
    Func_rho->SetParameter(1,1./3.);
    h_theta_star_RP[irho]->Fit(Func_rho,"NERQ");
    rhoreal[irho] = Func_rho->GetParameter(1);
    rhorealErr[irho] = Func_rho->GetParError(1);
    if(cut == 0) g_rho00outBefore->SetPoint(irho,rhoobs[irho],rhoreal[irho]);
    if(cut == 0) g_rho00outBefore->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));
    if(cut == 0) g_rho00outBeforeIn->SetPoint(irho,rho00in[irho],rhoreal[irho]);
    if(cut == 0) g_rho00outBeforeIn->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));


    TH1D *h_theta_clone = (TH1D*)h_theta[irho]->Clone("h_theta_clone");
    h_theta_clone->Sumw2();
    h_theta_before[irho]->Sumw2();
    h_theta_clone->Divide(h_theta_before[irho]);
    Func_Theta->SetParameter(0,h_theta_clone->GetBinContent(1));
    Func_Theta->SetParameter(1,0);
    h_theta_clone->Fit(Func_Theta,"ER");


    if(cut == 1 && method == 2){
      Fval[irho] = Func_Theta->GetParameter(1);
      Ferr[irho] = Func_Theta->GetParError(1);
      Gval[irho] = Func_Theta->GetParameter(2);
      Gerr[irho] = Func_Theta->GetParError(2);
      Hval[irho] = Func_Theta->GetParameter(3);
      Herr[irho] = Func_Theta->GetParError(3);
    }
    TCanvas *c0 = new TCanvas();
    c0->SetFillColor(0);
    c0->SetGrid(0,0);
    c0->SetTitle(0);
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    h_theta_clone->Draw("pE");
    Func_Theta->Draw("same"); 
    auto legend0 = new TLegend(0.1,0.7,0.48,0.9); 
    legend0->AddEntry(Func_Theta,"4th Order Fit","l");
    //legend->AddEntry(Func_Orig,"2nd Order Fit","l");
    legend0->Draw("same");
    c0->SaveAs(Form("FValues/yabs1%s/pt1_%d%s_method%d_FROMTHETA.pdf",cutoption.c_str(),irho,cutoption.c_str(),method));
    delete c0;


    TH1D *h_theta_star_clone = (TH1D*)h_theta_star[irho]->Clone("h_theta_star_clone");
    h_theta_star_clone->Sumw2();
    h_theta_star_before[irho]->Sumw2();
    h_theta_star_clone->Divide(h_theta_star_before[irho]); 
    Func_A->SetParameter(0,h_theta_star_clone->GetBinContent(1));
    Func_A->SetParameter(1,0);
    Func_A->SetParameter(2,0.0);
    //Func_A->FixParameter(2,0.0);
    Func_Orig->SetParameter(0,h_theta_star_clone->GetBinContent(1));
    Func_Orig->SetParameter(1,0);
    h_theta_star_clone->Fit(Func_Orig,"ER");

    if(cut == 1 && method == 0){
      h_theta_star_clone->Fit(Func_A,"ER");
      Fval[irho] = Func_A->GetParameter(1);
      Ferr[irho] = Func_A->GetParError(1);
      Gval[irho] = Func_A->GetParameter(2);
      Gerr[irho] = Func_A->GetParError(2);
      TCanvas *c1 = new TCanvas();
      c1->SetFillColor(0);
      c1->SetGrid(0,0);
      c1->SetTitle(0);
      c1->SetBottomMargin(0.15);
      c1->SetLeftMargin(0.15);
      h_theta_star_clone->Draw("pE");
      Func_A->Draw("same"); 
      auto legend = new TLegend(0.1,0.7,0.48,0.9); 
      legend->AddEntry(Func_A,"4th Order Fit","l");
      //legend->AddEntry(Func_Orig,"2nd Order Fit","l");
      legend->Draw("same");
      c1->SaveAs(Form("FValues/0p8y1/pt1_%d%s_method%d.pdf",irho,cutoption.c_str(),method));
      delete c1;
    }
 
    FvalOrig[irho] = Func_Orig->GetParameter(1);
    FerrOrig[irho] = Func_Orig->GetParError(1);

    TH1D *h_theta_star_clone_RP = (TH1D*)h_theta_star_RP[irho]->Clone("h_theta_star_clone_RP");
    h_theta_star_clone_RP->Sumw2();
    h_theta_star_before_RP[irho]->Sumw2();
    h_theta_star_clone_RP->Divide(h_theta_star_before_RP[irho]); 
    Func_A->SetParameter(0,h_theta_star_clone_RP->GetBinContent(1));
    Func_A->SetParameter(1,0);
    Func_A->SetParameter(2,0.0);
    //Func_A->FixParameter(2,0.0);
    Func_Orig->SetParameter(0,h_theta_star_clone_RP->GetBinContent(1));
    Func_Orig->SetParameter(1,0);
    h_theta_star_clone_RP->Fit(Func_Orig,"ER");
    h_theta_star_clone_RP->Fit(Func_A,"ER");

    if(cut == 1 && method == 1){
      Fval[irho] = Func_A->GetParameter(1);
      Ferr[irho] = Func_A->GetParError(1);
      Gval[irho] = Func_A->GetParameter(2);
      Gerr[irho] = Func_A->GetParError(2);
      TCanvas *c1 = new TCanvas();
      c1->SetFillColor(0);
      c1->SetGrid(0,0);
      c1->SetTitle(0);
      c1->SetBottomMargin(0.15);
      c1->SetLeftMargin(0.15);
      h_theta_star_clone_RP->Draw("pE");
      Func_A->Draw("same"); 
      auto legend = new TLegend(0.1,0.7,0.48,0.9); 
      legend->AddEntry(Func_A,"4th Order Fit","l");
      //legend->AddEntry(Func_Orig,"2nd Order Fit","l");
      legend->Draw("same");
      c1->SaveAs(Form("FValues/0p8y1/pt1_%d%s_method%d.pdf",irho,cutoption.c_str(),method));
      delete c1;
    }

    TF1 *Func_AFG = new TF1("Func_AFG",FuncAFG,0,1,5);
    //TF1 *Func_AFG = new TF1("Func_AFG",FuncAFGH,0,1,6);
    Func_AFG->SetParameter(0,h_theta_star_RP[irho]->GetBinContent(1));
    Func_AFG->SetParameter(1,1./3.);
    Func_AFG->FixParameter(2, Func_A->GetParameter(1));//Fval[irho]);
    Func_AFG->FixParameter(3, Func_A->GetParameter(2));//Gval[irho]);
    Func_AFG->FixParameter(4, 1.0);//Res_12);
    if(cut == 1 && method == 2) Func_AFG->FixParameter(2, Fval[irho]);
    if(cut == 1 && method == 2) Func_AFG->FixParameter(3, Gval[irho]);
    if(cut == 1 && method == 2) Func_AFG->FixParameter(4, Hval[irho]);
    if(cut == 1 && method == 2)  Func_AFG->FixParameter(5, 1.0);//Res_12);
    h_theta_star_RP[irho]->Fit(Func_AFG,"NERQ");
    rhoreal[irho] = Func_AFG->GetParameter(1);
    rhorealErr[irho] = Func_AFG->GetParError(1);
    cout << "rho00in = " << rho00in[irho] << ",    rho00 from RP with acceptance correction = " << rhoreal[irho] << endl;
    if(cut == 1) g_rho00outBefore->SetPoint(irho,rhoobs[irho],rhoreal[irho]);
    if(cut == 1) g_rho00outBefore->SetPointError(irho,0.0,0.0,Func_AFG->GetParError(1),Func_AFG->GetParError(1));
    if(cut == 1) g_rho00outBeforeIn->SetPoint(irho,rho00in[irho],rhoreal[irho]);
    if(cut == 1) g_rho00outBeforeIn->SetPointError(irho,0.0,0.0,Func_AFG->GetParError(1),Func_AFG->GetParError(1));


    /*TCanvas *c1 = new TCanvas();
    c1->SetFillColor(0);
    c1->SetGrid(0,0);
    c1->SetTitle(0);
    c1->SetBottomMargin(0.15);
    c1->SetLeftMargin(0.15);
    h_theta_star_clone->Draw("pE");
    Func_A->Draw("same"); 
    auto legend = new TLegend(0.1,0.7,0.48,0.9); 
    legend->AddEntry(Func_A,"4th Order Fit","l");
    legend->AddEntry(Func_Orig,"2nd Order Fit","l");
    legend->Draw("same");
    c1->SaveAs(Form("FValues/0p8y1/pt1_%d%s_method%d.pdf",irho,cutoption.c_str(),method));
    delete c1;
*/
    delete h_theta_clone;
    delete h_theta_star_clone;
    delete Func_rho;
    delete Func_A;
    delete Func_Theta;
    //MCFiles[irho]->Close();
  }
 
  double FfromFp[5] = {0.0};
  for(int irho = 0; irho < 5; irho++)
  {
    //FfromFp[irho] = -Fval[irho]/(2.+Fval[irho]);
    //double err = TMath::Abs((Fval[irho]/(2.-Fval[irho])/(2.-Fval[irho])-1./(2.-Fval[irho]))*Ferr[irho]);
    //cout << "F from F* = " << FfromFp[irho] << " +/- " << err << "      F from Orig Method = " << FvalOrig[irho] << " +/- " << FerrOrig[irho] << endl;
    //cout << "F from theta = " << Fpval[irho] << " +/- " << Fperr[ierr] << "      G from theta = " << Gpval[irho] << " +/- " << Gperr[irho] << endl;
    cout << "F = " << Fval[irho] << " +/- " << Ferr[irho] << endl;
    cout << "G = " << Gval[irho] << " +/- " << Gerr[irho] << endl;
    cout << "H = " << Hval[irho] << " +/- " << Herr[irho] << endl;
  }


  //TCanvas *c1[5];
  TGraphAsymmErrors *g_diff = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_diffI = new TGraphAsymmErrors();

  for(int irho = 0; irho < 5; irho++) 
  {
    //TF1 *Func_rho = new TF1("Func_rho",FuncAFG,0,1,5);
    TF1 *Func_rho = new TF1("Func_rho",FuncAFG,0,1,5);
    Func_rho->SetLineColor(kBlue);
    TF1 *Func_rhoOrig = new TF1("Func_rhoOrig",FuncAD,0,1,4);
    Func_rhoOrig->SetLineColor(kRed);
    //TF1 *Func_rho = new TF1("Func_rho",FuncAF_Updated,0,1,4);
    //TF1 *Func_rdl = new TF1("Func_rdl","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);

    h_theta_star[irho]->GetXaxis()->SetTitleOffset(1.2);
    h_theta_star[irho]->GetXaxis()->SetTitle("cos#theta*");
    h_theta_star[irho]->GetYaxis()->SetTitle("yield");
    h_theta_star[irho]->GetYaxis()->SetTitleOffset(1.0);
    h_theta_star[irho]->SetMarkerColor(2);
    h_theta_star[irho]->SetMarkerSize(1.8);
    h_theta_star[irho]->SetMarkerStyle(21);
    //h_theta_star[irho]->Draw();
    Func_rho->SetParameter(0, h_theta_star[irho]->GetBinContent(5));
    Func_rho->SetParameter(1, 0.333);
    //g_rho00outAfter->SetPoint(irho,rho00in[irho],Func_rho->GetParameter(1));
    //g_rho00outAfter->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));
    Func_rho->FixParameter(2, Fval[2]);
    Func_rho->FixParameter(3, Gval[2]);
    //Func_rho->FixParameter(4, Hval[irho]);
    Func_rho->FixParameter(4, resval);//Res_12);

    //Func_rho->FixParameter(3, 1.0);//Res_12);

    Func_rhoOrig->SetParameter(0, h_theta_star[irho]->GetBinContent(5));
    Func_rhoOrig->SetParameter(1, 0.333);
    Func_rhoOrig->FixParameter(2, FvalOrig[irho]);
    Func_rhoOrig->FixParameter(3, resval);//Res_12);


    //Func_rho->Draw("same");
    //delete Func_rdl;
    TString *Title = new TString(Form("fit/acceptanceQA/yield_pt_1_cent_5_rho%d%s_method%d.pdf",irho,cutoption.c_str(),method));
  
    TCanvas *c1 = new TCanvas(Form("c%d",irho),Form("c%d",irho),10,10,800,800);
    c1->cd();
    c1->SetFillColor(0);
    c1->SetGrid(0,0);
    c1->SetTitle(0);
    c1->SetBottomMargin(0.15);
    c1->SetLeftMargin(0.15);
    h_theta_star[irho]->Draw("pE");
    h_theta_star[irho]->Fit(Func_rhoOrig, "NMI");
    //Func_rhoOrig->Draw("same");
    h_theta_star[irho]->Fit(Func_rho, "NMI");
    Func_rho->Draw("same");
 
    auto legend = new TLegend(0.1,0.7,0.48,0.9); 
    legend->AddEntry(Func_rho,"4th Order Fit","l");
    //legend->AddEntry(Func_rhoOrig,"2nd Order Fit","l");
    legend->Draw("same");
    
    c1->SaveAs(Title->Data());
    //double rho00final = 4./(1.+3.*resval)*(Func_rho->GetParameter(1)-1/3.)+1./3.;
    double rho00final = Func_rho->GetParameter(1);
    cout << "rho after = " << rho00final << " +/- " << Func_rho->GetParError(1) << endl;
    g_rho00outAfter->SetPoint(irho,rhoobs[irho],rho00final);
    g_rho00outAfter->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));
    g_rho00outAfterIn->SetPoint(irho,rho00in[irho],rho00final);
    g_rho00outAfterIn->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));

    cout << "rho input = " << rho00in[irho] << ",    rho observed = " << rhoobs[irho] << endl;

    cout << "rho00in = " << rho00in[irho] << ",    rho00 from RP with acceptance correction = " << rhoreal[irho] << " +/- " << rhorealErr[irho] << endl;
    g_diff->SetPoint(irho,rho00in[irho],rho00final-rhoreal[irho]);
    g_diffI->SetPoint(irho,rho00in[irho],rho00final-rho00in[irho]);
    double uncorrErr = TMath::Sqrt(Func_rho->GetParError(1)*Func_rho->GetParError(1)+rhorealErr[irho]*rhorealErr[irho]);
    double corrErr = TMath::Sqrt(fabs(Func_rho->GetParError(1)*Func_rho->GetParError(1)-rhorealErr[irho]*rhorealErr[irho]));
    std::cout << "uncorrelated Error = " << uncorrErr << ",      correlated Error = " << corrErr << std::endl;
    g_diff->SetPointError(irho,0.0,0.0,corrErr,corrErr);
    g_diffI->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));

    delete c1;  
    delete Func_rho;
    delete Func_rhoOrig;
    delete Title;

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
  h_play->GetYaxis()->SetTitle("#rho_{00}^{reco}");
  h_play->GetXaxis()->SetTitle("#rho_{00}^{obs}");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetLimits(0.1,0.55);
  h_play->GetYaxis()->SetRangeUser(0.1,0.55);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");
  h_play->Draw("pE");

  //TF1 *linear = new TF1("linear","4./(1+3*[0])*(x-1./3.) + 1./3.",0.24,0.42);
  //linear->FixParameter(0,resval);
  //linear->Draw("same");

  g_rho00outBefore->SetMarkerStyle(24);
  g_rho00outBefore->SetMarkerColor(kRed);
  g_rho00outBefore->SetMarkerSize(3.0);
  g_rho00outBefore->SetLineWidth(3);
  g_rho00outBefore->Draw("pE Same");

  g_rho00outAfter->SetMarkerStyle(26);
  g_rho00outAfter->SetMarkerColor(kBlue);
  g_rho00outAfter->SetMarkerSize(3.0);
  g_rho00outAfter->SetLineWidth(3);
  g_rho00outAfter->Draw("pE Same");

  TLegend *leg1 = new TLegend(0.2,0.70,0.5,0.85);
  leg1->SetFillColor(10);
  leg1->SetBorderSize(0);
  leg1->AddEntry(g_rho00outBefore,"Fit to cos(#theta^{*})","p");
  leg1->AddEntry(g_rho00outAfter,"Fit to cos(#theta^{*'}) using R","p");
  //leg1->AddEntry(linear,"#frac{4}{1+3R} (#rho_{00}^{obs} - #frac{1}{3}) + #frac{1}{3}","l");
  leg1->Draw("same");

  string FigureName = Form("./figures/AcceptanceQA/AcceptanceQA_%s%s_method%d.pdf",vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method);
  c_play->SaveAs(FigureName.c_str());

  TCanvas *c_play2 = new TCanvas("c_play2","c_play2",10,10,800,800);
  c_play2->SetLeftMargin(0.15);
  c_play2->SetBottomMargin(0.15);
  c_play2->SetGrid(0,0);
  c_play2->SetTicks(1,1);
  c_play2->cd();

  h_play->GetXaxis()->SetTitle("#rho_{00} input");
  h_play->Draw("pE");

  TF1 *linear = new TF1("linear","x",0.24,0.42);
  //linear->FixParameter(0,resval);
  linear->Draw("same");


  g_rho00outBeforeIn->SetMarkerStyle(24);
  g_rho00outBeforeIn->SetMarkerColor(kRed);
  g_rho00outBeforeIn->SetMarkerSize(3.0);
  g_rho00outBeforeIn->SetLineWidth(3);
  g_rho00outBeforeIn->Draw("pE Same");

  g_rho00outAfterIn->SetMarkerStyle(26);
  g_rho00outAfterIn->SetMarkerColor(kBlue);
  g_rho00outAfterIn->SetMarkerSize(3.0);
  g_rho00outAfterIn->SetLineWidth(3);
  g_rho00outAfterIn->Draw("pE Same");

  TLegend *leg2 = new TLegend(0.2,0.70,0.5,0.85);
  leg2->SetFillColor(10);
  leg2->SetBorderSize(0);
  leg2->AddEntry(g_rho00outBeforeIn,"Fit to cos(#theta^{*})","p");
  leg2->AddEntry(g_rho00outAfterIn,"Fit to cos(#theta^{*'}) using R","p");
  //leg1->AddEntry(linear,"#frac{4}{1+3R} (#rho_{00}^{obs} - #frac{1}{3}) + #frac{1}{3}","l");
  leg2->Draw("same");

  FigureName = Form("./figures/AcceptanceQA/AcceptanceQAInput_%s%s_method%d.pdf",vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method);
  c_play2->SaveAs(FigureName.c_str());

  TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,800,800);
  c_diff->SetLeftMargin(0.15);
  c_diff->SetBottomMargin(0.15);
  c_diff->SetGrid(0,0);
  c_diff->SetTicks(1,1);
  c_diff->cd();

  g_diff->GetHistogram()->GetXaxis()->SetTitle("Input #rho_{00}");
  g_diff->GetHistogram()->GetYaxis()->SetTitle("#rho_{00}'-#rho_{00}");
  g_diff->GetHistogram()->SetTitle("Difference Plot");

  g_diff->SetMarkerStyle(20);
  g_diff->SetMarkerColor(kBlack);
  g_diff->Draw("APE");

  FigureName = Form("./figures/AcceptanceQA/AcceptanceQADifference_%s%s_method%d.pdf",vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method);
  c_diff->SaveAs(FigureName.c_str());

  TCanvas *c_diffI = new TCanvas("c_diffI","c_diffI",10,10,800,800);
  c_diffI->SetLeftMargin(0.15);
  c_diffI->SetBottomMargin(0.15);
  c_diffI->SetGrid(0,0);
  c_diffI->SetTicks(1,1);
  c_diffI->cd();

  g_diffI->GetHistogram()->GetXaxis()->SetTitle("Input #rho_{00}");
  g_diffI->GetHistogram()->GetYaxis()->SetTitle("#rho_{00}'-input #rho_{00}");
  g_diffI->GetHistogram()->SetTitle("Difference Plot");

  g_diffI->SetMarkerStyle(20);
  g_diffI->SetMarkerColor(kBlack);
  g_diffI->Draw("APE");

  FigureName = Form("./figures/AcceptanceQA/AcceptanceQADifferenceFromIdeal_%s%s_method%d.pdf",vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method);
  c_diffI->SaveAs(FigureName.c_str());
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

  double result = (1.+Bs*D/2.) + (As+D-Bs*D)*CosTheta*CosTheta + (As*D+Bs*D/2.)*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}

double FuncAF_Updated(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double R = par[3];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double result = (2. + F - Bs*F/2.) + (2.*As - F + As*F + Bs*F)*CosTheta*CosTheta - (As*F + Bs*F/2.)*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}

double Func4th(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double F = par[1];
  double G = par[2];

  double result = 1. + (4.*F+3.*G)/8. - (2.*F+3.*G)/4.*CosTheta*CosTheta + 3.*G/8.*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}

/*double Func4th(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double F = par[1];
  double G = par[2];

  double result = 2. + F + 3.*G/4. - F*CosTheta*CosTheta - 3.*G/2.*CosTheta*CosTheta + 3.*G/4.*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}*/

double FuncAFG(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double G = par[3];
  double R = par[4];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double order0 = 2. + F - Bs*F/2. + 3.*G/4. - Bs*G/2.;    
  double order2 = (2.*As - F + As*F + Bs*F - 3.*G/2. + 3.*As*G/4. + 3.*Bs*G/2.)*CosTheta*CosTheta;
  double order4 = (-1.*As*F - Bs*F/2. + 3.*G/4. - 3.*As*G/2. - 3.*Bs*G/2.)*CosTheta*CosTheta*CosTheta*CosTheta;
  double order6 = (3.*As*G/4. + Bs*G/2.)*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta;

  //double result = order0 + order2 + order4 + order6;

  return N*(order0 + order2 + order4 + order6);

}

double FuncAFGH(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double G = par[3];
  double H = par[4];
  double R = par[5];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double order0 = 64. - 16.*(-2. + Bs)*F + 8.*(3. - 2.*Bs)*G + 5.*(4. - 3.*Bs)*H;    
  double order2 = 4.*(As*(16. + 8.*F + 6.*G + 5.*H) + (-1. + Bs)*(8.*F + 12.*G + 15.*H))*CosTheta*CosTheta;
  double order4 = -2.*(8.*(2.*As + Bs)*F - 12.*G + 24.*(As + Bs)*G + 15.*(-2. + 2.*As + 3.*Bs)*H)*CosTheta*CosTheta*CosTheta*CosTheta;
  double order6 = 4.*(6.*As*G + 4.*Bs*G - 5.*H + 15.*(As + Bs)*H)*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta;
  double order8 = -5.*(4.*As + 3.*Bs)*H*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta;

  //double result = order0 + order2 + order4 + order6;

  return N*(order0 + order2 + order4 + order6 + order8);

}
