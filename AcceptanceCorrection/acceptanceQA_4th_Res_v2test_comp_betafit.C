#include "Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"

void acceptanceQA_4th_Res_v2test_comp_betafit(const int energy = 4, const int pid = 0, bool doall = true, std::string res = "0p4", double resval = 0.4, int cut = 1, int method = 1, int v2 = 0) {


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
 
  //std::string ytextfile = "0p2y0p4";
  //std::string ytext = "0.2<y<0.4";
  std::string ytextfile = "yabs1";
  std::string ytext = "-1<y<1";

  Double_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};
  Double_t pt_set[8] = {0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};
//  Double_t pt_set[8] = {3.0, 3.3, 3.6, 3.9, 4.3, 4.6, 4.9, 5.2};
//  Double_t pt_set[8] = {0.6, 1.4, 2.2, 3.0, 3.8, 4.6, 5.4, 7.2};

  TH1D *h_theta_star_before[6][3];
  TH1D *h_theta[6][3];
  TH1D *h_theta_before[6][3];
  TH1D *h_theta_star[6][3];
  TH1D *h_theta_star_RP[6][3];
  TH1D *h_theta_star_before_RP[6][3];
  TH1D *h_out1[6][3];
  TH1D *h_out2[6][3];


  TH1F*  h_cosbeta[6][3];
  TH1F* h_cos2beta[6][3];
  TH1F* h_cos4beta[6][3];
  TH1F*  h_cosbetaP[6][3];
  TH1F* h_cos2betaP[6][3];
  TH1F* h_cos4betaP[6][3];
 
  TProfile*  p_cosbeta[6][3];
  TProfile* p_cos2beta[6][3];
  TProfile* p_cos4beta[6][3];
  TProfile*  p_cosbetaP[6][3];
  TProfile* p_cos2betaP[6][3];
  TProfile* p_cos4betaP[6][3];

  TH1F* h_beta[6][3];
  TH1F* h_tstar[6][3];

  TH1F* h_betaP[6][3];
  TH1F* h_tstarP[6][3];

  TH2F* h_betatstar[6][3];
  TH2F* h_betatstarP[6][3];



  TFile *MCFiles[4]; 
  //int etamode[4] = {3,4,5,0};
  std::string etatext[6] = {"inf","1.0","0.8","0.6","0.4","0.2"};
  double rho00in[3] = {0.2500,0.3333,0.4000};
  int rho00val[3] = {2500,3333,4000};
  TGraphAsymmErrors *g_rho00outBefore[6] ;
  TGraphAsymmErrors *g_rho00outAfter[6]  ;
  TGraphAsymmErrors *g_rho00outBeforeIn[6];
  TGraphAsymmErrors *g_rho00outAfterIn[6] ;
  TGraphAsymmErrors *g_rho00RecoObs[6];

  double Fval[6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Ferr[6][6] = {0.0};
  double Gval[6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Gerr[6][6] = {0.0};
  double FvalB[6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double FerrB[6][6] = {0.0};
  double GvalB[6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double GerrB[6][6] = {0.0};
  double Fpval[5] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Fperr[5] = {0.0};
  double Gpval[5] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Gperr[5] = {0.0};

  double FvalOrig[5] = {0.0};
  double FerrOrig[5] = {0.0};

  double rhoobs[6][6] = {0.0};
  double rhoobsErr[6][6] = {0.0};
  double rhoreal[6][6] = {0.0};
  double rhorealErr[6][6] = {0.0};

  for(int irho = 0; irho < 3; irho++)
  {
    MCFiles[irho] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/ResolutionTesting/ResolutionTest_0p4_%s_v2%d_withBeta_MANYPLOTS/McAcceptanceOutput_pt1_energy%d_pid%d_cent4_EtaMode_0_nrho%d.root",vmsa::mPID[pid].c_str(),ytextfile.c_str(),v2,energy,pid,rho00val[irho]),"READ");

    for(int ieta = 0; ieta < 6; ieta++)
    {
      if(irho == 0){
        g_rho00outBefore[ieta] = new TGraphAsymmErrors();
        g_rho00outAfter[ieta]  = new TGraphAsymmErrors();
        g_rho00outBeforeIn[ieta] = new TGraphAsymmErrors();
        g_rho00outAfterIn[ieta]  = new TGraphAsymmErrors();
        g_rho00RecoObs[ieta]  = new TGraphAsymmErrors();
      }
      std::string histname;
      histname = Form("h_cosbeta_%d",ieta);
      h_cosbeta[ieta][irho]      =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_cosbeta_%d_%d",ieta,irho));
      histname = Form("h_cos2beta_%d",ieta);
      h_cos2beta[ieta][irho]     =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_cos2beta_%d_%d",ieta,irho));
      histname = Form("h_cos4beta_%d",ieta);
      h_cos4beta[ieta][irho]     =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_cos4beta_%d_%d",ieta,irho));
      
      histname = Form("h_cosbetaP_%d",ieta);
      h_cosbetaP[ieta][irho]     =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_cosbetaP_%d_%d",ieta,irho));
      histname = Form("h_cos2betaP_%d",ieta);
      h_cos2betaP[ieta][irho]    =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_cos2betaP_%d_%d",ieta,irho));
      histname = Form("h_cos4betaP_%d",ieta);
      h_cos4betaP[ieta][irho]    =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_cos4betaP_%d_%d",ieta,irho));

      histname = Form("p_cosbeta_%d",ieta);
      p_cosbeta[ieta][irho]      = (TProfile*)((TProfile*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("p_cosbeta_%d_%d",ieta,irho));
      histname = Form("p_cos2beta_%d",ieta);
      p_cos2beta[ieta][irho]     = (TProfile*)((TProfile*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("p_cos2beta_%d_%d",ieta,irho));
      histname = Form("p_cos4beta_%d",ieta);
      p_cos4beta[ieta][irho]     = (TProfile*)((TProfile*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("p_cos4beta_%d_%d",ieta,irho));
      
      histname = Form("p_cosbetaP_%d",ieta);
      p_cosbetaP[ieta][irho]     = (TProfile*)((TProfile*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("p_cosbetaP_%d_%d",ieta,irho));
      histname = Form("p_cos2betaP_%d",ieta);
      p_cos2betaP[ieta][irho]    = (TProfile*)((TProfile*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("p_cos2betaP_%d_%d",ieta,irho));
      histname = Form("p_cos4betaP_%d",ieta);
      p_cos4betaP[ieta][irho]    = (TProfile*)((TProfile*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("p_cos4betaP_%d_%d",ieta,irho));
      
      histname = Form("h_beta_%d",ieta);
      h_beta[ieta][irho]         =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_beta_%d_%d",ieta,irho));
      histname = Form("h_tstar_%d",ieta);
      h_tstar[ieta][irho]        =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_tstar_%d_%d",ieta,irho));
      
      histname = Form("h_betaP_%d",ieta);
      h_betaP[ieta][irho]        =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_betaP_%d_%d",ieta,irho));
      histname = Form("h_tstarP_%d",ieta);
      h_tstarP[ieta][irho]       =         (TH1F*)((TH1F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_tstarP_%d_%d",ieta,irho));
      
      histname = Form("h_betatstar_%d",ieta);
      h_betatstar[ieta][irho]    =         (TH2F*)((TH2F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_betatstar_%d_%d",ieta,irho));
      histname = Form("h_betatstarP_%d",ieta);
      h_betatstarP[ieta][irho]   =         (TH2F*)((TH2F*)MCFiles[irho]->Get(histname.c_str()))->Clone(Form("h_betatstarP_%d_%d",ieta,irho));


      h_theta_star[ieta][irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get(Form("h_theta_star_%d",ieta)))->Clone(Form("h_theta_star_%d_%d",irho,ieta));
      h_theta_star_RP[ieta][irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get(Form("h_theta_star_RP_%d",ieta)))->Clone(Form("h_theta_star_RP_%d_%d",irho,ieta));
      //h_out1[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_out1"))->Clone(Form("h_out1_%d",irho));
      //h_out2[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_out2"))->Clone(Form("h_out2_%d",irho));

      if(ieta > 0)
      { 
        TF1 *Func_rho = new TF1("Func_rho","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
        TF1 *Func_A = new TF1("Func_A",Func4th,0,1,3);
        Func_A->SetLineColor(kBlue);
        TF1 *Func_Orig = new TF1("Func_Orig","[0]*(1 + [1]*x*x)",0,1);
        //Func_Orig->SetLineColor(kRed);
        TF1 *Func_Theta = new TF1("Func_Theta","[0]*(1 + [1]*x*x + [2]*x*x*x*x + [3]*x*x*x*x*x*x)",0,1);
        //TF1 *Func_AFG = new TF1("Func_AFG",FuncAFG,0,1,3);

        TF2 *Func_beta = new TF2("Func_beta","[0] * (1. + [1]*(sin(x)*sin(y))^2 + [2]*(sin(x)*sin(y))^4)",0,TMath::Pi(),0.0,2.0*TMath::Pi());

        //Func_rho->SetParameter(0,h_theta_star[ieta][irho]->GetBinContent(1));
        //Func_rho->SetParameter(1,1./3.);
        //h_theta_star[ieta][irho]->Fit(Func_rho,"NERQ");
        //rhoobs[ieta][irho] = Func_rho->GetParameter(1);
        //rhoobsErr[ieta][irho] = Func_rho->GetParError(1);

        //Func_rho->SetParameter(0,h_theta_star_RP[ieta][irho]->GetBinContent(1));
        //Func_rho->SetParameter(1,1./3.);
        //h_theta_star_RP[ieta][irho]->Fit(Func_rho,"NERQ");
        //rhoreal[ieta][irho] = Func_rho->GetParameter(1);
        //rhorealErr[ieta][irho] = Func_rho->GetParError(1);
        //if(cut == 0) g_rho00outBefore[ieta]->SetPoint(irho,rhoobs[ieta][irho],rhoreal[ieta][irho]);
        //if(cut == 0) g_rho00outBefore[ieta]->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));
        //if(cut == 0) g_rho00outBeforeIn[ieta]->SetPoint(irho,rho00in[irho],rhoreal[ieta][irho]);
        //if(cut == 0) g_rho00outBeforeIn[ieta]->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));


        //TH1D *h_theta_clone = (TH1D*)h_theta[ieta][irho]->Clone("h_theta_clone");
        //h_theta_clone->Sumw2();
        //h_theta_before[ieta][irho]->Sumw2();
        //h_theta_clone->Divide(h_theta_before[ieta][irho]);
        //Func_Theta->SetParameter(0,h_theta_clone->GetBinContent(1));
        //Func_Theta->SetParameter(1,0);
        //h_theta_clone->Fit(Func_Theta,"ER");


        //if(cut == 1 && method == 2){
        //  Fval[ieta][irho] = Func_Theta->GetParameter(1);
        //  Ferr[ieta][irho] = Func_Theta->GetParError(1);
        //  Gval[ieta][irho] = Func_Theta->GetParameter(2);
        //  Gerr[ieta][irho] = Func_Theta->GetParError(2);
        //}
        //TCanvas *c0 = new TCanvas();
        //c0->SetFillColor(0);
        //c0->SetGrid(0,0);
        //c0->SetTitle(0);
        //c0->SetBottomMargin(0.15);
        //c0->SetLeftMargin(0.15);
        //h_theta_clone->Draw("pE");
        //Func_Theta->Draw("same"); 
        //auto legend0 = new TLegend(0.1,0.7,0.48,0.9); 
        //legend0->AddEntry(Func_Theta,"4th Order Fit","l");
        ////legend->AddEntry(Func_Orig,"2nd Order Fit","l");
        //legend0->Draw("same");
        //c0->SaveAs(Form("FValues/yabs1%s/pt1_%d%s_method%d_FROMTHETA.pdf",cutoption.c_str(),irho,cutoption.c_str(),method));
        //delete c0;

        TH2D *h_betatstarP_ratio = (TH2D*)h_betatstarP[ieta][irho]->Clone("h_betatstarP_ratio");
        h_betatstarP_ratio->Sumw2();
        h_betatstarP[0][irho]->Sumw2();
        h_betatstarP_ratio->Divide(h_betatstarP[0][irho]); 
        Func_beta->SetParameter(0,h_betatstarP_ratio->GetBinContent(1));
        Func_beta->SetParameter(1,0.0);
        Func_beta->SetParameter(2,0.0);
        //Func_A->FixParameter(2,0.0);
        //Func_beta->SetParameter(0,h_theta_star_clone->GetBinContent(1));
        //Func_beta->SetParameter(1,0);
        //h_theta_star_clone->Fit(Func_Orig,"ER");

        TH1D *h_theta_star_clone = (TH1D*)h_theta_star[ieta][irho]->Clone("h_theta_star_clone");
        h_theta_star_clone->Sumw2();
        h_theta_star[0][irho]->Sumw2();
        h_theta_star_clone->Divide(h_theta_star[0][irho]); 
        Func_A->SetParameter(0,h_theta_star_clone->GetBinContent(1));
        Func_A->SetParameter(1,0);
        Func_A->SetParameter(2,0.0);
        //Func_A->FixParameter(2,0.0);
        Func_Orig->SetParameter(0,h_theta_star_clone->GetBinContent(1));
        Func_Orig->SetParameter(1,0);
        h_theta_star_clone->Fit(Func_Orig,"ER");

        if(cut == 1 && method == 0){
          h_betatstarP_ratio->Fit(Func_beta,"ER");
          FvalB[ieta][irho] = Func_beta->GetParameter(1);
          FerrB[ieta][irho] = Func_beta->GetParError(1);
          GvalB[ieta][irho] = Func_beta->GetParameter(2);
          GerrB[ieta][irho] = Func_beta->GetParError(2);
          TCanvas *c1 = new TCanvas();
          c1->SetFillColor(0);
          c1->SetGrid(0,0);
          c1->SetTitle(0);
          c1->SetBottomMargin(0.15);
          c1->SetLeftMargin(0.15);
          h_betatstarP_ratio->Draw("colz");
          Func_beta->Draw("same"); 
          //legend->AddEntry(Func_A,"4th Order Fit","l");
          //legend->AddEntry(Func_Orig,"2nd Order Fit","l");
          c1->SaveAs(Form("FValues/yabs1/pt1_rho%d_eta%d%s_method%d_2D.pdf",irho,ieta,cutoption.c_str(),method));

          h_theta_star_clone->Fit(Func_A,"ER");
          Fval[ieta][irho] = Func_A->GetParameter(1);
          Ferr[ieta][irho] = Func_A->GetParError(1);
          Gval[ieta][irho] = Func_A->GetParameter(2);
          Gerr[ieta][irho] = Func_A->GetParError(2);
          c1->SetFillColor(0);
          c1->SetGrid(0,0);
          c1->SetTitle(0);
          c1->SetBottomMargin(0.15);
          c1->SetLeftMargin(0.15);
          h_theta_star_clone->Draw("pE");
          Func_A->Draw("same"); 
          auto legend = new TLegend(0.1,0.7,0.48,0.9); 
          //legend->AddEntry(Func_A,"4th Order Fit","l");
          //legend->AddEntry(Func_Orig,"2nd Order Fit","l");
          legend->Draw("same");
          c1->SaveAs(Form("FValues/yabs1/pt1_%d%s_method%d.pdf",irho,cutoption.c_str(),method));
          delete c1;
        }
 
        FvalOrig[irho] = Func_Orig->GetParameter(1);
        FerrOrig[irho] = Func_Orig->GetParError(1);

        TH2D *h_betatstar_ratio = (TH2D*)h_betatstar[ieta][irho]->Clone("h_betatstar_ratio");
        h_betatstar_ratio->Sumw2();
        h_betatstar[0][irho]->Sumw2();
        h_betatstar_ratio->Divide(h_betatstar[0][irho]); 
        Func_beta->SetParameter(0,h_betatstar_ratio->GetBinContent(1));
        Func_beta->SetParameter(1,0.0);
        Func_beta->SetParameter(2,0.0);

        TH1D *h_theta_star_clone_RP = (TH1D*)h_theta_star_RP[ieta][irho]->Clone("h_theta_star_clone_RP");
        h_theta_star_clone_RP->Sumw2();
        h_theta_star_RP[0][irho]->Sumw2();
        h_theta_star_clone_RP->Divide(h_theta_star_RP[0][irho]); 
        Func_A->SetParameter(0,h_theta_star_clone_RP->GetBinContent(1));
        Func_A->SetParameter(1,0);
        Func_A->SetParameter(2,0.0);
        //Func_A->FixParameter(2,0.0);
        Func_Orig->SetParameter(0,h_theta_star_clone_RP->GetBinContent(1));
        Func_Orig->SetParameter(1,0);
        h_theta_star_clone_RP->Fit(Func_Orig,"ER");
        h_theta_star_clone_RP->Fit(Func_A,"ER");
        h_betatstar_ratio->Fit(Func_beta,"ER");

        if(cut == 1 && method == 1){
          FvalB[ieta][irho] = Func_beta->GetParameter(1);
          FerrB[ieta][irho] = Func_beta->GetParError(1);
          GvalB[ieta][irho] = Func_beta->GetParameter(2);
          GerrB[ieta][irho] = Func_beta->GetParError(2);
          TCanvas *c1 = new TCanvas();
          c1->SetFillColor(0);
          c1->SetGrid(0,0);
          c1->SetTitle(0);
          c1->SetBottomMargin(0.15);
          c1->SetLeftMargin(0.15);
          h_betatstar_ratio->Draw("colz");
          Func_beta->Draw("same"); 
          //legend->AddEntry(Func_A,"4th Order Fit","l");
          //legend->AddEntry(Func_Orig,"2nd Order Fit","l");
          c1->SaveAs(Form("FValues/yabs1/pt1_rho%d_eta%d%s_method%d_2D.pdf",irho,ieta,cutoption.c_str(),method));

          Fval[ieta][irho] = Func_A->GetParameter(1);
          Ferr[ieta][irho] = Func_A->GetParError(1);
          Gval[ieta][irho] = Func_A->GetParameter(2);
          Gerr[ieta][irho] = Func_A->GetParError(2);
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
          c1->SaveAs(Form("FValues/yabs/pt1_%d%s_method%d.pdf",irho,cutoption.c_str(),method));
          delete c1;
        }

        TF1 *Func_AFG = new TF1("Func_AFG",FuncAFG,0,1,5);
        //TF1 *Func_AFG = new TF1("Func_AFG",FuncAFGH,0,1,6);
        Func_AFG->SetParameter(0,h_theta_star_RP[ieta][irho]->GetBinContent(1));
        Func_AFG->SetParameter(1,1./3.);
        Func_AFG->FixParameter(2, Func_A->GetParameter(1));//Fval[irho]);
        Func_AFG->FixParameter(3, Func_A->GetParameter(2));//Gval[irho]);
        Func_AFG->FixParameter(4, 1.0);//Res_12);
        //if(cut == 1 && method == 2) Func_AFG->FixParameter(2, Fval[irho]);
        //if(cut == 1 && method == 2) Func_AFG->FixParameter(3, Gval[irho]);
        //if(cut == 1 && method == 2) Func_AFG->FixParameter(4, Hval[irho]);
        //if(cut == 1 && method == 2)  Func_AFG->FixParameter(5, 1.0);//Res_12);
        h_theta_star_RP[ieta][irho]->Fit(Func_AFG,"NERQ");
        rhoreal[ieta][irho] = Func_AFG->GetParameter(1);
        rhorealErr[ieta][irho] = Func_AFG->GetParError(1);
        cout << "rho00in = " << rho00in[irho] << ",    rho00 from RP with acceptance correction = " << rhoreal[ieta][irho] << endl;
        if(cut == 1) g_rho00outBefore[ieta]->SetPoint(irho,rhoobs[ieta][irho],rhoreal[ieta][irho]);
        if(cut == 1) g_rho00outBefore[ieta]->SetPointError(irho,0.0,0.0,Func_AFG->GetParError(1),Func_AFG->GetParError(1));
        if(cut == 1) g_rho00outBeforeIn[ieta]->SetPoint(irho,rho00in[irho],rhoreal[ieta][irho]);
        if(cut == 1) g_rho00outBeforeIn[ieta]->SetPointError(irho,0.0,0.0,Func_AFG->GetParError(1),Func_AFG->GetParError(1));


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
        //delete h_theta_clone;
        delete h_theta_star_clone;
        delete Func_rho;
        delete Func_A;
        delete Func_Theta;
        //MCFiles[irho]->Close();
      }
    }
  }
 
  double FfromFp[5] = {0.0};
  for(int ieta = 1; ieta < 6; ieta++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
    //FfromFp[irho] = -Fval[irho]/(2.+Fval[irho]);
    //double err = TMath::Abs((Fval[irho]/(2.-Fval[irho])/(2.-Fval[irho])-1./(2.-Fval[irho]))*Ferr[irho]);
    //cout << "F from F* = " << FfromFp[irho] << " +/- " << err << "      F from Orig Method = " << FvalOrig[irho] << " +/- " << FerrOrig[irho] << endl;
    //cout << "F from theta = " << Fpval[irho] << " +/- " << Fperr[ierr] << "      G from theta = " << Gpval[irho] << " +/- " << Gperr[irho] << endl;
      cout << "irho = " << irho << ",   ieta = " << ieta << endl;
      cout << "F = " << Fval[ieta][irho] << " +/- " << Ferr[ieta][irho] << "           F from 2D = " << FvalB[ieta][irho] << " +/- " << FerrB[ieta][irho] << endl;
      cout << "G = " << Gval[ieta][irho] << " +/- " << Gerr[ieta][irho] << "           G from 2D = " << GvalB[ieta][irho] << " +/- " << GerrB[ieta][irho] << endl;
    } 
  }


  //TCanvas *c1[5];
  TGraphAsymmErrors *g_diff[6]; 
  TGraphAsymmErrors *g_diffI[6];

  for(int ieta = 1; ieta < 6; ieta++)
  {
    g_diff[ieta]  = new TGraphAsymmErrors();
    g_diffI[ieta] = new TGraphAsymmErrors();
    for(int irho = 0; irho < 3; irho++) 
    {
      cout << "ETA = " << ieta << ", IRHO = " << irho << endl;
      //TF1 *Func_rho = new TF1("Func_rho",FuncAFG,0,1,5);
      TF1 *Func_rho = new TF1("Func_rho",FuncAFG,0,1,5);
      Func_rho->SetLineColor(kBlue);
      TF1 *Func_rhoOrig = new TF1("Func_rhoOrig",FuncAD,0,1,4);
      Func_rhoOrig->SetLineColor(kRed);
      //TF1 *Func_rho = new TF1("Func_rho",FuncAF_Updated,0,1,4);
      //TF1 *Func_rdl = new TF1("Func_rdl","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);

      h_theta_star[ieta][irho]->GetXaxis()->SetTitleOffset(1.2);
      h_theta_star[ieta][irho]->GetXaxis()->SetTitle("cos#theta*");
      h_theta_star[ieta][irho]->GetYaxis()->SetTitle("yield");
      h_theta_star[ieta][irho]->GetYaxis()->SetTitleOffset(1.0);
      h_theta_star[ieta][irho]->SetMarkerColor(2);
      h_theta_star[ieta][irho]->SetMarkerSize(1.8);
      h_theta_star[ieta][irho]->SetMarkerStyle(21);
      //h_theta_star[irho]->Draw();
      Func_rho->SetParameter(0, h_theta_star[ieta][irho]->GetBinContent(5));
      Func_rho->SetParameter(1, 0.333);
      //g_rho00outAfter->SetPoint(irho,rho00in[irho],Func_rho->GetParameter(1));
      //g_rho00outAfter->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));
      //Func_rho->FixParameter(2, -1.57811);//Fval[1]);
      //Func_rho->FixParameter(3, 0.649864);//Gval[1]);
      Func_rho->FixParameter(2, Fval[ieta][1]);
      Func_rho->FixParameter(3, Gval[ieta][1]);
      //Func_rho->FixParameter(2, FvalB[ieta][1]);
      //Func_rho->FixParameter(3, GvalB[ieta][1]);
      //Func_rho->FixParameter(4, Hval[irho]);
      Func_rho->FixParameter(4, resval);//Res_12);

      //Func_rho->FixParameter(3, 1.0);//Res_12);

      Func_rhoOrig->SetParameter(0, h_theta_star[ieta][irho]->GetBinContent(5));
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
      h_theta_star[ieta][irho]->Draw("pE");
      h_theta_star[ieta][irho]->Fit(Func_rhoOrig, "NMI");
      //Func_rhoOrig->Draw("same");
      h_theta_star[ieta][irho]->Fit(Func_rho, "NMI");
      Func_rho->Draw("same");
 
      auto legend = new TLegend(0.1,0.7,0.48,0.9); 
      legend->AddEntry(Func_rho,"4th Order Fit","l");
      //legend->AddEntry(Func_rhoOrig,"2nd Order Fit","l");
      legend->Draw("same");
      
      c1->SaveAs(Title->Data());
      //double rho00final = 4./(1.+3.*resval)*(Func_rho->GetParameter(1)-1/3.)+1./3.;
      double rho00final = Func_rho->GetParameter(1);
      cout << "rho after = " << rho00final << " +/- " << Func_rho->GetParError(1) << endl;
      g_rho00outAfter[ieta]->SetPoint(irho,rhoobs[ieta][irho],rho00final);
      g_rho00outAfter[ieta]->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));
      g_rho00outAfterIn[ieta]->SetPoint(irho,rho00in[irho],rho00final);
      g_rho00outAfterIn[ieta]->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));
      g_rho00RecoObs[ieta]->SetPoint(irho,rho00in[irho],rhoobs[ieta][irho]);
      g_rho00RecoObs[ieta]->SetPointError(irho,0.0,0.0,rhoobsErr[ieta][irho],rhoobsErr[ieta][irho]);
      cout << "rho input = " << rho00in[irho] << ",    rho observed = " << rhoobs[ieta][irho] << endl;

      cout << "rho00in = " << rho00in[irho] << ",    rho00 from RP with acceptance correction = " << rhoreal[ieta][irho] << " +/- " << rhorealErr[ieta][irho] << endl;
      g_diff[ieta]->SetPoint(irho,rho00in[irho],rho00final-rhoreal[ieta][irho]);
      g_diffI[ieta]->SetPoint(irho,rho00in[irho]+0.0015*ieta,rho00final-rho00in[irho]);
      double uncorrErr = TMath::Sqrt(Func_rho->GetParError(1)*Func_rho->GetParError(1)+rhorealErr[ieta][irho]*rhorealErr[ieta][irho]);
      double corrErr = TMath::Sqrt(fabs(Func_rho->GetParError(1)*Func_rho->GetParError(1)-rhorealErr[ieta][irho]*rhorealErr[ieta][irho]));
      std::cout << "uncorrelated Error = " << uncorrErr << ",      correlated Error = " << corrErr << std::endl;
      g_diff[ieta]->SetPointError(irho,0.0,0.0,corrErr,corrErr);
      g_diffI[ieta]->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));

      delete c1;  
      delete Func_rho;
      delete Func_rhoOrig;
      delete Title;

    }
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
  h_play->GetXaxis()->SetLimits(0.22,0.42);
  h_play->GetYaxis()->SetRangeUser(0.18,0.58);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");
  //h_play->Draw("pE");

  //TF1 *linear = new TF1("linear","4./(1+3*[0])*(x-1./3.) + 1./3.",0.24,0.42);
  //linear->FixParameter(0,resval);
  //linear->Draw("same");

  //g_rho00outBefore->SetMarkerStyle(24);
  //g_rho00outBefore->SetMarkerColor(kRed);
  //g_rho00outBefore->SetMarkerSize(3.0);
  //g_rho00outBefore->SetLineWidth(3);
  //g_rho00outBefore->Draw("pE Same");

  //g_rho00outAfter->SetMarkerStyle(26);
  //g_rho00outAfter->SetMarkerColor(kBlue);
  //g_rho00outAfter->SetMarkerSize(3.0);
  //g_rho00outAfter->SetLineWidth(3);
  //g_rho00outAfter->Draw("pE Same");

  //TLegend *leg1 = new TLegend(0.2,0.70,0.5,0.85);
  //leg1->SetFillColor(10);
  //leg1->SetBorderSize(0);
  //leg1->AddEntry(g_rho00outBefore,"Fit to cos(#theta^{*})","p");
  //leg1->AddEntry(g_rho00outAfter,"Fit to cos(#theta^{*'}) using R","p");
  ////leg1->AddEntry(linear,"#frac{4}{1+3R} (#rho_{00}^{obs} - #frac{1}{3}) + #frac{1}{3}","l");
  //leg1->Draw("same");

  //string FigureName = Form("./figures/AcceptanceQA/AcceptanceQA_%s%s_method%d.pdf",vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method);
  //c_play->SaveAs(FigureName.c_str());

  //TCanvas *c_play2 = new TCanvas("c_play2","c_play2",10,10,800,800);
  //c_play2->SetLeftMargin(0.15);
  //c_play2->SetBottomMargin(0.15);
  //c_play2->SetGrid(0,0);
  //c_play2->SetTicks(1,1);
  //c_play2->cd();

  //h_play->GetXaxis()->SetTitle("#rho_{00} input");
  //h_play->Draw("pE");

  //TF1 *linear = new TF1("linear","x",0.24,0.42);
  ////linear->FixParameter(0,resval);
  //linear->Draw("same");


  //g_rho00outBeforeIn->SetMarkerStyle(24);
  //g_rho00outBeforeIn->SetMarkerColor(kRed);
  //g_rho00outBeforeIn->SetMarkerSize(3.0);
  //g_rho00outBeforeIn->SetLineWidth(3);
  //g_rho00outBeforeIn->Draw("pE Same");

  //g_rho00outAfterIn->SetMarkerStyle(26);
  //g_rho00outAfterIn->SetMarkerColor(kBlue);
  //g_rho00outAfterIn->SetMarkerSize(3.0);
  //g_rho00outAfterIn->SetLineWidth(3);
  //g_rho00outAfterIn->Draw("pE Same");

  //TLegend *leg2 = new TLegend(0.2,0.70,0.5,0.85);
  //leg2->SetFillColor(10);
  //leg2->SetBorderSize(0);
  //leg2->AddEntry(g_rho00outBeforeIn,"Fit to cos(#theta^{*})","p");
  //leg2->AddEntry(g_rho00outAfterIn,"Fit to cos(#theta^{*'}) using R","p");
  ////leg1->AddEntry(linear,"#frac{4}{1+3R} (#rho_{00}^{obs} - #frac{1}{3}) + #frac{1}{3}","l");
  //leg2->Draw("same");

  //FigureName = Form("./figures/AcceptanceQA/AcceptanceQAInput_%s%s_method%d.pdf",vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method);
  //c_play2->SaveAs(FigureName.c_str());

  TCanvas *c_play3 = new TCanvas("c_play3","c_play3",10,10,800,800);
  c_play3->SetLeftMargin(0.15);
  c_play3->SetBottomMargin(0.15);
  c_play3->SetGrid(0,0);
  c_play3->SetTicks(1,1);
  c_play3->cd();

  h_play->DrawCopy("pE");

  TF1 *linear = new TF1("linear","x",0.22,0.42);
  //linear->FixParameter(0,resval);
  linear->Draw("same");

  h_play->GetYaxis()->SetTitle("#rho_{00}^{reco}");
  h_play->GetXaxis()->SetTitle("#rho_{00}^{input}");
  TLegend *leg2 = new TLegend(0.6,0.175,0.9,0.425);
  leg2->SetFillColor(10);
  leg2->SetBorderSize(0);

  int styleb[6] = {0,29,21,34,20,33};
  int stylea[6] = {0,30,25,28,24,27};
  int color[6] = {0,3,4,6,2,9};

  for(int ieta = 1; ieta < 6; ieta++)
  {
    g_rho00RecoObs[ieta]->SetMarkerStyle(styleb[ieta]);
    g_rho00RecoObs[ieta]->SetMarkerColor(color[ieta]);
    g_rho00RecoObs[ieta]->SetMarkerSize(2.0);
    g_rho00RecoObs[ieta]->SetLineWidth(2);
    g_rho00RecoObs[ieta]->Draw("pE Same");

    g_rho00outAfterIn[ieta]->SetMarkerStyle(stylea[ieta]);
    g_rho00outAfterIn[ieta]->SetMarkerColor(color[ieta]);
    g_rho00outAfterIn[ieta]->SetMarkerSize(2.0);
    g_rho00outAfterIn[ieta]->SetLineWidth(2);
    g_rho00outAfterIn[ieta]->Draw("pE Same");

    leg2->AddEntry(g_rho00RecoObs[ieta],Form("|#eta|<%s, %s, Before Correction",etatext[ieta].c_str(),ytext.c_str()),"p");
    leg2->AddEntry(g_rho00outAfterIn[ieta],Form("|#eta|<%s, %s, After Correction",etatext[ieta].c_str(),ytext.c_str()),"p");
  }
  leg2->AddEntry(linear,"reco = input","l");
  leg2->Draw("same");

  std::string FigureName = Form("./figures/AcceptanceQA/AcceptanceQAInput_2D_v2%d_%s_RecoObs_%s%s_method%d.pdf",v2,ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method);
  c_play3->SaveAs(FigureName.c_str());


  //TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,800,800);
  //c_diff->SetLeftMargin(0.15);
  //c_diff->SetBottomMargin(0.15);
  //c_diff->SetGrid(0,0);
  //c_diff->SetTicks(1,1);
  //c_diff->cd();

  //g_diff->GetHistogram()->GetXaxis()->SetTitle("Input #rho_{00}");
  //g_diff->GetHistogram()->GetYaxis()->SetTitle("#rho_{00}'-#rho_{00}");
  //g_diff->GetHistogram()->SetTitle("Difference Plot");

  //g_diff->SetMarkerStyle(20);
  //g_diff->SetMarkerColor(kBlack);
  //g_diff->Draw("APE");

  //FigureName = Form("./figures/AcceptanceQA/AcceptanceQADifference_%s%s_method%d.pdf",vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method);
  //c_diff->SaveAs(FigureName.c_str());

  TCanvas *c_diffI = new TCanvas("c_diffI","c_diffI",10,10,800,800);
  c_diffI->SetLeftMargin(0.15);
  c_diffI->SetBottomMargin(0.15);
  c_diffI->SetGrid(0,0);
  c_diffI->SetTicks(1,1);
  c_diffI->cd();

  h_play->GetXaxis()->SetLimits(0.22,0.42);
  h_play->GetYaxis()->SetRangeUser(0.18,0.58);
  h_play->GetYaxis()->SetTitle("#rho_{00}^{reco}-#rho_{00}^{input}");

  TLegend *leg3 = new TLegend(0.65,0.65,0.9,0.85);
  leg3->SetFillColor(10);
  leg3->SetBorderSize(0);

  double max = -999;
  double min =  999;
  for(int ieta = 1; ieta < 6; ieta++)
  {
    g_diffI[ieta]->SetMarkerStyle(styleb[ieta]);
    g_diffI[ieta]->SetMarkerColor(color[ieta]);
    g_diffI[ieta]->SetLineColor(color[ieta]);
    g_diffI[ieta]->SetMarkerSize(2.0);
    g_diffI[ieta]->SetLineWidth(2);
  
    //double tmin = g_diffI[ieta]->GetMinimum();
    //double tmax = g_diffI[ieta]->GetMaximum();
    //if(tmin < min) min = tmin;
    //if(tmax > max) max = tmax;
   
    //cout << "min = " << tmin << endl;
    //cout << "max = " << tmax << endl;

    leg3->AddEntry(g_diffI[ieta],Form("|#eta|<%s, %s",etatext[ieta].c_str(),ytext.c_str()),"p");
  }
  //leg2->AddEntry(linear,"","l");
 
  h_play->GetYaxis()->SetRangeUser(-0.01,0.01);
  h_play->DrawCopy("pE");
  for(int ieta = 1; ieta < 6; ieta++)
  {  
    g_diffI[ieta]->Draw("pE Same");
  }
  leg3->Draw("same");
  TF1 *zero = new TF1("zero","0",0.22,0.42);
  //linear->FixParameter(0,resval);
  zero->SetLineStyle(2);
  zero->Draw("same");

  FigureName = Form("./figures/AcceptanceQA/AcceptanceQADifferenceFromIdeal_2D_v2%d_%s_%s%s_method%d.pdf",v2,ytextfile.c_str(),vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method);
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
