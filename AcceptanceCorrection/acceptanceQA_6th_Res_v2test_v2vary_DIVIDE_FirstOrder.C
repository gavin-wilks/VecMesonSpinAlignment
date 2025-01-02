#include "Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"

void acceptanceQA_6th_Res_v2test_v2vary_DIVIDE_FirstOrder(const int energy = 4, const int pid = 0, std::string res = "0p4", double resval = 0.250968, int cut = 1, int method = 0, int v2 = 1, int fit2d = 0, int EP = 1, int etamode = 0) {

  //int ieta = 1;
  //int etamode = 10;

  std::string cutoption = "";
  if(cut == 1) cutoption = "_EtaCut";

  gROOT->Reset();
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendFont(22);
  gStyle->SetLabelFont(22);
  gStyle->SetTitleFont(22);
 
  //std::string ytextfile = "0p2y0p4";
  //std::string ytext = "0.2<y<0.4";
  //std::string ytextfile = "yabs1";
  //std::string ytext = "0.8<y<1.0";
  //std::string ytextfile = "yabs1";
  std::string ytext = "-1.0<y<1.0";
  std::string ytextfile = "yabs1";
  //std::string ytext = "-1<y<1";

  std::string v2text[11] = {"0.0","0.025","0.050","0.075","0.10","0.125","0.15","0.175","0.2","0.225","0.25"};
  double      v2val[11] = {0.0,0.025,0.050,0.075,0.10,0.125,0.15,0.175,0.2,0.225,0.25};

  Double_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};
  Double_t pt_set[8] = {0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};
//  Double_t pt_set[8] = {3.0, 3.3, 3.6, 3.9, 4.3, 4.6, 4.9, 5.2};
//  Double_t pt_set[8] = {0.6, 1.4, 2.2, 3.0, 3.8, 4.6, 5.4, 7.2};

  TH1D *h_theta_star_before[11][6][3];
  TH1D *h_theta[11][6][3];
  TH1D *h_theta_before[11][6][3];
  TH1D *h_theta_star[11][6][3];
  TH1D *h_theta_star_RP[11][6][3];
  TH1D *h_theta_star_before_RP[11][6][3];

  //TH1D *h_theta_star_ratio[11][6][3];
  //TH1D *h_theta_star_ratio_RP[11][6][3];


  TH1F*  h_cosbeta[11][6][3];
  TH1F* h_cos2beta[11][6][3];
  TH1F* h_cos4beta[11][6][3];
  TH1F*  h_cosbetaP[11][6][3];
  TH1F* h_cos2betaP[11][6][3];
  TH1F* h_cos4betaP[11][6][3];
 
  TProfile*  p_cosbeta[11][6][3];
  TProfile* p_cos2beta[11][6][3];
  TProfile* p_cos4beta[11][6][3];
  TProfile*  p_cosbetaP[11][6][3];
  TProfile* p_cos2betaP[11][6][3];
  TProfile* p_cos4betaP[11][6][3];

  TH1F* h_beta[11][6][3];
  TH1F* h_tstar[11][6][3];

  TH1F* h_betaP[11][6][3];
  TH1F* h_tstarP[11][6][3];

  TH2F* h_betatstar[11][6][3];
  TH2F* h_betatstarP[11][6][3];



  TFile *MCFiles[11][6]; 
  std::string etatext[6] = {"inf","1.0","0.8","0.6","0.4","0.2"};
  double rho00in[3] = {0.2500,0.3333,0.4000};
  int rho00val[3] = {2500,3333,4000};
  
  TGraphAsymmErrors *g_rho00outBefore[6] ;
  TGraphAsymmErrors *g_rho00outAfter[6] ;
  TGraphAsymmErrors *g_rho00outBeforeIn[6];
  TGraphAsymmErrors *g_rho00outAfterIn[6];
  TGraphAsymmErrors *g_rho00RecoObs[6];

  double Fval[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Ferr[11][6][6] = {0.0};
  double Gval[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Gerr[11][6][6] = {0.0};
  double Hval[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Herr[11][6][6] = {0.0};
  double FvalRP[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double FerrRP[11][6][6] = {0.0};
  double GvalRP[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double GerrRP[11][6][6] = {0.0};
  double HvalRP[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double HerrRP[11][6][6] = {0.0};
  double FvalB[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double FerrB[11][6][6] = {0.0};
  double GvalB[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double GerrB[11][6][6] = {0.0};
  double FvalBRP[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double FerrBRP[11][6][6] = {0.0};
  double GvalBRP[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double GerrBRP[11][6][6] = {0.0};
  double Fpval[5] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Fperr[5] = {0.0};
  double Gpval[5] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Gperr[5] = {0.0};

  double FvalOrig[5] = {0.0};
  double FerrOrig[5] = {0.0};

  double rhoobs[11][6][6] = {0.0};
  double rhoobsErr[11][6][6] = {0.0};
  double rhoreal[11][6][6] = {0.0};
  double rhorealErr[11][6][6] = {0.0};

  double rhoBefore[11][4] = {0.0};
  double rhoBeforeErr[11][4] = {0.0};

  TGraphAsymmErrors *g_rhoResCorr[4];
  TGraphAsymmErrors *g_rhoResCorrDiff[4];
  
 int  rhoorder[3]={0,1,2};

  for(int ir = 0; ir < 3; ir++)
  {
    int irho = rhoorder[ir];
    g_rho00outBefore[irho] = new TGraphAsymmErrors();
    g_rho00outAfter[irho]  = new TGraphAsymmErrors();
    g_rho00outBeforeIn[irho] = new TGraphAsymmErrors();
    g_rho00outAfterIn[irho]  = new TGraphAsymmErrors();
    g_rho00RecoObs[irho]  = new TGraphAsymmErrors();
    g_rhoResCorr[irho] = new TGraphAsymmErrors();
    g_rhoResCorrDiff[irho] = new TGraphAsymmErrors();

    for(int iem = 0; iem < 1; iem++)
    {

      cout << "iem = " << iem << ",     irho = " << irho << endl;
      MCFiles[iem][irho] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/FirstOrder/McAcceptanceOutput_pt1_energy%d_pid%d_cent5_EtaMode_0_nrho%d.root",vmsa::mPID[pid].c_str(),energy,pid,rho00val[irho]),"READ");
     
      cout << "created graphs" << endl;

      for(int i = 0; i < 6; i++)
      {
        cout << "iem = " << iem << ",    ieta = " << i << ",     irho = " << irho << endl;
        std::string histname;
        histname = Form("h_cosbeta_%d",i);
        h_cosbeta[iem][i][irho]      =         (TH1F*)((TH1F*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("h_cosbeta_%d_%d_%d",iem,i,irho));
        histname = Form("h_cos2beta_%d",i);
        h_cos2beta[iem][i][irho]     =         (TH1F*)((TH1F*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("h_cos2beta_%d_%d_%d",iem,i,irho));
        histname = Form("h_cos4beta_%d",i);
        h_cos4beta[iem][i][irho]     =         (TH1F*)((TH1F*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("h_cos4beta_%d_%d_%d",iem,i,irho));
        
        histname = Form("h_cosbetaP_%d",i);
        h_cosbetaP[iem][i][irho]     =         (TH1F*)((TH1F*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("h_cosbetaP_%d_%d_%d",iem,i,irho));
        histname = Form("h_cos2betaP_%d",i);
        h_cos2betaP[iem][i][irho]    =         (TH1F*)((TH1F*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("h_cos2betaP_%d_%d_%d",iem,i,irho));
        histname = Form("h_cos4betaP_%d",i);
        h_cos4betaP[iem][i][irho]    =         (TH1F*)((TH1F*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("h_cos4betaP_%d_%d_%d",iem,i,irho));

        histname = Form("p_cosbeta_%d",i);
        p_cosbeta[iem][i][irho]      = (TProfile*)((TProfile*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("p_cosbeta_%d_%d_%d",iem,i,irho));
        histname = Form("p_cos2beta_%d",i);
        p_cos2beta[iem][i][irho]     = (TProfile*)((TProfile*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("p_cos2beta_%d_%d_%d",iem,i,irho));
        histname = Form("p_cos4beta_%d",i);
        p_cos4beta[iem][i][irho]     = (TProfile*)((TProfile*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("p_cos4beta_%d_%d_%d",iem,i,irho));
        
        histname = Form("p_cosbetaP_%d",i);
        p_cosbetaP[iem][i][irho]     = (TProfile*)((TProfile*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("p_cosbetaP_%d_%d_%d",iem,i,irho));
        histname = Form("p_cos2betaP_%d",i);
        p_cos2betaP[iem][i][irho]    = (TProfile*)((TProfile*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("p_cos2betaP_%d_%d_%d",iem,i,irho));
        histname = Form("p_cos4betaP_%d",i);
        p_cos4betaP[iem][i][irho]    = (TProfile*)((TProfile*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("p_cos4betaP_%d_%d_%d",iem,i,irho));
        
        histname = Form("h_beta_%d",i);
        h_beta[iem][i][irho]         =         (TH1F*)((TH1F*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("h_beta_%d_%d_%d",iem,i,irho));
        histname = Form("h_tstar_%d",i);
        h_tstar[iem][i][irho]        =         (TH1F*)((TH1F*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("h_tstar_%d_%d_%d",iem,i,irho));
        
        histname = Form("h_betaP_%d",i);
        h_betaP[iem][i][irho]        =         (TH1F*)((TH1F*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("h_betaP_%d_%d_%d",iem,i,irho));
        histname = Form("h_tstarP_%d",i);
        h_tstarP[iem][i][irho]       =         (TH1F*)((TH1F*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("h_tstarP_%d_%d_%d",iem,i,irho));
        
        histname = Form("h_betatstar_%d",i);
        h_betatstar[iem][i][irho]    =         (TH2F*)((TH2F*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("h_betatstar_%d_%d_%d",iem,i,irho));
        histname = Form("h_betatstarP_%d",i);
        h_betatstarP[iem][i][irho]   =         (TH2F*)((TH2F*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("h_betatstarP_%d_%d_%d",iem,i,irho));


        h_theta_star[iem][i][irho] = (TH1D*)((TH1D*)MCFiles[iem][irho]->Get(Form("h_theta_star_%d",i)))->Clone(Form("h_theta_star_%d_%d_%d",iem,i,irho));
        h_theta_star_RP[iem][i][irho] = (TH1D*)((TH1D*)MCFiles[iem][irho]->Get(Form("h_theta_star_RP_%d",i)))->Clone(Form("h_theta_star_RP_%d_%d_%d",iem,i,irho));

        h_theta_star[iem][i][irho]->Print();
        h_theta_star_RP[iem][i][irho]->Print();
     }
       
      cout << "Loaded all histograms" << endl;

    }
  }

  TGraphAsymmErrors *g_diff[6]; 
  TGraphAsymmErrors *g_diffI[6];

  //double ridx[3] = {1,0,2};
  for(int iem = 0; iem < 1/*11*/; iem++)
  {
    for(int ieta = 1; ieta < 4; ieta++)
    {
      g_rho00outAfterIn[ieta] = new TGraphAsymmErrors();
      g_rho00outBeforeIn[ieta] = new TGraphAsymmErrors();
      for(int irho = 0; irho < 3; irho++) 
      {
        //irho = ridx[ir];
        //TF1 *Func_rho = new TF1("Func_rho",FuncAFG,0,1,5);
        TF1 *Func_rho = new TF1("Func_rho",FuncAFGH,0,1,6);
        //TF2 *Func_2D = new TF2("Func_2D",Func2D_wB,0.0,TMath::Pi(),0.0,2.0*TMath::Pi(),6);
        TF2 *Func_2D = new TF2("Func_2D",Func2D,0.0,TMath::Pi(),0.0,2.0*TMath::Pi(),5);
        Func_rho->SetLineColor(kBlue);
        TF1 *Func_rhoOrig = new TF1("Func_rhoOrig",FuncAD,0,1,4);
        Func_rhoOrig->SetLineColor(kRed);
        //TF1 *Func_rho = new TF1("Func_rho",FuncAF_Updated,0,1,4);
        TF1 *Func_rdl = new TF1("Func_rdl","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
        TF1 *Func_rdlB = new TF1("Func_rdlB","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);

        h_theta_star[iem][ieta][irho]->GetXaxis()->SetTitleOffset(1.2);
        h_theta_star[iem][ieta][irho]->GetXaxis()->SetTitle("cos#theta*");
        h_theta_star[iem][ieta][irho]->GetYaxis()->SetTitle("yield");
        h_theta_star[iem][ieta][irho]->GetYaxis()->SetTitleOffset(1.0);
        h_theta_star[iem][ieta][irho]->SetMarkerColor(2);
        h_theta_star[iem][ieta][irho]->SetMarkerSize(1.8);
        h_theta_star[iem][ieta][irho]->SetMarkerStyle(21);
        //h_theta_star[irho]->Draw();
        TH1D* h_theta_star_ratio = (TH1D*)h_theta_star[iem][ieta][1]->Clone("ratio");
        TH1D* h_theta_star_clone = (TH1D*)h_theta_star[iem][ieta][irho]->Clone("clone");
        h_theta_star_ratio->Sumw2();
        h_theta_star_ratio->Divide(h_theta_star[iem][0][1]);
        h_theta_star_clone->Divide(h_theta_star_ratio);

        h_theta_star[iem][ieta][irho]->Fit(Func_rdlB, "NMI");
        double rhoB = Func_rdlB->GetParameter(1);
        double rhoBErr = Func_rdlB->GetParError(1);

        cout << "divided histograms" << endl;

        Func_rdl->SetParameter(0, h_theta_star_clone->GetBinContent(5));
        Func_rdl->SetParameter(1, 0.333);

        TString *Title = new TString(Form("fit/acceptanceQA/yield_pt_1_cent_5_rho%d%s_method%d.pdf",irho,cutoption.c_str(),method));
        
        TCanvas *c1 = new TCanvas(Form("c%d",irho),Form("c%d",irho),10,10,800,800);
        c1->cd();
        c1->SetFillColor(0);
        c1->SetGrid(0,0);
        c1->SetTitle(0);
        c1->SetBottomMargin(0.15);
        c1->SetLeftMargin(0.15);
        h_theta_star[iem][ieta][irho]->Draw("pE");
        h_theta_star[iem][ieta][irho]->Fit(Func_rhoOrig, "NMI");
        //Func_rhoOrig->Draw("same");
        //h_theta_star_clone->Fit(Func_rho, "WNMIER");
        h_theta_star_clone->Fit(Func_rdl, "NMI");
        if(fit2d == 1 && EP == 1) h_betatstarP[iem][ieta][irho]->Fit(Func_2D, "NMI");
        Func_rho->Draw("same");
 
        auto legend = new TLegend(0.1,0.7,0.48,0.9); 
        legend->AddEntry(Func_rho,"4th Order Fit","l");
        //legend->AddEntry(Func_rhoOrig,"2nd Order Fit","l");
        legend->Draw("same");
        
        c1->SaveAs(Title->Data());
        //double rho00final = 4./(1.+3.*resval)*(Func_rho->GetParameter(1)-1/3.)+1./3.;
         
        double rho_obs = Func_rdl->GetParameter(1);
        double rho_obs_err = Func_rdl->GetParError(1);
        
        // Error calculation
        double drhodobs = 4./(1.+3.*resval);  // calculate d(rho)/d(rho_obs)
        double drhodR = -12.*(rho_obs - 1./3.)/(1.+3.*resval)/(1.+3.*resval); // calculate d(rho)/d(R)

        double real_rho = 1./3. + 4./(1.+3.*resval)*(rho_obs - 1./3.);
        double real_rho_error = TMath::Sqrt((rho_obs_err*rho_obs_err)*(drhodobs*drhodobs));// + (resval*resval)*(drhodR*drhodR));
  
        double rho00final; 
        double rho00finalErr;
        if(EP == 1) { rho00final = real_rho;//Func_rho->GetParameter(1);
                      rho00finalErr = real_rho_error;//Func_rho->GetParError(1); 
                    }

        if(EP == 0) { rho00final = rhoreal[iem][ieta][irho];
                      rho00finalErr = rhorealErr[iem][ieta][irho]; }

        if(EP == 1 && fit2d == 1) { rho00final = Func_2D->GetParameter(1);
                                    rho00finalErr = Func_2D->GetParError(1); }

        if(EP == 0 && fit2d == 1) { rho00final = rhoreal[iem][ieta][irho];
                                    rho00finalErr = rhorealErr[iem][ieta][irho]; }


        cout << "rho after = " << rho00final << " +/- " << rho00finalErr << endl;
        cout << "about to fill After" << endl;
        g_rho00outAfterIn[ieta]->SetPoint(irho,rho00in[irho],rho00final);
        g_rho00outAfterIn[ieta]->SetPointError(irho,0.0,0.0,rho00finalErr,rho00finalErr);
        cout << "about to fill Before" << endl;
        g_rho00outBeforeIn[ieta]->SetPoint(irho,rho00in[irho],rhoB);
        g_rho00outBeforeIn[ieta]->SetPointError(irho,0.0,0.0,rhoBErr,rhoBErr);

        delete c1;  
        delete Func_rho;
        delete Func_rhoOrig;
        delete Title;
        delete h_theta_star_ratio;
        delete h_theta_star_clone;
      }
    }
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.20);
  c_play->SetGrid(0,0);
  c_play->SetTicks(1,1);
  c_play->cd();

  TH1F *h_play = new TH1F("h_play","h_play",100,-50,50);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetYaxis()->SetTitle("#rho_{00}");
  h_play->GetXaxis()->SetTitle("#rho_{00}^{input}");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetLimits(0.23,0.42);
  h_play->GetYaxis()->SetRangeUser(0.23,0.50);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");
  h_play->Draw("pE");

  TF1 *linear = new TF1("linear","x",0.23,0.42);
  //linear->FixParameter(0,resval);
  linear->Draw("same");

  int styleb[6] = {0,29,21,34,20,33};
  int stylea[6] = {0,30,25,28,24,27};
  int color[6] = {0,3,4,6,2,9};

  TLegend *leg1 = new TLegend(0.2,0.65,0.6,0.85);
  leg1->SetFillColor(10);
  leg1->SetBorderSize(0);
  for(int ieta = 1; ieta < 4; ieta++)
  {
    g_rho00outBeforeIn[ieta]->SetMarkerStyle(styleb[ieta]);
    g_rho00outBeforeIn[ieta]->SetMarkerColor(color[ieta]);
    g_rho00outBeforeIn[ieta]->SetLineColor(color[ieta]);
    g_rho00outBeforeIn[ieta]->SetMarkerSize(2.0);
    g_rho00outBeforeIn[ieta]->SetLineWidth(2);
    g_rho00outBeforeIn[ieta]->Draw("pE Same");

    g_rho00outAfterIn[ieta]->SetMarkerStyle(stylea[ieta]);
    g_rho00outAfterIn[ieta]->SetMarkerColor(color[ieta]);
    g_rho00outAfterIn[ieta]->SetLineColor(color[ieta]);
    g_rho00outAfterIn[ieta]->SetMarkerSize(2.0);
    g_rho00outAfterIn[ieta]->SetLineWidth(2);
    g_rho00outAfterIn[ieta]->Draw("pE Same");
    leg1->AddEntry(g_rho00outBeforeIn[ieta],Form("Before Correction |#eta|<%s",etatext[ieta].c_str()),"p");
    leg1->AddEntry(g_rho00outAfterIn[ieta],Form("After Correction |#eta|<%s",etatext[ieta].c_str()),"p");
  }

  leg1->AddEntry(linear,"#rho_{00} = #rho_{00}^{input}","l");
  leg1->Draw("same");

  string FigureName = Form("./figures/AcceptanceQA/AcceptanceQA_%s%s_method%d_BEFOREandAfter_FirstOrder.pdf",vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method);
  c_play->SaveAs(FigureName.c_str());


  //TCanvas *c_play3 = new TCanvas("c_play3","c_play3",10,10,800,800);
  //c_play3->SetLeftMargin(0.15);
  //c_play3->SetBottomMargin(0.15);
  //c_play3->SetGrid(0,0);
  //c_play3->SetTicks(1,1);
  //c_play3->cd();

  //h_play->SetTitle(Form("|#eta|<%s, %s",etatext[ieta].c_str(),ytext.c_str()));
  //h_play->DrawCopy("pE");

  //TF1 *linear = new TF1("linear","x",0.22,0.42);
  ////linear->FixParameter(0,resval);
  ////linear->Draw("same");

  //h_play->GetYaxis()->SetTitle("#rho_{00}^{reco}");
  //h_play->GetXaxis()->SetTitle("v_{2}^{input}");
  //TLegend *leg2 = new TLegend(0.6,0.175,0.9,0.425);
  //leg2->SetFillColor(10);
  //leg2->SetBorderSize(0);

  //int styleb[6] = {0,29,21,34,20,33};
  //int stylea[6] = {0,30,25,28,24,27};
  //int color[6] = {0,3,4,6,2,9}

  //for(int irho = 0; irho < 3; irho++)
  //{
  //  g_rho00RecoObs[irho]->SetMarkerStyle(styleb[irho+1]);
  //  g_rho00RecoObs[irho]->SetMarkerColor(color[irho+1]);
  //  g_rho00RecoObs[irho]->SetLineColor(color[irho+1]);
  //  g_rho00RecoObs[irho]->SetMarkerSize(2.0);
  //  g_rho00RecoObs[irho]->SetLineWidth(2);
  //  g_rho00RecoObs[irho]->Draw("pE Same");

  //  g_rho00outAfterIn[irho]->SetMarkerStyle(stylea[irho+1]);
  //  g_rho00outAfterIn[irho]->SetMarkerColor(color[irho+1]);
  //  g_rho00outAfterIn[irho]->SetLineColor(color[irho+1]);
  //  g_rho00outAfterIn[irho]->SetMarkerSize(2.0);
  //  g_rho00outAfterIn[irho]->SetLineWidth(2);
  //  g_rho00outAfterIn[irho]->Draw("pE Same");

  //  leg2->AddEntry(g_rho00RecoObs[irho],Form("#rho_{00}^{input} = %1.4f, Before Correction",rho00in[irho]),"p");
  //  leg2->AddEntry(g_rho00outAfterIn[irho],Form("#rho_{00}^{input} = %1.4f, After Correction",rho00in[irho]),"p");
  //}
  ////leg2->AddEntry(linear,"reco = input","l");
  //leg2->Draw("same");

  //std::string FigureName = Form("./figures/AcceptanceQA/V2Vary_6th_AcceptanceQAInput_EDGE_%s_RecoObs_v2%d_EPSmear%d_2d%d_%s%s_method%d_etamode%d_ieta%d.pdf",ytextfile.c_str(),v2,EP,fit2d,vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method,etamode,ieta);
  //c_play3->SaveAs(FigureName.c_str());


  //h_play->DrawCopy("pE");

  ////TF1 *linear = new TF1("linear","x",0.22,0.42);
  ////linear->FixParameter(0,resval);
  ////linear->Draw("same");

  //h_play->GetYaxis()->SetTitle("#rho_{00}^{reco}");
  //TLegend *leg2b = new TLegend(0.6,0.175,0.9,0.425);
  //leg2b->SetFillColor(10);
  //leg2b->SetBorderSize(0);

  //for(int irho = 0; irho < 3; irho++)
  //{
  //  g_rhoResCorr[irho]->SetMarkerStyle(styleb[irho+1]);
  //  g_rhoResCorr[irho]->SetMarkerColor(color[irho+1]);
  //  g_rhoResCorr[irho]->SetLineColor(color[irho+1]);
  //  g_rhoResCorr[irho]->SetMarkerSize(2.0);
  //  g_rhoResCorr[irho]->SetLineWidth(2);
  //  g_rhoResCorr[irho]->Draw("pE Same");

  //  //leg2b->AddEntry(g_rhoResCorr[iem],Form("%s, Before Res Correction",ytext[iem].c_str()),"p");
  //  leg2b->AddEntry(g_rhoResCorr[irho],Form("#rho_{00}^{input} = %1.4f, After Res Correction",rho00in[irho]),"p");
  //}
  ////leg2b->AddEntry(linear,"reco = input","l");
  //leg2b->Draw("same");

  //FigureName = Form("./figures/AcceptanceQA/VaryV2_6th_AcceptanceQAInput_ONLYRESCORR_EDGE_%s_RecoObs_v2%d_EPSmear%d_2d%d_%s%s_method%d_etamode%d_ieta%d.pdf",ytextfile.c_str(),v2,EP,fit2d,vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method,etamode,ieta);
  //c_play3->SaveAs(FigureName.c_str());

  ////TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,800,800);
  ////c_diff->SetLeftMargin(0.15);
  ////c_diff->SetBottomMargin(0.15);
  ////c_diff->SetGrid(0,0);
  ////c_diff->SetTicks(1,1);
  ////c_diff->cd();

  ////g_diff->GetHistogram()->GetXaxis()->SetTitle("Input #rho_{00}");
  ////g_diff->GetHistogram()->GetYaxis()->SetTitle("#rho_{00}'-#rho_{00}");
  ////g_diff->GetHistogram()->SetTitle("Difference Plot");

  ////g_diff->SetMarkerStyle(20);
  ////g_diff->SetMarkerColor(kBlack);
  ////g_diff->Draw("APE");

  ////FigureName = Form("./figures/AcceptanceQA/AcceptanceQADifference_%s%s_method%d.pdf",vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method);
  ////c_diff->SaveAs(FigureName.c_str());

  //TCanvas *c_diffI = new TCanvas("c_diffI","c_diffI",10,10,800,800);
  //c_diffI->SetLeftMargin(0.15);
  //c_diffI->SetBottomMargin(0.15);
  //c_diffI->SetGrid(0,0);
  //c_diffI->SetTicks(1,1);
  //c_diffI->cd();

  ////h_play->GetXaxis()->SetLimits(0.22,0.42);
  //h_play->GetYaxis()->SetRangeUser(0.18,0.58);
  //h_play->GetYaxis()->SetTitle("#rho_{00}^{reco}-#rho_{00}^{input}");

  //TLegend *leg3 = new TLegend(0.2,0.2,0.4,0.45);
  //leg3->SetFillColor(10);
  //leg3->SetBorderSize(0);

  //double max = -999;
  //double min =  999;
  //for(int irho = 0; irho < 3; irho++)
  //{
  //  cout << "YUP" << irho << endl;
  //  g_diffI[irho]->SetMarkerStyle(styleb[irho+1]);
  //  g_diffI[irho]->SetMarkerColor(color[irho+1]);
  //  g_diffI[irho]->SetLineColor(color[irho+1]);
  //  g_diffI[irho]->SetMarkerSize(2.0);
  //  g_diffI[irho]->SetLineWidth(2);
  //
  //  double tmin = g_diffI[irho]->GetMinimum();
  //  double tmax = g_diffI[irho]->GetMaximum();
  //  if(tmin < min) min = tmin;
  //  if(tmax > max) max = tmax;
  // 
  //  cout << "min = " << tmin << endl;
  //  cout << "max = " << tmax << endl;

  //  leg3->AddEntry(g_diffI[irho],Form("#rho_{00}^{input} = %1.4f",rho00in[irho]),"p");
  //  cout << "YEAH" << irho << endl;
  //}
  ////leg2->AddEntry(linear,"","l");
 
  //h_play->GetYaxis()->SetRangeUser(-0.15,0.05);
  //h_play->DrawCopy("pE");
  //for(int irho = 0; irho < 3; irho++)
  //{  
  //  cout << "YUP" << irho << endl;
  //  g_diffI[irho]->Draw("pE Same");
  //  cout << "YEAH" << irho << endl;
  //}
  //leg3->Draw("same");
  //TF1 *zero = new TF1("zero","0",-0.025,0.275);
  ////linear->FixParameter(0,resval);
  //zero->SetLineStyle(2);
  //zero->Draw("same");

  //FigureName = Form("./figures/AcceptanceQA/VaryV2_6th_AcceptanceQADifferenceFromIdeal_%s_EDGE_v2%d_EPSmear%d_2d%d_%s%s_method%d_etamode%d_ieta%d.pdf",ytextfile.c_str(),v2,EP,fit2d,vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method,etamode,ieta);
  //c_diffI->SaveAs(FigureName.c_str());


  //TLegend *leg3b = new TLegend(0.2,0.2,0.4,0.45);
  //leg3b->SetFillColor(10);
  //leg3b->SetBorderSize(0);

  //double max = -999;
  //double min =  999;
  //for(int irho = 0; irho < 3; irho++)
  //{
  //  cout << "YUP" << irho << endl;
  //  g_rhoResCorrDiff[irho]->SetMarkerStyle(styleb[irho+1]);
  //  g_rhoResCorrDiff[irho]->SetMarkerColor(color[irho+1]);
  //  g_rhoResCorrDiff[irho]->SetLineColor(color[irho+1]);
  //  g_rhoResCorrDiff[irho]->SetMarkerSize(2.0);
  //  g_rhoResCorrDiff[irho]->SetLineWidth(2);
  //
  //  double tmin = g_rhoResCorrDiff[irho]->GetMinimum();
  //  double tmax = g_rhoResCorrDiff[irho]->GetMaximum();
  //  if(tmin < min) min = tmin;
  //  if(tmax > max) max = tmax;
  // 
  //  cout << "min = " << tmin << endl;
  //  cout << "max = " << tmax << endl;

  //  leg3b->AddEntry(g_rhoResCorrDiff[irho],Form("#rho_{00}^{input} = %1.4f",rho00in[irho]),"p");
  //  cout << "YEAH" << irho << endl;
  //}
  ////leg2->AddEntry(linear,"","l");
 
  //h_play->GetYaxis()->SetRangeUser(-0.05,0.05);
  //h_play->DrawCopy("pE");
  //for(int irho = 0; irho < 3; irho++)
  //{  
  //  cout << "YUP" << irho << endl;
  //  g_rhoResCorrDiff[irho]->Draw("pE Same");
  //  cout << "YEAH" << irho << endl;
  //}
  //leg3b->Draw("same");
  ////TF1 *zero = new TF1("zero","0",0.22,0.42);
  ////linear->FixParameter(0,resval);
  //zero->SetLineStyle(2);
  //zero->Draw("same");

  //FigureName = Form("./figures/AcceptanceQA/VaryV2_6th_AcceptanceQADifferenceFromIdeal_ONLYRESCORR_EDGE_%s_v2%d_EPSmear%d_2d%d_%s%s_method%d_etamode%d_ieta%d.pdf",ytextfile.c_str(),v2,EP,fit2d,vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method,etamode,ieta);
  //c_diffI->SaveAs(FigureName.c_str());
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

double Func2D(double *x_val, double *par) {

  double ts = x_val[0];
  double b = x_val[1];

  double sts = TMath::Sin(ts);
  double cts = TMath::Cos(ts);
  double sb  = TMath::Sin(b);
  double c2b = TMath::Cos(2.*b);

  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double G = par[3];
  double R = par[4];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double acceptance = (1. + F*(sts*sts)*(sb*sb) + G*(sts*sts*sts*sts)*(sb*sb*sb*sb));
  double rawyield   = (1. + As*cts*cts + Bs*sts*sts*c2b);
 
  return N*rawyield*acceptance*(sts);

}

double Func2D_wB(double *x_val, double *par) {

  double ts = x_val[0];
  double b = x_val[1];

  double sts = TMath::Sin(ts);
  double cts = TMath::Cos(ts);
  double sb  = TMath::Sin(b);
  double c2b = TMath::Cos(2.*b);

  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double G = par[3];
  double R = par[4];
  double Bn = par[5];

  double B = -2.*TMath::Sqrt(2)*Bn/(1.-rho);
  double A = (3.*rho-1.)/(1.-rho);
  double As = (A*(1.+3.*R)+B*(3.-3.*R))/(4.+A*(1.-R)+B*(-1.+R));
  double Bs = (A*(1.-R)+B*(3.+R))/(4.+A*(1.-R)+B*(-1.+R));

  double acceptance = (1. + F*sts*sts*sb*sb + G*sts*sts*sts*sts*sb*sb*sb*sb);
  double rawyield   = (1. + As*cts*cts + Bs*sts*sts*c2b);
 
  return N*rawyield*acceptance*(-1.0/sts);
}

double Func6th(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double F = par[1];
  double G = par[2];
  double H = par[3];

  double result = 16. + (8.*F+6.*G+5.*H) - (8.*F+12.*G+15.*H)*CosTheta*CosTheta + (6.*G+15.*H)*CosTheta*CosTheta*CosTheta*CosTheta - 5.*H*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta;

  return N/16.*result;

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

  return N/64.*(order0 + order2 + order4 + order6 + order8);

}
