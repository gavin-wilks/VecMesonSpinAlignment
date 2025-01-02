#include "Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"


void acceptanceQA_includeBandBeta_4th_Res_v2test_v2vary(const int energy = 4, const int pid = 0, std::string res = "0p4", double resval = 0.4, int cut = 1, int method = 0, int v2 = 1, int fit2d = 0, int EP = 1, int etamode = 0, int ieta = 1) {
  
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

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
  std::string ytextfile = "yabs1";
  std::string ytext = "0.8<y<1.0";
  //std::string ytextfile = "yabs1";
  //std::string ytext = "-1.0<y<1.0";
  //std::string ytextfile = "yabs1";
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
  
  TGraphAsymmErrors *g_rho00outBefore[4] ;
  TGraphAsymmErrors *g_rho00outAfter[4] ;
  TGraphAsymmErrors *g_rho00outBeforeIn[4];
  TGraphAsymmErrors *g_rho00outAfterIn[4];
  TGraphAsymmErrors *g_rho00RecoObs[4];

  double Fval[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Ferr[11][6][6] = {0.0};
  double Gval[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double Gerr[11][6][6] = {0.0};
  double FvalRP[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double FerrRP[11][6][6] = {0.0};
  double GvalRP[11][6][6] = {0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  double GerrRP[11][6][6] = {0.0};
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

  double b2a[11][6][6] = {0.0};
  double b2b[11][6][6] = {0.0};
  double b2c[11][6][6] = {0.0};
  double b4a[11][6][6] = {0.0};
  double b4b[11][6][6] = {0.0};
  double b4c[11][6][6] = {0.0};

  double b2a_noEta[11][6][6] = {0.0};
  double b2b_noEta[11][6][6] = {0.0};
  double b2c_noEta[11][6][6] = {0.0};
  double b4a_noEta[11][6][6] = {0.0};
  double b4b_noEta[11][6][6] = {0.0};
  double b4c_noEta[11][6][6] = {0.0};

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
  
 int  rhoorder[3]={1,0,2};

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

    for(int iem = 0; iem < 11; iem++)
    {

      cout << "iem = " << iem << ",     irho = " << irho << endl;
      MCFiles[iem][irho] = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/%s/Acceptance/ResolutionTesting/ResolutionTest_0.4_EtaCut%d_%s_v2%d_withBeta_MANYPLOTS_EDGE_v2val%s/McAcceptanceOutput_pt1_energy%d_pid%d_cent4_EtaMode_0_nrho%d.root",vmsa::mPID[pid].c_str(),etamode,ytextfile.c_str(),v2,v2text[iem].c_str(),energy,pid,rho00val[irho]),"READ");
     
      MCFiles[iem][irho]->Print();

      cout << "created graphs" << endl;

      for(int i = 0; i < 6; i++)
      {
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
        p_cos2betaP[iem][i][irho]->Print();      
        histname = Form("p_cos4betaP_%d",i);
        p_cos4betaP[iem][i][irho]    = (TProfile*)((TProfile*)MCFiles[iem][irho]->Get(histname.c_str()))->Clone(Form("p_cos4betaP_%d_%d_%d",iem,i,irho));
        p_cos4betaP[iem][i][irho]->Print();      

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


        h_theta_star[iem][i][irho] = (TH1D*)((TH1D*)MCFiles[iem][irho]->Get(Form("h_theta_star_%d",i)))->Clone(Form("h_theta_star_%d_%d_%d",iem,irho,i));
        h_theta_star_RP[iem][i][irho] = (TH1D*)((TH1D*)MCFiles[iem][irho]->Get(Form("h_theta_star_RP_%d",i)))->Clone(Form("h_theta_star_RP_%d_%d_%d",iem,irho,i));
     }
       
      //h_out1[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_out1"))->Clone(Form("h_out1_%d",irho));
      //h_out2[irho] = (TH1D*)((TH1D*)MCFiles[irho]->Get("h_out2"))->Clone(Form("h_out2_%d",irho));
        
      cout << "Loaded all histograms" << endl;

      TF2 *Func_2Dbefore = new TF2("Func_2Dbefore",Func2D,0,TMath::Pi(),0,2.*TMath::Pi(),5);
      TF1 *Func_1D = new TF1("Func_1D",FuncAFG,0,1,5);
      //TF2 *Func_2Dbefore= new TF2("Func_2Dbefore","[0]*(1.-[1]+(3.*[1]-1)*(cos(x))^2)+y*0",0,TMath::Pi(),0,2.*TMath::Pi());
      //double rhoBefore;
      //double rhoBeforeErr;
      if(EP == 0 && fit2d == 1)
      {
        Func_2Dbefore->SetParameter(0,h_betatstar[iem][0][irho]->GetBinContent(10,10));
        Func_2Dbefore->SetParameter(1,1./3.);
        Func_2Dbefore->FixParameter(2,0.0);
        Func_2Dbefore->FixParameter(3,0.0);
        Func_2Dbefore->FixParameter(4,1.0);
        h_betatstar[iem][0][irho]->Fit(Func_2Dbefore,"NERQ");
        rhoBefore[iem][irho]    = Func_2Dbefore->GetParameter(1);
        rhoBeforeErr[iem][irho] = Func_2Dbefore->GetParError(1);
        TCanvas *c1 = new TCanvas();
        c1->SetFillColor(0);
        c1->SetGrid(0,0);
        c1->SetTitle(0);
        c1->SetBottomMargin(0.15);
        c1->SetLeftMargin(0.15);
        h_betatstar[iem][0][irho]->Draw("colz");
        //Func_beta->Draw("same"); 
        //legend->AddEntry(Func_A,"4th Order Fit","l");
        //legend->AddEntry(Func_Orig,"2nd Order Fit","l");
        c1->SaveAs(Form("FValues/yabs1/BEFORE_EP%d_iem%d_pt1_rho%d_eta%d%s_method%d_2D.pdf",EP,iem,irho,0,cutoption.c_str(),method));
      }
      if(EP == 1 && fit2d == 1)
      { 
        Func_2Dbefore->SetParameter(0,h_betatstarP[iem][0][irho]->GetBinContent(10,10));
        Func_2Dbefore->SetParameter(1,1./3.);
        Func_2Dbefore->FixParameter(2,0.0);
        Func_2Dbefore->FixParameter(3,0.0);
        Func_2Dbefore->FixParameter(4,resval);
        h_betatstarP[iem][0][irho]->Fit(Func_2Dbefore,"NERQ");
        rhoBefore[iem][irho]    = Func_2Dbefore->GetParameter(1);
        rhoBeforeErr[iem][irho] = Func_2Dbefore->GetParError(1);
        TCanvas *c1 = new TCanvas();
        c1->SetFillColor(0);
        c1->SetGrid(0,0);
        c1->SetTitle(0);
        c1->SetBottomMargin(0.15);
        c1->SetLeftMargin(0.15);
        h_betatstarP[iem][0][irho]->Draw("colz");
        //Func_beta->Draw("same"); 
        //legend->AddEntry(Func_A,"4th Order Fit","l");
        //legend->AddEntry(Func_Orig,"2nd Order Fit","l");
        c1->SaveAs(Form("FValues/yabs1/BEFORE_EP%d_iem%d_pt1_rho%d_eta%d%s_method%d_2D.pdf",EP,iem,irho,0,cutoption.c_str(),method));
      }  
      if(EP == 0 && fit2d == 0)
      {
        Func_1D->SetParameter(0,h_theta_star_RP[iem][0][irho]->GetBinContent(1));
        Func_1D->SetParameter(1,1./3.);
        Func_1D->FixParameter(2,0.0);
        Func_1D->FixParameter(3,0.0);
        Func_1D->FixParameter(4,1.0);
        h_theta_star_RP[iem][0][irho]->Fit(Func_1D,"NERQ");
        rhoBefore[iem][irho]    = Func_1D->GetParameter(1);
        rhoBeforeErr[iem][irho] = Func_1D->GetParError(1);
        TCanvas *c1 = new TCanvas();
        c1->SetFillColor(0);
        c1->SetGrid(0,0);
        c1->SetTitle(0);
        c1->SetBottomMargin(0.15);
        c1->SetLeftMargin(0.15);
        h_theta_star_RP[iem][0][irho]->Draw("pE");
        //Func_beta->Draw("same"); 
        //legend->AddEntry(Func_A,"4th Order Fit","l");
        //legend->AddEntry(Func_Orig,"2nd Order Fit","l");
        c1->SaveAs(Form("FValues/yabs1/BEFORE_EP%d_iem%d_pt1_rho%d_eta%d%s_method%d_1D.pdf",EP,iem,irho,0,cutoption.c_str(),method));
      }
      if(EP == 1 && fit2d == 0)
      {
        Func_1D->SetParameter(0,h_theta_star[iem][0][irho]->GetBinContent(1));
        Func_1D->SetParameter(1,1./3.);
        Func_1D->FixParameter(2,0.0);
        Func_1D->FixParameter(3,0.0);
        Func_1D->FixParameter(4,resval);
        h_theta_star[iem][0][irho]->Fit(Func_1D,"NERQ");
        rhoBefore[iem][irho]    = Func_1D->GetParameter(1);
        rhoBeforeErr[iem][irho] = Func_1D->GetParError(1);
        TCanvas *c1 = new TCanvas();
        c1->SetFillColor(0);
        c1->SetGrid(0,0);
        c1->SetTitle(0);
        c1->SetBottomMargin(0.15);
        c1->SetLeftMargin(0.15);
        h_theta_star[iem][0][irho]->Draw("pE");
        //Func_beta->Draw("same"); 
        //legend->AddEntry(Func_A,"4th Order Fit","l");
        //legend->AddEntry(Func_Orig,"2nd Order Fit","l");
        c1->SaveAs(Form("FValues/yabs1/BEFORE_EP%d_iem%d_pt1_rho%d_eta%d%s_method%d_1D.pdf",EP,iem,irho,0,cutoption.c_str(),method));
      }
      cout << "rhoBefore = " << rhoBefore[iem][irho] << " +/- " << rhoBeforeErr[iem][irho] << endl;
      g_rhoResCorr[irho]->SetPoint(iem,v2val[iem],rhoBefore[iem][irho]);
      g_rhoResCorr[irho]->SetPointError(iem,0.0,0.0,rhoBeforeErr[iem][irho],rhoBeforeErr[iem][irho]);
      g_rhoResCorrDiff[irho]->SetPoint(iem,v2val[iem],rhoBefore[iem][irho]-rho00in[irho]);
      g_rhoResCorrDiff[irho]->SetPointError(iem,0.0,0.0,rhoBeforeErr[iem][irho],rhoBeforeErr[iem][irho]);
      


      TF1 *Func_rho = new TF1("Func_rho","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);
      TF2 *Func_rho2D= new TF2("Func_rho2D",Func2D,0,TMath::Pi(),0,2.*TMath::Pi(),5);
      //TF2 *Func_rho2D= new TF2("Func_rho2D",Func2D_wB,0,TMath::Pi(),0,2.*TMath::Pi(),6);
      TF1 *Func_A = new TF1("Func_A",Func4th,0,1,9);
      Func_A->SetLineColor(kBlue);
      TF1 *Func_Orig = new TF1("Func_Orig","[0]*(1 + [1]*x*x)",0,1);
      //Func_Orig->SetLineColor(kRed);
      TF1 *Func_Theta = new TF1("Func_Theta","[0]*(1 + [1]*x*x + [2]*x*x*x*x + [3]*x*x*x*x*x*x)",0,1);
      //TF1 *Func_AFG = new TF1("Func_AFG",FuncAFG,0,1,3);

      TF2 *Func_beta = new TF2("Func_beta","[0] * (1. + [1]*(sin(x)*sin(y))^2 + [2]*(sin(x)*sin(y))^4)",0,TMath::Pi(),0.0,2.0*TMath::Pi());

      if(EP == 0 && fit2d == 0)
      {
        Func_rho->SetParameter(0,h_theta_star_RP[iem][ieta][irho]->GetBinContent(1));
        Func_rho->SetParameter(1,1./3.);
        h_theta_star_RP[iem][ieta][irho]->Fit(Func_rho,"NERQ");
        rhoobs[iem][ieta][irho] = Func_rho->GetParameter(1);
        rhoobsErr[iem][ieta][irho] = Func_rho->GetParError(1);
      }
      if(EP == 1 && fit2d == 0)
      { 
        Func_rho->SetParameter(0,h_theta_star[iem][ieta][irho]->GetBinContent(1));
        Func_rho->SetParameter(1,1./3.);
        h_theta_star[iem][ieta][irho]->Fit(Func_rho,"NERQ");
        rhoobs[iem][ieta][irho] = Func_rho->GetParameter(1);
        rhoobsErr[iem][ieta][irho] = Func_rho->GetParError(1);
      }
      if(EP == 0 && fit2d == 1)
      {
        Func_rho2D->SetParameter(0,h_betatstar[iem][ieta][irho]->GetBinContent(10,10));
        Func_rho2D->SetParameter(1,1./3.);
        Func_rho2D->FixParameter(2,0.0);
        Func_rho2D->FixParameter(3,0.0);
        Func_rho2D->FixParameter(4,1.0);
        h_betatstar[iem][ieta][irho]->Fit(Func_rho2D,"NERQ");
        rhoobs[iem][ieta][irho] = Func_rho2D->GetParameter(1);
        rhoobsErr[iem][ieta][irho] = Func_rho2D->GetParError(1);
      }
      if(EP == 1 && fit2d == 1)
      { 
        Func_rho2D->SetParameter(0,h_betatstarP[iem][ieta][irho]->GetBinContent(10,10));
        Func_rho2D->SetParameter(1,1./3.);
        Func_rho2D->FixParameter(2,0.0);
        Func_rho2D->FixParameter(3,0.0);
        Func_rho2D->FixParameter(4,1.0);
        h_betatstarP[iem][ieta][irho]->Fit(Func_rho2D,"NERQ");
        rhoobs[iem][ieta][irho] = Func_rho2D->GetParameter(1);
        rhoobsErr[iem][ieta][irho] = Func_rho2D->GetParError(1);
      }

      cout << "rho input = " << rho00in[irho] << ",    rho observed = " << rhoobs[iem][ieta][irho] << endl;
      //Func_rho->SetParameter(0,h_theta_star_RP[iem][ieta][irho]->GetBinContent(1));
      //Func_rho->SetParameter(1,1./3.);
      //h_theta_star_RP[iem][ieta][irho]->Fit(Func_rho,"NERQ");
      //rhoreal[iem][ieta][irho] = Func_rho->GetParameter(1);
      //rhorealErr[iem][ieta][irho] = Func_rho->GetParError(1);
      //if(cut == 0) g_rho00outBefore[ieta]->SetPoint(irho,rhoobs[iem][ieta][irho],rhoreal[iem][ieta][irho]);
      //if(cut == 0) g_rho00outBefore[ieta]->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));
      //if(cut == 0) g_rho00outBeforeIn[ieta]->SetPoint(irho,rho00in[irho],rhoreal[iem][ieta][irho]);
      //if(cut == 0) g_rho00outBeforeIn[ieta]->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));


      //TH1D *h_theta_clone = (TH1D*)h_theta[iem][ieta][irho]->Clone("h_theta_clone");
      //h_theta_clone->Sumw2();
      //h_theta_before[iem][ieta][irho]->Sumw2();
      //h_theta_clone->Divide(h_theta_before[iem][ieta][irho]);
      //Func_Theta->SetParameter(0,h_theta_clone->GetBinContent(1));
      //Func_Theta->SetParameter(1,0);
      //h_theta_clone->Fit(Func_Theta,"ER");


      //if(cut == 1 && method == 2){
      //  Fval[iem][ieta][irho] = Func_Theta->GetParameter(1);
      //  Ferr[iem][ieta][irho] = Func_Theta->GetParError(1);
      //  Gval[iem][ieta][irho] = Func_Theta->GetParameter(2);
      //  Gerr[iem][ieta][irho] = Func_Theta->GetParError(2);
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

      TH2D *h_betatstarP_ratio = (TH2D*)h_betatstarP[iem][ieta][irho]->Clone("h_betatstarP_ratio");
      h_betatstarP_ratio->Sumw2();
      h_betatstarP[iem][0][irho]->Sumw2();
      h_betatstarP_ratio->Divide(h_betatstarP[iem][0][irho]); 
      Func_beta->SetParameter(0,h_betatstarP_ratio->GetBinContent(1));
      Func_beta->SetParameter(1,0.0);
      Func_beta->SetParameter(2,0.0);
      //Func_A->FixParameter(2,0.0);
      //Func_beta->SetParameter(0,h_theta_star_clone->GetBinContent(1));
      //Func_beta->SetParameter(1,0);
      //h_theta_star_clone->Fit(Func_Orig,"ER");

      TH1D *h_theta_star_clone = (TH1D*)h_theta_star[iem][ieta][irho]->Clone("h_theta_star_clone");
      h_theta_star_clone->Sumw2();
      h_theta_star[iem][0][irho]->Sumw2();
      h_theta_star_clone->Divide(h_theta_star[iem][0][irho]); 
      Func_A->SetParameter(0,h_theta_star_clone->GetBinContent(1));
      Func_A->SetParameter(1,0);
      Func_A->SetParameter(2,0.0);
      //Func_A->FixParameter(2,0.0);
      Func_Orig->SetParameter(0,h_theta_star_clone->GetBinContent(1));
      Func_Orig->SetParameter(1,0);
      h_theta_star_clone->Fit(Func_Orig,"ER");

      TF1 *cos2bP = new TF1("cos2bP","[0]+[1]*x+[2]*x*x",0,1);
      TF1 *cos4bP = new TF1("cos4bP","[0]+[1]*x+[2]*x*x",0,1);
      TF1 *cos2bP_noEta = new TF1("cos2bP","[0]+[1]*x+[2]*x*x",0,1);
      TF1 *cos4bP_noEta = new TF1("cos4bP","[0]+[1]*x+[2]*x*x",0,1);
  
      p_cos2betaP[iem][ieta][irho]->Fit(cos2bP,"NMIER");
      p_cos4betaP[iem][ieta][irho]->Fit(cos4bP,"NMIER");
      p_cos2betaP[iem][0][irho]->Fit(cos2bP_noEta,"NMIER");
      p_cos4betaP[iem][0][irho]->Fit(cos4bP_noEta,"NMIER");
  
      b2a[iem][ieta][irho] = cos2bP->GetParameter(0);
      b2b[iem][ieta][irho] = cos2bP->GetParameter(1);
      b2c[iem][ieta][irho] = cos2bP->GetParameter(2);
      b4a[iem][ieta][irho] = cos4bP->GetParameter(0);
      b4b[iem][ieta][irho] = cos4bP->GetParameter(1);
      b4c[iem][ieta][irho] = cos4bP->GetParameter(2);

      b2a_noEta[iem][ieta][irho] = cos2bP_noEta->GetParameter(0);
      b2b_noEta[iem][ieta][irho] = cos2bP_noEta->GetParameter(1);
      b2c_noEta[iem][ieta][irho] = cos2bP_noEta->GetParameter(2);
      b4a_noEta[iem][ieta][irho] = cos4bP_noEta->GetParameter(0);
      b4b_noEta[iem][ieta][irho] = cos4bP_noEta->GetParameter(1);
      b4c_noEta[iem][ieta][irho] = cos4bP_noEta->GetParameter(2);

      if(cut == 1 && method == 0){
        if(fit2d){h_betatstarP_ratio->Fit(Func_beta,"ER");
        FvalB[iem][ieta][irho] = Func_beta->GetParameter(1);
        FerrB[iem][ieta][irho] = Func_beta->GetParError(1);
        GvalB[iem][ieta][irho] = Func_beta->GetParameter(2);
        GerrB[iem][ieta][irho] = Func_beta->GetParError(2);
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
       }
        Func_A->FixParameter(3,b2a_noEta[iem][ieta][irho]);
        Func_A->FixParameter(4,b2b_noEta[iem][ieta][irho]);
        Func_A->FixParameter(5,b2c_noEta[iem][ieta][irho]);
        Func_A->FixParameter(6,b4a_noEta[iem][ieta][irho]);
        Func_A->FixParameter(7,b4b_noEta[iem][ieta][irho]);
        Func_A->FixParameter(8,b4c_noEta[iem][ieta][irho]);
        h_theta_star_clone->Fit(Func_A,"NMIER");
        Fval[iem][ieta][irho] = Func_A->GetParameter(1);
        Ferr[iem][ieta][irho] = Func_A->GetParError(1);
        Gval[iem][ieta][irho] = Func_A->GetParameter(2);
        Gerr[iem][ieta][irho] = Func_A->GetParError(2);
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
        c1->SaveAs(Form("FValues/yabs1/pt1_rho%d_eta%d%s_method%d.pdf",irho,ieta,cutoption.c_str(),method));
        delete c1;
      }
 
      FvalOrig[irho] = Func_Orig->GetParameter(1);
      FerrOrig[irho] = Func_Orig->GetParError(1);

      TH2D *h_betatstar_ratio = (TH2D*)h_betatstar[iem][ieta][irho]->Clone("h_betatstar_ratio");
      h_betatstar_ratio->Sumw2();
      h_betatstar[iem][0][irho]->Sumw2();
      h_betatstar_ratio->Divide(h_betatstar[iem][0][irho]); 
      Func_beta->SetParameter(0,h_betatstar_ratio->GetBinContent(1));
      Func_beta->SetParameter(1,0.0);
      Func_beta->SetParameter(2,0.0);

      TH1D *h_theta_star_clone_RP = (TH1D*)h_theta_star_RP[iem][ieta][irho]->Clone("h_theta_star_clone_RP");
      //TH1D *h_theta_star_clone_RP = (TH1D*)h_theta_star[iem][ieta][irho]->Clone("h_theta_star_clone_RP");
      h_theta_star_clone_RP->Sumw2();
      h_theta_star_RP[iem][0][irho]->Sumw2();
      h_theta_star_clone_RP->Divide(h_theta_star_RP[iem][0][irho]); 
      Func_A->SetParameter(0,h_theta_star_clone_RP->GetBinContent(1));
      Func_A->SetParameter(1,0);
      Func_A->SetParameter(2,0.0);
      //Func_A->FixParameter(2,0.0);
      Func_Orig->SetParameter(0,h_theta_star_clone_RP->GetBinContent(1));
      Func_Orig->SetParameter(1,0);
      h_theta_star_clone_RP->Fit(Func_Orig,"ER");
      h_theta_star_clone_RP->Fit(Func_A,"NMIER");
      if(fit2d)h_betatstar_ratio->Fit(Func_beta,"ER");

      if(cut == 1 && method == 1){
         if(fit2d){
        FvalB[iem][ieta][irho] = Func_beta->GetParameter(1);
        FerrB[iem][ieta][irho] = Func_beta->GetParError(1);
        GvalB[iem][ieta][irho] = Func_beta->GetParameter(2);
        GerrB[iem][ieta][irho] = Func_beta->GetParError(2);
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
        }
        Fval[iem][ieta][irho] = Func_A->GetParameter(1);
        Ferr[iem][ieta][irho] = Func_A->GetParError(1);
        Gval[iem][ieta][irho] = Func_A->GetParameter(2);
        Gerr[iem][ieta][irho] = Func_A->GetParError(2);
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
        c1->SaveAs(Form("FValues/yabs1/pt1_rho%d_eta%d%s_method%d.pdf",irho,ieta,cutoption.c_str(),method));
        delete c1;
      }

      TF1 *Func_AFG = new TF1("Func_AFG",FuncAFG,0,1,5);
      TF2 *Func_2D = new TF2("Func_2D",Func2D,0.0,TMath::Pi(),0.0,2.0*TMath::Pi(),5);
      //TF2 *Func_2D = new TF2("Func_2D",Func2D_wB,0.0,TMath::Pi(),0.0,2.0*TMath::Pi(),6);
      //TF1 *Func_AFG = new TF1("Func_AFG",FuncAFGH,0,1,6);
      FvalRP[iem][ieta][irho] = Func_A->GetParameter(1);
      FerrRP[iem][ieta][irho] = Func_A->GetParError(1);
      GvalRP[iem][ieta][irho] = Func_A->GetParameter(2);
      GerrRP[iem][ieta][irho] = Func_A->GetParError(2);
      FvalBRP[iem][ieta][irho] = Func_beta->GetParameter(1);
      FerrBRP[iem][ieta][irho] = Func_beta->GetParError(1);
      GvalBRP[iem][ieta][irho] = Func_beta->GetParameter(2);
      GerrBRP[iem][ieta][irho] = Func_beta->GetParError(2);
      if(fit2d == 0)
      {
        Func_AFG->SetParameter(0,h_theta_star_RP[iem][ieta][irho]->GetBinContent(1));
        Func_AFG->SetParameter(1,1./3.);
        Func_AFG->FixParameter(2,FvalRP[iem][ieta][1]);//Fval[irho]);
        Func_AFG->FixParameter(3,GvalRP[iem][ieta][1]);//Gval[irho]);
        //Func_AFG->FixParameter(2,Func_A->GetParameter(1));//Fval[irho]);
        //Func_AFG->FixParameter(3,Func_A->GetParameter(2));//Gval[irho]);
        //Func_AFG->FixParameter(2,FvalRP[0][ieta][1]);//Fval[irho]);
        //Func_AFG->FixParameter(3,GvalRP[0][ieta][1]);//Gval[irho]);
        Func_AFG->FixParameter(4, 1.0);//Res_12);
        h_theta_star_RP[iem][ieta][irho]->Fit(Func_AFG,"NMIER");
        rhoreal[iem][ieta][irho] = Func_AFG->GetParameter(1);
        rhorealErr[iem][ieta][irho] = Func_AFG->GetParError(1);
      }
      if(fit2d == 1)
      {
        Func_2D->SetParameter(0,h_betatstar[iem][ieta][irho]->GetBinContent(10,10));
        Func_2D->SetParameter(1,1./3.);
        Func_2D->FixParameter(2,FvalBRP[iem][ieta][1]);//Fval[irho]);
        Func_2D->FixParameter(3,GvalBRP[iem][ieta][1]);//Gval[irho]);
        Func_2D->FixParameter(4, 1.0);//Res_12);
        h_betatstar[iem][ieta][irho]->Fit(Func_2D,"ER");
        rhoreal[iem][ieta][irho] = Func_2D->GetParameter(1);
        rhorealErr[iem][ieta][irho] = Func_2D->GetParError(1);
      }
      cout << "rho00in = " << rho00in[irho] << ",    rho00 from RP with acceptance correction = " << rhoreal[iem][ieta][irho] << endl;
      if(cut == 1) g_rho00outBefore[irho]->SetPoint(irho,rhoobs[iem][ieta][irho],rhoreal[iem][ieta][irho]);
      if(cut == 1) g_rho00outBefore[irho]->SetPointError(irho,0.0,0.0,rhorealErr[iem][ieta][irho],rhorealErr[iem][ieta][irho]);
      if(cut == 1) g_rho00outBeforeIn[irho]->SetPoint(irho,rho00in[irho],rhoreal[iem][ieta][irho]);
      if(cut == 1) g_rho00outBeforeIn[irho]->SetPointError(irho,0.0,0.0,rhorealErr[iem][ieta][irho],rhorealErr[iem][ieta][irho]);


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

  for(int iem = 0; iem < 11; iem++)
  {
    for(int irho = 0; irho < 3; irho++)
    {
    //FfromFp[irho] = -Fval[irho]/(2.+Fval[irho]);
    //double err = TMath::Abs((Fval[irho]/(2.-Fval[irho])/(2.-Fval[irho])-1./(2.-Fval[irho]))*Ferr[irho]);
    //cout << "F from F* = " << FfromFp[irho] << " +/- " << err << "      F from Orig Method = " << FvalOrig[irho] << " +/- " << FerrOrig[irho] << endl;
    //cout << "F from theta = " << Fpval[irho] << " +/- " << Fperr[ierr] << "      G from theta = " << Gpval[irho] << " +/- " << Gperr[irho] << endl;
      cout << "irho = " << irho << ",   iem = " << iem << endl;
      cout << "F = " << Fval[iem][ieta][irho] << " +/- " << Ferr[iem][ieta][irho] << "           F from 2D = " << FvalB[iem][ieta][irho] << " +/- " << FerrB[iem][ieta][irho] << endl;
      cout << "G = " << Gval[iem][ieta][irho] << " +/- " << Gerr[iem][ieta][irho] << "           G from 2D = " << GvalB[iem][ieta][irho] << " +/- " << GerrB[iem][ieta][irho] << endl;
    } 
  }

  TGraphAsymmErrors *g_diff[6]; 
  TGraphAsymmErrors *g_diffI[6];

  for(int irho = 0; irho < 3; irho++) 
  {
    g_diff[irho]  = new TGraphAsymmErrors();
    g_diffI[irho] = new TGraphAsymmErrors();
    for(int iem = 0; iem < 11; iem++)
    {
      //TF1 *Func_rho = new TF1("Func_rho",FuncAFG,0,1,5);
      TF1 *Func_rho = new TF1("Func_rho",FuncAFG,0,1,18);
      //TF2 *Func_2D = new TF2("Func_2D",Func2D_wB,0.0,TMath::Pi(),0.0,2.0*TMath::Pi(),6);
      TF2 *Func_2D = new TF2("Func_2D",Func2D,0.0,TMath::Pi(),0.0,2.0*TMath::Pi(),5);
      Func_rho->SetLineColor(kBlue);
      TF1 *Func_rhoOrig = new TF1("Func_rhoOrig",FuncAD,0,1,4);
      Func_rhoOrig->SetLineColor(kRed);
      //TF1 *Func_rho = new TF1("Func_rho",FuncAF_Updated,0,1,4);
      //TF1 *Func_rdl = new TF1("Func_rdl","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);

      h_theta_star[iem][ieta][irho]->GetXaxis()->SetTitleOffset(1.2);
      h_theta_star[iem][ieta][irho]->GetXaxis()->SetTitle("cos#theta*");
      h_theta_star[iem][ieta][irho]->GetYaxis()->SetTitle("yield");
      h_theta_star[iem][ieta][irho]->GetYaxis()->SetTitleOffset(1.0);
      h_theta_star[iem][ieta][irho]->SetMarkerColor(2);
      h_theta_star[iem][ieta][irho]->SetMarkerSize(1.8);
      h_theta_star[iem][ieta][irho]->SetMarkerStyle(21);
      //h_theta_star[irho]->Draw();
      Func_rho->SetParameter(0, h_theta_star[iem][ieta][irho]->GetBinContent(5));
      Func_rho->SetParameter(1, 0.333);
      //Func_rho->SetParLimits(1, 0.0, 1.0);
      //g_rho00outAfter->SetPoint(irho,rho00in[irho],Func_rho->GetParameter(1));
      //g_rho00outAfter->SetPointError(irho,0.0,0.0,Func_rho->GetParError(1),Func_rho->GetParError(1));
      //Func_rho->FixParameter(2, -1.57811);//Fval[1]);
      //Func_rho->FixParameter(3, 0.649864);//Gval[1]);
      Func_rho->FixParameter(2, Fval[iem][ieta][irho]);
      Func_rho->FixParameter(3, Gval[iem][ieta][irho]);
      //Func_rho->FixParameter(2, Fval[0][ieta][1]);
      //Func_rho->FixParameter(3, Gval[0][ieta][1]);

      //if(fit2d == 1) Func_rho->FixParameter(2, FvalB[iem][ieta][1]);
      //if(fit2d == 1) Func_rho->FixParameter(3, GvalB[iem][ieta][1]);
      //Func_rho->FixParameter(4, Hval[irho]);
      Func_rho->FixParameter(4, resval);//Res_12);
      Func_rho->FixParameter(5,  b2a_noEta[iem][ieta][irho]);
      Func_rho->FixParameter(6,  b2b_noEta[iem][ieta][irho]);
      Func_rho->FixParameter(7,  b2c_noEta[iem][ieta][irho]);
      Func_rho->FixParameter(8,  b4a_noEta[iem][ieta][irho]);
      Func_rho->FixParameter(9,  b4b_noEta[iem][ieta][irho]);
      Func_rho->FixParameter(10, b4c_noEta[iem][ieta][irho]);
      Func_rho->FixParameter(11, 0.0);
      Func_rho->FixParameter(12, 0.0);
      Func_rho->FixParameter(13, 0.0);
      Func_rho->FixParameter(14, 0.0);
      Func_rho->FixParameter(15, b2a[iem][ieta][irho]);
      Func_rho->FixParameter(16, b2b[iem][ieta][irho]);
      Func_rho->FixParameter(17, b2c[iem][ieta][irho]);

      //Func_rho->FixParameter(4, 1.0);//Res_12);

      if(fit2d == 1) 
      {
        Func_2D->SetParameter(0,h_betatstarP[iem][ieta][irho]->GetBinContent(10,10));
        Func_2D->SetParameter(1,0.3333);
        Func_2D->FixParameter(2,FvalB[iem][ieta][1]);
        Func_2D->FixParameter(3,GvalB[iem][ieta][1]);
        Func_2D->FixParameter(4,resval);
      }
      //Func_rho->FixParameter(3, 1.0);//Res_12);

      //Func_rhoOrig->SetParameter(0, h_theta_star[iem][ieta][irho]->GetBinContent(5));
      //Func_rhoOrig->SetParameter(1, 0.333);
      //Func_rhoOrig->FixParameter(2, FvalOrig[irho]);
      //Func_rhoOrig->FixParameter(3, resval);//Res_12);


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
      h_theta_star[iem][ieta][irho]->Draw("pE");
      //h_theta_star[iem][ieta][irho]->Fit(Func_rhoOrig, "NMI");
      //Func_rhoOrig->Draw("same");
      h_theta_star[iem][ieta][irho]->Fit(Func_rho, "NMIER");
      if(fit2d == 1 && EP == 1) h_betatstarP[iem][ieta][irho]->Fit(Func_2D, "NMI");
      Func_rho->Draw("same");
 
      auto legend = new TLegend(0.1,0.7,0.48,0.9); 
      legend->AddEntry(Func_rho,"4th Order Fit","l");
      //legend->AddEntry(Func_rhoOrig,"2nd Order Fit","l");
      legend->Draw("same");
      
      c1->SaveAs(Title->Data());
      //double rho00final = 4./(1.+3.*resval)*(Func_rho->GetParameter(1)-1/3.)+1./3.;
      double rho00final; 
      double rho00finalErr;
      if(EP == 1) { rho00final = Func_rho->GetParameter(1);
                    rho00finalErr = Func_rho->GetParError(1); }

      if(EP == 0) { rho00final = rhoreal[iem][ieta][irho];
                    rho00finalErr = rhorealErr[iem][ieta][irho]; }

      if(EP == 1 && fit2d == 1) { rho00final = Func_2D->GetParameter(1);
                                  rho00finalErr = Func_2D->GetParError(1); }

      if(EP == 0 && fit2d == 1) { rho00final = rhoreal[iem][ieta][irho];
                                  rho00finalErr = rhorealErr[iem][ieta][irho]; }


      cout << "rho after = " << rho00final << " +/- " << rho00finalErr << endl;
      g_rho00outAfter[irho]->SetPoint(iem,rhoobs[iem][ieta][irho],rho00final);
      g_rho00outAfter[irho]->SetPointError(iem,0.0,0.0,rho00finalErr,rho00finalErr);
      g_rho00outAfterIn[irho]->SetPoint(iem,v2val[iem],rho00final);
      g_rho00outAfterIn[irho]->SetPointError(iem,0.0,0.0,rho00finalErr,rho00finalErr);
      g_rho00RecoObs[irho]->SetPoint(iem,v2val[iem],rhoobs[iem][ieta][irho]);
      g_rho00RecoObs[irho]->SetPointError(iem,0.0,0.0,rhoobsErr[iem][ieta][irho],rhoobsErr[iem][ieta][irho]);
      cout << "rho input = " << rho00in[irho] << ",    rho observed = " << rhoobs[iem][ieta][irho] << endl;

      cout << "rho00in = " << rho00in[irho] << ",    rho00 from RP with acceptance correction = " << rhoreal[iem][ieta][irho] << " +/- " << rhorealErr[iem][ieta][irho] << endl;
      g_diff[irho]->SetPoint(iem,v2val[iem],rho00final-rhoreal[iem][ieta][irho]);
      g_diffI[irho]->SetPoint(iem,v2val[iem],rho00final-rho00in[irho]);
      double uncorrErr = TMath::Sqrt(rho00finalErr*rho00finalErr+rhorealErr[iem][ieta][irho]*rhorealErr[iem][ieta][irho]);
      double corrErr = TMath::Sqrt(fabs(rho00finalErr*rho00finalErr-rhorealErr[iem][ieta][irho]*rhorealErr[iem][ieta][irho]));
      std::cout << "uncorrelated Error = " << uncorrErr << ",      correlated Error = " << corrErr << std::endl;
      g_diff[irho]->SetPointError(iem,0.0,0.0,corrErr,corrErr);
      g_diffI[irho]->SetPointError(iem,0.0,0.0,rho00finalErr,rho00finalErr);

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

  TH1F *h_play = new TH1F("h_play","h_play",100,-50,50);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetYaxis()->SetTitle("#rho_{00}^{reco}");
  h_play->GetXaxis()->SetTitle("v_{2}^{input}");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetLimits(-0.025,0.275);
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

  h_play->SetTitle(Form("|#eta|<%s, %s",etatext[ieta].c_str(),ytext.c_str()));
  h_play->DrawCopy("pE");

  TF1 *linear = new TF1("linear","x",0.22,0.42);
  //linear->FixParameter(0,resval);
  //linear->Draw("same");

  h_play->GetYaxis()->SetTitle("#rho_{00}^{reco}");
  h_play->GetXaxis()->SetTitle("v_{2}^{input}");
  TLegend *leg2 = new TLegend(0.6,0.175,0.9,0.425);
  leg2->SetFillColor(10);
  leg2->SetBorderSize(0);

  int styleb[6] = {0,29,21,34,20,33};
  int stylea[6] = {0,30,25,28,24,27};
  int color[6] = {0,3,4,6,2,9};

  for(int irho = 0; irho < 3; irho++)
  {
    g_rho00RecoObs[irho]->SetMarkerStyle(styleb[irho+1]);
    g_rho00RecoObs[irho]->SetMarkerColor(color[irho+1]);
    g_rho00RecoObs[irho]->SetLineColor(color[irho+1]);
    g_rho00RecoObs[irho]->SetMarkerSize(2.0);
    g_rho00RecoObs[irho]->SetLineWidth(2);
    g_rho00RecoObs[irho]->Draw("pE Same");

    g_rho00outAfterIn[irho]->SetMarkerStyle(stylea[irho+1]);
    g_rho00outAfterIn[irho]->SetMarkerColor(color[irho+1]);
    g_rho00outAfterIn[irho]->SetLineColor(color[irho+1]);
    g_rho00outAfterIn[irho]->SetMarkerSize(2.0);
    g_rho00outAfterIn[irho]->SetLineWidth(2);
    g_rho00outAfterIn[irho]->Draw("pE Same");

    leg2->AddEntry(g_rho00RecoObs[irho],Form("#rho_{00}^{input} = %1.4f, Before Correction",rho00in[irho]),"p");
    leg2->AddEntry(g_rho00outAfterIn[irho],Form("#rho_{00}^{input} = %1.4f, After Correction",rho00in[irho]),"p");
  }
  //leg2->AddEntry(linear,"reco = input","l");
  leg2->Draw("same");

  std::string FigureName = Form("./figures/AcceptanceQA/V2Vary_AcceptanceQAInput_EDGE_%s_RecoObs_v2%d_EPSmear%d_2d%d_%s%s_method%d_etamode%d_ieta%d.pdf",ytextfile.c_str(),v2,EP,fit2d,vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method,etamode,ieta);
  c_play3->SaveAs(FigureName.c_str());


  h_play->DrawCopy("pE");

  //TF1 *linear = new TF1("linear","x",0.22,0.42);
  //linear->FixParameter(0,resval);
  //linear->Draw("same");

  h_play->GetYaxis()->SetTitle("#rho_{00}^{reco}");
  TLegend *leg2b = new TLegend(0.6,0.175,0.9,0.425);
  leg2b->SetFillColor(10);
  leg2b->SetBorderSize(0);

  for(int irho = 0; irho < 3; irho++)
  {
    g_rhoResCorr[irho]->SetMarkerStyle(styleb[irho+1]);
    g_rhoResCorr[irho]->SetMarkerColor(color[irho+1]);
    g_rhoResCorr[irho]->SetLineColor(color[irho+1]);
    g_rhoResCorr[irho]->SetMarkerSize(2.0);
    g_rhoResCorr[irho]->SetLineWidth(2);
    g_rhoResCorr[irho]->Draw("pE Same");

    //leg2b->AddEntry(g_rhoResCorr[iem],Form("%s, Before Res Correction",ytext[iem].c_str()),"p");
    leg2b->AddEntry(g_rhoResCorr[irho],Form("#rho_{00}^{input} = %1.4f, After Res Correction",rho00in[irho]),"p");
  }
  //leg2b->AddEntry(linear,"reco = input","l");
  leg2b->Draw("same");

  FigureName = Form("./figures/AcceptanceQA/VaryV2_AcceptanceQAInput_ONLYRESCORR_EDGE_%s_RecoObs_v2%d_EPSmear%d_2d%d_%s%s_method%d_etamode%d_ieta%d.pdf",ytextfile.c_str(),v2,EP,fit2d,vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method,etamode,ieta);
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

  //h_play->GetXaxis()->SetLimits(0.22,0.42);
  h_play->GetYaxis()->SetRangeUser(0.18,0.58);
  h_play->GetYaxis()->SetTitle("#rho_{00}^{reco}-#rho_{00}^{input}");

  TLegend *leg3 = new TLegend(0.2,0.2,0.4,0.45);
  leg3->SetFillColor(10);
  leg3->SetBorderSize(0);

  double max = -999;
  double min =  999;
  for(int irho = 0; irho < 3; irho++)
  {
    cout << "YUP" << irho << endl;
    g_diffI[irho]->SetMarkerStyle(styleb[irho+1]);
    g_diffI[irho]->SetMarkerColor(color[irho+1]);
    g_diffI[irho]->SetLineColor(color[irho+1]);
    g_diffI[irho]->SetMarkerSize(2.0);
    g_diffI[irho]->SetLineWidth(2);
  
    double tmin = g_diffI[irho]->GetMinimum();
    double tmax = g_diffI[irho]->GetMaximum();
    if(tmin < min) min = tmin;
    if(tmax > max) max = tmax;
   
    cout << "min = " << tmin << endl;
    cout << "max = " << tmax << endl;

    leg3->AddEntry(g_diffI[irho],Form("#rho_{00}^{input} = %1.4f",rho00in[irho]),"p");
    cout << "YEAH" << irho << endl;
  }
  //leg2->AddEntry(linear,"","l");
 
  h_play->GetYaxis()->SetRangeUser(-0.15,0.05);
  h_play->DrawCopy("pE");
  for(int irho = 0; irho < 3; irho++)
  {  
    cout << "YUP" << irho << endl;
    g_diffI[irho]->Draw("pE Same");
    cout << "YEAH" << irho << endl;
  }
  leg3->Draw("same");
  TF1 *zero = new TF1("zero","0",-0.025,0.275);
  //linear->FixParameter(0,resval);
  zero->SetLineStyle(2);
  zero->Draw("same");

  FigureName = Form("./figures/AcceptanceQA/VaryV2_AcceptanceQADifferenceFromIdeal_%s_EDGE_v2%d_EPSmear%d_2d%d_%s%s_method%d_etamode%d_ieta%d.pdf",ytextfile.c_str(),v2,EP,fit2d,vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method,etamode,ieta);
  c_diffI->SaveAs(FigureName.c_str());


  TLegend *leg3b = new TLegend(0.2,0.2,0.4,0.45);
  leg3b->SetFillColor(10);
  leg3b->SetBorderSize(0);

  double max = -999;
  double min =  999;
  for(int irho = 0; irho < 3; irho++)
  {
    cout << "YUP" << irho << endl;
    g_rhoResCorrDiff[irho]->SetMarkerStyle(styleb[irho+1]);
    g_rhoResCorrDiff[irho]->SetMarkerColor(color[irho+1]);
    g_rhoResCorrDiff[irho]->SetLineColor(color[irho+1]);
    g_rhoResCorrDiff[irho]->SetMarkerSize(2.0);
    g_rhoResCorrDiff[irho]->SetLineWidth(2);
  
    double tmin = g_rhoResCorrDiff[irho]->GetMinimum();
    double tmax = g_rhoResCorrDiff[irho]->GetMaximum();
    if(tmin < min) min = tmin;
    if(tmax > max) max = tmax;
   
    cout << "min = " << tmin << endl;
    cout << "max = " << tmax << endl;

    leg3b->AddEntry(g_rhoResCorrDiff[irho],Form("#rho_{00}^{input} = %1.4f",rho00in[irho]),"p");
    cout << "YEAH" << irho << endl;
  }
  //leg2->AddEntry(linear,"","l");
 
  h_play->GetYaxis()->SetRangeUser(-0.05,0.05);
  h_play->DrawCopy("pE");
  for(int irho = 0; irho < 3; irho++)
  {  
    cout << "YUP" << irho << endl;
    g_rhoResCorrDiff[irho]->Draw("pE Same");
    cout << "YEAH" << irho << endl;
  }
  leg3b->Draw("same");
  //TF1 *zero = new TF1("zero","0",0.22,0.42);
  //linear->FixParameter(0,resval);
  zero->SetLineStyle(2);
  zero->Draw("same");

  FigureName = Form("./figures/AcceptanceQA/VaryV2_AcceptanceQADifferenceFromIdeal_ONLYRESCORR_EDGE_%s_v2%d_EPSmear%d_2d%d_%s%s_method%d_etamode%d_ieta%d.pdf",ytextfile.c_str(),v2,EP,fit2d,vmsa::mBeamEnergy[energy].c_str(),cutoption.c_str(),method,etamode,ieta);
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

//double Func4th(double *x_val, double *par) {
//
//  double CosTheta = x_val[0];
//  double N = par[0];
//  double F = par[1];
//  double G = par[2];
//  double b2a = par[3];
//  double b2b = par[4];
//  double b2c = par[5];
//  double b4a = par[6];
//  double b4b = par[7];
//  double b4c = par[8];
//
//  double cos2b = b2a + b2b * CosTheta + b2c * CosTheta * CosTheta;
//  double cos4b = b4a + b4b * CosTheta + b4c * CosTheta * CosTheta;
//
//  double order0 = 1. + F/2. + 3.*G/8. - F/2.*cos2b - G/2.*cos2b + G/8.*cos4b;
//  double order2 = (-F/2. - 3.*G/4. + F/2.*cos2b + G*cos2b - G/4.*cos4b) * CosTheta * CosTheta;
//  double order4 = (3.*G/8. - G/2.*cos2b + G/8.*cos4b) * CosTheta * CosTheta * CosTheta * CosTheta;
//
//  return N * ( order0 + order2 + order4 );
//
//}

double Func4th(double *x_val, double *par) {

  double TwoPi = 1.0;//TMath::Pi()*2.0;
  double CosTheta = x_val[0];
  double N = par[0];
  double F = par[1];
  double G = par[2];
  double b2a = par[3];
  double b2b = par[4];
  double b2c = par[5];
  double b4a = par[6];
  double b4b = par[7];
  double b4c = par[8];

  double cos2b = TwoPi*(b2a + b2b * CosTheta + b2c * CosTheta * CosTheta);
  double cos4b = TwoPi*(b4a + b4b * CosTheta + b4c * CosTheta * CosTheta);

  double order0 = 1. + F/2. + 3.*G/8. - F/2.*cos2b - G/2.*cos2b + G/8.*cos4b;
  double order2 = (-F/2. - 3.*G/4. + F/2.*cos2b + G*cos2b - G/4.*cos4b) * CosTheta * CosTheta;
  double order4 = (3.*G/8. - G/2.*cos2b + G/8.*cos4b) * CosTheta * CosTheta * CosTheta * CosTheta;

  return N * ( order0 + order2 + order4 );

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

//double FuncAFG(double *x_val, double *par) {
//
//  double CosTheta = x_val[0];
//  double N = par[0];
//  double rho = par[1];
//  double F = par[2];
//  double G = par[3];
//  double R = par[4];
//  double b2a = par[5];
//  double b2b = par[6];
//  double b2c = par[7];
//  double b4a = par[8];
//  double b4b = par[9];
//  double b4c = par[10];
//  double b2b4a = par[11];
//  double b2b4b = par[12];
//  double b2b4c = par[13];
//  double rhop1n1 = par[14];
//
//  double cos2b = b2a + b2b * CosTheta + b2c * CosTheta * CosTheta;
//  double cos4b = b4a + b4b * CosTheta + b4c * CosTheta * CosTheta;
//  double cos2bcos4b = cos2b*cos4b; // temporary "solution" 
//  //double cos2bcos4b = b2b4a + b2b4b * CosTheta + b2b4c * CosTheta * CosTheta; 
//
//  double A = (3.*rho-1.)/(1.-rho);
//  double B = 2.*(-rhop1n1)/(1.-rho);
//  //double As = A*(1.+3.*R)/(4.+A*(1.-R));
//  //double Bs = A*(1.-R)/(4.+A*(1.-R));
//  double As = (A*(1.+3.*R)+B*(3.-3.*R))/(4.+A*(1.-R)+B*(-1.+R));
//  double Bs = (A*(1.-R)+B*(3.+R))/(4.+A*(1.-R)+B*(-1.+R));
//
//  double order0 = 1. + F/2. + 3.*G/8. + Bs*cos2b - F/2.*cos2b + Bs*F/2.*cos2b - G/2.*cos2b + 3.*Bs*G/8.*cos2b - Bs*F/4.*(1. + cos4b) 
//                  - Bs*G/4.*(1. + cos4b) + G/8.*cos4b + Bs*G/8.*cos2bcos4b;
//  double order2 = (As - F/2. + As*F/2. - 3.*G/4. + 3.*As*G/8. - Bs*cos2b + F/2.*cos2b - As*F/2.*cos2b - Bs*F*cos2b
//                   + G*cos2b - As*G/2.*cos2b - 9.*Bs*G/8.*cos2b + Bs*F/2.*(1. + cos4b) + 3.*Bs*G/4.*(1 + cos4b)
//                   - G/4.*cos4b + As*G/8.*cos4b - 3.*Bs*G/8.*cos2bcos4b) * CosTheta * CosTheta;
//  double order4 = (-As*F/2. + 3.*G/8. - 3.*As*G/4. + As*F/2.*cos2b + Bs*F/2.*cos2b - G/2.*cos2b + A*G*cos2b 
//                   + 9.*Bs*G/8.*cos2b - Bs*F/4.*(1 + cos4b) - 3.*Bs*G/4.*(1 + cos4b) + G/8.*cos4b - As*G/4.*cos4b 
//                   + 3.*Bs*G/8.*cos2bcos4b)* CosTheta * CosTheta * CosTheta * CosTheta;
//  double order6 = (3.*As*G/8. - As*G/2.*cos2b - 3.*Bs*G/8.*cos2b + Bs*G/4.*(1 + cos4b) + As*G/8.*cos4b 
//                   - Bs*G/8.*cos2bcos4b) * CosTheta * CosTheta * CosTheta * CosTheta * CosTheta * CosTheta;
//
//  double result = order0 + order2 + order4 + order6;
//
//  return N*result;
//
//}

double FuncAFG(double *x_val, double *par) {

  double TwoPi = 1.0;//TMath::Pi()*2.0;
  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double G = par[3];
  double R = par[4];
  double b2a = par[5];
  double b2b = par[6];
  double b2c = par[7];
  double b4a = par[8];
  double b4b = par[9];
  double b4c = par[10];
  double b2b4a = par[11];
  double b2b4b = par[12];
  double b2b4c = par[13];
  double rhop1n1 = par[14];
  double b02a = par[15];
  double b02b = par[16];
  double b02c = par[17];

  double cos2b0 = (b02a + b02b * CosTheta + b02c * CosTheta * CosTheta); // Multiply by 2pi to convert average to integral
  double cos2b = TwoPi*(b2a + b2b * CosTheta + b2c * CosTheta * CosTheta);
  double cos4b = TwoPi*(b4a + b4b * CosTheta + b4c * CosTheta * CosTheta);
  double cos2bcos4b = cos2b*cos4b; // temporary "solution" 
  //double cos2bcos4b = b2b4a + b2b4b * CosTheta + b2b4c * CosTheta * CosTheta; 

  double A = (3.*rho-1.)/(1.-rho);
  double B = 2.*(-rhop1n1)/(1.-rho);
  //double As = A*(1.+3.*R)/(4.+A*(1.-R));
  //double Bs = A*(1.-R)/(4.+A*(1.-R));
  double As = (A*(1.+3.*R)+B*(3.-3.*R))/(4.+A*(1.-R)+B*(-1.+R));
  double Bs = (A*(1.-R)+B*(3.+R))/(4.+A*(1.-R)+B*(-1.+R));

  double order0 = 1. + F/2. + 3.*G/8. - F/2.*cos2b - G/2.*cos2b + G/8.*cos4b + Bs*cos2b0 + Bs*F/2.*cos2b0 + 3.*Bs*G/8.*cos2b0 - Bs*F/2.*cos2b*cos2b0 
                  - Bs*G/2.*cos2b*cos2b0 + Bs*G/8.*cos2b0*cos4b;
  double order2 = (As - F/2. + As*F/2. - 3.*G/4. + 3.*As*G/8. + F/2.*cos2b - As*F/2.*cos2b
                   + G*cos2b - As*G/2.*cos2b - G/4.*cos4b + As*G/8.*cos4b - Bs*cos2b0 - Bs*F*cos2b0 - 9.*Bs*G/8.*cos2b0 + Bs*F*cos2b*cos2b0 
                   + 3.*Bs*G/2.*cos2b*cos2b0 - 3.*Bs*G/8.*cos2b0*cos4b) * CosTheta * CosTheta;
  double order4 = (-As*F/2. + 3.*G/8. - 3.*As*G/4. + As*F/2.*cos2b - G/2.*cos2b + A*G*cos2b + G/8.*cos4b - As*G/4.*cos4b  
                   + Bs*F/2.*cos2b0 + 9.*Bs*G/8.*cos2b0 - Bs*F/2.*cos2b*cos2b0 - 3.*Bs*G/2*cos2b*cos2b0
                   + 3.*Bs*G/8.*cos2b0*cos4b)* CosTheta * CosTheta * CosTheta * CosTheta;
  double order6 = (3.*As*G/8. - As*G/2.*cos2b + As*G/8.*cos4b - 3.*Bs*G/8.*cos2b0 + Bs*G/2.*cos2b*cos2b0
                   - Bs*G/8.*cos2b0*cos4b) * CosTheta * CosTheta * CosTheta * CosTheta * CosTheta * CosTheta;

  double result = order0 + order2 + order4 + order6;

  return N*result;

}
