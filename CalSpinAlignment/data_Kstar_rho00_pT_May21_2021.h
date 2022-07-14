
//pT dependence of kstar Rho00: Centrality 20-60%
// pair rapidity |y| < 1.0
namespace KStarBESI{

//11 GeV
const int kstar_Npt_11   = 2;
Float_t kstar_pt_11[kstar_Npt_11]    = {1.25, 3.25};
Float_t kstar_rho00_pT_11[kstar_Npt_11] = {0.302849, 0.373099};
Float_t kstar_rho00_pT_11_stat[kstar_Npt_11] = {0.030352, 0.0399644};
Float_t kstar_rho00_pT_11_syst[kstar_Npt_11] = {0.0132217, 0.0164683};


//14 GeV
const int kstar_Npt_14   = 3;
Float_t kstar_pt_14[kstar_Npt_14]    = {1.25, 3.25};
Float_t kstar_rho00_pT_14[kstar_Npt_14] = {0.311237, 0.304138};
Float_t kstar_rho00_pT_14_stat[kstar_Npt_14] = {0.0255322, 0.0354243};
Float_t kstar_rho00_pT_14_syst[kstar_Npt_14] = {0.017875, 0.0247494};


//19 GeV
const int kstar_Npt_19   = 3;
Float_t kstar_pt_19[kstar_Npt_19]    = {1.25, 3.25};
Float_t kstar_rho00_pT_19[kstar_Npt_19] = {0.349793, 0.334475};
Float_t kstar_rho00_pT_19_stat[kstar_Npt_19] = {0.0188304, 0.0199586};
Float_t kstar_rho00_pT_19_syst[kstar_Npt_19] = {0.00804442, 0.0160283};

//27 GeV
const int kstar_Npt_27   = 3;
Float_t kstar_pt_27[kstar_Npt_27]    = {1.25, 2.25, 4 };
Float_t kstar_rho00_pT_27[kstar_Npt_27] = {0.313043, 0.321985, 0.309001};
Float_t kstar_rho00_pT_27_stat[kstar_Npt_27] = {0.0134185, 0.0131507, 0.0691408};
Float_t kstar_rho00_pT_27_syst[kstar_Npt_27] = {0.00655509, 0.00697146, 0.0371412};

//39 GeV
const int kstar_Npt_39   = 3;
Float_t kstar_pt_39[kstar_Npt_39]    = {1.25, 2.25, 4};
Float_t kstar_rho00_pT_39[kstar_Npt_39] = {0.331466, 0.330098, 0.368276};
Float_t kstar_rho00_pT_39_stat[kstar_Npt_39] = {0.0090676, 0.00785387, 0.0376544};
Float_t kstar_rho00_pT_39_syst[kstar_Npt_39] = {0.00475486, 0.00639724, 0.0217396};

//54 GeV
const int kstar_Npt_54   = 6;
Float_t kstar_pt_54[kstar_Npt_54]    = {1.25, 1.75, 2.25, 2.75, 3.5, 4.5};
Float_t kstar_rho00_pT_54[kstar_Npt_54] = {0.320713, 0.343106, 0.354951, 0.36704, 0.326496, 0.303558};
Float_t kstar_rho00_pT_54_stat[kstar_Npt_54] = {0.00594582, 0.00588814, 0.00721296, 0.011111, 0.0155364, 0.0452151};
Float_t kstar_rho00_pT_54_syst[kstar_Npt_54] = {0.00494916, 0.00437556, 0.00791626, 0.00842355, 0.00973361, 0.0277798};


//200 GeV 
const int kstar_Npt_200   = 6;
Float_t kstar_pt_200[kstar_Npt_200]    = {1.25, 1.75, 2.25, 2.75, 3.5, 4.5};
Float_t kstar_rho00_pT_200[kstar_Npt_200] = {0.319872, 0.328844, 0.351253, 0.345452, 0.344547, 0.381185};
Float_t kstar_rho00_pT_200_stat[kstar_Npt_200] = {0.00617656, 0.00565801, 0.00594025, 0.00756617, 0.00844954, 0.0179879};
Float_t kstar_rho00_pT_200_syst[kstar_Npt_200] = {0.00891925, 0.00705815, 0.00664678, 0.00739795, 0.00760375, 0.0119649};



//x-error
Float_t kstar_pt_err_11[kstar_Npt_11] = {0.08, 0.08};
Float_t kstar_pt_err_14[kstar_Npt_14] = {0.08, 0.08, 0.08};
Float_t kstar_pt_err_19[kstar_Npt_19] = {0.08, 0.08, 0.08};
Float_t kstar_pt_err_27[kstar_Npt_27] = {0.08, 0.08, 0.08};
Float_t kstar_pt_err_39[kstar_Npt_39] = {0.08, 0.08, 0.08};
Float_t kstar_pt_err_54[kstar_Npt_54] = {0.08, 0.08, 0.08, 0.08, 0.08, 0.08};
Float_t kstar_pt_err_200[kstar_Npt_200] = {0.08, 0.08, 0.08, 0.08, 0.08, 0.08};

};
