#include "TMath.h"
#include "TLatex.h"
#include "TLine.h"

Double_t generalized_t_pdf(Double_t *x, Double_t *par) {
    Double_t A = par[0];       // Normalization factor
    Double_t nu = par[1];      // Degrees of freedom
    Double_t mu = par[2];      // Location parameter
    Double_t sigma = par[3];   // Scale parameter
    Double_t val = x[0];

    Double_t numerator = TMath::Gamma((nu + 1) / 2);
    Double_t denominator = TMath::Sqrt(nu * TMath::Pi()) * sigma * TMath::Gamma(nu / 2);
    Double_t term = TMath::Power(1 + TMath::Power((val - mu) / sigma, 2) / nu, - (nu + 1) / 2);

    return A * (numerator / denominator * term);
}

Double_t generalized_t_pdf_3part(Double_t *x, Double_t *par) {
    Double_t A1 = par[0];       // Normalization factor
    Double_t nu1 = par[1];      // Degrees of freedom
    Double_t mu1 = par[2];      // Location parameter
    Double_t sigma1 = par[3];   // Scale parameter

    Double_t A2 = par[4];       // Normalization factor
    Double_t nu2 = par[5];      // Degrees of freedom
    Double_t mu2 = par[6];      // Location parameter
    Double_t sigma2 = par[7];   // Scale parameter

    Double_t A3 = par[8];       // Normalization factor
    Double_t nu3 = par[9];      // Degrees of freedom
    Double_t mu3 = par[10];      // Location parameter
    Double_t sigma3 = par[11];   // Scale parameter

    Double_t val = x[0];

    Double_t numerator1 = TMath::Gamma((nu1 + 1) / 2);
    Double_t denominator1 = TMath::Sqrt(nu1 * TMath::Pi()) * sigma1 * TMath::Gamma(nu1 / 2);
    Double_t term1 = TMath::Power(1 + TMath::Power((val - mu1) / sigma1, 2) / nu1, - (nu1 + 1) / 2);

    Double_t numerator2 = TMath::Gamma((nu2 + 1) / 2);
    Double_t denominator2 = TMath::Sqrt(nu2 * TMath::Pi()) * sigma2 * TMath::Gamma(nu2 / 2);
    Double_t term2 = TMath::Power(1 + TMath::Power((val - mu2) / sigma2, 2) / nu2, - (nu2 + 1) / 2);

    Double_t numerator3 = TMath::Gamma((nu3 + 1) / 2);
    Double_t denominator3 = TMath::Sqrt(nu3 * TMath::Pi()) * sigma3 * TMath::Gamma(nu3 / 2);
    Double_t term3 = TMath::Power(1 + TMath::Power((val - mu3) / sigma3, 2) / nu3, - (nu3 + 1) / 2);

    return A1 * (numerator1 / denominator1 * term1) + A2 * (numerator2 / denominator2 * term2) + A3 * (numerator3 / denominator3 * term3) ;
}

double ThreeParticleGaussian2D(double *x, double *p)
{
  double GaussX1 = p[0] *TMath::Exp(-TMath::Power(x[0]-p[1], 2)/(2*p[2] *p[2]));
  double GaussX2 = p[3] *TMath::Exp(-TMath::Power(x[0]-p[4], 2)/(2*p[5] *p[5]));
  double GaussX3 = p[6] *TMath::Exp(-TMath::Power(x[0]-p[7], 2)/(2*p[8] *p[8]));
  double GaussY1 = p[9] *TMath::Exp(-TMath::Power(x[1]-p[10],2)/(2*p[11]*p[11]));
  double GaussY2 = p[12]*TMath::Exp(-TMath::Power(x[1]-p[13],2)/(2*p[14]*p[14]));
  double GaussY3 = p[15]*TMath::Exp(-TMath::Power(x[1]-p[16],2)/(2*p[17]*p[17]));

  return (GaussX1+GaussX2+GaussX3)*(GaussY1+GaussY2+GaussY3);

}

double SingleParticleGaussian2D(double *x, double *p)
{
  double GaussX1 = p[0]*TMath::Exp(-TMath::Power(x[0]-p[1],2)/(2*p[2]*p[2]));
  double GaussY1 = p[3]*TMath::Exp(-TMath::Power(x[1]-p[4],2)/(2*p[5]*p[5]));

  return (GaussX1)*(GaussY1);

}

double ThreeParticleGaussian(double *x, double *p)
{
  double GaussX1 = p[0]*TMath::Exp(-TMath::Power(x[0]-p[1],2)/(2*p[2]*p[2]));
  double GaussX2 = p[3]*TMath::Exp(-TMath::Power(x[0]-p[4],2)/(2*p[5]*p[5]));
  double GaussX3 = p[6]*TMath::Exp(-TMath::Power(x[0]-p[7],2)/(2*p[8]*p[8]));

  return (GaussX1+GaussX2+GaussX3);
}

double TwoParticleGaussian(double *x, double *p)
{
  double GaussX1 = p[0]*TMath::Exp(-TMath::Power(x[0]-p[1],2)/(2*p[2]*p[2]));
  double GaussX2 = p[3]*TMath::Exp(-TMath::Power(x[0]-p[4],2)/(2*p[5]*p[5]));

  return (GaussX1+GaussX2);
}

double SingleParticleGaussian(double *x, double *p)
{
  double GaussX1 = p[0]*TMath::Exp(-TMath::Power(x[0]-p[1],2)/(2*p[2]*p[2]));

  return (GaussX1);
}



double SpinDensity2Dcos(double *x, double *p)
{
  return p[5]*((1.-p[0])+(3.*p[0] -1.)*x[0]*x[0]
               -sqrt(2)*p[1]*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*cos(x[1])
               +sqrt(2)*p[2]*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*sin(x[1])
               -2.*p[3]*(1-x[0]*x[0])*cos(2.*x[1])
               +2.*p[4]*(1-x[0]*x[0])*sin(2.*x[1]));
}

double SpinDensity2Dcos_GlobalFromHelicity(double *x, double *p)
{
  return p[5]*((1.-p[0])+(3.*p[0] -1.)*(TMath::Sqrt(1-x[0]*x[0])*cos(x[1])+x[0])*(TMath::Sqrt(1-x[0]*x[0])*cos(x[1])+x[0])
               -sqrt(2)*p[1]*2*(x[0]-x[0]*x[0]*x[0])*cos(x[1])
               +sqrt(2)*p[2]*2*(x[0]-x[0]*x[0]*x[0])*sin(x[1])
               -2.*p[3]*(1-x[0]*x[0])*cos(x[1])*cos(x[1])
               +2.*p[4]*(1-x[0]*x[0])*sin(x[1])*cos(x[1]));
}

double SpinDensity2DcosSeparated(double *x, double *p)
{
  double lhs = (p[1]*p[1] + p[2]*p[2] + p[3]*p[3] + p[4]*p[4]);
  double rhs = (p[5]*p[5] + p[6]*p[6]);
 
  double penalty = 0.0;
  //if( lhs < rhs )
  //{
  //  penalty = 1e3 * (rhs - lhs) * (rhs - lhs);
  //}
  if( lhs < rhs )
  {
    return 1e20;
  }

  double funcValue =  p[5]*((1.-p[0])+(3.*p[0] -1.)*x[0]*x[0]
               -sqrt(2)*(p[1]-p[2])*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*cos(x[1])
               +sqrt(2)*(p[3]-p[4])*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*sin(x[1])
               -2.*p[5]*(1-x[0]*x[0])*cos(2.*x[1])
               +2.*p[6]*(1-x[0]*x[0])*sin(2.*x[1]));

  return funcValue + penalty;
}

double SpinDensity(double *x_val, double *par)
{
  double x = x_val[0];
  double rho00 = par[0];
  double Norm = par[1]; // Norm = 0.75 <= sampling | Norm from fitting <= extract rho00

  double dNdCosThetaStar = Norm*((1.0-rho00)+(3.0*rho00-1)*x*x);

  return dNdCosThetaStar;
}

double SpinDensityCosCosH(double *x_val, double *par)
{
  double Cos = x_val[0];
  double CosH = x_val[1];
  double rho00 = par[0];
  double rho00_h = par[1];
  double Norm = par[2]; // Norm = 0.75 <= sampling | Norm from fitting <= extract rho00

  double dNdCosThetaStar = Norm*((1.0-rho00)+(3.0*rho00-1)*Cos*Cos)*((1.0-rho00_h)+(3.0*rho00_h-1)*CosH*CosH);

  return dNdCosThetaStar;
}

double SpinDensityCosCosHCorr(double *x_val, double *par)
{
  double Cos = x_val[0];
  double CosH = x_val[1];
  double rho00 = par[0];
  double rho00_h = par[1];
  double Norm = par[2]; // Norm = 0.75 <= sampling | Norm from fitting <= extract rho00
  double alpha = par[3];

  double dNdCosThetaStar = Norm*((1.0-rho00)+(3.0*rho00-1)*Cos*Cos)*((1.0-rho00_h)+(3.0*rho00_h-1)*CosH*CosH)*(1.+alpha*Cos*Cos);

  return dNdCosThetaStar;
}

double m2_PID(double *x_val, double *par)
{
  double pT = x_val[0];
  double a = par[0];
  double b = par[1];
  double c = par[2];
  double d = par[3];
  double e = par[4];

  return 1./(a + b*TMath::Exp(-c*(pT-d))) + e;
}

double m2_PID_tanh(double *x_val, double *par)
{
  double pT = x_val[0];
  double a = par[0];
  double b = par[1];
  double c = par[2];
  double d = par[3];
  //double e = par[4];

  return a*TMath::TanH(b*(pT-c))+d;
}

double nsig_PID_tanh(double *x_val, double *par)
{
  double pT = x_val[0];
  double a = par[0];
  double b = par[1];
  double c = par[2];
  //double d = par[3];
  //double e = par[4];

  return a*TMath::TanH(b*(pT-c));
}

double PolyRes(double *x_val, double *par)
{
  double x = x_val[0];
  double y = par[0] + par[1]*(x-1.0/3.0) + 1.0/3.0;

  return y;
}

double flowSample(double *x_val, double *par) 
{
  double x, y, v2;
  x = x_val[0];
  v2 = par[0];
  y = 1.0 + 2.0*v2*cos(2.0*x);

  return y;
}

double flow(double *x_val, double *par) 
{
  double x, y, v2, Norm;
  x = x_val[0];
  v2 = par[0];
  Norm = par[1];
  y = Norm*(1.0 + 2.0*v2*cos(2.0*x));

  return y;
}

double EventPlaneGaus(double *var, double *par)
{
  double Psi2 = var[0];
  double res = par[0]; // res = <cos(2*(Psi2-Psi_RP))>

  double sigma = acos(res)/2.0;
  double sigmaSquare = sigma*sigma;
  double norm = 1.0/sqrt(2.0*sigmaSquare*TMath::Pi());
  double power = -1.0*Psi2*Psi2/(2.0*sigmaSquare);

  double y = norm*exp(power);

  return y;
}

double EventPlaneResolution(double *var, double *par)
{
  double chi = var[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  double y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

double EventPlaneDist(double *var, double *par)
{
  double DeltaPsi = var[0];
  double chi = par[0];
  double arg = chi/TMath::Sqrt(2.0);
  double arg2 = -0.5*chi*chi;
  double pi = TMath::Pi();
  double norm = par[1];

  double cos = TMath::Cos(2.0*DeltaPsi);
  double sin2 = TMath::Sin(2.0*DeltaPsi)*TMath::Sin(2.0*DeltaPsi);
  double y = norm*(TMath::Exp(arg2)+TMath::Sqrt(pi)*arg*cos*TMath::Exp(arg2*sin2)*(1.0+TMath::Erf(arg*cos)));

  return y;
}

double v2_pT_FitFunc(double *x_val, double *par)
{
  // Fit function for v2 vs. pT
  // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
  double v2, pT, a, b, c, d, n;
  pT = x_val[0];
  n  = par[0]; // number-of-constituent quarks
  a  = par[1];
  b  = par[2];
  c  = par[3];
  d  = par[4];

  if(c != 0.0)
  {
    v2 = a*n/(1.0 + exp(-(pT/n - b)/c)) - d*n;
  }
  else v2 = 0.0;

  return v2;
}

double Levy(double *var, double *par)
{
  double const m0 = 1.01940; // phi-meson mass
  double pT   = var[0];
  double mT   = sqrt(pT*pT+m0*m0);
  double dNdy = par[0];
  double n    = par[1];
  double T    = par[2];

  double numer = dNdy*(n-1)*(n-2);
  double denom = n*T*(n*T+m0*(n-2));
  double power = pow(1+(mT-m0)/(n*T),-1.0*n);

  double y = numer*power/denom;

  return y;
}

double pTLevy(double *var, double *par)
{
  double const m0 = 1.01940; // phi-meson mass
  double pT   = var[0];
  double mT   = sqrt(pT*pT+m0*m0);
  double dNdy = par[0];
  double n    = par[1];
  double T    = par[2];

  double numer = dNdy*(n-1)*(n-2);
  double denom = n*T*(n*T+m0*(n-2));
  double power = pow(1+(mT-m0)/(n*T),-1.0*n);

  double y = pT*numer*power/denom;

  return y;
}

double meanLevy(double *var, double *par)
{
  double const m0 = 1.01940; // phi-meson mass
  double pT   = var[0];
  double mT   = sqrt(pT*pT+m0*m0);
  double dNdy = par[0];
  double n    = par[1];
  double T    = par[2];

  double numer = dNdy*(n-1)*(n-2);
  double denom = n*T*(n*T+m0*(n-2));
  double power = pow(1+(mT-m0)/(n*T),-1.0*n);

  double y = pT*pT*numer*power/denom;

  return y;
}

//---------------Back Ground Subtraction-----------------------
double BreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  return BW;
}

double PolyBreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double Poly = par[3] + par[4]*x;

  double y = BW + Poly;

  return y;
}

double Poly(double *x_val, double *par)
{
  double x = x_val[0];
  double y = par[0] + par[1]*x;

  return y;
}

double Poly2BreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double Poly = par[3] + par[4]*x + par[5]*x*x;

  double y = BW + Poly;

  return y;
}

double Poly2(double *x_val, double *par)
{
  double x = x_val[0];
  double y = par[0] + par[1]*x + par[2]*x*x;

  return y;
}

double Poly3MBreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double Poly = par[3] + par[4]*(x-m0) + par[5]*(x-m0)*(x-m0) + par[6]*(x-m0)*(x-m0)*(x-m0);

  double y = BW + Poly;

  return y;
}

double Poly3M(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[4];
  double y = par[0] + par[1]*(x-m0) + par[2]*(x-m0)*(x-m0) + par[3]*(x-m0)*(x-m0)*(x-m0);

  return y;
}

double Poly3BreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double Poly = par[3] + par[4]*x + par[5]*x*x + par[6]*x*x*x;

  double y = BW + Poly;

  return y;
}

double Poly3(double *x_val, double *par)
{
  double x = x_val[0];
  double y = par[0] + par[1]*x + par[2]*x*x + par[3]*x*x*x;

  return y;
}

double Poly3BreitWignerMultiple(double *x_val, double *par)
{
  double x = x_val[0];
  double m01 = par[0];
  double Gamma1 = par[1];
  double Norm1 = par[2];
  double m02 = par[7];
  double Gamma2 = par[8];
  double Norm2 = par[9];
  double m03 = par[10];
  double Gamma3 = par[11];
  double Norm3 = par[12];

  double denom1 = 2.0*TMath::Pi()*((x-m01)*(x-m01)+Gamma1*Gamma1/4.0);
  double BW1 = Norm1*Gamma1/denom1;
  double denom2 = 2.0*TMath::Pi()*((x-m02)*(x-m02)+Gamma2*Gamma2/4.0);
  double BW2 = Norm2*Gamma2/denom2;
  double denom3 = 2.0*TMath::Pi()*((x-m03)*(x-m03)+Gamma3*Gamma3/4.0);
  double BW3 = Norm3*Gamma3/denom3;

  double Poly = par[3] + par[4]*x + par[5]*x*x + par[6]*x*x*x;

  double y = BW1 + BW2 + BW3 + Poly;

  return y;
}

double Poly3MultipleBW(double *x_val, double *par)
{
  double x = x_val[0];
  double m02 = par[4];
  double Gamma2 = par[5];
  double Norm2 = par[6];
  double m03 = par[7];
  double Gamma3 = par[8];
  double Norm3 = par[9];

  double denom2 = 2.0*TMath::Pi()*((x-m02)*(x-m02)+Gamma2*Gamma2/4.0);
  double BW2 = Norm2*Gamma2/denom2;
  double denom3 = 2.0*TMath::Pi()*((x-m03)*(x-m03)+Gamma3*Gamma3/4.0);
  double BW3 = Norm3*Gamma3/denom3;

  double Poly = par[0] + par[1]*x + par[2]*x*x + par[3]*x*x*x;

  double y = BW2 + BW3 + Poly;

  return y;
}

double Poly5BreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double Poly = par[3] + par[4]*x + par[5]*x*x + par[6]*x*x*x + par[7]*x*x*x*x + par[8]*x*x*x*x*x;

  double y = BW + Poly;

  return y;
}

double Poly5(double *x_val, double *par)
{
  double x = x_val[0];
  double y = par[0] + par[1]*x + par[2]*x*x + par[3]*x*x*x + par[4]*x*x*x*x + par[5]*x*x*x*x*x;

  return y;
}

double LegeBreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double p1 = x;
  double p2 = 0.5*(3.0*x*x-1.0);
  double Para = par[3] + par[4]*p1 + par[5]*p2;
  // double Para = par[3]*(1.0 + par[4]*p1 + par[5]*p2);

  double y = BW + Para;

  return y;
}

double Lege(double *x_val, double *par)
{
  double x = x_val[0];
  double p1 = x;
  double p2 = 0.5*(3.0*x*x-1.0);
  double y = par[0] + par[1]*p1 + par[2]*p2;
  // double y = par[0]*(1.0 + par[1]*p1 + par[2]*p2);

  return y;
}

double ChebBreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double t1 = x;
  double t2 = 2.0*x*x-1.0;
  double Cheb = par[3] + par[4]*t1 + par[5]*t2;
  // double Cheb = par[3]*(1.0 + par[4]*t1 + par[5]*t2);

  double y = BW + Cheb;

  return y;
}

double Cheb(double *x_val, double *par)
{
  double x = x_val[0];
  double t1 = x;
  double t2 = 2.0*x*x-1.0;
  double y = par[0] + par[1]*t1 + par[2]*t2;
  // double y = par[0]*(1.0 + par[1]*t1 + par[2]*t2);

  return y;
}

double LogaBreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double log = par[3] + par[4]*TMath::Log(par[5]*x+par[6]);

  double y = BW + log;

  return y;
}

double Loga(double *x_val, double *par)
{
  double x = x_val[0];
  double y = par[0] + par[1]*TMath::Log(par[2]*x+par[3]);

  return y;
}

double SqRtBreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double sqrt = par[3] + par[4]*TMath::Sqrt(x+par[5]);

  double y = BW + sqrt;

  return y;
}

double SqRt(double *x_val, double *par)
{
  double x = x_val[0];
  double y = par[0] + par[1]*TMath::Sqrt(x+par[2]);

  return y;
}

double DELPHIBreitWigner(double *var, double *par)
{
  double x = var[0];
  double m0 = par[0];
  double Gamma = par[1];
  double nsig = par[2];
  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double sig = nsig*Gamma/denom;

  double mMassKaon = 0.49368;
  double c0 = par[3];
  double c1 = par[4];
  double c2 = par[5];
  double c3 = par[6];
  double c4 = par[7];
  double nbkg = par[8];
  double x_var = x-mMassKaon*2.0;
  double bkg = nbkg*TMath::Power(x_var,c0)*TMath::Exp(c1*x+c2*x*x+c3*x*x*x+c4*x*x*x*x);

  return sig+bkg;
}

double DELPHI(double *var, double *par)
{
  double mMassKaon = 0.49368;
  double x = var[0];
  double c0 = par[0];
  double c1 = par[1];
  double c2 = par[2];
  double c3 = par[3];
  double c4 = par[4];
  double nbkg = par[5];
  double x_var = x-mMassKaon*2.0;
  double bkg = nbkg*TMath::Power(x_var,c0)*TMath::Exp(c1*x+c2*x*x+c3*x*x*x+c4*x*x*x*x);

  return bkg;
}

//---------------Back Ground Subtraction-----------------------

float ErrorAdd(float x, float y)
{
  return sqrt(x*x+y*y);
}

float ErrTimes(float x, float y, float dx, float dy)
{
  return x*y*ErrorAdd(dx/x,dy/y);
}

float ErrDiv(float x, float y, float dx, float dy)
{
  return x/y*ErrorAdd(dx/x,dy/y);
}


// full event plane resolution
double Resolution_Full(double *x_val, double *par)
{
  double y;
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

// tof matching efficiency
double tof_Kaon(double* x, double* par)
{
  return par[0]*(1.0 / (pow(x[0] - par[1], 2) + par[2]) - par[4] / (exp(x[0] - par[3]) + par[5]) + par[6]);
}


double rhoGrhoH(double *x, double*par)
{
  double v2 = x[0];
  double pT = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double rho = par[2];

  return 1./4.*(2.-2.*rho+pT*pT*(-1.-3.*rho*(-1.+v2)+v2)/p2)-1./3.;
}

double rhoGfromHv2(double *x, double*par)
{
  double v2 = x[0];
  double pT = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return 1./8.*(3. + 2.*rerho1n1 - rho00+v2 + 6.*rerho1n1*v2 - 3.*rho00*v2 
                - (-1. + 2.*rerho1n1 + 3.*rho00)*(-1. + v2)*cos2theta 
                + 2.*TMath::Sqrt(2)*real*(-1. + v2)*sin2theta) - 1./3.;
}

double realGfromHv2(double *x, double*par)
{
  double v2 = x[0];
  double pT = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return 0;
}

double imagGfromHv2(double *x, double*par)
{
  double v2 = x[0];
  double pT = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return v2*(imag*costheta + TMath::Sqrt(2.)*imrho1n1*sintheta);
}

double rerho1n1GfromHv2(double *x, double*par)
{
  double v2 = x[0];
  double pT = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return 1./16.*(- ((1.+6.*rerho1n1-3.*rho00)*(-1.+v2)) 
                 + (-1.+2.*rerho1n1+3.*rho00)*(3.+v2)*cos2theta 
                 - 2.*TMath::Sqrt(2.)*real*(3.+v2)*sin2theta);
}

double imrho1n1GfromHv2(double *x, double*par)
{
  double v2 = x[0];
  double pT = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return 0;
}

double rhoGfromHphipsi(double *x, double*par)
{
  double ppsi = x[0];
  double pT = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  double cosppsi = TMath::Cos(ppsi);
  double sinppsi = TMath::Sin(ppsi);
  double cos2ppsi = TMath::Cos(2*ppsi);
  double sin2ppsi = TMath::Sin(2*ppsi);

  return 1./8.*(3.+2.*rerho1n1-rho00+(1.+6.*rerho1n1-3.*rho00)*cos2ppsi
                +4.*sin2ppsi*(TMath::Sqrt(2.)*imag*costheta+2.*imrho1n1*sintheta)
                +2.*sinppsi*sinppsi*((-1.+2.*rerho1n1+3.*rho00)*cos2theta-2.*TMath::Sqrt(2.)*real*sin2theta))-1./3.;
}

double realGfromHphipsi(double *x, double*par)
{
  double ppsi = x[0];
  double pT = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  double cosppsi = TMath::Cos(ppsi);
  double sinppsi = TMath::Sin(ppsi);

  return cosppsi*(TMath::Sqrt(2.)*imrho1n1*costheta-imag*sintheta)
         +1./4.*sinppsi*(-4.*real*cos2theta+TMath::Sqrt(2.)*(1.-2.*rerho1n1-3.*rho00)*sin2theta);
}

double imagGfromHphipsi(double *x, double*par)
{
  double ppsi = x[0];
  double pT = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  double cosppsi = TMath::Cos(ppsi);
  double sinppsi = TMath::Sin(ppsi);
  double sin2ppsi = TMath::Sin(2*ppsi);

  return 1./4.*(-4.0*imag*costheta*sinppsi*sinppsi
                +cosppsi*cosppsi*(4.*imag*costheta+4.*TMath::Sqrt(2.)*imrho1n1*sintheta)
                -TMath::Sqrt(2.)*sintheta*(4.*imrho1n1*sinppsi*sinppsi+rerho1n1*sin2ppsi*sintheta)
                +cosppsi*sinppsi*(TMath::Sqrt(2.)*(-1.-5.*rerho1n1+3.*rho00)+TMath::Sqrt(2)*(-1.+rerho1n1+3.*rho00)*cos2theta-4.*real*sin2theta));
}

double rerho1n1GfromHphipsi(double *x, double*par)
{
  double ppsi = x[0];
  double pT = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  double cosppsi = TMath::Cos(ppsi);
  double sinppsi = TMath::Sin(ppsi);

  return 1./16.*(-((-1.+2.*rerho1n1+3.*rho00)*costheta*costheta*(-3.+sinppsi*sinppsi))
                 +(1.+6.*rerho1n1-3.*rho00)*(1.+sinppsi*sinppsi)
                 +4.*TMath::Sqrt(2)*real*costheta*(-3.+sinppsi*sinppsi)*sintheta
                 +(-1.+2.*rerho1n1+3.*rho00)*(-3.+sinppsi*sinppsi)*sintheta*sintheta
                 +8.*cosppsi*sinppsi*(-TMath::Sqrt(2.)*imag*costheta-2.*imrho1n1*sintheta)
                 +cosppsi*cosppsi*(-1.-6.*rerho1n1+3.*rho00+(-1.+2.*rerho1n1+3.*rho00)*costheta*costheta-4.*TMath::Sqrt(2.)*real*costheta*sintheta+(1.-2.*rerho1n1-3.*rho00)*sintheta*sintheta));
}

double imrho1n1GfromHphipsi(double *x, double*par)
{
  double ppsi = x[0];
  double pT = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  double cosppsi = TMath::Cos(ppsi);
  double sinppsi = TMath::Sin(ppsi);

  return 1./4.*( 2.*sinppsi*(2.*imrho1n1*costheta-TMath::Sqrt(2.)*imag*sintheta)
                +cosppsi*(2.*TMath::Sqrt(2.)*real*cos2theta+(-1.+2.*rerho1n1+3.*rho00)*sin2theta));
}

double rhoGfromHpt(double *x, double*par)
{
  double pT = x[0];
  double v2 = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return 1./8.*(3. + 2.*rerho1n1 - rho00+v2 + 6.*rerho1n1*v2 - 3.*rho00*v2 
		- (-1. + 2.*rerho1n1 + 3.*rho00)*(-1. + v2)*cos2theta 
		+ 2.*TMath::Sqrt(2)*real*(-1. + v2)*sin2theta) - 1./3.;
}

double realGfromHpt(double *x, double*par)
{
  double pT = x[0];
  double v2 = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return 0;
}

double imagGfromHpt(double *x, double*par)
{
  double pT = x[0];
  double v2 = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return v2*(imag*costheta + TMath::Sqrt(2.)*imrho1n1*sintheta);
}

double rerho1n1GfromHpt(double *x, double*par)
{
  double pT = x[0];
  double v2 = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return 1./16.*(- ((1.+6.*rerho1n1-3.*rho00)*(-1.+v2)) 
		 + (-1.+2.*rerho1n1+3.*rho00)*(3.+v2)*cos2theta 
		 - 2.*TMath::Sqrt(2.)*real*(3.+v2)*sin2theta);
}

double imrho1n1GfromHpt(double *x, double*par)
{
  double pT = x[0];
  double v2 = par[0];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(par[1]);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return 0;
}


double rhoGfromHy(double *x, double*par)
{
  double y = x[0];
  double v2 = par[0];
  double pT = par[1];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(y);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return 1./8.*(3. + 2.*rerho1n1 - rho00+v2 + 6.*rerho1n1*v2 - 3.*rho00*v2 
		- (-1. + 2.*rerho1n1 + 3.*rho00)*(-1. + v2)*cos2theta 
		+ 2.*TMath::Sqrt(2)*real*(-1. + v2)*sin2theta) - 1./3.;
}

double realGfromHy(double *x, double*par)
{
  double y = x[0];
  double v2 = par[0];
  double pT = par[1];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(y);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return 0;
}

double imagGfromHy(double *x, double*par)
{
  double y = x[0];
  double v2 = par[0];
  double pT = par[1];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(y);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return v2*(imag*costheta + TMath::Sqrt(2.)*imrho1n1*sintheta);
}

double rerho1n1GfromHy(double *x, double*par)
{
  double y = x[0];
  double v2 = par[0];
  double pT = par[1];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(y);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return 1./16.*(- ((1.+6.*rerho1n1-3.*rho00)*(-1.+v2)) 
		 + (-1.+2.*rerho1n1+3.*rho00)*(3.+v2)*cos2theta 
		 - 2.*TMath::Sqrt(2.)*real*(3.+v2)*sin2theta);
}

double imrho1n1GfromHy(double *x, double*par)
{
  double y = x[0];
  double v2 = par[0];
  double pT = par[1];
  double mT = TMath::Sqrt(1.019461*1.019461+pT*pT);
  double pz = mT*TMath::SinH(y);
  double p2 = pT*pT+pz*pz;
  double p  = TMath::Sqrt(p2);
  double rho00 = par[2];
  double real = par[3];
  double imag = par[4];
  double rerho1n1 = par[5];
  double imrho1n1 = par[6];

  double costheta = pT/p;
  double cos2theta = 2.*(pT*pT/p2)-1.;
  double sintheta = -pz/p;
  double sin2theta = -2.*pz*pT/p2;

  return 0;
}

double acceptance2D(double *x, double *par)
{
  double N = par[0];
  double F = par[1];
  double G = par[2];
  double H = par[3];
  double I = par[4];
  double cos = x[0];
  double sinbeta = x[1];

  return     N + F*TMath::Power((1.-cos*cos),1./2.)*TMath::Power(sinbeta,1.)
               + G*TMath::Power((1.-cos*cos),2./2.)*TMath::Power(sinbeta,2.) 
               + H*TMath::Power((1.-cos*cos),3./2.)*TMath::Power(sinbeta,3.)
               + I*TMath::Power((1.-cos*cos),4./2.)*TMath::Power(sinbeta,4.);
}

double acceptance2Dsimple(double *x, double *par)
{
  double N = par[0];
  double F = par[1];
  double G = par[2];
  double H = par[3];
  double I = par[4];
  double J = par[5];
  double K = par[6];
  //double L = par[7];
  //double M = par[8];
  //double N = par[9];
  //double O = par[10];
  //double P = par[11];
  //double Q = par[12];
  double cos = x[0];
  //double cos1beta = TMath::Cos(2.*x[1]);
  //double cos2beta = TMath::Cos(4.*x[1]);
  double sin1beta = TMath::Sin(x[1]);
  double sin2beta = TMath::Sin(2.*x[1]);
  //double sin3beta = TMath::Sin(2.*x[1]);
  //double sin4beta = TMath::Sin(4.*x[1]);

  //return N*(1. + F*TMath::Power(cos,2.) + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
  return N*(1. + F*TMath::Power(cos,2.) + G*TMath::Power(cos,4.) + H*TMath::Power(cos,6.)+ I*TMath::Power(cos,8.))*(1. + J*TMath::Power(sin1beta,2.) + K*TMath::Power(sin1beta,4.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
  //return A + F*TMath::Power(cos,2.)*TMath::Power(sin1beta,2.) 
  //         + G*TMath::Power(cos,2.)*TMath::Power(sin1beta,4.)
  //         + H*TMath::Power(cos,2.)*TMath::Power(sin1beta,6.)
  //         //+ I*TMath::Power(cos,2.)*TMath::Power(sin2beta,2.) 
  //         //+ J*TMath::Power(cos,2.)*TMath::Power(sin2beta,4.)
  //         //+ K*TMath::Power(cos,2.)*TMath::Power(sin2beta,6.)
  //         + I*TMath::Power(cos,4.)*TMath::Power(sin1beta,2.) 
  //         + J*TMath::Power(cos,4.)*TMath::Power(sin1beta,4.)
  //         + K*TMath::Power(cos,4.)*TMath::Power(sin1beta,6.);
  //         //+ O*TMath::Power(cos,4.)*TMath::Power(sin2beta,2.) 
  //         //+ P*TMath::Power(cos,4.)*TMath::Power(sin2beta,4.)
  //         //+ Q*TMath::Power(cos,4.)*TMath::Power(sin2beta,6.);

  //         + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
}
double acceptance2Dsimple2(double *x, double *par)
{
  double A = par[0];
  double F = par[1];
  double G = par[2];
  double H = par[3];
  double I = par[4];
  double J = par[5];
  double K = par[6];
  double L = par[7];
  double M = par[8];
  double N = par[9];
  //double O = par[10];
  //double P = par[11];
  //double Q = par[12];
  double cos = x[0];
  double cos1beta = TMath::Cos(x[1]);
  double cos2beta = TMath::Cos(2.*x[1]);
  double sin1beta = TMath::Sin(x[1]);
  double sin2beta = TMath::Sin(2.*x[1]);
  //double sin3beta = TMath::Sin(2.*x[1]);
  //double sin4beta = TMath::Sin(4.*x[1]);

  //return N*(1. + F*TMath::Power(cos,2.) + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
  return A  + F*(1.+ G*TMath::Power(sin1beta,2.) + H*TMath::Power(cos1beta,2.)) 
            + I*(1.+ J*TMath::Power(sin1beta,2.) + K*TMath::Power(cos1beta,2.))*TMath::Power(cos,2.)
            + L*(1.+ M*TMath::Power(sin1beta,2.) + N*TMath::Power(cos1beta,2.))*TMath::Power(cos,4.);
            //+ J*(1.+ K*TMath::Power(sin1beta,2.))*TMath::Power(cos,4.))*TMath::Power(sin1beta,2.);
           //+ I*TMath::Power(cos,2.)*TMath::Power(sin2beta,2.) 
           //+ J*TMath::Power(cos,2.)*TMath::Power(sin2beta,4.)
           //+ K*TMath::Power(cos,2.)*TMath::Power(sin2beta,6.)
           //+ I*TMath::Power(cos,4.)*TMath::Power(sin1beta,2.) 
           //+ J*TMath::Power(cos,4.)*TMath::Power(sin1beta,4.)
           //+ K*TMath::Power(cos,4.)*TMath::Power(sin1beta,6.);
           //+ O*TMath::Power(cos,4.)*TMath::Power(sin2beta,2.) 
           //+ P*TMath::Power(cos,4.)*TMath::Power(sin2beta,4.)
           //+ Q*TMath::Power(cos,4.)*TMath::Power(sin2beta,6.);

  //         + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
}

double acceptance2Dsimilar(double *x, double *p)
{
  return p[0]  +p[1]*x[0]*x[0]
               +p[2]*x[0]*x[0]*x[0]*x[0]
               +p[3]*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*cos(x[1])
               +p[4]*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*cos(2.*x[1])
               +p[5]*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*sin(x[1])
               +p[6]*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*sin(2.*x[1])
               +p[7]*(1-x[0]*x[0])*cos(x[1])
               +p[8]*(1-x[0]*x[0])*cos(2.*x[1])
               +p[9]*(1-x[0]*x[0])*sin(x[1])
               +p[10]*(1-x[0]*x[0])*sin(2.*x[1]);
}

double acceptance2Dsimple3(double *x, double *par)
{
  double F = par[0];
  //double F = par[0];
  double G = par[1];
  double H = par[2];
  double I = par[3];
  double J = par[4];
  double K = par[5];
  double L = par[6];
  double M = par[7];
  double N = par[8];
  double O = par[9];
  double P = par[10];
  double Q = par[11];
  double R = par[12];
  double S = par[13];
  double T = par[14];
  double cos = x[0];
  double cos1beta = TMath::Cos(x[1]);
  double cos2beta = TMath::Cos(2.*x[1]);
  double sin1beta = TMath::Sin(x[1]);
  double sin2beta = TMath::Sin(2.*x[1]);
  //double sin3beta = TMath::Sin(2.*x[1]);
  //double sin4beta = TMath::Sin(4.*x[1]);

  //return N*(1. + F*TMath::Power(cos,2.) + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
  return       F + G*TMath::Power(sin1beta,2.) + H*TMath::Power(cos1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(cos1beta,4.)
            + (K + L*TMath::Power(sin1beta,2.) + M*TMath::Power(cos1beta,2.) + N*TMath::Power(sin1beta,4.) + O*TMath::Power(cos1beta,4.))*TMath::Power(cos,2.)
            + (P + Q*TMath::Power(sin1beta,2.) + R*TMath::Power(cos1beta,2.) + S*TMath::Power(sin1beta,4.) + T*TMath::Power(cos1beta,4.))*TMath::Power(cos,4.);
            //+ J*(1.+ K*TMath::Power(sin1beta,2.))*TMath::Power(cos,4.))*TMath::Power(sin1beta,2.);
           //+ I*TMath::Power(cos,2.)*TMath::Power(sin2beta,2.) 
           //+ J*TMath::Power(cos,2.)*TMath::Power(sin2beta,4.)
           //+ K*TMath::Power(cos,2.)*TMath::Power(sin2beta,6.)
           //+ I*TMath::Power(cos,4.)*TMath::Power(sin1beta,2.) 
           //+ J*TMath::Power(cos,4.)*TMath::Power(sin1beta,4.)
           //+ K*TMath::Power(cos,4.)*TMath::Power(sin1beta,6.);
           //+ O*TMath::Power(cos,4.)*TMath::Power(sin2beta,2.) 
           //+ P*TMath::Power(cos,4.)*TMath::Power(sin2beta,4.)
           //+ Q*TMath::Power(cos,4.)*TMath::Power(sin2beta,6.);

  //         + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
}
double acceptance2Dsimple4(double *x, double *par)
{
  double F = par[0];
  //double F = par[0];
  double G = par[1];
  double H = par[2];
  double I = par[3];
  double J = par[4];
  double K = par[5];
  double L = par[6];
  double M = par[7];
  double N = par[8];
  double cos = x[0];
  double cos1beta = TMath::Cos(x[1]);
  double cos2beta = TMath::Cos(2.*x[1]);
  double sin1beta = TMath::Sin(x[1]);
  double sin2beta = TMath::Sin(2.*x[1]);
  //double sin3beta = TMath::Sin(2.*x[1]);
  //double sin4beta = TMath::Sin(4.*x[1]);

  //return N*(1. + F*TMath::Power(cos,2.) + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
  return    F*(1.+ G*TMath::Power(sin1beta,2.) + H*TMath::Power(sin1beta,4.)
            + (I + J*TMath::Power(sin1beta,2.) + K*TMath::Power(sin1beta,4.))*TMath::Power(cos,2.)
            + (L + M*TMath::Power(sin1beta,2.) + N*TMath::Power(sin1beta,4.))*TMath::Power(cos,4.));
            //+ J*(1.+ K*TMath::Power(sin1beta,2.))*TMath::Power(cos,4.))*TMath::Power(sin1beta,2.);
           //+ I*TMath::Power(cos,2.)*TMath::Power(sin2beta,2.) 
           //+ J*TMath::Power(cos,2.)*TMath::Power(sin2beta,4.)
           //+ K*TMath::Power(cos,2.)*TMath::Power(sin2beta,6.)
           //+ I*TMath::Power(cos,4.)*TMath::Power(sin1beta,2.) 
           //+ J*TMath::Power(cos,4.)*TMath::Power(sin1beta,4.)
           //+ K*TMath::Power(cos,4.)*TMath::Power(sin1beta,6.);
           //+ O*TMath::Power(cos,4.)*TMath::Power(sin2beta,2.) 
           //+ P*TMath::Power(cos,4.)*TMath::Power(sin2beta,4.)
           //+ Q*TMath::Power(cos,4.)*TMath::Power(sin2beta,6.);

  //         + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
}
double acceptance2Dsimple5(double *x, double *par)
{
  double F = par[0];
  //double F = par[0];
  double G = par[1];
  double H = par[2];
  double I = par[3];
  double J = par[4];
  double K = par[5];
  double L = par[6];
  double cos = x[0];
  double cos1beta = TMath::Cos(x[1]);
  double cos2beta = TMath::Cos(2.*x[1]);
  double sin1beta = TMath::Sin(x[1]);
  double sin2beta = TMath::Sin(2.*x[1]);
  //double sin3beta = TMath::Sin(2.*x[1]);
  //double sin4beta = TMath::Sin(4.*x[1]);

  //return N*(1. + F*TMath::Power(cos,2.) + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
  return       F 
            + (G + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.))*TMath::Power(cos,2.)
            + (J + K*TMath::Power(sin1beta,2.) + L*TMath::Power(sin1beta,4.))*TMath::Power(cos,4.);
            //+ J*(1.+ K*TMath::Power(sin1beta,2.))*TMath::Power(cos,4.))*TMath::Power(sin1beta,2.);
           //+ I*TMath::Power(cos,2.)*TMath::Power(sin2beta,2.) 
           //+ J*TMath::Power(cos,2.)*TMath::Power(sin2beta,4.)
           //+ K*TMath::Power(cos,2.)*TMath::Power(sin2beta,6.)
           //+ I*TMath::Power(cos,4.)*TMath::Power(sin1beta,2.) 
           //+ J*TMath::Power(cos,4.)*TMath::Power(sin1beta,4.)
           //+ K*TMath::Power(cos,4.)*TMath::Power(sin1beta,6.);
           //+ O*TMath::Power(cos,4.)*TMath::Power(sin2beta,2.) 
           //+ P*TMath::Power(cos,4.)*TMath::Power(sin2beta,4.)
           //+ Q*TMath::Power(cos,4.)*TMath::Power(sin2beta,6.);

  //         + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
}

double acceptance2Dsimple6(double *x, double *par)
{
  double F = par[0];
  //double F = par[0];
  double G = par[1];
  double H = par[2];
  double I = par[3];
  double J = par[4];
  double K = par[5];
  double L = par[6];
  double cos = x[0];
  double cos1beta = TMath::Cos(x[1]);
  double cos2beta = TMath::Cos(2.*x[1]);
  double sin1beta = TMath::Sin(x[1]);
  double sin2beta = TMath::Sin(2.*x[1]);
  //double sin3beta = TMath::Sin(2.*x[1]);
  //double sin4beta = TMath::Sin(4.*x[1]);

  //return N*(1. + F*TMath::Power(cos,2.) + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
  return       F + G*TMath::Power(sin1beta,2.) + H*TMath::Power(sin1beta,4.)
            + (I*TMath::Power(sin1beta,2.) + J*TMath::Power(sin1beta,4.))*TMath::Power(cos,2.)
            + (K*TMath::Power(sin1beta,2.) + L*TMath::Power(sin1beta,4.))*TMath::Power(cos,4.);
            //+ J*(1.+ K*TMath::Power(sin1beta,2.))*TMath::Power(cos,4.))*TMath::Power(sin1beta,2.);
           //+ I*TMath::Power(cos,2.)*TMath::Power(sin2beta,2.) 
           //+ J*TMath::Power(cos,2.)*TMath::Power(sin2beta,4.)
           //+ K*TMath::Power(cos,2.)*TMath::Power(sin2beta,6.)
           //+ I*TMath::Power(cos,4.)*TMath::Power(sin1beta,2.) 
           //+ J*TMath::Power(cos,4.)*TMath::Power(sin1beta,4.)
           //+ K*TMath::Power(cos,4.)*TMath::Power(sin1beta,6.);
           //+ O*TMath::Power(cos,4.)*TMath::Power(sin2beta,2.) 
           //+ P*TMath::Power(cos,4.)*TMath::Power(sin2beta,4.)
           //+ Q*TMath::Power(cos,4.)*TMath::Power(sin2beta,6.);

  //         + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
}
double acceptance2Dsimple7(double *x, double *par)
{
  double F = par[0];
  //double F = par[0];
  double G = par[1];
  double H = par[2];
  double I = par[3];
  double J = par[4];
  double K = par[5];
  double L = par[6];
  double M = par[7];
  double N = par[8];
  double O = par[9];
  double P = par[10];
  double Q = par[11];
  double cos = x[0];
  double cos1beta = TMath::Cos(x[1]);
  double cos2beta = TMath::Cos(2.*x[1]);
  double sin1beta = TMath::Sin(x[1]);
  double sin2beta = TMath::Sin(2.*x[1]);
  //double sin3beta = TMath::Sin(2.*x[1]);
  //double sin4beta = TMath::Sin(4.*x[1]);

  //return N*(1. + F*TMath::Power(cos,2.) + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
  return   F* (1.+ G*TMath::Power(sin1beta,2.) + H*TMath::Power(sin1beta,4.)
            + (I + J*TMath::Power(sin1beta,2.) + K*TMath::Power(sin1beta,4.))*TMath::Power(cos,2.)
            + (L + M*TMath::Power(sin1beta,2.) + N*TMath::Power(sin1beta,4.))*TMath::Power(cos,4.)
            + (O + P*TMath::Power(sin1beta,2.) + Q*TMath::Power(sin1beta,4.))*TMath::Power(cos,6.));
            //+ J*(1.+ K*TMath::Power(sin1beta,2.))*TMath::Power(cos,4.))*TMath::Power(sin1beta,2.);
           //+ I*TMath::Power(cos,2.)*TMath::Power(sin2beta,2.) 
           //+ J*TMath::Power(cos,2.)*TMath::Power(sin2beta,4.)
           //+ K*TMath::Power(cos,2.)*TMath::Power(sin2beta,6.)
           //+ I*TMath::Power(cos,4.)*TMath::Power(sin1beta,2.) 
           //+ J*TMath::Power(cos,4.)*TMath::Power(sin1beta,4.)
           //+ K*TMath::Power(cos,4.)*TMath::Power(sin1beta,6.);
           //+ O*TMath::Power(cos,4.)*TMath::Power(sin2beta,2.) 
           //+ P*TMath::Power(cos,4.)*TMath::Power(sin2beta,4.)
           //+ Q*TMath::Power(cos,4.)*TMath::Power(sin2beta,6.);

  //         + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
}
double acceptance2Dsimple8(double *x, double *par)
{
  double F = par[0];
  //double F = par[0];
  double G = par[1];
  double H = par[2];
  double I = par[3];
  double J = par[4];
  double K = par[5];
  double L = par[6];
  double M = par[7];
  double N = par[8];
  double O = par[9];
  double P = par[10];
  double Q = par[11];
  double R = par[12];
  double S = par[13];
  double T = par[14];
  double U = par[15];
  double cos = x[0];
  double cos1beta = TMath::Cos(x[1]);
  double cos2beta = TMath::Cos(2.*x[1]);
  double sin1beta = TMath::Sin(x[1]);
  double sin2beta = TMath::Sin(2.*x[1]);
  //double sin3beta = TMath::Sin(2.*x[1]);
  //double sin4beta = TMath::Sin(4.*x[1]);

  //return N*(1. + F*TMath::Power(cos,2.) + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
  return   F* (1.+ G*TMath::Power(sin1beta,2.) + H*TMath::Power(sin1beta,4.) + R*TMath::Power(sin1beta,6)
            + (I + J*TMath::Power(sin1beta,2.) + K*TMath::Power(sin1beta,4.) + S*TMath::Power(sin1beta,6))*TMath::Power(cos,2.)
            + (L + M*TMath::Power(sin1beta,2.) + N*TMath::Power(sin1beta,4.) + T*TMath::Power(sin1beta,6))*TMath::Power(cos,4.)
            + (O + P*TMath::Power(sin1beta,2.) + Q*TMath::Power(sin1beta,4.) + U*TMath::Power(sin1beta,6))*TMath::Power(cos,6.));
            //+ J*(1.+ K*TMath::Power(sin1beta,2.))*TMath::Power(cos,4.))*TMath::Power(sin1beta,2.);
           //+ I*TMath::Power(cos,2.)*TMath::Power(sin2beta,2.) 
           //+ J*TMath::Power(cos,2.)*TMath::Power(sin2beta,4.)
           //+ K*TMath::Power(cos,2.)*TMath::Power(sin2beta,6.)
           //+ I*TMath::Power(cos,4.)*TMath::Power(sin1beta,2.) 
           //+ J*TMath::Power(cos,4.)*TMath::Power(sin1beta,4.)
           //+ K*TMath::Power(cos,4.)*TMath::Power(sin1beta,6.);
           //+ O*TMath::Power(cos,4.)*TMath::Power(sin2beta,2.) 
           //+ P*TMath::Power(cos,4.)*TMath::Power(sin2beta,4.)
           //+ Q*TMath::Power(cos,4.)*TMath::Power(sin2beta,6.);

  //         + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
}
double SpinDensity2DcosAcc(double *x, double *par)
{
  //double F = par[0];
  double G = par[6];
  double H = par[7];
  double I = par[8];
  double J = par[9];
  double K = par[10];
  double L = par[11];
  double M = par[12];
  double N = par[13];
  double cos = x[0];
  double cos1beta = TMath::Cos(x[1]);
  double cos2beta = TMath::Cos(2.*x[1]);
  double sin1beta = TMath::Sin(x[1]);
  double sin2beta = TMath::Sin(2.*x[1]);
  //double sin3beta = TMath::Sin(2.*x[1]);
  //double sin4beta = TMath::Sin(4.*x[1]);
  double rhofunc =  par[5]*((1.-par[0])+(3.*par[0] -1.)*x[0]*x[0]
               -sqrt(2)*par[1]*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*cos1beta
               +sqrt(2)*par[2]*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*sin1beta
               -2.*par[3]*(1-x[0]*x[0])*cos2beta
               +2.*par[4]*(1-x[0]*x[0])*sin2beta);

  //return N*(1. + F*TMath::Power(cos,2.) + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
  return rhofunc * (1.+ G*TMath::Power(sin1beta,2.) + H*TMath::Power(sin1beta,4.)
                 + (I + J*TMath::Power(sin1beta,2.) + K*TMath::Power(sin1beta,4.))*TMath::Power(cos,2.)
                 + (L + M*TMath::Power(sin1beta,2.) + N*TMath::Power(sin1beta,4.))*TMath::Power(cos,4.));
            //+ J*(1.+ K*TMath::Power(sin1beta,2.))*TMath::Power(cos,4.))*TMath::Power(sin1beta,2.);
           //+ I*TMath::Power(cos,2.)*TMath::Power(sin2beta,2.) 
           //+ J*TMath::Power(cos,2.)*TMath::Power(sin2beta,4.)
           //+ K*TMath::Power(cos,2.)*TMath::Power(sin2beta,6.)
           //+ I*TMath::Power(cos,4.)*TMath::Power(sin1beta,2.) 
           //+ J*TMath::Power(cos,4.)*TMath::Power(sin1beta,4.)
           //+ K*TMath::Power(cos,4.)*TMath::Power(sin1beta,6.);
           //+ O*TMath::Power(cos,4.)*TMath::Power(sin2beta,2.) 
           //+ P*TMath::Power(cos,4.)*TMath::Power(sin2beta,4.)
           //+ Q*TMath::Power(cos,4.)*TMath::Power(sin2beta,6.);

  //         + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
}
double SpinDensity2DcosAccLonger(double *x, double *par)
{
  //double F = par[0];
  double G = par[6];
  double H = par[7];
  double I = par[8];
  double J = par[9];
  double K = par[10];
  double L = par[11];
  double M = par[12];
  double N = par[13];
  double O = par[14];
  double P = par[15];
  double Q = par[16];
  double cos = x[0];
  double cos1beta = TMath::Cos(x[1]);
  double cos2beta = TMath::Cos(2.*x[1]);
  double sin1beta = TMath::Sin(x[1]);
  double sin2beta = TMath::Sin(2.*x[1]);
  //double sin3beta = TMath::Sin(2.*x[1]);
  //double sin4beta = TMath::Sin(4.*x[1]);
  double rhofunc =  par[5]*((1.-par[0])+(3.*par[0] -1.)*x[0]*x[0]
               -sqrt(2)*par[1]*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*cos1beta
               +sqrt(2)*par[2]*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*sin1beta
               -2.*par[3]*(1-x[0]*x[0])*cos2beta
               +2.*par[4]*(1-x[0]*x[0])*sin2beta);

  //return N*(1. + F*TMath::Power(cos,2.) + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
  return rhofunc * (1.+ G*TMath::Power(sin1beta,2.) + H*TMath::Power(sin1beta,4.)
                 + (I + J*TMath::Power(sin1beta,2.) + K*TMath::Power(sin1beta,4.))*TMath::Power(cos,2.)
                 + (L + M*TMath::Power(sin1beta,2.) + N*TMath::Power(sin1beta,4.))*TMath::Power(cos,4.)
                 + (O + P*TMath::Power(sin1beta,2.) + Q*TMath::Power(sin1beta,4.))*TMath::Power(cos,6.));
            //+ J*(1.+ K*TMath::Power(sin1beta,2.))*TMath::Power(cos,4.))*TMath::Power(sin1beta,2.);
           //+ I*TMath::Power(cos,2.)*TMath::Power(sin2beta,2.) 
           //+ J*TMath::Power(cos,2.)*TMath::Power(sin2beta,4.)
           //+ K*TMath::Power(cos,2.)*TMath::Power(sin2beta,6.)
           //+ I*TMath::Power(cos,4.)*TMath::Power(sin1beta,2.) 
           //+ J*TMath::Power(cos,4.)*TMath::Power(sin1beta,4.)
           //+ K*TMath::Power(cos,4.)*TMath::Power(sin1beta,6.);
           //+ O*TMath::Power(cos,4.)*TMath::Power(sin2beta,2.) 
           //+ P*TMath::Power(cos,4.)*TMath::Power(sin2beta,4.)
           //+ Q*TMath::Power(cos,4.)*TMath::Power(sin2beta,6.);

  //         + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
}
double SpinDensity2DcosAccLongerer(double *x, double *par)
{
  //double F = par[0];
  double G = par[6];
  double H = par[7];
  double I = par[8];
  double J = par[9];
  double K = par[10];
  double L = par[11];
  double M = par[12];
  double N = par[13];
  double O = par[14];
  double P = par[15];
  double Q = par[16];
  double R = par[17];
  double S = par[18];
  double T = par[19];
  double U = par[20];
  double cos = x[0];
  double cos1beta = TMath::Cos(x[1]);
  double cos2beta = TMath::Cos(2.*x[1]);
  double sin1beta = TMath::Sin(x[1]);
  double sin2beta = TMath::Sin(2.*x[1]);
  //double sin3beta = TMath::Sin(2.*x[1]);
  //double sin4beta = TMath::Sin(4.*x[1]);
  double rhofunc =  par[5]*((1.-par[0])+(3.*par[0] -1.)*x[0]*x[0]
               -sqrt(2)*par[1]*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*cos1beta
               +sqrt(2)*par[2]*2*x[0]*TMath::Sqrt(1-x[0]*x[0])*sin1beta
               -2.*par[3]*(1-x[0]*x[0])*cos2beta
               +2.*par[4]*(1-x[0]*x[0])*sin2beta);

  //return N*(1. + F*TMath::Power(cos,2.) + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
  return rhofunc * (1.+ G*TMath::Power(sin1beta,2.) + H*TMath::Power(sin1beta,4.) + R*TMath::Power(sin1beta,6) 
                 + (I + J*TMath::Power(sin1beta,2.) + K*TMath::Power(sin1beta,4.) + S*TMath::Power(sin1beta,6))*TMath::Power(cos,2.)
                 + (L + M*TMath::Power(sin1beta,2.) + N*TMath::Power(sin1beta,4.) + T*TMath::Power(sin1beta,6))*TMath::Power(cos,4.)
                 + (O + P*TMath::Power(sin1beta,2.) + Q*TMath::Power(sin1beta,4.) + U*TMath::Power(sin1beta,6))*TMath::Power(cos,6.));
            //+ J*(1.+ K*TMath::Power(sin1beta,2.))*TMath::Power(cos,4.))*TMath::Power(sin1beta,2.);
           //+ I*TMath::Power(cos,2.)*TMath::Power(sin2beta,2.) 
           //+ J*TMath::Power(cos,2.)*TMath::Power(sin2beta,4.)
           //+ K*TMath::Power(cos,2.)*TMath::Power(sin2beta,6.)
           //+ I*TMath::Power(cos,4.)*TMath::Power(sin1beta,2.) 
           //+ J*TMath::Power(cos,4.)*TMath::Power(sin1beta,4.)
           //+ K*TMath::Power(cos,4.)*TMath::Power(sin1beta,6.);
           //+ O*TMath::Power(cos,4.)*TMath::Power(sin2beta,2.) 
           //+ P*TMath::Power(cos,4.)*TMath::Power(sin2beta,4.)
           //+ Q*TMath::Power(cos,4.)*TMath::Power(sin2beta,6.);

  //         + G*TMath::Power(cos,4.))*(1. + H*TMath::Power(sin1beta,2.) + I*TMath::Power(sin1beta,4.) + J*TMath::Power(sin2beta,2.) + K*TMath::Power(sin2beta,4.) + L*TMath::Power(sin1beta,6.) + M*TMath::Power(sin2beta,6.));// + L*TMath::Power(cos1beta,2.) + M*TMath::Power(cos1beta,4.) + N*TMath::Power(cos2beta,2.) + O*TMath::Power(cos2beta,4.));
}
