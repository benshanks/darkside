//Math Helper Functions for Beta Decay
// B. Shanks, 9/25/14

#include "TComplex.h"
#include <math.h>
#include "TMath.h"
#include "TF1.h"
#include "Math/GSLIntegrator.h"
#include "Math/WrappedTF1.h"

//copying some units from CLHEP in case this code is reused by someone without CLHEP
static const double fine_structure_const = 1/137.035989;
static const double electron_mass_c2 = 0.511; //MeV
static const double pi = TMath::Pi();

static const double krPeakCentroidMeV = 0.0415;//MeV

//------------------------------------------------------------------
//Math helpers
//------------------------------------------------------------------

//root apparently doesn't have any implementation of gamma for complex numbers.  Here are a few numerical implementations I had lying around.

TComplex ComplexLnGammaLanczos(const TComplex& z)
{
    static double lanczos_7_c[9] = {
        0.99999999999980993227684700473478,
        676.520368121885098567009190444019,
        -1259.13921672240287047156078755283,
        771.3234287776530788486528258894,
        -176.61502916214059906584551354,
        12.507343278686904814458936853,
        -0.13857109526572011689554707,
        9.984369578019570859563e-6,
        1.50563273514931155834e-7
    };
    double zr = z.Re();
    double zi = z.Im();
    int k;
    double Ag_r, Ag_i;
    zr -= 1.0; /* Lanczos writes z! instead of Gamma(z) */
    Ag_r = lanczos_7_c[0];
    Ag_i = 0.0;
    for(k=1; k<=8; k++) {
        double R = zr + k;
        double I = zi;
        double a = lanczos_7_c[k] / (R*R + I*I);
        Ag_r +=  a * R;
        Ag_i -=  a * I;
    }
    TComplex zDummy(zr+7.5, zi);
    zDummy = TComplex::Log(zDummy);
    double log1_r = zDummy.Re(), log1_i = zDummy.Im();
    zDummy = TComplex::Log(zDummy(Ag_r, Ag_i));
    double logAg_r = zDummy.Re(), logAg_i = zDummy.Im();
    
    static double logRootTwoPi = log(sqrt(2.*pi));
    double yr = (zr+0.5)*log1_r - zi*log1_i - (zr+7.5) + logRootTwoPi + logAg_r;
    double yi = zi*log1_r + (zr+0.5)*log1_i - zi + logAg_i;
    return TComplex(yr, yi);
}

TComplex ComplexLnGamma(const TComplex& z)
{
    double zr = z.Re();
    double zi = z.Im();
    if(zr <= 0.5) {
        /* Transform to right half plane using reflection;
         * in fact we do a little better by stopping at 1/2.
         */
        TComplex zDummy(1.-zr, -zi);
        zDummy = ComplexLnGammaLanczos(zDummy);
        double a = zDummy.Re(), b = zDummy.Im();
        zDummy = TComplex::Log(TComplex::Sin(zDummy(pi*zr, pi*zi)));
        double lnsin_r = zDummy.Re(), lnsin_i = zDummy.Im();
        static double lnpi = log(pi);
        double yr = lnpi - lnsin_r - a;
        double yi = -lnsin_i - b;
        return TComplex(yr, yi);
    }
    else {
        /* otherwise plain vanilla Lanczos */
        return ComplexLnGammaLanczos(z);
    }
}

TComplex GammaBKF(TComplex z)
{
    // complex gamma function
    // modified f2c of CERNLIB cgamma64.F
    
    Double_t zr=z.Re();
    
    if (z.Im()==0.0) {
        if (zr==0.0) {
            TComplex w(0.0,0.0);
            return w;
        } else if (-TMath::Abs(zr)==TMath::Floor(zr)) {
            TComplex w(0.0,0.0);
            return w;
        }
    }
    
    TComplex u;
    TComplex v;
    
    if (zr>=1.0) {
        u=1.0;
        v=z;
    } else if (zr>=0.0) {
        u=1.0/z;
        v=1.0+z;
    } else {
        u=1.0;
        v=1.0-z;
    }
    
    const Double_t c1=2.50662827463100050;
    const Double_t c2[16]={
        41.624436916439068,-51.224241022374774,
        11.338755813488977, -0.747732687772388,
        0.008782877493061, -0.000001899030264,
        0.000000001946335, -0.000000000199345,
        0.000000000008433,  0.000000000001486,
        -0.000000000000806,  0.000000000000293,
        -0.000000000000102,  0.000000000000037,
        -0.000000000000014,  0.000000000000006
    };
    
    TComplex w(1.0,0.0);
    TComplex s(c2[0],0.0);
    TComplex cpi(TMath::Pi(),0.0);
    
    for (int i=1;i<16;++i) {
        w*=(v-((Double_t)i))/(v+((Double_t)(i-1)));
        s+=c2[i]*w;
    }
    
    w=v+4.5;
    w=c1*s*TComplex::Exp((v-0.5)*TComplex::Log(w)-w);
    
    if (zr<0.0) w=cpi/(TComplex::Sin(cpi*z)*w);
    
    return w*u;
}


double AbsComplexGamma(double zr, double zi)
{
    //return TComplex::Abs(TComplex::Exp(ComplexLnGamma(TComplex(zr, zi))));
    return TComplex::Abs(GammaBKF(TComplex(zr, zi)));
}

double LogAbsComplexGamma(double zr, double zi)
{
    return log(TComplex::Abs(TComplex::Exp(ComplexLnGamma(TComplex(zr, zi)))));
}

//------------------------------------------------------------------
//Beta Decay stuff
//------------------------------------------------------------------

double Fermi(double E, double Z, double A)
{
    //uses log gamma
    
    if(E <= 0.) return 0.;
    double alphaZ = Z*fine_structure_const;
    double w = E/electron_mass_c2 + 1.;
    double p = sqrt(w*w-1.);
    double n = alphaZ*w/p;
    double k = sqrt(1.-alphaZ*alphaZ);
    double r = 0.426*fine_structure_const*pow(A,1./3.);
    double norm = 2.*(1.+k)/pow(TMath::Gamma(1.+2.*k), 2);
    return norm*pow(2.*p*r, 2.*k-2.)*exp(pi*n + 2.*log(AbsComplexGamma(k,n)));
}

double FermiBKF(double E, double Z, double A)
{
    //uses normal gamma (BKF)
    
    if(E <= 0.) return 0.;
    double alphaZ = Z*fine_structure_const;
    double w = E/electron_mass_c2 + 1.;
    double p = sqrt(w*w-1.);
    double n = alphaZ*w/p;
    double k = sqrt(1.-alphaZ*alphaZ);
    double r = 0.426*fine_structure_const*pow(A,1./3.);
    double norm = 2.*(1.+k)/pow(TMath::Gamma(1.+2.*k), 2);
    return norm*pow(2.*p*r, 2.*k-2.)*exp(pi*n)* pow(AbsComplexGamma(k,n),2);
}

double BetaSpec(double E, double Z, double A, double Q){
    //without forbiddenness correction
    if(E < 0 || E > Q) return 0;
    
    return Fermi(E,Z,A) * pow((Q - E), 2.) * (E + electron_mass_c2) * sqrt(E*E + 2.*electron_mass_c2*E);
}

double BetaSpec(double* E, double* Z, double* A, double* Q){
    //without forbiddenness correction
    
    return BetaSpec(*E, *Z, *A, *Q);
}

double ForbiddenBetaSpec(double E, double Z, double A, double Q){
    //Use the forbiddenness correction
    //if a==b, absorb a and b into overall fit, OK to set both to 1 here
    double a = 1;
    double b = 1;
    
    double correction = a*(pow(E+electron_mass_c2,2.)-pow(electron_mass_c2,2.)) + b*pow((Q-E),2.);
    return BetaSpec(E,Z,A,Q) * correction;
}

double ForbiddenBetaSpec(double* E, double* Z, double* A, double* Q){
    //Use the forbiddenness correction
    return ForbiddenBetaSpec(*E, *Z, *A, *Q);
}

double ForbiddenBetaSpecPEPdf(double PE, double Z, double A, double Q, double N, double LY){
    //converts from PE
    double E = PE /LY;
    return N*ForbiddenBetaSpec(E,Z,A,Q);
}


double ForbiddenBetaSpecPEPdf(double* x, double* p){
    //Convert forbidden beta spec to PE
    double PE = x[0]; //returns in PHOTOELECTRONS
    
    
    double Z = p[0];
    double A = p[1];
    double Q = p[2];
    double N = p[3]; //changes normalization. does nothing right now.
    double LY = p[4]; //light yield parameter
    
    double E = PE /LY;
    
    return N*ForbiddenBetaSpec(E,Z,A,Q);
}

//------------------------------------------------------------------
// Convolving in resolution, adding in the Kr peak
//------------------------------------------------------------------


double GaussResSig( double *x, double *p )
{
    double PE = x[0];
    
    double sig0 = p[0];
    double sig1 = p[1];
    double sig2 = p[2];
    
   return sqrt(sig0*sig0 + (1.+sig1*sig1)*TMath::Max(PE, 0.)  + sig2*sig2* pow(TMath::Max(PE, 0.),2) );
}


//double GaussRes( double *x, double *p )
//{
//    double PE = x[0];
//    
//    double sig0 = p[0];
//    double sig1 = p[1];
//    double sig2 = p[2];
//    
//    fGaussResSig->SetParameter(0, sig0);
//    fGaussResSig->SetParameter(1, sig1);
//    fGaussResSig->SetParameter(2, sig2);
//    
//    double sig = fGaussResSig->Eval(PE);
//    
//    return 1/sqrt(2.*pi) / sig * exp(-0.5*( pow(PE/sig),2)) );
//}


double GaussTimesBeta( double *x, double *p )
{
    //include here any noise contribution
    
    double PE = p[0]; //this is the PE you're evaluating at
    
    double Z = p[1];
    double A = p[2];
    double Q = p[3];
    double N = p[4]; //overall amplitude
    double LY = p[5]; //light yield (assume constant)
    
    double t = x[0]; //integration variable
    
    double sig = p[6];
    
    double bg = p[7];
    //return TMath::Gaus(PE,t, sig, 1);
    return (ForbiddenBetaSpecPEPdf(t, Z,A,Q,N,LY)+bg) * TMath::Gaus(PE,t, sig, 1);
}

double BetaSpecPdfRes(double* x, double* p)
{
    double PE = x[0]; //returns in PHOTOELECTRONS
    
    double Z = p[0];
    double A = p[1];
    double Q = p[2];
    
    double N = p[3]; //overall amplitude
    double LY = p[4]; //light yield (assume constant)
    
    double sig0 = p[5];
    double sig1 = p[6];
    double sig2 = p[7];
    
    double bg = p[8]; //flat background spectrum to add into convolution
    
    //resolution function
    //double sig = sqrt(sig0*sig0 + (1+sig1*sig1)*TMath::Max(PE, 0.)  + sig2*sig2*pow(TMath::Max(PE, 0.),2));
    
    
    double sig = GaussResSig(x, &p[5]);
    
    Double_t nSig = 5.; //number of sigma to convolve in
    Double_t xlow = PE - nSig*sig;
    Double_t xupp = PE + nSig*sig;
    
    //better use GSL fancy integration
    
    TF1 fGaussTimesBeta("fGaussTimesBeta", GaussTimesBeta, xlow, xupp, 8);
    fGaussTimesBeta.SetParameter(0,PE);
    fGaussTimesBeta.SetParameter(1,Z);
    fGaussTimesBeta.SetParameter(2,A);
    fGaussTimesBeta.SetParameter(3,Q);
    fGaussTimesBeta.SetParameter(4,N);
    fGaussTimesBeta.SetParameter(5,LY);
    fGaussTimesBeta.SetParameter(6,sig);
    fGaussTimesBeta.SetParameter(7,bg);
    
    ROOT::Math::GSLIntegrator ig( ROOT::Math::IntegrationOneDim::kADAPTIVE);

    const ROOT::Math::WrappedTF1 wf(fGaussTimesBeta);
    
    ig.SetFunction(wf);
    ig.SetRelTolerance(0.001);
    
    double conv = ig.Integral(xlow, xupp);
    return conv;
}

double Kr83Peak(double* x, double* p){
    //for an INDEPENDENT gaussian
    double PE = x[0]; //returns in PHOTOELECTRONS
   
    double LY = p[0];
    double KrAmp = p[1]; //photopeak amplitude
    double KrSigma = p[2];
    
    double KrCentroid = LY*krPeakCentroidMeV;
    
    return KrAmp*TMath::Gaus(PE,KrCentroid,KrSigma, 1);
    
}

double Kr83PeakSharedSigma(double* x, double* p){
    //for gaussian sigma tied to beta spectrum
    
    double PE = x[0]; //returns in PHOTOELECTRONS
    
    double LY = p[0];
    double sig0 = p[1];
    double sig1 = p[2];
    double sig2 = p[3];
    double KrAmp = p[4]; //photopeak amplitude
    
    double KrCentroid = LY*krPeakCentroidMeV;
        
    double sig = GaussResSig(x, &p[1]);
    //can't use root to normalize the gaussian, because the sigma changes as a function of E.
    //instead normalize to the sigma at median.
    
    double sigMedian = GaussResSig(&KrCentroid, &p[1]);
    
    return 1/sqrt(2*pi)/sigMedian * KrAmp*TMath::Gaus(PE,KrCentroid,sig, 0);
    
}

double Ar39andKr83Spec(double* x, double* p)
{
    //ar beta and kr photopeak fit separately
    
    double PE = x[0]; //returns in PHOTOELECTRONS

    //for the beta
    double Z = p[0];
    double A = p[1];
    double Q = p[2];
    
    double N = p[3]; //overall amplitude of the beta
    double LY = p[4]; //light yield (assume constant)
    
    double sig0 = p[5];
    double sig1 = p[6];
    double sig2 = p[7];
                  
    //flat background
    double bg = p[8]; //photopeak amplitude
    
    //for the gamma
    double KrAmp = p[9]; //photopeak amplitude
    double KrSigma = p[10];
    
    double krParam[3] = {LY, KrAmp, KrSigma};
    
    // Sum of beta and peak function
    double sum = BetaSpecPdfRes(x, p) + Kr83Peak(x, krParam);
                                                      
    return sum;
}

double Ar39andKr83SpecCombined(double* x, double* p)
{
    //uses a shared sigma parameter for kr and argon
    
    double PE = x[0]; //returns in PHOTOELECTRONS
    
    //for the beta
    double Z = p[0];
    double A = p[1];
    double Q = p[2];
    
    double N = p[3]; //overall amplitude of the beta
    double LY = p[4]; //light yield (assume constant)
    
    double sig0 = p[5];
    double sig1 = p[6];
    double sig2 = p[7];
                  
    //flat background
    double bg = p[8]; //photopeak amplitude
    
    //for the gamma
    double KrAmp = p[9]; //photopeak amplitude
    
    double krParam[5] = {LY, sig0, sig1, sig2, KrAmp};
    
    // Sum of beta and peak function
    double sum = BetaSpecPdfRes(x, p) + Kr83PeakSharedSigma(x,krParam);
    
    return sum;
}

