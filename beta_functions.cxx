//Math Helper Functions for Beta Decay
// B. Shanks, 9/25/14

#include "TComplex.h"
#include <cmath>
#include "TMath.h"

//copying some units from CLHEP in case this code is reused by someone without CLHEP
static const double fine_structure_const = 1/137.035989;
static const double electron_mass_c2 = 0.511;
static const double pi = TMath::Pi();

//------------------------------------------------------------------
//Math helpers
//------------------------------------------------------------------

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
    if(E <= 0.) return 0.;
    // Z = positive for electrons, negative for positrons
    //if(E <= 50.*eV) E = 50.*eV;
    double alphaZ = Z*fine_structure_const;
    double w = E/electron_mass_c2 + 1.;
    double p = sqrt(w*w-1.);
    double n = alphaZ*w/p;
    double k = sqrt(1.-alphaZ*alphaZ);
    double r = 0.426*fine_structure_const*pow(A,1./3.);
    double norm = 2.*(1.+k)/pow(TMath::Gamma(1.+2.*k), 2);
    return norm*pow(2.*p*r, 2.*k-2.)*exp(pi*n + 2.*log(AbsComplexGamma(k,n)));
}

double BetaSpec(double E, double Z, double A, double Q){
    //without forbiddenness correction
    return Fermi(E,Z,A) * pow((Q - E), 2.) * (E + electron_mass_c2) * sqrt(E + 2.*electron_mass_c2*E);
}

double BetaSpec(double* E, double* Z, double* A, double* Q){
    //without forbiddenness correction
    return BetaSpec(*E, *Z, *A, *Q);
}

double ForbiddenBetaSpec(double E, double Z, double A, double Q){
    //without forbiddenness correction
    //absorb a and b into overall fit, if a==b, OK to set both to 1
    double a = 1;
    double b = 1;
    
    double correction = a*(pow(E+electron_mass_c2,2.)-pow(electron_mass_c2,2.)) + b*pow((Q-E),2.);
    return BetaSpec(E,Z,A,Q) * correction;
}

double ForbiddenBetaSpec(double* E, double* Z, double* A, double* Q){
    //Use the forbiddenness correction
    return ForbiddenBetaSpec(*E, *Z, *A, *Q);
}
