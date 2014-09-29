//Trying to fit the whole beta spectrum with convolution for resolution
//Constant light yield, no e field

//B. Shanks, 9/26/14

#include "Riostream.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include <math.h>

#include "beta_functions.cxx"

void fitSpectrum(){
    gROOT->ProcessLine(".L beta_functions.cxx");
    gStyle->SetFillColor(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptStat(0);
    
  
    TFile *f = new TFile("kr.root");
//    
//    //hist is in photoelectrons
    TH1F *hpx = (TH1F*)f->Get("hS1_0");
    
    //-----Set up general fit values
    double fitMin = 200;
    double fitMax = 5000;
    double backgroundGuess = 0; //flat background spectrum

    
    //-----Set up the Beta Spectrum fit

    //Q value is at 0.565
    double ArZ = 19;
    double ArA = 39;
    double ArQ = 0.565; //in MeV

    double ArAmp = 1.5*pow(10, 4);
    double ArLY = 8000.;
    
    double ArSig0 = 175;
    double ArSig1 = 0.4;
    double ArSig2 = 0.0435;
    
    //-----Set up the Kr fit guesses
    
    //peak is at 41.5 keV.  looks here like 200-450 PE is fair, centered at 325
    
    double KrampGuess = 1.1*pow(10, 5);  ; //WITH the 1/sqrt(2pi)/sigma factor
    double KrcentroidGuess = 328;
    double KrsigmaGuess = 24.2; //converting from FWHM
    
    
    TF1* fAr39andKr83Spec = new TF1("fAr39andKr83Spec", Ar39andKr83Spec, fitMin, fitMax, 12);
    
//    fAr39andKr83Spec->SetParameter(0, ArZ);
//    fAr39andKr83Spec->SetParameter(1, ArA);
//    fAr39andKr83Spec->SetParameter(2, ArQ);
    fAr39andKr83Spec->SetParameter(3, ArAmp);
    fAr39andKr83Spec->SetParameter(4, ArLY);
    fAr39andKr83Spec->SetParameter(5, ArSig0);
//    fAr39andKr83Spec->SetParameter(6, ArSig1);
    fAr39andKr83Spec->SetParameter(7, ArSig2);
    fAr39andKr83Spec->SetParLimits(7, 0., 100.);
    
    fAr39andKr83Spec->SetParameter(8, KrampGuess);
    fAr39andKr83Spec->SetParameter(9, KrcentroidGuess);
    fAr39andKr83Spec->SetParameter(10, KrsigmaGuess);
//    fAr39andKr83Spec->SetParameter(11, backgroundGuess);
    
//    fAr39andKr83Spec->SetParNames("Z","A", "Q","Argon Amp.", "Light Yield", "Sig0", "Sig1", "Sig2", "Kr Amp.", "Kr Centroid","Kr Sigma", "BG const.");
    
    //fix the sigma1 from laser measurements.
    fAr39andKr83Spec->FixParameter(6, ArSig1);
    
    //don't float the Z, A, or Q.
    fAr39andKr83Spec->FixParameter(0, ArZ);
    fAr39andKr83Spec->FixParameter(1, ArA);
    fAr39andKr83Spec->FixParameter(2, ArQ);
    
    //not using this right now
    fAr39andKr83Spec->FixParameter(11, backgroundGuess);
    
    //    fBetaSpecPdfRes->SetParameter(0, ArZ);
    //    fBetaSpecPdfRes->SetParameter(1, ArA);
    //    fBetaSpecPdfRes->SetParameter(2, ArQ);
    //    fBetaSpecPdfRes->SetParameter(3, ArAmp);
    //    fBetaSpecPdfRes->SetParameter(4, ArLY);
    //    fBetaSpecPdfRes->SetParameter(5, ArSig0);
    //    fBetaSpecPdfRes->SetParameter(6, ArSig1);
    //    fBetaSpecPdfRes->SetParameter(7, ArSig2);
    
    //------------------------------------------------------------
    //-----Fit the sum
    //------------------------------------------------------------
    //hpx->Fit("fAr39andKr83Spec", "mern0");
    hpx->Fit("fAr39andKr83Spec", "mern0V");

    
    //Plot the fit result (redoing it explicitly because some funny business is happening w/ the hist color on load from .root file)
    
    TCanvas* c1 = new TCanvas();
    
    //TODO: this could be prettier...
    hpx->SetTitle("Fit to the ^{83}Kr Peak and ^{39}Ar Spectrum");
    hpx->GetXaxis()->SetTitle("PE");
    hpx->GetXaxis()->SetRange(0,fitMax);
    hpx->GetYaxis()->SetTitle("Counts");
    //hpx->GetYaxis()->CenterTitle();
    //hpx->GetXaxis()->SetTitleSize(0.06);
    //hpx->GetYaxis()->SetTitleSize(0.06);
    hpx->SetLineColor(kBlue);
    hpx->Draw();
    
    fAr39andKr83Spec->SetLineColor(kRed);
    fAr39andKr83Spec->SetNpx(5000);
    fAr39andKr83Spec->Draw("same");
//    
    c1->Update();
    
    c1->Print("spectrum.pdf");
    c1->Print("spectrum.root");
}

int main(){
    fitSpectrum();
}
