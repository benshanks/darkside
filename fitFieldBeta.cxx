//Fits the whole spectrum,  WITH e field, no Kr peak
//Constant light yield

//B. Shanks, 9/30/14

#include "Riostream.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include <math.h>

#include "beta_functions.cxx"

void fitField(){
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
    TH1F *hpx = (TH1F*)f->Get("hS1_3");
    
    //-----Set up general fit values
    double fitMin = 100;
    double fitMax = 5000;
    double backgroundGuess = 0; //flat background spectrum

    
    //-----Set up the Beta Spectrum fit

    //Q value is at 0.565
    double ArZ = 19;
    double ArA = 39;
    double ArQ = 0.565; //in MeV

    double ArAmp = 1.4*pow(10, 4);
    double ArLY = 5700;
    
    double ArSig0 = 150;
    double ArSig1 = 0.4;
    double ArSig2 = 0.016;
    
    //-----Set up the Kr fit guesses
    
    //peak is at 41.5 keV.  looks here like 200-450 PE is fair, centered at 325
    
    double KrampGuess = 6.45*pow(10, 4);  ; //WITH the 1/sqrt(2pi)/sigma factor
    double KrcentroidGuess = 291;
    double KrsigmaGuess = 17.5;
    
    
    //fit only to the beta spectrum
    TF1* fBetaSpec = new TF1("fBetaSpec", BetaSpecPdfRes, fitMin, fitMax, 8);
    
    //fix the sigma1 from laser measurements.
    fBetaSpec->FixParameter(6, ArSig1);
    
    //don't float the Z, A, or Q.
    fBetaSpec->FixParameter(0, ArZ);
    fBetaSpec->FixParameter(1, ArA);
    fBetaSpec->FixParameter(2, ArQ);
    
    fBetaSpec->SetParameter(3, ArAmp);
    fBetaSpec->SetParameter(4, ArLY);
    fBetaSpec->SetParameter(5, ArSig0);
    fBetaSpec->SetParameter(7, ArSig2);
    fBetaSpec->SetParLimits(7, 0., 100.);

    
    fBetaSpec->SetParNames("Z","A", "Q","Argon Amp.", "Light Yield", "Sig0", "Sig1", "Sig2");

    
    //------------------------------------------------------------
    //-----Fit the spectrum
    //------------------------------------------------------------
    //hpx->Fit("fAr39andKr83Spec", "mern0");
    hpx->Fit("fBetaSpec", "mern0V");

    
    //Plot the fit result (redoing it explicitly because some funny business is happening w/ the hist color on load from .root file)
    
    TCanvas* c1 = new TCanvas();
    
    //TODO: this could be prettier...
    hpx->SetTitle("Fit to the ^{39}Ar Spectrum");
    hpx->GetXaxis()->SetTitle("PE");
    hpx->GetXaxis()->SetRange(0,fitMax);
    hpx->GetYaxis()->SetTitle("Counts");
    //hpx->GetYaxis()->CenterTitle();
    //hpx->GetXaxis()->SetTitleSize(0.06);
    //hpx->GetYaxis()->SetTitleSize(0.06);
    hpx->SetLineColor(kBlue);
    hpx->Draw();
    
    fBetaSpec->SetLineColor(kRed);
    fBetaSpec->SetNpx(5000);
    fBetaSpec->Draw("same");
//    
    c1->Update();
    
    c1->Print("spectrum_field_betaonly.pdf");
    c1->Print("spectrum_field_betaonly.root");
}

int main(){
    fitField();
}
