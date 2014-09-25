#include <TROOT.h>
#include <TStyle.h>
#include "TF1.h"

void plot_forbidden_beta(){
    TCanvas* c1 = new TCanvas();
    
//    //pull the data file into an artificial spectrum...
//    TH1F *h1 = new TH1F("h1","^{39}Ar Spectrum",101,0,0.54);
//    ReadArSpectrum(h1);

//    h1->SetStats(0);
//    h1->Draw();
//    c1->Update();
    
    //now fit it...
    
    gROOT->ProcessLine(".L beta_functions.cxx");
    
    gStyle->SetFillColor(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    
    
    TF1* fBeta = new TF1("fBeta", "ForbiddenBetaSpec(x, 19, 39, 0.565)", 0, 0.565);
    fBeta->Draw();
    
    c1->Update();
}