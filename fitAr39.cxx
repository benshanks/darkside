#include "Riostream.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TH1.h"
#include "TF1.h"

void fitAr39(){
    gROOT->ProcessLine(".L beta_functions.cxx");
    gStyle->SetFillColor(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    
    
    TCanvas* c1 = new TCanvas();
    
    //pull the data file into an artificial spectrum...
    TH1F *h1 = new TH1F("h1","^{39}Ar Spectrum",101,0,0.54);
    ReadArSpectrum(h1);
    h1->SetStats(0);
    h1->Draw();
    c1->Update();
    
    //now fit it...
    
    TF1* fBeta = new TF1("fBeta", "[0]*ForbiddenBetaSpec(x, 19, 39, 0.565)", 0, 0.565);
    fBeta->SetParameter(0,50000);
    fBeta->Draw();
    h1->Fit("fBeta", "r");
    
    c1->Update();
}


//------------------------------------------------------------------
//read in the .dat file
//------------------------------------------------------------------

//borrowing heavily from root tutorial basic.c
void ReadArSpectrum(TH1F* h1) {
    
    ifstream in;
    in.open("Ar39.dat");
    Float_t energy,spec;
    
    int specMult = 10000; //whatever...
    
    while (1) {
        in >> energy >> spec;
        if (!in.good()) break;
        
        int specNum = TMath::Nint(spec*specMult);
        
        for(size_t iSpec=0;iSpec<=specNum;iSpec++){
            h1->Fill(energy);
        }
    }
    in.close();
}
