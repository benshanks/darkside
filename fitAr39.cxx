//Pulls in normalized spectrum from Ar39.dat, generates unnormalized hist,
//and fits forbidden & normal beta spectrum to it.
//B. Shanks, 9/24/14

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
    gStyle->SetOptStat(0);
    
    
    TCanvas* c1 = new TCanvas();
    
    //pull the data file into an artificial spectrum...
    TH1F *h1 = new TH1F("h1","^{39}Ar Spectrum",101,0,0.54);
    ReadArSpectrum(h1);
//    h1->SetStats(0);
    
    //now fit it with forbidden and non-forbidden...
    
    TF1* fBeta = new TF1("fBeta", "[0]*BetaSpec(x, 19, 39, 0.565)", 0, 0.565);
    fBeta->SetParameter(0,50000);
    fBeta->Draw();
    h1->Fit("fBeta", "0nr");
    c1->Update();
    
    TF1* fBetaForb = new TF1("fBetaForb", "[0]*ForbiddenBetaSpec(x, 19, 39, 0.565)", 0, 0.565);
    fBetaForb->SetParameter(0,50000);
    fBetaForb->SetLineColor(kMagenta);
    fBetaForb->Draw("SAME");
    h1->Fit("fBetaForb", "0nr");
    
    c1->Update();
    
    h1->Draw("Same");
    c1->Update();
    
    double histMean = h1->GetMean();
    double betaMean = fBeta->Mean(0,0.565);
    double betaForbMean = fBetaForb->Mean(0,0.565);
    
    char legHist[500];
    char legBeta[500];
    char legBetaMean[500];
    sprintf(legHist,".Dat File (Mean %.4f)", histMean);
    sprintf(legBeta,"Allowed Beta (Mean %.4f)", betaMean);
    sprintf(legBetaMean,"Forbidden Beta (Mean %.4f)", betaForbMean);
    
    TLegend* leg2 = new TLegend(.5, .75, .9, .9);
    leg2->SetTextFont(72);
    leg2->AddEntry(legHist,legHist,"l");
    leg2->AddEntry(fBeta,legBeta,"l");
    leg2->AddEntry(fBetaForb,legBetaMean,"l");
    leg2->Draw();
    
    c1->Print("Argon_Beta_Spectra.pdf");
    c1->Print("Argon_Beta_Spectra.root");
    
}


//------------------------------------------------------------------
//read in the .dat file
//------------------------------------------------------------------

//borrowing heavily from root tutorial basic.c
void ReadArSpectrum(TH1F* h1) {
    
    ifstream in;
    in.open("Ar39.dat");
    Float_t energy,spec;
    
    int specMult = 1000; //whatever...
    
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
