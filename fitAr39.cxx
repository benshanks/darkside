#include "Riostream.h"

#include "TH1.h"
#include "TComplex.h"
#include <cmath>
#include "TMath.h"


void fitAr39(){
    //pull the data file into an artificial spectrum...
    TH1F *h1 = new TH1F("h1","^{39}Ar Spectrum",101,0,0.54);
    ReadArSpectrum(h1);
    TCanvas* c1 = new TCanvas();
    h1->SetStats(0);
    h1->Draw();
    h1->Draw();
    c1->Update();
    
    //now fit it...
    
}


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
