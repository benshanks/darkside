//Just fits the Kr83 peak to a Gaussian
// B. Shanks, 9/23/14


void fitKr83()
{
    TFile *f = new TFile("kr.root");
    
    //hist is in photoelectrons
    TH1F *hpx = (TH1F*)f->Get("hS1_0");
    
    //peak is at 41.5 keV.  looks here like 200-450 PE is fair, centered at 325
    double PEFitmin=200;//in PE
    double PEFitmax=450;
    
    double PEcentroidGuess = 325;
    double PEsigmaGuess = 75/2.3548; //converting from FWHM
    
    double PEbgGuess = 440;
    double PEampGuess = 110000  ; //WITH the 1/sqrt(2pi)/sigma factor

    
    TF1 *gausFunc = new TF1("PeakKr83","[0] + [1]*1/sqrt(2.*TMath::Pi())/[3]*exp(-0.5*((x-[2])/[3])^2)",PEFitmin,PEFitmax);
    gausFunc->SetParameter(0,PEbgGuess);
    gausFunc->SetParameter(1,PEampGuess);
    gausFunc->SetParameter(2,PEcentroidGuess);
    gausFunc->SetParameter(3,PEsigmaGuess);
    
    gausFunc->SetParNames("Constant","Amplitude", "Centroid","Sigma");
    
    hpx->Fit("PeakKr83", "mern0");
    
    
    //Plot the fit result (redoing it explicitly because some funny business is happening w/ the hist color on load from .root file)
    TCanvas* c1 = new TCanvas();
    
    //TODO: this could be prettier...
    hpx->SetStats(0);
    hpx->SetTitle("Fit to the ^{83}Kr Peak");
    hpx->GetXaxis()->SetTitle("PE");
    //hpx->GetXaxis()->SetTitleSize(0.06);
    hpx->GetXaxis()->SetRange(PEFitmin,PEFitmax);
    hpx->GetYaxis()->SetTitle("Counts");
    //hpx->GetYaxis()->CenterTitle();
    //hpx->GetYaxis()->SetTitleSize(0.06);
    hpx->SetLineColor(kBlue);
    hpx->Draw();
    gausFunc->SetLineColor(kRed);
    gausFunc->Draw("SAME");
    c1->Update();
    
    
}