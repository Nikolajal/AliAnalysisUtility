//
//  Part of the AliAnalysisUtility package
//
//  General File to be included in the individual analysis
//
//  Author              Nicola Rubini
//  Created             22/11/2021
//  Last modified       22/11/2021
#ifndef ALIANALYSISUTILITY_H
#define ALIANALYSISUTILITY_H
//
//>>    Functions File
#include "AAU_GlobalUtilVariables.h"
#include "AAU_GlobalUtilFunctions.h"
#include "AAU_Style.h"
#include "AAU_Functions.h"
#include "AAU_Histograms.h"
#include "AAU_Resolution.h"
//
using namespace std;
using namespace RooFit;
//
//------------------------------//
//    FIT FUNCTIONS UTILITIES   //
//------------------------------//
//
//  --  --  FIT Utilities  --  --  //
//
//_____________________________________________________________________________
//                                                  //  Relates to global variable kFitScarseHisto, kScarseHistoDef, kScarseHistoMin
template < class Tclass >
bool                    fIsWorthFitting             ( Tclass * aTarget )    {
    if ( !aTarget ) {
        cout << "[ERROR] You are trying to fit a null histogram! Skipping this one..." << endl;
        return false;
    }
    if (  aTarget->GetEntries() == 0. ) {
        cout << "[ERROR] You are trying to fit an empty histogram! Skipping " << aTarget->GetName() << endl;
        return false;
    }
    if (  aTarget->GetEntries() <= kScarseHistoMin > 0 ? kScarseHistoMin : kScarseHistoDef*aTarget->GetNbinsX()*aTarget->GetNbinsY()*aTarget->GetNbinsZ() ) {
        cout << "[WARNING] You are trying to fit a scarsely populated histogram!" << endl;
        if ( kFitScarseHisto )  {
            cout << "[INFO] You chose to Fit it anyway, so I'm trying" << endl;
            return true;
        }
        else    {
            cout << "[INFO] You chose to skip this type of Fit, so I'm skipping" << endl;
            return false;
        }
    }
    return true;
}
//
//_____________________________________________________________________________
//                                                  //  Relates to global variable kLooseErrors
template < class Tclass >
Tclass                 *fLooseErrors                ( Tclass * aTarget, Float_t fLooseErrors = kLooseErrors )    {
    if ( fLooseErrors != kLooseErrors ) cout << "[INFO] Loosening factor for errors is different from the global set one, using the one provided in function call." << endl;
    if ( fLooseErrors < 1 )             cout << "[WARNING] Loosening factor for errors is less than 1, meaning your errors are less than what they really are. This might be causing problems to the fit function." << endl;
    if ( fLooseErrors == 1 )            cout << "[INFO] Loosening factor for errors is set to one, this is rendering this function useless. I'm returning a copy of the input." << endl;
    Tclass     *fResult =   (Tclass*)aTarget->Clone();
    if ( fLooseErrors == 1 )            return  fResult;
    Int_t       fNBins  =   aTarget->GetNbinsX()*aTarget->GetNbinsY()*aTarget->GetNbinsZ();
    for ( Int_t iBin = 1; iBin <= fNBins; iBin++ )  {
        fResult ->  SetBinError(iBin,fLooseErrors*(aTarget->GetBinError(iBin)));
    }
    return fResult;
}
//
//_____________________________________________________________________________
//                                                  //  Relates to global variable kLooseErrors
TGraphAsymmErrors      *fLooseErrors                ( TGraphAsymmErrors * aTarget, Float_t fLooseErrors = kLooseErrors )    {
    if ( fLooseErrors != kLooseErrors ) cout << "[INFO] Loosening factor for errors is different from the global set one, using the one provided in function call." << endl;
    if ( fLooseErrors < 1 )             cout << "[WARNING] Loosening factor for errors is less than 1, meaning your errors are less than what they really are. This might be causing problems to the fit function." << endl;
    if ( fLooseErrors == 1 )            cout << "[INFO] Loosening factor for errors is set to one, this is rendering this function useless. I'm returning a copy of the input." << endl;
    TGraphAsymmErrors     *fResult =   new TGraphAsymmErrors(*aTarget);
    if ( fLooseErrors == 1 )            return  fResult;
    Int_t       fNBins  =   aTarget->GetN();
    for ( Int_t iBin = 0; iBin < fNBins; iBin++ )  {
        fResult ->  SetPointEYhigh(iBin,fLooseErrors*(aTarget->GetErrorYhigh(iBin-1)));
        fResult ->  SetPointEYlow(iBin,fLooseErrors*(aTarget->GetErrorYlow(iBin-1)));
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
bool                    fIsResultAcceptable         ( Double_t fResult, Double_t fError, Double_t fTolerance = kMaximumError )  {
    return (fError/fResult >= fTolerance/100. ? false : true);
}
//
//_____________________________________________________________________________
//
void                    uTestGoodParameters         ( TGraphAsymmErrors *gTarget, TF1 *fFitFunction = fLevyTsallis )  {
    auto nIterations    =   10000;
    auto fMinLimit      =   5;
    auto fMaxLimit      =   25;
    if ( !gTarget ) { cout << "Invalid Graph to fit" << endl; return; }
    cout << "[INFO] Starting search for good parameters for function " << fFitFunction->GetName() << endl;
    cout << "[INFO] Maximum iterations: " << nIterations << endl;
    for ( Int_t iTer = 0; iTer < nIterations; iTer++ )  {
        cout << Form("[INFO] Starting iteration #%i",iTer) << endl;
        gTarget                             ->  Fit(fFitFunction,"IMREQ0SEX","");
        auto    fNDF    =   fFitFunction    ->  GetNDF();
        auto    fChi2   =   fFitFunction    ->  GetChisquare();
        auto    fCheck  =   fChi2/fNDF;
        cout << Form("[INFO] Currently we have %.1f Chi2/NDF",fCheck) << endl;
        if ( fCheck < fMinLimit )   { cout << Form("[INFO] The limit was set to %i, so I found a satisfying answer! bye bye...",fMinLimit) << endl; break; }
        else                { cout << Form("[INFO] The limit was set to %i, so not so great... Let's keep looking",fMinLimit) << endl;  }
        for ( Int_t iTer = 0; iTer < fFitFunction->GetNpar(); iTer++ )  {
            if ( fFitFunction->GetParError(iTer) == 0 ) continue;
            if ( fCheck > fMaxLimit )   {
                Double_t fMax, fMin;
                fFitFunction->GetParLimits(iTer,fMin,fMax);
                fFitFunction->SetParameter(iTer,uRandomGen->Uniform(fMin,fMax));
            }   else    {
                
            }
        }
    }
}
//
//_____________________________________________________________________________
//
//
//_____________________________________________________________________________
    
//!TODO: Sort things out, generalise where necessary
//
//----------------------------------------//
//    TO BE SORTED/COMPLETED UTILITIES    //
//----------------------------------------//
//
std::vector<TH1F*> _DUMP;
//
TH2F*
fInvertTH2F
 ( TH2F * hTarget ) {
    auto fResult = new TH2F(*hTarget);
    for ( Int_t iBin = 1; iBin <= hTarget->GetNbinsX(); iBin++ )    {
        for ( Int_t jBin = 1; jBin <= hTarget->GetNbinsY(); jBin++ )    {
            fResult->SetBinContent  ( iBin, jBin, hTarget->GetBinContent( jBin, iBin ) );
            fResult->SetBinError    ( iBin, jBin, hTarget->GetBinError  ( jBin, iBin ) );
        }
    }
    return fResult;
}
//_____________________________________________________________________________
TH1F*
uProjectX1D
( TH1F *hTarget )   {
    std::vector<Float_t> uAllPoints;
    for ( Int_t iBin = 0; iBin < hTarget->GetNbinsX(); iBin++ ) {
        auto fYValue    =   hTarget->GetBinContent(iBin+1);
        uAllPoints.push_back(fYValue);
    }
    return uBuildTH1F(uAllPoints,-1.,-1.);
}

//_____________________________________________________________________________
TH1F*
uProjectBarlow
( TH1F *hTarget, TH1F  *hTester )   {
    std::vector<Float_t> uAllPoints;
    for ( Int_t iBin = 0; iBin < hTarget->GetNbinsX(); iBin++ ) {
        auto fYValue    =   hTarget->GetBinContent  (iBin+1);
        auto fYError    =   hTarget->GetBinError    (iBin+1);
        auto fTValue    =   hTester->GetBinContent  (iBin+1);
        auto fTError    =   hTester->GetBinError    (iBin+1);
        uAllPoints.push_back(uBarlowPar(fYValue,fYError,fTValue,fTError));
    }
    return uBuildTH1F(uAllPoints,600,0.,-30.,30.);
}
TH1F*
uProjectBarlow
( TH2F *hTarget, TH2F  *hTester, Bool_t fDiag = false )   {
    std::vector<Float_t> uAllPoints;
    for ( Int_t iBin = 0; iBin < hTarget->GetNbinsX(); iBin++ ) {
        for ( Int_t jBin = 0; jBin < hTarget->GetNbinsY(); jBin++ ) {
            if ( ( iBin > jBin ) && fDiag ) continue;
            auto fYValue    =   hTarget->GetBinContent  (iBin+1,jBin+1);
            auto fYError    =   hTarget->GetBinError    (iBin+1,jBin+1);
            auto fTValue    =   hTester->GetBinContent  (iBin+1,jBin+1);
            auto fTError    =   hTester->GetBinError    (iBin+1,jBin+1);
            uAllPoints.push_back(uBarlowPar(fYValue,fYError,fTValue,fTError));
        }
    }
    return uBuildTH1F(uAllPoints,600,0.,-30.,30.);
}

//_____________________________________________________________________________
Bool_t
uIsRelevantVariation
( TH1F* hCheckHist, TString fFolder = "", TString fName = ""  )   {
    auto nPassedChecks  =   0;
    auto uMean          =   hCheckHist->GetMean();
    auto uSTDV          =   hCheckHist->GetRMS();
    auto uIntegralTot   =   hCheckHist->GetEntries();
    auto uIntegral1Sg   =   0.;
    auto uIntegral2Sg   =   0.;
    for ( Int_t iBin = 0; iBin <= hCheckHist->GetNbinsX(); iBin++ )    {
        auto fXValue    =   hCheckHist->GetBinCenter(iBin+1);
        if ( fabs( fXValue - uMean ) <= uSTDV   ) uIntegral1Sg += hCheckHist->GetBinContent(iBin+1);
        if ( fabs( fXValue - uMean ) <= 2*uSTDV ) uIntegral2Sg += hCheckHist->GetBinContent(iBin+1);
    }
    if ( fabs(uMean)                <= 0.10 ) nPassedChecks++;
    if ( uSTDV                      <= 1.10 ) nPassedChecks++;
    if ( uIntegral1Sg/uIntegralTot  >= 0.60 ) nPassedChecks++;
    if ( uIntegral2Sg/uIntegralTot  >= 0.88 ) nPassedChecks++;
    if ( !fFolder.IsNull() )   {
        gROOT->SetBatch();
        TCanvas    *cDrawResult =   new TCanvas();
        gStyle->SetOptStat(0);
        TLatex     *fLatex      =   new TLatex();
        hCheckHist  ->  SetTitle("");
        hCheckHist  ->  SetLineColor(kColors[2]);
        auto fMaxium = hCheckHist  ->  GetMaximum();
        hCheckHist  ->  SetMaximum(fMaxium*1.3);
        hCheckHist  ->  GetXaxis()->SetTitle("Barlow parameter");
        hCheckHist  ->  Draw("HIST");
        fLatex->SetTextFont(60);
        fLatex->SetTextSize(0.05);
        fLatex      ->  DrawLatexNDC(0.80,0.75,Form("%s",hCheckHist->GetName()));
        fLatex->SetTextFont(42);
        fLatex->SetTextSize(0.04);
        if ( fabs(uMean)                <= 0.10 )   fLatex  ->  DrawLatexNDC(0.18,0.84,Form("#color[8]{#mu:   %.2f < 0.10}",fabs(uMean)));
        else                                        fLatex  ->  DrawLatexNDC(0.18,0.84,Form("#color[2]{#mu:   %.2f < 0.10}",fabs(uMean)));
        if ( uSTDV                      <= 1.10 )   fLatex  ->  DrawLatexNDC(0.18,0.79,Form("#color[8]{#sigma:   %.2f < 1.10}",uSTDV));
        else                                        fLatex  ->  DrawLatexNDC(0.18,0.79,Form("#color[2]{#sigma:   %.2f < 1.10}",uSTDV));
        if ( uIntegral1Sg/uIntegralTot  >= 0.60 )   fLatex  ->  DrawLatexNDC(0.18,0.74,Form("#color[8]{INT_{1#sigma}: %.2f > 0.60}",uIntegral1Sg/uIntegralTot));
        else                                        fLatex  ->  DrawLatexNDC(0.18,0.74,Form("#color[2]{INT_{1#sigma}: %.2f > 0.60}",uIntegral1Sg/uIntegralTot));
        if ( uIntegral2Sg/uIntegralTot  >= 0.88 )   fLatex  ->  DrawLatexNDC(0.18,0.69,Form("#color[8]{INT_{2#sigma}: %.2f > 0.88}",uIntegral2Sg/uIntegralTot));
        else                                        fLatex  ->  DrawLatexNDC(0.18,0.69,Form("#color[2]{INT_{2#sigma}: %.2f > 0.88}",uIntegral2Sg/uIntegralTot));
        fLatex->SetTextSize(0.075);
        if ( nPassedChecks < 3 )    fLatex      ->  DrawLatexNDC(0.55,0.82,"#color[2]{CONSIDERED}");
        else                        fLatex      ->  DrawLatexNDC(0.60,0.82,"#color[8]{REJECTED}");
        gROOT       ->  ProcessLine (Form(".! mkdir -p %s",fFolder.Data()));
        cDrawResult ->  SaveAs((fFolder+fName+TString(Form("_%s.pdf",hCheckHist->GetName()))).Data());
        gROOT->SetBatch(kFALSE);
    }
    return nPassedChecks < 3;
}
std::vector<Bool_t>
uIsRelevantVariation
( TH1F* hStandard, std::vector<TH1F*> hVariations, TString fFolder = "", TString fName = "" )   {
    std::vector<Bool_t> uResults;
    for ( auto hVariation : hVariations ) {
        auto hTester    =   (TH1F*)hStandard->Clone();
        auto hTested    =   (TH1F*)hVariation->Clone();
        auto hCheckHist =   uProjectBarlow(hTester,hTested);
        hCheckHist      ->  SetName(hVariation->GetName());
        uResults.push_back( uIsRelevantVariation(hCheckHist,fFolder,fName) );
    }
    return uResults;
}
std::vector<Bool_t>
uIsRelevantVariation
( TH2F* hStandard, std::vector<TH2F*> hVariations, TString fFolder = "", TString fName = "", Bool_t fDiag = false  )   {
    std::vector<Bool_t> uResults;
    for ( auto hVariation : hVariations ) {
        auto hTester    =   (TH2F*)hStandard->Clone();
        auto hTested    =   (TH2F*)hVariation->Clone();
        auto hCheckHist =   uProjectBarlow(hTester,hTested,fDiag);
        hCheckHist      ->  SetName(hVariation->GetName());
        uResults.push_back( uIsRelevantVariation(hCheckHist,fFolder,fName) );
    }
    return uResults;
}

//_____________________________________________________________________________
std::vector<TH1F*>
uBuildVariationSpectra
( TH1F* hStandard, std::vector<TH1F*> hVariations, std::vector<Bool_t> uWhichIsToBeconsidered ) {
    std::vector<TH1F*> uResults;
    auto iTer = -1;
    for ( auto hVariation : hVariations ) {
        iTer++;
        if ( uWhichIsToBeconsidered.size() > iTer && !uWhichIsToBeconsidered.at(iTer) ) continue;
        auto hTester    =   (TH1F*)hStandard->Clone();
        hTester->Divide(hVariation,hStandard);
        uResults.push_back(hTester);
    }
    return uResults;
}
std::vector<TH1F*>
uBuildVariationSpectra
( TH1F* hStandard, std::vector<TH1F*> hVariations ) {
    std::vector<Bool_t> fDump;
    for ( auto fTested : hVariations ) fDump.push_back(kTRUE);
    return uBuildVariationSpectra(hStandard,hVariations,fDump);
}
std::vector<TH2F*>
uBuildVariationSpectra
( TH2F* hStandard, std::vector<TH2F*> hVariations, std::vector<Bool_t> uWhichIsToBeconsidered ) {
    std::vector<TH2F*> uResults;
    auto iTer = -1;
    for ( auto hVariation : hVariations ) {
        iTer++;
        if ( uWhichIsToBeconsidered.size() != 0 && !uWhichIsToBeconsidered.at(iTer) ) continue;
        auto hTester    =   (TH2F*)hStandard->Clone();
        hTester->Divide(hVariation,hStandard);
        uResults.push_back(hTester);
    }
    return uResults;
}
std::vector<TH2F*>
uBuildVariationSpectra
( TH2F* hStandard, std::vector<TH2F*> hVariations ) {
    std::vector<Bool_t> fDump;
    for ( auto fTested : hVariations ) fDump.push_back(kTRUE);
    return uBuildVariationSpectra(hStandard,hVariations,fDump);
}

//_____________________________________________________________________________
std::vector<TH1F*>
uBuildVariationBinByBin // Put bool to check or not for barlow, put bool to have point by point barlow check
( TH1F* hStandard, std::vector<TH1F*> hVariations, TString fFolder = "", TString fName = "", std::vector<Bool_t> fWhichToChoose = std::vector<Bool_t> () ) {
    std::vector<TH1F*>  hResults;
    auto    uSystVariationSpectra_  =   uBuildVariationSpectra(hStandard,hVariations,fWhichToChoose);
    for ( Int_t iPT = 0; iPT < hStandard->GetNbinsX(); iPT++ )    {
        hResults.push_back( new TH1F(Form("VarHist_%i",iPT),Form("VarHist_%i",iPT),2000,-1,+1) );
        for ( auto  hVariationSpectrum : uSystVariationSpectra_ ) {
            hResults.at(iPT)->Fill(hVariationSpectrum->GetBinContent(iPT+1)-1);
        }
        if ( !fFolder.IsNull() )    {
            gROOT->SetBatch();
            TCanvas    *cDrawResult =   new TCanvas();
            hResults.at(iPT)->Draw();
            cDrawResult->SaveAs((fFolder+fName+TString(Form("_%i.pdf",iPT))).Data());
            gROOT->SetBatch(kFALSE);
        }
    }
    return hResults;
}
std::vector<std::pair<TH1F*,Int_t>>
uBuildVariationBinByBin // Put bool to check or not for barlow, put bool to have point by point barlow check
( TH2F* hStandard, std::vector<TH2F*> hVariations, TString fFolder = "", TString fName = "", std::vector<Bool_t> fWhichToChoose = std::vector<Bool_t> (), Bool_t fDiag = false ) {
    std::vector<std::pair<TH1F*,Int_t>>  hResults;
    auto    uSystVariationSpectra_  =   uBuildVariationSpectra( hStandard, hVariations, fWhichToChoose );
    auto    iTH1    =   0;
    for ( Int_t iPT = 0; iPT < hStandard->GetNbinsX(); iPT++ )    {
        for ( Int_t jPT = 0; jPT < hStandard->GetNbinsY(); jPT++ )    {
            if ( ( iPT > jPT ) && fDiag ) continue;
            hResults.push_back( std::pair<TH1F*,Int_t>( new TH1F(Form("VarHist_%i_%i",iPT,jPT),Form("VarHist_%i_%i",iPT,jPT),2000,-1,+1) , hStandard->GetBin( iPT+1, jPT+1 ) ) );
            for ( auto  hVariationSpectrum : uSystVariationSpectra_ ) {
                hResults.at(iTH1).first->Fill(hVariationSpectrum->GetBinContent(iPT+1,jPT+1)-1);
            }
            if ( !fFolder.IsNull() )    {
                gROOT->SetBatch();
                TCanvas    *cDrawResult =   new TCanvas();
                gStyle->SetOptStat(1111);
                hResults.at(iTH1).first->Draw();
                cDrawResult->SaveAs((fFolder+fName+TString(Form("_%i.pdf",iTH1))).Data());
                gStyle->SetOptStat(0);
                gROOT->SetBatch(kFALSE);
            }
            iTH1++;
            if ( (iPT != jPT) && fDiag )    {
                hResults.push_back( std::pair<TH1F*,Int_t>( new TH1F(Form("VarHist_%i_%i",jPT,iPT),Form("VarHist_%i_%i",jPT,iPT),2000,-1,+1) , hStandard->GetBin( jPT+1, iPT+1 ) ) );
                for ( auto  hVariationSpectrum : uSystVariationSpectra_ ) {
                    hResults.at(iTH1).first->Fill(hVariationSpectrum->GetBinContent(iPT+1,jPT+1)-1);
                }
                iTH1++;
            }
        }
    }
    return hResults;
}


//_____________________________________________________________________________
THStack*
uBuildSystematicStack
( TH1F* hStandard, std::vector<TH1F*> hVariations, std::vector<Bool_t> fWhichToChoose = std::vector<Bool_t> () ) {
    auto uVariationsBinByBin = uBuildVariationBinByBin( hStandard, hVariations, "", "", fWhichToChoose );
    THStack    *hResults    =   new THStack ("uBuildSystematicStack","uBuildSystematicStack");
    TH1F       *hMean       =   (TH1F*)hStandard->Clone();
    hMean                   ->  SetFillColorAlpha(2,0.1);
    hMean                   ->  SetLineColorAlpha(2,1.0);
    hMean                   ->  SetLineWidth(1);
    TH1F       *hSTDV       =   (TH1F*)hStandard->Clone();
    hSTDV                   ->  SetFillColorAlpha(4,0.1);
    hSTDV                   ->  SetLineColorAlpha(4,1.0);
    hSTDV                   ->  SetLineWidth(1);
    auto    iBin    =   0;
    for ( auto uBinVariation : uVariationsBinByBin ) {
        iBin++;
        hMean   ->  SetBinContent   ( iBin, fabs(uBinVariation->GetMean()) );
        hMean   ->  SetBinError     ( iBin, 0   );
        hSTDV   ->  SetBinContent   ( iBin, uBinVariation->GetRMS()        );
        hSTDV   ->  SetBinError     ( iBin, 0   );
    }
    hMean           ->  Scale(100);
    hSTDV           ->  Scale(100);
    hResults        ->  Add(hMean);
    hResults        ->  Add(hSTDV);
    
    return  hResults;
}
std::vector<THStack*>
uBuildSystematicStack
( TH2F* hStandard, std::vector<TH2F*> hVariations, std::vector<Bool_t> fWhichToChoose = std::vector<Bool_t> ()  ) {
    std::vector<THStack*>   hResults;
    for ( Int_t iPT2D = 1; iPT2D <= hStandard->GetNbinsY(); iPT2D++ ) {
        TH1F * hTargetSTD = (TH1F*)(hStandard->ProjectionX(Form("%i",iPT2D),iPT2D,iPT2D))->Clone();
        std::vector<TH1F*>   hVariationSlice;
        for ( auto uSlice : hVariations )   {
            hVariationSlice.push_back((TH1F*)(uSlice->ProjectionX(Form("%i_2",iPT2D),iPT2D,iPT2D))->Clone());
        }
        hResults.push_back(uBuildSystematicStack( hTargetSTD, hVariationSlice, fWhichToChoose ));
    }
    return hResults;
}

//_____________________________________________________________________________
TH1F*
uBuildSystematicError
( TH1F* hStandard, std::vector<TH1F*> hVariations, TString fFolder = "", TString fName = "", std::vector<Bool_t> fWhichToChoose = std::vector<Bool_t> () ) {
    auto uVariationsBinByBin = uBuildVariationBinByBin( hStandard, hVariations, fFolder, fName, fWhichToChoose );
    TH1F   *hResults    =   (TH1F*)(hStandard->Clone());
    auto    iBin    =   0;
    for ( auto uBinVariation : uVariationsBinByBin ) {
        iBin++;
        hResults    ->  SetBinContent   ( iBin, fabs(uBinVariation->GetMean()) + uBinVariation->GetRMS() );
        hResults    ->  SetBinError     ( iBin, 0   );
    }
    return  hResults;
}
TH2F*
uBuildSystematicError
( TH2F* hStandard, std::vector<TH2F*> hVariations, TString fFolder = "", TString fName = "", std::vector<Bool_t> fWhichToChoose = std::vector<Bool_t> (), Bool_t fDiag = false ) {
    auto uVariationsBinByBin = uBuildVariationBinByBin( hStandard, hVariations, fFolder, fName, fWhichToChoose, fDiag );
    TH2F   *hResults    =   (TH2F*)(hStandard->Clone());
    hResults->Reset();
    for ( auto uBinVariation : uVariationsBinByBin ) {
        hResults    ->  SetBinContent   ( uBinVariation.second, fabs(uBinVariation.first->GetMean()) + uBinVariation.first->GetRMS() );
        hResults    ->  SetBinError     ( uBinVariation.second, 0   );
    }
    return  hResults;
}

//_____________________________________________________________________________
TCanvas*
uBuildSystAndStatStack
( TH1F* hStandard, THStack* hSystematical )  {
    TCanvas    *cResult =   new TCanvas();
    TGraphErrors   *gStatistical        =   new TGraphErrors();
    hSystematical->Draw("HIST F");
    gStatistical                        ->  SetFillColorAlpha(kGray,0.75);
    for ( Int_t iBin = 0; iBin < hStandard->GetNbinsX(); iBin++ )  {
        gStatistical->SetPoint      ( iBin, hStandard->GetBinCenter(iBin+1), 0. );
        gStatistical->SetPointError ( iBin, .5*hStandard->GetBinWidth(iBin+1), hStandard->GetBinError(iBin+1)/hStandard->GetBinContent(iBin+1) );
    }
    TMultiGraph    *cDrawAllGraphs      =   new TMultiGraph("cDrawAllGraphs","");
    //g1D_Syst_Err                              ->  SetFillColorAlpha(kGray+2,0.5);
    cDrawAllGraphs                      ->  Add         (gStatistical,      "AE2");
    cDrawAllGraphs                      ->  Add         (gStatistical,      "AE2");
    //cDrawAllGraphs                      ->  Add         (g1D_Syst_Err,      "AE2");
    cDrawAllGraphs->Draw("SAME ALP");
    cDrawAllGraphs->SetMinimum(0.);
    hSystematical->Draw("SAME HIST F");
    return  cResult;
}

//_____________________________________________________________________________ !TODO: !
TCanvas*
uBuildTotalSystematics
 ( std::vector<TH1F*> hSystematics ) {
    TCanvas    *cDrawResults    =   new TCanvas();
    TH1F*   hTotal  =   (TH1F*)hSystematics.at(0)->Clone();
    return  cDrawResults;
}

//_____________________________________________________________________________
void
uEvaluateRatioError
 ( TH1F* h1DStandard, std::vector<TH1F*> h1DVariations, TH2F* h2DStandard, std::vector<TH2F*> h2DVariations, TString fFolder = "", TH1F* h1DEfficiency = nullptr, TH2F* h2DEfficiency = nullptr ) {
    auto f1DRelevantVariations  =   uIsRelevantVariation(   h1DStandard,  h1DVariations,    (fFolder+TString("/plots/BarlowCheck/1D/")).Data(),"1D");
    auto f2DRelevantVariations  =   uIsRelevantVariation(   h2DStandard,  h2DVariations,    (fFolder+TString("/plots/BarlowCheck/2D/")).Data(),"2D",true);
    auto h1DFullSystematics     =   uBuildSystematicError(  h1DStandard,  h1DVariations,    (fFolder+TString("/plots/BinByBinCheck/1D/")).Data(),"1D",  f1DRelevantVariations);
    auto h2DFullSystematics     =   uBuildSystematicError(  h2DStandard,  h2DVariations,    (fFolder+TString("/plots/BinByBinCheck/2D/")).Data(),"2D",  f2DRelevantVariations,    true);
    //
    if ( h1DEfficiency ) h1DStandard->Divide(h1DEfficiency);
    if ( h2DEfficiency ) h2DStandard->Divide(h2DEfficiency);
    for ( Int_t iX = 0; iX < h1DStandard->GetNbinsX(); iX++ )   {
        h1DStandard->SetBinError(iX+1,h1DFullSystematics->GetBinContent(iX+1)*h1DStandard->GetBinContent(iX+1));
    }
    for ( Int_t iX = 0; iX < h2DStandard->GetNbinsX(); iX++ )   {
        for ( Int_t iY = 0; iY < h2DStandard->GetNbinsY(); iY++ )   {
            h2DStandard->SetBinError(iX+1,iY+1,h2DFullSystematics->GetBinContent(iX+1,iY+1)*h2DStandard->GetBinContent(iX+1,iY+1));
        }
    }
    //
    auto    iTer    =   0;
    std::vector<Float_t>   kSimpleRatio;
    std::vector<Float_t>   kSquareRatio;
    std::vector<Float_t>   k1Phi_Sng;
    std::vector<Float_t>   k2Phi_Sng;
    auto    k1DStdIntegralErr  =   0.;
    auto    k2DStdIntegralErr  =   0.;
    auto    k1DStdIntegral     =   h1DStandard->IntegralAndError(-1,10000,k1DStdIntegralErr,"width");
    auto    k2DStdIntegral     =   h2DStandard->IntegralAndError(-1,10000,-1,10000,k2DStdIntegralErr,"width");
            k1DStdIntegralErr /=   k1DStdIntegral;
            k2DStdIntegralErr /=   k2DStdIntegral;
    //
    for ( auto hVariation : h1DVariations )   {
        //
        if ( h1DEfficiency ) h1DVariations.at(iTer)->Divide(h1DEfficiency);
        if ( h2DEfficiency ) h2DVariations.at(iTer)->Divide(h2DEfficiency);
        auto    k1Dintegral     =   h1DVariations.at(iTer)->Integral("width");
        auto    k2Dintegral     =   h2DVariations.at(iTer)->Integral("width");
        //
        k1Phi_Sng.push_back(k1Dintegral/k1DStdIntegral -1);
        k2Phi_Sng.push_back(k2Dintegral/k2DStdIntegral -1);
        kSimpleRatio.push_back((k1DStdIntegral*k2Dintegral)/(k1Dintegral*k2DStdIntegral)-1);
        kSquareRatio.push_back((k1DStdIntegral*k1DStdIntegral*k2Dintegral)/(k1Dintegral*k1Dintegral*k2DStdIntegral)-1);
        iTer++;
    }
    //
    TH1F       *hSimpleRatioError   =   uBuildTH1F(kSimpleRatio,2000,0,-0.5,0.5);
    auto        fSimpleRatioError   =   0.;
    TH1F       *hSquareRatioError   =   uBuildTH1F(kSquareRatio,2000,0,-0.5,0.5);
    auto        fSquareRatioError   =   0.;
    TH1F       *h1Phi_Sng_Error     =   uBuildTH1F(k1Phi_Sng,2000,0,-0.5,0.5);
    auto        f1Phi_Sng_Error     =   0.;
    TH1F       *h2Phi_Sng_Error     =   uBuildTH1F(k2Phi_Sng,2000,0,-0.5,0.5);
    auto        f2Phi_Sng_Error     =   0.;
    fSimpleRatioError   +=  fabs(hSimpleRatioError->GetMean());
    fSimpleRatioError   +=  hSimpleRatioError->GetRMS();
    fSquareRatioError   +=  fabs(hSquareRatioError->GetMean());
    fSquareRatioError   +=  hSquareRatioError->GetRMS();
    f1Phi_Sng_Error     +=  fabs(h1Phi_Sng_Error->GetMean());
    f1Phi_Sng_Error     +=  h1Phi_Sng_Error->GetRMS();
    f2Phi_Sng_Error     +=  fabs(h2Phi_Sng_Error->GetMean());
    f2Phi_Sng_Error     +=  h2Phi_Sng_Error->GetRMS();
    //
    TH1F   *hStandard   =   new TH1F("hStandard",   "", 2,  0,  2);
    hStandard->GetYaxis()->SetTitle("Systematic uncertainty (%)");
    hStandard->GetXaxis()->SetNdivisions(2);
    hStandard->GetXaxis()->SetBinLabel(hStandard->GetXaxis()->FindBin(0.5),"#LT Y_{2#phi} #GT / #LT Y_{1#phi} #GT");
    hStandard->GetXaxis()->SetBinLabel(hStandard->GetXaxis()->FindBin(1.5),"#LT Y_{2#phi} #GT / #LT Y_{1#phi} #GT^{2}");
    hStandard->GetXaxis()->LabelsOption("h");
    hStandard->SetMarkerColor(kColors[3]);
    hStandard->SetLineWidth(3);
    hStandard->SetMarkerStyle(kMarkers[3]);
    hStandard->SetBinContent(1,fSimpleRatioError);
    hStandard->SetBinContent(2,fSquareRatioError);
    hStandard->SetBinError  (1,0);
    hStandard->SetBinError  (2,0);
    hStandard->Scale(100);
    TH1F   *hStandard_Sng   =   new TH1F("hStandard_Sng",   "", 2,  0,  2);
    hStandard_Sng->GetYaxis()->SetTitle("Systematic uncertainty (%)");
    hStandard_Sng->GetXaxis()->SetNdivisions(2);
    hStandard_Sng->GetXaxis()->SetBinLabel(hStandard->GetXaxis()->FindBin(0.5),"#LT Y_{#phi} #GT");
    hStandard_Sng->GetXaxis()->SetBinLabel(hStandard->GetXaxis()->FindBin(1.5),"#LT Y_{#phi#phi} #GT");
    hStandard_Sng->GetXaxis()->LabelsOption("h");
    hStandard_Sng->SetMarkerColor(kColors[3]);
    hStandard_Sng->SetLineWidth(3);
    hStandard_Sng->SetMarkerStyle(kMarkers[3]);
    hStandard_Sng->SetBinContent(1,f1Phi_Sng_Error);
    hStandard_Sng->SetBinContent(2,f2Phi_Sng_Error);
    hStandard_Sng->SetBinError  (1,0);
    hStandard_Sng->SetBinError  (2,0);
    hStandard_Sng->Scale(100);
    TH1F   *hSimple     =   new TH1F("hLinear",     "", 2,  0,  2);
    hSimple->SetMarkerColor(kColors[2]);
    hSimple->SetLineWidth(3);
    hSimple->SetMarkerStyle(kMarkers[4]);
    hSimple->SetBinContent(1,k1DStdIntegralErr+k2DStdIntegralErr);
    hSimple->SetBinContent(2,2*k1DStdIntegralErr+k2DStdIntegralErr);
    hSimple->SetBinError  (1,0);
    hSimple->SetBinError  (2,0);
    hSimple->Scale(100);
    TH1F   *hSimple_Sng     =   new TH1F("hLinear_Sng",     "", 2,  0,  2);
    hSimple_Sng->SetMarkerColor(kColors[2]);
    hSimple_Sng->SetLineWidth(3);
    hSimple_Sng->SetMarkerStyle(kMarkers[4]);
    hSimple_Sng->SetBinContent(1,k1DStdIntegralErr);
    hSimple_Sng->SetBinContent(2,k2DStdIntegralErr);
    hSimple_Sng->SetBinError  (1,0);
    hSimple_Sng->SetBinError  (2,0);
    hSimple_Sng->Scale(100);
    TH1F   *hSquare     =   new TH1F("hSquare",     "", 2,  0,  2);
    hSquare->SetMarkerColor(kColors[1]);
    hSquare->SetLineWidth(3);
    hSquare->SetMarkerStyle(kMarkers[5]);
    hSquare->SetBinContent(1,SquareSum( {k1DStdIntegralErr,k2DStdIntegralErr} ));
    hSquare->SetBinContent(2,SquareSum( {2*k1DStdIntegralErr,k2DStdIntegralErr} ));
    hSquare->SetBinError  (1,0);
    hSquare->SetBinError  (2,0);
    hSquare->Scale(100);
    //
    auto fMaximum  = 100 * max ( 2*k1DStdIntegralErr+k2DStdIntegralErr, max (fSimpleRatioError, fSquareRatioError ) );
    hStandard->SetMaximum(fMaximum*1.3);
    fMaximum  = 100 * max ( k2DStdIntegralErr, f2Phi_Sng_Error );
    hStandard_Sng->SetMaximum(fMaximum*1.3);
    //
    gROOT->SetBatch();
    TCanvas    *c1 = new TCanvas("Ratio");
    //
    TLegend    *lLegend = new TLegend(0.18,0.82,0.33,0.72);
    lLegend     ->  AddEntry( hStandard,    "Ratio err.", "P" );
    lLegend     ->  AddEntry( hSimple,      "Linear err.", "P" );
    lLegend     ->  AddEntry( hSquare,      "Square err.", "P" );
    //
    hStandard->Draw("][ EP MIN0");
    hSimple->Draw("SAME EP ][");
    hSquare->Draw("SAME EP ][");
    lLegend->Draw("same");
    //
    c1->SaveAs((fFolder+TString("/plots/Ratio.pdf")).Data());
    delete c1;
    delete lLegend;
    //
                c1      = new TCanvas("Ratio");
    //
                lLegend = new TLegend(0.18,0.82,0.33,0.72);
    lLegend     ->  AddEntry( hStandard,    "Yield err.", "P" );
    lLegend     ->  AddEntry( hSimple,      "Integral err.", "P" );
    //
    hStandard_Sng->Draw("][ EP MIN0");
    hSimple_Sng->Draw("SAME EP ][");
    lLegend->Draw("same");
    //
    c1->SaveAs((fFolder+TString("/plots/Ratio_Sng.pdf")).Data());
    delete c1;
    //
    TCanvas    *c2 = new TCanvas("Ratio_Simple");
    //
    hSimpleRatioError->Draw();
    //
    c2->SaveAs((fFolder+TString("/plots/Ratio_Simple.pdf")).Data());
    delete c2;
    //
    TCanvas    *c3 = new TCanvas("Ratio_Square");
    //
    hSquareRatioError->Draw();
    //
    c3->SaveAs((fFolder+TString("/plots/Ratio_Square.pdf")).Data());
    delete c3;
    //
    gROOT->SetBatch(kFALSE);
    //
    TFile      *fOutput =   new TFile   (Form("%s%s",fFolder.Data(),"/RT_Systematic.root"),"recreate");
    hStandard->Scale(0.01);
    hStandard->Write();
    fOutput->Close();
}

TCanvas*                fMultipleError      ( TGraphAsymmErrors *fGStat, TGraphAsymmErrors *fGSyst, TGraphAsymmErrors **fGVart, int initial, int stop, std::vector<TString> fNames )  {
    TCanvas        *cDrawComparison     =   new TCanvas("cDrawComparison","",1000,1000);
    gStyle                              ->  SetOptStat(0);
    //
    TMultiGraph    *cDrawAllGraphs      =   new TMultiGraph("cDrawAllGraphs","");
    //
    auto    kNColumms = 1+(int)((stop-initial)/(3.));
    TLegend        *cComparisonLegend   =   new TLegend(0.11,0.75,min(0.89,0.11+0.05*kNColumms),0.89);
    cComparisonLegend                   ->  SetLineColorAlpha(1,0.);
    cComparisonLegend                   ->  SetFillColorAlpha(1,0.);
    cComparisonLegend                   ->  SetNColumns(kNColumms);
    //
    fGStat                              ->  SetFillColorAlpha(kGray,0.75);
    fGSyst                              ->  SetFillColorAlpha(kGray+2,0.5);
    //
    cDrawAllGraphs                      ->  Add         (fGStat,      "AE2");
    cComparisonLegend                   ->  AddEntry    (fGStat,      "Stat. Err.",   "F");
    cDrawAllGraphs                      ->  Add         (fGSyst,      "AE2");
    cComparisonLegend                   ->  AddEntry    (fGSyst,      "Syst. Err.",   "F");
    for ( Int_t iTer = initial; iTer <= stop; iTer++ )    {
        if ( fNames.size() < iTer )  hTitle  =   Form("_%i",iTer);
        else            hTitle  =   fNames.at(iTer-1).Data();
        fGVart[iTer]                    ->  SetMarkerStyle(kRainbowMarker[2*(int)(iTer/6.)]);
        fGVart[iTer]                    ->  SetMarkerColor(fGetRainbowColor(iTer*2));
        fGVart[iTer]                    ->  SetLineColor(fGetRainbowColor(iTer*2));
        cDrawAllGraphs                  ->  Add         (fGVart[iTer],    "EPL");
        cComparisonLegend               ->  AddEntry    (fGVart[iTer],    hTitle,         "EP");
    }
    //
    cDrawAllGraphs                      ->  Draw("AC");
    cComparisonLegend                   ->  Draw("SAME");
    return                              cDrawComparison;
}
//
TH1F*  fMakeMeTH1F                          ( TGraphAsymmErrors   *fToBeTransformed )  {
    //  Checking the consistency of TH1*
    Int_t       fNPoints    =   fToBeTransformed ->  GetN();
    Double_t   *fBinArray   =   new Double_t    [fNPoints+1];
    for ( Int_t iFit = 0; iFit < fNPoints; iFit++ ) {
        auto    fXValue         =   ( fToBeTransformed ->  GetPointX(iFit) );
        auto    fXErrorELow     =   ( fToBeTransformed ->  GetErrorXlow(iFit) );
        fBinArray[iFit]         =   fXValue-fXErrorELow;
    }
    auto    fXValue         =   ( fToBeTransformed ->  GetPointX(fNPoints-1) );
    auto    fXErrorELow     =   ( fToBeTransformed ->  GetErrorXhigh(fNPoints-1) );
    fBinArray[fNPoints]     =   fXValue+fXErrorELow;
    //
    TH1F   *fResult =   new TH1F(Form("%s_toTH1F",fToBeTransformed->GetName()),fToBeTransformed->GetTitle(),fNPoints,fBinArray);
    //
    for ( Int_t iFit = 0; iFit < fNPoints; iFit++ ) {
        auto    fYValue         =   ( fToBeTransformed ->  GetPointY(iFit) );
        auto    fYErrorEHigh    =   ( fToBeTransformed ->  GetErrorYhigh(iFit) );
        auto    fYErrorELow     =   ( fToBeTransformed ->  GetErrorYlow(iFit) );
        fResult ->  SetBinContent   (iFit+1,fYValue);
        fResult ->  SetBinError     (iFit+1,max(fYErrorEHigh,fYErrorELow));
    }
    //
    return  fResult;
}
//
TCanvas*            fPublishSpectrum                ( TGraphAsymmErrors *gStat, TGraphAsymmErrors *gSyst, Int_t kColor = 2, Int_t kMarker = 8, TString fOption = "" )  {
    TCanvas        *fResult         =   new TCanvas();
    //fResult         ->  SetLogy();
    
    //  Prepping Stat Graph
    gStat           ->  SetMarkerStyle(kMarker);
    gStat           ->  SetMarkerColor(kColor);
    gStat           ->  SetFillColorAlpha(kColor,0.0);
    gStat           ->  SetLineColor(kColor);
    
    //  Prepping Syst Graph
    gSyst           ->  SetMarkerStyle(kMarker);
    gSyst           ->  SetMarkerColor(kColor);
    gSyst           ->  SetLineColor(kColor);
    gSyst           ->  SetFillColorAlpha(kColor,0.4);
    
    //  Prepping the multigraph
    TMultiGraph    *fPlotResult     =   new TMultiGraph();
    fPlotResult     ->  Add(gStat,"PE5");
    fPlotResult     ->  Add(gSyst,"PE2");
    fPlotResult     ->  Draw("A");
    
    //  Prepping Legend
    TLegend        *fLegendResult   =   new TLegend(0.7,0.7,0.9,0.9);
    fLegendResult   ->  SetFillStyle(0);
    fLegendResult   ->  SetBorderSize(0);
    fLegendResult   ->  AddEntry(gStat,"","EPF");
    fLegendResult   ->  AddEntry(gSyst,"","PF");
    fLegendResult   ->  Draw("same");
    
    fResult         ->  SaveAs("tt.pdf");
    fResult         ->  Write();
    return      fResult;
}
TGraphMultiErrors*  fGraphMultiErrors       ( TH1F *hStatistics, TH1F *hSystematics )  {
    auto    fNPoints    =   hStatistics->GetNbinsX();
    TGraphMultiErrors*  fResult =   new TGraphMultiErrors(fNPoints,3);
    fResult->SetMarkerColor(38);
    fResult->SetMarkerStyle(22);
    fResult->GetAttLine(0)->SetLineColorAlpha(kRed,1.);
    fResult->GetAttFill(0)->SetFillColorAlpha(kRed,0.);
    fResult->GetAttLine(1)->SetLineColorAlpha(kBlue,1.);
    fResult->GetAttFill(1)->SetFillColorAlpha(kBlue,0.3);
    for ( Int_t iPnt = 0; iPnt < fNPoints; iPnt++ ) {
        auto    fYValue =   hStatistics ->  GetBinContent(iPnt+1);
        auto    fCenter =   hStatistics ->  GetBinCenter(iPnt+1);
        auto    fXError =   hStatistics ->  GetBinLowEdge(iPnt+2) - fCenter;
        auto    fEStat  =   hStatistics ->  GetBinError(iPnt+1);
        auto    fESyst  =   hSystematics->  GetBinError(iPnt+1);
        fResult    ->  SetPoint        ( iPnt,   fCenter,    fYValue );
        fResult    ->  SetPointEX      ( iPnt,   fXError,    fXError     );
        fResult    ->  SetPointEY      ( iPnt,   0,          SquareSum({fESyst,fEStat}),     SquareSum({fESyst,fEStat}) );
        fResult    ->  SetPointEY      ( iPnt,   1,          fESyst,     fESyst );
        fResult    ->  SetPointEY      ( iPnt,   2,          fEStat,     fEStat );
    }
    return  fResult;
}
//
TCanvas*            fRatioPlot                  ( TH1   *fToBeCorrected, TH1    *fAccepted,  TH1    *fTotal,    Double_t fScale = 1. )  {
    return nullptr;
}
//
Double_t
uHistoIntegralAndError
( std::vector<TH1F*> fInputData, Double_t &fInputError ) {
    Double_t*   fResult     =   new Double_t[2];
    fResult[0] = 0;
    fResult[1] = 0;
    //
    auto    iHisto  =   0;
    for ( auto uTarget : fInputData )   {
        iHisto++;
        auto fPrevEr    =   fResult[1];
        auto fError     =   0.;
        auto fIntegral  =   uTarget->IntegralAndError(iHisto+1,10000,fError,"width");
        //
        auto fWidth     =   uTarget->GetBinWidth(iHisto);
        //
        fResult[0]     +=   fIntegral*fWidth;
        fResult[1]      =   SquareSum( { fError*fWidth, fPrevEr } );
    }
    iHisto  =   0;
    for ( auto uTarget : fInputData )   {
        iHisto++;
        auto fPrevEr    =   fResult[1];
        auto fError     =   0.;
        auto fIntegral  =   uTarget->IntegralAndError(iHisto,iHisto,fError,"width");
        //
        auto fWidth     =   uTarget->GetBinWidth(iHisto);
        //
        fResult[0]     +=   fIntegral*fWidth/(2.);
        fResult[1]      =   SquareSum( { fError*fWidth/(2.), fPrevEr } );
    }
    fInputError = fResult[1];
    return  fResult[0];
}

void
uCompareResultsTH1F
( std::vector<TString> fFileNames, std::vector<TString> fHistoName, std::vector<TString> fLegendEntries ) {
    TCanvas*    cDrawComparison =   new TCanvas("cDrawComparison","cDrawComparison",800,1100);
    cDrawComparison->SetTopMargin(0.3);
    TLegend*    lDrawLegend     =   new TLegend(0.1,0.75,0.9,0.92);
    lDrawLegend->SetNColumns(5);
    auto iTer = 0;
    for ( auto kFileName : fFileNames ) {
        auto    kHisto  =   (TH1F*)((new TFile((kFileName).Data()))->Get(fHistoName.at(iTer)));
        kHisto->SetMarkerStyle(fGetMarker( (int)((iTer)/8), 1 ));
        kHisto->SetMarkerColor(fGetColor(iTer));
        kHisto->SetLineColor(fGetColor(iTer));
        kHisto->Draw("same");
        lDrawLegend->AddEntry(kHisto,fLegendEntries.at(iTer).Data());
        iTer++;
    }
    lDrawLegend->Draw("same");
}
void
uCompareResultsTH1F
( std::vector<TString> fFileNames, TString fHistoName, std::vector<TString> fLegendEntries ) {
    std::vector<TString> fUtility;
    for ( auto kFile : fFileNames ) fUtility.push_back(fHistoName);
    uCompareResultsTH1F( fFileNames, fUtility, fLegendEntries );
}
//

#endif
