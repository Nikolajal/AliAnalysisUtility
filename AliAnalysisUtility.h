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
#include "AAU_Graphs.h"
#include "AAU_Resolution.h"
#include "AAU_Efficiency.h"
#include "AAU_Extrapolation.h"
#include "AAU_Systematics.h"
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
        fGVart[iTer]                    ->  SetMarkerColor(uGetRainbowColor(iTer*2));
        fGVart[iTer]                    ->  SetLineColor(uGetRainbowColor(iTer*2));
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
        kHisto->SetMarkerStyle(uGetMarker( (int)((iTer)/8), 1 ));
        kHisto->SetMarkerColor(uGetColor(iTer));
        kHisto->SetLineColor(uGetColor(iTer));
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
