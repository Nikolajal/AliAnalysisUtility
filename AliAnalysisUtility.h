// Global Values and constants file
//
// !TODO: 1. - (_!_) To be implemented (...)
// !TODO: 2. - (_!_) Make the Markerstyle in a kArray
// !TODO: 3. x (_!_) Low number clearer in colors and markers
//
#ifndef ALIANALYSISUTILITY_H
#define ALIANALYSISUTILITY_H
//
//>>    Functions File
#include "AliAnalysisUtility_Functions.h"
//
// C++
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <chrono>
//
// ROOT
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TBenchmark.h"
//
// RooFit
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooGlobalFunc.h"
#include "RooMsgService.h"
#include "RooClassFactory.h"
//
// RooFitFunction
#include "RooChebychev.h"
#include "RooArgusBG.h"
#include "RooBreitWigner.h"
#include "RooExponential.h"
#include "RooVoigtian.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooUniform.h"
//
using namespace std;
using namespace RooFit;
//
//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//
//
// Random Generator
TRandom            *uRandomGen                      =   new TRandom();
// Benchmark
TBenchmark         *uBenchmark                      =   new TBenchmark();
// TLatex
TLatex             *uLatex                          =   new TLatex();

// TO BE CLEANED UP
auto iBuilderTH1FCounter = 0;
// Title and Name for histograms
auto                hName                           =   "Name";
auto                hTitle                          =   "Title";
// General Parameters
const Bool_t        kFitScarseHisto                 =   kTRUE;      //  Skip the fit of histograms that are below the threshold set below
const Float_t       kScarseHistoDef                 =   0.;         //  % of entries w.r.t. total bin number
const Int_t         kScarseHistoMin                 =   1000.;      //  N of entries
const Float_t       kLooseErrors                    =   3.;         //  Multiplication factor for the Error looosening
const Double_t      kMaximumError                   =   2.5;        //  Maximum percentage error to accept fit result, will retry if higher
const Int_t         kRainbowColor[]                 =   {kRed+1,kOrange-1,kYellow+1,kSpring-1,kGreen+1,kTeal-1,kCyan+1,kAzure-1,kBlue+1,kViolet-1,kMagenta+1,kPink-1}; // Up to 12 spectra in Rainbow colours
const Int_t         kRainbowMarker[]                =   {20,21,22,23,24,25}; // Up to 12 Marker in Rainbow colours
// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};
// Canvas Style
bool kStyleSet = false;
void SetStyle(Bool_t graypalette = false) {
    if ( kStyleSet ) return;
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);


}
//
//>> Physics is by defualt in GeV, we define some constants to evaluate results in other units
//>> When inserting the values into the system use *MeV, when extracting from the system use /MeV
//>> To change this, use the constants below
auto const          KeV                             =   1e-6;
auto const          MeV                             =   1e-3;
auto const          GeV                             =   1;
auto const          TeV                             =   1e+3;
//
Double_t
SquareSum
 ( std::initializer_list<Double_t> list )  {
    Double_t    fResult =   0;
    for ( auto element : list )  {
        fResult +=  element*element;
    }
    return  TMath::Sqrt(fResult);
}
//
//--------------------------------//
//      BENCHMARK UTILITIES       //
//--------------------------------//
//
TString                 kMSG_PrintTimer             =   "\r[INFO] Event # %7.f %s | %3.1f %% | %7.2f %s events/s | Time: %02.0f:%02.0f:%02.0f | ETA: %02.0f:%02.0f:%02.0f";
//
void
fStartTimer
 ( TString fTimerName )    {
    uBenchmark->Start(fTimerName.Data());
    printf("[INFO] Starting %s \n", fTimerName.Data());
    fflush(stdout);
}
//
void
fStopTimer
 ( TString fTimerName )     {
    uBenchmark->Stop(fTimerName.Data());
    cout << endl;
    printf("[INFO] Stopping %s \n", fTimerName.Data());
    Float_t fElapsedS   = (float)(uBenchmark->GetRealTime(fTimerName.Data()));
    Float_t fElapsedM   = (Int_t)(fElapsedS/60.);
    printf("[INFO] It took %02.0f:%02.0f \n",   fElapsedM,  fElapsedS - 60.*fElapsedM);
    fflush(stdout);
}
//
void
fPrintLoopTimer
 ( TString fTimerName, Int_t iEvent, Int_t nEntries, Int_t iPrintInterval )   {
    if ( iEvent%iPrintInterval != 0 || iEvent == 0 ) return;
    
    // Suffix for events
    TString     fSuffix =   "";
    Int_t       fSfxCor =   iPrintInterval;
    if ( iPrintInterval/1e3 < 1e3         && iPrintInterval/1e3 >= 1 )       {
        fSuffix =   "k";
        fSfxCor =   (int)(iPrintInterval/1e3) + iPrintInterval%(int)1e3;
    }
    if ( iPrintInterval/1e6 < 1e6         && iPrintInterval/1e6 >= 1 )       {
        fSuffix =   "mln";
        fSfxCor =   (int)(iPrintInterval/1e6) + iPrintInterval%(int)1e6;
    }
    if ( iPrintInterval/1e9 < 1e9         && iPrintInterval/1e9 >= 1 )       {
        fSuffix =   "mld";
        fSfxCor =   (int)(iPrintInterval/1e9) + iPrintInterval%(int)1e9;
    }
    
    // Stopping timer
    uBenchmark->Stop(fTimerName.Data());

    //
    //  Elapsed Time
    Float_t fRealElapsedSec =   (float)(uBenchmark->GetRealTime(fTimerName.Data()));
    Float_t fRealElapsedMin =   (Int_t)(fRealElapsedSec/60.);
    Float_t fRealElapsedHor =   (Int_t)(fRealElapsedSec/3600.);
    //
    Float_t fShowElapsedSec =   fRealElapsedSec - fRealElapsedMin*60.;
    Float_t fShowElapsedMin =   fRealElapsedMin - fRealElapsedHor*60.;
    Float_t fShowElapsedHor =   fRealElapsedHor;
    
    //
    //  Event utilities
    Float_t fProcessedFrac  =   (float)iEvent/((float)nEntries);
    Float_t fShowPrintEvnt  =   (float)iEvent*(float)fSfxCor/((float)iPrintInterval);   //TODO: Clean-Up
    Float_t fShowSpeedEvnt  =   fShowPrintEvnt/fRealElapsedSec;
    
    //
    //  ETA
    Float_t fRealEstimatSec =   fRealElapsedSec/fProcessedFrac - fRealElapsedSec;
    Float_t fRealEstimatMin =   (Int_t)(fRealEstimatSec/60.);
    Float_t fRealEstimatHor =   (Int_t)(fRealEstimatSec/3600.);
    //
    Float_t fShowEstimatSec =   fRealEstimatSec - fRealEstimatMin*60.;
    Float_t fShowEstimatMin =   fRealEstimatMin - fRealEstimatHor*60.;
    Float_t fShowEstimatHor =   fRealEstimatHor;
    
    
    // Resuming timer
    uBenchmark->Start(fTimerName.Data());
    
    // Printing
    cout << "\33[2K" << flush;
    cout << Form(kMSG_PrintTimer.Data(),  fShowPrintEvnt,  fSuffix.Data(), 100.*fProcessedFrac, fShowSpeedEvnt,  fSuffix.Data(), fShowElapsedHor, fShowElapsedMin,  fShowElapsedSec, fShowEstimatHor, fShowEstimatMin, fShowEstimatSec) << flush;
}
//
//------------------------------//
//    HISTOGRAM UTILITIES       //
//------------------------------//
//
//  --  --  Generation Utilities  --  --  //
//
TH1F*
uBuildTH1F
 ( std::vector<Float_t> fInputData, Int_t fNofBins = -1, Float_t fOffset = 0., Float_t fLowBound = 0, Float_t fHigBound = 0 )    {
    auto    fMaxValue   =   *std::max_element(std::begin(fInputData),std::end(fInputData));
    auto    fMinValue   =   *std::min_element(std::begin(fInputData),std::end(fInputData));
    auto    fSizeOfAr   =   fInputData.size();
    if ( fLowBound == fHigBound )   {
        fLowBound   =   fMinValue - 0.2*(fMaxValue-fMinValue) + fOffset;
        fHigBound   =   fMaxValue + 0.2*(fMaxValue-fMinValue) + fOffset;
    }
    if ( fNofBins <= 0 )    {
        fNofBins = 12;
        if      ( fSizeOfAr >= 1.e2 )   fNofBins = (int)(fSizeOfAr/3.) + 2;
        if      ( fSizeOfAr >= 1.e3 )   fNofBins = (int)(fSizeOfAr/5.) + 2;
        if      ( fSizeOfAr >= 1.e3 )   fNofBins = 216;
    }
    TH1F   *fBuiltTH1F  =   new TH1F(Form("TH1F_from_vector_%i",iBuilderTH1FCounter),Form("TH1F_from_vector_%i",iBuilderTH1FCounter),fNofBins,fLowBound,fHigBound);
    for     ( auto iValue : fInputData )    fBuiltTH1F->Fill( iValue + fOffset );
    iBuilderTH1FCounter++;
    return  fBuiltTH1F;
}
//
//  --  --  Customisation Utilities  --  --  //
//                                                  //  Relates to global variable kRainbowColor
void
uCleanOutsiders
 ( std::vector<Float_t> &fInputData )    {
    auto    fCheckAgain = true;
    while ( fCheckAgain )   {
        fCheckAgain = false;
        auto    MeanOfDistribution  =   0.;
        auto    STDVOfDistribution  =   0.;
        for ( auto iValue : fInputData ) {
            MeanOfDistribution  +=  iValue;
        }
        MeanOfDistribution      /=  fInputData.size();
        for ( auto iValue : fInputData ) {
            STDVOfDistribution  +=  (iValue - MeanOfDistribution)*(iValue - MeanOfDistribution);
        }
        STDVOfDistribution      /=  fInputData.size();
        STDVOfDistribution      =   TMath::Sqrt(STDVOfDistribution);
        for ( auto iValue = fInputData.begin(); iValue != fInputData.end(); ) {
            if ( fabs(*iValue - MeanOfDistribution) >= 8*STDVOfDistribution )   {
                fCheckAgain = true;
                iValue = fInputData.erase(iValue);
            }   else    {
                ++iValue;
            }
        }
    }
}
//                                                  //  Relates to global variable kRainbowColor
Int_t
fGetRainbowColor
 ( Int_t iTerator, Bool_t iSparse = false )    {
    if  ( iSparse ) iTerator    *=  2;
    return kRainbowColor[iTerator%12];
}
//                                                  //  Relates to global variable kRainbowColor
TGraphAsymmErrors*
fTH1_to_TGAsymmErrors
 ( TH1*  hTarget )    {
    TGraphAsymmErrors  *fResult                     =   new TGraphAsymmErrors();
    Int_t       fNBins  =   hTarget->GetNbinsX();
    for ( Int_t iBin = 1; iBin <= fNBins; iBin++ )  {
        auto    fXValue =   hTarget         ->GetBinCenter(iBin);
        auto    fYValue =   hTarget         ->GetBinContent(iBin);
        auto    fXError =   0.5*( hTarget->GetBinLowEdge(iBin+1) - hTarget->GetBinLowEdge(iBin) );
        auto    fYError =   hTarget         ->GetBinError(iBin);
        fResult         ->  SetPoint        (iBin-1,fXValue,fYValue);
        fResult         ->  SetPointError   (iBin-1,fXError,fXError,fYError,fYError);
    }
    fResult             ->SetMaximum(hTarget->GetMaximum());
    fResult             ->SetMinimum(hTarget->GetMinimum());
    return fResult;
}
//                                                  //  Relates to global variable kRainbowColor
TGraphAsymmErrors**
fTH2_to_TGAsymmErrors
 ( TH2*  hTarget )    {
    Int_t       fNBinsY =   hTarget->GetNbinsY();
    TGraphAsymmErrors **fResult                     =   new TGraphAsymmErrors*[fNBinsY];
    for ( Int_t iBin = 1; iBin <= fNBinsY; iBin++ )  {
        fResult[iBin-1] =   fTH1_to_TGAsymmErrors(hTarget ->  ProjectionX(Form("%i",iBin),iBin,iBin));
    }
    return fResult;
}
//
//_____________________________________________________________________________
//                                                  //  !TODO: Make the Markerstyle in a kArray
//                                                  //  Relates to global variable kRainbowColor
template < class Tclass >
void                    fSetRainbowStyle            ( std::vector<Tclass**> fInputObjectLits, Int_t fStartIteratorAt = 0 )  {
    TCanvas*    fResult     =   new TCanvas();
    Int_t       fIterator   =   fStartIteratorAt;
    for ( Tclass*&fObjectToPlot   :   fInputObjectLits )  {
        if ( fIterator == 12 ) fIterator = 0;
        fObjectToPlot->SetMarkerStyle(20+fIterator >= 31 ? 20+fIterator+1 : 20+fIterator );
        fObjectToPlot->SetMarkerColor(kRainbowColor[fIterator]);
        fObjectToPlot->SetLineWidth(2.);
        fObjectToPlot->SetLineStyle(9);
        fObjectToPlot->SetLineColor(kRainbowColor[fIterator]);
        fIterator++;
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
template < class Tclass >
void                    fSetUniformBinning          ( Tclass *fArrBin, Tclass fMinBin, Tclass fMaxBin, Int_t fNPoints )  {
    for (int i = 0; i <= fNPoints; i++ )
    {
        fArrBin[i] = fMinBin+(i)*(fMaxBin - fMinBin)/(static_cast<Float_t>(fNPoints));
    }
}
////_____________________________________________________________________________
void
uOffset
( TH1* hTarget, Double_t kOffset ){
    for ( Int_t iBin = 1; iBin <= hTarget->GetNbinsX(); iBin++ ) {
        auto fYValue = hTarget->GetBinContent(iBin);
        hTarget->SetBinContent(iBin, fYValue+kOffset );
    }
}
//
void
uSetHisto
( TH1* hTarget, TString fOption ){
    hTarget->SetTitle("");
    if ( fOption.Contains("1D") || fOption.Contains("12D") )   {
        if ( fOption.Contains("EFF") ) {
            // Y-Axis limits
            hTarget->Scale(100.);
            hTarget->SetMaximum(100.);
            hTarget->SetMinimum(0.);
            // Y-Axis title
            hTarget->GetYaxis()->SetTitle("Efficiency #times Acceptance (%)");
            // Marker style
            hTarget->SetMarkerStyle(markers[5]);
            // Colour scheme
            hTarget->SetMarkerColor(colors[2]);
            hTarget->SetLineColor(colors[2]);
            hTarget->SetFillColorAlpha(fillColors[2],0.33);
            //
            hTarget->SetOption("EP L");
            if ( fOption.Contains("EFF2") ) {
                hTarget->SetMarkerStyle(markers[9]);
                hTarget->SetMarkerColor(colors[1]);
                hTarget->SetLineColor(colors[1]);
                hTarget->SetFillColorAlpha(fillColors[1],0.33);
                hTarget->SetOption("EP L SAME");
            }
            if ( fOption.Contains("SL") ) {
                // Y-Axis title
                hTarget->GetYaxis()->SetTitle("Signal Loss (%)");
                uOffset(hTarget,-100);
                hTarget->SetMaximum(10.);
                hTarget->SetMinimum(-5.);
            }
        }
        if ( fOption.Contains("SPT") ) {
            
            /*
             // Preferred colors and markers
             const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
             const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
             const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};
             */
            
            hTarget->SetMarkerStyle(markers[3]);
            hTarget->SetMarkerColor(colors[2]);
            hTarget->SetMarkerSize(1);
            hTarget->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hTarget->GetYaxis()->SetTitleOffset(1.5);
            hTarget->GetYaxis()->SetTitle("1/N_{ev}dN^{2}/(dydp_{T})");
            if ( fOption.Contains("12D") )  {
                hTarget->GetYaxis()->SetTitle("1/N_{ev}dN^{3}/(dydp_{T,#phi_{1}}dp_{T,#phi_{2}})");
                hTarget->GetXaxis()->SetTitle("#it{p}_{T,#phi_{1}} (GeV/#it{c})");
            }
            if ( fOption.Contains("STAT") ) {
                hTarget->SetLineColor(colors[2]);
                hTarget->SetFillColorAlpha(colors[2],0.33);
                hTarget->SetOption("PE2");
            }
            if ( fOption.Contains("SYST") ) {
                hTarget->SetLineColor(colors[3]);
                hTarget->SetFillColorAlpha(colors[3],0.);
                hTarget->SetOption("PE2");
            }
        }
    } else if ( fOption.Contains("2D") )   {
        
    } else if ( fOption.Contains("3D") )   {
        cout << " Buu 3D " << endl;
    } else  {
        cout << " CANT GUESS DIMENSION " << endl;
    }
}
//
TCanvas*
uPlotSpectrum
( TH1* hTarget, TH1* hTrSyst, TString fOption = "" ){
    //
    SetStyle();
    //
    TCanvas    *cDrawResult =   new TCanvas();
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    uSetHisto(hTarget,fOption + TString("SPT STAT"));
    uSetHisto(hTrSyst,fOption + TString("SPT SYST"));
    //
    TLegend    *lLegend     =   new TLegend(0.2,0.35,0.4,0.5);
    lLegend     ->  SetFillColorAlpha(0.,0.);
    lLegend     ->  AddEntry    (hTarget,"Data","P");
    lLegend     ->  AddEntry    (hTarget,"Stat","F");
    lLegend     ->  AddEntry    (hTrSyst,"Syst","F");
    //
    hTarget->Draw();
    hTrSyst->Draw("SAME E1");
    lLegend->Draw("SAME");
    //
    uLatex->SetTextFont(60);
    uLatex->SetTextSize(0.05);
    uLatex->DrawLatexNDC(0.20, 0.3,"ALICE");
    uLatex->SetTextFont(42);
    uLatex->SetTextSize(0.04);
    uLatex->DrawLatexNDC(0.20, 0.25,"pp #sqrt{#it{s}}= 7 TeV");
    uLatex->DrawLatexNDC(0.20, 0.2,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
    //
    return cDrawResult;
}
//
void
uSetHisto
( TGraphMultiErrors* hTarget, TString fOption = "" ){
    hTarget->SetTitle("");
    hTarget->GetYaxis()->SetTitle("1/N_{ev}dN/dy");
    //hTarget->GetXaxis()->SetNdivisions(2);
    //hTarget->GetXaxis()->SetBinLabel(hTarget->GetXaxis()->FindBin(1),"#LT Y_{1#phi} #GT");
    //hTarget->GetXaxis()->SetBinLabel(hTarget->GetXaxis()->FindBin(2),"#LT Y_{2#phi} #GT");
    //hTarget->GetXaxis()->LabelsOption("h");
    hTarget->GetXaxis()->SetLabelSize(0.08);
    hTarget->SetLineColorAlpha(0.,0.);
    hTarget->SetMarkerStyle(21);
    hTarget->SetMarkerColor(kRed);
    hTarget->GetAttLine(0)->SetLineColor(38);
    hTarget->GetAttLine(1)->SetLineColor(46);
    hTarget->GetAttFill(0)->SetFillColorAlpha(38,0.33);
    hTarget->GetAttFill(1)->SetFillColorAlpha(46,0.33);
}
//
//  --  --  Data Handling Utilities  --  --  //
//
//_____________________________________________________________________________
//
TGraphAsymmErrors      *fSumErrors                  ( TGraphAsymmErrors* gBasic, TGraphAsymmErrors* gAddition )    {
    //  Checking the consistency of TGraphs
    Int_t   fNPoints =   gBasic ->  GetN();
    if  ( fNPoints  != gAddition ->  GetN() )   {
        cout << "[ERROR] Systematics and Statistics do not have the same number of points! Skipping this one..." << endl;
        return nullptr;
    }
    //
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*gBasic);
    for ( Int_t iFit = 0; iFit < fNPoints; iFit++ ) {
        auto    fXErrBsicLow    =   ( gBasic ->  GetErrorXlow(iFit) );
        auto    fXErrBsicHigh   =   ( gBasic ->  GetErrorXhigh(iFit) );
        auto    fYErrBsicLow    =   ( gBasic ->  GetErrorYlow(iFit) );
        auto    fYErrBsicHigh   =   ( gBasic ->  GetErrorYhigh(iFit) );
        auto    fYErrAddtLow    =   ( gAddition ->  GetErrorYlow(iFit) );
        auto    fYErrAddtHigh   =   ( gAddition ->  GetErrorYhigh(iFit) );
        fResult ->  SetPointError(iFit,fXErrBsicLow,fXErrBsicHigh,sqrt(fYErrBsicLow*fYErrBsicLow+fYErrAddtLow*fYErrAddtLow),sqrt(fYErrBsicHigh*fYErrBsicHigh+fYErrAddtHigh*fYErrAddtHigh));
    }
    //
    return  fResult;
}
//
template < class Tclass >
Tclass                   *fSumErrors                  ( Tclass* gBasic, Tclass* gAddition )    {
    Tclass  *fResult =   new Tclass(*gBasic);
    for ( Int_t iBin = 0; iBin < gBasic->GetNbinsX(); iBin++ ) {
        fResult ->  SetBinError( iBin, SquareSum( { gBasic->GetBinError(iBin), gAddition->GetBinError(iBin) } ) );
    }
    //
    return  fResult;
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors      *fScaleWithError             ( TGraphAsymmErrors* gBasic, Double_t fScale, Double_t fScaleErrHigh = 0., Double_t fScaleErrLow = 0. )    {
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*gBasic);
    for ( Int_t iFit = 0; iFit < fResult->GetN(); iFit++ ) {
        auto    fYValue         =   ( gBasic ->  GetPointY(iFit) );
        auto    fYErrBsicLow    =   ( gBasic ->  GetErrorYlow(iFit) );
        auto    fYErrBsicHigh   =   ( gBasic ->  GetErrorYhigh(iFit) );
        fResult ->  SetPointY       ( iFit, fYValue*fScale);
        fResult ->  SetPointEYhigh  ( iFit, (fYValue*fScale)*sqrt(fYErrBsicHigh*fYErrBsicHigh/(fYValue*fYValue) + fScaleErrHigh*fScaleErrHigh/(fScale*fScale)));
        fResult ->  SetPointEYlow   ( iFit, (fYValue*fScale)*sqrt(fYErrBsicLow*fYErrBsicLow/(fYValue*fYValue) + fScaleErrLow*fScaleErrLow/(fScale*fScale)));
    }
    //
    return  fResult;
}
TH1F                   *fScaleWithError             ( TH1F* gBasic, Double_t fScale, Double_t fScaleError = 0. )    {
    TH1F  *fResult =   new TH1F(*gBasic);
    for ( Int_t iBin = 1; iBin <= fResult->GetNbinsX(); iBin++ ) {
        fResult->SetBinContent( iBin, gBasic->GetBinContent(iBin)/fScale );
        fResult->SetBinError( iBin, (gBasic->GetBinContent(iBin)/fScale)*SquareSum( { gBasic->GetBinError( iBin )/gBasic->GetBinContent(iBin), fScaleError/fScale } ) );
    }
    return  fResult;
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors      *uRandomizePoints            ( TGraphAsymmErrors* gStatic, TGraphAsymmErrors* gMoveable )    {
    //  Checking the consistency of TGraphs
    Int_t   fNPoints =   gStatic ->  GetN();
    if  ( fNPoints  != gMoveable ->  GetN() )
    {
        cout << "[ERROR] Systematics and Statistics do not have the same number of points! Skipping this one..." << endl;
        return nullptr;
    }
    //
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*gStatic);
    for ( Int_t iFit = 0; iFit < fNPoints; iFit++ ) {
        auto    fXValue         =   ( gStatic ->  GetPointX(iFit) );
        auto    fYValue         =   ( gStatic ->  GetPointY(iFit) );
        auto    fYErrMvblLow    =   ( gMoveable ->  GetErrorYlow(iFit) );
        auto    fYErrMvblHigh   =   ( gMoveable ->  GetErrorYhigh(iFit) );
        auto    fIsFluctLow     =   true;
        auto    fYNewValue      =   fYValue;
        ( uRandomGen  ->  Uniform (0.,1.) ) > 0.5? fIsFluctLow = true : fIsFluctLow = false;
        if ( fIsFluctLow )  {
            fYNewValue  -= fabs(uRandomGen  ->  Gaus(fYValue,fYErrMvblLow) - fYValue);
        }   else    {
            fYNewValue  += fabs(uRandomGen  ->  Gaus(fYValue,fYErrMvblHigh) - fYValue);
        }
        fResult->SetPoint(iFit,fXValue,fYNewValue);
    }
    return  fSumErrors(fResult,gMoveable);
}
//
template < class Tclass >
Tclass*                 uRandomizePoints            ( Tclass* gStatic, Tclass* gMoveable )    {
    Tclass* fResult =   (Tclass*)(gStatic->Clone());
    for ( Int_t iBin = 0; iBin < gStatic->GetNbinsX(); iBin++ ) {
        auto    fBinContent =   gStatic->GetBinContent(iBin+1);
        auto    fMoveError  =   gMoveable->GetBinError(iBin+1);
        auto fTest = uRandomGen -> Gaus(fBinContent,fMoveError);
        fResult->SetBinContent  ( iBin+1, fTest );
        fResult->SetBinError    ( iBin+1, fMoveError );
    }
    return  fSumErrors(fResult,gMoveable);
}
std::vector<TH1F*>      uRandomizePoints            ( std::vector<TH1F*>  gStatic, std::vector<TH1F*>  gMoveable )    {
    std::vector<TH1F*>  fResult;
    auto iTer = 0;
    for ( auto hStatic : gStatic ) {
        fResult.push_back(uRandomizePoints((TH1F*)hStatic->Clone(),(TH1F*)gMoveable.at(iTer)->Clone()));
        iTer++;
    }
    return  fResult;
}
std::vector<TH1F*>      uRandomizePointsSymm        ( std::vector<TH1F*>  gStatic, std::vector<TH1F*>  gMoveable )    {
    std::vector<TH1F*>  fResult;
    for ( auto hTarget : gStatic )  {
        fResult.push_back(new TH1F(*hTarget));
    }
    auto nHisto =   gStatic.size();
    auto nBins  =   gStatic.at(0)->GetNbinsX();
    for ( Int_t iTer = 1; iTer <= nHisto; iTer++ ) {
        for ( Int_t jTer = iTer; jTer <= nBins; jTer++ ) {
            auto    fBinContent =   gStatic.at(iTer-1)    ->GetBinContent (jTer);
            auto    fStatError  =   gStatic.at(iTer-1)    ->GetBinError   (jTer);
            auto    fMoveError  =   gMoveable.at(iTer-1)  ->GetBinError   (jTer);
            auto    fYNewValue  =   max( 0., uRandomGen -> Gaus(fBinContent,fMoveError) );
            auto    fYNewError  =   sqrt( fStatError*fStatError + fMoveError*fMoveError );
            fResult.at(iTer-1)->SetBinContent ( jTer, fYNewValue );
            fResult.at(iTer-1)->SetBinError   ( jTer, fYNewError );
            fResult.at(jTer-1)->SetBinContent ( iTer, fYNewValue );
            fResult.at(jTer-1)->SetBinError   ( iTer, fYNewError );
        }
    }
    return  fResult;
}
std::vector<TH1D*>      uRandomizePointsSymm        ( std::vector<TH1D*>  gStatic, std::vector<TH1D*>  gMoveable )    {
    std::vector<TH1D*>  fResult;
    for ( auto hTarget : gStatic )  {
        fResult.push_back(new TH1D(*hTarget));
    }
    auto nHisto =   gStatic.size();
    auto nBins  =   gStatic.at(0)->GetNbinsX();
    for ( Int_t iTer = 1; iTer <= nHisto; iTer++ ) {
        for ( Int_t jTer = iTer; jTer <= nBins; jTer++ ) {
            auto    fBinContent =   gStatic.at(iTer-1)    ->GetBinContent (jTer);
            auto    fStatError  =   gStatic.at(iTer-1)    ->GetBinError   (jTer);
            auto    fMoveError  =   gMoveable.at(iTer-1)  ->GetBinError   (jTer);
            auto    fYNewValue  =   max( 0., uRandomGen -> Gaus(fBinContent,fMoveError) );
            auto    fYNewError  =   sqrt( fStatError*fStatError + fMoveError*fMoveError );
            fResult.at(iTer-1)->SetBinContent ( jTer, fYNewValue );
            fResult.at(iTer-1)->SetBinError   ( jTer, fYNewError );
            fResult.at(jTer-1)->SetBinContent ( iTer, fYNewValue );
            fResult.at(jTer-1)->SetBinError   ( iTer, fYNewError );
        }
    }
    return  fResult;
}
//
//
//_____________________________________________________________________________
//
TH1F                    *fEfficiencycorrection       ( TH1   *fToBeCorrected, TH1    *fAccepted,  TH1   *fTotal,    Double_t fScale = 1. )  {
    TH1F   *fEfficiency =   (TH1F*)fAccepted->Clone();
    TH1F   *fResult     =   (TH1F*)fToBeCorrected->Clone();
    fEfficiency         ->  Divide(fAccepted,fTotal,1.,1.,"b");
    fResult             ->  Divide(fToBeCorrected,fEfficiency,fScale);
    return  fResult;
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors      *fEfficiencycorrection       ( TGraphAsymmErrors *fToBeCorrected, TH1    *fAccepted,  TH1    *fTotal,    Double_t fScale = 1. )  {
    //  Copying accordingly the TH1*
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*fToBeCorrected);
    //
    fResult->Divide(fAccepted,fTotal,"cl=0.683 b(1,1) mode");
    //
    Int_t   fNPoints =   fResult ->  GetN();
    for ( Int_t iFit = 0; iFit < fNPoints; iFit++ ) {
        auto    fYValueResult   =   ( fToBeCorrected ->  GetPointY(iFit) );
        auto    fYErrorRstHigh  =   ( fToBeCorrected ->  GetErrorYhigh(iFit) );
        auto    fYErrorRstLow   =   ( fToBeCorrected ->  GetErrorYlow(iFit) );
        auto    fYValueEffic    =   ( fResult ->  GetPointY(iFit) );
        auto    fYErrorEffHigh  =   ( fResult ->  GetErrorYhigh(iFit) );
        auto    fYErrorEffLow   =   ( fResult ->  GetErrorYlow(iFit) );
        fResult ->  SetPointY       (iFit,fScale*fYValueResult/fYValueEffic);
        fResult ->  SetPointEYlow   (iFit,(fScale*fYValueResult/fYValueEffic)*sqrt(fYErrorEffHigh*fYErrorEffHigh/(fYValueEffic*fYValueEffic) + fYErrorRstHigh*fYErrorRstHigh/(fYValueResult*fYValueResult) ));
        fResult ->  SetPointEYhigh  (iFit,(fScale*fYValueResult/fYValueEffic)*sqrt(fYErrorEffLow*fYErrorEffLow/(fYValueEffic*fYValueEffic) + fYErrorRstLow*fYErrorRstLow/(fYValueResult*fYValueResult) ));
    }
    //
    return  fResult;
}
//
//_____________________________________________________________________________
//                                                          // !TODO: To Be Implemented (...)
std::vector<TGraphAsymmErrors*> fEfficiencycorrection       ( TH2   *fToBeCorrected, TH2    *fAccepted,  TH2    *fTotal,    Double_t fScale = 1. )  {
    return std::vector<TGraphAsymmErrors*>();
}
//
std::vector<TH1F*> fEfficiencycorrection       ( TH2   *fToBeCorrected, TH1    *fAccepted,  TH1    *fTotal,    Double_t fScale = 1. )  {
    std::vector<TH1F*> fResult;
    if ( !fToBeCorrected )  { cout << "No fToBeCorrected" << endl; return fResult; }
    if ( !fAccepted )  { cout << "No fAccepted" << endl; return fResult; }
    if ( !fTotal )  { cout << "No fTotal" << endl; return fResult; }
    TH1F   *fEfficiency =   (TH1F*)fAccepted->Clone();
    fEfficiency         ->  Divide(fAccepted,fTotal,1.,1.,"b");
    for ( Int_t iHisto = 1; iHisto <= fToBeCorrected->GetNbinsY(); iHisto++ )    {
        auto    fConditional    =   fToBeCorrected->ProjectionY(Form("dd_%i",iHisto),iHisto,iHisto);
        TH1F*    fCorrCondit    =   fEfficiencycorrection( fConditional, fAccepted, fTotal, fScale );
        fResult.push_back( fScaleWithError( fCorrCondit, fEfficiency->GetBinContent(iHisto), fEfficiency->GetBinError(iHisto) ) );
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors*              fEfficiencycorrection       ( TGraphAsymmErrors    *fToBeCorrected,    TGraphAsymmErrors   *fEfficiency,     Double_t fScale = 1. )  {
    //  Checking the consistency of TGraphs
    Int_t   fNPoints =   fToBeCorrected ->  GetN();
    if  ( fNPoints  != fEfficiency ->  GetN() )   {
        cout << "[ERROR] Systematics and Statistics do not have the same number of points! Skipping this one..." << endl;
        return nullptr;
    }
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*fToBeCorrected);
    for ( Int_t iHisto = 1; iHisto <= fNPoints; iHisto++ )    {
        auto    fTargetY        =   fToBeCorrected->GetPointY(iHisto-1);
        auto    fEfficnY        =   fEfficiency->GetPointY(iHisto-1);
        auto    fTargetYlow     =   fToBeCorrected->GetErrorYlow(iHisto-1);
        auto    fEfficnYlow     =   fEfficiency->GetErrorYlow(iHisto-1);
        auto    fTargetYhigh    =   fToBeCorrected->GetErrorYlow(iHisto-1);
        auto    fEfficnYhigh    =   fEfficiency->GetErrorYlow(iHisto-1);
        fResult->SetPointY      (iHisto-1,  fTargetY*fEfficnY*fScale);
        fResult->SetPointEYlow  (iHisto-1,  sqrt(fTargetYlow*fTargetYlow+fEfficnYlow*fEfficnYlow));
        fResult->SetPointEYhigh (iHisto-1,  sqrt(fTargetYhigh*fTargetYhigh+fEfficnYhigh*fEfficnYhigh));
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
//  --  --  Specific purpose histogram generation  --  --  //
//
//_____________________________________________________________________________
//
TCanvas                *uPlotReferenceValue         ( TGraphAsymmErrors*  hMeasured,    Float_t fReference,   Float_t fRefError, TString fLabel = "PDG Value" )    {
    TCanvas    *fResult         =   new TCanvas();
    //
    auto        fNPoints        =   hMeasured->GetN();
    auto        fXLow           =   hMeasured->GetPointX(0);
    auto        fXLowErr        =   hMeasured->GetErrorXlow(0);
    auto        fXHig           =   hMeasured->GetPointX(fNPoints-1);
    auto        fXHigErr        =   hMeasured->GetErrorXhigh(fNPoints-1);
                hMeasured       ->  SetMaximum(max((float)(hMeasured->GetMaximum()+0.25*(hMeasured->GetMaximum()-hMeasured->GetMinimum())), fReference+fRefError*5));
                hMeasured       ->  SetMinimum(min((float)(hMeasured->GetMinimum()-0.25*(hMeasured->GetMaximum()-hMeasured->GetMinimum())), fReference-fRefError*5));
    //
                hMeasured       ->  Draw("APE");
    //
    TLegend    *fLegend         =   new TLegend(0.6,0.9,0.9,0.75);
                fLegend         ->  SetLineColorAlpha(1,0.);
                fLegend         ->  SetFillColorAlpha(1,0.);
                fLegend         ->  SetNColumns(2);
    //
    TH1F      **fPlotReference  =   new TH1F   *[7];
    for ( Int_t iTer = 0; iTer < 7; iTer ++ )   {
        fPlotReference[iTer]    =   new TH1F    (Form("%i",iTer),Form("%i",iTer),1,(fXLow-fXLowErr*1.5),(fXHig+fXHigErr*1.5));
        fPlotReference[iTer]    ->  SetBinContent(1, fReference + ( iTer -3 )*fRefError );
        fPlotReference[iTer]    ->  SetLineStyle(10);
        fPlotReference[iTer]    ->  SetLineWidth(3);
        fPlotReference[iTer]    ->  SetLineColorAlpha(920+(4-fabs(iTer-3)),0.75);
        fPlotReference[iTer]    ->  Draw("same ][");
    }
    fPlotReference[3]           ->  SetLineColorAlpha(2,0.75);
    fPlotReference[3]           ->  SetLineWidth(5);
    //
    fLegend                     ->  AddEntry( fPlotReference[3],    fLabel.Data(),  "L");
    fLegend                     ->  AddEntry( fPlotReference[2],    "#pm 1 #sigma", "L");
    fLegend                     ->  AddEntry( fPlotReference[1],    "#pm 2 #sigma", "L");
    fLegend                     ->  AddEntry( fPlotReference[0],    "#pm 3 #sigma", "L");
    fLegend                     ->  Draw("same");
    //
    return fResult;
}
TCanvas                *uPlotReferenceValue         ( TH1*  hMeasured,    Float_t fReference,   Float_t fRefError, TString fLabel = "PDG Value"  )    {
    return uPlotReferenceValue( fTH1_to_TGAsymmErrors(hMeasured), fReference, fRefError, fLabel );
}
//
//_____________________________________________________________________________
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
        cout << "[ERROR] You are trying to fit a null histogram! Skipping " << aTarget->GetName() << endl;
        return false;
    }
    if (  aTarget->GetEntries() == 0. ) {
        cout << "[ERROR] You are trying to fit an empty histogram! Skipping this one..." << endl;
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


//-------------------------------------//
//    GENERAL MULTIPURPOSE UTILITIES   //
//-------------------------------------//
//
//  --  --  BITMask Utilities  --  --  //
//
//_____________________________________________________________________________
//
bool                    fCheckMask                  ( Int_t fToBeChecked, Int_t iMaskCheck, Bool_t fOnlyThis = false )    {
    if ( fToBeChecked == 0  )   return false;
    if ( !fOnlyThis )   return  ((fToBeChecked & (int)(pow(2,iMaskCheck))) == (int)(pow(2,iMaskCheck)));
    else                return  ((fToBeChecked &~(int)(pow(2,iMaskCheck))) != 0 );
}
//
//_____________________________________________________________________________
//
//  --  --  Systematcis Utilities  --  --  //
//
//_____________________________________________________________________________
//
Double_t                uBarlowPar    ( Double_t  fStandard, Double_t fStdError, Double_t fVariatin, Double_t fVarError )  {
    auto    fSigmaStd       =   fStdError*fStdError;
    auto    fSigmaVar       =   fVarError*fVarError;
    auto    fSigmaDff       =   TMath::Sqrt( fabs( fSigmaStd - fSigmaVar ) );
    if      ( fSigmaDff    == 0 )  return false;
    auto    fParameter      =   ( fStandard - fVariatin )/fSigmaDff;
    return  fParameter;
}
bool                    fBarlowCheck    ( Double_t  fStandard, Double_t fStdError, Double_t fVariatin, Double_t fVarError )  {
    //return false;
    return  ( fabs ( uBarlowPar    ( fStandard, fStdError, fVariatin, fVarError ) )    <= 1 );
}
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
        hCheckHist  ->  SetLineColor(colors[2]);
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
    hStandard->SetMarkerColor(colors[3]);
    hStandard->SetLineWidth(3);
    hStandard->SetMarkerStyle(markers[3]);
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
    hStandard_Sng->SetMarkerColor(colors[3]);
    hStandard_Sng->SetLineWidth(3);
    hStandard_Sng->SetMarkerStyle(markers[3]);
    hStandard_Sng->SetBinContent(1,f1Phi_Sng_Error);
    hStandard_Sng->SetBinContent(2,f2Phi_Sng_Error);
    hStandard_Sng->SetBinError  (1,0);
    hStandard_Sng->SetBinError  (2,0);
    hStandard_Sng->Scale(100);
    TH1F   *hSimple     =   new TH1F("hLinear",     "", 2,  0,  2);
    hSimple->SetMarkerColor(colors[2]);
    hSimple->SetLineWidth(3);
    hSimple->SetMarkerStyle(markers[4]);
    hSimple->SetBinContent(1,k1DStdIntegralErr+k2DStdIntegralErr);
    hSimple->SetBinContent(2,2*k1DStdIntegralErr+k2DStdIntegralErr);
    hSimple->SetBinError  (1,0);
    hSimple->SetBinError  (2,0);
    hSimple->Scale(100);
    TH1F   *hSimple_Sng     =   new TH1F("hLinear_Sng",     "", 2,  0,  2);
    hSimple_Sng->SetMarkerColor(colors[2]);
    hSimple_Sng->SetLineWidth(3);
    hSimple_Sng->SetMarkerStyle(markers[4]);
    hSimple_Sng->SetBinContent(1,k1DStdIntegralErr);
    hSimple_Sng->SetBinContent(2,k2DStdIntegralErr);
    hSimple_Sng->SetBinError  (1,0);
    hSimple_Sng->SetBinError  (2,0);
    hSimple_Sng->Scale(100);
    TH1F   *hSquare     =   new TH1F("hSquare",     "", 2,  0,  2);
    hSquare->SetMarkerColor(colors[1]);
    hSquare->SetLineWidth(3);
    hSquare->SetMarkerStyle(markers[5]);
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

#endif
