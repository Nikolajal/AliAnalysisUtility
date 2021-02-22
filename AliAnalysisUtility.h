// Global Values and constants file
// !TODO: All set!

#ifndef ALIANALYSISUTILITY_H
#define ALIANALYSISUTILITY_H

// C++
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <chrono>

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

// RooFitFunction
#include "RooChebychev.h"
#include "RooArgusBG.h"
#include "RooBreitWigner.h"
#include "RooExponential.h"
#include "RooVoigtian.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooUniform.h"

using namespace std;
using namespace RooFit;

//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//
//
// Random Generator
//
TRandom    *fRandomGen      =   new TRandom();
//
// Benchmark
//
TBenchmark *fBenchmark      =   new TBenchmark();
//
// Title and Name for histograms
//
auto        hName           =   "Name";
auto        hTitle          =   "Title";
//
// Title and Name for histograms
//
const Bool_t    kFitScarseHisto =   kTRUE;      //  Skip the fit of histograms that are below the threshold set below
const Float_t   kScarseHistoDef =   0.;         //  % of entries w.r.t. total bin number
const Int_t     kScarseHistoMin =   1000.;      //  N of entries
const Float_t   kLooseErrors    =   3.;         //  Multiplication factor for the Error looosening
const Double_t  kMaximumError   =   100.;       //  Maximum percentage error to accept fit result, will retry if higher
const Int_t     kRainbowColor[] =   {kRed+1,kOrange-1,kYellow+1,kSpring-1,kGreen+1,kTeal-1,kCyan+1,kAzure-1,kBlue+1,kViolet-1,kMagenta+1,kPink-1}; // Up to 12 spectra in Rainbow colours
//
//---------------------------------------//
//- Physics is by defualt in GeV, we    -//
//- define some constants to evaluate   -//
//- results in other units              -//
//---------------------------------------//
//
auto const  KeV             =   1e6;
auto const  MeV             =   1e3;
auto const  GeV             =   1;
auto const  TeV             =   1e-3;
//
//--------------------------------//
//      BENCHMARK UTILITIES       //
//--------------------------------//
//
TString     fMSG_PrintTimer =   "[INFO] Event # %4.f %s | %02.0f %% | %2.2f %s events/s | Time: %02.0f:%02.0f | ETA: %02.0f:%02.0f \n";
//
//_____________________________________________________________________________
//
void    fStartTimer             ( TString fTimerName )    {
    fBenchmark->Start(fTimerName.Data());
    printf("[INFO] Starting %s \n", fTimerName.Data());
    fflush(stdout);
}
//
//_____________________________________________________________________________
//
void    fStopTimer              ( TString fTimerName )     {
    fBenchmark->Stop(fTimerName.Data());
    printf("[INFO] Stopping %s \n", fTimerName.Data());
    Float_t fElapsedS   = (float)(fBenchmark->GetRealTime(fTimerName.Data()));
    Float_t fElapsedM   = (Int_t)(fElapsedS/60.);
    printf("[INFO] It took %02.0f:%02.0f \n",   fElapsedM,  fElapsedS - 60.*fElapsedM);
    fflush(stdout);
}
//
//_____________________________________________________________________________
//
void    fPrintLoopTimer         ( TString fTimerName, Int_t iEvent, Int_t nEntries, Int_t iPrintInterval )   {
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
    fBenchmark->Stop(fTimerName.Data());
    
    // Evaluating informations
    Float_t fFraction   =   (float)iEvent/((float)nEntries);
    Float_t fElapsedS   =   (float)(fBenchmark->GetRealTime(fTimerName.Data()));
    Float_t fElapsedM   =   (Int_t)(fElapsedS/60.);
    Float_t fPrintEvt   =   (float)iEvent*(float)fSfxCor/((float)iPrintInterval);
    Float_t fSpeedvsS   =   fPrintEvt/fElapsedS;
    Float_t fEta____S   =   (Int_t)(fElapsedS*((float)nEntries/((float)iEvent) -1));
    Float_t fEta____M   =   (Int_t)(fEta____S/60.);
    
    // Printing
    "[INFO] Event # %4.f %s | %02.0f %% | %1.2f %s events/s | Time: %02.0f:%02.0f | ETA: %02.0f:%02.0f \n";
    printf(fMSG_PrintTimer.Data(),  fPrintEvt,  fSuffix.Data(), 100.*fFraction, fSpeedvsS,  fSuffix.Data(), fElapsedM,  fElapsedS -60.*fElapsedM,  fEta____M,  fEta____S -60.*fEta____M);
    fflush(stdout);
    
    // Resuming timer
    fBenchmark->Start(fTimerName.Data());
}
//
//------------------------------//
//    HISTOGRAM UTILITIES       //
//------------------------------//
//
template < class Tclass >
bool                fIsWorthFitting             ( Tclass * aTarget )    {
    if ( !aTarget ) {
        cout << "[ERROR] You are trying to fit a null histogram! Skipping this one..." << endl;
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
//
template < class Tclass >
Tclass *            fLooseErrors                ( Tclass * aTarget, Float_t fLooseErrors = kLooseErrors )    {
    if ( fLooseErrors != kLooseErrors ) cout << "[INFO] Loosening factor for errors is different from the global set one, using the one provided in function call." << endl;
    if ( fLooseErrors < 1 )             cout << "[WARNING] Loosening factor for errors is less than 1, meaning your errors are less than what they really are. This might be causing problems to the fit function." << endl;
    if ( fLooseErrors == 1 )            cout << "[INFO] Loosening factor for errors is set to one, this is rendering this function useless. I'm returning a copy of the input." << endl;
    Tclass     *fResult =   new Tclass(*aTarget);
    if ( fLooseErrors == 1 )            return  fResult;
    Int_t       fNBins  =   aTarget->GetNbinsX()*aTarget->GetNbinsY()*aTarget->GetNbinsZ();
    for ( Int_t iBin = 1; iBin <= fNBins; iBin++ )  {
        fResult ->  SetBinError(iBin,fLooseErrors*(aTarget->GetBinError(iBin)));
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
bool                fCheckMask                  ( UChar_t fMask, Int_t iMaskBit )               {
    return  ( ((int)pow(2,iMaskBit) & fMask)       ==  (int)pow(2,iMaskBit));
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors*  fSumGraphErrors             ( TGraphAsymmErrors* gBasic, TGraphAsymmErrors* gAddition )    {
    //  Checking the consistency of TGraphs
    Int_t   fNPoints =   gBasic ->  GetN();
    if  ( fNPoints  != gAddition ->  GetN() )
    {
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
//_____________________________________________________________________________
//
TGraphAsymmErrors*  fScaleWithError             ( TGraphAsymmErrors* gBasic, Double_t fScale, Double_t fScaleErrHigh = 0., Double_t fScaleErrLow = 0. )    {
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
//
//_____________________________________________________________________________
//
TGraphAsymmErrors*  fRandomizePoints            ( TGraphAsymmErrors* gStatic, TGraphAsymmErrors* gMoveable )    {
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
        ( fRandomGen  ->  Uniform (0.,1.) ) > 0.5? fIsFluctLow = true : fIsFluctLow = false;
        if ( fIsFluctLow )  {
            fYNewValue  -= fabs(fRandomGen  ->  Gaus(fYValue,fYErrMvblLow) - fYValue);
        }   else    {
            fYNewValue  += fabs(fRandomGen  ->  Gaus(fYValue,fYErrMvblHigh) - fYValue);
        }
        fResult->SetPoint(iFit,fXValue,fYNewValue);
    }
    return  fSumGraphErrors(fResult,gMoveable);
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors*  fEfficiencycorrection       ( TH1   *fToBeCorrected, TH1    *fAccepted,  TH1    *fTotal,    Double_t fScale = 1. )  {
    //  Checking the consistency of TH1*
    //
    //  Copying accordingly the TH1*
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(fAccepted);
    //
    fResult->Divide(fAccepted,fTotal,"cl=0.683 b(1,1) mode");
    //
    Int_t   fNPoints =   fResult ->  GetN();
    for ( Int_t iFit = 0; iFit < fNPoints; iFit++ ) {
        auto    fYValueResult   =   ( fToBeCorrected ->  GetBinContent(iFit+1) );
        auto    fYErrorRstUnif  =   ( fToBeCorrected ->  GetBinError(iFit+1) );
        auto    fYValueEffic    =   ( fResult ->  GetPointY(iFit) );
        auto    fYErrorEffHigh  =   ( fResult ->  GetErrorYhigh(iFit) );
        auto    fYErrorEffLow   =   ( fResult ->  GetErrorYlow(iFit) );
        fResult ->  SetPointY       (iFit,fScale*fYValueResult/fYValueEffic);
        fResult ->  SetPointEYlow   (iFit,(fScale*fYValueResult/fYValueEffic)*sqrt(fYErrorEffHigh*fYErrorEffHigh/(fYValueEffic*fYValueEffic) + fYErrorRstUnif*fYErrorRstUnif/(fYValueResult*fYValueResult) ));
        fResult ->  SetPointEYhigh  (iFit,(fScale*fYValueResult/fYValueEffic)*sqrt(fYErrorEffLow*fYErrorEffLow/(fYValueEffic*fYValueEffic) + fYErrorRstUnif*fYErrorRstUnif/(fYValueResult*fYValueResult) ));
    }
    //
    return  fResult;
}
//
//_____________________________________________________________________________// To Be Implemented ...
//
std::vector<TGraphAsymmErrors*> fEfficiencycorrection       ( TH2   *fToBeCorrected, TH2    *fAccepted,  TH2    *fTotal,    Double_t fScale = 1. )  {
    return std::vector<TGraphAsymmErrors*>();
}
//
//_____________________________________________________________________________
//
std::vector<TGraphAsymmErrors*> fEfficiencycorrection       ( TH2   *fToBeCorrected, TH1    *fAccepted,  TH1    *fTotal,    Double_t fScale = 1. )  {
    std::vector<TGraphAsymmErrors*> fResult;
    if ( !fToBeCorrected )  { cout << "No fToBeCorrected" << endl; return fResult; }
    if ( !fAccepted )  { cout << "No fAccepted" << endl; return fResult; }
    if ( !fTotal )  { cout << "No fTotal" << endl; return fResult; }
    TGraphAsymmErrors  *fEfficiency =   new TGraphAsymmErrors(fAccepted);
    fEfficiency ->  Divide(fAccepted,fTotal,"cl=0.683 b(1,1) mode");
    for ( Int_t iHisto = 1; iHisto <= fToBeCorrected->GetNbinsX(); iHisto++ )    {
        auto    fConditional    =   fToBeCorrected->ProjectionX(Form("dd_%i",iHisto),iHisto,iHisto);
        if ( fConditional->GetEntries() == 0 ) continue;
        auto    fAddition       =   fEfficiencycorrection(fConditional,fAccepted,fTotal,fScale);
        auto    fConditionalEff =   fEfficiency->GetPointY(iHisto-1);
        auto    fConditEffHigh  =   fEfficiency->GetErrorYhigh(iHisto-1);
        auto    fConditEffLow   =   fEfficiency->GetErrorYlow(iHisto-1);
        fResult.push_back(fScaleWithError(fAddition, 1./fConditionalEff,fConditEffHigh/(fConditionalEff*fConditionalEff),fConditEffLow/(fConditionalEff*fConditionalEff)));
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
template < class Tclass >
void                fSetRainbowStyle            ( std::vector<Tclass*> fInputObjectLits, Int_t fStartIteratorAt = 0 )  {
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
bool                fIsResultAcceptable         ( Double_t fResult, Double_t fError )  {
    return (fError/fResult >= kMaximumError/100. ? false : true);
}
//
//_____________________________________________________________________________
//
//TCanvas*            fRatioPlot                  ( TH1   *fToBeCorrected, TH1    *fAccepted,  TH1    *fTotal,    Double_t fScale = 1. )  {
//    return nullptr;
//}
//
//------------------------------//
//    GENERAL FUNCTIONS         //
//------------------------------//
//
Double_t                fLevyFunc1D     ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fPT     = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fEnne   = fParams[1];
    Double_t    fSlop   = fParams[2];
    Double_t    fdNdY   = fParams[3];
    
    Double_t    fNum1   = (fEnne-1)*(fEnne-2);
    Double_t    fDen1   = fEnne*fSlop*(fEnne*fSlop+fMass*(fEnne-2));
    Double_t    fFac1   = fNum1/fDen1;
    
    Double_t    fMasT   = sqrt(fMass*fMass+fPT*fPT);
    Double_t    fNum2   = fMasT - fMass;
    Double_t    fDen2   = fEnne*fSlop;
    Double_t    fFac2   = TMath::Power((1 + fNum2/fDen2),(-fEnne));
    
    return      fPT*fdNdY*fFac1*fFac2;
}
//
//_____________________________________________________________________________
//
TF1 *                   fLevyFit1D      = new TF1 ("fLevyFunc1D",fLevyFunc1D,0.,100.,4);
//
//_____________________________________________________________________________
//
bool                    fCheckFlag      (Int_t fFlag, Int_t fToBeChecked )  {
    int fFlagCheck = (int)pow(2,fFlag);
    return (fToBeChecked & fFlagCheck) == fFlagCheck;
}
//
//_____________________________________________________________________________
//
bool                    fCheckMask      ( Int_t fToBeChecked, Int_t iMaskCheck )   {
    return  ((fToBeChecked & (int)(pow(2,iMaskCheck))) == (int)(pow(2,iMaskCheck)));
}
//
#endif
