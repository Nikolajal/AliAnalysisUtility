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
const Bool_t    kFitScarseHisto =   kTRUE;     // Skip the fit of histograms that are below the threshold set below
const Float_t   kScarseHistoDef =   0.;        // % of entries w.r.t. total bin number
const Int_t     kScarseHistoMin =   1000.;     // N of entries
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
bool    fIsWorthFitting                         ( Tclass * aTarget )    {
    if ( !aTarget ) {
        cout << "[ERROR] You are trying to fit a null histogram! Skipping this one..." << endl;
        return false;
    }
    if (  aTarget->GetEntries() == 0. ) {
        cout << "[ERROR] You are trying to fit an empty histogram! Skipping this one..." << endl;
        return false;
    }
    if (  aTarget->GetEntries() <= kScarseHistoDef*aTarget->GetNbinsX()*aTarget->GetNbinsY()*aTarget->GetNbinsZ() + kScarseHistoMin ) {
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
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors();
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
        fResult ->  SetPointEYlow   (iFit,fScale*sqrt(fYErrorEffHigh*fYErrorEffHigh + fYErrorRstUnif*fYErrorRstUnif ));
        fResult ->  SetPointEYhigh  (iFit,fScale*sqrt(fYErrorEffLow*fYErrorEffLow + fYErrorRstUnif*fYErrorRstUnif ));
    }
    //
    return  fResult;
}
//
//_____________________________________________________________________________
//
//TCanvas*            fRatioPlot                  ( TH1   *fToBeCorrected, TH1    *fAccepted,  TH1    *fTotal,    Double_t fScale = 1. )  {
    return nullptr;
}
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

#endif
