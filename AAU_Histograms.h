//
//  Part of the AliAnalysisUtility package
//
//  Utilities for Histogram handling
//
//  Author              Nicola Rubini
//  Created             22/11/2021
//  Last modified       24/11/2021
#ifndef ALIANALYSISUTILITY_HISTOGRAMS_H
#define ALIANALYSISUTILITY_HISTOGRAMS_H
//  TODO: Implement histogram check for all functions
//
//  Global File w/ Constants and other functions
#include "AliAnalysisUtility.h"
//
//  --- GLOBAL VARIABLES
//
Int_t       iBuilderTH1_TypeCounter =   0;
//
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//  --- --- --- --- --- --- //! GENERAL UTILITY FUNCTIONS
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
template <  typename THXTarget_Type >
Int_t
uGetTHDimension
 ( THXTarget_Type*   fTarget ) {
    TObject* kTObj1DTestTarget  =   dynamic_cast< TObject* >( fTarget );
    TH1* kHist1DTestTarget  =   dynamic_cast< TH1* >( fTarget );
    if ( !kHist1DTestTarget || !kTObj1DTestTarget || !fTarget )   {
        if ( !fTarget ) cout << "[ERROR] Target is a nullptr" << endl;
        if ( !kHist1DTestTarget &&  kTObj1DTestTarget ) cout << "[ERROR] Target " << fTarget->GetName() << " is not a histogram!" << endl;
        if ( !kHist1DTestTarget && !kTObj1DTestTarget ) cout << "[ERROR] Target is not a TObject!" << endl;
        return -1;
    }
    TH2* kHist2DTestTarget  =   dynamic_cast< TH2* >( fTarget );
    TH3* kHist3DTestTarget  =   dynamic_cast< TH3* >( fTarget );
    if ( !kHist2DTestTarget && !kHist3DTestTarget ) return 1;
    if (  kHist2DTestTarget && !kHist3DTestTarget ) return 2;
    if ( !kHist2DTestTarget &&  kHist3DTestTarget ) return 3;
    return -1;
}
//
template <  typename THXTarget_Type,
            typename THXSource_Type >
Int_t
uGetTHPairDimension
( THXTarget_Type*   fTarget_1,    THXSource_Type*   fTarget_2 )  {
    auto kTarget_1_Dim  =   uGetTHDimension(fTarget_1);
    auto kTarget_2_Dim  =   uGetTHDimension(fTarget_2);
    if ( kTarget_1_Dim == kTarget_2_Dim )   return  kTarget_1_Dim;
    cout << "[ERROR] Dimension of histograms are not coherent!" << endl;
    cout << "[INFO] Dimension of Target 1 is " << kTarget_1_Dim << endl;
    cout << "[INFO] Dimension of Target 2 is " << kTarget_2_Dim << endl;
    return -1;
}
//
template <  typename THXTarget_Type,
            typename THXSource_Type >
Bool_t
uIsTHPairConsistent
( THXTarget_Type*   fTarget_1,    THXSource_Type*   fTarget_2 )  {
    auto    nDimension  = uGetTHPairDimension( fTarget_1, fTarget_2 );
    if ( nDimension < 0 ) return false;
    if ( fTarget_1->GetNbinsX() != fTarget_2->GetNbinsX() ) return false;
    if ( fTarget_1->GetNbinsY() != fTarget_2->GetNbinsY() ) return false;
    if ( fTarget_1->GetNbinsZ() != fTarget_2->GetNbinsZ() ) return false;
    for ( Int_t iBin = 1; iBin <= fTarget_1->GetNbinsX(); iBin++ ) {
        for ( Int_t jBin = 1; jBin <= fTarget_1->GetNbinsY(); jBin++ ) {
            for ( Int_t kBin = 1; kBin <= fTarget_1->GetNbinsZ(); kBin++ ) {
                if ( fTarget_1->GetXaxis()->GetBinLowEdge( iBin ) != fTarget_2->GetXaxis()->GetBinLowEdge( iBin ) ) return false;
                if ( fTarget_1->GetYaxis()->GetBinLowEdge( jBin ) != fTarget_2->GetYaxis()->GetBinLowEdge( jBin ) ) return false;
                if ( fTarget_1->GetZaxis()->GetBinLowEdge( kBin ) != fTarget_2->GetZaxis()->GetBinLowEdge( kBin ) ) return false;
            }
        }
    }
    return true;
}
//
template <  Bool_t   TCheckBoundaries = true,
            typename THXTarget_Type,
            typename THXSource_Type >
Bool_t
uIsTHPairRebinnable
( THXTarget_Type*   fTarget_1,    THXSource_Type*   fTarget_2 )  {
    auto    nDimension  = uGetTHPairDimension( fTarget_1, fTarget_2 );
    if ( nDimension < 0 ) return false;
    if ( TCheckBoundaries ) {
        if ( fTarget_1->GetXaxis()->GetXmax()() != fTarget_2->GetXaxis()->GetXmax()() ) return false;
        if ( fTarget_1->GetYaxis()->GetXmax()() != fTarget_2->GetYaxis()->GetXmax()() ) return false;
        if ( fTarget_1->GetZaxis()->GetXmax()() != fTarget_2->GetZaxis()->GetXmax()() ) return false;
        if ( fTarget_1->GetXaxis()->GetXmin()() != fTarget_2->GetXaxis()->GetXmin()() ) return false;
        if ( fTarget_1->GetYaxis()->GetXmin()() != fTarget_2->GetYaxis()->GetXmin()() ) return false;
        if ( fTarget_1->GetZaxis()->GetXmin()() != fTarget_2->GetZaxis()->GetXmin()() ) return false;
    }
    for ( Int_t iBin = 1; iBin <= fTarget_1->GetNbinsX(); iBin++ ) {
        auto    kLowBinEdge_1   =   fTarget_1->GetXaxis()->GetBinLowEdge(iBin);
        auto    kLowBinEdge_2   =   fTarget_2->GetXaxis()->GetBinLowEdge(fTarget_2->GetXaxis()->FindBin(kLowBinEdge_1));
        auto    kLowBinEdge_3   =   fTarget_1->GetXaxis()->GetBinLowEdge(fTarget_1->GetXaxis()->FindBin(kLowBinEdge_2));
        if ( kLowBinEdge_3 != kLowBinEdge_2 ) return false;
    }
    for ( Int_t jBin = 1; jBin <= fTarget_1->GetNbinsY(); jBin++ ) {
        auto    kLowBinEdge_1   =   fTarget_1->GetYaxis()->GetBinLowEdge(jBin);
        auto    kLowBinEdge_2   =   fTarget_2->GetYaxis()->GetBinLowEdge(fTarget_2->GetYaxis()->FindBin(kLowBinEdge_1));
        auto    kLowBinEdge_3   =   fTarget_1->GetYaxis()->GetBinLowEdge(fTarget_1->GetYaxis()->FindBin(kLowBinEdge_2));
        if ( kLowBinEdge_3 != kLowBinEdge_2 ) return false;
    }
    for ( Int_t kBin = 1; kBin <= fTarget_1->GetNbinsZ(); kBin++ ) {
        auto    kLowBinEdge_1   =   fTarget_1->GetZaxis()->GetBinLowEdge(kBin);
        auto    kLowBinEdge_2   =   fTarget_2->GetZaxis()->GetBinLowEdge(fTarget_2->GetZaxis()->FindBin(kLowBinEdge_1));
        auto    kLowBinEdge_3   =   fTarget_1->GetZaxis()->GetBinLowEdge(fTarget_1->GetZaxis()->FindBin(kLowBinEdge_2));
        if ( kLowBinEdge_3 != kLowBinEdge_2 ) return false;
    }
    return true;
}
//
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//  --- --- --- --- --- --- //! HISTOGRAM BUILDING FUNCTIONS
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
//  --  --  --  TODO: Make Build TH2 and TH3 and THn
template<   typename TH1_Type = TH1F,
            typename stdVec_Type,
            typename = typename std::enable_if<std::is_arithmetic<stdVec_Type>::value, stdVec_Type>::type >
TH1_Type*
uBuildTH1
 ( std::vector<stdVec_Type> fInputData, Int_t fNofBins = -1, Float_t fOffset = 0., Float_t fLowBound = 0, Float_t fHigBound = 0 )    {
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
    TH1_Type   *fBuiltTH1_Type  =   new TH1_Type(Form("TH1_Type_from_vector_%i",iBuilderTH1_TypeCounter),Form("TH1_Type_from_vector_%i",iBuilderTH1_TypeCounter),fNofBins,fLowBound,fHigBound);
    for     ( auto iValue : fInputData )    fBuiltTH1_Type->Fill( iValue + fOffset );
    iBuilderTH1_TypeCounter++;
    return  fBuiltTH1_Type;
}
//
//  --  --  --  TODO: Implement
template<   typename TH1_Type = TH2F,
            typename stdVec_Type,
            typename = typename std::enable_if<std::is_arithmetic<stdVec_Type>::value, stdVec_Type>::type >
TH1_Type*
uBuildTH2
 ( std::vector<stdVec_Type> fInputData, Int_t fNofBins = -1, Float_t fOffset = 0., Float_t fLowBound = 0, Float_t fHigBound = 0 )    {
    return  new TH1_Type;
}
//
//  --  --  --  TODO: Implement
template<   typename TH1_Type = TH3F,
            typename stdVec_Type,
            typename = typename std::enable_if<std::is_arithmetic<stdVec_Type>::value, stdVec_Type>::type >
TH1_Type*
uBuildTH3
 ( std::vector<stdVec_Type> fInputData, Int_t fNofBins = -1, Float_t fOffset = 0., Float_t fLowBound = 0, Float_t fHigBound = 0 )    {
    return  new TH1_Type;
}
//
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//  --- --- --- --- --- --- //! HISTOGRAM MANIPULATION FUNCTIONS
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
//
//  --- --- --- --- --- --- BINNING FUNCTIONS
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
//  --  --  --  TODO: Generalise
template <  typename TH1Target_Type,
            typename TH1Source_Type >
void
uRebin1D
( TH1Target_Type*   fTarget,    TH1Source_Type*   fSource )  {
    auto    kPairDimension  =   uGetPairDimension( fTarget, fSource );
    if      ( kPairDimension != 1 ) {
        cout << "[ERROR] Histograms are inconsistent, skipping rebinning" << endl;
        return;
    }
    auto fCurrentBin        =   1;
    auto fCurrentMergeN     =   0;
    auto fCurrentContent    =   0.;
    auto fCurrentError      =   0.;
    for ( Int_t iBin = 0; iBin < fSource->GetNbinsX(); iBin++ ) {
        auto    fYvalue     =   fSource ->  GetBinContent   ( iBin );
        auto    fYerror     =   fSource ->  GetBinError     ( iBin );
        auto    fXvalue     =   fSource ->  GetBinCenter    ( iBin );
        auto    fLowEdge    =   fTarget ->  GetBinLowEdge   ( fCurrentBin + 1 );
        auto    fMinEdge    =   fTarget ->  GetBinLowEdge   ( fCurrentBin );
        if ( fXvalue < fMinEdge )   continue;
        if ( fXvalue > fLowEdge )   {
            fCurrentBin++;
            iBin--;
            fCurrentMergeN  =   0;
            fCurrentContent =   0;
            fCurrentError   =   0;
            continue;
        }
        fCurrentMergeN++;
        fCurrentContent +=  fYvalue;
        fCurrentError   +=  fYerror*fYerror;
        fTarget ->  SetBinContent   ( fCurrentBin, fCurrentContent / ( fCurrentMergeN*1. ) );
        fTarget ->  SetBinError     ( fCurrentBin, sqrt(fCurrentError) / ( fCurrentMergeN*1. ) );
    }
}
//
//  --  --  --  TODO: Implement
template <  typename TH2Taregt_Type,
            typename TH2Source_Type >
void
uRebin2D
( TH2Taregt_Type*   fTarget,    TH2Source_Type*   fSource )  {
    return;
}
//
//  --  --  --  TODO: Implement
template <  typename TH3Target_Type,
            typename TH3Source_Type >
void
uRebin3D
( TH3Target_Type*   fTarget,    TH3Source_Type*   fSource )  {
    return;
}
//
//  --  --  --  TODO: Avoid in generalisation
template <  typename THXTarget_Type,
            typename THXSource_Type >
void
uRebin
( THXTarget_Type*   fTarget,    THXSource_Type*   fSource )  {
    auto    kPairDimension  =   uGetPairDimension( fTarget, fSource );
    if      ( kPairDimension == 1 ) uRebin1D(fTarget,fSource);
    else if ( kPairDimension == 2 ) uRebin2D(fTarget,fSource);
    else if ( kPairDimension == 3 ) uRebin3D(fTarget,fSource);
    else    {
        cout << "[ERROR] Histograms are inconsistent, skipping rebinning" << endl;
    }
}
//
template<   typename stdArr_Type,
            typename = typename std::enable_if<std::is_arithmetic<stdArr_Type>::value, stdArr_Type>::type >
Float_t*
uGetUniformBinningArray
( stdArr_Type fMinBin, stdArr_Type fMaxBin, Int_t fNBins ) {
    stdArr_Type*    fResult =   new stdArr_Type[fNBins+1];
    for ( Int_t iBin = 0; iBin < fNBins+1; iBin++ )    fResult[iBin]    =   fMinBin + ( iBin )*( fMaxBin - fMinBin ) / ( static_cast<Float_t>( fNBins ) );
    return  fResult;
}
//
template<   typename TInput1,
            typename TInput2,
            typename TInput3,
            typename = typename std::enable_if<std::is_arithmetic<TInput1>::value, TInput1>::type,
            typename = typename std::enable_if<std::is_arithmetic<TInput2>::value, TInput2>::type,
            typename = typename std::enable_if<std::is_arithmetic<TInput3>::value, TInput3>::type >
std::pair<Int_t,Float_t*>
uUniformBinning
( TInput1 fBinWidth, TInput2 fLowEdge, TInput3 fHighEdge )  {
    Int_t       nBinNumber  =   ( fHighEdge - fLowEdge ) / fBinWidth;
    if ( ( ( ( fHighEdge - fLowEdge ) / fBinWidth ) - nBinNumber ) >= fBinWidth ) nBinNumber++;
    Float_t*    fResult     =   new Float_t [ nBinNumber+1 ];
    for ( Int_t iBin = 0; iBin <= nBinNumber; iBin++ )  fResult[iBin]   =   fLowEdge + iBin * fBinWidth;
    if ( fResult[ nBinNumber ] != fHighEdge )  cout << "[WARNING] [uUniformBinning] Changed High Edge from " << fHighEdge << " to " << fResult[ nBinNumber ] << endl;
    return  std::pair<Int_t,Float_t*> ( nBinNumber, fResult );
}
//
//  --- --- --- --- --- --- DATA HANDLING FUNCTIONS
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
template<   typename TH1_Type >
void
uOffset
( TH1_Type* hTarget, Double_t kOffset, Bool_t kAbsolute = false ){
    for ( Int_t iBin = 1; iBin <= hTarget->GetNbinsX(); iBin++ ) {
        auto fYValue = hTarget->GetBinContent(iBin);
        if ( kAbsolute )    hTarget->SetBinContent(iBin, fabs( fYValue+kOffset ) );
        else                hTarget->SetBinContent(iBin, fYValue+kOffset );
    }
}
//
template<   typename TH1_Type >
void
uAbsolute
( TH1_Type* hTarget ){
    uOffset( hTarget, 0, kTRUE );
}
//
template<   typename THXTarget_Type,
            typename TInput1 = Double_t,
            typename TInput2 = Double_t,
            typename = typename std::enable_if<std::is_arithmetic<TInput1>::value, TInput1>::type,
            typename = typename std::enable_if<std::is_arithmetic<TInput2>::value, TInput2>::type >
THXTarget_Type*
uScale
 ( THXTarget_Type* hTarget, TInput1 fScaleFactor, TInput2 fScaleError = -1. )  {
    auto    nDimension  =   uGetTHDimension( hTarget );
    auto    fResult     =   (THXTarget_Type*)(hTarget->Clone());
    if ( nDimension < 0 ) return fResult;
    for ( Int_t iBin = 1; iBin <= fResult->GetNbinsX(); iBin++ ) {
        for ( Int_t jBin = 1; jBin <= fResult->GetNbinsY(); jBin++ ) {
            for ( Int_t kBin = 1; kBin <= fResult->GetNbinsZ(); kBin++ ) {
                auto    kGlobalBin  =   fResult->GetBin( iBin, jBin, kBin );
                fResult ->  SetBinContent   ( kGlobalBin, fScaleFactor  * hTarget   ->  GetBinContent   ( iBin ) );
                if ( fScaleError <= -1 )fResult ->  SetBinError     ( kGlobalBin,   fScaleFactor    *   hTarget   ->  GetBinError     ( iBin ) );
                else                    fResult ->  SetBinError     ( kGlobalBin,   ( 1 + fScaleError ) * fScaleError     *   hTarget   ->  GetBinError     ( iBin ) );
            }
        }
    }
    return fResult;
}
//
// TODO: Generalise w/ THXTarget_Type_3 as return
template <  Bool_t      TSquareSum          = kTRUE,
            typename    THXTarget_Type_1    = TH1F,
            typename    THXTarget_Type_2    = TH1F  >
THXTarget_Type_1*
uSumErrors
 ( THXTarget_Type_1* hTarget_1, THXTarget_Type_2* hTarget_2 )    {
    THXTarget_Type_1*   fResult =   (THXTarget_Type_1*)(hTarget_1->Clone());
    if ( !uIsTHPairConsistent( hTarget_1, hTarget_2 ) ) return fResult;
    for ( Int_t iBin = 1; iBin <= fResult->GetNbinsX(); iBin++ ) {
        for ( Int_t jBin = 1; jBin <= fResult->GetNbinsY(); jBin++ ) {
            for ( Int_t kBin = 1; kBin <= fResult->GetNbinsZ(); kBin++ ) {
                auto    kGlobalBin  =   fResult->GetBin( iBin, jBin, kBin );
                auto    kNewError   =   0.;
                if ( TSquareSum )   kNewError   =   SquareSum( { hTarget_1->GetBinError(iBin), hTarget_2->GetBinError(iBin) } );
                else                kNewError   =   hTarget_1->GetBinError(iBin) + hTarget_2->GetBinError(iBin);
                fResult ->  SetBinError( iBin, kNewError );
            }
        }
    }
    return  fResult;
}
//
template<   Bool_t TSquareSum = true,
            typename THXTarget_Type_1,
            typename THXTarget_Type_2 = THXTarget_Type_1 >
THXTarget_Type_2*
uRandomisePoints
( THXTarget_Type_1* hTarget ) {
    auto    nDimension  =   uGetTHDimension( hTarget );
    THXTarget_Type_2*  fResult =   (THXTarget_Type_2*)(hTarget->Clone());
    if ( nDimension < 0 )   return  fResult;
    for ( Int_t iBin = 1; iBin <= fResult->GetNbinsX(); iBin++ ) {
        for ( Int_t jBin = 1; jBin <= fResult->GetNbinsY(); jBin++ ) {
            for ( Int_t kBin = 1; kBin <= fResult->GetNbinsZ(); kBin++ ) {
                auto    kGlobalBin  =   fResult->GetBin         ( iBin, jBin, kBin );
                auto    kBinContent =   fResult->GetBinContent  ( kGlobalBin );
                auto    kBinError   =   fResult->GetBinError    ( kGlobalBin );
                auto    kNewBinCont =   uRandomGen -> Gaus( kBinContent, kBinError );
                fResult ->  SetBinContent( kGlobalBin, kNewBinCont );
            }
        }
    }
    return fResult;
}
//
template<   Bool_t TSquareSum = true,
            typename THXTarget_Type_1,
            typename THXTarget_Type_2,
            typename THXTarget_Type_3 = THXTarget_Type_1 >
THXTarget_Type_3*
uRandomisePoints
( THXTarget_Type_1* hTarget_1, THXTarget_Type_2* hTarget_2  ) {
    auto    nDimension  =   uGetTHDimension( hTarget_1 );
    THXTarget_Type_3*  fResult =   (THXTarget_Type_3*)(hTarget_1->Clone());
    if ( nDimension < 0 )   return  fResult;
    fResult =   uRandomisePoints( hTarget_1 );
    fResult =   uSumErrors      ( hTarget_1, hTarget_2 );
    return fResult;
}
//
//  --- --- --- --- --- --- WRITE/READ FILE FUNCTIONS
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
// TODO: Generalise for N cycles, for labels, for any TH1* if input std::vector& you can desume the TH1 type
template<   Int_t TDimension = 0,
            typename THXTarget_Type = TH1F,
            Bool_t   TResverse  =   false,
            typename TReturn_Type = typename TNVector< TDimension, THXTarget_Type* >::type >
TReturn_Type
uLoadHistograms
( TFile* kDataFile, TString kHistogramName, TString kNewName = ""  ) {
    TReturn_Type fResult;
    uInternalLoadHisto( fResult, kDataFile, kHistogramName, kNewName );
    return fResult;
}
template<   typename THXTarget_Type = TH1F >
void
uInternalLoadHisto
( THXTarget_Type*& fResult, TFile* kDataFile, TString kHistogramName, TString kNewName = "" ) {
    if (  kDataFile->Get(kHistogramName) ) fResult = ( THXTarget_Type* )( ( kDataFile->Get( kHistogramName ) )->Clone() );
    if ( !kNewName.IsNull() )   fResult->SetName((kNewName));
    if ( !fResult ) cout << "[ERROR] No Histogram match found for " << kHistogramName.Data() << endl;
    return fResult;
}
template<   typename THXTarget_Type = TH1F >
void
uInternalLoadHisto
( std::vector<THXTarget_Type*>& fResult, TFile* kDataFile, TString kHistogramName, TString kNewName = "" ) {
    auto iTer = 0;
    while ( true ) {
        if ( !kDataFile->Get(Form(kHistogramName,iTer)) ) break;
        fResult.push_back( (THXTarget_Type*)( ( kDataFile->Get( Form( kHistogramName, iTer ) ) )->Clone() ) );
        if ( !kNewName.IsNull() )    fResult.at(iTer)->SetName(Form(kNewName,iTer));
        iTer++;
    }
    if ( fResult.size() == 0 ) cout << "[ERROR] No Histogram match found for " << kHistogramName.Data() << endl;
}
template<   typename THXTarget_Type = TH1F >
void
uInternalLoadHisto
( std::vector<std::vector<THXTarget_Type*>>& fResult, TFile* kDataFile, TString kHistogramName, TString kNewName = "" ) {
    auto iTer = 0;
    while ( true ) {
        std::vector<THXTarget_Type*>   kUtility;
        auto jTer = 0;
        while ( true ) {
            if ( !kDataFile->Get(Form(kHistogramName,iTer,jTer)) ) break;
            kUtility.push_back( (THXTarget_Type*)( ( kDataFile->Get( Form( kHistogramName, iTer, jTer ) ) )->Clone() ) );
            if ( !kNewName.IsNull() )    kUtility.at(jTer)->SetName(Form(kNewName,iTer,jTer));
            jTer++;
        }
        if ( jTer == 0 ) break;
        fResult.push_back( kUtility );
        iTer++;
    }
    if ( fResult.size() == 0 ) cout << "[ERROR] No Histogram match found for " << kHistogramName.Data() << endl;
}
template<   typename THXTarget_Type = TH1F >
void
uInternalLoadHisto
( std::vector<std::vector<std::vector<THXTarget_Type*>>>& fResult, TFile* kDataFile, TString kHistogramName, TString kNewName = "" ) {
    auto iTer = 0;
    while ( true ) {
        std::vector<std::vector<THXTarget_Type*>>  kUtilit2;
        auto jTer = 0;
        while ( true ) {
            std::vector<THXTarget_Type*>   kUtility;
            auto kTer = 0;
            while ( true ) {
                if ( !kDataFile->Get(Form(kHistogramName,iTer,jTer,kTer)) ) break;
                kUtility.push_back( (THXTarget_Type*)( ( kDataFile->Get( Form( kHistogramName, iTer, jTer, kTer ) ) )->Clone() ) );
                if ( !kNewName.IsNull() )    kUtility.at(kTer)->SetName(Form(kNewName,iTer,jTer,kTer));
                kTer++;
            }
            if ( kTer == 0 ) break;
            kUtilit2.push_back( kUtility );
            jTer++;
        }
        if ( jTer == 0 ) break;
        fResult.push_back( kUtilit2 );
        iTer++;
    }
    if ( fResult.size() == 0 ) cout << "[ERROR] No Histogram match found for " << kHistogramName.Data() << endl;
}
//
template<   typename THXTarget_Type,
            typename TInput1 = Float_t,
            typename = typename std::enable_if<std::is_arithmetic<TInput1>::value > >
void
uAddSumHistogram
 ( std::vector<THXTarget_Type*> &hTarget, TString kNewName = "", std::vector<TInput1> kWeights = {} ) {
    // TODO: null vec, warning weitghs less than histos
    //if ( size == 0 ) warnign erorr
    auto    iTer = 0;
    auto    hResult =   (THXTarget_Type*)(hTarget.at(0)->Clone());
    hResult ->  Reset();
    for ( auto kHisto : hTarget )    {
        if ( kWeights.size() < iTer+1 ) hResult -> Add( kHisto );
        else                            hResult -> Add( kHisto, kWeights.at(iTer) );
        iTer++;
    }
    if ( kNewName.IsNull() )    hResult ->  SetName( "SumHisto_from_uAddSumHistogram" );
    else                        hResult ->  SetName( kNewName );
    hTarget                     .insert( hTarget.begin(), hResult );
}
//
template<   typename THXTarget_Type,
            typename TInput1,
            typename = typename std::enable_if<std::is_arithmetic<TInput1>::value > >
void
uAddSumHistogram
 ( std::vector<std::vector<THXTarget_Type*>> &hTarget, TString kNewName = "", std::vector<TInput1> kWeights = {} ) {
    std::vector<THXTarget_Type*>  fResult;
    for ( Int_t jTer = 0; jTer < hTarget.at(0).size(); jTer++ )   {
        std::vector<THXTarget_Type*> kUtility;
        for ( Int_t iTer = 0; iTer < hTarget.size(); iTer++ )   {
            kUtility.push_back( hTarget.at(iTer).at(jTer) );
        }
        uAddSumHistogram( kUtility, Form( kNewName, jTer), kWeights );
        fResult.push_back( kUtility.at(0) );
    }
    hTarget.push_back( fResult );
}
//
//TODO: Implement
template<   typename THXTarget_Type = TH1F >
std::vector<std::vector<THXTarget_Type*>>
uReverseStructure
( std::vector<std::vector<THXTarget_Type*>> hTarget ) {
    std::vector<std::vector<THXTarget_Type*>> fResult;
    for ( Int_t iTer = 0; iTer < hTarget.at(0).size(); iTer++ )  {
        std::vector<THXTarget_Type*> kUtility;
        for ( Int_t jTer = 0; jTer < hTarget.size(); jTer++ ) {
            kUtility.push_back( (THXTarget_Type*)hTarget[jTer][iTer]->Clone() );
        }
        fResult.push_back(  *(new std::vector<THXTarget_Type*> (kUtility)) );
    }
    return fResult;
}
//
template<   typename THXTarget_Type = TH1F >
std::vector<std::vector<std::vector<THXTarget_Type*>>>
uReverseStructure
( std::vector<std::vector<std::vector<THXTarget_Type*>>> hTarget ) {
    std::vector<std::vector<std::vector<THXTarget_Type*>>> fResult;
    for ( Int_t iTer = 0; iTer < hTarget.at(0).size(); iTer++ )  {
        std::vector<std::vector<THXTarget_Type*>> kUtilit2;
        for ( Int_t jTer = 0; jTer < hTarget.size(); jTer++ ) {
            std::vector<THXTarget_Type*> kUtility;
            for ( Int_t kTer = 0; kTer < hTarget.size(); kTer++ ) {
                kUtility.push_back( (THXTarget_Type*)hTarget[jTer][kTer][iTer]->Clone() );
            }
            kUtilit2.push_back(  *(new std::vector<THXTarget_Type*> (kUtility)) );
        }
        fResult.push_back(  *(new std::vector<std::vector<THXTarget_Type*>> (kUtilit2)) );
    }
    return fResult;
}
//
template<   Int_t TDimension = 2,
            typename THXTarget_Type = TH1F,
            typename TTarget_Type = typename TNVector< TDimension, THXTarget_Type* >::type >
TTarget_Type
uReverseStructure
( TTarget_Type hVecTarget ) {
    return uReverseStructure(hVecTarget);
}
//
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//  --- --- --- --- --- --- //! TODO: LEGACY, TO BE CHECKED AGAIN
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
//  --  --  --  RETRO COMPATIBILITY, TO BE CLEANED
TH1F*
uBuildTH1F
 ( std::vector<Float_t> fInputData, Int_t fNofBins = -1, Float_t fOffset = 0., Float_t fLowBound = 0, Float_t fHigBound = 0 )    {
    return uBuildTH1<TH1F>(fInputData,fNofBins,fOffset,fLowBound,fHigBound);
}
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
template < class Tclass >
Tclass*                 fSumErrors              ( Tclass* gBasic, Tclass* gAddition )    {
    Tclass  *fResult =   new Tclass(*gBasic);
    for ( Int_t iBin = 0; iBin < gBasic->GetNbinsX(); iBin++ ) {
        fResult ->  SetBinError( iBin, SquareSum( { gBasic->GetBinError(iBin), gAddition->GetBinError(iBin) } ) );
    }
    //
    return  fResult;
}
TH1F*                   fEfficiencycorrection   ( TH1   *fToBeCorrected, TH1    *fAccepted,  TH1   *fTotal,    Double_t fScale = 1. )  {
    TH1F   *fEfficiency =   (TH1F*)fAccepted->Clone();
    TH1F   *fResult     =   (TH1F*)fToBeCorrected->Clone();
    fEfficiency         ->  Divide(fAccepted,fTotal,1.,1.,"b");
    fResult             ->  Divide(fToBeCorrected,fEfficiency,fScale);
    return  fResult;
}
std::vector<TH1F*>      fEfficiencycorrection   ( TH2   *fToBeCorrected, TH1    *fAccepted,  TH1    *fTotal,    Double_t fScale = 1. )  {
    std::vector<TH1F*> fResult;
    if ( !fToBeCorrected )  { cout << "No fToBeCorrected" << endl; return fResult; }
    if ( !fAccepted )  { cout << "No fAccepted" << endl; return fResult; }
    if ( !fTotal )  { cout << "No fTotal" << endl; return fResult; }
    TH1F   *fEfficiency =   (TH1F*)fAccepted->Clone();
    fEfficiency         ->  Divide(fAccepted,fTotal,1.,1.,"b");
    for ( Int_t iHisto = 1; iHisto <= fToBeCorrected->GetNbinsY(); iHisto++ )    {
        auto    fConditional    =   fToBeCorrected->ProjectionY(Form("dd_%i",iHisto),iHisto,iHisto);
        TH1F*    fCorrCondit    =   fEfficiencycorrection( fConditional, fAccepted, fTotal, fScale );
        fResult.push_back( uScale( fCorrCondit, fEfficiency->GetBinContent(iHisto), fEfficiency->GetBinError(iHisto) ) );
    }
    return fResult;
}
std::vector<TH1F*>      uRandomizePointsSymm    ( std::vector<TH1F*>  gStatic, std::vector<TH1F*>  gMoveable )  {
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
std::vector<TH1D*>      uRandomizePointsSymm    ( std::vector<TH1D*>  gStatic, std::vector<TH1D*>  gMoveable )  {
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
            hTarget->SetMarkerStyle(kMarkers[5]);
            // Colour scheme
            hTarget->SetMarkerColor(kColors[2]);
            hTarget->SetLineColor(kColors[2]);
            hTarget->SetFillColorAlpha(kFillColors[2],0.33);
            //
            hTarget->SetOption("EP L");
            if ( fOption.Contains("EFF2") ) {
                hTarget->SetMarkerStyle(kMarkers[9]);
                hTarget->SetMarkerColor(kColors[1]);
                hTarget->SetLineColor(kColors[1]);
                hTarget->SetFillColorAlpha(kFillColors[1],0.33);
                hTarget->SetOption("EP L SAME");
            }
            if ( fOption.Contains("SL") ) {
                // Y-Axis title
                hTarget->GetYaxis()->SetTitle("Signal Loss (%)");
                uOffset(hTarget,-100);
                hTarget->SetMaximum(10.);
                hTarget->SetMinimum(-1.);
            }
        }
        if ( fOption.Contains("SPT") ) {
            
             // Preferred kColors and kMarkers
             //const Int_t kFillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
             //const Int_t kColors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
             //const Int_t kMarkers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};
            
            hTarget->SetMarkerStyle(kMarkers[3]);
            hTarget->SetMarkerColor(kColors[2]);
            hTarget->SetMarkerSize(1);
            hTarget->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            hTarget->GetYaxis()->SetTitleOffset(1.5);
            hTarget->GetYaxis()->SetTitle("1/N_{ev}dN^{2}/(dydp_{T})");
            if ( fOption.Contains("12D") )  {
                hTarget->GetYaxis()->SetTitle("1/N_{ev}dN^{3}/(dydp_{T,#phi_{1}}dp_{T,#phi_{2}})");
                hTarget->GetXaxis()->SetTitle("#it{p}_{T,#phi_{1}} (GeV/#it{c})");
            }
        }
        if ( fOption.Contains("MPT") ) {
            hTarget->SetMarkerStyle(kMarkers[3]);
            hTarget->SetMarkerColor(kColors[2]);
            hTarget->SetMarkerSize(1);
            hTarget->GetXaxis()->SetTitleOffset(1.1);
            hTarget->GetXaxis()->SetTitle("#it{p}_{T,#phi_{1}} (GeV/#it{c})");
            hTarget->GetYaxis()->SetTitleOffset(1.5);
            hTarget->GetYaxis()->SetTitle("#LT #it{p}_{T,#phi_{2}} #GT (GeV/#it{c})");
        }
        if ( fOption.Contains("STAT") ) {
            hTarget->SetLineColor(kColors[2]);
            hTarget->SetFillColorAlpha(kColors[2],0.33);
            hTarget->SetOption("PE2");
        }
        if ( fOption.Contains("SYST") ) {
            hTarget->SetLineColor(kColors[3]);
            hTarget->SetFillColorAlpha(kColors[3],0.);
            hTarget->SetOption("PE2");
        }
    } else if ( fOption.Contains("2D") )   {
        
    } else if ( fOption.Contains("3D") )   {
        cout << " Buu 3D " << endl;
    } else  {
        cout << " CANT GUESS DIMENSION " << endl;
    }
}

template<   typename THXTarget_Type = TH1F >
Double_t
uGetFWHM
 ( TH1* hTarget ) {
    Double_t    kMaxVal =   0;
    Double_t    kMaxX   =   0;
    Double_t    kMinX   =   0;
    for ( Int_t iBin = 1; iBin <= hTarget->GetNbinsX(); iBin++ )    {
        kMaxVal =   max( kMaxVal, hTarget->GetBinContent(iBin) );
    }
    for ( Int_t iBin = 1; iBin <= hTarget->GetNbinsX(); iBin++ )    {
        cout << "Searching: " << iBin << " Val: " << hTarget->GetBinContent(iBin) << "  " << kMaxVal/2. << endl;
        if ( hTarget->GetBinContent(iBin) > kMaxVal/2. ) {
            kMinX = hTarget->GetBinCenter(iBin);
            break;
        }
    }
    cout << "dddd" << endl;
    for ( Int_t iBin = hTarget->GetNbinsX(); iBin > 0; iBin-- )    {
        cout << "Searching: " << iBin << " Val: " << hTarget->GetBinContent(iBin) << "  " << kMaxVal/2. << endl;
        if ( hTarget->GetBinContent(iBin) > kMaxVal/2. ) {
            kMaxX = hTarget->GetBinCenter(iBin-1);
            break;
        }
    }
    cout << "ddrrdd" << endl;
    return kMaxX - kMinX;
    /*
    int bin1 = h1->FindFirstBinAbove(h1->GetMaximum()/2);
    int bin2 = h1->FindLastBinAbove(h1->GetMaximum()/2);
    double fwhm = h1->GetBinCenter(bin2) - h1->GetBinCenter(bin1);
     */
}
//  --  --  --  RETRO COMPATIBILITY, TO BE CLEANED
//
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//  --- --- --- --- --- --- //! TODO: NOT PROPER, TO BE CHECKED AGAIN
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
std::vector<TH1F*>
uMakeRatio
 ( std::vector<TH1F*> hTargets )   {
    std::vector<TH1F*> fResult;
    for ( auto kHisto : hTargets )  {
        auto    kNewHisto   =   (TH1F*)(kHisto->Clone());
        kNewHisto->Divide(kHisto,hTargets.at(0));
        fResult.push_back( kNewHisto );
    }
    return fResult;
}
//
#endif  /* AAU_Histograms_h */
