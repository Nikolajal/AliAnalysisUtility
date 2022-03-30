//
//  Part of the AliAnalysisUtility package
//
//  Global Functions
//
//  Author              Nicola Rubini
//  Created             22/11/2021
//  Last modified       22/11/2021
#ifndef ALIANALYSISUTILITY_GLOBALUTILFUNCTIONS_H
#define ALIANALYSISUTILITY_GLOBALUTILFUNCTIONS_H
//
//  --- Global Functions
//
//  --- --- Type definition helpers
//
template< Int_t TDimension, typename TContent_Type >    struct TNVector                         { typedef std::vector< typename TNVector< TDimension - 1, TContent_Type >::type > type; };
template< typename TContent_Type >                      struct TNVector < 0, TContent_Type >    { typedef TContent_Type type; };
//
//---------------------------------------------------------------------------------------------------------------------------------------------------   Square Sum of Inputs
//  !TODO: Name Change  -> +u...
Double_t
SquareSum
 ( std::initializer_list<Double_t> list )  {
    Double_t    fResult =   0;
    for ( auto element : list ) fResult +=  element*element;
    return  TMath::Sqrt(fResult);
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
template< typename stdVec_Type, typename = typename std::enable_if<std::is_arithmetic<stdVec_Type>::value, stdVec_Type>::type >
void
uCleanOutsiders
 ( std::vector<stdVec_Type> &fInputData, Float_t nSigmaCut = 10. )    {
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
            if ( fabs(*iValue - MeanOfDistribution) >= nSigmaCut*STDVOfDistribution )   {
                fCheckAgain = true;
                iValue = fInputData.erase(iValue);
            }   else    {
                ++iValue;
            }
        }
    }
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
//  TODO: pass stringfromat from outside
TString
uCurrentDateTimeOffset
( Int_t fOffset ) {
    time_t     now = time(0) + fOffset;
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%H:%M:%S", &tstruct);
    return TString(buf);
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
//  TODO: pass stringfromat from outside
TString
uCurrentDateTime
() {
    return uCurrentDateTimeOffset(0);
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
//  !TODO: Name Change f -> u
bool
fCheckMask
 ( Int_t fToBeChecked, Int_t iMaskCheck, Bool_t fOnlyThis = false )    {
    if ( ( fToBeChecked == 0 ) && ( iMaskCheck == -1 ) )    return true;
    if ( ( fToBeChecked <= 0 ) || ( iMaskCheck <= 0 ) )     return false;
    if ( !fOnlyThis )   return  ( fToBeChecked & BIT(iMaskCheck) );
    else                return  ( fToBeChecked ==   BIT(iMaskCheck) );
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
Double_t
uBarlowPar
 ( Double_t  fStandard, Double_t fStdError, Double_t fVariatin, Double_t fVarError )  {
    auto    fSigmaStd       =   fStdError*fStdError;
    auto    fSigmaVar       =   fVarError*fVarError;
    auto    fSigmaDff       =   TMath::Sqrt( fabs( fSigmaStd - fSigmaVar ) );
    if      ( fSigmaDff    == 0 )  return false;
    auto    fParameter      =   ( fStandard - fVariatin )/fSigmaDff;
    return  fParameter;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
//  !TODO: Name Change f -> u
bool
fBarlowCheck
 ( Double_t  fStandard, Double_t fStdError, Double_t fVariatin, Double_t fVarError )  {
    //return false;
    return  ( fabs ( uBarlowPar    ( fStandard, fStdError, fVariatin, fVarError ) )    <= 1 );
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
//                                                                                                                                                  --- BENCHMARK UTILITIES
//---------------------------------------------------------------------------------------------------------------------------------------------------
void
fStartTimer
 ( TString fTimerName )    {
    uBenchmark->Reset();
    uBenchmark->Start(fTimerName.Data());
    printf("[INFO] Starting %s \n", fTimerName.Data());
    fflush(stdout);
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
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
//---------------------------------------------------------------------------------------------------------------------------------------------------
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
    cout << Form(kMSG_PrintTimer.Data(), uCurrentDateTime().Data(),  fShowPrintEvnt,  fSuffix.Data(), 100.*fProcessedFrac, fShowSpeedEvnt,  fSuffix.Data(), fShowElapsedHor, fShowElapsedMin,  fShowElapsedSec, fShowEstimatHor, fShowEstimatMin, fShowEstimatSec, uCurrentDateTimeOffset( fRealEstimatSec ).Data()) << flush;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
//
template < typename TTarget_Type >
void
push_to_front
( std::vector<TTarget_Type>& vTarget, TTarget_Type tAddition )  {
    vTarget.insert( vTarget.begin(), tAddition );
}
//
//---------------------------------------------------------------------------------------------------LEGACY-----------------------------------------------
Double_t
fGammaPhiValue
 ( Double_t fYieldPhi, Double_t fYieldPhiPhi )  {
    return  2*fYieldPhiPhi/fYieldPhi -fYieldPhi;
}
Double_t
fGammaPhiError
 ( Double_t fYieldPhi, Double_t fYieldPhiPhi, Double_t fErrorPhi, Double_t fErrorPhiPhi)  {
    auto    fPar1   =   2*fErrorPhiPhi/fYieldPhi;
    auto    fPar2   =   (2*fYieldPhiPhi/(fYieldPhi*fYieldPhi)+1)*fErrorPhi;
    return  fPar1 + fPar2;
}
Double_t
fSigmaPhiValue
 ( Double_t fYieldPhi, Double_t fYieldPhiPhi )  {
    return  2*fYieldPhiPhi + fYieldPhi - fYieldPhi*fYieldPhi;
}
Double_t
fSigmaPhiError
 ( Double_t fYieldPhi, Double_t fYieldPhiPhi, Double_t fErrorPhi, Double_t fErrorPhiPhi)  {
    return SquareSum( { 2*fErrorPhiPhi, (-1+2*fYieldPhi)*fErrorPhi } );
}
Double_t
uFindPercentile
 ( TH1D* hReference, Double_t kCurrent )  {
    auto kBin = hReference->GetXaxis()->FindBin(kCurrent);
    auto kTot = hReference->Integral(-1,-1);
    auto kPrt = hReference->Integral(1,kBin);
    return kPrt/kTot;
}
#endif
