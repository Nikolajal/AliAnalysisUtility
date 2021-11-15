//
//  Re-Include protection
#ifndef ALIANALYSISUTILITY_GLOBALUTILFUNCTIONS_H
#define ALIANALYSISUTILITY_GLOBALUTILFUNCTIONS_H
//
//>>    GLOBAL FUNCTIONS
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
uCurrentDateTime
() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%H:%M:%S", &tstruct);
    return TString(buf);
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
//  !TODO: Name Change f -> u
bool
fCheckMask
 ( Int_t fToBeChecked, Int_t iMaskCheck, Bool_t fOnlyThis = false )    {
    if ( fToBeChecked == 0  )   return false;
    if ( !fOnlyThis )   return  ( fToBeChecked &    BIT(iMaskCheck) );
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
//
#endif
