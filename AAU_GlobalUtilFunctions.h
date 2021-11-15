//
//  Re-Include protection
#ifndef ALIANALYSISUTILITY_GLOBALUTILFUNCTIONS_H
#define ALIANALYSISUTILITY_GLOBALUTILFUNCTIONS_H
//
//>>    GLOBAL FUNCTIONS
//
Double_t
SquareSum
 ( std::initializer_list<Double_t> list )  {
    Double_t    fResult =   0;
    for ( auto element : list ) fResult +=  element*element;
    return  TMath::Sqrt(fResult);
}
//  --  --  --  TODO: Add Nsigma cut from outside
template< typename stdVec_Type, typename = typename std::enable_if<std::is_arithmetic<stdVec_Type>::value, stdVec_Type>::type >                                            //  Relates to global variable kRainbowColor
void
uCleanOutsiders
 ( std::vector<stdVec_Type> &fInputData )    {
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
//
#endif
