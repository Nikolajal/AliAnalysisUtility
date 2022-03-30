//
//  Part of the AliAnalysisUtility package
//
//  Utilities for Graphs handling
//
//  Author              Nicola Rubini
//  Created             23/03/2022
//  Last modified       23/03/2022
//
#ifndef ALIANALYSISUTILITY_GRAPHS_H
#define ALIANALYSISUTILITY_GRAPHS_H
//
//  TODO: Implement histogram check for all functions
//
//  Global File w/ Constants and other functions
#include "AliAnalysisUtility.h"
//
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//  --- --- --- --- --- --- //! GENERAL UTILITY VARIABLES
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
//  ---
//
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//  --- --- --- --- --- --- //! GENERAL UTILITY FUNCTIONS
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
//  --- --- --- --- --- --- CONVERSION
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
// !TODO: Generalise function
template<   typename THXTarget_Type = TH1F >
TGraphErrors*
uMakeMeTGE
( THXTarget_Type* hTarget ) {
    auto fResult = new TGraphErrors();
    for ( Int_t iBin = 1; iBin <= hTarget->GetNbinsX(); iBin++ ) {
        auto xValue = hTarget->GetBinCenter(iBin);
        auto yValue = hTarget->GetBinContent(iBin);
        auto xError = 0.5*hTarget->GetBinWidth(iBin);
        auto yError = hTarget->GetBinError(iBin);
        fResult -> SetPoint     ( iBin-1, xValue, yValue );
        fResult -> SetPointError( iBin-1, xError, yError );
        cout << "xValue: " << xValue << endl;
        cout << "yValue: " << yValue << endl;
        cout << "xError: " << xError << endl;
        cout << "yError: " << yError << endl;
    }
    uMimicStyle( fResult, hTarget );
    return fResult;
}
//
#endif  /* AAU_Histograms_h */
