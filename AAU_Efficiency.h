//
//  Part of the AliAnalysisUtility package
//
//  Utilities for Efficiency calculation and handling
//
//  Author              Nicola Rubini
//  Created             24/11/2021
//  Last modified       24/11/2021
#ifndef AAU_Efficiency_h
#define AAU_Efficiency_h
//
//  Global File w/ Constants and other functions
#include "AliAnalysisUtility.h"
//
template<   typename THXTarget_Type >
void
//  Author:     Anders G. Knospe, The University of Texas at Austin
//  Created:    26/01/2014
uReweightEfficiency
( THXTarget_Type* hTarget, THXTarget_Type* hGenerated, THXTarget_Type* hReconstructed, TF1* hFitFunction ) {
    if ( uIsTHPairConsistent( hGenerated, hReconstructed ) ) { cout << "hGenerated is not compatible w/ hReconstructed" << endl; return; }
    
}
            
//
#endif /* AAU_Efficiency_h */
