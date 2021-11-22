//
//  Part of the AliAnalysisUtility package
//
//  Style Functions
//
//  Author              Nicola Rubini
//  Created             22/11/2021
//  Last modified       22/11/2021
#ifndef ALIANALYSISUTILITY_STYLE_H
#define ALIANALYSISUTILITY_STYLE_H
//
//  Global File w/ Constants and other functions
#include "AliAnalysisUtility.h"
//
//>>    GLOBAL VARIABLES
//
const Int_t         kRainbowColor[]                 =   {kRed+1,kOrange-1,kYellow+1,kSpring-1,kGreen+1,kTeal-1,kCyan+1,kAzure-1,kBlue+1,kViolet-1,kMagenta+1,kPink-1}; // Up to 12 spectra in Rainbow colours
const Int_t         kRainbowMarker[]                =   {20,21,22,23,24,25}; // Up to 12 Marker in Rainbow colours
// Preferred kColors and kMarkers
const Int_t         kFillColors[]                   =   {kGray+1,       kRed-10,    kBlue-9,        kGreen-8,   kMagenta-9,     kOrange-9,  kCyan-8,    kYellow-7}; // for syst bands
const Int_t         kColors[]                       =   {kBlack,        kRed+1 ,    kBlue+1,        kGreen+3,   kMagenta+1,     kOrange-1,  kCyan+2,    kYellow+2};
const Int_t         kMarkers[]                      =   {kOpenCircle,   kOpenSquare,kOpenDiamond,   kOpenCross, kOpenStar,      kFullCircle,kFullSquare,kFullCross, kFullDiamond,   kFullStar};
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
//>>    STYLE FUNCTIONS
//
Int_t
uGetRainbowColor
 ( Int_t iTerator, Bool_t iSparse = false )    {
    if  ( iSparse ) iTerator    *=  2;
    return kRainbowColor[iTerator%12];
}                                               //  Relates to global variable kRainbowColor
Int_t
uGetColor
 ( Int_t iTerator, Bool_t iSparse = false )    {
    if  ( iSparse ) iTerator    *=  2;
    return kColors[iTerator%8];
}                                               //  Relates to global variable kRainbowColor
Int_t
uGetMarker
 ( Int_t iTerator, Int_t OpenOnly = -1 )    {
    if  ( OpenOnly == 0 )   {
        return kMarkers[iTerator%5];
    } else if  ( OpenOnly == 1 )   {
        iTerator    +=  5;
        return kMarkers[iTerator%5];
    } else   {
        return kMarkers[iTerator%10];
    }
}                                               //  Relates to global variable kRainbowColor
Int_t
uGetFillColor
 ( Int_t iTerator, Bool_t iSparse = false )    {
    if  ( iSparse ) iTerator    *=  2;
    return kFillColors[iTerator%8];
}

//------------------------LEGACY
Int_t
fGetRainbowColor
 ( Int_t iTerator, Bool_t iSparse = false )    {
    if  ( iSparse ) iTerator    *=  2;
    return kRainbowColor[iTerator%12];
}                                               //  Relates to global variable kRainbowColor
Int_t
fGetColor
 ( Int_t iTerator, Bool_t iSparse = false )    {
    if  ( iSparse ) iTerator    *=  2;
    return kColors[iTerator%8];
}                                               //  Relates to global variable kRainbowColor
Int_t
fGetMarker
( Int_t iTerator, Int_t OpenOnly = -1 )    {
    if  ( OpenOnly == 0 )   {
        return kMarkers[iTerator%5];
    } else if  ( OpenOnly == 1 )   {
        iTerator    +=  5;
        return kMarkers[iTerator%5];
    } else   {
        return kMarkers[iTerator%10];
    }
}                                               //  Relates to global variable kRainbowColor
Int_t
fGetFillColor
 ( Int_t iTerator, Bool_t iSparse = false )    {
    if  ( iSparse ) iTerator    *=  2;
    return kFillColors[iTerator%8];
}

#endif
