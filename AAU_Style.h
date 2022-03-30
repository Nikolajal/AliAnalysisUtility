//
//  Part of the AliAnalysisUtility package
//
//  Style Functions
//
//  Author              Nicola Rubini
//  Created             22/11/2021
//  Last modified       23/03/2022
//
#ifndef ALIANALYSISUTILITY_STYLE_H
#define ALIANALYSISUTILITY_STYLE_H
//
//  Global File w/ Constants and other functions
#include "AliAnalysisUtility.h"
//
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//  --- --- --- --- --- --- //! GENERAL UTILITY VARIABLES
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
const Int_t kRainbowColor[]                 =   {kRed+1,kOrange-1,kYellow+1,kSpring-1,kGreen+1,kTeal-1,kCyan+1,kAzure-1,kBlue+1,kViolet-1,kMagenta+1,kPink-1}; // Up to 12 spectra in Rainbow colours
const Int_t kRainbowMarker[]                =   {20,21,22,23,24,25}; // Up to 12 Marker in Rainbow colours
// Preferred kColors and kMarkers
const Int_t kFillColors[]                   =   {kGray+1,       kRed-10,    kBlue-9,        kGreen-8,   kMagenta-9,     kOrange-9,  kCyan-8,    kYellow-7}; // for syst bands
const Int_t kColors[]                       =   {kBlack,        kRed+1 ,    kBlue+1,        kGreen+3,   kMagenta+1,     kOrange-1,  kCyan+2,    kYellow+2};
const Int_t kMarkers[]                      =   {kOpenCircle,   kOpenSquare,kOpenDiamond,   kOpenCross, kOpenStar,      kFullCircle,kFullSquare,kFullCross, kFullDiamond,   kFullStar};
//
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//  --- --- --- --- --- --- //! GENERAL UTILITY FUNCTIONS
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
void
SetStyle
( TString kStyleOption = "ALICE" ) {
    if ( kStyleOption.Contains("ALICE") )   {
        gStyle->Reset("Plain");
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetPalette(1);
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
        gStyle->SetLegendBorderSize(0);
        gStyle->SetLegendFillColor(kWhite);
        gStyle->SetLegendFont(42);
    }
}
//
Int_t
uGetRainbowColor
 ( Int_t iTerator, Bool_t iSparse = false )    {
    if  ( iSparse ) iTerator    *=  2;
    return kRainbowColor[iTerator%12];
}                                               //  Relates to global variable kRainbowColor
//
Int_t
uGetColor
 ( Int_t iTerator, Bool_t iSparse = false )    {
    if  ( iSparse ) iTerator    *=  2;
    return kColors[iTerator%8];
}                                               //  Relates to global variable kRainbowColor
//
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
//
Int_t
uGetFillColor
 ( Int_t iTerator, Bool_t iSparse = false )    {
    if  ( iSparse ) iTerator    *=  2;
    return kFillColors[iTerator%8];
}
//
//
template<   typename THXTarget_Type1,
            typename THXTarget_Type2 >
void
uMimicStyle
 ( THXTarget_Type1*& fTarget, THXTarget_Type2*& fTemplate )   {
    //
    //  --- Marker
    fTarget->SetMarkerColor( fTemplate->GetMarkerColor() );
    fTarget->SetMarkerStyle( fTemplate->GetMarkerStyle() );
    fTarget->SetMarkerSize( fTemplate->GetMarkerSize() );
    //
    //  --- Fill
    fTarget->SetFillColor( fTemplate->GetFillColor() );
    fTarget->SetFillStyle( fTemplate->GetFillStyle() );
    //
    //  --- Line
    fTarget->SetLineColor( fTemplate->GetLineColor() );
    fTarget->SetLineStyle( fTemplate->GetLineStyle() );
    fTarget->SetLineWidth( fTemplate->GetLineWidth() );
}
//
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//  --- --- --- --- --- --- //! SPECIFIC UTILITY FUNCTIONS
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
//
//  --- --- --- --- --- --- HISTOGRAMS FUNCTIONS
//  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
//
template<   typename THXTarget_Type >
void
uSetHisto
( THXTarget_Type* hTarget, TString fOption ) {
    //
    //  Null Title for plotting
    hTarget->SetTitle("");
    //
    //  Set Final quantitites plotting
    if ( fOption.Contains("FNL") ){
        hTarget->GetXaxis()->SetBinLabel(1,"#LT Y_{#phi} #GT");
        hTarget->GetXaxis()->SetBinLabel(2,"#LT Y_{#phi#phi} #GT");
        hTarget->GetXaxis()->SetBinLabel(3,"#frac{ #LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT  } ");
        hTarget->GetXaxis()->SetBinLabel(4,"#frac{ #LT Y_{#phi#phi} #GT }{ #LT Y_{#phi} #GT^{2} }");
        hTarget->GetXaxis()->SetBinLabel(5,"#sigma^{2}_{#phi}");
        hTarget->GetXaxis()->SetBinLabel(6,"#gamma_{#phi}");
    }
    //
    // Preferred kColors and kMarkers
    //const Int_t kFillColors[] =   {kGray+1,       kRed-10,    kBlue-9,        kGreen-8,   kMagenta-9,     kOrange-9,  kCyan-8,    kYellow-7};
    //const Int_t kColors[]     =   {kBlack,        kRed+1 ,    kBlue+1,        kGreen+3,   kMagenta+1,     kOrange-1,  kCyan+2,    kYellow+2};
    //const Int_t kMarkers[]    =   {kOpenCircle,   kOpenSquare,kOpenDiamond,   kOpenCross, kOpenStar,      kFullCircle,kFullSquare,kFullCross, kFullDiamond,   kFullStar};
    //
    //  Set Stat / Syst style
    if ( fOption.Contains("STAT") ) {
        hTarget->SetMarkerStyle     ( uGetMarker(5) );
        hTarget->SetMarkerColor     ( uGetColor(1)  );
        hTarget->SetMarkerSize      ( 2 );
        hTarget->SetLineColor       ( uGetColor(1) );
        hTarget->SetFillColorAlpha  ( uGetColor(1), 0. );
    }
    if ( fOption.Contains("SYST") ) {
        hTarget->SetMarkerStyle     ( uGetMarker(5) );
        hTarget->SetMarkerColor     ( uGetColor(1)  );
        hTarget->SetMarkerSize      ( 2 );
        hTarget->SetLineColor       ( uGetColor(1) );
        hTarget->SetFillColorAlpha  ( uGetColor(1), 0.25 );
    }
    if ( fOption.Contains("SYST") || fOption.Contains("STAT") ) {
        if ( fOption.Contains("VAR3") ) hTarget->SetMarkerStyle     ( uGetMarker(1) );
        if ( fOption.Contains("VAR2") ) hTarget->SetMarkerStyle     ( uGetMarker(6) );
        if ( fOption.Contains("VAR1") ) hTarget->SetMarkerStyle     ( uGetMarker(0) );
    }
    //
    //  Set Axes for Spectra
    if ( fOption.Contains("SPT") ) {
        if ( fOption.Contains("1D") )   {
            //  --- X Axis
            hTarget->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c}) ");
            //  --- Y Axis
            hTarget->GetYaxis()->SetTitle("1/#it{N}_{ev}d#it{N}^{2}/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");
        }
        if ( fOption.Contains("2D") )  {
            //  --- X Axis
            hTarget->GetXaxis()->SetTitle("#it{p}_{T,#phi_{1}} (GeV/#it{c}) ");
            //  --- Y Axis
            hTarget->GetYaxis()->SetTitle("#it{p}_{T,#phi_{2}} (GeV/#it{c}) ");
            //  --- Z Axis
            hTarget->GetZaxis()->SetTitle("1/#it{N}_{ev}d#it{N}^{3}/(d#it{y}d#it{p}_{T,#phi_{1}}d#it{p}_{T,#phi_{2}}) (GeV/#it{c})^{-2}");
        }
        if ( fOption.Contains("12D") )  {
            //  --- X Axis
            hTarget->GetXaxis()->SetTitle("#it{p}_{T,#phi_{1}} (GeV/#it{c}) ");
            //  --- Y Axis
            hTarget->GetYaxis()->SetTitle("1/#it{N}_{ev}d#it{N}^{3}/(d#it{y}d#it{p}_{T,#phi_{1}}d#it{p}_{T,#phi_{2}}) (GeV/#it{c})^{-2}");
        }
    }
    //
    //  Set Axes for Mean PT Spectra
    if ( fOption.Contains("MPT") ) {
        //  --- X Axis
        hTarget->GetXaxis()->SetTitle("#it{p}_{T,#phi_{1}} (GeV/#it{c}) ");
        //  --- Y Axis
        hTarget->GetYaxis()->SetTitle("#LT #it{p}_{T,#phi_{2}} #GT (GeV/#it{c})");
        hTarget->GetYaxis()->SetTitleOffset(1.35);
        hTarget->GetYaxis()->SetTitleSize(0.05);
    }
    //
    //  Set Axes for Efficiency and Signal Loss
    if ( fOption.Contains("EFF") ) {
        if ( fOption.Contains("1D") )   {
            //  --- X Axis
            hTarget->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        }
        if ( fOption.Contains("12D") )  {
            //  --- X Axis
            hTarget->GetXaxis()->SetTitle("#it{p}_{T,#phi_{1}} (GeV/#it{c})");
        }
        //  --- Y Axis
        hTarget ->  GetYaxis    ()      ->SetTitle("Efficiency #times Acceptance (%)");
        hTarget ->  Scale       (100.);
        hTarget ->  SetMaximum  (100.);
        hTarget ->  SetMinimum  (0.);
    }
    if ( fOption.Contains("SYSSTACK") ) {
        hTarget ->  SetMinimum(0);
        hTarget ->  SetMaximum( hTarget -> GetMaximum()*1.3 );
        hTarget ->  GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hTarget ->  GetYaxis()->SetTitleOffset(1.5);
        hTarget ->  GetYaxis()->SetTitle("Systematic Uncertainty (%)");
    }
    //
    //  !TODO: LEGACY TO BE CHECKED AGAIN
    if ( fOption.Contains("1D") || fOption.Contains("12D") )   {
        if ( fOption.Contains("EFF") ) {
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
                //uOffset(hTarget,-100);
                hTarget->SetMaximum(10.);
                hTarget->SetMinimum(-1.);
            }
        }
    }
}
//
void
uSetHisto
( THStack* hTarget, TString fOption ) {
    //
    //  Null Title for plotting
    hTarget->SetTitle("");
    //
    if ( fOption.Contains("SYSSTACK") ) {
        hTarget ->  SetMinimum(0);
        hTarget ->  SetMaximum( hTarget -> GetMaximum()*1.3 );
        hTarget ->  GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hTarget ->  GetYaxis()->SetTitleOffset(1.5);
        hTarget ->  GetYaxis()->SetTitle("Systematic Uncertainty (%)");
    }
}
//
#endif
