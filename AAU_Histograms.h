//
//  Re-Include protection
#ifndef ALIANALYSISUTILITY_HISTOGRAMS_H
#define ALIANALYSISUTILITY_HISTOGRAMS_H
//
//  Global File w/ Constants and other functions
#include "AliAnalysisUtility.h"
//
//>>    GLOBAL VARIABLES
//
Int_t       iBuilderTH1_TypeCounter =   0;
//
//>>
//>>    GENERAL UTILITY FUNCTIONS
//>>
//
template < class THXTarget_Type, class THXSource_Type >
Int_t
uGetPairDimension
( THXTarget_Type*   fTarget,    THXSource_Type*   fSource )  {
    //  Check is a histogram
    TH1* kHist1DTestTarget  =   dynamic_cast< TH1* >( fTarget );
    TH1* kHist1DTestSource  =   dynamic_cast< TH1* >( fSource );
    if ( !kHist1DTestTarget || !kHist1DTestSource )   {
        if ( !kHist1DTestTarget )   cout << "[ERROR] Target is not a histogram!" << endl;
        if ( !kHist1DTestSource )   cout << "[ERROR] Source is not a histogram!" << endl;
        return -1;
    }
    TH2* kHist2DTestTarget  =   dynamic_cast< TH2* >( fTarget );
    TH2* kHist2DTestSource  =   dynamic_cast< TH2* >( fSource );
    TH3* kHist3DTestTarget  =   dynamic_cast< TH3* >( fTarget );
    TH3* kHist3DTestSource  =   dynamic_cast< TH3* >( fSource );
    //
    auto kTargetIs1D        =   false;
    auto kSourceIs1D        =   false;
    auto kTargetIs2D        =   false;
    auto kSourceIs2D        =   false;
    auto kTargetIs3D        =   false;
    auto kSourceIs3D        =   false;
    //
    if ( !kHist2DTestTarget && !kHist3DTestTarget ) kTargetIs1D = true;
    if ( !kHist2DTestSource && !kHist3DTestSource ) kSourceIs1D = true;
    if (  kHist2DTestTarget && !kHist3DTestTarget ) kTargetIs2D = true;
    if (  kHist2DTestSource && !kHist3DTestSource ) kSourceIs2D = true;
    if ( !kHist2DTestTarget &&  kHist3DTestTarget ) kTargetIs3D = true;
    if ( !kHist2DTestSource &&  kHist3DTestSource ) kSourceIs3D = true;
    //
    if      ( kTargetIs1D && kSourceIs1D ) return 1;
    else if ( kTargetIs2D && kSourceIs2D ) return 2;
    else if ( kTargetIs3D && kSourceIs3D ) return 3;
    else    {
        // !TODO: Signal the dimension of targt source
        cout << "[ERROR] Dimension of histograms are not coherent!" << endl;
        return -1;
    }
}
//
//>>
//>>    FUNCTIONS TO BUILD HISTOGRAMS
//>>
//
template< typename TH1_Type = TH1F, typename stdVec_Type, typename = typename std::enable_if<std::is_arithmetic<stdVec_Type>::value, stdVec_Type>::type >
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
template< typename TH1_Type = TH2F, typename stdVec_Type, typename = typename std::enable_if<std::is_arithmetic<stdVec_Type>::value, stdVec_Type>::type >
TH1_Type*
uBuildTH2
 ( std::vector<stdVec_Type> fInputData, Int_t fNofBins = -1, Float_t fOffset = 0., Float_t fLowBound = 0, Float_t fHigBound = 0 )    {
    return  new TH1_Type;
}
//
//  --  --  --  TODO: Implement
template< typename TH1_Type = TH3F, typename stdVec_Type, typename = typename std::enable_if<std::is_arithmetic<stdVec_Type>::value, stdVec_Type>::type >
TH1_Type*
uBuildTH3
 ( std::vector<stdVec_Type> fInputData, Int_t fNofBins = -1, Float_t fOffset = 0., Float_t fLowBound = 0, Float_t fHigBound = 0 )    {
    return  new TH1_Type;
}
//
//>>
//>>    FUNCTIONS TO MANIPULATE HISTOGRAMS
//>>
//
//>>    >>  BINNING FUNCTIONS
//
template < class TH1Target_Type, class TH1Source_Type >
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
        // !TODO: Generalise to 2D
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
template < class TH2Taregt_Type, class TH2Source_Type >
void
uRebin2D
( TH2Taregt_Type*   fTarget,    TH2Source_Type*   fSource )  {
    return;
}
//
//  --  --  --  TODO: Implement
template < class TH3Target_Type, class TH3Source_Type >
void
uRebin3D
( TH3Target_Type*   fTarget,    TH3Source_Type*   fSource )  {
    return;
}
//
template < class THXTarget_Type, class THXSource_Type >
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
template< typename stdArr_Type, typename = typename std::enable_if<std::is_arithmetic<stdArr_Type>::value, stdArr_Type>::type >
Float_t*
uGetUniformBinningArray
( stdArr_Type fMinBin, stdArr_Type fMaxBin, Int_t fNBins ) {
    stdArr_Type*    fResult =   new stdArr_Type[fNBins+1];
    for ( Int_t iBin = 0; iBin < fNBins+1; iBin++ )    fResult[iBin]    =   fMinBin + ( iBin )*( fMaxBin - fMinBin ) / ( static_cast<Float_t>( fNBins ) );
    return  fResult;
}
//
//>>
//>>    --  --  --  TODO: Clean
//>>    LEGACY, TO BE CHECKED AGAIN
//>>
//
//  --  --  --  RETRO COMPATIBILITY, TO BE CLEANED
TH1F*
uBuildTH1F
 ( std::vector<Float_t> fInputData, Int_t fNofBins = -1, Float_t fOffset = 0., Float_t fLowBound = 0, Float_t fHigBound = 0 )    {
    return uBuildTH1<TH1F>(fInputData,fNofBins,fOffset,fLowBound,fHigBound);
}
//  --  --  --  RETRO COMPATIBILITY, TO BE CLEANED
//_____________________________________________________________________________
//                                                  //  !TODO: Make the kMarkerstyle in a kArray
//                                                  //  Relates to global variable kRainbowColor
template < class Tclass >
void                    fSetRainbowStyle            ( std::vector<Tclass**> fInputObjectLits, Int_t fStartIteratorAt = 0 )  {
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
template < class Tclass >
void                    fSetUniformBinning          ( Tclass *fArrBin, Tclass fMinBin, Tclass fMaxBin, Int_t fNPoints )  {
    for (int i = 0; i <= fNPoints; i++ )
    {
        fArrBin[i] = fMinBin+(i)*(fMaxBin - fMinBin)/(static_cast<Float_t>(fNPoints));
    }
}
////_____________________________________________________________________________
void
uOffset
( TH1* hTarget, Double_t kOffset, Bool_t kAbsolute = false ){
    for ( Int_t iBin = 1; iBin <= hTarget->GetNbinsX(); iBin++ ) {
        auto fYValue = hTarget->GetBinContent(iBin);
        if ( kAbsolute )    hTarget->SetBinContent(iBin, fabs( fYValue+kOffset ) );
        else                hTarget->SetBinContent(iBin, fYValue+kOffset );
    }
}
//
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
            
            /*
             // Preferred kColors and kMarkers
             const Int_t kFillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
             const Int_t kColors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
             const Int_t kMarkers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};
             */
            
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
//
TCanvas*
uPlotSpectrum
( TH1* hTarget, TH1* hTrSyst, TString fOption = "" ){
    //
    SetStyle();
    //
    TCanvas    *cDrawResult =   new TCanvas();
    gStyle->SetOptStat(0);
    if ( fOption.Contains("SPT") ) gPad->SetLogy();
    uSetHisto(hTarget,fOption + TString(" STAT"));
    uSetHisto(hTrSyst,fOption + TString(" SYST"));
    hTarget->SetMaximum(1.25*max(hTarget->GetMaximum(),hTrSyst->GetMaximum()));
    hTarget->SetMinimum(0.75*min(hTarget->GetMinimum(),hTrSyst->GetMinimum()));
    if ( fOption.Contains("SPT") )  hTarget->SetMaximum(2.0*max(hTarget->GetMaximum(),hTrSyst->GetMaximum()));
    if ( fOption.Contains("SPT") )  hTarget->SetMinimum(0.5*min(hTarget->GetMinimum(),hTrSyst->GetMinimum()));
    //
    TLegend    *lLegend;
    if ( fOption.Contains("R") )    lLegend =   new TLegend(0.65,0.35,0.85,0.5);
        else                        lLegend =   new TLegend(0.2,0.35,0.4,0.5);
    lLegend     ->  SetFillColorAlpha(0.,0.);
    lLegend     ->  AddEntry    (hTarget,"Data","P");
    lLegend     ->  AddEntry    (hTarget,"Stat","F");
    lLegend     ->  AddEntry    (hTrSyst,"Syst","F");
    //
    hTarget->Draw();
    hTrSyst->Draw("SAME E1");
    lLegend->Draw("SAME");
    //
    if ( fOption.Contains("R") )    {
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.65, 0.3,"ALICE");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.65, 0.25,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.65, 0.2,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
    } else  {
        uLatex->SetTextFont(60);
        uLatex->SetTextSize(0.05);
        uLatex->DrawLatexNDC(0.20, 0.3,"ALICE");
        uLatex->SetTextFont(42);
        uLatex->SetTextSize(0.04);
        uLatex->DrawLatexNDC(0.20, 0.25,"pp #sqrt{#it{s}}= 7 TeV");
        uLatex->DrawLatexNDC(0.20, 0.2,"#phi #rightarrow K^{+}K^{-}, |#it{y}|<0.5");
    }
    //
    return cDrawResult;
}
//
void
uSetHisto
( TGraphMultiErrors* hTarget, TString fOption = "" ){
    hTarget->SetTitle("");
    hTarget->GetYaxis()->SetTitle("1/N_{ev}dN/dy");
    //hTarget->GetXaxis()->SetNdivisions(2);
    //hTarget->GetXaxis()->SetBinLabel(hTarget->GetXaxis()->FindBin(1),"#LT Y_{1#phi} #GT");
    //hTarget->GetXaxis()->SetBinLabel(hTarget->GetXaxis()->FindBin(2),"#LT Y_{2#phi} #GT");
    //hTarget->GetXaxis()->LabelsOption("h");
    hTarget->GetXaxis()->SetLabelSize(0.08);
    hTarget->SetLineColorAlpha(0.,0.);
    hTarget->SetMarkerStyle(21);
    hTarget->SetMarkerColor(kRed);
    hTarget->GetAttLine(0)->SetLineColor(38);
    hTarget->GetAttLine(1)->SetLineColor(46);
    hTarget->GetAttFill(0)->SetFillColorAlpha(38,0.33);
    hTarget->GetAttFill(1)->SetFillColorAlpha(46,0.33);
}
//
//  --  --  Data Handling Utilities  --  --  //
//
//_____________________________________________________________________________
//
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
//
template < class Tclass >
Tclass                   *fSumErrors                  ( Tclass* gBasic, Tclass* gAddition )    {
    Tclass  *fResult =   new Tclass(*gBasic);
    for ( Int_t iBin = 0; iBin < gBasic->GetNbinsX(); iBin++ ) {
        fResult ->  SetBinError( iBin, SquareSum( { gBasic->GetBinError(iBin), gAddition->GetBinError(iBin) } ) );
    }
    //
    return  fResult;
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors      *fScaleWithError             ( TGraphAsymmErrors* gBasic, Double_t fScale, Double_t fScaleErrHigh = 0., Double_t fScaleErrLow = 0. )    {
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
TH1F                   *fScaleWithError             ( TH1F* gBasic, Double_t fScale, Double_t fScaleError = 0. )    {
    TH1F  *fResult =   new TH1F(*gBasic);
    for ( Int_t iBin = 1; iBin <= fResult->GetNbinsX(); iBin++ ) {
        fResult->SetBinContent( iBin, gBasic->GetBinContent(iBin)/fScale );
        fResult->SetBinError( iBin, (gBasic->GetBinContent(iBin)/fScale)*SquareSum( { gBasic->GetBinError( iBin )/gBasic->GetBinContent(iBin), fScaleError/fScale } ) );
    }
    return  fResult;
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors      *uRandomizePoints            ( TGraphAsymmErrors* gStatic, TGraphAsymmErrors* gMoveable )    {
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
        ( uRandomGen  ->  Uniform (0.,1.) ) > 0.5? fIsFluctLow = true : fIsFluctLow = false;
        if ( fIsFluctLow )  {
            fYNewValue  -= fabs(uRandomGen  ->  Gaus(fYValue,fYErrMvblLow) - fYValue);
        }   else    {
            fYNewValue  += fabs(uRandomGen  ->  Gaus(fYValue,fYErrMvblHigh) - fYValue);
        }
        fResult->SetPoint(iFit,fXValue,fYNewValue);
    }
    return  fSumErrors(fResult,gMoveable);
}
//
template < class Tclass >
Tclass*                 uRandomizePoints            ( Tclass* gStatic, Tclass* gMoveable )    {
    Tclass* fResult =   (Tclass*)(gStatic->Clone());
    for ( Int_t iBin = 0; iBin < gStatic->GetNbinsX(); iBin++ ) {
        auto    fBinContent =   gStatic->GetBinContent(iBin+1);
        auto    fMoveError  =   gMoveable->GetBinError(iBin+1);
        auto fTest = uRandomGen -> Gaus(fBinContent,fMoveError);
        fResult->SetBinContent  ( iBin+1, fTest );
        fResult->SetBinError    ( iBin+1, fMoveError );
    }
    return  fSumErrors(fResult,gMoveable);
}
std::vector<TH1F*>      uRandomizePoints            ( std::vector<TH1F*>  gStatic, std::vector<TH1F*>  gMoveable )    {
    std::vector<TH1F*>  fResult;
    auto iTer = 0;
    for ( auto hStatic : gStatic ) {
        fResult.push_back(uRandomizePoints((TH1F*)hStatic->Clone(),(TH1F*)gMoveable.at(iTer)->Clone()));
        iTer++;
    }
    return  fResult;
}
std::vector<TH1F*>      uRandomizePointsSymm        ( std::vector<TH1F*>  gStatic, std::vector<TH1F*>  gMoveable )    {
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
std::vector<TH1D*>      uRandomizePointsSymm        ( std::vector<TH1D*>  gStatic, std::vector<TH1D*>  gMoveable )    {
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
//
//
//_____________________________________________________________________________
//
TH1F                    *fEfficiencycorrection       ( TH1   *fToBeCorrected, TH1    *fAccepted,  TH1   *fTotal,    Double_t fScale = 1. )  {
    TH1F   *fEfficiency =   (TH1F*)fAccepted->Clone();
    TH1F   *fResult     =   (TH1F*)fToBeCorrected->Clone();
    fEfficiency         ->  Divide(fAccepted,fTotal,1.,1.,"b");
    fResult             ->  Divide(fToBeCorrected,fEfficiency,fScale);
    return  fResult;
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors      *fEfficiencycorrection       ( TGraphAsymmErrors *fToBeCorrected, TH1    *fAccepted,  TH1    *fTotal,    Double_t fScale = 1. )  {
    //  Copying accordingly the TH1*
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*fToBeCorrected);
    //
    fResult->Divide(fAccepted,fTotal,"cl=0.683 b(1,1) mode");
    //
    Int_t   fNPoints =   fResult ->  GetN();
    for ( Int_t iFit = 0; iFit < fNPoints; iFit++ ) {
        auto    fYValueResult   =   ( fToBeCorrected ->  GetPointY(iFit) );
        auto    fYErrorRstHigh  =   ( fToBeCorrected ->  GetErrorYhigh(iFit) );
        auto    fYErrorRstLow   =   ( fToBeCorrected ->  GetErrorYlow(iFit) );
        auto    fYValueEffic    =   ( fResult ->  GetPointY(iFit) );
        auto    fYErrorEffHigh  =   ( fResult ->  GetErrorYhigh(iFit) );
        auto    fYErrorEffLow   =   ( fResult ->  GetErrorYlow(iFit) );
        fResult ->  SetPointY       (iFit,fScale*fYValueResult/fYValueEffic);
        fResult ->  SetPointEYlow   (iFit,(fScale*fYValueResult/fYValueEffic)*sqrt(fYErrorEffHigh*fYErrorEffHigh/(fYValueEffic*fYValueEffic) + fYErrorRstHigh*fYErrorRstHigh/(fYValueResult*fYValueResult) ));
        fResult ->  SetPointEYhigh  (iFit,(fScale*fYValueResult/fYValueEffic)*sqrt(fYErrorEffLow*fYErrorEffLow/(fYValueEffic*fYValueEffic) + fYErrorRstLow*fYErrorRstLow/(fYValueResult*fYValueResult) ));
    }
    //
    return  fResult;
}
//
//_____________________________________________________________________________
//                                                          // !TODO: To Be Implemented (...)
std::vector<TGraphAsymmErrors*> fEfficiencycorrection       ( TH2   *fToBeCorrected, TH2    *fAccepted,  TH2    *fTotal,    Double_t fScale = 1. )  {
    return std::vector<TGraphAsymmErrors*>();
}
//
std::vector<TH1F*> fEfficiencycorrection       ( TH2   *fToBeCorrected, TH1    *fAccepted,  TH1    *fTotal,    Double_t fScale = 1. )  {
    std::vector<TH1F*> fResult;
    if ( !fToBeCorrected )  { cout << "No fToBeCorrected" << endl; return fResult; }
    if ( !fAccepted )  { cout << "No fAccepted" << endl; return fResult; }
    if ( !fTotal )  { cout << "No fTotal" << endl; return fResult; }
    TH1F   *fEfficiency =   (TH1F*)fAccepted->Clone();
    fEfficiency         ->  Divide(fAccepted,fTotal,1.,1.,"b");
    for ( Int_t iHisto = 1; iHisto <= fToBeCorrected->GetNbinsY(); iHisto++ )    {
        auto    fConditional    =   fToBeCorrected->ProjectionY(Form("dd_%i",iHisto),iHisto,iHisto);
        TH1F*    fCorrCondit    =   fEfficiencycorrection( fConditional, fAccepted, fTotal, fScale );
        fResult.push_back( fScaleWithError( fCorrCondit, fEfficiency->GetBinContent(iHisto), fEfficiency->GetBinError(iHisto) ) );
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
TGraphAsymmErrors*              fEfficiencycorrection       ( TGraphAsymmErrors    *fToBeCorrected,    TGraphAsymmErrors   *fEfficiency,     Double_t fScale = 1. )  {
    //  Checking the consistency of TGraphs
    Int_t   fNPoints =   fToBeCorrected ->  GetN();
    if  ( fNPoints  != fEfficiency ->  GetN() )   {
        cout << "[ERROR] Systematics and Statistics do not have the same number of points! Skipping this one..." << endl;
        return nullptr;
    }
    TGraphAsymmErrors  *fResult =   new TGraphAsymmErrors(*fToBeCorrected);
    for ( Int_t iHisto = 1; iHisto <= fNPoints; iHisto++ )    {
        auto    fTargetY        =   fToBeCorrected->GetPointY(iHisto-1);
        auto    fEfficnY        =   fEfficiency->GetPointY(iHisto-1);
        auto    fTargetYlow     =   fToBeCorrected->GetErrorYlow(iHisto-1);
        auto    fEfficnYlow     =   fEfficiency->GetErrorYlow(iHisto-1);
        auto    fTargetYhigh    =   fToBeCorrected->GetErrorYlow(iHisto-1);
        auto    fEfficnYhigh    =   fEfficiency->GetErrorYlow(iHisto-1);
        fResult->SetPointY      (iHisto-1,  fTargetY*fEfficnY*fScale);
        fResult->SetPointEYlow  (iHisto-1,  sqrt(fTargetYlow*fTargetYlow+fEfficnYlow*fEfficnYlow));
        fResult->SetPointEYhigh (iHisto-1,  sqrt(fTargetYhigh*fTargetYhigh+fEfficnYhigh*fEfficnYhigh));
    }
    return fResult;
}
//
//_____________________________________________________________________________
//
//  --  --  Specific purpose histogram generation  --  --  //
//
//_____________________________________________________________________________
//>>    LEGACY, TO BE CHECKED AGAIN
//
TGraphAsymmErrors*
fTH1_to_TGAsymmErrors
 ( TH1*  hTarget )    {
    TGraphAsymmErrors  *fResult                     =   new TGraphAsymmErrors();
    Int_t       fNBins  =   hTarget->GetNbinsX();
    for ( Int_t iBin = 1; iBin <= fNBins; iBin++ )  {
        auto    fXValue =   hTarget         ->GetBinCenter(iBin);
        auto    fYValue =   hTarget         ->GetBinContent(iBin);
        auto    fXError =   0.5*( hTarget->GetBinLowEdge(iBin+1) - hTarget->GetBinLowEdge(iBin) );
        auto    fYError =   hTarget         ->GetBinError(iBin);
        fResult         ->  SetPoint        (iBin-1,fXValue,fYValue);
        fResult         ->  SetPointError   (iBin-1,fXError,fXError,fYError,fYError);
    }
    fResult             ->SetMaximum(hTarget->GetMaximum());
    fResult             ->SetMinimum(hTarget->GetMinimum());
    return fResult;
}
//
TGraphAsymmErrors**
fTH2_to_TGAsymmErrors
 ( TH2*  hTarget )    {
    Int_t       fNBinsY =   hTarget->GetNbinsY();
    TGraphAsymmErrors **fResult                     =   new TGraphAsymmErrors*[fNBinsY];
    for ( Int_t iBin = 1; iBin <= fNBinsY; iBin++ )  {
        fResult[iBin-1] =   fTH1_to_TGAsymmErrors(hTarget ->  ProjectionX(Form("%i",iBin),iBin,iBin));
    }
    return fResult;
}
//
TCanvas                *uPlotReferenceValue         ( TGraphAsymmErrors*  hMeasured,    Float_t fReference,   Float_t fRefError, TString fLabel = "PDG Value" )    {
    TCanvas    *fResult         =   new TCanvas();
    //
    auto        fNPoints        =   hMeasured->GetN();
    auto        fXLow           =   hMeasured->GetPointX(0);
    auto        fXLowErr        =   hMeasured->GetErrorXlow(0);
    auto        fXHig           =   hMeasured->GetPointX(fNPoints-1);
    auto        fXHigErr        =   hMeasured->GetErrorXhigh(fNPoints-1);
                hMeasured       ->  SetMaximum(max((float)(hMeasured->GetMaximum()+0.25*(hMeasured->GetMaximum()-hMeasured->GetMinimum())), fReference+fRefError*5));
                hMeasured       ->  SetMinimum(min((float)(hMeasured->GetMinimum()-0.25*(hMeasured->GetMaximum()-hMeasured->GetMinimum())), fReference-fRefError*5));
    //
                hMeasured       ->  Draw("APE");
    //
    TLegend    *fLegend         =   new TLegend(0.6,0.9,0.9,0.75);
                fLegend         ->  SetLineColorAlpha(1,0.);
                fLegend         ->  SetFillColorAlpha(1,0.);
                fLegend         ->  SetNColumns(2);
    //
    TH1F      **fPlotReference  =   new TH1F   *[7];
    for ( Int_t iTer = 0; iTer < 7; iTer ++ )   {
        fPlotReference[iTer]    =   new TH1F    (Form("%i",iTer),Form("%i",iTer),1,(fXLow-fXLowErr*1.5),(fXHig+fXHigErr*1.5));
        fPlotReference[iTer]    ->  SetBinContent(1, fReference + ( iTer -3 )*fRefError );
        fPlotReference[iTer]    ->  SetLineStyle(10);
        fPlotReference[iTer]    ->  SetLineWidth(3);
        fPlotReference[iTer]    ->  SetLineColorAlpha(920+(4-fabs(iTer-3)),0.75);
        fPlotReference[iTer]    ->  Draw("same ][");
    }
    fPlotReference[3]           ->  SetLineColorAlpha(2,0.75);
    fPlotReference[3]           ->  SetLineWidth(5);
    //
    fLegend                     ->  AddEntry( fPlotReference[3],    fLabel.Data(),  "L");
    fLegend                     ->  AddEntry( fPlotReference[2],    "#pm 1 #sigma", "L");
    fLegend                     ->  AddEntry( fPlotReference[1],    "#pm 2 #sigma", "L");
    fLegend                     ->  AddEntry( fPlotReference[0],    "#pm 3 #sigma", "L");
    fLegend                     ->  Draw("same");
    //
    return fResult;
}
TCanvas                *uPlotReferenceValue         ( TH1*  hMeasured,    Float_t fReference,   Float_t fRefError, TString fLabel = "PDG Value"  )    {
    return uPlotReferenceValue( fTH1_to_TGAsymmErrors(hMeasured), fReference, fRefError, fLabel );
}
//
#endif
