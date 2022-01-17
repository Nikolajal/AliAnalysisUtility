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
template<   typename THXTarget_Type,
            typename THXGenRec_Type >
void
//  Author:     Anders G. Knospe, The University of Texas at Austin
//  Created:    26/01/2014
uReweightEfficiency
( THXTarget_Type* hMeasured, THXGenRec_Type* hGenerated, THXGenRec_Type* hReconstructed, TF1* hFitFunction ) {
    //  --- Initial consistency checks
    if ( uIsTHPairConsistent( hGenerated, hReconstructed ) )    { cout << "hGenerated is not compatible w/ hReconstructed" << endl; return; }
    if ( uIsTHPairRebinnable( hGenerated, hMeasured ) )         { cout << "hGenerated is not compatible w/ hMeasured" << endl; return; }
    if ( hGenerated->GetNbinsX() < hMeasured->GetNbinsX() || hGenerated->GetNbinsY() < hMeasured->GetNbinsY() || hGenerated->GetNbinsZ() < hMeasured->GetNbinsZ() ) { cout << "hGenerated must have finer binning w.r.t. hMeasured" << endl; return; }
    //
    //  --- Starting iterations
    // TODO: Parametrise max iterations
    std::vector<THXGenRec_Type*>    hGenerated_ReWork;
    std::vector<THXGenRec_Type*>    hReconstrd_ReWork;
    std::vector<THXGenRec_Type*>    hGenerated_ReBin_ReWork;
    std::vector<THXGenRec_Type*>    hReconstrd_ReBin_ReWork;
    std::vector<THXTarget_Type*>    hMeasured__ReWork;
    /*
    for ( Int_t iTer = 0; iTer < 10; iTer++ )   {
        hMeasured__ReWork.push_back();
        hGenerated_ReWork.push_back();
        hReconstrd_ReWork.push_back();
        
    }
     */
}
//
// TODO: ADD DIMENSION CHECK
template<   typename THXTarget_Type,
            typename THXGenRec_Type >
THXTarget_Type*
uEfficiencyCorrection1D
( THXTarget_Type* hTarget, THXGenRec_Type* hRec, THXGenRec_Type* hGen, Double_t fScale = 1.){
    auto    hEfficiency =   ( THXGenRec_Type* )( hRec->Clone() );
    auto    fResult     =   ( THXTarget_Type* )( hTarget->Clone() );
    if ( !hTarget ) { cout << "No hTrg" << endl; return fResult; }
    if ( !hRec )    { cout << "No hRec" << endl; return fResult; }
    if ( !hGen )    { cout << "No hGen" << endl; return fResult; }
    if ( !uIsTHPairConsistent(hRec,hGen) )      { cout << "hRec and hGen inconsistent" << endl; return fResult; }
    if ( !uIsTHPairConsistent(hRec,hTarget) )   { cout << "hRec and hTarget inconsistent" << endl; return fResult; }
    hEfficiency         ->  Divide( hRec, hGen, 1., 1., "b" );
    fResult             ->  Divide( hTarget, hEfficiency, fScale );
    return fResult;
}
//
// TODO: Generalise to n dimension correction for 1D input in rec and gen
template<   typename THXTarget_Type,
            typename THXGenRec_Type >
std::vector<THXGenRec_Type*>
uEfficiencyCorrection2D_std
( THXTarget_Type* hTarget, THXGenRec_Type* hRec, THXGenRec_Type* hGen, Double_t fScale = 1. )    {
    std::vector<THXGenRec_Type*> fResult;
    auto    hEfficiency =   ( THXGenRec_Type* )( hRec->Clone() );
    if ( !hTarget ) { cout << "No hTrg" << endl; return fResult; }
    if ( !hRec )    { cout << "No hRec" << endl; return fResult; }
    if ( !hGen )    { cout << "No hGen" << endl; return fResult; }
    if ( !uIsTHPairConsistent(hRec,hGen) )      { cout << "hRec and hGen inconsistent" << endl; return fResult; }
    cout << "ss" << endl;
    hEfficiency         ->  Divide( hRec, hGen, 1., 1., "b" );
    for ( Int_t xBin = 1; xBin <= hTarget->GetNbinsX(); xBin++ )    {
        auto    hCurrent_Slice  =   ( THXGenRec_Type* )( hRec->Clone() );
        for ( Int_t yBin = 1; yBin <= hTarget->GetNbinsY(); yBin++ )    {
            auto    kglobalBin  =   hTarget->GetBin( xBin, yBin );
            hCurrent_Slice  ->  SetBinContent( yBin, hTarget->GetBinContent(kglobalBin) / ( hEfficiency->GetBinContent(xBin) * hEfficiency->GetBinContent(yBin) ) );
            hCurrent_Slice  ->  SetBinError( yBin, hCurrent_Slice  ->  GetBinContent( yBin ) * SquareSum({ hTarget->GetBinError(kglobalBin)/hTarget->GetBinContent(kglobalBin) , hEfficiency->GetBinError(xBin)/hEfficiency->GetBinContent(xBin) , hEfficiency->GetBinError(yBin)/hEfficiency->GetBinContent(yBin) }) );
        }
        fResult.push_back( uScale( hCurrent_Slice, fScale ) );
    }
    return fResult;
}
//
// TODO: Generalise to n dimension correction for 1D input in rec and gen
template<   typename THXTarget_Type,
            typename THXGenRec_Type >
THXTarget_Type*
uEfficiencyCorrection2D
( THXTarget_Type* hTarget, THXGenRec_Type* hRec, THXGenRec_Type* hGen, Double_t fScale = 1. )    {
    if ( !hTarget ) { cout << "No hTrg" << endl; return nullptr; }
    if ( !hRec )    { cout << "No hRec" << endl; return nullptr; }
    if ( !hGen )    { cout << "No hGen" << endl; return nullptr; }
    auto    fResult     =   ( THXTarget_Type* )( hTarget->Clone() );
    auto    hEfficiency =   ( THXGenRec_Type* )( hRec->Clone() );
    if ( !uIsTHPairConsistent(hRec,hGen) )      { cout << "hRec and hGen inconsistent" << endl; return fResult; }
    hEfficiency         ->  Divide( hRec, hGen, 1., 1., "b" );
    for ( Int_t xBin = 1; xBin <= hTarget->GetNbinsX(); xBin++ )    {
        for ( Int_t yBin = 1; yBin <= hTarget->GetNbinsY(); yBin++ )    {
            auto    kGlobalBin  =   hTarget->GetBin( xBin, yBin );
            auto    kBinContent =   fScale * hTarget->GetBinContent(kGlobalBin) / ( hEfficiency->GetBinContent(xBin) * hEfficiency->GetBinContent(yBin) );
            auto    kBinError   =   SquareSum({ hTarget->GetBinError(kGlobalBin)/hTarget->GetBinContent(kGlobalBin) , hEfficiency->GetBinError(xBin)/hEfficiency->GetBinContent(xBin) , hEfficiency->GetBinError(yBin)/hEfficiency->GetBinContent(yBin) });
            fResult  ->  SetBinContent  ( kGlobalBin, kBinContent );
            fResult  ->  SetBinError    ( kGlobalBin, kBinContent * kBinError );
        }
    }
    return fResult;
}
//
template<   typename THXTarget_Type >
TCanvas*
uPlotEfficiencies
 ( std::vector<THXTarget_Type*> hTarget, std::vector<TString> fLegend = {} )  {
    TCanvas*        cDrawEfficiencies   =   new TCanvas("cDrawEfficiencies","cDrawEfficiencies",1200,1500);
    //
    TLegend*        lEfficiencies   =   new TLegend(0.625,0.88,0.88,0.7);
    lEfficiencies   ->  SetNColumns(2);
    lEfficiencies   ->  SetFillColorAlpha(0.,0.);
    lEfficiencies   ->  SetLineColorAlpha(0.,0.);
    //
    auto iTer = 0;
    TPad*   kUpperPlot  =   new TPad("kUpperPlot", "kUpperPlot", 0, 0.3, 1, 1.0);
    kUpperPlot      ->  SetLogx();
    kUpperPlot      ->  SetGridy();
    gStyle          ->  SetOptStat(0);
    kUpperPlot->SetBottomMargin(0);
    kUpperPlot->Draw();
    kUpperPlot->cd();
    for ( auto kSinglePeriodEff : hTarget )  {
        uSetHisto( kSinglePeriodEff, "EFF 1D" );
        if ( iTer != 0 )    kSinglePeriodEff ->  SetMarkerStyle ( uGetMarker(4) );
        kSinglePeriodEff    ->  SetMarkerColor ( uGetColor(iTer) );
        kSinglePeriodEff    ->  SetLineColor ( uGetColor(iTer) );
        kSinglePeriodEff    ->  Draw( "SAME" );
        if ( iTer+1 > fLegend.size() )          lEfficiencies->AddEntry( kSinglePeriodEff, kSinglePeriodEff->GetName(),  "EP" );
        else if ( !fLegend.at(iTer).IsNull() )  lEfficiencies->AddEntry( kSinglePeriodEff, fLegend.at(iTer),             "EP" );
        else                                    lEfficiencies->AddEntry( kSinglePeriodEff, kSinglePeriodEff->GetName(),  "EP" );
        iTer++;
    }
    lEfficiencies->Draw("SAME");
    //
    cDrawEfficiencies-> cd();
    TPad*   kLowerPlot  =   new TPad("kLowerPlot", "kLowerPlot", 0, 0.0, 1, 0.3);
    kLowerPlot      ->  SetLogx();
    kLowerPlot      ->  SetGridy();
    gStyle          ->  SetOptStat(0);
    gPad            ->  SetLogx();
    gPad            ->  SetGridy();
    kLowerPlot->SetTopMargin(0);
    kLowerPlot->Draw();
    kLowerPlot->cd();
    auto    kInclusiveReference =   ( THXTarget_Type* )( hTarget.at(0)->Clone() );
    iTer = 0;
    for ( auto kSinglePeriodEff : hTarget )  {
        auto    kPlotUtility        =   ( THXTarget_Type* )( kSinglePeriodEff->Clone() );
        kPlotUtility    ->  Divide( kInclusiveReference );
        kPlotUtility    ->  SetMaximum( 1.25 );
        kPlotUtility    ->  SetMinimum( 0.75 );
        kPlotUtility    ->  GetXaxis()  ->  SetTitleOffset(1.3);
        kPlotUtility    ->  GetXaxis()  ->  SetTitleSize(0.045);
        kPlotUtility    ->  GetYaxis()  ->  SetTitle("Ratio to Inclusive");
        kPlotUtility    ->  GetYaxis()  ->  SetTitleOffset(1.3);
        kPlotUtility    ->  GetYaxis()  ->  SetTitleSize(0.045);
        if ( iTer != 0 ) kPlotUtility->Draw("SAME");
        iTer++;
    }
    //
    kUpperPlot->cd();
    //
    return cDrawEfficiencies;
}
//
//
/*
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
 */
//
#endif /* AAU_Efficiency_h */
