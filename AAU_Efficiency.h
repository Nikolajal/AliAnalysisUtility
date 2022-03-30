//
//  Part of the AliAnalysisUtility package
//
//  Utilities for Efficiency calculation and handling
//
//  Author              Nicola Rubini
//  Mail                nicola.rubini@cern.ch
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
//  --- Author:         Anders G. Knospe, The University of Texas at Austin
//  --- Created:        26/01/2014
//  --- Last Modified:  25/04/2019
// !TODO: To be implemented
uReweightEfficiency
( THXTarget_Type* hMeasured, THXGenRec_Type* hGenerated, THXGenRec_Type* hReconstructed, TF1* hFitFunction ) {
    //  --- Initial consistency checks
    if ( uIsTHPairConsistent( hGenerated, hReconstructed ) )    { cout << "hGenerated is not compatible w/ hReconstructed" << endl; return; }
    if ( uIsTHPairRebinnable( hGenerated, hMeasured ) )         { cout << "hGenerated is not compatible w/ hMeasured" << endl;      return; }
    if ( hGenerated->GetNbinsX() < hMeasured->GetNbinsX() || hGenerated->GetNbinsY() < hMeasured->GetNbinsY() || hGenerated->GetNbinsZ() < hMeasured->GetNbinsZ() ) { cout << "hGenerated must have finer binning w.r.t. hMeasured" << endl; return; }
    //
    //  --- Starting iterations
    // TODO: Parametrise max iterations
    std::vector<THXGenRec_Type*>    hGenerated_ReWork;
    std::vector<THXGenRec_Type*>    hReconstrd_ReWork;
    std::vector<THXGenRec_Type*>    hGenerated_ReBin_ReWork;
    std::vector<THXGenRec_Type*>    hReconstrd_ReBin_ReWork;
    std::vector<THXTarget_Type*>    hMeasured__ReWork;
    //
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
 ( std::vector<THXTarget_Type*> hTarget, std::vector<TString> fLegend = {}, TString fNewName = "", Bool_t kSignalLoss = false )  {
    TCanvas*        cDrawEfficiencies   =   new TCanvas(fNewName,fNewName,1200,1500);
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
        if ( kSignalLoss )  uSetHisto( kSinglePeriodEff, "EFF SL 1D" );
        else                uSetHisto( kSinglePeriodEff, "EFF 1D" );
        if ( iTer != 0 )    kSinglePeriodEff ->  SetMarkerStyle ( uGetMarker(4) );
        kSinglePeriodEff    ->  SetMarkerColor ( uGetColor(iTer) );
        kSinglePeriodEff    ->  SetLineColor ( uGetColor(iTer) );
        if ( iTer == 0 && kSignalLoss ) kSinglePeriodEff    ->  SetMaximum(10);
        if ( iTer == 0 && kSignalLoss ) kSinglePeriodEff    ->  GetXaxis()  ->  SetTitle("Signal Loss (%)");
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
        kPlotUtility    ->  SetMaximum( 1.35 );
        kPlotUtility    ->  SetMinimum( 0.65 );
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
template<   typename THXGenRec_Type1D,
            typename THXGenRec_Type2D >
TCanvas*
uEfficiencyCompare_1D_2D
// !TODO: Solve temporary HACK
( THXGenRec_Type1D* hRec_1D, THXGenRec_Type1D* hGen_1D, THXGenRec_Type2D* hRec_2D, THXGenRec_Type2D* hGen_2D, Bool_t kSigLoss = false )    {
    //
    //  --- Result TCanvas
    TCanvas*    cDrawEfficiencyCompare  =   new TCanvas(Form("cDrawEfficiencyCompare_%i",iBuilderTH1_TypeCounter),Form("cDrawEfficiencyCompare_%i",iBuilderTH1_TypeCounter),1000,1000);
    cDrawEfficiencyCompare->Divide(2,4);
    //
    //  --- Canonical Checks on Input
    if ( !hRec_1D ) { cout << "No hRec_1D" << endl; return new TCanvas; }
    if ( !hGen_1D ) { cout << "No hGen_1D" << endl; return new TCanvas; }
    if ( !hRec_2D ) { cout << "No hRec_2D" << endl; return new TCanvas; }
    if ( !hGen_2D ) { cout << "No hGen_2D" << endl; return new TCanvas; }
    if ( !uIsTHPairConsistent(hRec_1D,hGen_1D) ) { cout << "hRec_1D and hGen_1D inconsistent" << endl; return new TCanvas; }
    if ( !uIsTHPairConsistent(hRec_2D,hGen_2D) ) { cout << "hRec_2D and hGen_2D inconsistent" << endl; return new TCanvas; }
    if ( !uIsTHPairConsistent(hRec_1D,hRec_2D->ProjectionX("tmp",1,1)) ) { cout << "hRec_1D and hRec_2D inconsistent" << endl; return new TCanvas; }
    //
    //  --- Legend
    TLegend*    lDrawLegend =   new TLegend(0.65,0.72,0.95,0.95);
    lDrawLegend->SetLineColorAlpha(kWhite,0.0);
    lDrawLegend->SetFillColorAlpha(kWhite,0.0);
    //
    //  --- Generate Efficiencies
    auto    hEff_1D     =    ( THXGenRec_Type1D* )( hRec_1D->Clone(Form("%s_copy_%i",hRec_1D->GetName(),iBuilderTH1_TypeCounter) ) );
    auto    hEff_2D     =    ( THXGenRec_Type2D* )( hRec_2D->Clone(Form("%s_copy_%i",hRec_1D->GetName(),iBuilderTH1_TypeCounter) ) );
    hEff_1D ->  Divide( hRec_1D, hGen_1D, 1., 1., "b" );
    hEff_2D ->  Divide( hRec_2D, hGen_2D, 1., 1., "b" );
    for ( Int_t xBin = 1; xBin <= hRec_1D->GetNbinsX(); xBin++ )    {
        cDrawEfficiencyCompare->cd( xBin );
        //
        gPad->SetTopMargin      (0.02);
        gPad->SetBottomMargin   (0.20);
        gPad->SetRightMargin    (0.02);
        gPad->SetLeftMargin     (0.10);
        //
        auto    hTemporary_1D   =   ( THXGenRec_Type1D* )( hEff_1D->Clone() );
        hTemporary_1D = uScale( hTemporary_1D, hEff_1D->GetBinContent(xBin), hEff_1D->GetBinError(xBin) );
        auto    hTemporary_2D   =   hEff_2D->ProjectionX( Form("tmp_%i_%i",xBin,iBuilderTH1_TypeCounter), xBin, xBin );
        //
        auto    kMaximum    =   max ( hTemporary_1D->GetMaximum(), hTemporary_2D->GetMaximum() );
        auto    kMinimum    =   max ( hTemporary_1D->GetMinimum(), hTemporary_2D->GetMinimum() );
        //
        if ( !kSigLoss ) {
            uSetHisto( hTemporary_1D, "EFF  12D " );
            uSetHisto( hTemporary_2D, "EFF2 12D " );
        } else {
            uSetHisto( hTemporary_1D, "EFF  SL  12D " );
            uSetHisto( hTemporary_2D, "EFF2 SL  12D " );
        }
        //
        hTemporary_1D->SetMaximum(160*kMaximum);
        if ( kSigLoss ) hTemporary_1D->SetMinimum(0 - kMinimum);
        hTemporary_1D->GetXaxis()->SetTitleSize(0.075);
        hTemporary_1D->GetYaxis()->SetTitleSize(0.062);
        hTemporary_1D->GetYaxis()->SetTitleOffset(0.75);
        hTemporary_1D->Draw("SAME");
        hTemporary_2D->Draw("SAME");
        //
        if ( xBin == 1 )   {
            lDrawLegend->AddEntry( hTemporary_1D, "#varepsilon_{1Dx1D}",    "EP" );
            lDrawLegend->AddEntry( hTemporary_2D, "#varepsilon_{2D}",       "EP" );
            lDrawLegend->Draw("same");
        }
        uLatex      ->  SetTextSize(0.075);
        uLatex      ->  DrawLatexNDC(0.18,0.07,Form("#it{p}_{T,#phi_{2}} (GeV/#it{c}) #in [%.1f;%.1f]",hTemporary_1D->GetBinLowEdge(xBin),hTemporary_1D->GetBinLowEdge(xBin+1)));
    }
    //
    iBuilderTH1_TypeCounter++;
    //
    return cDrawEfficiencyCompare;
}
//
#endif /* AAU_Efficiency_h */
