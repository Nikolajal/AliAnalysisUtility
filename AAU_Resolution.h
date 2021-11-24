//
//  Part of the AliAnalysisUtility package
//
//  Utilities for evaluation of detector resolution in HEP analysis
//
//  Author              Nicola Rubini
//  Created             22/11/2021
//  Last modified       22/11/2021
#ifndef AAU_Resolution_h
#define AAU_Resolution_h
//
//  Global File w/ Constants and other functions
#include "AliAnalysisUtility.h"
//
//  --- RESONANCES BASED CALCULATION
//
template < Int_t TEvalMethod = 0, typename TH1VecType, typename TH1Template >
std::vector<TH1Template*>
uCalculateResolution
 ( std::vector<TH1VecType*> hInput, TH1Template* hTemplate ) {
    std::vector<TH1Template*>   fResult;
    for ( Int_t iTer = 0; iTer < 5; iTer++ )    {
        auto    hCurrent_Resolution =   (TH1Template*)(hTemplate->Clone());
        switch ( TEvalMethod )  {
            default:
                hCurrent_Resolution->SetName(Form("hRes_RMS_%i",iTer+1));
                break;
            case 1:
                hCurrent_Resolution->SetName(Form("hRes_GAU_%i",iTer+1));
                break;
        }
        fResult.push_back( hCurrent_Resolution );
    }
    auto iHist = 1;
    for ( auto kHisto : hInput ) {
        //
        //  --- Recovering Full Histogram Mean and STDV
        auto    kFull_Mean  =   kHisto->GetMean();
        auto    kFull_STDV  =   kHisto->GetRMS();
        //
        //  --- Restricting in +- Nsigma
        for ( Int_t iTer = 0; iTer < 5; iTer++ ) {
            switch ( TEvalMethod )  {
                default:
                    kHisto  ->  GetXaxis()  ->  SetRangeUser( kFull_Mean - ( iTer + 1 )*kFull_STDV, kFull_Mean + ( iTer + 1 )*kFull_STDV  );
                    fResult.at( iTer )  ->  SetBinContent   ( iHist, kHisto  ->  GetRMS()       * ( kGaussStndDevtScale[iTer] ) );
                    fResult.at( iTer )  ->  SetBinError     ( iHist, kHisto  ->  GetRMSError()  * ( kGaussStndDevtScale[iTer] ) );
                    break;
                case 1:
                    fGauss  ->  SetParameter( 1, kFull_Mean );
                    fGauss  ->  SetParameter( 2, kFull_STDV );
                    kHisto  ->  Fit( fGauss, "IMESQ", "R", kFull_Mean - ( iTer + 1 )*kFull_STDV, kFull_Mean + ( iTer + 1 )*kFull_STDV );
                    fResult.at( iTer )  ->  SetBinContent   ( iHist, fGauss -> GetParameter( 2 ) );
                    fResult.at( iTer )  ->  SetBinError     ( iHist, fGauss -> GetParError ( 2 ) );
                    break;
            }
        }
        kHisto  ->  GetXaxis()  ->  SetRangeUser( kFull_Mean - ( 15 )*kFull_STDV, kFull_Mean + ( 15 )*kFull_STDV  );
        //
        iHist++;
    }
    return fResult;
}
//
template < typename TH1VecType_1, typename TH1Template >
TH1Template*
uCalculateResolutionTrueMassFIT
 ( std::vector<TH1VecType_1*> hInput, TH1Template* hFill  );
//
// TODO: Add check plots with bands where RMS +- nSigma is located, w/ superimpoxition of gauss fit resutls.
template < typename TH1VecType_1, typename TH1VecType_2, typename TH1Template >
std::vector<TH1F*>
uCalculateResolution
( std::vector<TH1VecType_1*> hDeltaMass, std::vector<TH1VecType_2*> hTrueInvMass, TH1Template* hTemplate  ) {
    auto    fResult     =   uCalculateResolution<0>( hDeltaMass, hTemplate );
    auto    fAppendVect =   uCalculateResolution<1>( hDeltaMass, hTemplate );
    auto    fAppendHist =   uCalculateResolutionTrueMassFIT( hTrueInvMass, hTemplate );
    for ( auto kHist : fAppendVect ) fResult.push_back( kHist );
    fResult.push_back( fAppendHist );
    return fResult;
}
//
template < typename TH1VecType_1 >
std::vector<TCanvas*>
uPlotResolution
 ( std::vector<TH1VecType_1*> hResolutionPlots ) {
    std::vector<TCanvas*>   fResult;
    TCanvas*    cFullEfficiencies    =   new TCanvas("cFullEfficiencies","",1200,1200);
    gStyle  ->  SetOptStat(0);
    gPad    ->  SetLogx();
    gPad    ->  SetGridy();
    //
    TLegend*        lFullLegend     =   new TLegend(0.625,0.88,0.88,0.7);
    lFullLegend ->  SetNColumns(2);
    lFullLegend ->  SetFillColorAlpha(0.,0.);
    lFullLegend ->  SetLineColorAlpha(0.,0.);
    //
    for ( Int_t iTer = 0; iTer < 5; iTer++ )    {
        hResolutionPlots.at( iTer ) ->  SetMaximum  ( hResolutionPlots.at( iTer )->GetMaximum() * 1.4 );
        hResolutionPlots.at( iTer ) ->  SetMinimum  ( 0 );
        hResolutionPlots.at( iTer ) ->  SetLineColor( uGetColor(2) );
        hResolutionPlots.at( iTer ) ->  SetMarkerColor( uGetColor(2) );
        hResolutionPlots.at( iTer ) ->  SetMarkerStyle( uGetMarker(iTer) );
        if ( iTer+1 > kResolutionLegend.size() )            lFullLegend ->  AddEntry( hResolutionPlots.at(iTer), hResolutionPlots.at(iTer)->GetName(),   "EP" );
        else if ( !kResolutionLegend.at(iTer).IsNull() )    lFullLegend ->  AddEntry( hResolutionPlots.at(iTer), kResolutionLegend.at(iTer),             "EP" );
        else                                                lFullLegend ->  AddEntry( hResolutionPlots.at(iTer), hResolutionPlots.at(iTer)->GetName(),   "EP" );
        hResolutionPlots.at( iTer ) ->  Draw("SAME EP");
    }
    for ( Int_t iTer = 5; iTer < 10; iTer++ )    {
        hResolutionPlots.at( iTer ) ->  SetLineColor( uGetColor(3) );
        hResolutionPlots.at( iTer ) ->  SetMarkerColor( uGetColor(3) );
        hResolutionPlots.at( iTer ) ->  SetMarkerStyle( uGetMarker(iTer-5) );
        if ( iTer+1 > kResolutionLegend.size() )            lFullLegend ->  AddEntry( hResolutionPlots.at(iTer), hResolutionPlots.at(iTer)->GetName(),   "EP" );
        else if ( !kResolutionLegend.at(iTer).IsNull() )    lFullLegend ->  AddEntry( hResolutionPlots.at(iTer), kResolutionLegend.at(iTer),             "EP" );
        else                                                lFullLegend ->  AddEntry( hResolutionPlots.at(iTer), hResolutionPlots.at(iTer)->GetName(),   "EP" );
        hResolutionPlots.at( iTer ) ->  Draw("SAME EP");
    }
    hResolutionPlots.at( 10 ) ->  SetLineColor( uGetColor(4) );
    hResolutionPlots.at( 10 ) ->  SetMarkerColor( uGetColor(4) );
    hResolutionPlots.at( 10 ) ->  SetMarkerStyle( uGetMarker(0) );
    if ( 11 > kResolutionLegend.size() )                lFullLegend ->  AddEntry( hResolutionPlots.at(10), hResolutionPlots.at(10)->GetName(),   "EP" );
    else if ( !kResolutionLegend.at(10).IsNull() )      lFullLegend ->  AddEntry( hResolutionPlots.at(10), kResolutionLegend.at(10),             "EP" );
    else                                                lFullLegend ->  AddEntry( hResolutionPlots.at(10), hResolutionPlots.at(10)->GetName(),   "EP" );
    hResolutionPlots.at( 10 ) ->  Draw("SAME EP");
    lFullLegend ->  Draw("SAME");
    fResult.push_back( cFullEfficiencies );
    //
    TCanvas*    cPartEfficiencies    =   new TCanvas("cPartEfficiencies","",1200,1200);
    gStyle  ->  SetOptStat(0);
    gPad    ->  SetLogx();
    gPad    ->  SetGridy();
    //
    TLegend*        lPartLegend     =   new TLegend(0.625,0.88,0.88,0.7);
    lPartLegend ->  SetNColumns(2);
    lPartLegend ->  SetFillColorAlpha(0.,0.);
    lPartLegend ->  SetLineColorAlpha(0.,0.);
    //
    if ( 3 > kResolutionLegend.size() )                 lPartLegend ->  AddEntry( hResolutionPlots.at(2), hResolutionPlots.at(2)->GetName(),   "EP" );
    else if ( !kResolutionLegend.at(2).IsNull() )       lPartLegend ->  AddEntry( hResolutionPlots.at(2), kResolutionLegend.at(2),             "EP" );
    else                                                lPartLegend ->  AddEntry( hResolutionPlots.at(2), hResolutionPlots.at(2)->GetName(),   "EP" );
    hResolutionPlots.at( 2 ) ->  Draw("SAME EP");
    //
    if ( 7 > kResolutionLegend.size() )                 lPartLegend ->  AddEntry( hResolutionPlots.at(6), hResolutionPlots.at(6)->GetName(),   "EP" );
    else if ( !kResolutionLegend.at(6).IsNull() )       lPartLegend ->  AddEntry( hResolutionPlots.at(6), kResolutionLegend.at(6),             "EP" );
    else                                                lPartLegend ->  AddEntry( hResolutionPlots.at(6), hResolutionPlots.at(6)->GetName(),   "EP" );
    hResolutionPlots.at( 6 ) ->  Draw("SAME EP");
    //
    if ( 11 > kResolutionLegend.size() )                lPartLegend ->  AddEntry( hResolutionPlots.at(10), hResolutionPlots.at(10)->GetName(),   "EP" );
    else if ( !kResolutionLegend.at(10).IsNull() )      lPartLegend ->  AddEntry( hResolutionPlots.at(10), kResolutionLegend.at(10),             "EP" );
    else                                                lPartLegend ->  AddEntry( hResolutionPlots.at(10), hResolutionPlots.at(10)->GetName(),   "EP" );
    hResolutionPlots.at( 10 ) ->  Draw("SAME EP");
    lPartLegend ->  Draw("SAME");
    //
    fResult.push_back( cPartEfficiencies );
    //
    TCanvas*    cPartEfficienciesNorm    =   new TCanvas("cPartEfficienciesNorm","",1200,1200);
    gStyle  ->  SetOptStat(0);
    gPad    ->  SetLogx();
    gPad    ->  SetGridy();
    //
    std::vector<TH1F*> kUtility;
    kUtility.push_back( hResolutionPlots.at( 2  ) );
    kUtility.push_back( hResolutionPlots.at( 6  ) );
    kUtility.push_back( hResolutionPlots.at( 10 ) );
    auto kDrawUtility = uMakeRatio( kUtility );
    //
    kDrawUtility.at(0)  ->  SetMaximum(1.4);
    kDrawUtility.at(0)  ->  SetMinimum(0.6);
    kDrawUtility.at(0)  ->  Draw("SAME EP");
    kDrawUtility.at(1)  ->  Draw("SAME EP");
    kDrawUtility.at(2)  ->  Draw("SAME EP");
    lPartLegend     ->  Draw("SAME");
    //
    fResult.push_back( cPartEfficienciesNorm );
    //
    return fResult;
}
//
#endif /* AAU_Resolution_h */
