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
 ( std::vector<TH1VecType*> hInput, TH1Template* hTemplate, TString kFolder = "" ) {
    std::vector<TH1Template*>   fResult;
    auto    kNSigmas    =   6;
    for ( Int_t iTer = 0; iTer < kNSigmas; iTer++ )    {
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
    if ( !kFolder.IsNull() )    {
        gStyle->SetPadTopMargin     (0.15);
        gStyle->SetPadBottomMargin  (0.20);
        gStyle->SetPadRightMargin   (0.15);
        gStyle->SetPadLeftMargin    (0.15);
    }
    TCanvas*    cDrawFitRMS =   new TCanvas("cDrawRMS","",9000,6000);
    TLine*      kUtilLine   =   new TLine();
    if ( !kFolder.IsNull() )    {
        cDrawFitRMS ->  Divide(3,2);
        cDrawFitRMS ->  cd(1);
    }
    auto iHist = 1;
    for ( auto kHisto : hInput ) {
        //
        //  --- Recovering Full Histogram Mean and STDV
        auto    kFull_Mean  =   kHisto  ->  GetMean();
        auto    kFull_STDV  =   kHisto  ->  GetRMS();
        auto    kFull_MAX_  =   kHisto  ->  GetXaxis()  ->  GetXmax();
        auto    kFull_MIN_  =   kHisto  ->  GetXaxis()  ->  GetXmin();
        if ( !kFolder.IsNull() ) {
            //kHisto  ->  SetTitle("");
            kHisto  ->  GetXaxis()  ->  SetTitle("m_{rec} - m_{gen} (GeV/c^{2})");
            kHisto  ->  SetMarkerStyle  ( uGetMarker(2) );
            kHisto  ->  SetMarkerColor  ( uGetColor(2) );
            kHisto  ->  SetLineColor    ( uGetColor(0) );
            kHisto  ->  SetMinimum      ( 0 );
            kHisto  ->  SetMaximum      ( 1.4 * kHisto    ->  GetMaximum() );
        }
        //
        //  --- Restricting in +- Nsigma
        for ( Int_t iTer = 0; iTer < kNSigmas; iTer++ ) {
            switch ( TEvalMethod )  {
                default:
                    kHisto  ->  GetXaxis()  ->  SetRangeUser( kFull_Mean - ( iTer + 1 )*kFull_STDV, kFull_Mean + ( iTer + 1 )*kFull_STDV  );
                    fResult.at( iTer )  ->  SetBinContent   ( iHist, kHisto  ->  GetRMS()       * ( kGaussStndDevtScale[iTer] ) );
                    fResult.at( iTer )  ->  SetBinError     ( iHist, kHisto  ->  GetRMSError()  * ( kGaussStndDevtScale[iTer] ) );
                    if ( !kFolder.IsNull() ) {
                        cDrawFitRMS->cd(iTer+1);
                        kHisto->Draw("PE1 MIN0");
                        //
                        //  --- Mean
                        auto    kXPos_mean  =   0.5 + 0.6*( kFull_Mean/(kFull_MAX_-kFull_MIN_) );
                        kUtilLine   ->  SetLineWidth( 2 );
                        kUtilLine   ->  SetLineStyle( 1 );
                        kUtilLine   ->  SetLineColor( uGetColor(1) );
                        kUtilLine   ->  DrawLineNDC(kXPos_mean,0.20,kXPos_mean,0.70);
                        //
                        //  --- Stdv
                        auto    kXPos_plus  =   kXPos_mean + 0.6*( ( iTer + 1 )*kFull_STDV/(kFull_MAX_-kFull_MIN_) );
                        auto    kXPos_mnus  =   kXPos_mean - 0.6*( ( iTer + 1 )*kFull_STDV/(kFull_MAX_-kFull_MIN_) );
                        kUtilLine   ->  SetLineWidth( 2 );
                        kUtilLine   ->  SetLineStyle( 2 );
                        kUtilLine   ->  SetLineColor( uGetColor(1) );
                        kUtilLine   ->  DrawLineNDC(kXPos_plus,0.20,kXPos_plus,0.51);
                        kUtilLine   ->  DrawLineNDC(kXPos_mnus,0.20,kXPos_mnus,0.51);
                        //
                        uLatex      ->  DrawLatexNDC( 0.2, 0.8, Form("#sigma_{res} = %5.5f", kHisto  ->  GetRMS()  * ( kGaussStndDevtScale[iTer] ) ) );
                    }
                    break;
                case 1:
                    fGauss  ->  SetParameter( 1, kFull_Mean );
                    fGauss  ->  SetParameter( 2, kFull_STDV );
                    kHisto  ->  Fit( fGauss, "IMEQN0", "", kFull_Mean - ( iTer + 1 )*kFull_STDV, kFull_Mean + ( iTer + 1 )*kFull_STDV );
                    fResult.at( iTer )  ->  SetBinContent   ( iHist, fGauss -> GetParameter( 2 ) );
                    fResult.at( iTer )  ->  SetBinError     ( iHist, fGauss -> GetParError ( 2 ) );
                    if ( !kFolder.IsNull() ) {
                        cDrawFitRMS->cd(iTer+1);
                        kHisto->Draw("PE1 MIN0");
                        fGauss->SetRange( kFull_Mean - ( iTer + 1 )*kFull_STDV, kFull_Mean + ( iTer + 1 )*kFull_STDV );
                        fGauss->DrawCopy("same");
                        //
                        //  --- Mean
                        auto    kXPos_mean  =   0.5 + 0.6*( fGauss -> GetParameter( 1 )/(kFull_MAX_-kFull_MIN_) );
                        kUtilLine   ->  SetLineWidth( 2 );
                        kUtilLine   ->  SetLineStyle( 1 );
                        kUtilLine   ->  SetLineColor( uGetColor(1) );
                        kUtilLine   ->  DrawLineNDC(kXPos_mean,0.20,kXPos_mean,0.70);
                        //
                        //  --- Stdv
                        auto    kXPos_plus  =   kXPos_mean + 0.6*( fGauss -> GetParameter( 2 )/(kFull_MAX_-kFull_MIN_) );
                        auto    kXPos_mnus  =   kXPos_mean - 0.6*( fGauss -> GetParameter( 2 )/(kFull_MAX_-kFull_MIN_) );
                        kUtilLine   ->  SetLineWidth( 2 );
                        kUtilLine   ->  SetLineStyle( 2 );
                        kUtilLine   ->  SetLineColor( uGetColor(1) );
                        kUtilLine   ->  DrawLineNDC(kXPos_plus,0.20,kXPos_plus,0.51);
                        kUtilLine   ->  DrawLineNDC(kXPos_mnus,0.20,kXPos_mnus,0.51);
                        //
                        uLatex      ->  DrawLatexNDC( 0.2, 0.8, Form("#sigma_{res} = %5.5f", fGauss -> GetParameter( 2 ) ) );
                    }
                    break;
            }
        }
        kHisto  ->  GetXaxis()  ->  SetRangeUser( kFull_Mean - ( 15 )*kFull_STDV, kFull_Mean + ( 15 )*kFull_STDV  );
        //
        if ( !kFolder.IsNull() )    {
            if ( TEvalMethod == 0 ) cDrawFitRMS ->  SaveAs(kFolder+TString(Form("Res_RMS_Check_%i.pdf",iHist)));
            if ( TEvalMethod == 1 ) cDrawFitRMS ->  SaveAs(kFolder+TString(Form("Res_Fit_Check_%i.pdf",iHist)));
        }
        //
        iHist++;
    }
    delete  cDrawFitRMS;
    return fResult;
}
//
template < typename TH1VecType_1, typename TH1VecType_2, typename TH1Template >
TH1Template*
uCalculateResolutionTrueMassFIT
( std::vector<TH1VecType_1*> hRecInvMass, std::vector<TH1VecType_1*> hTrueInvMass, std::vector<TH1VecType_2*> hDeltaMass, TH1Template* hTemplate, std::tuple<Float_t,Float_t,Float_t> kMass, std::tuple<Float_t,Float_t,Float_t> kWidth, Float_t kSigma, TString kFolder = "" )   {
    auto    fResult =   (TH1Template*)( hTemplate->Clone() );
    auto    iHist   =   1;
    for ( auto kHisto : hRecInvMass ) {
        //
        //  --- Setting up the Fit
        gROOT   ->  SetBatch( kTRUE );
        //
        auto            hCurrent_TrueInvMass    =   hTrueInvMass.at(iHist-1);
        //
        RooRealVar      InvMass =   RooRealVar      ("InvMass", "m_{REC}",          get<1>( kMass ),    get<2>( kMass ) );
        RooRealVar      InvMas2 =   RooRealVar      ("InvMas2", "m_{REC}",          1.0194,    1.0195 );
        RooDataHist*    dPft    =   new RooDataHist ("dPft",    "dPft",             InvMas2,            RooFit::Import(*hCurrent_TrueInvMass) );
        RooDataHist*    data    =   new RooDataHist ("Data",    "Data",             InvMass,            RooFit::Import(*kHisto) );
        //
        hCurrent_TrueInvMass->GetXaxis()->SetRangeUser(1.0194,    1.0195 );
        RooRealVar      sMass, sMPFt, sWPFt, sWidt, sSlop;
                        sMass   =   RooRealVar      ("bMass",   "bMass",    get<0>( kMass ),    get<1>( kMass ),    get<2>( kMass )     );
                        sMPFt   =   RooRealVar      ("bMPFt",   "bMPFt",    hCurrent_TrueInvMass->GetMean());
                        sWPFt   =   RooRealVar      ("bWPFt",   "bWPFt",    get<0>( kWidth ),   get<1>( kWidth ),   get<2>( kWidth )    );
                        sWidt   =   RooRealVar      ("bWidt",   "bWidt",    get<0>( kWidth )    );
                        sSlop   =   RooRealVar      ("bSlop",   "bSlop",    kSigma,             0.0,                1000.0              );
        RooVoigtian     fSig    =   RooVoigtian     ("fSig",    "fSig",     InvMass,            sMass,              sWidt,              sSlop);
        RooBreitWigner  fPFt    =   RooBreitWigner  ("fPFt",    "fPFt",     InvMas2,            sMPFt,              sWPFt);
        //
        auto fFitResults    =   fPFt.fitTo( *data, RooFit::Save() );
        auto N_Raw          =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().find("bWPFt"));
        //
        auto RMS_Width      =   hCurrent_TrueInvMass->GetMean();
        //cout << N_Raw->getVal() << endl;
        //sWidt               .   setVal      ( N_Raw->getVal() );
        //
        TCanvas *cDrawPlot  =   new TCanvas();
        auto fSaveToFrame   =   InvMas2.frame( RooFit::Name(""), RooFit::Title("") );
        dPft                ->  plotOn  ( fSaveToFrame, RooFit::MarkerColor(38),    RooFit::MarkerStyle(26),    RooFit::Name("RooData") );
        fPFt                .   plotOn  ( fSaveToFrame, RooFit::LineColor(4),       RooFit::LineStyle(kSolid),  RooFit::Name("RooPFt")  );
        fSaveToFrame        ->  SetTitle("");
        fSaveToFrame        ->  Draw();
        uLatex              ->  DrawLatexNDC( 0.2, 0.8, Form("#Gamma = %5.5f", 1.e3*N_Raw -> getVal()  ) );
        cDrawPlot           ->  SaveAs( kFolder+TString( Form( "/Res_TFR_True_Check_%i.pdf", iHist ) ) );
        delete              fSaveToFrame;
        delete              cDrawPlot;
        //
        sMass                   .   setVal      ( get<0>( kMass ) - hDeltaMass.at(iHist-1)->GetMean() );
        //
        fFitResults = fSig.fitTo( *data, RooFit::Save(), RooFit::InitialHesse( kTRUE ), RooFit::Minos( kTRUE ) );
        N_Raw  =   static_cast<RooRealVar*>(fFitResults ->floatParsFinal().find("bSlop"));
        //
        fResult->SetBinContent  ( iHist, N_Raw -> getVal() );
        fResult->SetBinError    ( iHist, N_Raw -> getError() );
        //
        cDrawPlot  =   new TCanvas();
        fSaveToFrame        =   InvMass.frame( RooFit::Name(""), RooFit::Title("") );
        data                ->  plotOn  ( fSaveToFrame, RooFit::MarkerColor(38),    RooFit::MarkerStyle(26),    RooFit::Name("RooData") );
        fSig                .   plotOn  ( fSaveToFrame, RooFit::LineColor(4),       RooFit::LineStyle(kSolid),  RooFit::Name("RooMod")  );
        fSaveToFrame        ->  SetTitle("");
        fSaveToFrame        ->  Draw();
        uLatex              ->  DrawLatexNDC( 0.2, 0.8, Form("#sigma_{res} = %5.5f", N_Raw -> getVal()  ) );
        cDrawPlot           ->  SaveAs( kFolder+TString( Form( "/Res_TFR_Check_%i.pdf", iHist ) ) );
        delete              fSaveToFrame;
        delete              cDrawPlot;
        //
        iHist++;
    }
    return fResult;
}
//
// TODO: Add check plots with bands where RMS +- nSigma is located, w/ superimpoxition of gauss fit resutls.
template < typename TH1VecType_1, typename TH1VecType_2, typename TH1Template >
std::vector<TH1F*>
uCalculateResolution
( std::vector<TH1VecType_1*> hDeltaMass, std::vector<TH1VecType_2*> hRecInvMass, std::vector<TH1VecType_2*> hTrueInvMass, TH1Template* hTemplate, std::tuple<Float_t,Float_t,Float_t> kMass, std::tuple<Float_t,Float_t,Float_t> kWidth, Float_t kSigma, TString kFolder = ""  ) {
    auto    fResult     =   uCalculateResolution<0>( hDeltaMass, hTemplate, kFolder );
    auto    fAppendVect =   uCalculateResolution<1>( hDeltaMass, hTemplate, kFolder );
    auto    fAppendHist =   uCalculateResolutionTrueMassFIT( hRecInvMass, hTrueInvMass, hDeltaMass, hTemplate, kMass, kWidth, kSigma, kFolder );
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
    TLegend*        lFullLegend     =   new TLegend(0.20,0.85,0.82,0.71);
    lFullLegend ->  SetNColumns(5);
    lFullLegend ->  SetFillColorAlpha(0.,0.);
    lFullLegend ->  SetLineColorAlpha(0.,0.);
    //
    for ( Int_t iTer = 0; iTer < 6; iTer++ )    {
        hResolutionPlots.at( iTer ) ->  SetMaximum( 1.50*hResolutionPlots.at( iTer ) ->GetMaximum() );
        hResolutionPlots.at( iTer ) ->  SetMinimum( 0.50*hResolutionPlots.at( iTer ) ->GetMinimum() );
        hResolutionPlots.at( iTer ) ->  SetLineColor  ( uGetColor(2) );
        hResolutionPlots.at( iTer ) ->  SetMarkerColor( uGetColor(2) );
        hResolutionPlots.at( iTer ) ->  SetMarkerStyle( uGetMarker(iTer) );
        if ( iTer+1 > kResolutionLegend.size() )            lFullLegend ->  AddEntry( hResolutionPlots.at(iTer), hResolutionPlots.at(iTer)->GetName(),   "EP" );
        else if ( !kResolutionLegend.at(iTer).IsNull() )    lFullLegend ->  AddEntry( hResolutionPlots.at(iTer), kResolutionLegend.at(iTer),             "EP" );
        else                                                lFullLegend ->  AddEntry( hResolutionPlots.at(iTer), hResolutionPlots.at(iTer)->GetName(),   "EP" );
        hResolutionPlots.at( iTer ) ->  Draw("SAME EP");
    }
    for ( Int_t iTer = 6; iTer < 12; iTer++ )    {
        hResolutionPlots.at( iTer ) ->  SetLineColor( uGetColor(3) );
        hResolutionPlots.at( iTer ) ->  SetMarkerColor( uGetColor(3) );
        hResolutionPlots.at( iTer ) ->  SetMarkerStyle( uGetMarker(iTer-5) );
        if ( iTer+1 > kResolutionLegend.size() )            lFullLegend ->  AddEntry( hResolutionPlots.at(iTer), hResolutionPlots.at(iTer)->GetName(),   "EP" );
        else if ( !kResolutionLegend.at(iTer).IsNull() )    lFullLegend ->  AddEntry( hResolutionPlots.at(iTer), kResolutionLegend.at(iTer),             "EP" );
        else                                                lFullLegend ->  AddEntry( hResolutionPlots.at(iTer), hResolutionPlots.at(iTer)->GetName(),   "EP" );
        hResolutionPlots.at( iTer ) ->  Draw("SAME EP");
    }
    hResolutionPlots.at( 12 ) ->  SetLineColor( uGetColor(4) );
    hResolutionPlots.at( 12 ) ->  SetMarkerColor( uGetColor(4) );
    hResolutionPlots.at( 12 ) ->  SetMarkerStyle( uGetMarker(0) );
    if ( 13 > kResolutionLegend.size() )                lFullLegend ->  AddEntry( hResolutionPlots.at(12), hResolutionPlots.at(12)->GetName(),   "EP" );
    else if ( !kResolutionLegend.at(12).IsNull() )      lFullLegend ->  AddEntry( hResolutionPlots.at(12), kResolutionLegend.at(12),             "EP" );
    else                                                lFullLegend ->  AddEntry( hResolutionPlots.at(12), hResolutionPlots.at(12)->GetName(),   "EP" );
    hResolutionPlots.at( 12 ) ->  Draw("SAME EP");
    lFullLegend ->  Draw("SAME");
    fResult.push_back( cFullEfficiencies );
    //
    TCanvas*    cPartEfficiencies    =   new TCanvas("cPartEfficiencies","",1200,1200);
    gStyle  ->  SetOptStat(0);
    gPad    ->  SetLogx();
    gPad    ->  SetGridy();
    //
    TLegend*        lPartLegend     =   new TLegend(0.20,0.85,0.82,0.71);
    lPartLegend ->  SetNColumns(3);
    lPartLegend ->  SetFillColorAlpha(0.,0.);
    lPartLegend ->  SetLineColorAlpha(0.,0.);
    //
    if ( 3 > kResolutionLegend.size() )                 lPartLegend ->  AddEntry( hResolutionPlots.at(2), hResolutionPlots.at(2)->GetName(),   "EP" );
    else if ( !kResolutionLegend.at(2).IsNull() )       lPartLegend ->  AddEntry( hResolutionPlots.at(2), kResolutionLegend.at(2),             "EP" );
    else                                                lPartLegend ->  AddEntry( hResolutionPlots.at(2), hResolutionPlots.at(2)->GetName(),   "EP" );
    hResolutionPlots.at( 2 ) ->  Draw("SAME EP");
    //
    if ( 8 > kResolutionLegend.size() )                 lPartLegend ->  AddEntry( hResolutionPlots.at(7), hResolutionPlots.at(7)->GetName(),   "EP" );
    else if ( !kResolutionLegend.at(7).IsNull() )       lPartLegend ->  AddEntry( hResolutionPlots.at(7), kResolutionLegend.at(7),             "EP" );
    else                                                lPartLegend ->  AddEntry( hResolutionPlots.at(7), hResolutionPlots.at(7)->GetName(),   "EP" );
    hResolutionPlots.at( 7 ) ->  Draw("SAME EP");
    //
    if ( 13 > kResolutionLegend.size() )                lPartLegend ->  AddEntry( hResolutionPlots.at(12), hResolutionPlots.at(12)->GetName(),   "EP" );
    else if ( !kResolutionLegend.at(12).IsNull() )      lPartLegend ->  AddEntry( hResolutionPlots.at(12), kResolutionLegend.at(12),             "EP" );
    else                                                lPartLegend ->  AddEntry( hResolutionPlots.at(12), hResolutionPlots.at(12)->GetName(),   "EP" );
    hResolutionPlots.at( 12 ) ->  Draw("SAME EP");
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
    kUtility.push_back( hResolutionPlots.at( 7  ) );
    kUtility.push_back( hResolutionPlots.at( 12 ) );
    auto kDrawUtility = uMakeRatio( kUtility );
    //
    kDrawUtility.at(0)  ->  SetMaximum( 1.75*kDrawUtility.at(0)->GetMaximum() );
    kDrawUtility.at(0)  ->  SetMinimum( 0.75*kDrawUtility.at(0)->GetMinimum() );
    kDrawUtility.at(0)  ->  Draw("SAME EP");
    kDrawUtility.at(1)  ->  Draw("SAME EP");
    kDrawUtility.at(2)  ->  Draw("SAME EP");
    lPartLegend         ->  Draw("SAME");
    //
    fResult.push_back( cPartEfficienciesNorm );
    //
    return fResult;
}
//
#endif /* AAU_Resolution_h */
