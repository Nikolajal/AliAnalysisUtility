//
//  Part of the AliAnalysisUtility package
//
//  Utilities for evaluation of detector resolution in HEP analysis
//
//  Author              Nicola Rubini
//  Created             01/12/2021
//  Last modified       01/12/2021
#ifndef AAU_Extrapolation_h
#define AAU_Extrapolation_h
//
//  Global File w/ Constants and other functions
#include "AliAnalysisUtility.h"
//
template <  typename THXTarget_Type >
std::map<TString,std::tuple<Float_t,Float_t,Float_t>>
uMeasureFullYield
( THXTarget_Type* hTarget_stat, THXTarget_Type* hTarget_syst, std::vector<std::tuple<TF1*,Float_t,Float_t,TString>> fExtrapolationFunc_Array, Int_t kIterations, TString kSaveFolder = "", TString kSaveName = "", Float_t kScale = 1. ) {
    std::map<TString,std::tuple<Float_t,Float_t,Float_t>>   fResult;
    //
    //  --- Implementation for TH1 of 1D
    //  TODO: Get dimension as specxialisation variable and overload (?)
    if ( !uIsTHPairConsistent( hTarget_stat, hTarget_syst ) )   { cout << "Stat and Syst Err Histos must be consistent! " << endl; return fResult; }
    //  TODO: Implement this (rvalue)
    //if ( uGetTHDimension( hTarget_stat ) == 2 )                 return uMeasureFullYield2D( hTarget_stat, hTarget_syst, fExtrapolationFunc_Array, kIterations, kSaveFolder, kSaveName );
    if ( uGetTHDimension( hTarget_stat ) != 1 )                 return fResult;
    //
    //  --- Optimisation mode
    gROOT->SetBatch(kTRUE);
    //
    //  --- --- Integral Definition w/ errors
    auto        fExtrapolationFunc  =   fExtrapolationFunc_Array.at(0);
    Double_t    kIntegral_Yield_stat, kIntegral_Yield_syst, kIntegral_TMom_stat, kIntegral_TMom_syst, kExtrapol_Yield_stat, kExtrapol_Yield_syst, kExtrapol_TMom_stat, kExtrapol_TMom_syst;
    auto        kIntegral_Yield     =   hTarget_syst->IntegralAndError(-1,-1,kIntegral_Yield_syst,"width");
                kIntegral_Yield     =   hTarget_stat->IntegralAndError(-1,-1,kIntegral_Yield_stat,"width");
                kIntegral_TMom_stat =   uEvaluateMeanPTError( hTarget_stat );
                kIntegral_TMom_syst =   uEvaluateMeanPTError( hTarget_syst );
    //
    auto        hFullUncertainties  =   uSumErrors( hTarget_stat, hTarget_syst );
                hFullUncertainties  ->  Fit( get<0>(fExtrapolationFunc), get<3>(fExtrapolationFunc), "", get<1>(fExtrapolationFunc), get<2>(fExtrapolationFunc) );
    auto        kExtrapol_Yield     =   get<0>(fExtrapolationFunc)->Integral( 0., hFullUncertainties->GetBinLowEdge(1) );
    auto        kFull_TMom          =   uEvaluateMeanPT( hFullUncertainties, get<0>(fExtrapolationFunc) );
    //
    //  --- --- --- Assigning Integral Results
    fResult["YL_INT"] = { kScale*kIntegral_Yield,   kScale*kIntegral_Yield_stat,    kScale*kIntegral_Yield_syst };
    fResult["PT_INT"] = { kFull_TMom,               kIntegral_TMom_stat,            kIntegral_TMom_syst         };
    //
    //  TODO: Check this can be taken care of by an external function
    if ( !kSaveFolder.IsNull() ) {
        TCanvas*    cSaveMeasureFit =   new TCanvas( "cSaveMeasureFit", "cSaveMeasureFit", 2000, 600 );
        cSaveMeasureFit             ->  Divide(3,1);
        //
        hFullUncertainties          ->  SetTitle("");
        hFullUncertainties          ->  SetMarkerStyle  ( uGetMarker(0) );
        hFullUncertainties          ->  SetMarkerColor  ( uGetColor (2) );
        hFullUncertainties          ->  SetLineColor    ( uGetColor (2) );
        get<0>(fExtrapolationFunc)  ->  SetLineColor    ( uGetColor (1) );
        get<0>(fExtrapolationFunc)  ->  SetLineWidth    ( 2             );
        //
        cSaveMeasureFit             ->  cd(1);
        gPad                        ->  SetLogy();
        hFullUncertainties          ->  DrawCopy();
        get<0>(fExtrapolationFunc)  ->  DrawCopy("SAME");
        //
        auto    kLowLimit           =   get<1>(fExtrapolationFunc);
        auto    kHigLimit           =   get<1>(fExtrapolationFunc) + 0.15*( get<2>(fExtrapolationFunc) - get<1>(fExtrapolationFunc) );
        hFullUncertainties          ->  GetXaxis()  ->  SetRangeUser( kLowLimit, kHigLimit );
        cSaveMeasureFit             ->  cd(2);
        gPad                        ->  SetLogy();
        hFullUncertainties          ->  DrawCopy();
        get<0>(fExtrapolationFunc)  ->  DrawCopy("SAME");
        //
                kLowLimit           =   get<2>(fExtrapolationFunc) - 0.5*( get<2>(fExtrapolationFunc) - get<1>(fExtrapolationFunc) );
                kHigLimit           =   get<2>(fExtrapolationFunc);
        hFullUncertainties          ->  GetXaxis()  ->  SetRangeUser( kLowLimit, kHigLimit );
        cSaveMeasureFit             ->  cd(3);
        gPad                        ->  SetLogy();
        hFullUncertainties          ->  DrawCopy();
        get<0>(fExtrapolationFunc)  ->  DrawCopy("SAME");
        //
        cSaveMeasureFit  ->  SaveAs( ( kSaveFolder + TString("/FullYield") + kSaveName + TString(".pdf") ) );
        cSaveMeasureFit  ->  SaveAs( ( kSaveFolder + TString("/FullYield") + kSaveName + TString(".eps") ) );
        delete cSaveMeasureFit;
    }
    //
    auto    kStatEval   =   uExtrapolateToLowPT( hTarget_stat, hTarget_syst, fExtrapolationFunc, kIterations,   kSaveFolder, TString("StatEvaluation") + kSaveName );
    auto    kSystEval   =   uExtrapolateToLowPT( hTarget_syst, hTarget_stat, fExtrapolationFunc, kIterations,   kSaveFolder, TString("SystEvaluation") + kSaveName );
    auto    kFFitEval   =   uExtrapolateToLowPT( hTarget_syst, hTarget_stat, fExtrapolationFunc_Array,          kSaveFolder, TString("MFitEvaluation") + kSaveName );
    //
    kExtrapol_Yield_stat    =   kStatEval.first.second;
    kExtrapol_Yield_syst    =   SquareSum( { kSystEval.first.second, kFFitEval.first.second } );
    kExtrapol_TMom_stat     =   kStatEval.second.second;
    kExtrapol_TMom_syst     =   SquareSum( { kSystEval.second.second, kFFitEval.second.second } );
    //
    fResult["YL_EXT"] = { kScale*kExtrapol_Yield,     kScale*kExtrapol_Yield_stat,       kScale*kExtrapol_Yield_syst        };
    fResult["PT_EXT"] = { kFull_TMom,        kExtrapol_TMom_stat,    kExtrapol_TMom_syst    };
    fResult["YL_FLL"] = { get<0>(fResult["YL_INT"]) + get<0>(fResult["YL_EXT"]), SquareSum( { get<1>(fResult["YL_INT"]), get<1>(fResult["YL_EXT"]) } ), SquareSum( { get<2>(fResult["YL_INT"]), get<2>(fResult["YL_EXT"]) } ) };
    fResult["PT_FLL"] = { get<0>(fResult["PT_INT"]) + get<0>(fResult["PT_EXT"]), SquareSum( { get<1>(fResult["PT_INT"]), get<1>(fResult["PT_EXT"]) } ), SquareSum( { get<2>(fResult["PT_INT"]), get<2>(fResult["PT_EXT"]) } ) };
    //
    //  --- Optimisation off
    gROOT->SetBatch(kFALSE);
    //
    return fResult;
}
//
void
add_ext_results
( std::vector<std::map<TString,std::tuple<Float_t,Float_t,Float_t>>> &fTarget  )  {
    std::map<TString,std::tuple<Float_t,Float_t,Float_t>>   kToAdd;
    get<0>(kToAdd["YL_FLL"])    =   0.;
    get<1>(kToAdd["YL_FLL"])    =   0.;
    get<2>(kToAdd["YL_FLL"])    =   0.;
    get<0>(kToAdd["YL_INT"])    =   0.;
    get<1>(kToAdd["YL_INT"])    =   0.;
    get<2>(kToAdd["YL_INT"])    =   0.;
    get<0>(kToAdd["YL_EXT"])    =   0.;
    get<1>(kToAdd["YL_EXT"])    =   0.;
    get<2>(kToAdd["YL_EXT"])    =   0.;
    for ( auto kResult : fTarget )    {
        auto kBinValue              =   get<0>(kResult["YL_FLL"]);
        auto kBin_stat              =   get<1>(kResult["YL_FLL"]);
        auto kBin_syst              =   get<2>(kResult["YL_FLL"]);
        get<0>(kToAdd["YL_FLL"])   +=   kBinValue;
        get<1>(kToAdd["YL_FLL"])   +=   kBin_stat*kBin_stat;
        get<2>(kToAdd["YL_FLL"])   +=   kBin_stat*kBin_stat;
        kBinValue                   =   get<0>(kResult["YL_INT"]);
        kBin_stat                   =   get<1>(kResult["YL_INT"]);
        kBin_syst                   =   get<2>(kResult["YL_INT"]);
        get<0>(kToAdd["YL_INT"])   +=   kBinValue;
        get<1>(kToAdd["YL_INT"])   +=   kBin_stat*kBin_stat;
        get<2>(kToAdd["YL_INT"])   +=   kBin_stat*kBin_stat;
        kBinValue                   =   get<0>(kResult["YL_EXT"]);
        kBin_stat                   =   get<1>(kResult["YL_EXT"]);
        kBin_syst                   =   get<2>(kResult["YL_EXT"]);
        get<0>(kToAdd["YL_EXT"])   +=   kBinValue;
        get<1>(kToAdd["YL_EXT"])   +=   kBin_stat*kBin_stat;
        get<2>(kToAdd["YL_EXT"])   +=   kBin_stat*kBin_stat;
    }
    get<1>(kToAdd["YL_FLL"])    =   TMath::Sqrt( get<1>(kToAdd["YL_FLL"]) );
    get<2>(kToAdd["YL_FLL"])    =   TMath::Sqrt( get<2>(kToAdd["YL_FLL"]) );
    get<1>(kToAdd["YL_INT"])    =   TMath::Sqrt( get<1>(kToAdd["YL_INT"]) );
    get<2>(kToAdd["YL_INT"])    =   TMath::Sqrt( get<2>(kToAdd["YL_INT"]) );
    get<1>(kToAdd["YL_EXT"])    =   TMath::Sqrt( get<1>(kToAdd["YL_EXT"]) );
    get<2>(kToAdd["YL_EXT"])    =   TMath::Sqrt( get<2>(kToAdd["YL_EXT"]) );
    push_to_front( fTarget, kToAdd );
}
//
// TODO: Implement checks
template <  typename THXTarget_Type >
std::vector<std::map<TString,std::tuple<Float_t,Float_t,Float_t>>>
uMeasureFullYield2D
( THXTarget_Type* hTarget_stat, THXTarget_Type* hTarget_syst, std::vector<std::tuple<TF1*,Float_t,Float_t,TString>> fExtrapolationFunc_Array, Int_t kIterations, TString kSaveFolder = "", TString kSaveName = "", Float_t kScale = 2 ) {
    std::vector<std::map<TString,std::tuple<Float_t,Float_t,Float_t>>>   fResult;
    //
    //  --- Implementation for TH2 of 2D
    //
    std::vector<TH1D*>  kStat_Array;
    std::vector<TH1D*>  kSyst_Array;
    for ( Int_t iBin = 1; iBin <= hTarget_stat->GetNbinsX(); iBin++ )   kStat_Array.push_back( (TH1D*)(hTarget_stat->ProjectionX(Form("%s_%i",hTarget_stat->GetName(),iBin),iBin,iBin))->Clone() );
    for ( Int_t iBin = 1; iBin <= hTarget_syst->GetNbinsX(); iBin++ )   kSyst_Array.push_back( (TH1D*)(hTarget_syst->ProjectionX(Form("%s_%i",hTarget_syst->GetName(),iBin),iBin,iBin))->Clone() );
    for ( Int_t iBin = 1; iBin <= hTarget_stat->GetNbinsX(); iBin++ )   kStat_Array.at( iBin-1 )->Scale( 1./(kStat_Array.at(0)->GetBinWidth(iBin)) );
    for ( Int_t iBin = 1; iBin <= hTarget_syst->GetNbinsX(); iBin++ )   kSyst_Array.at( iBin-1 )->Scale( 1./(kSyst_Array.at(0)->GetBinWidth(iBin)) );
    //
    fStartTimer("2D Full Yield Evaluation");
    //
    for ( Int_t iTer = 0; iTer <  kStat_Array.size(); iTer++ )  {
        auto    kCurrent_Scale  =   0.5*(kStat_Array.at(0)->GetBinWidth(iTer+1))*(kStat_Array.at(0)->GetBinWidth(iTer+1));
        fResult.push_back( uMeasureFullYield(kStat_Array.at(iTer),kSyst_Array.at(iTer),fExtrapolationFunc_Array,kIterations,kSaveFolder,Form(kSaveName,iTer),kCurrent_Scale) );
        fPrintLoopTimer("2D Full Yield Evaluation",iTer+1,hTarget_stat->GetNbinsX()+2,1);
    }
    //
    auto    iBin    = 0;
    TH1D*   h2D_LowPT_Extrap_stat   =   (TH1D*)(hTarget_stat->ProjectionX(Form("%s_%i",hTarget_stat->GetName(),-1),1,1))->Clone();
    TH1D*   h2D_LowPT_Extrap_syst   =   (TH1D*)(hTarget_syst->ProjectionX(Form("%s_%i",hTarget_syst->GetName(),-1),1,1))->Clone();
    for ( auto kResult : fResult )    {
        iBin++;
        auto kBinValue  =   get<0>(kResult["YL_EXT"]);
        auto kBin_stat  =   get<1>(kResult["YL_EXT"]);
        auto kBin_syst  =   get<2>(kResult["YL_EXT"]);
        h2D_LowPT_Extrap_stat   ->  SetBinContent   ( iBin, kBinValue );
        h2D_LowPT_Extrap_stat   ->  SetBinError     ( iBin, kBin_stat );
        h2D_LowPT_Extrap_syst   ->  SetBinContent   ( iBin, kBinValue );
        h2D_LowPT_Extrap_syst   ->  SetBinError     ( iBin, kBin_syst );
    }
    //
    push_to_front( fResult, uMeasureFullYield(h2D_LowPT_Extrap_stat,h2D_LowPT_Extrap_syst,fExtrapolationFunc_Array,kIterations,kSaveFolder,Form(kSaveName,-1)) );
    //
    fStopTimer("2D Full Yield Evaluation");
    //
    add_ext_results(fResult);
    //
    return fResult;
}
//
template <  typename THXTarget_Type >
std::pair<std::pair<Float_t,Float_t>,std::pair<Float_t,Float_t>>
uExtrapolateToLowPT
( THXTarget_Type* hTarget_move, THXTarget_Type* hTarget_stat, std::tuple<TF1*,Float_t,Float_t,TString> fExtrapolationFunc, Int_t kIterations, TString kSaveFolder = "", TString kSaveName = "" ) {
    std::pair<std::pair<Float_t,Float_t>,std::pair<Float_t,Float_t>>    fNullResult;
    //
    //  --- Implementation for TH1 of 1D
    //  TODO: Get dimension as specxialisation variable and overload (?)
    if ( !uIsTHPairConsistent( hTarget_move, hTarget_stat ) )   return fNullResult;
    if ( uGetTHDimension( hTarget_move ) != 1 )                 return fNullResult;
    //
    //  --- Optimisation mode
    gROOT->SetBatch(kTRUE);
    //
    auto        hFullUncertainties      =   uSumErrors( hTarget_move, hTarget_stat );
    auto        hFullUncertainties_Draw =   ( THXTarget_Type* )( hFullUncertainties->Clone() );
                hFullUncertainties      ->  Fit( get<0>(fExtrapolationFunc), get<3>(fExtrapolationFunc), "", get<1>(fExtrapolationFunc), get<2>(fExtrapolationFunc) );
    auto        kStandardExtrapolation  =   get<0>(fExtrapolationFunc)->Integral( 0, hFullUncertainties->GetBinLowEdge(1) );
    auto        kStandardIntegral       =   hFullUncertainties->Integral("width");
    auto        kStandardYield          =   kStandardExtrapolation+kStandardIntegral;
    auto        kStandardMeanPT         =   uEvaluateMeanPT( hFullUncertainties, get<0>(fExtrapolationFunc) );
    //
    TCanvas*    cSaveMeasureFit         =   new TCanvas( "cSaveMeasureFit", "cSaveMeasureFit", 2000, 1200 );
    if ( !kSaveFolder.IsNull() ) {
        cSaveMeasureFit             ->  Divide(3,2);
        //
        hFullUncertainties_Draw     ->  SetTitle("");
        hFullUncertainties_Draw     ->  SetMarkerStyle  ( uGetMarker(0) );
        hFullUncertainties_Draw     ->  SetMarkerColor  ( uGetColor (2) );
        hFullUncertainties_Draw     ->  SetLineColor    ( uGetColor (2) );
        get<0>(fExtrapolationFunc)  ->  SetLineColor    ( uGetColor (1) );
        get<0>(fExtrapolationFunc)  ->  SetLineWidth    ( 2             );
        get<0>(fExtrapolationFunc)  ->  SetLineStyle    ( 1             );
        //
        auto    kHigLimit           =   hFullUncertainties  ->  GetXaxis()  ->  GetXmax();
        auto    kLowLimit           =   hFullUncertainties  ->  GetXaxis()  ->  GetXmin();
        auto    kMaximum            =   hFullUncertainties  ->  GetMaximum();
        auto    kMinimum            =   hFullUncertainties  ->  GetMinimum();
        auto    kBinFocusEdge       =   hFullUncertainties  ->  GetXaxis()  ->  FindBin( 0.25*( kLowLimit + kHigLimit ) );
        auto    kBinFocusContent    =   hFullUncertainties  ->  GetBinContent( kBinFocusEdge );
        hFullUncertainties_Draw     ->  SetMaximum( (50.)*kMaximum );
        hFullUncertainties_Draw     ->  SetMinimum( (0.2)*kMinimum );
        cSaveMeasureFit             ->  cd(1);
        gPad                        ->  SetLogy();
        hFullUncertainties_Draw     ->  GetXaxis()  ->  SetRangeUser( 0., kHigLimit*1.2 );
        hFullUncertainties_Draw     ->  DrawCopy();
        //
        cSaveMeasureFit             ->  cd(3);
        gPad                        ->  SetLogy();
        hFullUncertainties_Draw     ->  GetXaxis()  ->  SetRangeUser( 0.25*( kLowLimit + kHigLimit ), kHigLimit*1.2  );
        hFullUncertainties_Draw     ->  DrawCopy();
        //
        cSaveMeasureFit             ->  cd(2);
        gPad                        ->  SetLogy();
        hFullUncertainties_Draw     ->  GetXaxis()  ->  SetRangeUser( 0., 0.25*( kLowLimit + kHigLimit ) );
        hFullUncertainties_Draw     ->  SetMaximum( 10.*kMaximum );
        hFullUncertainties_Draw     ->  SetMinimum( 0.2*kBinFocusContent );
        hFullUncertainties_Draw     ->  DrawCopy();
        //
        cSaveMeasureFit             ->  cd(4);
        gPad                        ->  SetLogx();
        auto    hRatioReference_Draw=   ( THXTarget_Type* )( hFullUncertainties_Draw->Clone() );
        hRatioReference_Draw        ->  Divide( hFullUncertainties_Draw, hFullUncertainties_Draw );
        hRatioReference_Draw        ->  GetYaxis()  ->  SetTitle( "FIT/DATA" );
        hRatioReference_Draw        ->  SetMaximum( 1.5 );
        hRatioReference_Draw        ->  SetMinimum( 0.5 );
        hRatioReference_Draw        ->  SetFillColorAlpha( 0., 0. );
        hRatioReference_Draw        ->  DrawCopy("HIST L");
        //
        hFullUncertainties_Draw     ->  GetXaxis()  ->  SetRangeUser( 0., kHigLimit*1.2  );
    }
    //
    //  --- Evaluation Loop
    std::vector<Float_t>    kFullYieldVariations;
    std::vector<Float_t>    kFullMeanPTVariations;
    for ( Int_t iTer = 0; iTer < kIterations; ++iTer )  {
        TCanvas*    kDump   =   new TCanvas();
        auto    kRandomised =   uSumErrors( uRandomisePoints( hTarget_move ), hTarget_stat );
        kRandomised         ->  Fit( get<0>(fExtrapolationFunc), get<3>(fExtrapolationFunc), "", get<1>(fExtrapolationFunc), get<2>(fExtrapolationFunc) );
        delete      kDump;
        //
        kFullYieldVariations    .push_back  ( kStandardIntegral + get<0>(fExtrapolationFunc)->Integral( 0, get<1>(fExtrapolationFunc) ) );
        kFullMeanPTVariations   .push_back  ( uEvaluateMeanPT( hFullUncertainties, get<0>(fExtrapolationFunc) ) );
        //
        if ( !kSaveFolder.IsNull() ) {
            cSaveMeasureFit->cd(1);
            get<0>(fExtrapolationFunc)->DrawCopy("SAME");
            cSaveMeasureFit->cd(2);
            get<0>(fExtrapolationFunc)->DrawCopy("SAME");
            cSaveMeasureFit->cd(3);
            get<0>(fExtrapolationFunc)->DrawCopy("SAME");
            cSaveMeasureFit->cd(4);
            auto    hRatioReference_Draw=   ( THXTarget_Type* )( hFullUncertainties_Draw->Clone() );
            hRatioReference_Draw    ->  Divide( get<0>(fExtrapolationFunc) );
            hRatioReference_Draw    ->  SetLineColor( uGetColor (1) );
            hRatioReference_Draw    ->  SetLineWidth( 2 );
            hRatioReference_Draw    ->  SetFillColorAlpha( 0., 0. );
            hRatioReference_Draw    ->  DrawCopy("SAME HIST L");
        }
    }
    //
    //  --- Building Error Evaluation Histograms
    auto    hFullYieldUncertEval    =   uBuildTH1(  kFullYieldVariations,   max( 40., kIterations/5. ) );
    auto    hFullMeanPTUncertEval   =   uBuildTH1(  kFullMeanPTVariations,  max( 40., kIterations/5. ) );
    //
    auto    kYieldMeanContrib       =   fabs( hFullYieldUncertEval->GetMean() - kStandardYield );
    auto    kYieldStdvContrib       =   hFullYieldUncertEval->GetRMS();
    auto    kMnPT_MeanContrib       =   fabs( hFullMeanPTUncertEval->GetMean() - kStandardMeanPT );
    auto    kMnPT_StdvContrib       =   hFullMeanPTUncertEval->GetRMS();
    //
    if ( !kSaveFolder.IsNull() ) {
        cSaveMeasureFit         ->  cd(1);
        hFullUncertainties_Draw ->  Draw("SAME");
        uLatex                  ->  DrawLatexNDC( 0.61, 0.81, Form("#frac{dN}{dy}: %.6f",kStandardExtrapolation+kStandardIntegral) );
        uLatex                  ->  DrawLatexNDC( 0.60, 0.71, Form("#LT#it{p}_{T}#GT: %.6f",kStandardMeanPT) );
        //
        cSaveMeasureFit         ->  cd(2);
        hFullUncertainties_Draw ->  Draw("SAME");
        uLatex                  ->  DrawLatexNDC( 0.61, 0.81, Form("#frac{dN}{dy}: %.6f",kStandardExtrapolation+kStandardIntegral) );
        uLatex                  ->  DrawLatexNDC( 0.60, 0.71, Form("#LT#it{p}_{T}#GT: %.6f",kStandardMeanPT) );
        //
        cSaveMeasureFit         ->  cd(3);
        hFullUncertainties_Draw ->  Draw("SAME");
        uLatex                  ->  DrawLatexNDC( 0.61, 0.81, Form("#frac{dN}{dy}: %.6f",kStandardExtrapolation+kStandardIntegral) );
        uLatex                  ->  DrawLatexNDC( 0.60, 0.71, Form("#LT#it{p}_{T}#GT: %.6f",kStandardMeanPT) );
        //
        cSaveMeasureFit             ->  cd(4);
        auto    hRatioReference_Draw=   ( THXTarget_Type* )( hFullUncertainties_Draw->Clone() );
        hRatioReference_Draw        ->  Divide( hFullUncertainties_Draw, hFullUncertainties_Draw );
        hRatioReference_Draw        ->  SetFillColorAlpha( 0., 0. );
        hRatioReference_Draw        ->  DrawCopy("SAME HIST L");
        for ( Int_t iBin = 1; iBin <= hTarget_move->GetNbinsX(); iBin++ ) hRatioReference_Draw->SetBinError( iBin, hTarget_move->GetBinError(iBin) / hTarget_move->GetBinContent(iBin)  );
        hRatioReference_Draw        ->  SetFillColorAlpha( kGray, 0.5 );
        hRatioReference_Draw        ->  DrawCopy("SAME E3");
        //
        cSaveMeasureFit->cd(5);
        hFullYieldUncertEval    ->  SetLineColor    ( uGetColor (2) );
        hFullYieldUncertEval    ->  SetMaximum( hFullYieldUncertEval->GetMaximum()*1.4 );
        hFullYieldUncertEval    ->  GetXaxis()  ->  SetNdivisions(5);
        hFullYieldUncertEval    ->  GetXaxis()  ->  SetTitle("#frac{dN}{dy}");
        hFullYieldUncertEval    ->  Draw();
        uLatex                  ->  DrawLatexNDC( 0.180, 0.82, Form("#sigma_{mean}: %2.2f%s",100*(kYieldMeanContrib)/(kStandardYield), "%") );
        uLatex                  ->  DrawLatexNDC( 0.195, 0.76, Form("#sigma_{stdv}: %2.2f%s",100*(kYieldStdvContrib)/(kStandardYield), "%") );
        uLatex                  ->  DrawLatexNDC( 0.650, 0.80, Form("#sigma_{totl}: %2.2f%s",100*(kYieldMeanContrib+kYieldStdvContrib)/(kStandardYield), "%") );
        //
        cSaveMeasureFit->cd(6);
        hFullMeanPTUncertEval   ->  SetLineColor    ( uGetColor (2) );
        hFullMeanPTUncertEval   ->  SetMaximum( hFullMeanPTUncertEval->GetMaximum()*1.4 );
        hFullMeanPTUncertEval   ->  GetXaxis()  ->  SetNdivisions(5);
        hFullMeanPTUncertEval   ->  GetXaxis()  ->  SetTitle("#LT#it{p}_{T}#GT (GeV/#it{c})");
        hFullMeanPTUncertEval   ->  Draw();
        uLatex                  ->  DrawLatexNDC( 0.180, 0.82, Form("#sigma_{mean}: %2.2f%s",100*(kMnPT_MeanContrib)/(kStandardMeanPT), "%") );
        uLatex                  ->  DrawLatexNDC( 0.195, 0.76, Form("#sigma_{stdv}: %2.2f%s",100*(kMnPT_StdvContrib)/(kStandardMeanPT), "%") );
        uLatex                  ->  DrawLatexNDC( 0.650, 0.80, Form("#sigma_{totl}: %2.2f%s",100*(kMnPT_MeanContrib+kMnPT_StdvContrib)/(kStandardMeanPT), "%") );
        //
        cSaveMeasureFit         ->  SaveAs( ( kSaveFolder + TString("/") + kSaveName + TString(".pdf") ) );
        cSaveMeasureFit         ->  SaveAs( ( kSaveFolder + TString("/") + kSaveName + TString(".eps") ) );
    }
    delete cSaveMeasureFit;
    //
    //  --- Optimisation off
    gROOT->SetBatch(kFALSE);
    //
    return  std::pair<std::pair<Float_t,Float_t>,std::pair<Float_t,Float_t>> {{kStandardYield,(kYieldMeanContrib+kYieldStdvContrib)},{kStandardMeanPT,(kMnPT_MeanContrib+kMnPT_StdvContrib)}};;
}
//
template <  typename THXTarget_Type >
std::pair<std::pair<Float_t,Float_t>,std::pair<Float_t,Float_t>>
uExtrapolateToLowPT
( THXTarget_Type* hTarget_move, THXTarget_Type* hTarget_stat, std::vector<std::tuple<TF1*,Float_t,Float_t,TString>> fExtrapolationFunc = kStandardSystematicFitFunctions, TString kSaveFolder = "", TString kSaveName = "" ) {
    std::pair<std::pair<Float_t,Float_t>,std::pair<Float_t,Float_t>>    fNullResult;
    //
    //  --- Implementation for TH1 of 1D
    //  TODO: Get dimension as specxialisation variable and overload (?)
    if ( !uIsTHPairConsistent( hTarget_move, hTarget_stat ) )   return fNullResult;
    if ( uGetTHDimension( hTarget_move ) != 1 )                 return fNullResult;
    //
    //  --- Optimisation mode
    gROOT->SetBatch(kTRUE);
    //
    auto        hFullUncertainties      =   uSumErrors( hTarget_move, hTarget_stat );
    auto        hFullUncertainties_Draw =   ( THXTarget_Type* )( hFullUncertainties->Clone() );
                hFullUncertainties      ->  Fit( get<0>(fExtrapolationFunc.at(0)), get<3>(fExtrapolationFunc.at(0)), "", get<1>(fExtrapolationFunc.at(0)), get<2>(fExtrapolationFunc.at(0)) );
    auto        kStandardExtrapolation  =   get<0>(fExtrapolationFunc.at(0))->Integral( 0, get<1>(fExtrapolationFunc.at(0)) );
    auto        kStandardIntegral       =   hFullUncertainties->Integral("width");
    auto        kStandardYield          =   kStandardExtrapolation+kStandardIntegral;
    auto        kStandardMeanPT         =   uEvaluateMeanPT( hFullUncertainties, get<0>(fExtrapolationFunc.at(0)) );
    //
    TCanvas*    cSaveMeasureFit         =   new TCanvas( "cSaveMeasureFit", "cSaveMeasureFit", 2000, 1200 );
    TLegend*    cAllFunctions           =   new TLegend( 0.18, 0.87, 0.6, 0.7 );
    if ( !kSaveFolder.IsNull() ) {
        cSaveMeasureFit             ->  Divide(3,2);
        //
        hFullUncertainties_Draw     ->  SetTitle("");
        hFullUncertainties_Draw     ->  SetMarkerStyle  ( uGetMarker(0) );
        hFullUncertainties_Draw     ->  SetMarkerColor  ( kBlack );
        hFullUncertainties_Draw     ->  SetLineColor    ( kBlack );
        hFullUncertainties_Draw     ->  SetFillColor    ( kBlack );
        //
        auto    kHigLimit           =   hFullUncertainties  ->  GetXaxis()  ->  GetXmax();
        auto    kLowLimit           =   hFullUncertainties  ->  GetXaxis()  ->  GetXmin();
        auto    kMaximum            =   hFullUncertainties  ->  GetMaximum();
        auto    kMinimum            =   hFullUncertainties  ->  GetMinimum();
        auto    kBinFocusEdge       =   hFullUncertainties  ->  GetXaxis()  ->  FindBin( 0.25*( kLowLimit + kHigLimit ) );
        auto    kBinFocusContent    =   hFullUncertainties  ->  GetBinContent( kBinFocusEdge );
        hFullUncertainties_Draw     ->  SetMaximum( (50.)*kMaximum );
        hFullUncertainties_Draw     ->  SetMinimum( (0.2)*kMinimum );
        cSaveMeasureFit             ->  cd(1);
        gPad                        ->  SetLogy();
        hFullUncertainties_Draw     ->  GetXaxis()  ->  SetRangeUser( 0., kHigLimit*1.2 );
        hFullUncertainties_Draw     ->  DrawCopy();
        //
        cSaveMeasureFit             ->  cd(3);
        gPad                        ->  SetLogy();
        hFullUncertainties_Draw     ->  GetXaxis()  ->  SetRangeUser( 0.25*( kLowLimit + kHigLimit ), kHigLimit*1.2  );
        hFullUncertainties_Draw     ->  DrawCopy();
        //
        cSaveMeasureFit             ->  cd(2);
        gPad                        ->  SetLogy();
        hFullUncertainties_Draw     ->  GetXaxis()  ->  SetRangeUser( 0., 0.25*( kLowLimit + kHigLimit ) );
        hFullUncertainties_Draw     ->  SetMaximum( 10.*kMaximum );
        hFullUncertainties_Draw     ->  SetMinimum( 0.2*kBinFocusContent );
        hFullUncertainties_Draw     ->  DrawCopy();
        //
        cSaveMeasureFit             ->  cd(4);
        gPad                        ->  SetLogx();
        auto    hRatioReference_Draw=   ( THXTarget_Type* )( hFullUncertainties_Draw->Clone() );
        hRatioReference_Draw        ->  Divide( hFullUncertainties_Draw, hFullUncertainties_Draw );
        hRatioReference_Draw        ->  GetYaxis()  ->  SetTitle( "FIT/DATA" );
        hRatioReference_Draw        ->  SetMaximum( 1.5 );
        hRatioReference_Draw        ->  SetMinimum( 0.5 );
        hRatioReference_Draw        ->  SetFillColorAlpha( 0., 0. );
        hRatioReference_Draw        ->  DrawCopy("HIST L");
        //
        hFullUncertainties_Draw     ->  GetXaxis()  ->  SetRangeUser( 0., kHigLimit*1.2  );
    }
    //
    //  --- Evaluation Loop
    std::vector<Float_t>    kFullYieldVariations;
    std::vector<Float_t>    kFullMeanPTVariations;
    std::map<TString,Int_t> kDrawColor;
    Int_t                   kDrawColorIterator  = 1;
    for ( auto kCurrent_Function : fExtrapolationFunc )  {
        TCanvas*    kDump   =   new TCanvas();
        hFullUncertainties  ->  Fit( get<0>(kCurrent_Function), get<3>(kCurrent_Function), "", get<1>(kCurrent_Function), get<2>(kCurrent_Function) );
        delete      kDump;
        //
        kFullYieldVariations    .push_back  ( kStandardIntegral + get<0>(kCurrent_Function)->Integral( 0, hFullUncertainties->GetBinLowEdge(1) ) );
        kFullMeanPTVariations   .push_back  ( uEvaluateMeanPT( hFullUncertainties, get<0>(kCurrent_Function) ) );
        //
        if ( !kSaveFolder.IsNull() ) {
            if ( kDrawColor[get<0>(kCurrent_Function)->GetName()] == 0 )    {
                kDrawColor[get<0>(kCurrent_Function)->GetName()] = kDrawColorIterator;
                get<0>(kCurrent_Function)->SetLineColor( uGetColor( kDrawColorIterator ) );
                get<0>(kCurrent_Function)->SetLineStyle( kDrawColorIterator%2 + 3 );
                get<0>(kCurrent_Function)->SetLineWidth( 2 );
                cAllFunctions->AddEntry( get<0>(kCurrent_Function), get<0>(kCurrent_Function)->GetName() );
                kDrawColorIterator++;
            }
            cSaveMeasureFit->cd(1);
            get<0>(kCurrent_Function)-> DrawCopy("SAME");
            cSaveMeasureFit->cd(2);
            get<0>(kCurrent_Function)-> DrawCopy("SAME");
            cSaveMeasureFit->cd(3);
            get<0>(kCurrent_Function)-> DrawCopy("SAME");
            cSaveMeasureFit->cd(4);
            auto    hRatioReference_Draw=   ( THXTarget_Type* )( hFullUncertainties_Draw->Clone() );
            hRatioReference_Draw    ->  Divide( get<0>(kCurrent_Function) );
            hRatioReference_Draw    ->  SetLineColor( uGetColor( kDrawColor[get<0>(kCurrent_Function)->GetName()] ) );
            hRatioReference_Draw    ->  SetLineStyle( kDrawColor[get<0>(kCurrent_Function)->GetName()]%2 + 3 );
            hRatioReference_Draw    ->  SetFillColorAlpha( 0., 0. );
            hRatioReference_Draw    ->  DrawCopy("SAME HIST L");
        }
    }
    //
    //  --- Building Error Evaluation Histograms
    auto    hFullYieldUncertEval    =   uBuildTH1(  kFullYieldVariations,   max( 40., fExtrapolationFunc.size()/5. ) );
    auto    hFullMeanPTUncertEval   =   uBuildTH1(  kFullMeanPTVariations,  max( 40., fExtrapolationFunc.size()/5. ) );
    //
    auto    kYieldMeanContrib       =   fabs( hFullYieldUncertEval->GetMean() - kStandardYield );
    auto    kYieldStdvContrib       =   hFullYieldUncertEval->GetRMS();
    auto    kMnPT_MeanContrib       =   fabs( hFullMeanPTUncertEval->GetMean() - kStandardMeanPT );
    auto    kMnPT_StdvContrib       =   hFullMeanPTUncertEval->GetRMS();
    //
    if ( !kSaveFolder.IsNull() ) {
        cSaveMeasureFit         ->  cd(1);
        hFullUncertainties_Draw ->  Draw("SAME");
        cAllFunctions           ->  SetNColumns(2);
        cAllFunctions           ->  Draw("SAME");
        uLatex                  ->  DrawLatexNDC( 0.61, 0.81, Form("#frac{dN}{dy}: %.6f",kStandardExtrapolation+kStandardIntegral) );
        uLatex                  ->  DrawLatexNDC( 0.60, 0.71, Form("#LT#it{p}_{T}#GT: %.6f",kStandardMeanPT) );
        //
        cSaveMeasureFit         ->  cd(2);
        hFullUncertainties_Draw ->  Draw("SAME");
        cAllFunctions           ->  SetNColumns(2);
        cAllFunctions           ->  Draw("SAME");
        uLatex                  ->  DrawLatexNDC( 0.61, 0.81, Form("#frac{dN}{dy}: %.6f",kStandardExtrapolation+kStandardIntegral) );
        uLatex                  ->  DrawLatexNDC( 0.60, 0.71, Form("#LT#it{p}_{T}#GT: %.6f",kStandardMeanPT) );
        //
        cSaveMeasureFit         ->  cd(3);
        hFullUncertainties_Draw ->  Draw("SAME");
        cAllFunctions           ->  SetNColumns(2);
        cAllFunctions           ->  Draw("SAME");
        uLatex                  ->  DrawLatexNDC( 0.61, 0.81, Form("#frac{dN}{dy}: %.6f",kStandardExtrapolation+kStandardIntegral) );
        uLatex                  ->  DrawLatexNDC( 0.60, 0.71, Form("#LT#it{p}_{T}#GT: %.6f",kStandardMeanPT) );
        //
        cSaveMeasureFit             ->  cd(4);
        auto    hRatioReference_Draw=   ( THXTarget_Type* )( hFullUncertainties_Draw->Clone() );
        hRatioReference_Draw        ->  Divide( hFullUncertainties_Draw, hFullUncertainties_Draw );
        hRatioReference_Draw        ->  SetFillColorAlpha( 0., 0. );
        hRatioReference_Draw        ->  DrawCopy("SAME HIST L");
        for ( Int_t iBin = 1; iBin <= hTarget_move->GetNbinsX(); iBin++ ) hRatioReference_Draw->SetBinError( iBin, hTarget_move->GetBinError(iBin) / hTarget_move->GetBinContent(iBin) );
        hRatioReference_Draw        ->  SetFillColorAlpha( kGray, 0.5 );
        hRatioReference_Draw        ->  DrawCopy("SAME E3");
        //
        cSaveMeasureFit         ->  cd(5);
        hFullYieldUncertEval    ->  SetLineColor    ( uGetColor (2) );
        hFullYieldUncertEval    ->  SetMaximum( hFullYieldUncertEval->GetMaximum()*1.4 );
        hFullYieldUncertEval    ->  GetXaxis()  ->  SetNdivisions(5);
        hFullYieldUncertEval    ->  GetXaxis()  ->  SetTitle("#frac{dN}{dy}");
        hFullYieldUncertEval    ->  Draw();
        uLatex                  ->  DrawLatexNDC( 0.180, 0.82, Form("#sigma_{mean}: %2.2f%s",100*(kYieldMeanContrib)/(kStandardYield), "%") );
        uLatex                  ->  DrawLatexNDC( 0.195, 0.76, Form("#sigma_{stdv}: %2.2f%s",100*(kYieldStdvContrib)/(kStandardYield), "%") );
        uLatex                  ->  DrawLatexNDC( 0.650, 0.80, Form("#sigma_{totl}: %2.2f%s",100*(kYieldMeanContrib+kYieldStdvContrib)/(kStandardYield), "%") );
        //
        cSaveMeasureFit         ->  cd(6);
        hFullMeanPTUncertEval   ->  SetLineColor    ( uGetColor (2) );
        hFullMeanPTUncertEval   ->  SetMaximum( hFullMeanPTUncertEval->GetMaximum()*1.4 );
        hFullMeanPTUncertEval   ->  GetXaxis()  ->  SetNdivisions(5);
        hFullMeanPTUncertEval   ->  GetXaxis()  ->  SetTitle("#LT#it{p}_{T}#GT (GeV/#it{c})");
        hFullMeanPTUncertEval   ->  Draw();
        uLatex                  ->  DrawLatexNDC( 0.180, 0.82, Form("#sigma_{mean}: %2.2f%s",100*(kMnPT_MeanContrib)/(kStandardMeanPT), "%") );
        uLatex                  ->  DrawLatexNDC( 0.195, 0.76, Form("#sigma_{stdv}: %2.2f%s",100*(kMnPT_StdvContrib)/(kStandardMeanPT), "%") );
        uLatex                  ->  DrawLatexNDC( 0.650, 0.80, Form("#sigma_{totl}: %2.2f%s",100*(kMnPT_MeanContrib+kMnPT_StdvContrib)/(kStandardMeanPT), "%") );
        //
        cSaveMeasureFit         ->  SaveAs( ( kSaveFolder + TString("/") + kSaveName + TString(".pdf") ) );
        cSaveMeasureFit         ->  SaveAs( ( kSaveFolder + TString("/") + kSaveName + TString(".eps") ) );
    }
    delete cSaveMeasureFit;
    //
    //  --- Optimisation off
    gROOT->SetBatch(kFALSE);
    //
    return  std::pair<std::pair<Float_t,Float_t>,std::pair<Float_t,Float_t>> {{kStandardYield,(kYieldMeanContrib+kYieldStdvContrib)},{kStandardMeanPT,(kMnPT_MeanContrib+kMnPT_StdvContrib)}};
}
//
template <  typename THXTarget_Type >
Float_t
uEvaluateMeanPT
( THXTarget_Type* hTarget, TF1* fLowPtModel )  {
    Float_t fResult         = 0;
    Float_t fDenominator    = 0;
    //
    //  --- Utility vectors
    std::vector<float>  kBinContent;
    std::vector<float>  kBinWitdh;
    std::vector<float>  kBinMeanPt;
    //
    //  Extrapolated Contribution
    if ( fLowPtModel )     {
        kBinContent.    push_back( fLowPtModel  ->  Integral(   0, hTarget->GetBinLowEdge(1) ) );
        kBinWitdh.      push_back( 1. );
        kBinMeanPt.     push_back( fLowPtModel  ->  Moment( 1,  0, hTarget->GetBinLowEdge(1) ) );
    }
    //
    //  Integrated Contribution
    for ( Int_t iBin = 1; iBin <= hTarget->GetNbinsX(); iBin++ )   {
        kBinContent.    push_back( hTarget->GetBinContent   (iBin) );
        kBinWitdh.      push_back( hTarget->GetBinWidth     (iBin) );
        kBinMeanPt.     push_back( hTarget->GetBinCenter    (iBin) );
    }
    //
    //  Calculate Mean PT
    for ( Int_t iBin = 0; iBin < kBinContent.size(); iBin++ )   {
        fResult         +=  kBinContent.at(iBin) * kBinWitdh.at(iBin) * kBinMeanPt.at(iBin);
        fDenominator    +=  kBinContent.at(iBin) * kBinWitdh.at(iBin);
    }
    return fResult/fDenominator;
}
//
template <  typename THXTarget_Type >
Double_t
uEvaluateMeanPTError
( THXTarget_Type* hTarget  )  {
    Float_t fHighLimit      = 0;
    Float_t fLowLimit       = 0;
    Float_t fHighDenom      = 0;
    Float_t fLowDenom       = 0;
    //
    //  --- Utility vectors
    std::vector<float>  kBinContent_High;
    std::vector<float>  kBinContent_Low;
    std::vector<float>  kBinWitdh;
    std::vector<float>  kBinMeanPt;
    //
    //  Integrated Contribution
    for ( Int_t iBin = 1; iBin <= hTarget->GetNbinsX(); iBin++ )   {
        kBinContent_High.   push_back( hTarget->GetBinContent   (iBin) + hTarget->GetBinError   (iBin) );
        kBinContent_Low.    push_back( hTarget->GetBinContent   (iBin) - hTarget->GetBinError   (iBin) );
        kBinWitdh.          push_back( hTarget->GetBinWidth     (iBin) );
        kBinMeanPt.         push_back( hTarget->GetBinCenter    (iBin) );
    }
    //
    //  Calculate Mean PT
    for ( Int_t iBin = 0; iBin < kBinWitdh.size(); iBin++ )   {
        fHighLimit      +=  kBinContent_High.at(iBin)   * kBinWitdh.at(iBin)    * kBinMeanPt.at(iBin);
        fLowLimit       +=  kBinContent_Low.at(iBin)    * kBinWitdh.at(iBin)    * kBinMeanPt.at(iBin);
        fHighDenom      +=  kBinContent_High.at(iBin)   * kBinWitdh.at(iBin);
        fLowDenom       +=  kBinContent_Low.at(iBin)    * kBinWitdh.at(iBin);
    }
    return 0.5*fabs( ( fLowLimit / fLowDenom ) - ( fHighLimit / fHighDenom ) );
}
//
#endif /* AAU_Extrapolation_h */
