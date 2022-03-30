//
//  Part of the AliAnalysisUtility package
//
//  Utilities for Systematics calculation and handling
//
//  Author              Nicola Rubini
//  Mail                nicola.rubini@cern.ch
//  Created             01/02/2022
//  Last modified       01/02/2022
#ifndef AAU_Systematics_h
#define AAU_Systematics_h
//
//  Global File w/ Constants and other functions
#include "AliAnalysisUtility.h"
//

//_____________________________________________________________________________
TH1F*
uProjectBarlow
( TH1F *hTarget, TH1F  *hTester )   {
    std::vector<Float_t> uAllPoints;
    for ( Int_t iBin = 0; iBin < hTarget->GetNbinsX(); iBin++ ) {
        auto fYValue    =   hTarget->GetBinContent  (iBin+1);
        auto fYError    =   hTarget->GetBinError    (iBin+1);
        auto fTValue    =   hTester->GetBinContent  (iBin+1);
        auto fTError    =   hTester->GetBinError    (iBin+1);
        uAllPoints.push_back(uBarlowPar(fYValue,fYError,fTValue,fTError));
    }
    return uBuildTH1F(uAllPoints,600,0.,-30.,30.);
}
TH1F*
uProjectBarlow
( TH2F *hTarget, TH2F  *hTester, Bool_t fDiag = false )   {
    std::vector<Float_t> uAllPoints;
    for ( Int_t iBin = 0; iBin < hTarget->GetNbinsX(); iBin++ ) {
        for ( Int_t jBin = 0; jBin < hTarget->GetNbinsY(); jBin++ ) {
            if ( ( iBin > jBin ) && fDiag ) continue;
            auto fYValue    =   hTarget->GetBinContent  (iBin+1,jBin+1);
            auto fYError    =   hTarget->GetBinError    (iBin+1,jBin+1);
            auto fTValue    =   hTester->GetBinContent  (iBin+1,jBin+1);
            auto fTError    =   hTester->GetBinError    (iBin+1,jBin+1);
            uAllPoints.push_back(uBarlowPar(fYValue,fYError,fTValue,fTError));
        }
    }
    return uBuildTH1F(uAllPoints,600,0.,-30.,30.);
}

//_____________________________________________________________________________
Bool_t
uIsRelevantVariation
( TH1F* hCheckHist, TString fFolder = "", TString fName = ""  )   {
    //  --- General Settings
    auto kMaxMean       =   0.80;    // Maximum absolute value of Barlow parameter mean
    auto kMaxStdv       =   1.30;    // Maximum absolute value of Barlow parameter stdv
    auto kMinInt1       =   0.55;    // Maximum absolute value of Barlow parameter Integral within +-1stdv
    auto kMinInt2       =   0.75;    // Maximum absolute value of Barlow parameter Integral within +-2stdv
    //
    auto nPassedChecks  =   0;
    auto uMean          =   hCheckHist->GetMean();
    auto uSTDV          =   hCheckHist->GetRMS();
    auto uIntegralTot   =   hCheckHist->GetEntries();
    auto uIntegral1Sg   =   0.;
    auto uIntegral2Sg   =   0.;
    for ( Int_t iBin = 0; iBin <= hCheckHist->GetNbinsX(); iBin++ )    {
        auto fXValue    =   hCheckHist->GetBinCenter(iBin+1);
        if ( fabs( fXValue - uMean ) <= uSTDV   ) uIntegral1Sg += hCheckHist->GetBinContent(iBin+1);
        if ( fabs( fXValue - uMean ) <= 2*uSTDV ) uIntegral2Sg += hCheckHist->GetBinContent(iBin+1);
    }
    if ( fabs(uMean)                <= kMaxMean ) nPassedChecks++;
    if ( uSTDV                      <= kMaxStdv ) nPassedChecks++;
    if ( uIntegral1Sg/uIntegralTot  >= kMinInt1 ) nPassedChecks++;
    if ( uIntegral2Sg/uIntegralTot  >= kMinInt2 ) nPassedChecks++;
    if ( !fFolder.IsNull() )   {
        gROOT->SetBatch();
        TCanvas    *cDrawResult =   new TCanvas();
        gStyle->SetOptStat(0);
        TLatex     *fLatex      =   new TLatex();
        hCheckHist  ->  SetTitle("");
        hCheckHist  ->  SetLineColor(kColors[2]);
        auto fMaxium = hCheckHist  ->  GetMaximum();
        hCheckHist  ->  SetMaximum(fMaxium*1.3);
        hCheckHist  ->  GetXaxis()->SetTitle("Barlow parameter");
        hCheckHist  ->  Draw("HIST");
        fLatex->SetTextFont(60);
        fLatex->SetTextSize(0.05);
        fLatex      ->  DrawLatexNDC(0.80,0.75,Form("%s",hCheckHist->GetName()));
        fLatex->SetTextFont(42);
        fLatex->SetTextSize(0.04);
        if ( fabs(uMean)                <= kMaxMean )   fLatex  ->  DrawLatexNDC(0.18,0.84,Form("#color[8]{#mu:   %.2f < %.2f}",fabs(uMean),kMaxMean));
        else                                            fLatex  ->  DrawLatexNDC(0.18,0.84,Form("#color[2]{#mu:   %.2f < %.2f}",fabs(uMean),kMaxMean));
        if ( uSTDV                      <= kMaxStdv )   fLatex  ->  DrawLatexNDC(0.18,0.79,Form("#color[8]{#sigma:   %.2f < %.2f}",uSTDV,kMaxStdv));
        else                                            fLatex  ->  DrawLatexNDC(0.18,0.79,Form("#color[2]{#sigma:   %.2f < %.2f}",uSTDV,kMaxStdv));
        if ( uIntegral1Sg/uIntegralTot  >= kMinInt1 )   fLatex  ->  DrawLatexNDC(0.18,0.74,Form("#color[8]{INT_{1#sigma}: %.2f > %.2f}",uIntegral1Sg/uIntegralTot,kMinInt1));
        else                                            fLatex  ->  DrawLatexNDC(0.18,0.74,Form("#color[2]{INT_{1#sigma}: %.2f > %.2f}",uIntegral1Sg/uIntegralTot,kMinInt1));
        if ( uIntegral2Sg/uIntegralTot  >= kMinInt2 )   fLatex  ->  DrawLatexNDC(0.18,0.69,Form("#color[8]{INT_{2#sigma}: %.2f > %.2f}",uIntegral2Sg/uIntegralTot,kMinInt2));
        else                                            fLatex  ->  DrawLatexNDC(0.18,0.69,Form("#color[2]{INT_{2#sigma}: %.2f > %.2f}",uIntegral2Sg/uIntegralTot,kMinInt2));
        fLatex->SetTextSize(0.075);
        if ( nPassedChecks < 3 )    fLatex      ->  DrawLatexNDC(0.55,0.82,"#color[2]{CONSIDERED}");
        else                        fLatex      ->  DrawLatexNDC(0.60,0.82,"#color[8]{REJECTED}");
        gROOT       ->  ProcessLine (Form(".! mkdir -p %s",fFolder.Data()));
        cDrawResult ->  SaveAs((fFolder+fName+TString(Form("_%s.pdf",hCheckHist->GetName()))).Data());
        gROOT->SetBatch(kFALSE);
    }
    return nPassedChecks < 3;
}
std::vector<Bool_t>
uIsRelevantVariation
( TH1F* hStandard, std::vector<TH1F*> hVariations, TString fFolder = "", TString fName = "" )   {
    std::vector<Bool_t> uResults;
    for ( auto hVariation : hVariations ) {
        auto hTester    =   (TH1F*)hStandard->Clone();
        auto hTested    =   (TH1F*)hVariation->Clone();
        auto hCheckHist =   uProjectBarlow(hTester,hTested);
        hCheckHist      ->  SetName(hVariation->GetName());
        uResults.push_back( uIsRelevantVariation(hCheckHist,fFolder,fName) );
    }
    return uResults;
}
std::vector<Bool_t>
uIsRelevantVariation
( TH2F* hStandard, std::vector<TH2F*> hVariations, TString fFolder = "", TString fName = "", Bool_t fDiag = false  )   {
    std::vector<Bool_t> uResults;
    for ( auto hVariation : hVariations ) {
        auto hTester    =   (TH2F*)hStandard->Clone();
        auto hTested    =   (TH2F*)hVariation->Clone();
        auto hCheckHist =   uProjectBarlow(hTester,hTested,fDiag);
        hCheckHist      ->  SetName(hVariation->GetName());
        uResults.push_back( uIsRelevantVariation(hCheckHist,fFolder,fName) );
    }
    return uResults;
}

//_____________________________________________________________________________
std::vector<TH1F*>
uBuildVariationSpectra
( TH1F* hStandard, std::vector<TH1F*> hVariations, std::vector<Bool_t> uWhichIsToBeconsidered ) {
    std::vector<TH1F*> uResults;
    auto iTer = -1;
    for ( auto hVariation : hVariations ) {
        iTer++;
        if ( uWhichIsToBeconsidered.size() > iTer && !uWhichIsToBeconsidered.at(iTer) ) continue;
        auto hTester    =   (TH1F*)hStandard->Clone();
        hTester->Divide(hVariation,hStandard);
        uResults.push_back(hTester);
    }
    return uResults;
}
std::vector<TH1F*>
uBuildVariationSpectra
( TH1F* hStandard, std::vector<TH1F*> hVariations ) {
    std::vector<Bool_t> fDump;
    for ( auto fTested : hVariations ) fDump.push_back(kTRUE);
    return uBuildVariationSpectra(hStandard,hVariations,fDump);
}
std::vector<TH2F*>
uBuildVariationSpectra
( TH2F* hStandard, std::vector<TH2F*> hVariations, std::vector<Bool_t> uWhichIsToBeconsidered ) {
    std::vector<TH2F*> uResults;
    auto iTer = -1;
    for ( auto hVariation : hVariations ) {
        iTer++;
        if ( uWhichIsToBeconsidered.size() != 0 && !uWhichIsToBeconsidered.at(iTer) ) continue;
        auto hTester    =   (TH2F*)hStandard->Clone();
        hTester->Divide(hVariation,hStandard);
        uResults.push_back(hTester);
    }
    return uResults;
}
std::vector<TH2F*>
uBuildVariationSpectra
( TH2F* hStandard, std::vector<TH2F*> hVariations ) {
    std::vector<Bool_t> fDump;
    for ( auto fTested : hVariations ) fDump.push_back(kTRUE);
    return uBuildVariationSpectra(hStandard,hVariations,fDump);
}

//_____________________________________________________________________________
std::vector<TH1F*>
uBuildVariationBinByBin // Put bool to check or not for barlow, put bool to have point by point barlow check
( TH1F* hStandard, std::vector<TH1F*> hVariations, TString fFolder = "", TString fName = "", std::vector<Bool_t> fWhichToChoose = std::vector<Bool_t> () ) {
    std::vector<TH1F*>  hResults;
    auto    uSystVariationSpectra_  =   uBuildVariationSpectra(hStandard,hVariations,fWhichToChoose);
    for ( Int_t iPT = 0; iPT < hStandard->GetNbinsX(); iPT++ )    {
        hResults.push_back( new TH1F(Form("VarHist_%i",iPT),Form("VarHist_%i",iPT),2000,-1,+1) );
        for ( auto  hVariationSpectrum : uSystVariationSpectra_ ) {
            hResults.at(iPT)->Fill(hVariationSpectrum->GetBinContent(iPT+1)-1);
        }
        if ( !fFolder.IsNull() )    {
            gROOT->SetBatch();
            TCanvas    *cDrawResult =   new TCanvas();
            hResults.at(iPT)->Draw();
            cDrawResult->SaveAs((fFolder+fName+TString(Form("_%i.pdf",iPT))).Data());
            gROOT->SetBatch(kFALSE);
        }
    }
    return hResults;
}
std::vector<std::pair<TH1F*,Int_t>>
uBuildVariationBinByBin // Put bool to check or not for barlow, put bool to have point by point barlow check
( TH2F* hStandard, std::vector<TH2F*> hVariations, TString fFolder = "", TString fName = "", std::vector<Bool_t> fWhichToChoose = std::vector<Bool_t> (), Bool_t fDiag = false ) {
    std::vector<std::pair<TH1F*,Int_t>>  hResults;
    auto    uSystVariationSpectra_  =   uBuildVariationSpectra( hStandard, hVariations, fWhichToChoose );
    auto    iTH1    =   0;
    for ( Int_t iPT = 0; iPT < hStandard->GetNbinsX(); iPT++ )    {
        for ( Int_t jPT = 0; jPT < hStandard->GetNbinsY(); jPT++ )    {
            if ( ( iPT > jPT ) && fDiag ) continue;
            hResults.push_back( std::pair<TH1F*,Int_t>( new TH1F(Form("VarHist_%i_%i",iPT,jPT),Form("VarHist_%i_%i",iPT,jPT),2000,-1,+1) , hStandard->GetBin( iPT+1, jPT+1 ) ) );
            for ( auto  hVariationSpectrum : uSystVariationSpectra_ ) {
                hResults.at(iTH1).first->Fill(hVariationSpectrum->GetBinContent(iPT+1,jPT+1)-1);
            }
            if ( !fFolder.IsNull() )    {
                gROOT->SetBatch();
                TCanvas    *cDrawResult =   new TCanvas();
                gStyle->SetOptStat(1111);
                hResults.at(iTH1).first->Draw();
                cDrawResult->SaveAs((fFolder+fName+TString(Form("_%i.pdf",iTH1))).Data());
                gStyle->SetOptStat(0);
                gROOT->SetBatch(kFALSE);
            }
            iTH1++;
            if ( (iPT != jPT) && fDiag )    {
                hResults.push_back( std::pair<TH1F*,Int_t>( new TH1F(Form("VarHist_%i_%i",jPT,iPT),Form("VarHist_%i_%i",jPT,iPT),2000,-1,+1) , hStandard->GetBin( jPT+1, iPT+1 ) ) );
                for ( auto  hVariationSpectrum : uSystVariationSpectra_ ) {
                    hResults.at(iTH1).first->Fill(hVariationSpectrum->GetBinContent(iPT+1,jPT+1)-1);
                }
                iTH1++;
            }
        }
    }
    return hResults;
}

//_____________________________________________________________________________
THStack*
uBuildSystematicStack
( TH1F* hStandard, std::vector<TH1F*> hVariations, std::vector<Bool_t> fWhichToChoose = std::vector<Bool_t> () ) {
    auto uVariationsBinByBin = uBuildVariationBinByBin( hStandard, hVariations, "", "", fWhichToChoose );
    THStack    *hResults    =   new THStack ("uBuildSystematicStack","uBuildSystematicStack");
    TH1F       *hMean       =   (TH1F*)hStandard->Clone();
    hMean                   ->  SetFillColorAlpha(2,0.1);
    hMean                   ->  SetLineColorAlpha(2,1.0);
    hMean                   ->  SetLineWidth(1);
    TH1F       *hSTDV       =   (TH1F*)hStandard->Clone();
    hSTDV                   ->  SetFillColorAlpha(4,0.1);
    hSTDV                   ->  SetLineColorAlpha(4,1.0);
    hSTDV                   ->  SetLineWidth(1);
    auto    iBin    =   0;
    for ( auto uBinVariation : uVariationsBinByBin ) {
        iBin++;
        hMean   ->  SetBinContent   ( iBin, fabs(uBinVariation->GetMean()) );
        hMean   ->  SetBinError     ( iBin, 0   );
        hSTDV   ->  SetBinContent   ( iBin, uBinVariation->GetRMS()        );
        hSTDV   ->  SetBinError     ( iBin, 0   );
    }
    hMean           ->  Scale(100);
    hSTDV           ->  Scale(100);
    hResults        ->  Add(hMean);
    hResults        ->  Add(hSTDV);
    
    return  hResults;
}
std::vector<THStack*>
uBuildSystematicStack
( TH2F* hStandard, std::vector<TH2F*> hVariations, std::vector<Bool_t> fWhichToChoose = std::vector<Bool_t> ()  ) {
    std::vector<THStack*>   hResults;
    for ( Int_t iPT2D = 1; iPT2D <= hStandard->GetNbinsY(); iPT2D++ ) {
        TH1F * hTargetSTD = (TH1F*)(hStandard->ProjectionX(Form("%i",iPT2D),iPT2D,iPT2D))->Clone();
        std::vector<TH1F*>   hVariationSlice;
        for ( auto uSlice : hVariations )   {
            hVariationSlice.push_back((TH1F*)(uSlice->ProjectionX(Form("%i_2",iPT2D),iPT2D,iPT2D))->Clone());
        }
        hResults.push_back(uBuildSystematicStack( hTargetSTD, hVariationSlice, fWhichToChoose ));
    }
    return hResults;
}

//_____________________________________________________________________________
TH1F*
uBuildSystematicError
( TH1F* hStandard, std::vector<TH1F*> hVariations, TString fFolder = "", TString fName = "", std::vector<Bool_t> fWhichToChoose = std::vector<Bool_t> () ) {
    auto uVariationsBinByBin = uBuildVariationBinByBin( hStandard, hVariations, fFolder, fName, fWhichToChoose );
    TH1F   *hResults    =   (TH1F*)(hStandard->Clone());
    auto    iBin    =   0;
    for ( auto uBinVariation : uVariationsBinByBin ) {
        iBin++;
        hResults    ->  SetBinContent   ( iBin, fabs(uBinVariation->GetMean()) + uBinVariation->GetRMS() );
        hResults    ->  SetBinError     ( iBin, 0   );
    }
    return  hResults;
}
TH2F*
uBuildSystematicError
( TH2F* hStandard, std::vector<TH2F*> hVariations, TString fFolder = "", TString fName = "", std::vector<Bool_t> fWhichToChoose = std::vector<Bool_t> (), Bool_t fDiag = false ) {
    auto uVariationsBinByBin = uBuildVariationBinByBin( hStandard, hVariations, fFolder, fName, fWhichToChoose, fDiag );
    TH2F   *hResults    =   (TH2F*)(hStandard->Clone());
    hResults->Reset();
    for ( auto uBinVariation : uVariationsBinByBin ) {
        hResults    ->  SetBinContent   ( uBinVariation.second, fabs(uBinVariation.first->GetMean()) + uBinVariation.first->GetRMS() );
        hResults    ->  SetBinError     ( uBinVariation.second, 0   );
    }
    return  hResults;
}

//_____________________________________________________________________________
TCanvas*
uBuildSystAndStatStack
( TH1F* hStandard, THStack* hSystematical )  {
    TCanvas    *cResult =   new TCanvas();
    TGraphErrors   *gStatistical        =   new TGraphErrors();
    hSystematical->Draw("HIST F");
    gStatistical                        ->  SetFillColorAlpha(kGray,0.75);
    for ( Int_t iBin = 0; iBin < hStandard->GetNbinsX(); iBin++ )  {
        gStatistical->SetPoint      ( iBin, hStandard->GetBinCenter(iBin+1), 0. );
        gStatistical->SetPointError ( iBin, .5*hStandard->GetBinWidth(iBin+1), hStandard->GetBinError(iBin+1)/hStandard->GetBinContent(iBin+1) );
    }
    TMultiGraph    *cDrawAllGraphs      =   new TMultiGraph("cDrawAllGraphs","");
    //g1D_Syst_Err                              ->  SetFillColorAlpha(kGray+2,0.5);
    cDrawAllGraphs                      ->  Add         (gStatistical,      "AE2");
    cDrawAllGraphs                      ->  Add         (gStatistical,      "AE2");
    //cDrawAllGraphs                      ->  Add         (g1D_Syst_Err,      "AE2");
    cDrawAllGraphs->Draw("SAME ALP");
    cDrawAllGraphs->SetMinimum(0.);
    hSystematical->Draw("SAME HIST F");
    return  cResult;
}

//_____________________________________________________________________________ !TODO: !
TCanvas*
uBuildTotalSystematics
 ( std::vector<TH1F*> hSystematics ) {
    TCanvas    *cDrawResults    =   new TCanvas();
    TH1F*   hTotal  =   (TH1F*)hSystematics.at(0)->Clone();
    return  cDrawResults;
}

//_____________________________________________________________________________
void
uEvaluateRatioError
 ( TH1F* h1DStandard, std::vector<TH1F*> h1DVariations, TH2F* h2DStandard, std::vector<TH2F*> h2DVariations, TString fFolder = "", TH1F* h1DEfficiency = nullptr, TH2F* h2DEfficiency = nullptr ) {
    auto f1DRelevantVariations  =   uIsRelevantVariation(   h1DStandard,  h1DVariations,    (fFolder+TString("/plots/BarlowCheck/1D/")).Data(),"1D");
    auto f2DRelevantVariations  =   uIsRelevantVariation(   h2DStandard,  h2DVariations,    (fFolder+TString("/plots/BarlowCheck/2D/")).Data(),"2D",true);
    auto h1DFullSystematics     =   uBuildSystematicError(  h1DStandard,  h1DVariations,    (fFolder+TString("/plots/BinByBinCheck/1D/")).Data(),"1D",  f1DRelevantVariations);
    auto h2DFullSystematics     =   uBuildSystematicError(  h2DStandard,  h2DVariations,    (fFolder+TString("/plots/BinByBinCheck/2D/")).Data(),"2D",  f2DRelevantVariations,    true);
    //
    if ( h1DEfficiency ) h1DStandard->Divide(h1DEfficiency);
    if ( h2DEfficiency ) h2DStandard->Divide(h2DEfficiency);
    for ( Int_t iX = 0; iX < h1DStandard->GetNbinsX(); iX++ )   {
        h1DStandard->SetBinError(iX+1,h1DFullSystematics->GetBinContent(iX+1)*h1DStandard->GetBinContent(iX+1));
    }
    for ( Int_t iX = 0; iX < h2DStandard->GetNbinsX(); iX++ )   {
        for ( Int_t iY = 0; iY < h2DStandard->GetNbinsY(); iY++ )   {
            h2DStandard->SetBinError(iX+1,iY+1,h2DFullSystematics->GetBinContent(iX+1,iY+1)*h2DStandard->GetBinContent(iX+1,iY+1));
        }
    }
    //
    auto    iTer    =   0;
    std::vector<Float_t>   kSimpleRatio;
    std::vector<Float_t>   kSquareRatio;
    std::vector<Float_t>   k1Phi_Sng;
    std::vector<Float_t>   k2Phi_Sng;
    auto    k1DStdIntegralErr  =   0.;
    auto    k2DStdIntegralErr  =   0.;
    auto    k1DStdIntegral     =   h1DStandard->IntegralAndError(-1,10000,k1DStdIntegralErr,"width");
    auto    k2DStdIntegral     =   h2DStandard->IntegralAndError(-1,10000,-1,10000,k2DStdIntegralErr,"width");
            k1DStdIntegralErr /=   k1DStdIntegral;
            k2DStdIntegralErr /=   k2DStdIntegral;
    //
    for ( auto hVariation : h1DVariations )   {
        //
        if ( h1DEfficiency ) h1DVariations.at(iTer)->Divide(h1DEfficiency);
        if ( h2DEfficiency ) h2DVariations.at(iTer)->Divide(h2DEfficiency);
        auto    k1Dintegral     =   h1DVariations.at(iTer)->Integral("width");
        auto    k2Dintegral     =   h2DVariations.at(iTer)->Integral("width");
        //
        k1Phi_Sng.push_back(k1Dintegral/k1DStdIntegral -1);
        k2Phi_Sng.push_back(k2Dintegral/k2DStdIntegral -1);
        kSimpleRatio.push_back((k1DStdIntegral*k2Dintegral)/(k1Dintegral*k2DStdIntegral)-1);
        kSquareRatio.push_back((k1DStdIntegral*k1DStdIntegral*k2Dintegral)/(k1Dintegral*k1Dintegral*k2DStdIntegral)-1);
        iTer++;
    }
    //
    TH1F       *hSimpleRatioError   =   uBuildTH1F(kSimpleRatio,2000,0,-0.5,0.5);
    auto        fSimpleRatioError   =   0.;
    TH1F       *hSquareRatioError   =   uBuildTH1F(kSquareRatio,2000,0,-0.5,0.5);
    auto        fSquareRatioError   =   0.;
    TH1F       *h1Phi_Sng_Error     =   uBuildTH1F(k1Phi_Sng,2000,0,-0.5,0.5);
    auto        f1Phi_Sng_Error     =   0.;
    TH1F       *h2Phi_Sng_Error     =   uBuildTH1F(k2Phi_Sng,2000,0,-0.5,0.5);
    auto        f2Phi_Sng_Error     =   0.;
    fSimpleRatioError   +=  fabs(hSimpleRatioError->GetMean());
    fSimpleRatioError   +=  hSimpleRatioError->GetRMS();
    fSquareRatioError   +=  fabs(hSquareRatioError->GetMean());
    fSquareRatioError   +=  hSquareRatioError->GetRMS();
    f1Phi_Sng_Error     +=  fabs(h1Phi_Sng_Error->GetMean());
    f1Phi_Sng_Error     +=  h1Phi_Sng_Error->GetRMS();
    f2Phi_Sng_Error     +=  fabs(h2Phi_Sng_Error->GetMean());
    f2Phi_Sng_Error     +=  h2Phi_Sng_Error->GetRMS();
    //
    TH1F   *hStandard   =   new TH1F("hStandard",   "", 2,  0,  2);
    hStandard->GetYaxis()->SetTitle("Systematic uncertainty (%)");
    hStandard->GetXaxis()->SetNdivisions(2);
    hStandard->GetXaxis()->SetBinLabel(hStandard->GetXaxis()->FindBin(0.5),"#LT Y_{2#phi} #GT / #LT Y_{1#phi} #GT");
    hStandard->GetXaxis()->SetBinLabel(hStandard->GetXaxis()->FindBin(1.5),"#LT Y_{2#phi} #GT / #LT Y_{1#phi} #GT^{2}");
    hStandard->GetXaxis()->LabelsOption("h");
    hStandard->SetMarkerColor(kColors[3]);
    hStandard->SetLineWidth(3);
    hStandard->SetMarkerStyle(kMarkers[3]);
    hStandard->SetBinContent(1,fSimpleRatioError);
    hStandard->SetBinContent(2,fSquareRatioError);
    hStandard->SetBinError  (1,0);
    hStandard->SetBinError  (2,0);
    hStandard->Scale(100);
    TH1F   *hStandard_Sng   =   new TH1F("hStandard_Sng",   "", 2,  0,  2);
    hStandard_Sng->GetYaxis()->SetTitle("Systematic uncertainty (%)");
    hStandard_Sng->GetXaxis()->SetNdivisions(2);
    hStandard_Sng->GetXaxis()->SetBinLabel(hStandard->GetXaxis()->FindBin(0.5),"#LT Y_{#phi} #GT");
    hStandard_Sng->GetXaxis()->SetBinLabel(hStandard->GetXaxis()->FindBin(1.5),"#LT Y_{#phi#phi} #GT");
    hStandard_Sng->GetXaxis()->LabelsOption("h");
    hStandard_Sng->SetMarkerColor(kColors[3]);
    hStandard_Sng->SetLineWidth(3);
    hStandard_Sng->SetMarkerStyle(kMarkers[3]);
    hStandard_Sng->SetBinContent(1,f1Phi_Sng_Error);
    hStandard_Sng->SetBinContent(2,f2Phi_Sng_Error);
    hStandard_Sng->SetBinError  (1,0);
    hStandard_Sng->SetBinError  (2,0);
    hStandard_Sng->Scale(100);
    TH1F   *hSimple     =   new TH1F("hLinear",     "", 2,  0,  2);
    hSimple->SetMarkerColor(kColors[2]);
    hSimple->SetLineWidth(3);
    hSimple->SetMarkerStyle(kMarkers[4]);
    hSimple->SetBinContent(1,k1DStdIntegralErr+k2DStdIntegralErr);
    hSimple->SetBinContent(2,2*k1DStdIntegralErr+k2DStdIntegralErr);
    hSimple->SetBinError  (1,0);
    hSimple->SetBinError  (2,0);
    hSimple->Scale(100);
    TH1F   *hSimple_Sng     =   new TH1F("hLinear_Sng",     "", 2,  0,  2);
    hSimple_Sng->SetMarkerColor(kColors[2]);
    hSimple_Sng->SetLineWidth(3);
    hSimple_Sng->SetMarkerStyle(kMarkers[4]);
    hSimple_Sng->SetBinContent(1,k1DStdIntegralErr);
    hSimple_Sng->SetBinContent(2,k2DStdIntegralErr);
    hSimple_Sng->SetBinError  (1,0);
    hSimple_Sng->SetBinError  (2,0);
    hSimple_Sng->Scale(100);
    TH1F   *hSquare     =   new TH1F("hSquare",     "", 2,  0,  2);
    hSquare->SetMarkerColor(kColors[1]);
    hSquare->SetLineWidth(3);
    hSquare->SetMarkerStyle(kMarkers[5]);
    hSquare->SetBinContent(1,SquareSum( {k1DStdIntegralErr,k2DStdIntegralErr} ));
    hSquare->SetBinContent(2,SquareSum( {2*k1DStdIntegralErr,k2DStdIntegralErr} ));
    hSquare->SetBinError  (1,0);
    hSquare->SetBinError  (2,0);
    hSquare->Scale(100);
    //
    auto fMaximum  = 100 * max ( 2*k1DStdIntegralErr+k2DStdIntegralErr, max (fSimpleRatioError, fSquareRatioError ) );
    hStandard->SetMaximum(fMaximum*1.3);
    fMaximum  = 100 * max ( k2DStdIntegralErr, f2Phi_Sng_Error );
    hStandard_Sng->SetMaximum(fMaximum*1.3);
    //
    gROOT->SetBatch();
    TCanvas    *c1 = new TCanvas("Ratio");
    //
    TLegend    *lLegend = new TLegend(0.18,0.82,0.33,0.72);
    lLegend     ->  AddEntry( hStandard,    "Ratio err.", "P" );
    lLegend     ->  AddEntry( hSimple,      "Linear err.", "P" );
    lLegend     ->  AddEntry( hSquare,      "Square err.", "P" );
    //
    hStandard->Draw("][ EP MIN0");
    hSimple->Draw("SAME EP ][");
    hSquare->Draw("SAME EP ][");
    lLegend->Draw("same");
    //
    c1->SaveAs((fFolder+TString("/plots/Ratio.pdf")).Data());
    delete c1;
    delete lLegend;
    //
                c1      = new TCanvas("Ratio");
    //
                lLegend = new TLegend(0.18,0.82,0.33,0.72);
    lLegend     ->  AddEntry( hStandard,    "Yield err.", "P" );
    lLegend     ->  AddEntry( hSimple,      "Integral err.", "P" );
    //
    hStandard_Sng->Draw("][ EP MIN0");
    hSimple_Sng->Draw("SAME EP ][");
    lLegend->Draw("same");
    //
    c1->SaveAs((fFolder+TString("/plots/Ratio_Sng.pdf")).Data());
    delete c1;
    //
    TCanvas    *c2 = new TCanvas("Ratio_Simple");
    //
    hSimpleRatioError->Draw();
    //
    c2->SaveAs((fFolder+TString("/plots/Ratio_Simple.pdf")).Data());
    delete c2;
    //
    TCanvas    *c3 = new TCanvas("Ratio_Square");
    //
    hSquareRatioError->Draw();
    //
    c3->SaveAs((fFolder+TString("/plots/Ratio_Square.pdf")).Data());
    delete c3;
    //
    gROOT->SetBatch(kFALSE);
    //
    TFile      *fOutput =   new TFile   (Form("%s%s",fFolder.Data(),"/RT_Systematic.root"),"recreate");
    hStandard->Scale(0.01);
    hStandard->Write();
    fOutput->Close();
}

//
void
uPlotStack
 ( THStack* hTarget, TString fFolder, TString fSubFolder = "", TString fName = "" )   {
    gROOT->SetBatch();
    //
    TCanvas*    cDrawStackSystematics   =   new TCanvas();
    gStyle      ->  SetOptStat(0);
    gPad        ->  SetLogx();
    //
    TLegend    *lLegend =   new TLegend(0.2,0.85,0.4,0.7);
    hTarget    ->  Draw();
    uSetHisto   ( hTarget, "SYSSTACK" );
    //
    lLegend->AddEntry((TH1F*)hTarget->GetHists()->At(0),"#mu contr.","F");
    lLegend->AddEntry((TH1F*)hTarget->GetHists()->At(1),"#sigma contr.","F");
    //
    hTarget    ->  Draw();
    lLegend             ->  Draw("SAME");
    //
    cDrawStackSystematics   ->  SaveAs((fFolder+TString("/plots/")+fSubFolder+TString("/Full_")+ fName + TString("_Sys.pdf")).Data());
    gROOT->SetBatch(kFALSE);
}
//
void
uPlotStack
 ( std::vector<THStack*> hTarget, TString fFolder, TString fSubFolder = "", TString fName = "" )   {
    gROOT->SetBatch();
    //
    TCanvas*    cDrawStackSystematics   =   new TCanvas();
    gStyle      ->  SetOptStat(0);
    gPad        ->  SetLogx();
    //
    TLegend    *lLegend =   new TLegend(0.2,0.85,0.4,0.7);
    lLegend     ->  AddEntry((TH1F*)hTarget.at(0)->GetHists()->At(0),"#mu contr.","F");
    lLegend     ->  AddEntry((TH1F*)hTarget.at(0)->GetHists()->At(1),"#sigma contr.","F");
    auto    iStack = 0;
    for ( auto kCurrent_Stack : hTarget )  {
        kCurrent_Stack      ->  Draw();
        uSetHisto   ( kCurrent_Stack, "SYSSTACK" );
        //
        kCurrent_Stack    ->  Draw();
        lLegend             ->  Draw("SAME");
        //
        cDrawStackSystematics   ->  SaveAs((fFolder+TString("/plots/")+fSubFolder+TString("/Full_")+ fName + TString(Form("_%i_Sys.pdf",iStack))).Data());
        iStack++;
    }
    gROOT->SetBatch(kFALSE);
}
//
template <  typename THXTarget_Type >
THXTarget_Type*
uSysEvaluate_BinByBin
( THXTarget_Type* hStandard, std::vector<THXTarget_Type*> hVariations, TString fFolder, TString fSubFolder = "", TString fName = "", Bool_t kNoBarlowCheck = false  ) {
    //
    //  --- Silence Warnings
    gErrorIgnoreLevel   =   kWarning;
    //
    //  --- Create Substructure
    if ( !fSubFolder.IsNull() ) {
        gROOT   ->  ProcessLine ( Form(".! mkdir -p %s/",(fFolder+TString("/plots/")+fSubFolder+TString("/BarlowCheck/")+fName+TString("/")).Data()) );
        gROOT   ->  ProcessLine ( Form(".! mkdir -p %s/",(fFolder+TString("/plots/")+fSubFolder+TString("/BinByBinCheck/")+fName+TString("/")).Data()) );
        gROOT   ->  ProcessLine ( Form(".! mkdir -p %s/",(fFolder+TString("/plots/")+fSubFolder+fName+TString("/")).Data()) );
    } else {
        gROOT   ->  ProcessLine ( Form(".! mkdir -p %s/",(fFolder+TString("/plots/BarlowCheck/")+fName+TString("/")).Data()) );
        gROOT   ->  ProcessLine ( Form(".! mkdir -p %s/",(fFolder+TString("/plots/BinByBinCheck/")+fName+TString("/")).Data()) );
    }
    //
    //  --- Check for relevant variations if requested
    std::vector<Bool_t> kIsRelevantVariation;
    if ( kNoBarlowCheck )   for ( auto kNull : hVariations )  kIsRelevantVariation.push_back(true);
    else                    kIsRelevantVariation = uIsRelevantVariation(hStandard,hVariations,(fFolder+TString("/plots/")+fSubFolder+TString("/BarlowCheck/")+fName+TString("/")).Data(),fName);
    //
    //  --- Build Systematical Stacks
    auto    uStackSystematic    =   uBuildSystematicStack( hStandard, hVariations, kIsRelevantVariation );
    //
    //  --- Plot Results
    //
    uPlotStack( uStackSystematic, fFolder, fSubFolder, fName );
    //
    TFile      *fOutput =   new TFile   (Form("%s/%s",fFolder.Data(),Form("/%s_Systematic.root",fName.Data())),"recreate");
    uBuildSystematicError(hStandard,hVariations,(fFolder+TString("/plots/")+fSubFolder+TString("/BinByBinCheck/")+fName+TString("/")).Data(),fName,kIsRelevantVariation)->Write();
    fOutput->Close();
    //
    gErrorIgnoreLevel   =   kInfo;
    //
    return uBuildSystematicError(hStandard,hVariations,(fFolder+TString("/plots/")+fSubFolder+TString("/BinByBinCheck/")+fName+TString("/")).Data(),fName,kIsRelevantVariation);
}
//
template <  typename TH1Target_Type,
            typename TH2Target_Type >
void
uSysEvaluate_Extrapolation_Custom1D
( TH1Target_Type* hStandard_1D, std::vector<TH1Target_Type*> hVariations_1D, TH2Target_Type* hStandard_2D, std::vector<TH2Target_Type*> hVariations_2D, TString fFolder = "", TString fSubFolder = "", TString fName = "", Bool_t kNoBarlowCheck = false  ) {
    auto    k1DSystematics  =   uSysEvaluate_BinByBin( hStandard_1D, hVariations_1D, fFolder, fSubFolder, "1D", kNoBarlowCheck );
    auto    k2DSystematics  =   uSysEvaluate_BinByBin( hStandard_2D, hVariations_2D, fFolder, fSubFolder, "2D", kNoBarlowCheck );
    //
    //  --- Silence Warnings
    gErrorIgnoreLevel   =   kWarning;
    //gROOT->SetBatch();
    //
    //  --- Create Substructure
    TString kPlotDir    =   fFolder+TString("/plots/")+fSubFolder+TString("/");
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s/",(kPlotDir).Data()) );
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s/",(kPlotDir+TString("/1D/")).Data()) );
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s/",(kPlotDir+TString("/2D/")).Data()) );
    gROOT   ->  ProcessLine ( Form(".! mkdir -p %s/",(kPlotDir+TString("/Full/")).Data()) );
    //
    //  --- Check for relevant variations if requested
    std::vector<Bool_t> kIsRelevantVariation1D;
    std::vector<Bool_t> kIsRelevantVariation2D;
    if ( kNoBarlowCheck )   for ( auto kNull : hVariations_1D )  kIsRelevantVariation1D.push_back(true);
    else                    kIsRelevantVariation1D = uIsRelevantVariation(hStandard_1D,hVariations_1D);
    if ( kNoBarlowCheck )   for ( auto kNull : hVariations_2D )  kIsRelevantVariation2D.push_back(true);
    else                    kIsRelevantVariation2D = uIsRelevantVariation(hStandard_2D,hVariations_2D);
    //
    std::vector<Float_t>    kVariation_1D_Yield;
    std::vector<Float_t>    kVariation_2D_Yield;
    std::vector<Float_t>    kVariation_XD_YY_over_Y;
    std::vector<Float_t>    kVariation_XD_YY_over_Y_Y;
    std::vector<Float_t>    kVariation_XD_Sigma;
    std::vector<Float_t>    kVariation_XD_Gamma;
    //  --- Start Inclusive Yield Systematics
    fStartTimer("Inclusive Yield Systematics Evaluation");
    //
    auto    iTer = 0;
    std::map<TString,std::tuple<Float_t,Float_t,Float_t>>               kStandard_1D    =   uMeasureFullYield( hStandard_1D, uScale(hStandard_1D,1.,-2.), {kStandardSystematicFitFunctions.at(0)}, 2, kPlotDir+TString("/1D/"), "STD_1D" );
    std::vector<std::map<TString,std::tuple<Float_t,Float_t,Float_t>>>  kVariations_1D;
    for ( auto kCurrent_Variation : hVariations_1D )    {
        iTer++;
        fPrintLoopTimer( "Inclusive Yield Systematics Evaluation", iTer, hVariations_1D.size()+1, 1 );
        kVariations_1D.push_back( uMeasureFullYield( kCurrent_Variation, uScale(kCurrent_Variation,1.,-2.), {kStandardSystematicFitFunctions.at(0)}, 8, kPlotDir+TString("/1D/"), TString(kCurrent_Variation->GetName())+TString("_1D") ) );
    }
    //
    fStopTimer("Inclusive Yield Systematics Evaluation");
    //
    //  --- Start Inclusive Pair Yield Systematics
    fStartTimer("Inclusive Pair Yield Systematics Evaluation");
    //
    iTer = 0;
    std::vector<std::map<TString,std::tuple<Float_t,Float_t,Float_t>>>              kStandard_2D    =   uMeasureFullYield2D( hStandard_2D, uScale(hStandard_2D,1.,-2.), {kStandardSystematicFitFunctions.at(0)}, 8, kPlotDir+TString("/2D/"), "STD_2D" );
    std::vector<std::vector<std::map<TString,std::tuple<Float_t,Float_t,Float_t>>>> kVariations_2D;
    auto                hMPT_Standard   =   uBuildMeanPT<TH1F>(kStandard_2D);
    std::vector<TH1F*>  hMPT_Variations;
    for ( auto kCurrent_Variation : hVariations_2D )    {
        iTer++;
        fPrintLoopTimer( "Inclusive Pair Yield Systematics Evaluation", iTer, hVariations_2D.size()+1, 1 );
        kVariations_2D.push_back( uMeasureFullYield2D( kCurrent_Variation, uScale(kCurrent_Variation,1.,-2.), {kStandardSystematicFitFunctions.at(0)}, 8, kPlotDir+TString("/2D/"), TString(kCurrent_Variation->GetName())+TString("_2D_%i") ) );
        auto    hCurrent_MPT_Histo = uBuildMeanPT<TH1F>(kVariations_2D.at(iTer-1));
        hCurrent_MPT_Histo  ->  SetName( kCurrent_Variation->GetName() );
        hMPT_Variations.push_back( hCurrent_MPT_Histo );
    }
    //
    fStopTimer("Inclusive Pair Yield Systematics Evaluation");
    //
    auto    kStandard_1D_Yield  =   get<0>(kStandard_1D["YL_FLL"]);
    auto    kStandard_2D_Yield  =   get<0>(kStandard_2D.at(0)["YL_FLL"]);
    for ( Int_t iVar = 0; iVar < min( hVariations_1D.size(), hVariations_2D.size() ); iVar++ )  {
        if ( ( !kIsRelevantVariation1D.at(iVar) ) && ( !kIsRelevantVariation2D.at(iVar) ) ) continue;
        auto    kCurrent_1D_Yield   =   get<0>((kVariations_1D.at(iVar))["YL_FLL"]);
        auto    kCurrent_2D_Yield   =   get<0>((kVariations_2D.at(iVar)).at(0)["YL_FLL"]);
        if ( kIsRelevantVariation1D.at(iVar) )  {
            kVariation_1D_Yield.push_back( kCurrent_1D_Yield / kStandard_1D_Yield );
        }
        if ( kIsRelevantVariation2D.at(iVar) )  {
            kVariation_2D_Yield.push_back( kCurrent_2D_Yield / kStandard_2D_Yield );
        }
        kVariation_XD_YY_over_Y     .push_back( ( kCurrent_2D_Yield * kStandard_1D_Yield ) / ( kStandard_2D_Yield * kCurrent_1D_Yield ) );
        kVariation_XD_YY_over_Y_Y   .push_back( ( kCurrent_2D_Yield * kStandard_1D_Yield * kStandard_1D_Yield ) / ( kStandard_2D_Yield * kCurrent_1D_Yield * kCurrent_1D_Yield ) );
        kVariation_XD_Sigma         .push_back( fSigmaPhiValue( kCurrent_1D_Yield, kCurrent_2D_Yield ) / fSigmaPhiValue( kStandard_1D_Yield, kStandard_2D_Yield ) );
        kVariation_XD_Gamma         .push_back( fGammaPhiValue( kCurrent_1D_Yield, kCurrent_2D_Yield ) / fGammaPhiValue( kStandard_1D_Yield, kStandard_2D_Yield ) );
    }
    //
    //  --- Start Comparison Systematics
    fStartTimer("Comparison Systematics Evaluation");
    //
    auto    h1D_Stat    =   (TH1F*)hStandard_1D->Clone("h1D_Stat");
    auto    h2D_Stat    =   (TH2F*)hStandard_2D->Clone("h2D_Stat");
    auto    h1D_Syst    =   (TH1F*)hStandard_1D->Clone("h1D_Syst");
    auto    h2D_Syst    =   (TH2F*)hStandard_2D->Clone("h2D_Syst");
    for ( Int_t iBin = 1; iBin <= hStandard_1D->GetNbinsX(); iBin++ )    {
        h1D_Syst    ->  SetBinError( iBin, k1DSystematics->GetBinContent(iBin)*hStandard_1D->GetBinContent(iBin) );
    }
    auto    k1DPropagated   =   uMeasureFullYield( h1D_Stat, h1D_Syst, {kStandardSystematicFitFunctions.at(0)}, 8, kPlotDir+TString("/1D/"), "PRG_1D" );
    for ( Int_t iBin = 1; iBin <= hStandard_2D->GetNbinsX(); iBin++ )    {
        for ( Int_t jBin = 1; jBin <= hStandard_2D->GetNbinsY(); jBin++ ) {
            h2D_Syst    ->  SetBinError( iBin, jBin, k2DSystematics->GetBinContent(iBin,jBin)*hStandard_2D->GetBinContent(iBin,jBin) );
        }
    }
    auto k2DPropagated  = uMeasureFullYield2D( h2D_Stat, h2D_Syst, {kStandardSystematicFitFunctions.at(0)}, 8, kPlotDir+TString("/2D/"), "_PRG_2D_%i" );
    //
    fStopTimer("Comparison Systematics Evaluation");
    //
    TH1F*   kRatioPropagated_   =   new TH1F( "kRatioPropagated_", "kRatioPropagated_", 6, 0.5, 6.5 );
    uSetHisto(kRatioPropagated_,"FNL");
    auto    k1DError    =   get<2>(k1DPropagated["YL_FLL"]) / get<0>(k1DPropagated["YL_FLL"]);
    auto    k2DError    =   get<2>(k2DPropagated.at(0)["YL_FLL"]) / get<0>(k2DPropagated.at(0)["YL_FLL"]);
    auto    kR1Error    =   SquareSum( { k2DError, k1DError } );
    auto    kR2Error    =   SquareSum( { k2DError, k1DError, k1DError } );
    auto    kP1Error    =   fSigmaPhiError( get<0>(k1DPropagated["YL_FLL"]), get<0>(k2DPropagated.at(0)["YL_FLL"]), get<2>(k1DPropagated["YL_FLL"]), get<2>(k2DPropagated.at(0)["YL_FLL"]) ) / fSigmaPhiValue( get<0>(k1DPropagated["YL_FLL"]), get<0>(k2DPropagated.at(0)["YL_FLL"]) );
    auto    kP2Error    =   fGammaPhiError( get<0>(k1DPropagated["YL_FLL"]), get<0>(k2DPropagated.at(0)["YL_FLL"]), get<2>(k1DPropagated["YL_FLL"]), get<2>(k2DPropagated.at(0)["YL_FLL"]) ) / fGammaPhiValue( get<0>(k1DPropagated["YL_FLL"]), get<0>(k2DPropagated.at(0)["YL_FLL"]) );
    kRatioPropagated_   ->  SetBinContent( 1, k1DError );
    kRatioPropagated_   ->  SetBinContent( 2, k2DError );
    kRatioPropagated_   ->  SetBinContent( 3, kR1Error );
    kRatioPropagated_   ->  SetBinContent( 4, kR2Error );
    kRatioPropagated_   ->  SetBinContent( 5, kP1Error );
    kRatioPropagated_   ->  SetBinContent( 6, kP2Error );
    //
    TH1F*   kRatioSystematics   =   new TH1F( "kRatioSystematics", "kRatioSystematics", 6, 0.5, 6.5 );
    uSetHisto(kRatioSystematics,"FNL");
    auto    h1DYield    =   uBuildTH1<TH1F>( kVariation_1D_Yield,       250, -1 );
    auto    h2DYield    =   uBuildTH1<TH1F>( kVariation_2D_Yield,       250, -1 );
    auto    hXDRati1    =   uBuildTH1<TH1F>( kVariation_XD_YY_over_Y,   250, -1 );
    auto    hXDRati2    =   uBuildTH1<TH1F>( kVariation_XD_YY_over_Y_Y, 250, -1 );
    auto    hXDPara1    =   uBuildTH1<TH1F>( kVariation_XD_Sigma,       250, -1 );
    auto    hXDPara2    =   uBuildTH1<TH1F>( kVariation_XD_Gamma,       250, -1 );
    //
    kRatioSystematics   ->  SetBinContent( 1, fabs(h1DYield->GetMean()) + h1DYield->GetRMS() );
    kRatioSystematics   ->  SetBinContent( 2, fabs(h2DYield->GetMean()) + h2DYield->GetRMS() );
    kRatioSystematics   ->  SetBinContent( 3, fabs(hXDRati1->GetMean()) + hXDRati1->GetRMS() );
    kRatioSystematics   ->  SetBinContent( 4, fabs(hXDRati2->GetMean()) + hXDRati2->GetRMS() );
    kRatioSystematics   ->  SetBinContent( 5, fabs(hXDPara1->GetMean()) + hXDPara1->GetRMS() );
    kRatioSystematics   ->  SetBinContent( 6, fabs(hXDPara2->GetMean()) + hXDPara2->GetRMS() );
    //
    TCanvas*    cDrawComparison =   new TCanvas( "cDrawComparison", "cDrawComparison", 1500, 1500 );
    //
    TLegend*    cLegendComparison   =   new TLegend( 0.3,0.75,0.7,0.85 );
    cLegendComparison   ->  SetNColumns(2);
    cLegendComparison   ->  AddEntry( kRatioPropagated_, "Propagated" );
    cLegendComparison   ->  AddEntry( kRatioSystematics, "Calculated" );
    //
    kRatioPropagated_->SetLineColor(kBlue);
    kRatioPropagated_->SetLineWidth(2);
    kRatioPropagated_->SetMaximum( 1.4*max( 1.*kRatioPropagated_->GetMaximum(), 1.*kRatioSystematics->GetMaximum() ) );
    kRatioPropagated_->Draw("SAME");
    //
    kRatioSystematics->SetLineColor(kRed);
    kRatioSystematics->SetLineWidth(2);
    kRatioSystematics->Draw("SAME");
    //
    cLegendComparison->Draw("SAME");
    //
    cDrawComparison ->  SaveAs( kPlotDir+TString("/Full/ComparePropagated.pdf") );
    //
    delete cDrawComparison;
    //
    auto kPTSystematics =   uSysEvaluate_BinByBin( hMPT_Standard, hMPT_Variations, fFolder, fSubFolder, "MPT", kNoBarlowCheck );
    kPTSystematics->SetName("kPTSystematics");
    //
    gROOT->SetBatch(kFALSE);
    //
    TFile*  fOutFile    =   new TFile( fFolder+TString("/FullSystematics.root"), "recreate" );
    //
    k1DSystematics      ->  SetName("k1DSystematics");
    k2DSystematics      ->  SetName("k2DSystematics");
    kRatioSystematics   ->  SetName("kRatioSystematics");
    kPTSystematics      ->  SetName("kPTSystematics");
    k1DSystematics      ->  Write();
    k2DSystematics      ->  Write();
    kRatioSystematics   ->  Write();
    kPTSystematics      ->  Write();
    h1D_Syst            ->  Write();
    h2D_Syst            ->  Write();
    //
    fOutFile->Close();
}
//
template <  Bool_t      TSquareSum          = kTRUE,
            typename    THXTarget_Type_1    = TH1F,
            typename    THXTarget_Type_2    = TH1F  >
THXTarget_Type_1*
uSumSystErrors
 ( THXTarget_Type_1* hTarget_1, THXTarget_Type_2* hTarget_2 )    {
    THXTarget_Type_1*   fResult =   (THXTarget_Type_1*)(hTarget_1->Clone());
    if ( !uIsTHPairConsistent( hTarget_1, hTarget_2 ) ) return fResult;
    for ( Int_t iBin = 1; iBin <= fResult->GetNbinsX(); iBin++ ) {
        for ( Int_t jBin = 1; jBin <= fResult->GetNbinsY(); jBin++ ) {
            for ( Int_t kBin = 1; kBin <= fResult->GetNbinsZ(); kBin++ ) {
                auto    kGlobalBin  =   fResult->GetBin( iBin, jBin, kBin );
                auto    kNewBinCn   =   0.;
                if ( TSquareSum )   kNewBinCn   =   SquareSum( { hTarget_1->GetBinContent(kGlobalBin), hTarget_2->GetBinContent(kGlobalBin) } );
                else                kNewBinCn   =   hTarget_1->GetBinContent(kGlobalBin) + hTarget_2->GetBinContent(kGlobalBin);
                fResult ->  SetBinContent( kGlobalBin, kNewBinCn );
            }
        }
    }
    return  fResult;
}
//
template <  Bool_t      TSquareSum      = kTRUE,
            typename    THXTarget_Type  = TH1F >
THXTarget_Type*
uSumSystErrors
 ( std::vector<THXTarget_Type*> hTarget )    {
    THXTarget_Type*   fResult =   (THXTarget_Type*)(hTarget.at(0)->Clone());
    Bool_t  bSkipFirst = true;
    for ( auto kCurrentHisto : hTarget )    {
        if ( bSkipFirst ) {
            bSkipFirst = false;
            continue;
        }
        fResult =   uSumSystErrors<TSquareSum,THXTarget_Type>( fResult, kCurrentHisto );
    }
    return  fResult;
}
//
#endif /* AAU_Systematics_h */
