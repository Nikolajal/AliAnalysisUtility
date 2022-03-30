//
//  Part of the AliAnalysisUtility package
//
//  Global Variables
//
//  Author              Nicola Rubini
//  Created             22/11/2021
//  Last modified       22/11/2021
#ifndef ALIANALYSISUTILITY_GLOBALUTILVARIABLES_H
#define ALIANALYSISUTILITY_GLOBALUTILVARIABLES_H
//
//  --- GLOBAL VARIABLES
//
//  --- --- Random Generator
TRandom            *uRandomGen                      =   new TRandom();
//  --- --- Benchmark
TBenchmark         *uBenchmark                      =   new TBenchmark();
//  --- --- TLatex
TLatex             *uLatex                          =   new TLatex();
//
//  --- --- Title and Name for histograms
auto                hName                           =   "Name";
auto                hTitle                          =   "Title";
//
//  --- --- General Parameters
const Bool_t        kFitScarseHisto                 =   kTRUE;     //  Skip the fit of histograms that are below the threshold set below
const Float_t       kScarseHistoDef                 =   0.;         //  % of entries w.r.t. total bin number
const Int_t         kScarseHistoMin                 =   1000.;      //  N of entries
const Float_t       kLooseErrors                    =   3.;         //  Multiplication factor for the Error Loosening
const Double_t      kMaximumError                   =   2.5;        //  Maximum percentage error to accept fit result, will retry if higher
//
//  --- --- Physics is by defualt in GeV, we define some constants to evaluate results in other units
//  --- --- When inserting the values into the system use *MeV, when extracting from the system use /MeV
//  --- --- To change this, use the constants below
auto const          KeV                             =   1e-6;
auto const          MeV                             =   1e-3;
auto const          GeV                             =   1;
auto const          TeV                             =   1e+3;
//
//  --- --- BenchMark Utility
TString                 kMSG_PrintTimer             =   "\r[INFO] [%s] Event # %7.f %s | %3.1f %% | %7.2f %s events/s | Time: %02.0f:%02.0f:%02.0f | ETA: %02.0f:%02.0f:%02.0f [%s]";
//
//  --- --- Scale for restricted Gauss Integrals
Float_t const       kGaussIntegralScale []          =   {0.68268949,0.95449974,0.99730020,0.99993666,0.99999943,1.00000000};
Float_t const       kGaussStndDevtScale []          =   {1.85349570,1.13687290,1.01350180,1.00048960,1.00000000,1.00000000};
//
//  --- --- Standard Legend for Resolution control plot
std::vector<TString>    kResolutionLegend           =   {"RMS 1 #sigma","RMS 2 #sigma","RMS 3 #sigma","RMS 4 #sigma","RMS 5 #sigma","RMS 6 #sigma","Gaus 1 #sigma","Gaus 2 #sigma","Gaus 3 #sigma","Gaus 4 #sigma","Gaus 5 #sigma","Gaus 6 #sigma","True Mass Fit"};
//
#endif
