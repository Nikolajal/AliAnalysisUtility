//
//  Re-Include protection
#ifndef ALIANALYSISUTILITY_GLOBALUTILVARIABLES_H
#define ALIANALYSISUTILITY_GLOBALUTILVARIABLES_H
//
//>>    GLOBAL VARIABLES
//
// Random Generator
TRandom            *uRandomGen                      =   new TRandom();
// Benchmark
TBenchmark         *uBenchmark                      =   new TBenchmark();
// TLatex
TLatex             *uLatex                          =   new TLatex();

// Title and Name for histograms
auto                hName                           =   "Name";
auto                hTitle                          =   "Title";
// General Parameters
const Bool_t        kFitScarseHisto                 =   kTRUE;      //  Skip the fit of histograms that are below the threshold set below
const Float_t       kScarseHistoDef                 =   0.;         //  % of entries w.r.t. total bin number
const Int_t         kScarseHistoMin                 =   1000.;      //  N of entries
const Float_t       kLooseErrors                    =   3.;         //  Multiplication factor for the Error looosening
const Double_t      kMaximumError                   =   2.5;        //  Maximum percentage error to accept fit result, will retry if higher
//
//>> Physics is by defualt in GeV, we define some constants to evaluate results in other units
//>> When inserting the values into the system use *MeV, when extracting from the system use /MeV
//>> To change this, use the constants below
auto const          KeV                             =   1e-6;
auto const          MeV                             =   1e-3;
auto const          GeV                             =   1;
auto const          TeV                             =   1e+3;
#endif
