//
//  Part of the AliAnalysisUtility package
//
//  Fit Functions
//
//  Author              Nicola Rubini
//  Created             22/11/2021
//  Last modified       22/11/2021
#ifndef ALIANALYSISUTILITY_FUNCTIONS_H
#define ALIANALYSISUTILITY_FUNCTIONS_H
//
//>>    Functions File
#include "AliAnalysisUtility.h"
//
//  --  --  Global Variables  --  --  //
//
const Double_t  integralPrecision   =   1.e-12;
//
//  --  --  FIT Custom Functions  --  --  //
//
//
//>>-->>    PT-Dependent
//
//_______________________________________________________________>> LevyTsallis
//_____________________________________________________________________________
//
Double_t                _LevyTsallis                ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fPT     = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fEnne   = fParams[1];
    Double_t    fSlop   = fParams[2];
    Double_t    fdNdY   = fParams[3];
    
    Double_t    fNum1   = (fEnne-1)*(fEnne-2);
    Double_t    fDen1   = fEnne*fSlop*(fEnne*fSlop+fMass*(fEnne-2));
    Double_t    fFac1   = fNum1/fDen1;
    
    Double_t    fMasT   = sqrt(fMass*fMass+fPT*fPT);
    Double_t    fNum2   = fMasT - fMass;
    Double_t    fDen2   = fEnne*fSlop;
    Double_t    fFac2   = TMath::Power((1 + fNum2/fDen2),(-fEnne));
    
    return      fPT*fdNdY*fFac1*fFac2;
}
TF1                    *fLevyTsallis                = new TF1 ("LevyTsallis",_LevyTsallis,0.,100.,4);
void                    fSetLevyTsallis             ()  {
    fLevyTsallis        ->  SetParNames("Mass","n","T","dN_dy");
}
//
//_____________________________________________________________>> PTExponential
//_____________________________________________________________________________
//
Double_t                _PTExponential              ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fPT     = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fSlop   = fParams[0];
    Double_t    fdNdY   = fParams[1];
    
    return      fPT*(fdNdY/(fSlop*fSlop))*TMath::Exp(-fPT / fSlop);
}
TF1                    *fPTExponential              = new TF1 ("PTExponential",_PTExponential,0.,100.,3);
void                    fSetPTExponential           ()  {
    fPTExponential      ->  SetParNames("Mass","T","dN_dy");
}
//
//_____________________________________________________________>> MTExponential
//_____________________________________________________________________________
//
Double_t                _MTExponential              ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fPT     = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fSlop   = fParams[1];
    Double_t    fdNdY   = fParams[2];
    
    Double_t    fMassT  = TMath::Sqrt( fPT*fPT + fMass*fMass );
    
    return      fdNdY*fPT*TMath::Exp( (fMass-fMassT) / fSlop )/( fSlop* ( fSlop + fMass ) );
}
TF1                    *fMTExponential              = new TF1 ("MTExponential",_MTExponential,0.,100.,3);
void                    fSetMTExponential           ()  {
    fMTExponential      ->  SetParNames("Mass","T","dN_dy");
}
//
//________________________________________________________________>> FermiDirac
//_____________________________________________________________________________
//
Double_t                _FermiDirac                 ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fPT     = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fSlop   = fParams[1];
    Double_t    fdNdY   = fParams[2];
    
    Double_t    fMassT  = TMath::Sqrt( fPT*fPT + fMass*fMass );
    
    return      fPT*fdNdY*(1./( TMath::Exp( fMassT / fSlop) +1 ));
}
TF1                    *fFermiDirac                 = new TF1 ("FermiDirac",_FermiDirac,0.,100.,3);
void                    fSetFermiDirac              ()  {
    fFermiDirac         ->  SetParNames("Mass","T","dN_dy");
}
//
//_________________________________________________________________>> Boltzmann
//_____________________________________________________________________________
//
Double_t                _Boltzmann                  ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fPT     = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fSlop   = fParams[1];
    Double_t    fdNdY   = fParams[2];
    
    Double_t    fMassT  = TMath::Sqrt( fPT*fPT + fMass*fMass );
    
    return      fdNdY*fPT*fMassT*TMath::Exp( (fMass-fMassT) / fSlop )/( fSlop* ( 2*fSlop*fSlop + 2*fSlop*fMass + fMass*fMass ) );
}
TF1                    *fBoltzmann                  = new TF1 ("Boltzmann",_Boltzmann,0.,100.,3);
void                    fSetBoltzmann               ()  {
    fBoltzmann          ->  SetParNames("Mass","T","dN_dy");
}
//
//_____________________________________________________________>> Bose-Einstein
//_____________________________________________________________________________
//
Double_t                _BoseEinstein               ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fPT     = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fSlop   = fParams[1];
    Double_t    fdNdY   = fParams[2];
    
    Double_t    fMassT  = TMath::Sqrt( fPT*fPT + fMass*fMass );
    
    return      fPT*fdNdY*(1./( TMath::Exp( fMassT / fSlop) -1 ))*(TMath::Exp( fMass / fSlop ) -1 );
}
TF1                    *fBoseEinstein               = new TF1 ("BoseEinstein",_BoseEinstein,0.,100.,3);
void                    fSetBoseEinstein            ()  {
    fBoseEinstein       ->  SetParNames("Mass","T","dN_dy");
}
//
//_____________________________________________________________>> Bose-Einstein
//_____________________________________________________________________________
//
Double_t                _PowerLaw                   ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fPT     = fVar[0];
    Double_t    fEnne   = fParams[0];
    Double_t    fSlop   = fParams[1];
    Double_t    fdNdY   = fParams[2];
    
    return      fPT*fdNdY*TMath::Power( (1. + fPT/fSlop) , -fEnne);
}
TF1                    *fPowerLaw                   = new TF1 ("PowerLaw",_PowerLaw,0.,100.,3);
void                    fSetPowerLaw                ()  {
    fPowerLaw       ->  SetParNames("n","T","dN_dy");
}
//
Double_t                _fTest                      ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fPT     = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fEnne   = fParams[1];
    Double_t    fSlop   = fParams[2];
    Double_t    fdNdY   = fParams[3];
    
    Double_t    fNum1   = (fEnne-1)*(fEnne-2);
    Double_t    fDen1   = fEnne*fSlop*(fEnne*fSlop+fMass*(fEnne-2));
    Double_t    fFac1   = fNum1/fDen1;
    
    Double_t    fMasT   = sqrt(fMass*fMass+fPT*fPT);
    Double_t    fNum2   = fMasT - fMass;
    Double_t    fDen2   = fEnne*fSlop;
    Double_t    fFac2   = TMath::Power((1 + fNum2/fDen2),(-fEnne));
    
    return      fPT*fdNdY*fFac1*fFac2*fPT*fFac1*fFac2;
}
TF1                    *fLevyTsallis2                = new TF1 ("LevyTsallis2",_fTest,0.,100.,4);
//
Double_t                _fTest2D                    ( Double_t * fVar, Double_t * fParams ) {
    Double_t    *   fVarX   =  new Double_t [1];
    fVarX[0] = fVar[0];
    Double_t    *   fVarY   =  new Double_t [1];
    fVarY[0] = fVar[1];
    Double_t    *   fParX   =  new Double_t [4];
    fParX[0] = fParams[0];
    fParX[1] = fParams[1];
    fParX[2] = fParams[2];
    fParX[3] = TMath::Sqrt(fParams[6]);
    Double_t    *   fParY   =  new Double_t [4];
    fParY[0] = fParams[3];
    fParY[1] = fParams[4];
    fParY[2] = fParams[5];
    fParY[3] = TMath::Sqrt(fParams[6]);
    return      _LevyTsallis(fVarX,fParX)*_LevyTsallis(fVarY,fParY);
}
TF2                    *fLevyTsallis2D               = new TF2 ("LevyTsallis2D",_fTest2D,0.,100.,0.,100.,7,2);
void                    fSetLevyTsallis2D           ()  {
    fLevyTsallis2D      ->  SetParNames("XMass","Xn","XT","YMass","Yn","YT","dN_dy");
}
//
Double_t                _BreitWigner                    ( Double_t * fVar, Double_t * fParams ) {
    return      TMath::BreitWigner(fVar[0],fParams[0],fParams[1]);
}
TF1                    *fBreitWigner                   = new TF1 ("BreitWigner",_BreitWigner,0.,100.,2);
void                    fSetBreitWigner           ()  {
    fBreitWigner      ->  SetParNames("Mean","Sigma");
}
//
Double_t                _BreitWigner2D                    ( Double_t * fVar, Double_t * fParams ) {
    return      TMath::BreitWigner(fVar[0],fParams[0],fParams[1])*TMath::BreitWigner(fVar[1],fParams[0],fParams[1]);
}
TF2                    *fBreitWigner2D                   = new TF2 ("BreitWigner2D",_BreitWigner2D,0.,100.,0.,100.,2,2);
void                    fSetBreitWigner2D           ()  {
    fBreitWigner2D      ->  SetParNames("Mean","Sigma");
}
//
//
std::vector<TF1*>   kAllFunctions   =   { fLevyTsallis, fMTExponential, fBoseEinstein, fBoltzmann, fPowerLaw };
std::vector<std::tuple<TF1*,Float_t,Float_t,TString>>   kStandardSystematicFitFunctions =
    {   {   fLevyTsallis,   0.4,    10.,    "EMRQSI"},
        {   fLevyTsallis,   0.4,    5.0,    "EMRQSI"},
        {   fLevyTsallis,   0.4,    3.0,    "EMRQSI"},
        {   fLevyTsallis,   0.4,    2.0,    "EMRQSI"},
        {   fLevyTsallis,   0.4,    1.5,    "EMRQSI"},
        {   fLevyTsallis,   0.4,    1.0,    "EMRQSI"},
        {   fMTExponential, 0.4,    2.0,    "EMRQSI"},
        {   fMTExponential, 0.4,    1.5,    "EMRQSI"},
        {   fMTExponential, 0.4,    1.2,    "EMRQSI"},
        {   fMTExponential, 0.4,    1.0,    "EMRQSI"},
        {   fBoseEinstein,  0.4,    2.0,    "EMRQSI"},
        {   fBoseEinstein,  0.4,    1.5,    "EMRQSI"},
        {   fBoseEinstein,  0.4,    1.2,    "EMRQSI"},
        {   fBoseEinstein,  0.4,    1.0,    "EMRQSI"},
        {   fBoltzmann,     0.4,    2.0,    "EMRQSI"},
        {   fBoltzmann,     0.4,    1.5,    "EMRQSI"},
        {   fBoltzmann,     0.4,    1.2,    "EMRQSI"},
        {   fBoltzmann,     0.4,    1.0,    "EMRQSI"}
    };


/*
std::vector<std::pair<TF1*,std::vector<float>>> fSystFitFunctions  =   {  {fMTExponential,{1.2,1.4,1.6,2.0}}, {fBoseEinstein,{1.2,1.4,1.6,2.0}}, {fBoltzmann,{1.2,1.4,1.6,2.0}}, {fPowerLaw,{2.0,2.8,4.0}}, {fLevyTsallis,{1.2,1.4,1.6,2.0,2.8,4.0}} };//, {fBGBlastWave,{1.2,1.4,1.6,2.0,2.8,4.0}} };
 */




//________________________________________________>> Boltzmann-Gibbs Blast-Wave
//_____________________________________________________________________________
/*
//
Double_t                _BGBlastWaveIntg            ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fR      = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fPT     = fParams[1];
    Double_t    fMxBeta = fParams[2];
    Double_t    fSlop   = fParams[3];
    Double_t    fEnne   = fParams[4];
    
    Double_t    fMassT  = TMath::Sqrt( fPT*fPT + fMass*fMass );
    
    Double_t    fBeta   =   min( fMxBeta*TMath::Power( fR,fEnne ) , ( 1-1.e-16 ) ); // Just a degeneration protection
    Double_t    fRho    =   TMath::ATanH( fBeta );
    Double_t    fFac1   =   min( ( fPT/fSlop )*TMath::SinH( fRho ) , 700. ); // Just a degeneration protection
    Double_t    fFac2   =   ( fMassT/fSlop )*TMath::CosH( fRho );
    
    return      fR*fMassT*TMath::BesselI0( fFac1 )*TMath::BesselK1( fFac2 );
}
TF1                    *fBGBlastWaveIntg            = new TF1 ("BGBlastWaveIntg",_BGBlastWaveIntg,0.,1.,5);
Double_t                _BGBlastWaveNorm            ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fPT     = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fMxBeta = fParams[1];
    Double_t    fTemp   = fParams[2];
    Double_t    fEnne   = fParams[3];
    
    if ( !fBGBlastWaveIntg )    fBGBlastWaveIntg    = new TF1 ("BGBlastWaveIntg",_BGBlastWaveIntg,0.,1.,5);
    fBGBlastWaveIntg->SetParameters(fMass,fPT,fMxBeta,fTemp,fEnne);
    
    return      fPT*fBGBlastWaveIntg->Integral(0., 1., integralPrecision);
}
TF1                    *fBGBlastWaveNorm            = new TF1 ("BGBlastWaveNorm",_BGBlastWaveIntg,0.,100.,4);
Double_t                _BGBlastWave                ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fPT     = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fMxBeta = fParams[1];
    Double_t    fTemp   = fParams[2];
    Double_t    fEnne   = fParams[3];
    Double_t    fdNdY   = fParams[4];
    
    if ( !fBGBlastWaveNorm )    fBGBlastWaveNorm    = new TF1 ("BGBlastWaveNorm",_BGBlastWaveIntg,0.,100.,4);
    fBGBlastWaveNorm->SetParameters(fMass,fMxBeta,fTemp,fEnne);
    
    Double_t    fIntg   = fBGBlastWaveNorm->Integral(0., 10., integralPrecision) ;
    Double_t    fVal    = fBGBlastWaveNorm->Eval(fPT);
    
    return      fdNdY*fVal/fIntg;
}
TF1                    *fBGBlastWave                = new TF1 ("BGBlastWave",_BGBlastWave,0.,100.,5);
void                    fSetBGBlastWave             ()  {
    fBGBlastWave        ->  SetParNames("Mass","MaxBeta","Temp","n","dN_dY");
}
//
*/
//_____________________________________________________________________________
//________________________________________________>> Asymmetric Gaussian
//_____________________________________________________________________________
//
Double_t                _AsymmGauss                 ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fVarX   = fVar[0];
    Double_t    fNorm   = fParams[0];
    Double_t    fMean   = fParams[1];
    Double_t    fSig1   = fParams[2];
    Double_t    fSig2   = fParams[3];
    
    if ( fVarX <= fMean )   return  fNorm/(TMath::Sqrt(2*TMath::Pi()) + TMath::Sqrt(fSig1) + TMath::Sqrt(fSig2) )*TMath::Gaus(fVarX,fMean,fSig1);
    else                    return  fNorm/(TMath::Sqrt(2*TMath::Pi()) + TMath::Sqrt(fSig1) + TMath::Sqrt(fSig2) )*TMath::Gaus(fVarX,fMean,fSig2);
}
TF1 *                   fAsymmGauss                 = new TF1 ("AsymmGauss",_AsymmGauss,-100.,100.,4);
void                    fSetAsymmGauss                 ()  {
    fAsymmGauss         ->  SetParNames("Norm","Mean","SigmaL","SigmaR");
    
}
//
//_____________________________________________________________________________
//
Double_t                _Gauss                      ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fVarX   = fVar[0];
    Double_t    fNorm   = fParams[0];
    Double_t    fMean   = fParams[1];
    Double_t    fSig_   = fParams[2];
    
    return  fNorm*TMath::Gaus(fVarX,fMean,fSig_,true);
}
TF1 *                   fGauss                      = new TF1 ("Gauss",_Gauss,-100.,100.,3);
void                    fSetGauss                   ()  {
    fGauss              ->  SetParNames("Norm","Mean","Sigma");
}
//
//_____________________________________________________________________________
//

//
//_____________________________________________________________________________


/* To be cleanse */

/*****************************************************************/
/* BOLTZMANN-GIBBS BLAST-WAVE */
/*****************************************************************/

static TF1 *fBGBlastWave_Integrand = NULL;
static TF1 *fBGBlastWave_Integrand_num = NULL;
static TF1 *fBGBlastWave_Integrand_den = NULL;
Double_t
BGBlastWave_Integrand(const Double_t *x, const Double_t *p)
{
  
  /*
     x[0] -> r (radius)
     p[0] -> mT (transverse mass)
     p[1] -> pT (transverse momentum)
     p[2] -> beta_max (surface velocity)
     p[3] -> T (freezout temperature)
     p[4] -> n (velocity profile)
  */
  
  Double_t r = x[0];
  Double_t mt = p[0];
  Double_t pt = p[1];
  Double_t beta_max = p[2];
  Double_t temp_1 = 1. / p[3];
  Double_t n = p[4];

  Double_t beta = beta_max * TMath::Power(r, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;
  Double_t rho = TMath::ATanH(beta);
  Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
  if (argI0 > 700.) argI0 = 700.;
  Double_t argK1 = mt * TMath::CosH(rho) * temp_1;
  //  if (argI0 > 100 || argI0 < -100)
  //    printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", r, pt, beta_max, 1. / temp_1, n, mt, beta, rho, argI0, argK1);
  return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);
  
}

Double_t
BGBlastWave_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  Double_t norm = p[4];
  
  if (!fBGBlastWave_Integrand)
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
  Double_t integral = fBGBlastWave_Integrand->Integral(0., 1., 1.e-6);
  return norm * pt * integral;
}

Double_t
BGBlastWaveRatio_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max_num = p[1];
  Double_t temp_num = p[2];
  Double_t n_num = p[3];
  Double_t norm_num = p[4];
  Double_t beta_max_den = p[5];
  Double_t temp_den = p[6];
  Double_t n_den = p[7];
  Double_t norm_den = p[8];
  
  if (!fBGBlastWave_Integrand_num)
    fBGBlastWave_Integrand_num = new TF1("fBGBlastWave_Integrand_num", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand_num->SetParameters(mt, pt, beta_max_num, temp_num, n_num);
  Double_t integral_num = fBGBlastWave_Integrand_num->Integral(0., 1.);

  if (!fBGBlastWave_Integrand_den)
    fBGBlastWave_Integrand_den = new TF1("fBGBlastWave_Integrand_den", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand_den->SetParameters(mt, pt, beta_max_den, temp_den, n_den);
  Double_t integral_den = fBGBlastWave_Integrand_den->Integral(0., 1.);

  return (norm_num / norm_den) * (integral_num / integral_den);
}

Double_t
BGBlastWaveParticleRatio_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass_num = p[0];
  Double_t mass_den = p[1];
  Double_t mt_num = TMath::Sqrt(pt * pt + mass_num * mass_num);
  Double_t mt_den = TMath::Sqrt(pt * pt + mass_den * mass_den);
  Double_t beta_max = p[2];
  Double_t temp = p[3];
  Double_t n = p[4];
  Double_t norm_num = p[5];
  Double_t norm_den = p[6];
  
  if (!fBGBlastWave_Integrand_num)
    fBGBlastWave_Integrand_num = new TF1("fBGBlastWave_Integrand_num", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand_num->SetParameters(mt_num, pt, beta_max, temp, n);
  Double_t integral_num = fBGBlastWave_Integrand_num->Integral(0., 1.);
  
  if (!fBGBlastWave_Integrand_den)
    fBGBlastWave_Integrand_den = new TF1("fBGBlastWave_Integrand_den", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand_den->SetParameters(mt_den, pt, beta_max, temp, n);
  Double_t integral_den = fBGBlastWave_Integrand_den->Integral(0., 1.);

  return (norm_num / norm_den) * (integral_num / integral_den);
}

Double_t
BGBlastWave_Func_OneOverPt(const Double_t *x, const Double_t *p)
{
  /* 1/pt dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  Double_t norm = p[4];
  
  if (!fBGBlastWave_Integrand)
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
  Double_t integral = fBGBlastWave_Integrand->Integral(0., 1., 1.e-3);

  return norm * integral;
}

TF1 *
BGBlastWave(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{
  
  TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func, 0., 10., 5);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("Mass", "beta_max", "T", "n", "dN_dy");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0.01, 0.99);
  fBGBlastWave->SetParLimits(2, 0.01, 1.);
  fBGBlastWave->SetParLimits(3, 0.01, 50.);
  return fBGBlastWave;
}

TF1 *
BGBlastWaveRatio(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{
  
  TF1 *fBGBlastWave = new TF1(name, BGBlastWaveRatio_Func, 0., 10., 9);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("mass", "beta_max_num", "T_num", "n_num", "norm_num", "beta_max_den", "T_den", "n_den", "norm_den");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0.01, 0.99);
  fBGBlastWave->SetParLimits(2, 0.01, 1.);
  fBGBlastWave->SetParLimits(3, 0.01, 10.);
  fBGBlastWave->SetParLimits(5, 0.01, 0.99);
  fBGBlastWave->SetParLimits(6, 0.01, 1.);
  fBGBlastWave->SetParLimits(7, 0.01, 10.);
  return fBGBlastWave;
}

TF1 *
BGBlastWaveParticleRatio(const Char_t *name, Double_t mass_num, Double_t mass_den, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm_num = 1.e6, Double_t norm_den = 1.e6)
{
  
  TF1 *fBGBlastWave = new TF1(name, BGBlastWaveParticleRatio_Func, 0., 10., 7);
  fBGBlastWave->SetParameters(mass_num, mass_den, beta_max, temp, n, norm_num, norm_den);
  fBGBlastWave->SetParNames("mass_num", "mass_den", "beta_max", "T", "n", "norm_num", "norm_den");
  fBGBlastWave->FixParameter(0, mass_num);
  fBGBlastWave->FixParameter(1, mass_den);
  fBGBlastWave->SetParLimits(2, 0.01, 0.99);
  fBGBlastWave->SetParLimits(3, 0.01, 1.);
  fBGBlastWave->SetParLimits(4, 0.01, 10.);
  return fBGBlastWave;
}

TF1 *BGBlastWave_OneOverPT(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{
  
  TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func_OneOverPt, 0., 10., 5);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0.01, 0.99);
  fBGBlastWave->SetParLimits(2, 0.01, 1.);
  fBGBlastWave->SetParLimits(3, 0.01, 50.);
  return fBGBlastWave;
}

/*****************************************************************/
/* TSALLIS BLAST-WAVE */
/*****************************************************************/

static TF1 *fTsallisBlastWave_Integrand_r = NULL;
Double_t
TsallisBlastWave_Integrand_r(const Double_t *x, const Double_t *p)
{
  /*
     x[0] -> r (radius)
     p[0] -> mT (transverse mass)
     p[1] -> pT (transverse momentum)
     p[2] -> beta_max (surface velocity)
     p[3] -> T (freezout temperature)
     p[4] -> n (velocity profile)
     p[5] -> q
     p[6] -> y (rapidity)
     p[7] -> phi (azimuthal angle)
  */
  
  Double_t r = x[0];
  Double_t mt = p[0];
  Double_t pt = p[1];
  Double_t beta_max = p[2];
  Double_t temp_1 = 1. / p[3];
  Double_t n = p[4];
  Double_t q = p[5];
  Double_t y = p[6];
  Double_t phi = p[7];

  if (q <= 1.) return r;

  Double_t beta = beta_max * TMath::Power(r, n);
  Double_t rho = TMath::ATanH(beta);
  
  Double_t part1 = mt * TMath::CosH(y) * TMath::CosH(rho);
  Double_t part2 = pt * TMath::SinH(rho) * TMath::Cos(phi);
  Double_t part3 = part1 - part2;
  Double_t part4 = 1 + (q - 1.) * temp_1 * part3;
  Double_t expo = -1. / (q - 1.);
  //  printf("part1=%f, part2=%f, part3=%f, part4=%f, expo=%f\n", part1, part2, part3, part4, expo);
  Double_t part5 = TMath::Power(part4, expo);

  return r * part5;
}

static TF1 *fTsallisBlastWave_Integrand_phi = NULL;
Double_t
TsallisBlastWave_Integrand_phi(const Double_t *x, const Double_t *p)
{
  /*
     x[0] -> phi (azimuthal angle)
  */
  
  Double_t phi = x[0];
  fTsallisBlastWave_Integrand_r->SetParameter(7, phi);
  Double_t integral = fTsallisBlastWave_Integrand_r->Integral(0., 1., integralPrecision);
  return integral;
}

static TF1 *fTsallisBlastWave_Integrand_y = NULL;
Double_t
TsallisBlastWave_Integrand_y(const Double_t *x, const Double_t *p)
{
  /*
     x[0] -> y (rapidity)
  */

  Double_t y = x[0];
  fTsallisBlastWave_Integrand_r->SetParameter(6, y);
  Double_t integral = fTsallisBlastWave_Integrand_phi->Integral(-TMath::Pi(), TMath::Pi(), integralPrecision);
  return TMath::CosH(y) * integral;
}

Double_t
TsallisBlastWave_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  Double_t q = p[4];
  Double_t norm = p[5];

  if (!fTsallisBlastWave_Integrand_r)
    fTsallisBlastWave_Integrand_r = new TF1("fTsallisBlastWave_Integrand_r", TsallisBlastWave_Integrand_r, 0., 1., 8);
  if (!fTsallisBlastWave_Integrand_phi)
    fTsallisBlastWave_Integrand_phi = new TF1("fTsallisBlastWave_Integrand_phi", TsallisBlastWave_Integrand_phi, -TMath::Pi(), TMath::Pi(), 0);
  if (!fTsallisBlastWave_Integrand_y)
    fTsallisBlastWave_Integrand_y = new TF1("fTsallisBlastWave_Integrand_y", TsallisBlastWave_Integrand_y, -0.5, 0.5, 0);

  fTsallisBlastWave_Integrand_r->SetParameters(mt, pt, beta_max, temp, n, q, 0., 0.);
  Double_t integral = fTsallisBlastWave_Integrand_y->Integral(-0.5, 0.5, integralPrecision);
  return norm * pt * integral;
}

TF1 *
TsallisBlastWave(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t q = 1.1, Double_t norm = 1.e6)
{
  
  TF1 *fTsallisBlastWave = new TF1(name, TsallisBlastWave_Func, 0., 10., 6);
  fTsallisBlastWave->SetParameters(mass, beta_max, temp, n, q, norm);
  fTsallisBlastWave->SetParNames("mass", "beta_max", "T", "n", "q", "norm");
  fTsallisBlastWave->FixParameter(0, mass);
  fTsallisBlastWave->SetParLimits(1, 0.01, 0.99);
  fTsallisBlastWave->SetParLimits(2, 0.01, 1.);
  fTsallisBlastWave->SetParLimits(3, 0.1, 10.);
  fTsallisBlastWave->SetParLimits(4, 1.001, 1.2);
  return fTsallisBlastWave;
}


//_____________________________________________________________________________
//
Double_t                fFlat2D                     ( Double_t * fVar, Double_t * fParams ) {
    Double_t    fX      = fVar[0];
    Double_t    fY      = fVar[1];
    
    Double_t    fFlatPar= fParams[0];
    
    return      fFlatPar;
}
TF1 *                   fFitFlat2D                  = new TF2 ("fFitFlat2D",fFlat2D,0.,100.,0.,100.,1);
//
//_____________________________________________________________________________

void                    fSetAllFunctions            ()  {
    fSetLevyTsallis();
    fSetPTExponential();
    fSetMTExponential();
    fSetFermiDirac();
    fSetBoltzmann();
    //fSetBGBlastWave();
    fSetBoseEinstein();
    fSetPowerLaw();
    fSetAsymmGauss();
    fSetGauss();
    fSetLevyTsallis2D();
    fSetBreitWigner();
    fSetBreitWigner2D();
}
TF1                    *fBGBlastWave                = BGBlastWave("BGBlastWave",1.019455);

#endif
