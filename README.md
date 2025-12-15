# Description

This is repository is the copy of ```https://github.com/QSSRHelper/QSSRHelper```, which is unaccessible due to the account is flagged for unknown reason

## Install&Usage

This package require the lateset version of FeynCalc.

Install this package by download the whole repository into default directory of Mathematica package.

e.g. download the repository as a zip file and unpackage it into ``` ...\Mathematica\Applications\ ```, which should looks like
```
.../Mathematica
    |
    └--
    └--Applications
    |   |
    |   └-- FeynCalc
    |   |   |...
    |   |   
    |   └-- QSSRHelper
    |...    |
            └-- Auxiliary
            └-- Documentation
            |...
```

Load ```FeynCalc``` before load this package:

```
$LoadFeynArts = True;
Global`$LoadAddOns = {"TARCER", "FeynHelpers"};
(*DistributeDefinitions[$LoadFeynArts,Global`$LoadAddOns];*)

<< FeynCalc`
$FAVerbose = 0;
<< QSSRHelper`
```

The simplest example, run below code in Mathematica Notebook, two outputs are identical  

```
(*propagator in momentum space*)QuarkP[p_] = I   GSD[p]  FAD[p];
GluonP[p_, \[Mu]_, \[Nu]_] = -I   MTD[\[Mu], \[Nu]]  FAD[p];

(*propagator in coordinate space*)
QuarkX[x_] = FourierPX[QuarkP[p], {p, x}];
GluonX[x_, \[Mu]_, \[Nu]_] = 
  FourierPX[GluonP[p, \[Mu], \[Nu]], {p, x}];
(*trace over dirac and su(n) indices*)
TTR[expr_] := SUNSimplify[TR[SUNTrace[expr, SUNTraceEvaluate -> True]]]

(*massless sunrise diagram for hybrid currents J^\[Mu]=g u^bar \
G^\[Mu]\[Alpha] \[Gamma]_\[Alpha] u *)

(*momentum integral for <J^\[Mu]J^\[Nu]>*)
dia = g^2   TTR[
    GAD[\[Alpha]] . SUNT[n] . QuarkP[p + l + k] . GAD[\[Beta]] . 
      SUNT[n] . 
      QuarkP[k]  (-I   FVD[-l, \[Mu]]  MTD[\[Alpha], \[Eta]] + 
       I   FVD[-l, \[Alpha]]  MTD[\[Mu], \[Eta]])  GluonP[
      l, \[Eta], \[Rho]]  (I   FVD[l, \[Nu]]  MTD[\[Beta], \[Rho]] - 
       I   FVD[l, \[Beta]]  MTD[\[Nu], \[Rho]])];
dia = IntegrateP[dia, {k, l}];
QEvaluate[I   dia   ScaleMu^(2  (4 - D)), p]


(*Fourier Transformation for<J^\[Mu](x)J^\[Nu](0)>*)
(* the propagator connecting two G^\[Mu]\[Nu] *)
GGX = FourierPX[(-I   FVD[-l, \[Mu]]  MTD[\[Alpha], \[Eta]] + 
      I   FVD[-l, \[Alpha]]  MTD[\[Mu], \[Eta]])  GluonP[
     l, \[Eta], \[Rho]]  (I   FVD[l, \[Nu]]  MTD[\[Beta], \[Rho]] - 
      I   FVD[l, \[Beta]]  MTD[\[Nu], \[Rho]]), {l, x}];

dia = g^2   TTR[
    GAD[\[Alpha]] . SUNT[n] . QuarkX[x] . GAD[\[Beta]] . SUNT[n] . 
      QuarkX[-x]   GGX];
dia = FourierXP[dia, {x, p}];
QEvaluate[I   dia   ScaleMu^(2  (4 - D)), p]
```

The Documentation is out of date and may useless at this time, the development of this package is far from complete due to the lack of time.

## to do:

- enhance the QSimplify so that it can recongize the color indicex in the Tetraquark current
- SU(3) flavor decomposition
- generate diagrams

