(* ::Package:: *)

(* Wolfram Language Package *)


(*
	This software is covered by the GNU General Public License 3.
	Copyright (C) 2021 Shuang-Hong Li
*)


ShowasTable::usage="Whether show the result as a table for some functions."


Parallelized::usage="Whether parallelize the evaluation."


Renormalization::usage="Whether renormalize the current automatically."


EpsOrder::usage="To which order of \[Epsilon], EpsOrder->0 by default."


HoldFlavor::usage="Whether evaluate FlavorDelta, since there are many hardon currents merry differ by flavors, it's convient to hold the flavor when evaluate the OPE and release FlavorDelta later after replace the flavors."


NLO::usage="Whether to evaluate the perturbative contribution to next leading order."


Condensates::usage="Set the condensates which should be calculated."


AutoNDR::usage="Whether automatically choose \!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\)-scheme. If true, for even number of \!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\), DiracGammaScheme will be set to \"BMHV\"; for old number of \!\(\*SuperscriptBox[\(\[Gamma]\), \(5\)]\), DiracGammaScheme will be set to \"NDR\"."


Pole::usage="When evaluating leading perturbative contribution, this Option sets the \!\(\*FractionBox[\(1\), \(\[Epsilon]\)]\)-pole of renormalized current."


SimplifyGamma::usage="Whether simplify the Gamma function appear in qfact1."


(*Begin["`Private`"];*)


(*ShowasTable       =True;
Parallelized      =False;
Renormalization   =True;
EpsOrder          =0;
HoldFlavor        =True;
NLO               ="NLOonly";*)
(*Condensates       ={};*)



(*End[];*)
