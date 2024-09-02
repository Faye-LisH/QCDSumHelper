(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)






QSimplify::usage= "QSimplify[expr_] is a specific Simplify function that simplify the dummy indices and Gamma functions";


Begin["`Private`QSimplify`"]


Options[QSimplify] = {
	NonCommutative->{},
	Contract->False,
	Lorentz->"Auto",
	Color->"Auto",
	Symmetry->"None",
	Separate->{"null"},
	SimplifyFlavor->False,
	SimplifyGamma->True}
(*---------------------------------------------------------------*)

QSimplify[expr_List,opts___]:=QSimplify[#,opts]&/@expr


QSimplify[expr_,OptionsPattern[]]:=Block[{tmp,flaor=OptionValue[SimplifyFlavor]},




tmp=IndexSimplify[
	If[flavor===True,
		expr/.FlavorDelta[i_,j_]:>If[i===j,1,0]
	,
		expr
	],
	QSimplify2->True,
	NonCommutative->OptionValue[NonCommutative],
	Contract->OptionValue[Contract],
	Lorentz->OptionValue[Lorentz],
	Color->OptionValue[Color],
	Symmetry->OptionValue[Symmetry],
	QSimplify2Rules->OptionValue[Separate],
	SimplifyGamma->OptionValue[SimplifyGamma]]/.qfact2[aa_]:>qfact2[aa//Simplify];
	
If[FreeQ[expr,qfact1]&&FreeQ[expr,qfact2],
	tmp/.{qfact1->Identity,qfact2->Identity}
,
	tmp
]
	
]
	


End[]
