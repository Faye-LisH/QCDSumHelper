(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



qfactor2::usage =
	"qfactor2[] is just a symbol for combine constant factors of Gamma function symbols in qfact1"

Begin["`Private`qfactor2`"]	


qfactor2/: MakeBoxes[qfactor2[expr_],TraditionalForm]:=
	RowBox[{"(",ToBoxes[expr,TraditionalForm],")"}];
	
qfactor2[xx_Integer]:=xx
qfactor2[xx_Rational]:=xx



End[]
