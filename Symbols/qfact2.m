(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



qfact2::usage =
	"qfact2[] is just a symbol for combine constant factors of Gamma function symbols in qfact1"
Begin["`Private`qfact2`"]


qfact2/: MakeBoxes[qfact2[expr_],TraditionalForm]:=
	RowBox[{"(",ToBoxes[expr,TraditionalForm],")"}];
	
qfact2[xx_Integer]:=xx
qfact2[xx_Rational]:=xx


End[]
