(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



qfact1::usage =
	"qfact1[] is just a symbol for combine constant factors involved in Feynman integral and Fourier transformation"
	
Begin["`Private`qfact1`"]


qfact1 /: MakeBoxes[qfact1[expr_],TraditionalForm]:=
	RowBox[{"(",ToBoxes[expr,TraditionalForm],")"}];
	
qfact1[xx_Integer]:=xx
qfact1[xx_Rational]:=xx


End[]
