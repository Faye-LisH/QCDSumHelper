(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



qfactor1::usage =
	"qfactor1[] is just a symbol for combine constant factors involved in Feynman integral and Fourier transformation"

Begin["`Private`qfactor1`"]	


qfactor1 /: MakeBoxes[qfactor1[expr_],TraditionalForm]:=
	RowBox[{"(",ToBoxes[expr,TraditionalForm],")"}];
	
qfactor1[xx_Integer]:=xx
qfactor1[xx_Rational]:=xx



End[]
