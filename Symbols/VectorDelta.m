(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



VectorDelta::usage =
	"VectorDelta[vector_] is a Dirac-delta function, \[Integral]\!\(\*
StyleBox[\"dv\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"f\",\nFontSize->14]\)\!\(\*
StyleBox[\"[\",\nFontSize->14]\)\!\(\*
StyleBox[\"v\",\nFontSize->14]\)\!\(\*
StyleBox[\"]\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"VectorDelta\",\nFontSize->14]\)\!\(\*
StyleBox[\"[\",\nFontSize->14]\)\!\(\*
StyleBox[\"v\",\nFontSize->14]\)\!\(\*
StyleBox[\"]\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"=\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"f\",\nFontSize->14]\)\!\(\*
StyleBox[\"[\",\nFontSize->14]\)\!\(\*
StyleBox[\"0\",\nFontSize->14]\)\!\(\*
StyleBox[\"]\",\nFontSize->14]\)."	

Begin["`Private`VectorDelta`"]	
(*Options[VectorDelta] = {}*)



VectorDelta/:MakeBoxes[VectorDelta[vector_],TraditionalForm]:=RowBox[{"\[Delta]","(",ToBoxes[vector,TraditionalForm],")"}]


End[]
