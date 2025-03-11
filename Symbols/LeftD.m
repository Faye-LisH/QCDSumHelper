(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



LeftD::usage="Left partial derivative symbol"


Begin["`Private`LeftD`"]	


LeftD/:MakeBoxes[LeftD[pv_],TraditionalForm]:=Block[{tmp},
tmp=FCI[pv];
SubsuperscriptBox[OverscriptBox["\[PartialD]","\[LeftArrow]"],ToBoxes[tmp[[2,1]],TraditionalForm],ToBoxes[tmp[[1,1]],TraditionalForm]]
]


End[]
