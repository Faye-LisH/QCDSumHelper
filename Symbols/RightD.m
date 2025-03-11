(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



RightD::usage="Partial derivative symbol"


Begin["`Private`RightD`"]	


RightD/:MakeBoxes[RightD[pv_],TraditionalForm]:=Block[{tmp},
tmp=FCI[pv];
SubsuperscriptBox["\[PartialD]",ToBoxes[tmp[[2,1]],TraditionalForm],ToBoxes[tmp[[1,1]],TraditionalForm]]
]


End[]
