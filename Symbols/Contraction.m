(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



Contraction::usage="Symbol to denote the Wick-Contraction for two fields."


Begin["`Private`Contraction`"]	


Contraction/:MakeBoxes[Contraction[f1:Except[_List],f2_,{x1_,x2_},labels___],TraditionalForm]:=Block[{tmp},
	OverscriptBox[RowBox[{ToBoxes[f1,TraditionalForm],"(",ToBoxes[x1,TraditionalForm],")",ToBoxes[f2,TraditionalForm],"(",ToBoxes[x2,TraditionalForm],")"}],"\[OverBracket]"]
]
Contraction/:MakeBoxes[Contraction[{f1_,f2_},k_,{x1_,x2_},labels___],TraditionalForm]:=Block[{tmp},
	OverscriptBox[RowBox[{ToBoxes[f1,TraditionalForm],"(",ToBoxes[x1,TraditionalForm],")",ToBoxes[f2,TraditionalForm],"(",ToBoxes[x2,TraditionalForm],")"}],"\[OverBracket]"]
]


End[]
