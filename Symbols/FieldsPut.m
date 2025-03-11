(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



FieldsPut::usage="FieldsPut[O_,x_] show an operator O as O(x), means put it at x"
Begin["`Private`FieldsPut`"]	


FieldsPut/:MakeBoxes[FieldsPut[f_,x_],TraditionalForm]:=RowBox[{ToBoxes[f,TraditionalForm],"(",ToBoxes[x],")"}]


End[]
