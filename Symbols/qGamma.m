(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



qGamma::usage = 
"qGamma[x] is just qGamma symble and doesn't do any evaluation";
Begin["`Private`qGamma`"]


qGamma /: MakeBoxes[qGamma[x_],TraditionalForm]:=
	RowBox[{"\[CapitalGamma]","(",ToBoxes[ExpandAll[x],TraditionalForm],")"}]
	
	(*qGamma[1]=1;
qGamma[2]=1;
qGamma[3]=2;
*)
(* keep all \[CapitalGamma] unevaluated to avoid some error/bug be hidden in function like FourierXP, Qdx : e.g. \[CapitalGamma][2]\[CapitalGamma][0]-\[CapitalGamma][1]\[CapitalGamma][0] = 0  
  although the equality hold, but the \[CapitalGamma][0] indicate there have something wrong in evaluation/input *)



End[]
