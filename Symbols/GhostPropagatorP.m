(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



GhostPropagatorP::usage="Symbol of Free Ghost propagator in momentum space"


Begin["`Private`GhostPropagatorP`"]	


GhostPropagatorP/:MakeBoxes[GhostPropagatorP[p_],TraditionalForm]:=RowBox[{"\!\(\*OverscriptBox[\(\[EmptyCircle]\), \(~\)]\)","(",ToBoxes[p],")"}]
GhostPropagatorP/:MakeBoxes[GhostPropagatorP[p_,{a_,b_}],TraditionalForm]:=SuperscriptBox[RowBox[{"\!\(\*OverscriptBox[\(\[EmptyCircle]\), \(~\)]\)","(",ToBoxes[p],")"}],RowBox[{ToBoxes[a],ToBoxes[b]}]]


End[]
