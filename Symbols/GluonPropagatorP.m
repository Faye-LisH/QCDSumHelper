(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



GluonPropagatorP::usage="Symbol of Free Gluon propagator in momentum space"


Begin["`Private`GluonPropagatorP`"]	


GluonPropagatorP/:MakeBoxes[GluonPropagatorP[p_],TraditionalForm]:=RowBox[{"\!\(\*OverscriptBox[\(\[CapitalDelta]\), \(~\)]\)","(",ToBoxes[p],")"}]
GluonPropagatorP/:MakeBoxes[GluonPropagatorP[p_,{u_,v_},{a_,b_}],TraditionalForm]:=SubsuperscriptBox[RowBox[{"\!\(\*OverscriptBox[\(\[CapitalDelta]\), \(~\)]\)","(",ToBoxes[p],")"}],RowBox[{ToBoxes[u],ToBoxes[v]}],RowBox[{ToBoxes[a],ToBoxes[b]}]]


End[]
