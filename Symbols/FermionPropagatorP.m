(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



FermionPropagatorP::usage="Symbol of Free Fermion propagator in momentum space"


Begin["`Private`FermionPropagatorP`"]	


FermionPropagatorP/:MakeBoxes[FermionPropagatorP[p_],TraditionalForm]:=RowBox[{"\!\(\*OverscriptBox[\(S\), \(~\)]\)","(",ToBoxes[p],")"}]
FermionPropagatorP/:MakeBoxes[FermionPropagatorP[p_,f_],TraditionalForm]:=RowBox[{SubscriptBox["\!\(\*OverscriptBox[\(S\), \(~\)]\)",ToBoxes[f]],"(",ToBoxes[p],")"}]

FermionPropagatorP/:MakeBoxes[FermionPropagatorP[p_,f_,{i_,j_},{a_,b_}],TraditionalForm]:=SubsuperscriptBox[RowBox[{SubscriptBox["\!\(\*OverscriptBox[\(S\), \(~\)]\)",ToBoxes[f]],"(",ToBoxes[p],")"}],RowBox[{ToBoxes[i],ToBoxes[j]}],RowBox[{ToBoxes[a],ToBoxes[b]}]]


End[]
