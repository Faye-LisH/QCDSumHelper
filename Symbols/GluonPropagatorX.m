(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



GluonPropagatorX::usage="Symbol of Free Gluon propagator in coordinate space"


Begin["`Private`GluonPropagatorX`"]	


GluonPropagatorX/:MakeBoxes[GluonPropagatorX[x:Except[_List]],TraditionalForm]:=RowBox[{"\[CapitalDelta]","(",ToBoxes[x],")"}]
GluonPropagatorX/:MakeBoxes[GluonPropagatorX[x:Except[_List],{u_,v_},{a_,b_}],TraditionalForm]:=SubsuperscriptBox[RowBox[{"\[CapitalDelta]","(",ToBoxes[x],")"}],RowBox[{ToBoxes[u],ToBoxes[v]}],RowBox[{ToBoxes[a],ToBoxes[b]}]]


GluonPropagatorX/:MakeBoxes[GluonPropagatorX[{x_,y_}],TraditionalForm]:=RowBox[{"\[CapitalDelta]","(",ToBoxes[x],",",ToBoxes[y],")"}]
GluonPropagatorX/:MakeBoxes[GluonPropagatorX[{x_,y_},{u_,v_},{a_,b_}],TraditionalForm]:=SubsuperscriptBox[RowBox[{"\[CapitalDelta]","(",ToBoxes[x],",",ToBoxes[y],")"}],RowBox[{ToBoxes[u],ToBoxes[v]}],RowBox[{ToBoxes[a],ToBoxes[b]}]]


End[]
