(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



GhostPropagatorX::usage="Symbol of Free Ghost propagator in coordinate space"


Begin["`Private`GhostPropagatorX`"]	


GhostPropagatorX/:MakeBoxes[GhostPropagatorX[x:Except[_List]],TraditionalForm]:=RowBox[{"\[EmptyCircle]","(",ToBoxes[x],")"}]
GhostPropagatorX/:MakeBoxes[GhostPropagatorX[x:Except[_List],{a_,b_}],TraditionalForm]:=SuperscriptBox[RowBox[{"\[EmptyCircle]","(",ToBoxes[x],")"}],RowBox[{ToBoxes[a],ToBoxes[b]}]]


GhostPropagatorX/:MakeBoxes[GhostPropagatorX[{x_,y_}],TraditionalForm]:=RowBox[{"\[EmptyCircle]","(",ToBoxes[x],",",ToBoxes[y],")"}]
GhostPropagatorX/:MakeBoxes[GhostPropagatorX[{x_,y_},{a_,b_}],TraditionalForm]:=SuperscriptBox[RowBox[{"\[EmptyCircle]","(",ToBoxes[x],",",ToBoxes[y],")"}],RowBox[{ToBoxes[a],ToBoxes[b]}]]


End[]
