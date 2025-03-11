(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



FermionPropagatorX::usage="Symbol of Free Fermion propagator in coordinate space"


Begin["`Private`FermionPropagatorX`"]	


FermionPropagatorX/:MakeBoxes[FermionPropagatorX[x:Except[_List]],TraditionalForm]:=RowBox[{"S","(",ToBoxes[x],")"}]
FermionPropagatorX/:MakeBoxes[FermionPropagatorX[x:Except[_List],f_],TraditionalForm]:=RowBox[{SubscriptBox["S",ToBoxes[f]],"(",ToBoxes[x],")"}]

FermionPropagatorX/:MakeBoxes[FermionPropagatorX[x:Except[_List],f_,{i_,j_},{a_,b_}],TraditionalForm]:=SubsuperscriptBox[RowBox[{SubscriptBox["S",ToBoxes[f]],"(",ToBoxes[x],")"}],RowBox[{ToBoxes[i],ToBoxes[j]}],RowBox[{ToBoxes[a],ToBoxes[b]}]]


FermionPropagatorX/:MakeBoxes[FermionPropagatorX[{x_,y_}],TraditionalForm]:=RowBox[{"S","(",ToBoxes[x],",",ToBoxes[y],")"}]
FermionPropagatorX/:MakeBoxes[FermionPropagatorX[{x_,y_},f_],TraditionalForm]:=RowBox[{SubscriptBox["S",ToBoxes[f]],"(",ToBoxes[x],",",ToBoxes[y],")"}]

FermionPropagatorX/:MakeBoxes[FermionPropagatorX[{x_,y_},f_,{i_,j_},{a_,b_}],TraditionalForm]:=SubsuperscriptBox[RowBox[{SubscriptBox["S",ToBoxes[f]],"(",ToBoxes[x],",",ToBoxes[y],")"}],RowBox[{ToBoxes[i],ToBoxes[j]}],RowBox[{ToBoxes[a],ToBoxes[b]}]]


End[]
