(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



Quark::usage = "Quark[f_] gives the quark symbol with flavor f."
Begin["`Private`Quark`"]


Quark/: MakeBoxes[Quark[i_],TraditionalForm]:=ToBoxes[i]
Quark[f_,1]:=Quark[f]
(*Quark[f_]:=Quark[ToString[f]]/;Head[f]=!=String
Quark[f_,-1]:=Quark[ToString[f],-1]/;Head[f]=!=String*)
Quark[{f_,-1}]:=Quark[f,-1]
Quark/: MakeBoxes[Quark[i_,-1],TraditionalForm]:=OverscriptBox[ToBoxes[i],"_"]
Quark[Quark[qs__]]:=Quark[qs]
Quark/:MakeBoxes[Quark[qs__Quark],TraditionalForm]:=RowBox[{qs}/.{Quark[q_]:>ToBoxes[q], Quark[q_,-1]:>OverscriptBox[ToBoxes[q],"_"]}]


End[]
