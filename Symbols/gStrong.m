(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



gStrong::usage = "Symbol of QCD copuling constant g."
Begin["`Private`gStrong`"]


gStrong/: MakeBoxes[gStrong,TraditionalForm]:=SubscriptBox["g", "s"]


End[]
