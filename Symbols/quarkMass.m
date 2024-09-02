(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



quarkMass::usage = "quarkMass[a_] give the symbol of quark mass with flavor a."
Begin["`Private`quarkMass`"]


quarkMass/: MakeBoxes[quarkMass[i_],TraditionalForm]:=SubscriptBox["m", ToBoxes[i]]
quarkMass/: MakeBoxes[quarkMass[],TraditionalForm]:="m"



End[]
