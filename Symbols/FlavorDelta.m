(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



FlavorDelta::usage ="FlavorDelta[i_,j_] is a Kroncker-delta for flavor i, j."
	
Begin["`Private`FlavorDelta`"]

Options[FlavorDelta] = {HoldFlavor->True}


(*SetAttributes[{FlavorDelta}, Orderless];*)

(*FlavorDelta/:Times[FlavorDelta[a_,b_,ops:OptionsPattern[]],FlavorDelta[b_,c_,ops:OptionsPattern[]]]:=FlavorDelta[a,c,ops]/;(OptionValue[HoldFlavor]===False)
FlavorDelta/:Times[FlavorDelta[a_,b_,ops:OptionsPattern[]],FlavorDelta[c_,b_,ops:OptionsPattern[]]]:=FlavorDelta[a,c,ops]/;(OptionValue[HoldFlavor]===False)
FlavorDelta/:Times[FlavorDelta[b_,a_,ops:OptionsPattern[]],FlavorDelta[b_,c_,ops:OptionsPattern[]]]:=FlavorDelta[a,c,ops]/;(OptionValue[HoldFlavor]===False)
FlavorDelta/:Times[FlavorDelta[b_,a_,ops:OptionsPattern[]],FlavorDelta[c_,b_,ops:OptionsPattern[]]]:=FlavorDelta[a,c,ops]/;(OptionValue[HoldFlavor]===False)*)



FlavorDelta[i_,j_,ops:OptionsPattern[]]:=FlavorDelta[j,i,ops]/;(i=!=j)&&(!OrderedQ[{i,j}])
FlavorDelta[i_,j_]:=FlavorDelta[i,j,HoldFlavor->True]

FlavorDelta/: MakeBoxes[FlavorDelta[i_,j_,OptionsPattern[]],TraditionalForm]:=SubscriptBox[ToBoxes["\[Delta]",TraditionalForm],RowBox[{ToBoxes[i,TraditionalForm],ToBoxes[j,TraditionalForm]}]]/;OrderedQ[{i,j}]

FlavorDelta[i_,j_,OptionsPattern[]]:=If[i===j,1,0]/;(OptionValue[HoldFlavor]===False)




End[]
