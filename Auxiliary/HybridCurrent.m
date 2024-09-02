(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


HybirdCurrent::usage=""
HybirdCurrent::inderr="Not a Color Singlet Current!"



Begin["`Private`HybirdCurrent`"]



(*----------------------------------------------------------*)


HybirdCurrent[{f1_},{lora:Except[_List],lorb_,lorr___},{f2_}]:=HybirdCurrent[{f1},{lora,lorb,lorr},1,{f2}]

HybirdCurrent[{f1_},{lora:Except[_List],lorb_,lorr___},-1*v_,{f2_}]:=-HybirdCurrent[{f1},{lora,lorb,lorr},v,{f2}]

HybirdCurrent[{f1_},{lora:Except[_List],lorb_,lorr___},v_,{f2_}]:=Signature[{lora,lorb}]HybirdCurrent[{f1},{##,lorr},v,{f2}]&@@Sort[{lora,lorb}]/;!OrderedQ[{lora,lorb}]

HybirdCurrent/:MakeBoxes[HybirdCurrent[{f1_},{lora:Except[_List],lorb_,lorr___},v_,{f2_}],TraditionalForm]:=RowBox[{SubscriptBox[OverscriptBox["\[CapitalPsi]","_"],ToBoxes[f1,TraditionalForm]],ToBoxes[GluonStrength[lora,lorb,lorr]v,TraditionalForm],SubscriptBox["\[CapitalPsi]",ToBoxes[f2,TraditionalForm]]}]/;(OrderedQ[{lora,lorb}]&&FreeQ[v,SUNTF])




(*----------------------------------------------------------*)


HybirdCurrent[{f1_},{{lora_,lorb_,sun_},lorr___},{f2_}]:=HybirdCurrent[{f1},{{lora,lorb,sun},lorr},1,{f2}]

HybirdCurrent[{f1_},{{lora_,lorb_,sun_},lorr___},-1*v_,{f2_}]:=-HybirdCurrent[{f1},{{lora,lorb,sun},lorr},v,{f2}]

HybirdCurrent[{f1_},{{lora_,lorb_,sun_},lorr___},v_,{f2_}]:=Signature[{lora,lorb}]HybirdCurrent[{f1},{{##,sun},lorr},v,{f2}]&@@Sort[{lora,lorb}]/;!OrderedQ[{lora,lorb}]

HybirdCurrent/:MakeBoxes[HybirdCurrent[{f1_},{{lora_,lorb_,sun_},lorr___},v_,{f2_}],TraditionalForm]:=Block[{cola=$ac[Unique[]],colb=$ac[Unique[]]},RowBox[{SubsuperscriptBox[OverscriptBox["\[CapitalPsi]","_"],ToBoxes[f1,TraditionalForm],ToBoxes[cola,TraditionalForm]],ToBoxes[SUNTF[sun,cola,colb]GluonStrength[{lora,lorb,sun},lorr]v ,TraditionalForm],SubsuperscriptBox["\[CapitalPsi]",ToBoxes[f2,TraditionalForm],ToBoxes[colb,TraditionalForm]]}]]/;(OrderedQ[{lora,lorb}]&&FreeQ[v,SUNTF])




(*----------------------------------------------------------*)


HybirdCurrent[{f1_,c1_},{lora:Except[_List],lorb_,lorr___},{f2_,c2_}]:=HybirdCurrent[{f1,c1},{lora,lorb,lorr},1,{f2,c2}]

HybirdCurrent[{f1_,c1_},{lora:Except[_List],lorb_,lorr___},-1*v_,{f2_,c2_}]:=-HybirdCurrent[{f1,c1},{lora,lorb,lorr},v,{f2,c2}]

HybirdCurrent[{f1_,c1_},{lora:Except[_List],lorb_,lorr___},v_,{f2_,c2_}]:=Signature[{lora,lorb}]HybirdCurrent[{f1,c1},{##,lorr},v,{f2,c2}]&@@Sort[{lora,lorb}]/;!OrderedQ[{lora,lorb}]

HybirdCurrent/:MakeBoxes[HybirdCurrent[{f1_,c1_},{lora:Except[_List],lorb_,lorr___},v_,{f2_,c2_}],TraditionalForm]:=Block[{sun=$AC[Unique[]]},RowBox[{SubsuperscriptBox[OverscriptBox["\[CapitalPsi]","_"],ToBoxes[f1,TraditionalForm],ToBoxes[c1,TraditionalForm]],ToBoxes[SUNTF[sun,c1,c2]GluonStrength[{lora,lorb,sun},lorr]v,TraditionalForm],SubsuperscriptBox["\[CapitalPsi]",ToBoxes[f2,TraditionalForm],ToBoxes[c1,TraditionalForm]]}]]/;(OrderedQ[{lora,lorb}]&&FreeQ[v,SUNTF])




(*----------------------------------------------------------*)


HybirdCurrent[{f1_,c1_},{{lora_,lorb_,sun_},lorr___},{f2_,c2_}]:=HybirdCurrent[{f1,c1},{{lora,lorb,sun},lorr},1,{f2,c2}]

HybirdCurrent[{f1_,c1_},{{lora_,lorb_,sun_},lorr___},-1*v_,{f2_,c2_}]:=-HybirdCurrent[{f1,c1},{{lora,lorb,sun},lorr},v,{f2,c2}]

HybirdCurrent[{f1_,c1_},{{lora_,lorb_,sun_},lorr___},v_,{f2_,c2_}]:=Signature[{lora,lorb}]HybirdCurrent[{f1,c1},{{##,sun},lorr},v,{f2,c2}]&@@Sort[{lora,lorb}]/;!OrderedQ[{lora,lorb}]

HybirdCurrent/:MakeBoxes[HybirdCurrent[{f1_,c1_},{{lora_,lorb_,sun_},lorr___},v_,{f2_,c2_}],TraditionalForm]:=RowBox[{SubsuperscriptBox[OverscriptBox["\[CapitalPsi]","_"],ToBoxes[f1,TraditionalForm],ToBoxes[c1,TraditionalForm]],ToBoxes[SUNTF[sun,c1,c2]GluonStrength[{lora,lorb,sun},lorr]v,TraditionalForm],SubsuperscriptBox["\[CapitalPsi]",ToBoxes[f2,TraditionalForm],ToBoxes[c1,TraditionalForm]]}]/;(OrderedQ[{lora,lorb}]&&FreeQ[v,SUNTF])

HybirdCurrent/:MakeBoxes[HybirdCurrent[{f1_,c1_},{{lora_,lorb_,sun_},lorr___},v_,{f2_,c2_}],TraditionalForm]:=Block[{},
If[FreeQ[v,SUNTF[sun,c1,c2]],
	Message[HybirdCurrent::inderr];
	0
,
	RowBox[{SubsuperscriptBox[OverscriptBox["\[CapitalPsi]","_"],ToBoxes[f1,TraditionalForm],ToBoxes[c1,TraditionalForm]],ToBoxes[SUNTF[sun,c1,c2]GluonStrength[{lora,lorb,sun},lorr]v,TraditionalForm],SubsuperscriptBox["\[CapitalPsi]",ToBoxes[f2,TraditionalForm],ToBoxes[c1,TraditionalForm]]}]]
]/;(OrderedQ[{lora,lorb}]&&!FreeQ[v,SUNTF])


(*--------------------------------------------------------------------------------------------------*)
(*--------------------------------------------------------------------------------------------------*)
(*--------------------------------------------------------------------------------------------------*)


HybirdCurrent[f1:Except[_List],{lora:Except[_List],lorb_,lorr___},f2:Except[_List]]:=HybirdCurrent[f1,{lora,lorb,lorr},1,f2]

HybirdCurrent[f1:Except[_List],{lora:Except[_List],lorb_,lorr___},-1*v_,f2:Except[_List]]:=-HybirdCurrent[f1,{lora,lorb,lorr},v,f2]

HybirdCurrent[f1:Except[_List],{lora:Except[_List],lorb_,lorr___},v_,f2:Except[_List]]:=Signature[{lora,lorb}]HybirdCurrent[f1,{##,lorr},v,f2]&@@Sort[{lora,lorb}]/;!OrderedQ[{lora,lorb}]

HybirdCurrent/:MakeBoxes[HybirdCurrent[f1:Except[_List],{lora:Except[_List],lorb_,lorr___},v_,f2:Except[_List]],TraditionalForm]:=RowBox[{ToBoxes[f1,TraditionalForm],ToBoxes[GluonStrength[lora,lorb,lorr]v,TraditionalForm],ToBoxes[f2,TraditionalForm]}]/;(OrderedQ[{lora,lorb}]&&FreeQ[v,SUNTF])



(*----------------------------------------------------------*)



HybirdCurrent[f1:Except[_List],{{lora_,lorb_,sun_},lorr___},f2:Except[_List]]:=HybirdCurrent[f1,{{lora,lorb,sun},lorr},1,f2]

HybirdCurrent[f1:Except[_List],{{lora_,lorb_,sun_},lorr___},-1*v_,f2:Except[_List]]:=-HybirdCurrent[f1,{{lora,lorb,sun},lorr},v,f2]

HybirdCurrent[f1:Except[_List],{{lora_,lorb_,sun_},lorr___},v_,f2:Except[_List]]:=Signature[{lora,lorb}]HybirdCurrent[f1,{{##,sun},lorr},v,f2]&@@Sort[{lora,lorb}]/;!OrderedQ[{lora,lorb}]

HybirdCurrent/:MakeBoxes[HybirdCurrent[f1:Except[_List],{{lora_,lorb_,sun_},lorr___},v_,f2:Except[_List]],TraditionalForm]:=Block[{cola=$ac[Unique[]],colb=$ac[Unique[]]},RowBox[{SubscriptBox[ToBoxes[f1,TraditionalForm],ToBoxes[cola,TraditionalForm]],ToBoxes[SUNTF[sun,cola,colb]GluonStrength[{lora,lorb,sun},lorr]v ,TraditionalForm],SubscriptBox[ToBoxes[f2,TraditionalForm],ToBoxes[colb,TraditionalForm]]}]]/;(OrderedQ[{lora,lorb}]&&FreeQ[v,SUNTF])



(*----------------------------------------------------------*)



HybirdCurrent[f1:Except[_List],c1_,{lora:Except[_List],lorb_,lorr___},f2:Except[_List],c2_]:=HybirdCurrent[f1,c1,{lora,lorb,lorr},1,f2,c2]

HybirdCurrent[f1:Except[_List],c1_,{lora:Except[_List],lorb_,lorr___},-1*v_,f2:Except[_List],c2_]:=-HybirdCurrent[f1,c1,{lora,lorb,lorr},v,f2,c2]

HybirdCurrent[f1:Except[_List],c1_,{lora:Except[_List],lorb_,lorr___},v_,f2:Except[_List],c2_]:=Signature[{lora,lorb}]HybirdCurrent[f1,c1,{##,lorr},v,f2,c2]&@@Sort[{lora,lorb}]/;(!OrderedQ[{lora,lorb}]&&FreeQ[v,SUNTF])

HybirdCurrent/:MakeBoxes[HybirdCurrent[f1:Except[_List],c1_,{lora:Except[_List],lorb_,lorr___},v_,f2:Except[_List],c2_],TraditionalForm]:=Block[{sun=$AC[Unique[]]},RowBox[{SubscriptBox[ToBoxes[f1,TraditionalForm],ToBoxes[c1,TraditionalForm]],ToBoxes[SUNTF[sun,c1,c2]GluonStrength[{lora,lorb,sun},lorr]v,TraditionalForm],SubscriptBox[ToBoxes[f2,TraditionalForm],ToBoxes[c1,TraditionalForm]]}]]/;(OrderedQ[{lora,lorb}]&&FreeQ[v,SUNTF])



(*----------------------------------------------------------*)


HybirdCurrent[f1_,c1_,{{lora_,lorb_,sun_},lorr___},f2_,c2_]:=HybirdCurrent[f1,c1,{{lora,lorb,sun},lorr},1,f2,c2]

HybirdCurrent[f1_,c1_,{{lora_,lorb_,sun_},lorr___},-1*v_,f2_,c2_]:=-HybirdCurrent[f1,c1,{{lora,lorb,sun},lorr},v,f2,c2]

HybirdCurrent[f1_,c1_,{{lora_,lorb_,sun_},lorr___},v_,f2_,c2_]:=Signature[{lora,lorb}]HybirdCurrent[f1,c1,{{##,sun},lorr},v,f2,c2]&@@Sort[{lora,lorb}]/;!OrderedQ[{lora,lorb}]

HybirdCurrent/:MakeBoxes[HybirdCurrent[f1_,c1_,{{lora_,lorb_,sun_},lorr___},v_,f2_,c2_],TraditionalForm]:=RowBox[{SubscriptBox[ToBoxes[f1,TraditionalForm],ToBoxes[c1,TraditionalForm]],ToBoxes[SUNTF[sun,c1,c2]GluonStrength[{lora,lorb,sun},lorr]v,TraditionalForm],SubscriptBox[ToBoxes[f2,TraditionalForm],ToBoxes[c1,TraditionalForm]]}]/;(OrderedQ[{lora,lorb}]&&FreeQ[v,SUNTF])

HybirdCurrent/:MakeBoxes[HybirdCurrent[f1_,c1_,{{lora_,lorb_,sun_},lorr___},v_,f2_,c2_],TraditionalForm]:=Block[{},If[FreeQ[v,SUNTF[sun,c1,c2]],
	Message[HybirdCurrent::inderr];
	0
,
	RowBox[{SubscriptBox[ToBoxes[f1,TraditionalForm],ToBoxes[c1,TraditionalForm]],ToBoxes[GluonStrength[{lora,lorb,sun},lorr]v,TraditionalForm],SubscriptBox[ToBoxes[f2,TraditionalForm],ToBoxes[c1,TraditionalForm]]}]]
]/;(OrderedQ[{lora,lorb}]&&!FreeQ[v,SUNTF])


End[]
(*EndPackage[]*)
