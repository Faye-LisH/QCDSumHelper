(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


Current::usage=""
Current::inderr="Dummy indices conflict!"



Begin["`Private`Current`"]





Current[aa_,-1*gg_,bb_]:=-Current[aa,gg,bb]

Current[f1_,col1_,-1*g_,f2_,col2_]:=-Current[f1,col1,g,f2,col2]

Current[-1,f1_,col1_,-1*g_,-1,f2_,col2_]:=-Current[-1,f1,col1,g,-1,f2,col2]



(*Current[aa_List]:=Current@@aa;*)

(*------------------------------------------*)
(* Meson current *)
Current/:MakeBoxes[Current[{f1_},{f2_}],TraditionalForm]:=RowBox[{SubscriptBox[OverscriptBox["\[CapitalPsi]","_"],ToBoxes[f1,TraditionalForm]],SubscriptBox["\[CapitalPsi]",ToBoxes[f2,TraditionalForm]]}]

Current/:MakeBoxes[Current[f1:Except[_List],f2:Except[_List]],TraditionalForm]:=RowBox[{OverscriptBox[ToBoxes[f1,TraditionalForm],"_"],ToBoxes[f2,TraditionalForm]}]

Current/:MakeBoxes[Current[{f1_},g_,{f2_}],TraditionalForm]:=RowBox[{SubscriptBox[OverscriptBox["\[CapitalPsi]","_"],ToBoxes[f1,TraditionalForm]],ToBoxes[g,TraditionalForm],SubscriptBox["\[CapitalPsi]",ToBoxes[f2,TraditionalForm]]}]

Current/:MakeBoxes[Current[f1:Except[_List],g_,f2:Except[_List]],TraditionalForm]:=RowBox[{OverscriptBox[ToBoxes[f1,TraditionalForm],"_"],ToBoxes[g,TraditionalForm],ToBoxes[f2,TraditionalForm]}]



(*------------------------------------------*)
(* quark-quark current *)
Current/:MakeBoxes[Current[{f1_,col1_},{f2_,col2_}]/;col1=!=col2,TraditionalForm]:=RowBox[{SuperscriptBox[SubsuperscriptBox["\[CapitalPsi]",ToBoxes[f1,TraditionalForm],ToBoxes[col1,TraditionalForm]]," T"],"C",SubsuperscriptBox["\[CapitalPsi]",ToBoxes[f2,TraditionalForm],ToBoxes[col2,TraditionalForm]]}]

Current/:MakeBoxes[Current[f1_,col1_,f2_,col2_]/;col1=!=col2,TraditionalForm]:=RowBox[{SubsuperscriptBox[ToBoxes[f1,TraditionalForm],ToBoxes[col1,TraditionalForm],"T"],"C",SubscriptBox[ToBoxes[f2,TraditionalForm],ToBoxes[col2,TraditionalForm]]}]

Current/:MakeBoxes[Current[{f1_,col1_},g_,{f2_,col2_}]/;col1=!=col2,TraditionalForm]:=RowBox[{SuperscriptBox[SubsuperscriptBox["\[CapitalPsi]",ToBoxes[f1,TraditionalForm],ToBoxes[col1,TraditionalForm]]," T"],"C",ToBoxes[g,TraditionalForm],SubsuperscriptBox["\[CapitalPsi]",ToBoxes[f2,TraditionalForm],ToBoxes[col2,TraditionalForm]]}]

Current/:MakeBoxes[Current[f1_,col1_,g_,f2_,col2_]/;col1=!=col2,TraditionalForm]:=RowBox[{SubsuperscriptBox[ToBoxes[f1,TraditionalForm],ToBoxes[col1,TraditionalForm],"T"],"C",ToBoxes[g,TraditionalForm],SubscriptBox[ToBoxes[f2,TraditionalForm],ToBoxes[col2,TraditionalForm]]}]




(*------------------------------------------*)
(* antiquark-antiquark current *)
Current/:MakeBoxes[Current[{-1,f1_,col1_},{-1,f2_,col2_}]/;col1=!=col2,TraditionalForm]:=RowBox[{SubsuperscriptBox[OverscriptBox["\[CapitalPsi]","_"],ToBoxes[f1,TraditionalForm],ToBoxes[col1,TraditionalForm]],"C",SuperscriptBox[SubsuperscriptBox[OverscriptBox["\[CapitalPsi]","_"],ToBoxes[f2,TraditionalForm],ToBoxes[col2,TraditionalForm]],"T"]}]

Current/:MakeBoxes[Current[-1,f1_,col1_,-1,f2_,col2_]/;col1=!=col2,TraditionalForm]:=RowBox[{SubscriptBox[OverscriptBox[ToBoxes[f1,TraditionalForm],"_"],ToBoxes[col1,TraditionalForm]],"C",SuperscriptBox[SubscriptBox[OverscriptBox[ToBoxes[f2,TraditionalForm],"_"],ToBoxes[col2,TraditionalForm]],"T"]}]

Current/:MakeBoxes[Current[{-1,f1_,col1_},g_,{-1,f2_,col2_}]/;col1=!=col2,TraditionalForm]:=RowBox[{SubsuperscriptBox[OverscriptBox["\[CapitalPsi]","_"],ToBoxes[f1,TraditionalForm],ToBoxes[col1,TraditionalForm]],ToBoxes[g,TraditionalForm],"C",SuperscriptBox[SubsuperscriptBox[OverscriptBox["\[CapitalPsi]","_"],ToBoxes[f2,TraditionalForm],ToBoxes[col2,TraditionalForm]],"T"]}]

Current/:MakeBoxes[Current[-1,f1_,col1_,g_,-1,f2_,col2_]/;col1=!=col2,TraditionalForm]:=RowBox[{SubscriptBox[OverscriptBox[ToBoxes[f1,TraditionalForm],"_"],ToBoxes[col1,TraditionalForm]],ToBoxes[g,TraditionalForm],"C",SuperscriptBox[SubscriptBox[OverscriptBox[ToBoxes[f2,TraditionalForm],"_"],ToBoxes[col2,TraditionalForm]],"T"]}]





End[]
(*EndPackage[]*)
