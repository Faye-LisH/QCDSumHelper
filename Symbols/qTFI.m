(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



qTFI::usage = "qTFI[{u_,v_,r_,s_,t_},{v1_,v2_,v3_,v4_,v5_},{k1_,k2_,p_,{Rules___}}] not do any evalutaiton;
 it keep the {k1_,k2_,p_,{Rules___}} to indicate where the 2-loop comes form."
Begin["`Private`qTFI`"]


(*-------------------------------------------------------------------------------------------*)
(* show the TFI with momentum k1, k2, p and perhaps the rules of shift and/or rescale the k1, k2, to indicate how this 2-loop integral getted *)
qTFI[list1_,list2_,{k1_,k2_,p_,rules_List}]/;!FreeQ[rules,List,{2}]:=
	If[Plus@@list1===0,
		qTFI[Tarcer`TFI[D,(p/.{-aa_:>aa,aa_+bb_/;!FreeQ[aa,-1]&&!FreeQ[bb,-1]:>-aa-bb})^2,list2],
				{k1,k2,p/.{-aa_:>aa,aa_+bb_/;!FreeQ[aa,-1]&&!FreeQ[bb,-1]:>-aa-bb},DeleteCases[Flatten[rules],Rule[a_,a_]]}]
																	
	,
		qTFI[Tarcer`TFI[D,(p/.{-aa_:>aa,aa_+bb_/;!FreeQ[aa,-1]&&!FreeQ[bb,-1]:>-aa-bb})^2,list1,list2],
				{k1,k2,p/.{-aa_:>aa,aa_+bb_/;!FreeQ[aa,-1]&&!FreeQ[bb,-1]:>-aa-bb},DeleteCases[Flatten[rules],Rule[a_,a_]]}]
	]

																																				
(*-------------------*)

qTFI[aa_,{k1_,k2_,p_,{}}]:=qTFI[aa,{k1,k2,p}]


End[]
