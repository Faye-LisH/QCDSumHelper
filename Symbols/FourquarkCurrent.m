(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


FourquarkCurrent::usage=""
FourquarkCurrent::colerr="Not a color signlet current!"



Begin["`Private`FourquarkCurrent`"]

Options[FourquarkCurrent] = {
	};




(*SetAttributes[FourquarkCurrent,Orderless];*)
(*FourquarkCurrent[aa_,bb_]:=Block[{factor=1,tmpa,tmpb,rule,dot,null},

If[simplify===True,
	tmpa=(aa//FCI)/.Dot\[Rule]dot;
	tmpb=(bb//FCI)/.Dot\[Rule]dot;

	rule=Flatten[Cases[#,DiracSigma[DiracGamma[LorentzIndex[lora_,dim___],dim___],DiracGamma[LorentzIndex[lorb_,dim___],dim___]]\[RuleDelayed]
				{dot[cc___,DiracGamma[LorentzIndex[lora,dim2___],dim2___],DiracGamma[LorentzIndex[lorb,dim2___],dim2___],dd___]/;MatchQ[null[dim2],null[D]|null[]]\[RuleDelayed](factor=-I factor;dot[cc,DiracSigma[DiracGamma[LorentzIndex[lora,dim2],dim2],DiracGamma[LorentzIndex[lorb,dim2],dim2]],dd]),
				dot[cc___,DiracGamma[LorentzIndex[lorb,dim2___],dim2___],DiracGamma[LorentzIndex[lora,dim2___],dim2___],dd___]/;MatchQ[null[dim2],null[D]|null[]]\[RuleDelayed](factor=I factor;dot[cc,DiracSigma[DiracGamma[LorentzIndex[lora,dim2],dim2],DiracGamma[LorentzIndex[lorb,dim2],dim2]],dd])},Infinity]]&/@{tmpa,tmpb};

(* Replace, e.g Current[_, GAD[u,v], _] Current[_, Sigma[u,v], _] -> -i Current[_, Sigma[u,v], _] Current[_, Sigma[u,v], _] *)

	tmpa=tmpa/.rule[[2]]/.dot\[Rule]Dot;
	tmpb=tmpb/.rule[[1]]/.dot\[Rule]Dot;

	factor FourquarkCurrent[tmpa,tmpb]
,

	FourquarkCurrent[aa,bb]
]

]*)




FourquarkCurrent[-aa_Current,bb_Current]:=-FourquarkCurrent[aa,bb]
FourquarkCurrent[aa_Current,-bb_Current]:=-FourquarkCurrent[aa,bb]
FourquarkCurrent[-aa_Current,-bb_Current]:=FourquarkCurrent[aa,bb]


FourquarkCurrent/:MakeBoxes[FourquarkCurrent[aa_Current,bb_Current],TraditionalForm]:=Block[{tmpcol},



If[And@@(MatchQ[#,Current[{_},{_}]|Current[Except[_List],Except[_List]]|Current[{_},_,{_}]|Current[Except[_List],_,Except[_List]]]&/@{aa,bb}),
(* meson-meson current *)

	
	RowBox[{ToBoxes[aa,TraditionalForm],ToBoxes[bb,TraditionalForm]}]

,
(* diquark-antidiquark current *)
  (* map color-index in quark to color-index, color-index in antiquark to 1/(color-index); for color signlect, the product of them = 1 *)
	tmpcol=aa bb/.{Current[{_,col1_},{_,col2_}]:>col1 col2,Current[{_,col1_},_,{_,col2_}]:>col1 col2,Current[_,col1_,_,col2_]:>col1 col2,Current[_,col1_,_,_,col2_]:>col1 col2,
					Current[{-1,_,col1_},{-1,_,col2_}]:>1/(col1 col2),Current[{-1,_,col1_},_,{-1,_,col2_}]:>1/(col1 col2),Current[-1,_,col1_,-1,_,col2_]:>1/(col1 col2),
					Current[-1,_,col1_,_,-1,_,col2_]:>1/(col1 col2)};
	
	If[tmpcol=!=1,
		Message[FourquarkCurrent::colerr];Abort[]
	,
		Which[
			MatchQ[aa,Current[{_,_},{_,_}]|Current[_,_,_,_]|Current[{_,_},_,{_,_}]|Current[_,_,_,_,_]]&&MatchQ[bb,Current[{-1,_,_},{-1,_,_}]|Current[-1,_,_,-1,_,_]|Current[{-1,_,_},_,{-1,_,_}]|Current[-1,_,_,_,-1,_,_]],

			RowBox[{ToBoxes[aa,TraditionalForm],ToBoxes[bb,TraditionalForm]}]
		,
			MatchQ[bb,Current[{_,_},{_,_}]|Current[_,_,_,_]|Current[{_,_},_,{_,_}]|Current[_,_,_,_,_]]&&MatchQ[aa,Current[{-1,_,_},{-1,_,_}]|Current[-1,_,_,-1,_,_]|Current[{-1,_,_},_,{-1,_,_}]|Current[-1,_,_,_,-1,_,_]],
		
			FourquarkCurrent[bb,aa]
			(*RowBox[{ToBoxes[bb,TraditionalForm],ToBoxes[aa,TraditionalForm]}]*)
		]

	]
]

]






End[]
(*EndPackage[]*)
