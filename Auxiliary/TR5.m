(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



TR5::usage = "TR5[expr_,DiracGammaScheme->\"BMHV\"] is a spesific version of TR, which set DiracGammaScheme to \"NDR\" if even number of GA[5] involved, 
and set DiracGammaScheme to \"BMHV\" by default if even number of GA[5] in volved."

TR5::seperr = "Fail to separate the terms with even and odd number \!\(\*SuperscriptBox[OverscriptBox[\(\[Gamma]\), \(_\)], \(5\)]\), some DiracTrace will not been evaluated."


Begin["`Private`TR5`"]
Options[TR5] = {
	DiracGammaScheme->"BMHV",
	ToD4->False
}


TR5[expr_,OptionsPattern[]]:=Block[{null0,tmp,tmp1,tmp2,odd,tmp3,tmp01,scheme,setscheme=OptionValue[DiracGammaScheme],dot,plus2},


tmp=(expr//FCI);

If[OptionValue[ToD4]===True,
	tmp=tmp/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]};
	tmp=TR[tmp]

,(* only at D-dimension the gamma5 should be dealing seriously  *)
(* store the DiracGammaScheme *)
	scheme=FCGetDiracGammaScheme[];

(* to counter the number of gamma^5, the experssion should be expanded to the form: factor1*term1[dot product of gamma] + ... *)


(*--- split the SUNTrace[...] to 2 terms by the number of gamma^5 ---*)
	tmp=tmp/.SUNTrace[aa_/;!FreeQ[aa,DiracGamma|DiracSigma]]:>(tmp3=aa//Expand;
															If[Head[tmp3]===Plus,
																tmp3=List@@tmp3;
																tmp01=evenodd[#]&/@tmp3;
																tmp3=Plus@@(tmp01 tmp3);
																tmp01=Plus@@((1-tmp01)tmp3);

																SUNTrace[tmp3]+SUNTrace[tmp01]
															,
																SUNTrace[tmp3]]
															);


	If[!FreeQ[tmp,DiracTrace],tmp=tmp/.DiracTrace[aa_]:>TR5[aa,DiracGammaScheme->setscheme]];


	If[FreeQ[tmp,DiracGamma|DiracSigma],

		4 tmp

	,
	(*--- expand ---*)
		tmp=tmp/.aa_Plus/;FreeQ[aa,DiracGamma[5]]:>plus2@@aa;
		tmp=tmp//QExpand2;

		tmp=tmp/.{plus2->Plus};


		If[Head[tmp]===Plus,tmp=List@@tmp,tmp={tmp}];
		
		odd=evenodd[#]&/@tmp;
		
	(* ---------------------------- *)
		FCSetDiracGammaScheme[setscheme];
		tmp1=TR[Plus@@(odd tmp),DiracTraceEvaluate->True];

		FCSetDiracGammaScheme["NDR"];
		tmp2=TR[Plus@@((1-odd) tmp),DiracTraceEvaluate->True];

	(* recover the DiracGammaScheme *)
		FCSetDiracGammaScheme[scheme];
			
		tmp1+tmp2

		]
]

]





evenodd[expr_]:=Block[{tmp,dot,eo},
If[!FreeQ[expr,Dot],
	tmp=expr/.Dot->dot;
	tmp=tmp/.dot[ss__SUNT]->1;(* discard SUN matrices *)

	(*--- count the number of gamma^5 ---*)
	tmp=tmp/.tr_dot:>dot[Cases[tr,DiracGamma[5]->1,Infinity]];
	eo=Boole[OddQ[Cases[{tmp},dot[{n__Integer}]:>n,Infinity]]];(* must add braces, otherwise it doesn't work for tmp=dot[{_Integer}] *)

	
	eo=DeleteDuplicates[eo]/.{{}->0,{0}->0,{1}->1,{_,_}:>-1};
	
	(*--- if involve terms have even and odd number of gamma^5, keep the term with odd number of gamma^5 ("NDR" Scheme) ---*)
	If[eo===-1,Message[tr5::seperr];eo=0];
	
	eo
,	
	0
]		
]


End[]
