(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


FourquarkLeading::usage="FourquarkLeading[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}}] give leading perturbative contribution of <J1 J2>, with J1=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a1\)]\)\!\(\*SubscriptBox[\(v1\[CapitalPsi]\), \(b1\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c1\)]\)\!\(\*SubscriptBox[\(v2\[CapitalPsi]\), \(d1\)]\) and J2=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(b2\)]\)\!\(\*SubscriptBox[\(v3\[CapitalPsi]\), \(a2\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(d2\)]\)\!\(\*SubscriptBox[\(v4\[CapitalPsi]\), \(c2\)]\) "

FourquarkLeading::inderr="Dummy indices conflict!"



Begin["`Private`FourquarkLeading`"]

Options[FourquarkLeading] = {
	HoldFlavor->False,
	EpsOrder->0,
	Pole->1,
	AutoNDR->False,
	Parallelized->False
}


(*------------------------------------------------------------------*)
(* propagators *)

prop[q_]=I GSD[q]FAD[q];
prop[q_,a_]=I (GSD[q]+quarkMass[a])FAD[q];

(*xprop[x_]=FourierPX[prop[q],{q,x}];*)
xprop[x_]=1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];

(*xprop[x_,a_]=FourierPX[prop[q,a],{q,x}];*)
xprop[x_,a_]=(Pair[Momentum[x,D],Momentum[x,D]]^(1-D/2) quarkMass[a]qfact1[1/4 I^(2-D) \[Pi]^(-D/2) qGamma[-1+D/2]] 
			+ DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qfact1[1/2 I^(1-D) \[Pi]^(-D/2) qGamma[D/2]]);
			


(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
(* leading diagram *)



FourquarkLeading[qq_,{v1_,v2_,{f11_,f12_,f13_,f14_}},{v3_,v4_,{f21_,f22_,f23_,f24_}},OptionsPattern[]]:=Block[
{null,tmp,tmp1,tmp2,tra1,dia1,tra2,dia2,k1,k2,k3,x,q,hv3,hv4,pole,holdf=OptionValue[HoldFlavor],atr,fdir,pall,ord=OptionValue[EpsOrder],files,dot},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[FourquarkLeading::inderr];
	Abort[]
];



pole=If[!FreeQ[OptionValue[Pole],Epsilon],OptionValue[Pole],1];



(* B A^+ B *)
hv3=v3//ComplexConjugate//FCI;
hv4=v4//ComplexConjugate//FCI;
(* ComplexConjugate[sigma] will expand the DiracSigma to gamma matrices, recover it to DiracSigma *)
hv4=(hv4/.Dot->dot)/.f1_ dot[aa___,DiracGamma[LorentzIndex[lor1_,dim___],dim___],DiracGamma[LorentzIndex[lor2_,dim___],dim___],bb___]+f2_ dot[aa___,DiracGamma[LorentzIndex[lor2_,dim___],dim___],DiracGamma[LorentzIndex[lor1_,dim___],dim___],bb___]/;(f1+f2===0):>-2I f1 dot[aa,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],bb];
hv4=hv4/.{null->1,dot->Dot};


If[OptionValue[AutoNDR]===True,atr=TR5,atr=TR];

(*----------------------------------------------------------------*)

If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];


(*------------------------------------------------------------------*)

If[pall===True,

	DistributeDefinitions[v1,v2,f11,f12,f13,f14,f21,f22,f23,f24];
	tmp= WaitAll[ParallelSubmit[{hv3,hv4,holdf,atr,ord,pole},
									leading[qq,{v1,v2,{f11,f12,f13,f14}},{hv3,hv4,{f21,f22,f23,f24}},holdf,atr,ord,pole]]
				];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,
	If[pall==="External",

		DistributeDefinitions[v1,v2,f11,f12,f13,f14,f21,f22,f23,f24];

		If[fdir=="None",
		(* evaluation, no import and export *)
			{ParallelSubmit[{hv3,hv4,holdf,atr,ord,pole},
									leading[qq,{v1,v2,{f11,f12,f13,f14}},{hv3,hv4,{f21,f22,f23,f24}},holdf,atr,ord,pole]]
			}

		,
		(* evaluation, import and export *)
			files=(StringSplit[ToString[#],"`"][[-1]])&/@{leading};
			files=("FourquarkLeading_"<>#)&/@files;
			
			
			ImExport[fdir,
						files,
						{{hv3,hv4,holdf,atr,ord,pole},
						{qq,{v1,v2,{f11,f12,f13,f14}},{hv3,hv4,{f21,f22,f23,f24}},holdf,atr,ord,pole},
						{leading}}
						]
		](* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
	
	,
	
		
		leading[qq,{v1,v2,{f11,f12,f13,f14}},{hv3,hv4,{f21,f22,f23,f24}},holdf,atr,ord,pole]/.CA-2CF->1/CA

	]
]

]

(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)


leading[qq_,{v1_,v2_,{f11_,f12_,f13_,f14_}},{hv3_,hv4_,{f21_,f22_,f23_,f24_}},holdf_,atr_,ord_,pole_]:=Block[{x,q,tra1,tra2,dia1,dia2},


(*------------------------------------------------------------------*)
(* figure *)

tra1=(FlavorDelta[f11,f21]FlavorDelta[f12,f22]FlavorDelta[f13,f23]FlavorDelta[f14,f24]atr[SUNTrace[v1 . xprop[x] . hv3 . xprop[-x]]]atr[SUNTrace[v2 . xprop[x] . hv4 . xprop[-x]]]
+
(* \[CapitalGamma]_3 \[TwoWayRule] \[CapitalGamma]_4] *)
FlavorDelta[f11,f23]FlavorDelta[f12,f24]FlavorDelta[f13,f21]FlavorDelta[f14,f22]atr[SUNTrace[v1 . xprop[x] . hv4 . xprop[-x]]]atr[SUNTrace[v2 . xprop[x] . hv3 . xprop[-x]]]);


dia1=FourierXP[tra1//SUNSimplify,{x,q}];
dia1=QEvaluate[I ScaleMu^(3(4-D))dia1 pole,q,EpsOrder->ord,Parallelized->False,HoldFlavor->holdf];


(*------------------------------------------------------------------*)
(* figure *)
tra2=(FlavorDelta[f12,f22]FlavorDelta[f21,f13]FlavorDelta[f14,f24]FlavorDelta[f23,f11]atr[SUNTrace[v1 . xprop[x] . hv3 . xprop[-x] . v2 . xprop[x] . hv4 . xprop[-x]]]
+
(* \[CapitalGamma]_3 \[TwoWayRule] \[CapitalGamma]_4 *)
FlavorDelta[f12,f24]FlavorDelta[f23,f13]FlavorDelta[f14,f22]FlavorDelta[f21,f11]atr[SUNTrace[v1 . xprop[x] . hv4 . xprop[-x] . v2 . xprop[x] . hv3 . xprop[-x]]]);


dia2=FourierXP[tra2//SUNSimplify,{x,q}];
dia2=QEvaluate[-I ScaleMu^(3(4-D))dia2 pole,q,EpsOrder->ord,Parallelized->False,HoldFlavor->holdf];(* -1 for odd number of fermion loops *)


(*------------------------------------------------------------------*)
QGather[dia1+dia2//SUNSimplify,q,ShowasTable->False]/.q->qq


]



End[]
(*EndPackage[]*)
