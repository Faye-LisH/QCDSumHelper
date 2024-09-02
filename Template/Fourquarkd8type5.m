(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
Fourquarkd8type5::usage="Fourquarkd8type5[q_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}}] give the \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)Gq\[RightAngleBracket]\[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)q\[RightAngleBracket] contribution of OPE result about <J1 J2>, with J1=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a1\)]\)\!\(\*SubscriptBox[\(v1\[CapitalPsi]\), \(b1\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c1\)]\)\!\(\*SubscriptBox[\(v2\[CapitalPsi]\), \(d1\)]\) and J2=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(b2\)]\)\!\(\*SubscriptBox[\(v3\[CapitalPsi]\), \(a2\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(d2\)]\)\!\(\*SubscriptBox[\(v4\[CapitalPsi]\), \(c2\)]\), here the background gluon comes from quark propatator. "

Fourquarkd8type5::inderr="Dummy indices conflict!"

Fourquarkd8type5::curerr="Unknow current structure!"



(* ::Code::Initialization::Plain:: *)
Begin["`Private`Fourquarkd8type5`"]




(* ::Code::Initialization::Plain:: *)
Options[Fourquarkd8type5] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True
}



(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(* propagators *)
(*
prop[q_]= I GSD[q]FAD[q];*)


(* ::Code::Initialization::Plain:: *)
(* xprop[x_]=FourierPX[prop[q],{q,x}]; *)
xprop[x_]= 1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];



(* ::Code::Initialization::Plain:: *)
(*  quark propagator with background gluon  *)
(*  ---<--- q *)
(*prog[q_,lora_,lorb_,cola_]=-1/2 FAD[{q}] GSD[q].GAD[lorb].SUNT[cola].FourDivergence[GSD[q]FAD[{q}],FVD[q,lora]];*)
(*prog[q_,lora_,lorb_,cola_]= -1/2 FAD[{q}]GSD[q] . GAD[lorb] . SUNT[cola] . (FAD[q] GAD[lora]- 2FAD[q,q] FVD[q,lora] GSD[q]);*)


(* ::Code::Initialization::Plain:: *)
(* x ---<--- 0 *)
(*xprog[x_,lora_,lorb_,cola_]=FourierPX[prog[q,lora,lorb,cola],{q,x}];*)
(*xprog[x_,lora_,lorb_,cola_]= (GAD[lora] . GAD[lorb] . GSD[x] qfact1[1/32 I^(-D) \[Pi]^(-D/2) qGamma[-1+D/2]] SPD[x,x]^(1-D/2) SUNT[cola]
							+GSD[x] . GAD[lorb] . GAD[lora] qfact1[qfact2[-(1/32) I^(-D) \[Pi]^(-D/2)] qGamma[-1+D/2]] SPD[x,x]^(1-D/2) SUNT[cola]
							+2 FVD[x,lora] FVD[x,lorb] GSD[x] qfact1[-(1/16) I^(-D) \[Pi]^(-D/2) qGamma[D/2]] SPD[x,x]^(-D/2) SUNT[cola]);
							*)


(* ::Code::Initialization::Plain:: *)
(*--- Condensate ---*)
(*qgq[lora_,lorb_,cola_]=  SUNT[cola] . DiracSigma[GAD[lora,lorb]];*)

cqgq[a_,b_]=-1/(CA CF 4(D-1)D)Condensate[{a,"G",b}];
cqq[f1_,f2_]=-1/(4CA)Condensate[{f1,f2}];		

(*-------------------------------------------------*)
(* gamma^uva in Eq.18 in arxiv 2111.13897; when contract with x^u x^v and antisymmertrize about a <-> b *)
ggv = 1/(2*2(D+2))(SPD[x] DiracSigma[GAD[lora,lorb]] + I(FVD[x,lora]GSD[x] . GAD[lorb]-FVD[x,lorb]GSD[x] . GAD[lora]));
ggv2 = 1/(2(D+2))(SPD[x] DiracSigma[GAD[lora,lorb]] + I(FVD[x,lora]GSD[x] . GAD[lorb]-FVD[x,lorb]GSD[x] . GAD[lora]));
(*ggv=1/2 I/(2(D+2))(SPD[x](GAD[lora,lorb]-GAD[lorb,lora])/2 + (FVD[x,lora]GSD[x] . GAD[lorb]-FVD[x,lorb]GSD[x] . GAD[lora]));
ggv2=1/2 I/(D+2)(SPD[x](GAD[lora,lorb]-GAD[lorb,lora])/2 + (FVD[x,lora]GSD[x] . GAD[lorb]-FVD[x,lorb]GSD[x] . GAD[lora]));
*)


Fourquarkd8type5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmp,tmp1,tmp2,hv3,hv4,holdf=OptionValue[HoldFlavor],atr,diagrams,fdir,files,pall,dot},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[Fourquarkd8type5::inderr];
	Abort[]
];




If[OptionValue[AutoNDR]===True,
	atr=TR5
	,
	atr=TR
];


(* B A^+ B *)
hv3=v3//ComplexConjugate;
hv4=v4//ComplexConjugate;
(* ComplexConjugate[sigma] will expand the DiracSigma to gamma matrices, recover it to DiracSigma *)
hv4=(hv4/.Dot->dot)/.f1_ dot[aa___,DiracGamma[LorentzIndex[lor1_,dim___],dim___],DiracGamma[LorentzIndex[lor2_,dim___],dim___],bb___]+f2_ dot[aa___,DiracGamma[LorentzIndex[lor2_,dim___],dim___],DiracGamma[LorentzIndex[lor1_,dim___],dim___],bb___]/;(f1+f2===0):>-2I f1 dot[aa,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],bb];
hv4=hv4 /. {null->1,dot->Dot};


(*--------------------------------------------*)		

					
If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];


diagrams={xtype1,xtype2,xtype3,xtype4,xtype5,xtype6,xtype7,xtype8,xtype9,xtype10,xtype11,xtype12};
(*{type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14};*)


(*---------------------------------------------------*)
If[pall===True,

	DistributeDefinitions[qq,v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
	tmp=Plus@@WaitAll[ParallelSubmit[{hv3,hv4,holdf,atr},#[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr]]&/@diagrams];
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.{CA-2CF->1/CA,2CA^3CF-CA^2(4CF^2+1)+1->0}
	
,

	If[pall==="External",
	
		DistributeDefinitions[qq,v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
		
		If[fdir==="None",
		
			ParallelSubmit[{hv3,hv4,holdf,atr},#[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr]]&/@diagrams
		,
		
			files=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files=("Fourquarkd8type5_"<>#)&/@files;
			
			ImExport[fdir,
						files,
						{{hv3,hv4,holdf,atr},
						{qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr},
						diagrams}
						]
		](* !!! do not touch tmp before WaitAll[tmp], any function may modifiy the things in EvaluationObject[...] should not be used here. !!! *)
	,
	

	
		tmp=Plus@@Map[#[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr]&,diagrams];
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.{CA-2CF->1/CA,2CA^3CF-CA^2(4CF^2+1)+1->0}
	]
]

]



(* ::Code::Initialization::Plain:: *)
(*--- The following are generated by algorithm ---*)


(* ::Input::Initialization:: *)
xtype1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,cola,lora,lorb,str,tr,dot,contract},


dia1=cqgq[b1, b2]*cqq[d1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v2]]]*tr[str[dot[hv3, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[d1, b2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, v1, xprop[x]]]]*tr[str[dot[hv3, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[a2, c1]*cqq[d1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]] + cqgq[a2, a1]*cqq[b1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, ggv . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]] - cqgq[d1, d2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[d1, b2]*cqq[b1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, ggv . SUNT[cola], hv4, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[c2, c1]*cqq[d1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, ggv . SUNT[cola]]]] - cqgq[d1, d2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, ggv . SUNT[cola], v1, xprop[x], hv3, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia1=(FourierXP[dia1,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(4-D) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,cola,lora,lorb,str,tr,dot,contract},


dia2=cqgq[c2, a1]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, v2, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x]]]] + cqgq[c2, a1]*cqq[d1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x]]]] + cqgq[b1, d2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[b1, b2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, ggv . SUNT[cola], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]] - cqgq[d1, b2]*cqq[b1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, hv4, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[b1, b2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]], hv3, xprop[-x], v2, xprop[x]]]] - cqgq[a2, c1]*cqq[d1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2]]] - cqgq[d1, b2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, ggv . SUNT[cola], v1, xprop[x], hv4, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia2=(FourierXP[dia2,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(4-D) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,cola,lora,lorb,str,tr,dot,contract},


dia3=cqgq[a2, c1]*cqq[b1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, xprop[-x], v1]]]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] + cqgq[b1, b2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, ggv . SUNT[cola], v2, xprop[x]]]] + cqgq[d1, b2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, ggv . SUNT[cola], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]] - cqgq[c2, c1]*cqq[b1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] - cqgq[b1, d2]*cqq[d1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]], hv4, xprop[-x], v2]]] - cqgq[c2, c1]*cqq[d1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2]]] - cqgq[c2, c1]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, ggv . SUNT[cola], v1, xprop[x], hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]] - cqgq[b1, b2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola], hv3, ggv . SUNT[cola], v2, xprop[x]]]];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia3=(FourierXP[dia3,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=QEvaluate[I ScaleMu^(4-D) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,cola,lora,lorb,str,tr,dot,contract},


dia4=cqgq[c2, c1]*cqq[d1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2]]]*tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]] + cqgq[b1, d2]*cqq[d1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2, ggv . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[a2, a1]*cqq[b1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, hv4, xprop[-x], v2, xprop[x]]]] - cqgq[c2, a1]*cqq[b1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, hv3, xprop[-x], v2, xprop[x]]]] - cqgq[d1, d2]*cqq[b1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, hv3, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[a2, c1]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, ggv . SUNT[cola], v1, xprop[x], hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]] - cqgq[c2, a1]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x], hv3, ggv . SUNT[cola], v2, xprop[x]]]] - cqgq[c2, a1]*cqq[d1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x], hv3, xprop[-x], v2, ggv . SUNT[cola]]]];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia4=(FourierXP[dia4,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=QEvaluate[I ScaleMu^(4-D) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,cola,lora,lorb,str,tr,dot,contract},


dia5=cqgq[d1, b2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, ggv . SUNT[cola], v1, xprop[x]]]] + cqgq[c2, c1]*cqq[b1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] + cqgq[c2, a1]*cqq[d1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2, ggv . SUNT[cola]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x]]]] + cqgq[b1, d2]*cqq[d1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2]]]*tr[str[dot[hv4, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[d1, d2]*cqq[b1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, ggv . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[b1, d2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]], hv4, v2, xprop[x]]]] - cqgq[a2, a1]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x], hv4, ggv . SUNT[cola], v2, xprop[x]]]] - cqgq[b1, d2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola], hv4, ggv . SUNT[cola], v2, xprop[x]]]];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia5=(FourierXP[dia5,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=QEvaluate[I ScaleMu^(4-D) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,cola,lora,lorb,str,tr,dot,contract},


dia6=cqgq[a2, a1]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, v2, xprop[x]]]]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x]]]] + cqgq[d1, d2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv4, ggv . SUNT[cola], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] + cqgq[b1, d2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, ggv . SUNT[cola], v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] + cqgq[a2, a1]*cqq[d1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, ggv . SUNT[cola]]]] + cqgq[d1, d2]*cqq[b1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1]]]*tr[str[dot[hv4, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[d1, d2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, ggv . SUNT[cola], v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[c2, a1]*cqq[d1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x], hv3, xprop[-x], v2]]] - cqgq[b1, b2]*cqq[d1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]], hv3, xprop[-x], v2]]];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia6=(FourierXP[dia6,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia6=QEvaluate[I ScaleMu^(4-D) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,cola,lora,lorb,str,tr,dot,contract},


dia7=cqgq[a2, c1]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]]*tr[str[dot[hv4, ggv . SUNT[cola], v1, xprop[x]]]] + cqgq[a2, a1]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x]]]]*tr[str[dot[hv4, ggv . SUNT[cola], v2, xprop[x]]]] + cqgq[c2, a1]*cqq[b1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, ggv . SUNT[cola]]]] - cqgq[a2, a1]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x], hv4, v2, xprop[x]]]] - cqgq[d1, b2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[b1, b2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, ggv . SUNT[cola], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola], hv3, xprop[-x], v2, xprop[x]]]] - cqgq[c2, a1]*cqq[b1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, ggv . SUNT[cola], hv3, xprop[-x], v2, xprop[x]]]] - cqgq[a2, c1]*cqq[d1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, ggv . SUNT[cola]]]];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia7=(FourierXP[dia7,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia7=QEvaluate[I ScaleMu^(4-D) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,cola,lora,lorb,str,tr,dot,contract},


dia8=cqgq[a2, a1]*cqq[d1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v2]]]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x]]]] + cqgq[b1, d2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]]*tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]] + cqgq[c2, c1]*cqq[d1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, ggv . SUNT[cola]]]] + cqgq[d1, b2]*cqq[b1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v1, ggv . SUNT[cola]]]] - cqgq[c2, a1]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x], hv3, v2, xprop[x]]]] - cqgq[a2, a1]*cqq[b1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, ggv . SUNT[cola], hv4, xprop[-x], v2, xprop[x]]]] - cqgq[c2, c1]*cqq[b1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, ggv . SUNT[cola], hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]] - cqgq[b1, d2]*cqq[d1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola], hv4, xprop[-x], v2, ggv . SUNT[cola]]]];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia8=(FourierXP[dia8,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia8=QEvaluate[I ScaleMu^(4-D) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,cola,lora,lorb,str,tr,dot,contract},


dia9=cqgq[c2, a1]*cqq[b1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1]]]*tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]] + cqgq[b1, d2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv4, ggv . SUNT[cola], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] + cqgq[c2, a1]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, ggv . SUNT[cola], v2, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x]]]] + cqgq[d1, d2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[c2, c1]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, v1, xprop[x], hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] - cqgq[a2, a1]*cqq[d1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x], hv4, xprop[-x], v2]]] - cqgq[d1, b2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, ggv . SUNT[cola], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[b1, b2]*cqq[d1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola], hv3, xprop[-x], v2, ggv . SUNT[cola]]]];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia9=(FourierXP[dia9,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia9=QEvaluate[I ScaleMu^(4-D) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,cola,lora,lorb,str,tr,dot,contract},


dia10=cqgq[a2, c1]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, v1, xprop[x]]]]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] + cqgq[d1, b2]*cqq[b1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, xprop[-x], v1]]]*tr[str[dot[hv3, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[c2, c1]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, v1, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] + cqgq[a2, c1]*cqq[b1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, ggv . SUNT[cola]]]] + cqgq[a2, c1]*cqq[d1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, ggv . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]] + cqgq[b1, b2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]] - cqgq[b1, d2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]], hv4, xprop[-x], v2, xprop[x]]]] - cqgq[a2, a1]*cqq[d1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x], hv4, xprop[-x], v2, ggv . SUNT[cola]]]];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia10=(FourierXP[dia10,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia10=QEvaluate[I ScaleMu^(4-D) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,cola,lora,lorb,str,tr,dot,contract},


dia11=cqgq[b1, b2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, v2, xprop[x]]]]*tr[str[dot[hv3, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[d1, d2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]]*tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]] + cqgq[c2, c1]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, ggv . SUNT[cola], v1, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]] - cqgq[d1, b2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, v1, xprop[x], hv4, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[d1, d2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, v1, xprop[x], hv3, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[a2, c1]*cqq[b1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] - cqgq[b1, b2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]], hv3, v2, xprop[x]]]] - cqgq[a2, c1]*cqq[b1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, ggv . SUNT[cola], hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia11=(FourierXP[dia11,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia11=QEvaluate[I ScaleMu^(4-D) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,cola,lora,lorb,str,tr,dot,contract},


dia12=cqgq[c2, c1]*cqq[b1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, ggv . SUNT[cola]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]] + cqgq[d1, b2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]] + cqgq[b1, b2]*cqq[d1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v2, ggv . SUNT[cola]]]] + cqgq[a2, a1]*cqq[b1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]] - cqgq[a2, c1]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, v1, xprop[x], hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] - cqgq[b1, d2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, ggv . SUNT[cola], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola], hv4, xprop[-x], v2, xprop[x]]]] - cqgq[d1, d2]*cqq[b1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, ggv . SUNT[cola], hv3, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[d1, d2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, ggv . SUNT[cola], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia12=(FourierXP[dia12,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia12=QEvaluate[I ScaleMu^(4-D) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
