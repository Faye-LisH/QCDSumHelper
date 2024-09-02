(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


Fourquarkd8type2::usage="Fourquarkd8type2[q_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}}] give the \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)Gq\[RightAngleBracket]\[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)q\[RightAngleBracket] contribution of OPE result about <J1 J2>, with J1=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a1\)]\)\!\(\*SubscriptBox[\(v1\[CapitalPsi]\), \(b1\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c1\)]\)\!\(\*SubscriptBox[\(v2\[CapitalPsi]\), \(d1\)]\) and J2=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(b2\)]\)\!\(\*SubscriptBox[\(v3\[CapitalPsi]\), \(a2\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(d2\)]\)\!\(\*SubscriptBox[\(v4\[CapitalPsi]\), \(c2\)]\), here the background gluon comes from quark propatator. "

Fourquarkd8type2::inderr="Dummy indices conflict!"

Fourquarkd8type2::curerr="Unknow current structure!"



Begin["`Private`Fourquarkd8type2`"]


(* EM-first, contains <q_bar(x)q(x)> ~ <q_bar Gq> *)


Options[Fourquarkd8type2] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True
}


(*------------------------------------------------------------------*)
(* propagators *)
(*prop[q_]= I GSD[q]FAD[q];*)


(* xprop[x_]=FourierPX[prop[q],{q,x}]; *)
xprop[x_]= 1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];



(*  quark propagator with background gluon  *)
(*  ---<--- q *)
(*prog[q_,lora_,lorb_,cola_]=-1/2 FAD[{q}] GSD[q].GAD[lorb].SUNT[cola].FourDivergence[GSD[q]FAD[{q}],FVD[q,lora]];*)
(*prog[q_,lora_,lorb_,cola_]= -1/2 FAD[{q}]GSD[q] . GAD[lorb] . SUNT[cola] . (FAD[q] GAD[lora]- 2FAD[q,q] FVD[q,lora] GSD[q]);*)


(* x ---<--- 0 *)
(*xprog[x_,lora_,lorb_,cola_]=FourierPX[prog[q,lora,lorb,cola],{q,x}];*)
(*xprog[x_,lora_,lorb_,cola_]= (GAD[lora] . GAD[lorb] . GSD[x] qfact1[1/32 I^(-D) \[Pi]^(-D/2) qGamma[-1+D/2]] SPD[x,x]^(1-D/2) SUNT[cola]
							+GSD[x] . GAD[lorb] . GAD[lora] qfact1[qfact2[-(1/32) I^(-D) \[Pi]^(-D/2)] qGamma[-1+D/2]] SPD[x,x]^(1-D/2) SUNT[cola]
							+2 FVD[x,lora] FVD[x,lorb] GSD[x] qfact1[-(1/16) I^(-D) \[Pi]^(-D/2) qGamma[D/2]] SPD[x,x]^(-D/2) SUNT[cola]);
*)							


(*--- Condensate ---*)
(*
qgq[lora_,lorb_,cola_]=  SUNT[cola] . DiracSigma[GAD[lora,lorb]];*)

cqgq[a_,b_]:=-1/(CA CF 4(D-1)D)Condensate[{a,"G",b}]
cqq[f1_,f2_]:=-1/(4CA)Condensate[{f1,f2}];		

(*-------------------------------------------------*)

(* gamma^uva in Eq.18 in arxiv 2111.13897; when contract with x^u x^v and antisymmertrize about a <-> b *)

(*ggv=1/2 I/(2(D+2))(SPD[x](GAD[lora,lorb]-GAD[lorb,lora])/2 + (FVD[x,lora]GSD[x] . GAD[lorb]-FVD[x,lorb]GSD[x] . GAD[lora]));
ggv2=1/2 I/(D+2)(SPD[x](GAD[lora,lorb]-GAD[lorb,lora])/2 + (FVD[x,lora]GSD[x] . GAD[lorb]-FVD[x,lorb]GSD[x] . GAD[lora]));
*)
ggv = 1/(2*2(D+2))(SPD[x] DiracSigma[GAD[lora,lorb]] + I(FVD[x,lora]GSD[x] . GAD[lorb]-FVD[x,lorb]GSD[x] . GAD[lora]));
ggv2 = 1/(2(D+2))(SPD[x] DiracSigma[GAD[lora,lorb]] + I(FVD[x,lora]GSD[x] . GAD[lorb]-FVD[x,lorb]GSD[x] . GAD[lora]));



Fourquarkd8type2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmp,tmp1,tmp2,hv3,hv4,holdf=OptionValue[HoldFlavor],atr,diagrams,fdir,files,pall,dot},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[Fourquarkd8type2::inderr];
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



diagrams={xtype1,xtype2,xtype3,xtype4,xtype5,xtype6,xtype7,xtype8,xtype9,xtype10,xtype11,xtype12,xtype13,xtype14,xtype15,xtype16};
(*{type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16};*)


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
			files=("Fourquarkd8type2_"<>#)&/@files;
			
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



(*--- The following are generated by algorithm ---*)


(* ::Input::Initialization:: *)
xtype1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,cola,lora,lorb,str,tr,dot,contract},


dia1=cqgq[d1, d2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]]*tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]] + cqgq[c2, c1]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, ggv . SUNT[cola], v1, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]] + cqgq[d1, d2]*cqq[b1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1]]]*tr[str[dot[hv4, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[a2, b2]*cqq[d1, a1]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[v2, ggv2 . SUNT[cola], v1, xprop[x], hv4, xprop[-x]]]] - cqgq[a2, c1]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, v1, xprop[x], hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] - cqgq[d1, d2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[d1, b2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, ggv . SUNT[cola], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[c2, b2]*cqq[d1, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, a1]*tr[str[dot[hv3, xprop[-x], v2, ggv2 . SUNT[cola], v1, xprop[x], hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia1=(FourierXP[dia1,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(4-D) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,cola,lora,lorb,str,tr,dot,contract},


dia2=cqgq[d1, b2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, v1, xprop[x]]]]*tr[str[dot[hv3, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[b1, d2]*cqq[d1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2, ggv . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[a2, b2]*cqq[d1, c1]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[v2, ggv2 . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]] + cqgq[b1, b2]*cqq[d1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v2, ggv . SUNT[cola]]]] + (CF*cqgq[d1, c1]*cqq[a2, d2]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1]*qfact2[-1 + D]*SPD[x]*tr[str[dot[v2]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3]]])/2 - (CF*cqgq[b1, c1]*cqq[a2, d2]*FlavorDelta[a2, d2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*qfact2[-1 + D]*SPD[x]*tr[str[dot[hv4, xprop[-x], v1, v2, xprop[x], hv3]]])/2 - cqgq[c2, c1]*cqq[d1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2]]] - cqgq[a2, d2]*cqq[d1, a1]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1]*tr[str[dot[hv4, xprop[-x], v2, ggv2 . SUNT[cola], v1, xprop[x], hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia2=(FourierXP[dia2,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(4-D) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,cola,lora,lorb,str,tr,dot,contract},


dia3=cqgq[a2, c1]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, v1, xprop[x]]]]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] + cqgq[d1, b2]*cqq[b1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, xprop[-x], v1]]]*tr[str[dot[hv3, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[c2, c1]*cqq[b1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] - (CF*cqgq[d1, a1]*cqq[c2, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, a1]*qfact2[-1 + D]*SPD[x]*tr[str[dot[hv3, xprop[-x], v2, v1, xprop[x], hv4]]])/2 - cqgq[c2, c1]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, v1, xprop[x], hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] - cqgq[b1, d2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]], hv4, v2, xprop[x]]]] - cqgq[b1, b2]*cqq[d1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]], hv3, xprop[-x], v2]]] - cqgq[a2, c1]*cqq[d1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2]]];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia3=(FourierXP[dia3,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=QEvaluate[I ScaleMu^(4-D) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,cola,lora,lorb,str,tr,dot,contract},


dia4=-(cqgq[c2, d2]*cqq[d1, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, c1]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[v2, ggv2 . SUNT[cola]]]]*tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]) + cqgq[d1, d2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv4, ggv . SUNT[cola], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] + cqgq[b1, b2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, ggv . SUNT[cola], v2, xprop[x]]]] + cqgq[c2, a1]*cqq[d1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x]]]] - cqgq[a2, b2]*cqq[b1, a1]*FlavorDelta[a2, b2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[v1, ggv2 . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]] - cqgq[a2, a1]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x], hv4, v2, xprop[x]]]] - cqgq[c2, b2]*cqq[b1, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, ggv2 . SUNT[cola], v2, xprop[x], hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[a2, d2]*cqq[b1, c1]*FlavorDelta[a2, d2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, xprop[-x], v1, ggv2 . SUNT[cola], v2, xprop[x], hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia4=(FourierXP[dia4,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=QEvaluate[I ScaleMu^(4-D) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,cola,lora,lorb,str,tr,dot,contract},


dia5=cqgq[c2, c1]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, v1, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] + cqgq[d1, d2]*cqq[b1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, ggv . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] + cqgq[b1, b2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]] + (CF*cqgq[d1, c1]*cqq[c2, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, c1]*qfact2[-1 + D]*SPD[x]*tr[str[dot[v2]]]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4]]])/2 + cqgq[c2, b2]*cqq[d1, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, c1]*tr[str[dot[v2, ggv2 . SUNT[cola]]]]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[c2, c1]*cqq[b1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, ggv . SUNT[cola], hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]] - cqgq[c2, c1]*cqq[d1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, ggv . SUNT[cola]]]] - cqgq[c2, a1]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x], hv3, ggv . SUNT[cola], v2, xprop[x]]]];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia5=(FourierXP[dia5,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=QEvaluate[I ScaleMu^(4-D) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,cola,lora,lorb,str,tr,dot,contract},


dia6=cqgq[b1, d2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv4, ggv . SUNT[cola], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] + cqgq[c2, a1]*cqq[b1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, ggv . SUNT[cola]]]] + cqgq[b1, d2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[d1, d2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - (CF*cqgq[b1, a1]*cqq[a2, b2]*FlavorDelta[a2, b2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*qfact2[-1 + D]*SPD[x]*tr[str[dot[hv3]]]*tr[str[dot[v1]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]])/2 + cqgq[c2, d2]*cqq[b1, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv3, xprop[-x], v1, ggv2 . SUNT[cola], v2, xprop[x]]]] - cqgq[c2, a1]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x], hv3, v2, xprop[x]]]] - cqgq[a2, c1]*cqq[d1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, ggv . SUNT[cola]]]];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia6=(FourierXP[dia6,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia6=QEvaluate[I ScaleMu^(4-D) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,cola,lora,lorb,str,tr,dot,contract},


dia7=cqgq[a2, c1]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]]*tr[str[dot[hv4, ggv . SUNT[cola], v1, xprop[x]]]] + cqgq[c2, a1]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, v2, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x]]]] + cqgq[b1, d2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, ggv . SUNT[cola], v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] + cqgq[a2, a1]*cqq[b1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]] - cqgq[d1, b2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[d1, d2]*cqq[b1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, hv3, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[a2, a1]*cqq[b1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, ggv . SUNT[cola], hv4, xprop[-x], v2, xprop[x]]]] - cqgq[b1, b2]*cqq[d1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola], hv3, xprop[-x], v2, ggv . SUNT[cola]]]];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia7=(FourierXP[dia7,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia7=QEvaluate[I ScaleMu^(4-D) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,cola,lora,lorb,str,tr,dot,contract},


dia8=-1/2*(CF*cqgq[b1, a1]*cqq[c2, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2]*qfact2[-1 + D]*SPD[x]*tr[str[dot[hv4]]]*tr[str[dot[v1]]]*tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]) + cqgq[d1, b2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, ggv . SUNT[cola], v1, xprop[x]]]] + (CF*cqgq[d1, a1]*cqq[c2, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, a1]*qfact2[-1 + D]*SPD[x]*tr[str[dot[hv4]]]*tr[str[dot[v2, v1, xprop[x], hv3, xprop[-x]]]])/2 + cqgq[c2, b2]*cqq[b1, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2]*tr[str[dot[v1, ggv2 . SUNT[cola]]]]*tr[str[dot[hv3, xprop[-x], v2, xprop[x], hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] + cqgq[c2, d2]*cqq[d1, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, a1]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[v2, ggv2 . SUNT[cola], v1, xprop[x], hv3, xprop[-x]]]] - cqgq[b1, b2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]], hv3, xprop[-x], v2, xprop[x]]]] - cqgq[b1, b2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]], hv3, v2, xprop[x]]]] - cqgq[d1, b2]*cqq[b1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, ggv . SUNT[cola], hv4, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia8=(FourierXP[dia8,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia8=QEvaluate[I ScaleMu^(4-D) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,cola,lora,lorb,str,tr,dot,contract},


dia9=cqgq[a2, a1]*cqq[d1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v2]]]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x]]]] - cqgq[c2, d2]*cqq[b1, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[v1, ggv2 . SUNT[cola]]]]*tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]] - (CF*cqgq[d1, c1]*cqq[a2, b2]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1]*qfact2[-1 + D]*SPD[x]*tr[str[dot[hv3]]]*tr[str[dot[v2]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]])/2 + cqgq[a2, c1]*cqq[d1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, ggv . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]] + (CF*cqgq[b1, c1]*cqq[c2, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2]*qfact2[-1 + D]*SPD[x]*tr[str[dot[hv4]]]*tr[str[dot[hv3, xprop[-x], v1, v2, xprop[x]]]])/2 - (CF*cqgq[b1, c1]*cqq[c2, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2]*qfact2[-1 + D]*SPD[x]*tr[str[dot[hv3, xprop[-x], v1, v2, xprop[x], hv4]]])/2 - cqgq[d1, b2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, v1, xprop[x], hv4, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[b1, d2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola], hv4, ggv . SUNT[cola], v2, xprop[x]]]];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia9=(FourierXP[dia9,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia9=QEvaluate[I ScaleMu^(4-D) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,cola,lora,lorb,str,tr,dot,contract},


dia10=cqgq[b1, b2]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, v2, xprop[x]]]]*tr[str[dot[hv3, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[c2, c1]*cqq[d1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, ggv . SUNT[cola]]]] + (CF*cqgq[b1, a1]*cqq[a2, d2]*FlavorDelta[a2, d2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*qfact2[-1 + D]*SPD[x]*tr[str[dot[v1]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x], hv3]]])/2 - cqgq[a2, a1]*cqq[d1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x], hv4, xprop[-x], v2]]] - cqgq[c2, a1]*cqq[d1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x], hv3, xprop[-x], v2]]] - cqgq[b1, b2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, ggv . SUNT[cola], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola], hv3, xprop[-x], v2, xprop[x]]]] - cqgq[d1, d2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, ggv . SUNT[cola], v1, xprop[x], hv3, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[c2, a1]*cqq[d1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x], hv3, xprop[-x], v2, ggv . SUNT[cola]]]];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia10=(FourierXP[dia10,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia10=QEvaluate[I ScaleMu^(4-D) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,cola,lora,lorb,str,tr,dot,contract},


dia11=cqgq[a2, a1]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, v2, xprop[x]]]]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, xprop[x]]]] - (CF*cqgq[d1, c1]*cqq[c2, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, c1]*qfact2[-1 + D]*SPD[x]*tr[str[dot[hv4]]]*tr[str[dot[v2]]]*tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]])/2 + cqgq[b1, d2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]]*tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]] + cqgq[a2, c1]*cqq[b1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, ggv . SUNT[cola]]]] + cqgq[d1, b2]*cqq[b1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v1, ggv . SUNT[cola]]]] + cqgq[d1, b2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]] + (CF*cqgq[b1, c1]*cqq[a2, b2]*FlavorDelta[a2, b2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*qfact2[-1 + D]*SPD[x]*tr[str[dot[hv3]]]*tr[str[dot[hv4, xprop[-x], v1, v2, xprop[x]]]])/2 - cqgq[b1, d2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, ggv . SUNT[cola], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola], hv4, xprop[-x], v2, xprop[x]]]];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia11=(FourierXP[dia11,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia11=QEvaluate[I ScaleMu^(4-D) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,cola,lora,lorb,str,tr,dot,contract},


dia12=cqgq[a2, c1]*cqq[b1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, xprop[-x], v1]]]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] + cqgq[a2, a1]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x]]]]*tr[str[dot[hv4, ggv . SUNT[cola], v2, xprop[x]]]] + cqgq[c2, c1]*cqq[b1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, xprop[-x], v1, ggv . SUNT[cola]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]] + cqgq[d1, b2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, ggv . SUNT[cola], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]] + cqgq[d1, d2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, ggv . SUNT[cola], v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] + cqgq[a2, d2]*cqq[d1, c1]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1]*tr[str[dot[v2, ggv2 . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[b1, d2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]], hv4, xprop[-x], v2, xprop[x]]]] - cqgq[d1, d2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, ggv . SUNT[cola], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia12=(FourierXP[dia12,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia12=QEvaluate[I ScaleMu^(4-D) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13,cola,lora,lorb,str,tr,dot,contract},


dia13=cqgq[b1, b2]*cqq[d1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v2]]]*tr[str[dot[hv3, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + cqgq[c2, a1]*cqq[b1, d2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1]]]*tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]] + cqgq[b1, b2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, ggv . SUNT[cola], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]] + cqgq[a2, b2]*cqq[b1, c1]*FlavorDelta[a2, b2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v1, ggv2 . SUNT[cola], v2, xprop[x]]]] - cqgq[d1, b2]*cqq[b1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, hv4, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[d1, d2]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, v1, xprop[x], hv3, xprop[-x], v2, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] - cqgq[a2, c1]*cqq[b1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] - cqgq[d1, b2]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, ggv . SUNT[cola], v1, xprop[x], hv4, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]];


dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia13=(FourierXP[dia13,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia13=QEvaluate[I ScaleMu^(4-D) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14,cola,lora,lorb,str,tr,dot,contract},


dia14=cqgq[a2, c1]*cqq[d1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]] + cqgq[a2, a1]*cqq[b1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, ggv . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]] + cqgq[a2, d2]*cqq[b1, a1]*FlavorDelta[a2, d2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[v1, ggv2 . SUNT[cola]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x], hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]] - cqgq[c2, a1]*cqq[b1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, hv3, xprop[-x], v2, xprop[x]]]] - cqgq[a2, a1]*cqq[c2, c1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x], hv4, ggv . SUNT[cola], v2, xprop[x]]]] - cqgq[b1, d2]*cqq[d1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola], hv4, xprop[-x], v2, ggv . SUNT[cola]]]] - cqgq[c2, a1]*cqq[b1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, ggv . SUNT[cola], hv3, xprop[-x], v2, xprop[x]]]] - cqgq[d1, d2]*cqq[b1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, ggv . SUNT[cola], hv3, xprop[-x], v2, DiracSigma[GAD[lora, lorb]] . SUNT[cola]]]];


dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia14=(FourierXP[dia14,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia14=QEvaluate[I ScaleMu^(4-D) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15,cola,lora,lorb,str,tr,dot,contract},


dia15=cqgq[b1, d2]*cqq[d1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2]]]*tr[str[dot[hv4, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]]]]] + (CF*cqgq[d1, a1]*cqq[a2, b2]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1]*qfact2[-1 + D]*SPD[x]*tr[str[dot[hv3]]]*tr[str[dot[v2, v1, xprop[x], hv4, xprop[-x]]]])/2 - (CF*cqgq[d1, a1]*cqq[a2, d2]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1]*qfact2[-1 + D]*SPD[x]*tr[str[dot[hv4, xprop[-x], v2, v1, xprop[x], hv3]]])/2 - cqgq[a2, a1]*cqq[b1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v1, hv4, xprop[-x], v2, xprop[x]]]] - cqgq[c2, c1]*cqq[b1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2, xprop[x]]]] - cqgq[b1, d2]*cqq[d1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, (CF*ggv) . DiracSigma[GAD[lora, lorb]], hv4, xprop[-x], v2]]] - cqgq[c2, c1]*cqq[a2, a1]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, ggv . SUNT[cola], v1, xprop[x], hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]] - cqgq[b1, b2]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, DiracSigma[GAD[lora, lorb]] . SUNT[cola], hv3, ggv . SUNT[cola], v2, xprop[x]]]];


dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia15=(FourierXP[dia15,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia15=QEvaluate[I ScaleMu^(4-D) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16,cola,lora,lorb,str,tr,dot,contract},


dia16=cqgq[c2, c1]*cqq[d1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . (CF*ggv), v2]]]*tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]] + cqgq[c2, a1]*cqq[a2, c1]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, ggv . SUNT[cola], v2, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x]]]] + cqgq[c2, a1]*cqq[d1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v2, ggv . SUNT[cola]]]]*tr[str[dot[hv4, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x]]]] + cqgq[a2, a1]*cqq[d1, d2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, ggv . SUNT[cola]]]] + (CF*cqgq[b1, a1]*cqq[c2, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2]*qfact2[-1 + D]*SPD[x]*tr[str[dot[v1]]]*tr[str[dot[hv3, xprop[-x], v2, xprop[x], hv4]]])/2 - cqgq[a2, a1]*cqq[d1, b2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v1, xprop[x], hv4, xprop[-x], v2, ggv . SUNT[cola]]]] - cqgq[a2, c1]*cqq[c2, a1]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, ggv . SUNT[cola], v1, xprop[x], hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]] - cqgq[a2, c1]*cqq[b1, b2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, ggv . SUNT[cola], hv3, DiracSigma[GAD[lora, lorb]] . SUNT[cola], v2, xprop[x]]]];


dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia16=(FourierXP[dia16,{x,q}]//SUNSimplify)/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia16=QEvaluate[I ScaleMu^(4-D) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]






(* ::Input::Initialization:: *)
(*type1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,coa,lora,lorb},


dia1=atr[SUNTrace[hv3 . xprop[-x] . v2]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v1 . xprop[x]]] cqgq[c2, a1] cqq[d1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[v1 . ggv2 . SUNT[cola]]] atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x] . hv3 . qgq[lora, lorb, cola]]] cqgq[a2, d2] cqq[b1, a1] FlavorDelta[a2, d2] FlavorDelta[b1, a1] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . v1 . xprop[x] . hv4 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[c2, c1] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . qgq[lora, lorb, cola] . hv4 . xprop[-x] . v2 . ggv . SUNT[cola]]] cqgq[b1, d2] cqq[d1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . hv3 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[a2, c1] cqq[b1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]] atr[SUNTrace[hv4 . ggv . SUNT[cola] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + CF atr[SUNTrace[v1]] atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x] . hv4]] cqgq[b1, a1] cqq[c2, b2] FlavorDelta[a2, c1] FlavorDelta[b1, a1] FlavorDelta[c2, b2] FlavorDelta[d1, d2] qfact2[-1 + D] SPD[x] - CF atr[SUNTrace[hv3 . xprop[-x] . v1 . v2 . xprop[x] . hv4]] cqgq[b1, c1] cqq[c2, b2] FlavorDelta[a2, a1] FlavorDelta[b1, c1] FlavorDelta[c2, b2] FlavorDelta[d1, d2] qfact2[-1 + D] SPD[x];


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(4-D) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,coa,lora,lorb},


dia2=atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . ggv . SUNT[cola]]] cqgq[c2, a1] cqq[b1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . v1 . ggv . SUNT[cola] . qgq[lora, lorb, cola] . hv4 . xprop[-x] . v2 . xprop[x]]] cqgq[b1, d2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . hv4 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[c2, c1] cqq[b1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . v2 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] + atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . xprop[x]]] atr[SUNTrace[hv4 . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[a2, a1] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . xprop[-x] . v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[b1, b2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . v2 . xprop[x]]] atr[SUNTrace[hv3 . xprop[-x] . v1 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[b1, b2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] - CF atr[SUNTrace[hv4 . xprop[-x] . v1 . v2 . xprop[x] . hv3]] cqgq[b1, c1] cqq[a2, d2] FlavorDelta[a2, d2] FlavorDelta[b1, c1] FlavorDelta[c2, a1] FlavorDelta[d1, b2] qfact2[-1 + D] SPD[x];


dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(4-D) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,coa,lora,lorb},


dia3=atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v1]] cqgq[c2, a1] cqq[b1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv3 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v1 . ggv . SUNT[cola]]] cqgq[d1, b2] cqq[b1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . xprop[-x] . v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v2 . ggv . SUNT[cola]]] cqgq[b1, d2] cqq[d1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . ggv . SUNT[cola] . hv3 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[b1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . v1 . ggv . SUNT[cola] . qgq[lora, lorb, cola] . hv3 . xprop[-x] . v2 . xprop[x]]] cqgq[b1, b2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . v1 . xprop[x]]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[c2, c1] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v1]] cqgq[a2, a1] cqq[b1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] - CF atr[SUNTrace[hv3]] atr[SUNTrace[v1]] atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]] cqgq[b1, a1] cqq[a2, b2] FlavorDelta[a2, b2] FlavorDelta[b1, a1] FlavorDelta[c2, c1] FlavorDelta[d1, d2] qfact2[-1 + D] SPD[x];


dia3=FourierXP[dia3,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=QEvaluate[I ScaleMu^(4-D) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,coa,lora,lorb},


dia4=atr[SUNTrace[hv4 . xprop[-x] . v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[b1, d2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . ggv . SUNT[cola] . hv4 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[b1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] + atr[SUNTrace[v2 . ggv2 . SUNT[cola]]] atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . qgq[lora, lorb, cola]]] cqgq[a2, d2] cqq[d1, c1] FlavorDelta[a2, d2] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, c1] + atr[SUNTrace[v2 . ggv2 . SUNT[cola]]] atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . qgq[lora, lorb, cola]]] cqgq[c2, b2] cqq[d1, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, b2] FlavorDelta[d1, c1] - atr[SUNTrace[hv4 . xprop[-x] . v1 . ggv . SUNT[cola] . qgq[lora, lorb, cola] . hv3 . xprop[-x] . v2]] cqgq[b1, b2] cqq[d1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v2 . xprop[x]]] atr[SUNTrace[hv3 . ggv . SUNT[cola] . v1 . xprop[x]]] cqgq[c2, c1] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . v1 . xprop[x]]] atr[SUNTrace[hv4 . xprop[-x] . v2 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + CF atr[SUNTrace[hv4]] atr[SUNTrace[v2 . v1 . xprop[x] . hv3 . xprop[-x]]] cqgq[d1, a1] cqq[c2, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, d2] FlavorDelta[d1, a1] qfact2[-1 + D] SPD[x];


dia4=FourierXP[dia4,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=QEvaluate[I ScaleMu^(4-D) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,coa,lora,lorb},


dia5=atr[SUNTrace[hv4 . qgq[lora, lorb, cola]]] atr[SUNTrace[v2 . ggv2 . SUNT[cola] . v1 . xprop[x] . hv3 . xprop[-x]]] cqgq[c2, d2] cqq[d1, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, d2] FlavorDelta[d1, a1] + atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]] atr[SUNTrace[hv3 . v2 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]] atr[SUNTrace[hv4 . ggv . SUNT[cola] . v1 . qgq[lora, lorb, cola]]] cqgq[b1, d2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v2 . ggv . SUNT[cola]]] cqgq[a2, c1] cqq[d1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . ggv . SUNT[cola] . qgq[lora, lorb, cola] . hv4 . xprop[-x] . v2]] cqgq[b1, d2] cqq[d1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v1 . hv3 . xprop[-x] . v2 . xprop[x]]] cqgq[c2, a1] cqq[b1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . xprop[-x] . v2]] atr[SUNTrace[hv3 . xprop[-x] . v1 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[b1, b2] cqq[d1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] - CF atr[SUNTrace[hv3]] atr[SUNTrace[v2]] atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]] cqgq[d1, c1] cqq[a2, b2] FlavorDelta[a2, b2] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, c1] qfact2[-1 + D] SPD[x];


dia5=FourierXP[dia5,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=QEvaluate[I ScaleMu^(4-D) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,coa,lora,lorb},


dia6=atr[SUNTrace[hv4 . v1 . xprop[x]]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[a2, c1] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv3 . xprop[-x] . v2]] atr[SUNTrace[hv4 . xprop[-x] . v1 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[b1, d2] cqq[d1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . ggv . SUNT[cola] . v1 . xprop[x] . hv4 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[c2, c1] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . xprop[x] . hv4 . xprop[-x] . v2 . ggv . SUNT[cola]]] cqgq[a2, a1] cqq[d1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola]]] atr[SUNTrace[v2 . ggv2 . SUNT[cola]]] atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]] cqgq[a2, b2] cqq[d1, c1] FlavorDelta[a2, b2] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, c1] - atr[SUNTrace[hv3 . xprop[-x] . v1 . ggv2 . SUNT[cola] . v2 . xprop[x] . hv4 . qgq[lora, lorb, cola]]] cqgq[c2, b2] cqq[b1, c1] FlavorDelta[a2, a1] FlavorDelta[b1, c1] FlavorDelta[c2, b2] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . xprop[-x] . v1]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[c2, c1] cqq[b1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . v2 . xprop[x]]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v1 . xprop[x]]] cqgq[a2, a1] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia6=FourierXP[dia6,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia6=QEvaluate[I ScaleMu^(4-D) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,coa,lora,lorb},


dia7=atr[SUNTrace[hv3 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . ggv . SUNT[cola] . v1 . xprop[x]]] cqgq[d1, b2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . ggv . SUNT[cola] . v1 . qgq[lora, lorb, cola] . hv4 . xprop[-x] . v2 . xprop[x]]] cqgq[b1, d2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . ggv . SUNT[cola] . hv4 . xprop[-x] . v2 . xprop[x]]] cqgq[a2, a1] cqq[b1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . qgq[lora, lorb, cola] . hv4 . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[b1, d2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . v2 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v1 . xprop[x] . hv3 . xprop[-x] . v2]] cqgq[c2, a1] cqq[d1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]] atr[SUNTrace[hv4 . v2 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] - CF atr[SUNTrace[hv4 . xprop[-x] . v2 . v1 . xprop[x] . hv3]] cqgq[d1, a1] cqq[a2, d2] FlavorDelta[a2, d2] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, a1] qfact2[-1 + D] SPD[x];


dia7=FourierXP[dia7,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia7=QEvaluate[I ScaleMu^(4-D) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,coa,lora,lorb},


dia8=atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]] atr[SUNTrace[hv3 . ggv . SUNT[cola] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . qgq[lora, lorb, cola] . v2 . ggv . SUNT[cola]]] cqgq[c2, c1] cqq[d1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . ggv . SUNT[cola] . v1 . xprop[x] . hv3 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[a2, c1] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]] atr[SUNTrace[hv3 . v1 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[b1, b2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . xprop[x]]] atr[SUNTrace[hv4 . xprop[-x] . v2 . ggv . SUNT[cola]]] cqgq[a2, a1] cqq[d1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . xprop[-x] . v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v2 . ggv . SUNT[cola]]] cqgq[b1, b2] cqq[d1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] - CF atr[SUNTrace[hv4]] atr[SUNTrace[v1]] atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]] cqgq[b1, a1] cqq[c2, d2] FlavorDelta[a2, c1] FlavorDelta[b1, a1] FlavorDelta[c2, d2] FlavorDelta[d1, b2] qfact2[-1 + D] SPD[x] + CF atr[SUNTrace[v2]] atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4]] cqgq[d1, c1] cqq[c2, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, b2] FlavorDelta[d1, c1] qfact2[-1 + D] SPD[x];


dia8=FourierXP[dia8,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia8=QEvaluate[I ScaleMu^(4-D) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,coa,lora,lorb},


dia9=atr[SUNTrace[hv3 . qgq[lora, lorb, cola]]] atr[SUNTrace[v2 . ggv2 . SUNT[cola] . v1 . xprop[x] . hv4 . xprop[-x]]] cqgq[a2, b2] cqq[d1, a1] FlavorDelta[a2, b2] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, a1] + atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v2 . xprop[x]]] atr[SUNTrace[hv4 . ggv . SUNT[cola] . v1 . xprop[x]]] cqgq[a2, c1] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v1 . hv4 . xprop[-x] . v2 . xprop[x]]] cqgq[a2, a1] cqq[b1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v1 . xprop[x] . hv4 . xprop[-x] . v2]] cqgq[a2, a1] cqq[d1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v1 . ggv2 . SUNT[cola] . v2 . xprop[x]]] cqgq[c2, d2] cqq[b1, c1] FlavorDelta[a2, a1] FlavorDelta[b1, c1] FlavorDelta[c2, d2] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . ggv . SUNT[cola] . v1 . qgq[lora, lorb, cola] . hv3 . xprop[-x] . v2 . xprop[x]]] cqgq[b1, b2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v2 . xprop[x]]] atr[SUNTrace[hv3 . xprop[-x] . v1 . ggv . SUNT[cola]]] cqgq[c2, c1] cqq[b1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v2 . ggv . SUNT[cola]]] cqgq[c2, c1] cqq[d1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia9=FourierXP[dia9,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia9=QEvaluate[I ScaleMu^(4-D) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,coa,lora,lorb},


dia10=-(atr[SUNTrace[hv4 . xprop[-x] . v1 . ggv2 . SUNT[cola] . v2 . xprop[x] . hv3 . qgq[lora, lorb, cola]]] cqgq[a2, d2] cqq[b1, c1] FlavorDelta[a2, d2] FlavorDelta[b1, c1] FlavorDelta[c2, a1] FlavorDelta[d1, b2]) + atr[SUNTrace[hv4 . v1 . xprop[x]]] atr[SUNTrace[hv3 . xprop[-x] . v2 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . xprop[x]]] atr[SUNTrace[hv3 . xprop[-x] . v2 . ggv . SUNT[cola]]] cqgq[c2, a1] cqq[d1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . ggv . SUNT[cola] . v1 . xprop[x] . hv4 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . ggv . SUNT[cola] . qgq[lora, lorb, cola] . hv4 . v2 . xprop[x]]] cqgq[b1, d2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . ggv . SUNT[cola] . hv3 . xprop[-x] . v2 . xprop[x]]] cqgq[c2, a1] cqq[b1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . qgq[lora, lorb, cola] . v2 . ggv . SUNT[cola]]] cqgq[a2, c1] cqq[d1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . xprop[x] . hv3 . xprop[-x] . v2 . ggv . SUNT[cola]]] cqgq[c2, a1] cqq[d1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2];


dia10=FourierXP[dia10,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia10=QEvaluate[I ScaleMu^(4-D) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,coa,lora,lorb},


dia11=-(atr[SUNTrace[hv3 . xprop[-x] . v1 . ggv . SUNT[cola] . hv4 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[c2, c1] cqq[b1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2]) - atr[SUNTrace[hv3 . xprop[-x] . v1 . hv4 . xprop[-x] . v2 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[b1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . xprop[x] . hv4 . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[a2, a1] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . hv3 . xprop[-x] . v2 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[b1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v2]] cqgq[a2, c1] cqq[d1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[v1 . ggv2 . SUNT[cola]]] atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x] . hv4 . qgq[lora, lorb, cola]]] cqgq[c2, b2] cqq[b1, a1] FlavorDelta[a2, c1] FlavorDelta[b1, a1] FlavorDelta[c2, b2] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . xprop[-x] . v2]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v1 . xprop[x]]] cqgq[a2, a1] cqq[d1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + CF atr[SUNTrace[v1]] atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x] . hv3]] cqgq[b1, a1] cqq[a2, d2] FlavorDelta[a2, d2] FlavorDelta[b1, a1] FlavorDelta[c2, c1] FlavorDelta[d1, b2] qfact2[-1 + D] SPD[x];


dia11=FourierXP[dia11,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia11=QEvaluate[I ScaleMu^(4-D) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,coa,lora,lorb},


dia12=-(atr[SUNTrace[hv4 . xprop[-x] . v2 . ggv2 . SUNT[cola] . v1 . xprop[x] . hv3 . qgq[lora, lorb, cola]]] cqgq[a2, d2] cqq[d1, a1] FlavorDelta[a2, d2] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, a1]) + atr[SUNTrace[hv4 . xprop[-x] . v1]] atr[SUNTrace[hv3 . xprop[-x] . v2 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[b1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v1 . xprop[x] . hv4 . v2 . xprop[x]]] cqgq[a2, a1] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola]]] atr[SUNTrace[v2 . ggv2 . SUNT[cola]]] atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]] cqgq[c2, d2] cqq[d1, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, d2] FlavorDelta[d1, c1] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola]]] atr[SUNTrace[v1 . ggv2 . SUNT[cola]]] atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]] cqgq[a2, b2] cqq[b1, a1] FlavorDelta[a2, b2] FlavorDelta[b1, a1] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]] atr[SUNTrace[hv3 . ggv . SUNT[cola] . v1 . qgq[lora, lorb, cola]]] cqgq[b1, b2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v1 . ggv . SUNT[cola]]] cqgq[d1, d2] cqq[b1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + CF atr[SUNTrace[hv3]] atr[SUNTrace[v2 . v1 . xprop[x] . hv4 . xprop[-x]]] cqgq[d1, a1] cqq[a2, b2] FlavorDelta[a2, b2] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, a1] qfact2[-1 + D] SPD[x];


dia12=FourierXP[dia12,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia12=QEvaluate[I ScaleMu^(4-D) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13,coa,lora,lorb},


dia13=atr[SUNTrace[hv3 . v2 . xprop[x]]] atr[SUNTrace[hv4 . xprop[-x] . v1 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[b1, d2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . xprop[x]]] atr[SUNTrace[hv3 . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[c2, a1] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . xprop[-x] . v1]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[a2, c1] cqq[b1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . ggv . SUNT[cola] . qgq[lora, lorb, cola] . hv3 . v2 . xprop[x]]] cqgq[b1, b2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . ggv . SUNT[cola] . hv3 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[a2, c1] cqq[b1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . v1 . xprop[x] . hv3 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[a2, c1] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - CF atr[SUNTrace[hv3 . xprop[-x] . v2 . v1 . xprop[x] . hv4]] cqgq[d1, a1] cqq[c2, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, b2] FlavorDelta[d1, a1] qfact2[-1 + D] SPD[x] + CF atr[SUNTrace[hv3]] atr[SUNTrace[hv4 . xprop[-x] . v1 . v2 . xprop[x]]] cqgq[b1, c1] cqq[a2, b2] FlavorDelta[a2, b2] FlavorDelta[b1, c1] FlavorDelta[c2, a1] FlavorDelta[d1, d2] qfact2[-1 + D] SPD[x];


dia13=FourierXP[dia13,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia13=QEvaluate[I ScaleMu^(4-D) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14,coa,lora,lorb},


dia14=-(atr[SUNTrace[hv3 . xprop[-x] . v2 . ggv2 . SUNT[cola] . v1 . xprop[x] . hv4 . qgq[lora, lorb, cola]]] cqgq[c2, b2] cqq[d1, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, b2] FlavorDelta[d1, a1]) + atr[SUNTrace[hv3 . v2 . xprop[x]]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v1 . xprop[x]]] cqgq[c2, a1] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]] atr[SUNTrace[hv4 . v1 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[b1, d2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v2]] cqgq[a2, c1] cqq[d1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . v1 . xprop[x] . hv4 . xprop[-x] . v2 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . ggv . SUNT[cola] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . ggv . SUNT[cola] . v1 . xprop[x]]] cqgq[d1, d2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] - CF atr[SUNTrace[hv4]] atr[SUNTrace[v2]] atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]] cqgq[d1, c1] cqq[c2, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, d2] FlavorDelta[d1, c1] qfact2[-1 + D] SPD[x];


dia14=FourierXP[dia14,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia14=QEvaluate[I ScaleMu^(4-D) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15,coa,lora,lorb},


dia15=atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v2 . xprop[x]]] atr[SUNTrace[hv4 . xprop[-x] . v1 . ggv . SUNT[cola]]] cqgq[a2, c1] cqq[b1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . ggv . SUNT[cola] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v2]] cqgq[c2, c1] cqq[d1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v1 . xprop[x] . hv3 . v2 . xprop[x]]] cqgq[c2, a1] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . xprop[x] . hv3 . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[c2, a1] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v1 . ggv2 . SUNT[cola] . v2 . xprop[x]]] cqgq[a2, b2] cqq[b1, c1] FlavorDelta[a2, b2] FlavorDelta[b1, c1] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + CF atr[SUNTrace[hv4]] atr[SUNTrace[hv3 . xprop[-x] . v1 . v2 . xprop[x]]] cqgq[b1, c1] cqq[c2, d2] FlavorDelta[a2, a1] FlavorDelta[b1, c1] FlavorDelta[c2, d2] FlavorDelta[d1, b2] qfact2[-1 + D] SPD[x] + CF atr[SUNTrace[v2]] atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3]] cqgq[d1, c1] cqq[a2, d2] FlavorDelta[a2, d2] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, c1] qfact2[-1 + D] SPD[x];


dia15=FourierXP[dia15,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia15=QEvaluate[I ScaleMu^(4-D) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16,coa,lora,lorb},


dia16=-(atr[SUNTrace[hv4 . qgq[lora, lorb, cola]]] atr[SUNTrace[v1 . ggv2 . SUNT[cola]]] atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]] cqgq[c2, d2] cqq[b1, a1] FlavorDelta[a2, c1] FlavorDelta[b1, a1] FlavorDelta[c2, d2] FlavorDelta[d1, b2]) - atr[SUNTrace[hv4 . xprop[-x] . v1 . qgq[lora, lorb, cola] . hv3 . ggv . SUNT[cola] . v2 . xprop[x]]] cqgq[b1, b2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . ggv . SUNT[cola] . v1 . xprop[x] . hv3 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . v1 . xprop[x] . hv3 . xprop[-x] . v2 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . qgq[lora, lorb, cola] . hv3 . xprop[-x] . v2 . ggv . SUNT[cola]]] cqgq[b1, b2] cqq[d1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . ggv . SUNT[cola]]] cqgq[a2, a1] cqq[b1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . xprop[-x] . v1]] atr[SUNTrace[hv4 . xprop[-x] . v2 . ggv . SUNT[cola] . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[b1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . ggv . SUNT[cola] . v2]] cqgq[c2, c1] cqq[d1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia16=FourierXP[dia16,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia16=QEvaluate[I ScaleMu^(4-D) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


End[]
(*EndPackage[]*)
