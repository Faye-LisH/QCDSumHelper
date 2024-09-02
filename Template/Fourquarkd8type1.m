(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
Fourquarkd8type1::usage="Fourquarkd8type1[q_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}}] give the \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)Gq\[RightAngleBracket]\[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)q\[RightAngleBracket] contribution of OPE result about <J1 J2>, with J1=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a1\)]\)\!\(\*SubscriptBox[\(v1\[CapitalPsi]\), \(b1\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c1\)]\)\!\(\*SubscriptBox[\(v2\[CapitalPsi]\), \(d1\)]\) and J2=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(b2\)]\)\!\(\*SubscriptBox[\(v3\[CapitalPsi]\), \(a2\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(d2\)]\)\!\(\*SubscriptBox[\(v4\[CapitalPsi]\), \(c2\)]\), here the background gluon comes from quark propatator. "

Fourquarkd8type1::inderr="Dummy indices conflict!"

Fourquarkd8type1::curerr="Unknow current structure!"



(* ::Code::Initialization::Plain:: *)
Begin["`Private`Fourquarkd8type1`"]


(* ::Code::Initialization::Plain:: *)
Options[Fourquarkd8type1] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True
}



(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(* propagators *)

prop[q_]= I GSD[q]FAD[q];


(* ::Code::Initialization::Plain:: *)
(* xprop[x_]=FourierPX[prop[q],{q,x}]; *)
xprop[x_]= 1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];

(*  quark propagator with background gluon  *)
(*  ---<--- q *)

(* x ---<--- 0 *)
(*xprog[x_,lora_,lorb_,cola_]=FourierPX[prog[q,lora,lorb,cola],{q,x}];*)
xprog[x_,lora_,lorb_,cola_]= (GAD[lora] . GAD[lorb] . GSD[x] qfact1[1/32 I^(-D) \[Pi]^(-D/2) qGamma[-1+D/2]] SPD[x,x]^(1-D/2) SUNT[cola]
							+GSD[x] . GAD[lorb] . GAD[lora] qfact1[qfact2[-(1/32) I^(-D) \[Pi]^(-D/2)] qGamma[-1+D/2]] SPD[x,x]^(1-D/2) SUNT[cola]
							+2 FVD[x,lora] FVD[x,lorb] GSD[x] qfact1[-(1/16) I^(-D) \[Pi]^(-D/2) qGamma[D/2]] SPD[x,x]^(-D/2) SUNT[cola]);
							


(* ::Code::Initialization::Plain:: *)
(*--- Condensate ---*)

qgq[lora_,lorb_,cola_]=  SUNT[cola] . DiracSigma[GAD[lora,lorb]];

cqgq[a_,b_]:=-1/(CA CF 4(D-1)D)Condensate[{a,"G",b}]
cqq[f1_,f2_]:=-1/(4CA)Condensate[{f1,f2}];		

(*-------------------------------------------------*)


Fourquarkd8type1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmp,tmp1,tmp2,hv3,hv4,holdf=OptionValue[HoldFlavor],atr,diagrams,fdir,pall,files,dot},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[Fourquarkd8type1::inderr];
	Abort[]
];


If[OptionValue[AutoNDR]===True,
	atr=TR5
	,
	atr=TR
];


(* B A^+ B *)
hv3=v3//ComplexConjugate;
hv4=v4//ComplexConjugate//FCI;
(* ComplexConjugate[sigma] will expand the DiracSigma to gamma matrices, recover it to DiracSigma *)
hv4=(hv4/.Dot->dot)/.(f1_ dot[aa___,DiracGamma[LorentzIndex[lor1_,dim___],dim___],DiracGamma[LorentzIndex[lor2_,dim___],dim___],bb___]+f2_ dot[aa___,DiracGamma[LorentzIndex[lor2_,dim___],dim___],DiracGamma[LorentzIndex[lor1_,dim___],dim___],bb___])/;(f1+f2===0):>-2I f1 dot[aa,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],bb];
hv4=hv4/. {null->1,dot->Dot};


diagrams={type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16};

(*--------------------------------------------*)		

					
If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];	


(*---------------------------------------------------*)
If[pall===True,

	DistributeDefinitions[qq,v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
	
	tmp=Plus@@WaitAll[ParallelSubmit[{hv3,hv4,holdf,atr},#[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr]]&/@diagrams];
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	
,

	If[pall==="External",
	
		DistributeDefinitions[qq,v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
		
		If[fdir==="None",
		
			ParallelSubmit[{hv3,hv4,holdf,atr},#[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr]]&/@diagrams
		,
		
			files=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files=("Fourquarkd8type1_"<>#)&/@files;
			
			ImExport[fdir,
						files,
						{{hv3,hv4,holdf,atr},
						{qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr},
						diagrams}
						]
		]
		(* !!! do not touch tmp before WaitAll[tmp], any function may modifiy the things in EvaluationObject[...] should not be used here. !!! *)
	,
	
	
		tmp=Plus@@Map[#[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr]&,diagrams];
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	]
]

]



(* ::Code::Initialization::Plain:: *)
(*--- The following are generated by algorithm ---*)


(* ::Input::Initialization:: *)
type1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,cola,lora,lorb},


dia1=-(atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v2 . v1 . xprop[x] . hv3 . qgq[lora, lorb, cola]]] cqgq[a2, d2] cqq[d1, a1] FlavorDelta[a2, d2] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, a1]) + atr[SUNTrace[hv3 . v2 . xprop[x]]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . qgq[lora, lorb, cola]]] cqgq[b1, d2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1]] atr[SUNTrace[hv3 . xprop[-x] . v2 . xprog[x, lora, lorb, cola]]] cqgq[c2, a1] cqq[b1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . xprop[-x] . v1]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[b1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v2 . xprop[x] . hv3]] cqgq[b1, a1] cqq[a2, d2] FlavorDelta[a2, d2] FlavorDelta[b1, a1] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . hv4 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[c2, c1] cqq[b1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[v2]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v1 . xprog[x, lora, lorb, cola]]] cqgq[c2, d2] cqq[d1, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, d2] FlavorDelta[d1, c1] - atr[SUNTrace[hv4 . xprop[-x] . v1 . hv3 . xprog[-x, lora, lorb, cola] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[b1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3]] atr[SUNTrace[hv4 . xprop[-x] . v1 . qgq[lora, lorb, cola] . v2 . xprog[x, lora, lorb, cola]]] cqgq[b1, c1] cqq[a2, b2] FlavorDelta[a2, b2] FlavorDelta[b1, c1] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . xprop[-x] . v2]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . xprog[x, lora, lorb, cola]]] cqgq[a2, a1] cqq[d1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia1=FourierXP[dia1,{x,q}]/.SUNN->CA;

dia1=QEvaluate[I ScaleMu^(4-D) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,cola,lora,lorb},


dia2=-(atr[SUNTrace[hv4 . xprop[-x] . v2 . qgq[lora, lorb, cola] . v1 . xprog[x, lora, lorb, cola] . hv3]] cqgq[d1, a1] cqq[a2, d2] FlavorDelta[a2, d2] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, a1]) + atr[SUNTrace[hv3 . qgq[lora, lorb, cola]]] atr[SUNTrace[v2 . v1 . xprog[x, lora, lorb, cola] . hv4 . xprop[-x]]] cqgq[a2, b2] cqq[d1, a1] FlavorDelta[a2, b2] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, a1] + atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v2]] atr[SUNTrace[hv4 . xprop[-x] . v1 . qgq[lora, lorb, cola]]] cqgq[b1, d2] cqq[d1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[v1]] atr[SUNTrace[hv4 . xprop[-x] . v2 . xprog[x, lora, lorb, cola] . hv3 . qgq[lora, lorb, cola]]] cqgq[a2, d2] cqq[b1, a1] FlavorDelta[a2, d2] FlavorDelta[b1, a1] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . hv4 . qgq[lora, lorb, cola] . v2 . xprog[x, lora, lorb, cola]]] cqgq[c2, c1] cqq[b1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] + atr[SUNTrace[v2]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . xprop[x] . hv3 . qgq[lora, lorb, cola]]] cqgq[a2, d2] cqq[d1, c1] FlavorDelta[a2, d2] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, c1] - atr[SUNTrace[hv3]] atr[SUNTrace[v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v1 . xprog[x, lora, lorb, cola]]] cqgq[d1, c1] cqq[a2, b2] FlavorDelta[a2, b2] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, c1] - atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . v2 . xprop[x] . hv4 . qgq[lora, lorb, cola]]] cqgq[c2, b2] cqq[b1, c1] FlavorDelta[a2, a1] FlavorDelta[b1, c1] FlavorDelta[c2, b2] FlavorDelta[d1, d2] - atr[SUNTrace[hv3]] atr[SUNTrace[v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v2 . xprog[x, lora, lorb, cola]]] cqgq[b1, a1] cqq[a2, b2] FlavorDelta[a2, b2] FlavorDelta[b1, a1] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . xprop[-x] . v1]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v2 . xprog[x, lora, lorb, cola]]] cqgq[c2, c1] cqq[b1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia2=FourierXP[dia2,{x,q}]/.SUNN->CA;

dia2=QEvaluate[I ScaleMu^(4-D) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,cola,lora,lorb},


dia3=atr[SUNTrace[hv3 . qgq[lora, lorb, cola]]] atr[SUNTrace[v2 . v1 . xprop[x] . hv4 . xprog[-x, lora, lorb, cola]]] cqgq[a2, b2] cqq[d1, a1] FlavorDelta[a2, b2] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, a1] + atr[SUNTrace[hv4]] atr[SUNTrace[v2 . qgq[lora, lorb, cola] . v1 . xprog[x, lora, lorb, cola] . hv3 . xprop[-x]]] cqgq[d1, a1] cqq[c2, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, d2] FlavorDelta[d1, a1] - atr[SUNTrace[hv3 . xprop[-x] . v1 . xprog[x, lora, lorb, cola] . hv4 . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4]] atr[SUNTrace[v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v2 . xprop[x]]] cqgq[b1, a1] cqq[c2, d2] FlavorDelta[a2, c1] FlavorDelta[b1, a1] FlavorDelta[c2, d2] FlavorDelta[d1, b2] + atr[SUNTrace[v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . xprop[x] . hv3]] cqgq[d1, c1] cqq[a2, d2] FlavorDelta[a2, d2] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, c1] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . xprog[x, lora, lorb, cola] . hv3 . v2 . xprop[x]]] cqgq[c2, a1] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . hv3 . xprop[-x] . v2 . xprog[x, lora, lorb, cola]]] cqgq[c2, a1] cqq[b1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . hv3 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[b1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . qgq[lora, lorb, cola] . hv3 . xprog[-x, lora, lorb, cola] . v2]] cqgq[b1, b2] cqq[d1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v2 . xprog[x, lora, lorb, cola]]] cqgq[b1, b2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia3=FourierXP[dia3,{x,q}]/.SUNN->CA;

dia3=QEvaluate[I ScaleMu^(4-D) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,cola,lora,lorb},


dia4=atr[SUNTrace[hv3 . v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . xprop[x]]] cqgq[d1, b2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v2]] atr[SUNTrace[hv4 . xprop[-x] . v1 . xprog[x, lora, lorb, cola]]] cqgq[a2, c1] cqq[d1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . hv4 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[b1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . xprog[x, lora, lorb, cola] . hv4 . v2 . xprop[x]]] cqgq[a2, a1] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . qgq[lora, lorb, cola] . hv4 . xprog[-x, lora, lorb, cola] . v2]] cqgq[b1, d2] cqq[d1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] + atr[SUNTrace[v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v1 . xprog[x, lora, lorb, cola] . hv3]] cqgq[d1, c1] cqq[a2, d2] FlavorDelta[a2, d2] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, c1] + atr[SUNTrace[v2]] atr[SUNTrace[hv3 . xprop[-x] . v1 . xprog[x, lora, lorb, cola] . hv4 . qgq[lora, lorb, cola]]] cqgq[c2, b2] cqq[d1, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, b2] FlavorDelta[d1, c1] + atr[SUNTrace[v1]] atr[SUNTrace[hv3 . xprop[-x] . v2 . xprog[x, lora, lorb, cola] . hv4 . qgq[lora, lorb, cola]]] cqgq[c2, b2] cqq[b1, a1] FlavorDelta[a2, c1] FlavorDelta[b1, a1] FlavorDelta[c2, b2] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . v1 . xprop[x]]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v1 . xprog[x, lora, lorb, cola]]] cqgq[d1, d2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia4=FourierXP[dia4,{x,q}]/.SUNN->CA;

dia4=QEvaluate[I ScaleMu^(4-D) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,cola,lora,lorb},


dia5=-(atr[SUNTrace[hv3 . xprop[-x] . v2 . qgq[lora, lorb, cola] . v1 . xprog[x, lora, lorb, cola] . hv4]] cqgq[d1, a1] cqq[c2, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, b2] FlavorDelta[d1, a1]) + atr[SUNTrace[hv4]] atr[SUNTrace[v2 . qgq[lora, lorb, cola] . v1 . xprop[x] . hv3 . xprog[-x, lora, lorb, cola]]] cqgq[d1, a1] cqq[c2, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, d2] FlavorDelta[d1, a1] + atr[SUNTrace[hv3 . v2 . xprog[x, lora, lorb, cola]]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . xprop[x]]] cqgq[c2, a1] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . v1 . xprop[x] . hv4 . qgq[lora, lorb, cola] . v2 . xprog[x, lora, lorb, cola]]] cqgq[c2, c1] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . v1 . xprog[x, lora, lorb, cola] . hv4 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . xprop[x] . hv4 . qgq[lora, lorb, cola] . v2]] cqgq[c2, c1] cqq[d1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[v1]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v2 . xprog[x, lora, lorb, cola]]] cqgq[c2, d2] cqq[b1, a1] FlavorDelta[a2, c1] FlavorDelta[b1, a1] FlavorDelta[c2, d2] FlavorDelta[d1, b2] + atr[SUNTrace[v2]] atr[SUNTrace[hv4 . xprop[-x] . v1 . xprog[x, lora, lorb, cola] . hv3 . qgq[lora, lorb, cola]]] cqgq[a2, d2] cqq[d1, c1] FlavorDelta[a2, d2] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, c1] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . xprop[x] . hv3 . v2 . xprog[x, lora, lorb, cola]]] cqgq[c2, a1] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[b1, c1] cqq[a2, b2] FlavorDelta[a2, b2] FlavorDelta[b1, c1] FlavorDelta[c2, a1] FlavorDelta[d1, d2];


dia5=FourierXP[dia5,{x,q}]/.SUNN->CA;

dia5=QEvaluate[I ScaleMu^(4-D) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,cola,lora,lorb},


dia6=atr[SUNTrace[hv3]] atr[SUNTrace[v2 . qgq[lora, lorb, cola] . v1 . xprog[x, lora, lorb, cola] . hv4 . xprop[-x]]] cqgq[d1, a1] cqq[a2, b2] FlavorDelta[a2, b2] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, a1] + atr[SUNTrace[hv3]] atr[SUNTrace[v2 . qgq[lora, lorb, cola] . v1 . xprop[x] . hv4 . xprog[-x, lora, lorb, cola]]] cqgq[d1, a1] cqq[a2, b2] FlavorDelta[a2, b2] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, a1] + atr[SUNTrace[hv4 . qgq[lora, lorb, cola]]] atr[SUNTrace[v2 . v1 . xprog[x, lora, lorb, cola] . hv3 . xprop[-x]]] cqgq[c2, d2] cqq[d1, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, d2] FlavorDelta[d1, a1] + atr[SUNTrace[hv4 . qgq[lora, lorb, cola]]] atr[SUNTrace[v2 . v1 . xprop[x] . hv3 . xprog[-x, lora, lorb, cola]]] cqgq[c2, d2] cqq[d1, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, d2] FlavorDelta[d1, a1] - atr[SUNTrace[hv4 . xprop[-x] . v1 . qgq[lora, lorb, cola] . v2 . xprog[x, lora, lorb, cola] . hv3]] cqgq[b1, c1] cqq[a2, d2] FlavorDelta[a2, d2] FlavorDelta[b1, c1] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . hv4 . xprog[-x, lora, lorb, cola] . v2 . xprop[x]]] cqgq[a2, a1] cqq[b1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] + atr[SUNTrace[v2]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . xprop[x] . hv4 . qgq[lora, lorb, cola]]] cqgq[c2, b2] cqq[d1, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, b2] FlavorDelta[d1, c1] - atr[SUNTrace[hv4]] atr[SUNTrace[v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . xprop[x]]] cqgq[d1, c1] cqq[c2, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, d2] FlavorDelta[d1, c1] - atr[SUNTrace[hv4]] atr[SUNTrace[v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v1 . xprog[x, lora, lorb, cola]]] cqgq[d1, c1] cqq[c2, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, d2] FlavorDelta[d1, c1] - atr[SUNTrace[hv4 . v1 . xprog[x, lora, lorb, cola] . hv3 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[a2, c1] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2];


dia6=FourierXP[dia6,{x,q}]/.SUNN->CA;

dia6=QEvaluate[I ScaleMu^(4-D) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,cola,lora,lorb},


dia7=-(atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v2 . qgq[lora, lorb, cola] . v1 . xprop[x] . hv4]] cqgq[d1, a1] cqq[c2, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, b2] FlavorDelta[d1, a1]) + atr[SUNTrace[hv4 . v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v2 . xprop[x]]] cqgq[b1, d2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . v1 . xprop[x]]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . v1 . qgq[lora, lorb, cola] . hv4 . xprop[-x] . v2 . xprog[x, lora, lorb, cola]]] cqgq[b1, d2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . xprop[x] . hv4 . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . hv3 . xprog[-x, lora, lorb, cola] . v2 . xprop[x]]] cqgq[c2, a1] cqq[b1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . v1 . xprop[x] . hv3 . xprog[-x, lora, lorb, cola] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . xprop[-x] . v1]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[b1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . xprop[x]]] cqgq[d1, d2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v2]] atr[SUNTrace[hv3 . xprop[-x] . v1 . qgq[lora, lorb, cola]]] cqgq[b1, b2] cqq[d1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia7=FourierXP[dia7,{x,q}]/.SUNN->CA;

dia7=QEvaluate[I ScaleMu^(4-D) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,cola,lora,lorb},


dia8=atr[SUNTrace[v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v2 . xprog[x, lora, lorb, cola] . hv3]] cqgq[b1, a1] cqq[a2, d2] FlavorDelta[a2, d2] FlavorDelta[b1, a1] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . xprop[x] . hv4 . xprog[-x, lora, lorb, cola] . v2]] cqgq[a2, a1] cqq[d1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . xprog[x, lora, lorb, cola] . hv4 . qgq[lora, lorb, cola] . v2]] cqgq[c2, c1] cqq[d1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . v1 . qgq[lora, lorb, cola] . hv3 . xprop[-x] . v2 . xprog[x, lora, lorb, cola]]] cqgq[b1, b2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . xprop[x] . hv3 . qgq[lora, lorb, cola] . v2]] cqgq[a2, c1] cqq[d1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . xprop[x] . hv3 . xprog[-x, lora, lorb, cola] . v2]] cqgq[c2, a1] cqq[d1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . v2 . xprop[x]]] cqgq[a2, b2] cqq[b1, c1] FlavorDelta[a2, b2] FlavorDelta[b1, c1] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v2 . xprog[x, lora, lorb, cola] . hv4]] cqgq[b1, a1] cqq[c2, b2] FlavorDelta[a2, c1] FlavorDelta[b1, a1] FlavorDelta[c2, b2] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . v1 . xprog[x, lora, lorb, cola]]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[c2, c1] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v2 . xprop[x]]] cqgq[a2, a1] cqq[b1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia8=FourierXP[dia8,{x,q}]/.SUNN->CA;

dia8=QEvaluate[I ScaleMu^(4-D) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,cola,lora,lorb},


dia9=atr[SUNTrace[hv3 . v2 . xprog[x, lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v1 . qgq[lora, lorb, cola]]] cqgq[b1, d2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv3 . xprop[-x] . v2]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . qgq[lora, lorb, cola]]] cqgq[b1, d2] cqq[d1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . hv4 . xprop[-x] . v2 . xprog[x, lora, lorb, cola]]] cqgq[a2, a1] cqq[b1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . qgq[lora, lorb, cola] . hv4 . v2 . xprog[x, lora, lorb, cola]]] cqgq[b1, d2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . xprog[x, lora, lorb, cola] . hv4 . xprop[-x] . v2]] cqgq[a2, a1] cqq[d1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . v2 . xprop[x]]] cqgq[c2, d2] cqq[b1, c1] FlavorDelta[a2, a1] FlavorDelta[b1, c1] FlavorDelta[c2, d2] FlavorDelta[d1, b2] - atr[SUNTrace[v2]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v1 . xprog[x, lora, lorb, cola]]] cqgq[a2, b2] cqq[d1, c1] FlavorDelta[a2, b2] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, c1] - atr[SUNTrace[hv4 . v1 . xprop[x] . hv3 . qgq[lora, lorb, cola] . v2 . xprog[x, lora, lorb, cola]]] cqgq[a2, c1] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . v2 . xprog[x, lora, lorb, cola] . hv4 . qgq[lora, lorb, cola]]] cqgq[c2, b2] cqq[b1, c1] FlavorDelta[a2, a1] FlavorDelta[b1, c1] FlavorDelta[c2, b2] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1]] atr[SUNTrace[hv4 . xprop[-x] . v2 . xprog[x, lora, lorb, cola]]] cqgq[a2, a1] cqq[b1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia9=FourierXP[dia9,{x,q}]/.SUNN->CA;

dia9=QEvaluate[I ScaleMu^(4-D) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,cola,lora,lorb},


dia10=-(atr[SUNTrace[hv4 . xprop[-x] . v1 . v2 . xprog[x, lora, lorb, cola] . hv3 . qgq[lora, lorb, cola]]] cqgq[a2, d2] cqq[b1, c1] FlavorDelta[a2, d2] FlavorDelta[b1, c1] FlavorDelta[c2, a1] FlavorDelta[d1, b2]) + atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v2 . xprop[x]]] cqgq[c2, a1] cqq[b1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v2]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . xprop[x]]] cqgq[c2, a1] cqq[d1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . qgq[lora, lorb, cola] . hv4 . v2 . xprop[x]]] cqgq[b1, d2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[b1, c1] cqq[c2, d2] FlavorDelta[a2, a1] FlavorDelta[b1, c1] FlavorDelta[c2, d2] FlavorDelta[d1, b2] + atr[SUNTrace[hv4]] atr[SUNTrace[hv3 . xprop[-x] . v1 . qgq[lora, lorb, cola] . v2 . xprog[x, lora, lorb, cola]]] cqgq[b1, c1] cqq[c2, d2] FlavorDelta[a2, a1] FlavorDelta[b1, c1] FlavorDelta[c2, d2] FlavorDelta[d1, b2] - atr[SUNTrace[v2]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . xprop[x]]] cqgq[c2, d2] cqq[d1, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, d2] FlavorDelta[d1, c1] - atr[SUNTrace[hv4 . v1 . qgq[lora, lorb, cola] . hv3 . xprog[-x, lora, lorb, cola] . v2 . xprop[x]]] cqgq[b1, b2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . v1 . xprog[x, lora, lorb, cola] . hv3 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v2]] atr[SUNTrace[hv3 . xprop[-x] . v1 . xprog[x, lora, lorb, cola]]] cqgq[c2, c1] cqq[d1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia10=FourierXP[dia10,{x,q}]/.SUNN->CA;

dia10=QEvaluate[I ScaleMu^(4-D) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,cola,lora,lorb},


dia11=-(atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . v2 . xprop[x] . hv3 . qgq[lora, lorb, cola]]] cqgq[a2, d2] cqq[b1, c1] FlavorDelta[a2, d2] FlavorDelta[b1, c1] FlavorDelta[c2, a1] FlavorDelta[d1, b2]) + atr[SUNTrace[hv3 . v2 . xprop[x]]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . xprog[x, lora, lorb, cola]]] cqgq[c2, a1] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . v1 . xprog[x, lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . xprop[x] . hv4]] cqgq[d1, c1] cqq[c2, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, b2] FlavorDelta[d1, c1] - atr[SUNTrace[hv4 . xprop[-x] . v1 . hv3 . qgq[lora, lorb, cola] . v2 . xprog[x, lora, lorb, cola]]] cqgq[a2, c1] cqq[b1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . qgq[lora, lorb, cola] . hv3 . xprop[-x] . v2]] cqgq[b1, b2] cqq[d1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v1 . v2 . xprog[x, lora, lorb, cola]]] cqgq[a2, b2] cqq[b1, c1] FlavorDelta[a2, b2] FlavorDelta[b1, c1] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . qgq[lora, lorb, cola] . v2 . xprop[x] . hv4]] cqgq[b1, c1] cqq[c2, b2] FlavorDelta[a2, a1] FlavorDelta[b1, c1] FlavorDelta[c2, b2] FlavorDelta[d1, d2] - atr[SUNTrace[hv3]] atr[SUNTrace[v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v2 . xprop[x]]] cqgq[b1, a1] cqq[a2, b2] FlavorDelta[a2, b2] FlavorDelta[b1, a1] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . v1 . xprog[x, lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia11=FourierXP[dia11,{x,q}]/.SUNN->CA;

dia11=QEvaluate[I ScaleMu^(4-D) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,cola,lora,lorb},


dia12=-(atr[SUNTrace[hv3 . xprop[-x] . v2 . v1 . xprog[x, lora, lorb, cola] . hv4 . qgq[lora, lorb, cola]]] cqgq[c2, b2] cqq[d1, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, b2] FlavorDelta[d1, a1]) + atr[SUNTrace[hv4 . v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v2 . xprog[x, lora, lorb, cola]]] cqgq[b1, d2] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . v1 . xprop[x] . hv4 . xprog[-x, lora, lorb, cola] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v1 . v2 . xprog[x, lora, lorb, cola]]] cqgq[c2, d2] cqq[b1, c1] FlavorDelta[a2, a1] FlavorDelta[b1, c1] FlavorDelta[c2, d2] FlavorDelta[d1, b2] - atr[SUNTrace[hv3]] atr[SUNTrace[v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . xprop[x]]] cqgq[d1, c1] cqq[a2, b2] FlavorDelta[a2, b2] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, c1] - atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . hv3 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[a2, c1] cqq[b1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . qgq[lora, lorb, cola] . v2 . xprog[x, lora, lorb, cola] . hv4]] cqgq[b1, c1] cqq[c2, b2] FlavorDelta[a2, a1] FlavorDelta[b1, c1] FlavorDelta[c2, b2] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . v1 . xprop[x]]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v2 . xprog[x, lora, lorb, cola]]] cqgq[c2, c1] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1]] atr[SUNTrace[hv4 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[b1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . v2 . xprog[x, lora, lorb, cola]]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . xprop[x]]] cqgq[a2, a1] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia12=FourierXP[dia12,{x,q}]/.SUNN->CA;

dia12=QEvaluate[I ScaleMu^(4-D) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13,cola,lora,lorb},


dia13=-(atr[SUNTrace[hv4 . xprop[-x] . v2 . v1 . xprog[x, lora, lorb, cola] . hv3 . qgq[lora, lorb, cola]]] cqgq[a2, d2] cqq[d1, a1] FlavorDelta[a2, d2] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, a1]) + atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[a2, c1] cqq[b1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v2]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . xprop[x]]] cqgq[a2, c1] cqq[d1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[v1]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v2 . xprop[x]]] cqgq[c2, d2] cqq[b1, a1] FlavorDelta[a2, c1] FlavorDelta[b1, a1] FlavorDelta[c2, d2] FlavorDelta[d1, b2] + atr[SUNTrace[v1]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v2 . xprop[x] . hv4 . qgq[lora, lorb, cola]]] cqgq[c2, b2] cqq[b1, a1] FlavorDelta[a2, c1] FlavorDelta[b1, a1] FlavorDelta[c2, b2] FlavorDelta[d1, d2] + atr[SUNTrace[v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v2 . xprop[x] . hv4]] cqgq[b1, a1] cqq[c2, b2] FlavorDelta[a2, c1] FlavorDelta[b1, a1] FlavorDelta[c2, b2] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v2 . xprop[x]]] cqgq[b1, b2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . v2 . xprog[x, lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v1 . qgq[lora, lorb, cola]]] cqgq[b1, b2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . xprop[-x] . v2]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . qgq[lora, lorb, cola]]] cqgq[b1, b2] cqq[d1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v2]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . xprop[x]]] cqgq[c2, c1] cqq[d1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia13=FourierXP[dia13,{x,q}]/.SUNN->CA;

dia13=QEvaluate[I ScaleMu^(4-D) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14,cola,lora,lorb},


dia14=atr[SUNTrace[v1]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v2 . xprop[x] . hv3 . qgq[lora, lorb, cola]]] cqgq[a2, d2] cqq[b1, a1] FlavorDelta[a2, d2] FlavorDelta[b1, a1] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . v1 . qgq[lora, lorb, cola] . hv4 . xprog[-x, lora, lorb, cola] . v2 . xprop[x]]] cqgq[b1, d2] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . xprop[x] . hv4 . v2 . xprog[x, lora, lorb, cola]]] cqgq[a2, a1] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . qgq[lora, lorb, cola] . hv4 . xprop[-x] . v2]] cqgq[b1, d2] cqq[d1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4]] atr[SUNTrace[v1 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v2 . xprog[x, lora, lorb, cola]]] cqgq[b1, a1] cqq[c2, d2] FlavorDelta[a2, c1] FlavorDelta[b1, a1] FlavorDelta[c2, d2] FlavorDelta[d1, b2] - atr[SUNTrace[v2]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . xprop[x]]] cqgq[a2, b2] cqq[d1, c1] FlavorDelta[a2, b2] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, c1] - atr[SUNTrace[hv4 . xprop[-x] . v1 . qgq[lora, lorb, cola] . hv3 . v2 . xprog[x, lora, lorb, cola]]] cqgq[b1, b2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . xprog[x, lora, lorb, cola] . hv3 . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[v1]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v2 . xprop[x]]] cqgq[a2, b2] cqq[b1, a1] FlavorDelta[a2, b2] FlavorDelta[b1, a1] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . v2 . xprop[x]]] atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1 . qgq[lora, lorb, cola]]] cqgq[b1, b2] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia14=FourierXP[dia14,{x,q}]/.SUNN->CA;

dia14=QEvaluate[I ScaleMu^(4-D) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15,cola,lora,lorb},


dia15=-(atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v2 . qgq[lora, lorb, cola] . v1 . xprop[x] . hv3]] cqgq[d1, a1] cqq[a2, d2] FlavorDelta[a2, d2] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, a1]) + atr[SUNTrace[hv3 . v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v1 . xprog[x, lora, lorb, cola]]] cqgq[d1, b2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . xprop[-x] . v1]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v2 . xprog[x, lora, lorb, cola]]] cqgq[a2, c1] cqq[b1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . v1 . xprop[x]]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v2 . xprog[x, lora, lorb, cola]]] cqgq[a2, c1] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv3 . xprop[-x] . v2]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . xprog[x, lora, lorb, cola]]] cqgq[c2, a1] cqq[d1, b2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . v1 . xprog[x, lora, lorb, cola] . hv4 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[c2, c1] cqq[a2, a1] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] + atr[SUNTrace[v2 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv3 . xprop[-x] . v1 . xprog[x, lora, lorb, cola] . hv4]] cqgq[d1, c1] cqq[c2, b2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, b2] FlavorDelta[d1, c1] - atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . qgq[lora, lorb, cola] . hv3 . v2 . xprop[x]]] cqgq[b1, b2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . xprop[-x] . v1 . xprog[x, lora, lorb, cola] . hv3 . qgq[lora, lorb, cola] . v2]] cqgq[a2, c1] cqq[d1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] + atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v1]] atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[c2, c1] cqq[b1, b2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia15=FourierXP[dia15,{x,q}]/.SUNN->CA;

dia15=QEvaluate[I ScaleMu^(4-D) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16,cola,lora,lorb},


dia16=-(atr[SUNTrace[hv3 . xprog[-x, lora, lorb, cola] . v2 . v1 . xprop[x] . hv4 . qgq[lora, lorb, cola]]] cqgq[c2, b2] cqq[d1, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, b2] FlavorDelta[d1, a1]) - atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . qgq[lora, lorb, cola] . v2 . xprop[x] . hv3]] cqgq[b1, c1] cqq[a2, d2] FlavorDelta[a2, d2] FlavorDelta[b1, c1] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1]] atr[SUNTrace[hv3 . xprop[-x] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[b1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] + atr[SUNTrace[hv4 . v1 . xprog[x, lora, lorb, cola]]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v2 . xprop[x]]] cqgq[a2, c1] cqq[c2, a1] FlavorDelta[a2, c1] FlavorDelta[b1, d2] FlavorDelta[c2, a1] FlavorDelta[d1, b2] - atr[SUNTrace[hv3 . xprop[-x] . v1 . hv4 . xprog[-x, lora, lorb, cola] . v2 . qgq[lora, lorb, cola]]] cqgq[d1, b2] cqq[b1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, d2] FlavorDelta[c2, c1] FlavorDelta[d1, b2] - atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v1 . xprop[x] . hv3 . v2 . qgq[lora, lorb, cola]]] cqgq[d1, d2] cqq[a2, c1] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[hv4 . qgq[lora, lorb, cola] . v1 . xprog[x, lora, lorb, cola] . hv3 . xprop[-x] . v2]] cqgq[c2, a1] cqq[d1, d2] FlavorDelta[a2, c1] FlavorDelta[b1, b2] FlavorDelta[c2, a1] FlavorDelta[d1, d2] - atr[SUNTrace[v1]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola]]] atr[SUNTrace[hv4 . xprop[-x] . v2 . xprog[x, lora, lorb, cola]]] cqgq[a2, b2] cqq[b1, a1] FlavorDelta[a2, b2] FlavorDelta[b1, a1] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . v2 . xprop[x]]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . xprog[x, lora, lorb, cola]]] cqgq[a2, a1] cqq[c2, c1] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2] + atr[SUNTrace[hv4 . xprog[-x, lora, lorb, cola] . v2]] atr[SUNTrace[hv3 . qgq[lora, lorb, cola] . v1 . xprop[x]]] cqgq[a2, a1] cqq[d1, d2] FlavorDelta[a2, a1] FlavorDelta[b1, b2] FlavorDelta[c2, c1] FlavorDelta[d1, d2];


dia16=FourierXP[dia16,{x,q}]/.SUNN->CA;

dia16=QEvaluate[I ScaleMu^(4-D) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
