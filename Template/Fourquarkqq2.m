(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



Fourquarkqq2::usage="Fourquarkqq2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},atr_] give the \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)q\!\(\*SuperscriptBox[\(\[RightAngleBracket]\), \(2\)]\) contribution of OPE result about <J1 J2>, with J1=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a1\)]\)\!\(\*SubscriptBox[\(v1\[CapitalPsi]\), \(b1\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c1\)]\)\!\(\*SubscriptBox[\(v2\[CapitalPsi]\), \(d1\)]\) and J2=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(b2\)]\)\!\(\*SubscriptBox[\(v3\[CapitalPsi]\), \(a2\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(d2\)]\)\!\(\*SubscriptBox[\(v4\[CapitalPsi]\), \(c2\)]\) "

Fourquarkqq2::inderr="Dummy indices conflict!"

Begin["`Private`Fourquarkqq2`"]



Options[Fourquarkqq2]={
	Parallelized->True,
	AutoNDR->True,
	HoldFlavor->False
}



Fourquarkqq2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmp,tmp1,tmp2,hv3,hv4,atr,holdf=OptionValue[HoldFlavor],fdir,pall,files,diagrams,dot},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[Fourquarkqq2::inderr];
	Abort[]
];


(*-------------------*)


If[OptionValue[AutoNDR]===True,
	atr=TR5
	,
	atr=TR
];


(* B A^+ B *)
hv3=v3//ComplexConjugate;
hv4=v4//ComplexConjugate//FCI;
(* ComplexConjugate[sigma] will expand the DiracSigma to gamma matrices, recover it to DiracSigma *)
hv4=(hv4/.Dot->dot)/.f1_ dot[aa___,DiracGamma[LorentzIndex[lor1_,dim___],dim___],DiracGamma[LorentzIndex[lor2_,dim___],dim___],bb___]+f2_ dot[aa___,DiracGamma[LorentzIndex[lor2_,dim___],dim___],DiracGamma[LorentzIndex[lor1_,dim___],dim___],bb___]/;(f1+f2===0):>-2I f1 dot[aa,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],bb];
hv4=hv4/.{null->1,dot->Dot};

(*--------------------------------------------*)		

					
If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];


diagrams={type1,type2,type3,type4,type5};


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
			files=("Fourquarkqq2_"<>#)&/@files;
			
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



(*------------------------------------------------------------------*)
(* propagators *)

(*prop[q_]=I GSD[q]FAD[q];*)


(*xprop[x_]=FourierPX[prop[q],{q,x}];*)
xprop[x_]=1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];



(* the follwoing are generated by algorithem *)


(* ::Input::Initialization:: *)
type1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1},


dia1=(Condensate[{a2, b2}]*Condensate[{d1, a1}]*Contract[atr[SUNTrace[hv3]]*atr[SUNTrace[v2 . v1 . xprop[x] . hv4 . xprop[-x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1])/(16*CA^2) + (Condensate[{a2, c1}]*Condensate[{c2, a1}]*Contract[atr[SUNTrace[hv3 . v2 . xprop[x]]]*atr[SUNTrace[hv4 . v1 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2])/(16*CA^2) + (Condensate[{c2, a1}]*Condensate[{d1, b2}]*Contract[atr[SUNTrace[hv3 . xprop[-x] . v2]]*atr[SUNTrace[hv4 . v1 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2])/(16*CA^2) - (Condensate[{b1, d2}]*Condensate[{c2, c1}]*Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . hv4 . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2])/(16*CA^2) - (Condensate[{c2, c1}]*Condensate[{d1, b2}]*Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . v2]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2])/(16*CA^2) + (Condensate[{b1, c1}]*Condensate[{c2, d2}]*Contract[atr[SUNTrace[hv4]]*atr[SUNTrace[hv3 . xprop[-x] . v1 . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2])/(16*CA^2) - (Condensate[{a2, b2}]*Condensate[{d1, c1}]*Contract[atr[SUNTrace[hv3]]*atr[SUNTrace[v2]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1])/(16*CA^2) - (Condensate[{a2, c1}]*Condensate[{c2, a1}]*Contract[atr[SUNTrace[hv4 . v1 . xprop[x] . hv3 . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2])/(16*CA^2);


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(4-D) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2},


dia2=(Condensate[{a2, c1}]*Condensate[{b1, d2}]*Contract[atr[SUNTrace[hv3 . v2 . xprop[x]]]*atr[SUNTrace[hv4 . xprop[-x] . v1]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2])/(16*CA^2) + (Condensate[{b1, d2}]*Condensate[{d1, b2}]*Contract[atr[SUNTrace[hv3 . xprop[-x] . v2]]*atr[SUNTrace[hv4 . xprop[-x] . v1]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2])/(16*CA^2) - (Condensate[{a2, a1}]*Condensate[{d1, b2}]*Contract[atr[SUNTrace[hv3 . v1 . xprop[x] . hv4 . xprop[-x] . v2]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2])/(16*CA^2) - (Condensate[{b1, d2}]*Condensate[{d1, b2}]*Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . hv4 . xprop[-x] . v2]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2])/(16*CA^2) - (Condensate[{b1, a1}]*Condensate[{c2, d2}]*Contract[atr[SUNTrace[hv4]]*atr[SUNTrace[v1]]*atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2])/(16*CA^2) - (Condensate[{b1, b2}]*Condensate[{c2, a1}]*Contract[atr[SUNTrace[hv4 . v1 . hv3 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2])/(16*CA^2) - (Condensate[{a2, c1}]*Condensate[{d1, d2}]*Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . v2]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2])/(16*CA^2) + (Condensate[{a2, a1}]*Condensate[{d1, d2}]*Contract[atr[SUNTrace[hv3 . v1 . xprop[x]]]*atr[SUNTrace[hv4 . xprop[-x] . v2]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2])/(16*CA^2);


dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(4-D) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3},


dia3=-(Condensate[{a2, d2}]*Condensate[{d1, a1}]*Contract[atr[SUNTrace[hv4 . xprop[-x] . v2 . v1 . xprop[x] . hv3]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1])/(16*CA^2) + (Condensate[{a2, d2}]*Condensate[{b1, a1}]*Contract[atr[SUNTrace[v1]]*atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x] . hv3]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2])/(16*CA^2) + (Condensate[{a2, d2}]*Condensate[{d1, c1}]*Contract[atr[SUNTrace[v2]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1])/(16*CA^2) - (Condensate[{c2, a1}]*Condensate[{d1, d2}]*Contract[atr[SUNTrace[hv4 . v1 . xprop[x] . hv3 . xprop[-x] . v2]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2])/(16*CA^2) + (Condensate[{a2, b2}]*Condensate[{b1, c1}]*Contract[atr[SUNTrace[hv3]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . v2 . xprop[x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2])/(16*CA^2) + (Condensate[{b1, a1}]*Condensate[{c2, b2}]*Contract[atr[SUNTrace[v1]]*atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x] . hv4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2])/(16*CA^2) + (Condensate[{b1, b2}]*Condensate[{c2, c1}]*Contract[atr[SUNTrace[hv3 . xprop[-x] . v1]]*atr[SUNTrace[hv4 . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2])/(16*CA^2) + (Condensate[{a2, a1}]*Condensate[{b1, b2}]*Contract[atr[SUNTrace[hv3 . v1]]*atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2])/(16*CA^2);


dia3=FourierXP[dia3,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=QEvaluate[I ScaleMu^(4-D) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4},


dia4=-(Condensate[{c2, b2}]*Condensate[{d1, a1}]*Contract[atr[SUNTrace[hv3 . xprop[-x] . v2 . v1 . xprop[x] . hv4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, a1])/(16*CA^2) + (Condensate[{c2, d2}]*Condensate[{d1, a1}]*Contract[atr[SUNTrace[hv4]]*atr[SUNTrace[v2 . v1 . xprop[x] . hv3 . xprop[-x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, a1])/(16*CA^2) - (Condensate[{a2, d2}]*Condensate[{b1, c1}]*Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . v2 . xprop[x] . hv3]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2])/(16*CA^2) + (Condensate[{b1, d2}]*Condensate[{c2, a1}]*Contract[atr[SUNTrace[hv4 . v1]]*atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2])/(16*CA^2) + (Condensate[{a2, c1}]*Condensate[{d1, b2}]*Contract[atr[SUNTrace[hv3 . v2]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2])/(16*CA^2) + (Condensate[{c2, b2}]*Condensate[{d1, c1}]*Contract[atr[SUNTrace[v2]]*atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, c1])/(16*CA^2) - (Condensate[{c2, d2}]*Condensate[{d1, c1}]*Contract[atr[SUNTrace[hv4]]*atr[SUNTrace[v2]]*atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, c1])/(16*CA^2) - (Condensate[{a2, c1}]*Condensate[{b1, b2}]*Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . hv3 . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2])/(16*CA^2);


dia4=FourierXP[dia4,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=QEvaluate[I ScaleMu^(4-D) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5},


dia5=-(Condensate[{a2, a1}]*Condensate[{b1, d2}]*Contract[atr[SUNTrace[hv3 . v1 . hv4 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2])/(16*CA^2) - (Condensate[{a2, a1}]*Condensate[{c2, c1}]*Contract[atr[SUNTrace[hv3 . v1 . xprop[x] . hv4 . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2])/(16*CA^2) - (Condensate[{b1, b2}]*Condensate[{d1, d2}]*Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . hv3 . xprop[-x] . v2]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2])/(16*CA^2) - (Condensate[{b1, c1}]*Condensate[{c2, b2}]*Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . v2 . xprop[x] . hv4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2])/(16*CA^2) - (Condensate[{a2, b2}]*Condensate[{b1, a1}]*Contract[atr[SUNTrace[hv3]]*atr[SUNTrace[v1]]*atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2])/(16*CA^2) + (Condensate[{a2, a1}]*Condensate[{c2, c1}]*Contract[atr[SUNTrace[hv3 . v1 . xprop[x]]]*atr[SUNTrace[hv4 . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2])/(16*CA^2) + (Condensate[{b1, b2}]*Condensate[{d1, d2}]*Contract[atr[SUNTrace[hv3 . xprop[-x] . v1]]*atr[SUNTrace[hv4 . xprop[-x] . v2]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2])/(16*CA^2) + (Condensate[{c2, c1}]*Condensate[{d1, d2}]*Contract[atr[SUNTrace[hv4 . v2]]*atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2])/(16*CA^2);


dia5=FourierXP[dia5,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=QEvaluate[I ScaleMu^(4-D) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_1 *)
type1[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},atr_]:=Block[{tmp0,x,q,lora,lorb,cola,colb,k1,k2,k3,l,d1a,d1b,d1c,d1d,dia1},



d1a=If[b1===b2&&a2===a1&&d1===d2&&c1===c2,atr[SUNTrace[g1.xprop[x].g3.xprop[-x]]]atr[SUNTrace[g2.g4]]qq2 cqq[d1,d2]cqq[c1,c2],0]+If[b1===b2&&a2===a1&&c1===d1&&c2===d2,-atr[SUNTrace[g1.xprop[x].g3.xprop[-x]]]atr[SUNTrace[g2]]atr[SUNTrace[g4]]qq2 cqq[c1,d1]cqq[c2,d2],0];

d1b=If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[g1.xprop[x].g3.g2.g4.xprop[-x]]]qq2 cqq[a2,c1]cqq[d1,d2],0]+If[b1===b2&&a2===d2&&c2===a1&&c1===d1,atr[SUNTrace[g1.xprop[x].g3.g4.xprop[-x]]]atr[SUNTrace[g2]]qq2 cqq[a2,d2]cqq[c1,d1],0];

d1c=If[b1===d2&&c2===c1&&d1===b2&&a2===a1,-atr[SUNTrace[g1.xprop[x].g4.g2.g3.xprop[-x]]]qq2 cqq[c2,c1]cqq[d1,b2],0]+If[b1===d2&&c2===b2&&a2===a1&&c1===d1,atr[SUNTrace[g1.xprop[x].g4.g3.xprop[-x]]]atr[SUNTrace[g2]]qq2 cqq[c2,b2]cqq[c1,d1],0];

d1d=If[b1===d2&&c2===a1&&d1===b2&&a2===c1,atr[SUNTrace[g1.xprop[x].g4.xprop[-x]]]atr[SUNTrace[g2.g3]]qq2 cqq[d1,b2]cqq[c1,a2],0]+If[b1===d2&&c2===a1&&c1===d1&&a2===b2,-atr[SUNTrace[g1.xprop[x].g4.xprop[-x]]]atr[SUNTrace[g2]]atr[SUNTrace[g3]]qq2 cqq[c1,d1]cqq[a2,b2],0];


(*-----------------------------*)
dia1=FourierXP[d1a+d1b+d1c+d1d,{x,q}];

QEvaluate[I dia1 ScaleMu^(4-D),q,Parallelized->False]/.q->qq

]
*)


(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_2 *)
type2[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},atr_]:=Block[{tmp0,x,q,lora,lorb,cola,colb,k1,k2,k3,l,d2a,d2b,d2c,d2d,dia2},


d2a=If[b1===b2&&a2===a1&&d1===d2&&c2===c1,atr[SUNTrace[g1.xprop[x].g3]]atr[SUNTrace[g2.g4.xprop[-x]]]qq2 cqq[a1,a2]cqq[d1,d2],0]+If[b1===b2&&a2===d2&&c2===c1&&d1===a1,-atr[SUNTrace[g1.xprop[x].g3.g4.xprop[-x].g2]]qq2 cqq[a2,d2]cqq[a1,d1],0];

d2b=If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[g1.xprop[x].g3.xprop[-x].g2.g4]]qq2 cqq[d1,d2]cqq[c2,a1],0]+If[b1===b2&&a2===c1&&d1===a1&&c2===d2,atr[SUNTrace[g1.xprop[x].g3.xprop[-x].g2]]atr[SUNTrace[g4]]qq2 cqq[a1,d1]cqq[c2,d2],0];

d2c=If[b1===d2&&c2===c1&&d1===b2&&a2===a1,-atr[SUNTrace[g1.xprop[x].g4.xprop[-x].g2.g3]]qq2 cqq[d1,b2]cqq[a2,a1],0]+If[b1===d2&&c2===c1&&d1===a1&&b2===a2,atr[SUNTrace[g1.xprop[x].g4.xprop[-x].g2]]atr[SUNTrace[g3]]qq2 cqq[a1,d1]cqq[a2,b2],0];

d2d=If[b1===d2&&c2===a1&&a2===c1&&d1===b2,atr[SUNTrace[g1.xprop[x].g4]]atr[SUNTrace[g3.xprop[-x].g2]]qq2 cqq[b1,c2]cqq[d1,b2],0]+If[b1===d2&&c2===b2&&a2===c1&&d1===a1,-atr[SUNTrace[g1.xprop[x].g4.g3.xprop[-x].g2]]qq2 cqq[c2,b2]cqq[d1,a1],0];



(*-----------------------------*)
dia2=FourierXP[d2a+d2b+d2c+d2d,{x,q}];

QEvaluate[I dia2 ScaleMu^(4-D),q,Parallelized->False]/.q->qq
]
*)


(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_3 *)
type3[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},atr_]:=Block[{tmp0,x,q,lora,lorb,cola,colb,k1,k2,k3,l,d3a,d3b,d3c,d3d,dia3},


d3a=If[b1===b2&&a2===a1&&d1===d2&&c2===c1,atr[SUNTrace[g1.g3.xprop[-x]]]atr[SUNTrace[g2.xprop[x].g4]]qq2 cqq[b1,b2]cqq[c2,c1],0]+If[a2===a1&&b1===c1&&d1===d2&&c2===b2,-atr[SUNTrace[g3.xprop[-x].g1.g2.xprop[x].g4]]qq2 cqq[b1,c1]cqq[c2,b2],0];

d3b=If[b1===b2&&a2===c1&&d1===d2&&c2===b1,-atr[SUNTrace[g1.g3.g2.xprop[x].g4.xprop[-x]]]qq2 cqq[b1,b2]cqq[a2,c1],0]+If[b1===c1&&d1===d2&&c2===a1&&a2===b2,atr[SUNTrace[g1.g2.xprop[x].g4.xprop[-x]]]atr[SUNTrace[g3]]qq2 cqq[a2,b2]cqq[b1,c1],0];

d3c=If[b1===d2&&c2===c1&&d1===b2&&a2===a1,-atr[SUNTrace[g1.g4.g2.xprop[x].g3.xprop[-x]]]qq2 cqq[b1,d2]cqq[c2,c1],0]+If[b1===c1&&d1===b2&&a2===a1&&c2===d2,atr[SUNTrace[g1.g2.xprop[x].g3.xprop[-x]]]atr[SUNTrace[g4]]qq2 cqq[b1,c1]cqq[c2,d2],0];

d3d=If[b1===d2&&c2===a1&&d1===b2&&a2===c1,atr[SUNTrace[g1.g4.xprop[-x]]]atr[SUNTrace[g2.xprop[x].g3]]qq2 cqq[b1,d2]cqq[a2,c1],0]+If[b1===c1&&d1===b2&&a2===d2&&c2===a1,-atr[SUNTrace[g1.g2.xprop[x].g3.g4.xprop[-x]]]qq2 cqq[b1,c1]cqq[a2,d2],0];


(*-----------------------------*)
dia3=FourierXP[d3a+d3b+d3c+d3d,{x,q}];

QEvaluate[I dia3 ScaleMu^(4-D),q,Parallelized->False]/.q->qq

]*)


(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_4 *)
type4[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},atr_]:=Block[{tmp0,x,q,lora,lorb,cola,colb,k1,k2,k3,l,d4a,d4b,d4c,d4d,dia4},

d4a=If[b1===b2&&a2===a1&&d1===d2&&c2===c1,atr[SUNTrace[g1.g3]]atr[SUNTrace[g2.xprop[x].g4.xprop[-x]]]qq2 cqq[b1,b2]cqq[a1,a2],0]+If[a1===b1&&a2===b2&&c1===c2&&d2===d1,-atr[SUNTrace[g1]]atr[SUNTrace[g3]]atr[SUNTrace[g2.xprop[x].g4.xprop[-x]]]qq2 cqq[a1,b1]cqq[a2,b2],0];

d4b=If[a2===a1&&b1===d2&&c2===c1&&d1===b2,-atr[SUNTrace[g3.g1.g4.xprop[-x].g2.xprop[x]]]qq2 cqq[a2,a1]cqq[b1,d2],0]+If[a1===b1&&a2===d2&&c2===c1&&d1===b2,atr[SUNTrace[g1]]atr[SUNTrace[g2.xprop[x].g3.g4.xprop[-x]]]qq2 cqq[a1,b1]cqq[a2,d2],0];

d4c=If[b1===b2&&a2===c1&&d1===d2&&c2===a1,atr[SUNTrace[g1.g3.xprop[-x].g2.xprop[x].g4]]qq2 cqq[b1,b2]cqq[c2,a1],0]+If[a1===b1&&a2===c1&&d1===d2&&c2===b2,atr[SUNTrace[g1]]atr[SUNTrace[g3.xprop[-x].g2.xprop[x].g4]]qq2 cqq[a1,b1]cqq[c2,b2],0];

d4d=If[d1===d2&&a2===c1&&b1===d2&&c2===a1,atr[SUNTrace[g2.xprop[x].g3.xprop[-x]]]atr[SUNTrace[g1.g4]]qq2 cqq[b1,d2]cqq[c2,a1],0]+If[a1===b1&&c2===d2&&d1===b2&&a2===c1,-atr[SUNTrace[g1]]atr[SUNTrace[g4]]atr[SUNTrace[g2.xprop[x].g3.xprop[-x]]]qq2 cqq[a1,b1]cqq[c2,d2],0];



(*-----------------------------*)
dia4=FourierXP[d4a+d4b+d4c+d4d,{x,q}];

QEvaluate[I dia4 ScaleMu^(4-D),q,Parallelized->False]/.q->qq
]

*)


(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_5 *)
type5[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},atr_]:=Block[{tmp0,x,q,lora,lorb,cola,colb,k1,k2,k3,l,d5a,d5b,d5c,d5d,dia5},

d5a=If[b1===b2&&a2===a1&&d1===d2&&c2===c1,atr[SUNTrace[g1.xprop[x].g3]]atr[SUNTrace[g2.xprop[x].g4]]qq2 cqq[a2,a1]cqq[c1,c2],0]+If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[g1.xprop[x].g3.g2.xprop[x].g4]]qq2 cqq[a2,c1]cqq[c2,a1],0];

d5b=If[b1===d2&&c2===a1&&d1===b2&&a2===c1,atr[SUNTrace[g1.xprop[x].g4]]atr[SUNTrace[g2.xprop[x].g3]]qq2 cqq[c2,a1]cqq[a2,c1],0]+If[b1===d2&&c2===c1&&d1===b2&&a2===a1,-atr[SUNTrace[g1.xprop[x].g4.g2.xprop[x].g3]]qq2 cqq[c2,c1]cqq[a2,a1],0];

d5c=If[a1===a2&&b1===b2&&d1===d2&&c2===c1,atr[SUNTrace[g3.xprop[-x].g1]]atr[SUNTrace[g4.xprop[-x].g2]]qq2 cqq[a1,a2]cqq[c1,c2],0]+If[b1===d2&&c2===c1&&d1===b2&&a2===a1,-atr[SUNTrace[g3.xprop[-x].g1.g4.xprop[-x].g2]]qq2 cqq[b1,d2]cqq[d1,b2],0];

d5d=If[b1===d2&&c2===a1&&a2===c1&&d1===b2,atr[SUNTrace[g4.xprop[-x].g1]]atr[SUNTrace[g3.xprop[-x].g2]]qq2 cqq[b1,d2]cqq[d1,b2],0]+If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[g1.g3.xprop[-x].g2.g4.xprop[-x]]]qq2 cqq[b1,b2]cqq[d1,d2],0];


(*-----------------------------*)
dia5=FourierXP[d5a+d5b+d5c+d5d,{x,q}];

QEvaluate[I dia5 ScaleMu^(4-D),q,Parallelized->False]/.q->qq

]


*)


(*-------------------------------------------------------------------------------------------------------*)


(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_1 *)
type1[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},atr_]:=Block[{tmp0,x,q,lora,lorb,cola,colb,k1,k2,k3,l,d1a,d1b,d1c,d1d,dia1},



(* first diagram *)
d1a=If[{a1,b1}==={a2,b2},-atr[SUNTrace[g1.prop[k1+q].g3.prop[k1]]],0]*(
If[{c1,d1}==={c2,d2},-atr[SUNTrace[g2.g4 qq2]],0]cqq[c1,c2]cqq[d1,d2]
+If[c1===d1&&c2===d2,atr[SUNTrace[g2]]qq2 atr[SUNTrace[g4]]cqq[d1,c1]cqq[d2,c2],0]);

(* second diagram *)
d1b=If[{a1,b1}==={c2,b2}&&{c1,d1}==={a2,d2},-atr[SUNTrace[g1.prop[k1+q].g3.g2.g4.prop[k1] qq2]],0]cqq[d1,d2]cqq[c1,a2]
+If[{a1,b1}==={c2,b2}&&c1===d1&&a2===d2,atr[SUNTrace[g1.prop[k1+q].g3.g4.prop[k1] ]]qq2 atr[SUNTrace[g2]],0]cqq[d1,c1]cqq[a2,d2];

(* third diagram *)
d1c=If[a1===a2&&b1===d2&&d1===b2&&c1===c2,-atr[SUNTrace[g3.prop[k1].g1.prop[k1+q].g4.g2 qq2]],0]cqq[d1,b2]cqq[c1,c2]
+If[a1===a2&&b1===d2&&d1===c1&&b2===c2,atr[SUNTrace[g3.prop[k1].g1.prop[k1+q].g4 ]]qq2 atr[SUNTrace[g2]],0]cqq[d1,c1]cqq[b2,c2];

(* fourth diagram *)
d1d=If[{a1,b1}==={c2,d2},-atr[SUNTrace[g1.prop[k1+q].g4.prop[k1]]],0]*(
If[{c1,d1}==={a2,b2},-atr[SUNTrace[g2.g3 qq2]],0]cqq[d1,b2]cqq[c1,a2]
+If[c1===d1&&a2===b2,atr[SUNTrace[g2]]atr[SUNTrace[g3]]qq2,0]cqq[d1,c1]cqq[b2,a2]);






(*(*-----------------------------*)
(* first diagram *)
d1a=If[{a1,b1}==={a2,b2},-atr[SUNTrace[g1 . prop[k1+q] . g3 . prop[k1]]],0](If[{c1,d1}==={c2,d2},-atr[SUNTrace[g2 . g4 qq2]],0]cqq[c1,c2]cqq[d1,d2]+If[c1===d1&&c2===d2,atr[SUNTrace[g2]]qq2 atr[SUNTrace[g4]]cqq[d1,c1]cqq[d2,c2],0]);


(*-----------------------------*)
(* second diagram *)
d1b=If[{a1,b1}==={c2,b2}&&{c1,d1}==={a2,d2},-atr[SUNTrace[g1 . prop[k1+q] . g3 . g2 . g4 . prop[k1] qq2]],0]cqq[d1,d2]cqq[c1,a2]+If[{a1,b1}==={c2,b2}&&c1===d1&&a2===d2,atr[SUNTrace[g1 . prop[k1+q] . g3 . g4 . prop[k1] ]]qq2 atr[SUNTrace[g2]],0]cqq[d1,c1]cqq[a2,d2];


(*-----------------------------*)
(* third diagram *)
d1c=If[a1===a2&&b1===d2&&d1===b2&&c1===c2,-atr[SUNTrace[g3 . prop[k1] . g1 . prop[k1+q] . g4 . g2 qq2]],0]cqq[d1,b2]cqq[c1,c2]+If[a1===a2&&b1===d2&&d1===c1&&b2===c2,atr[SUNTrace[g3 . prop[k1] . g1 . prop[k1+q] . g4 ]]qq2 atr[SUNTrace[g2]],0]cqq[d1,c1]cqq[b2,c2];


(*-----------------------------*)
(* fourth diagram *)
d1d=If[{a1,b1}==={c2,d2},-atr[SUNTrace[g1 . prop[k1+q] . g4 . prop[k1]]],0](If[{c1,d1}==={a2,b2},-atr[SUNTrace[g2 . g3 qq2]],0]cqq[d1,b2]cqq[c1,a2]+If[c1===d1&&a2===b2,atr[SUNTrace[g2]]atr[SUNTrace[g3]]qq2,0]cqq[d1,c1]cqq[b2,a2]);

*)
(*-----------------------------*)

dia1=IntegrateP[d1a+d1b+d1c+d1d,k1]//SUNSimplify;
QEvaluate[I dia1 ScaleMu^(4-D),q,Parallelized->False]/.q->qq

]
(*-------------------------------------------------*)
(*-------------------------------------------------*)*)


(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_2 *)
type2[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},atr_]:=Block[{tmp0,x,q,lora,lorb,cola,colb,k1,k2,k3,l,d2a,d2b,d2c,d2d,dia2},

(* first diagram *)
d2a=If[{a1,b1}==={a2,b2}&&{c1,d1}==={c2,d2},atr[SUNTrace[prop[k1+q].g3.g1]]qq2 atr[SUNTrace[g2.g4.prop[k1]]],0]cqq[a1,a2]cqq[d1,d2]
+If[b1===b2&&a2===d2&&c1===c2&&d1===a1,-atr[SUNTrace[g1.prop[k1+q].g3.g4.prop[k1].g2 qq2]],0]cqq[a1,d1]cqq[a2,d2];

(* second diagram *)
d2b=If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[prop[k1+q].g3.prop[k1].g2.g4.g1 qq2]],0]cqq[d1,d2]cqq[a1,c2]
+If[b1===b2&&a2===c1&&d1===a1&&d2===c2,atr[SUNTrace[g1.prop[k1+q].g3.prop[k1].g2 ]]qq2 atr[SUNTrace[g4]],0]cqq[a1,d1]cqq[c2,d2];

(* third diagram *)
d2c=If[b1===d2&&c2===c1&&d1===b2&&a1===a2,-atr[SUNTrace[g1.prop[k1+q].g4.prop[k1].g2.g3 qq2]],0]cqq[a1,a2]cqq[d1,b2]
+If[b1===d2&&c2===c1&&d1===a1&&b2===a2,atr[SUNTrace[g1.prop[k1+q].g4.prop[k1].g2 ]]qq2 atr[SUNTrace[g3]],0]cqq[a1,d1]cqq[a2,b2];

(* fourth diagram *)
d2d=If[b1===d2&&c2===a1&&a2===c1&&d1===b2,atr[SUNTrace[g1.prop[k1+q].g4]]qq2 atr[SUNTrace[g3.prop[k1].g2]],0]cqq[a1,c2]cqq[d1,b2]
+If[b1===d2&&c2===b2&&a2===c1&&d1===a1,-atr[SUNTrace[g1.prop[k1+q].g4.g3.prop[k1].g2 qq2]],0]cqq[a1,d1]cqq[b2,c2];



(*(*-----------------------------*)
(* first diagram *)
d2a=If[{a1,b1}==={a2,b2}&&{c1,d1}==={c2,d2},atr[SUNTrace[prop[k1+q] . g3 . g1]]qq2 atr[SUNTrace[g2 . g4 . prop[k1]]],0]cqq[a1,a2]cqq[d1,d2]+If[b1===b2&&a2===d2&&c1===c2&&d1===a1,-atr[SUNTrace[g1 . prop[k1+q] . g3 . g4 . prop[k1] . g2 qq2]],0]cqq[a1,d1]cqq[a2,d2];

(*-----------------------------*)
(* second diagram *)
d2b=If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[prop[k1+q] . g3 . prop[k1] . g2 . g4 . g1 qq2]],0]cqq[d1,d2]cqq[a1,c2]+If[b1===b2&&a2===c1&&d1===a1&&d2===c2,atr[SUNTrace[g1 . prop[k1+q] . g3 . prop[k1] . g2 ]]qq2 atr[SUNTrace[g4]],0]cqq[a1,d1]cqq[c2,d2];

(*-----------------------------*)
(* third diagram *)
d2c=If[b1===d2&&c2===c1&&d1===b2&&a1===a2,-atr[SUNTrace[g1 . prop[k1+q] . g4 . prop[k1] . g2 . g3 qq2]],0]cqq[a1,a2]cqq[d1,b2]+If[b1===d2&&c2===c1&&d1===a1&&b2===a2,atr[SUNTrace[g1 . prop[k1+q] . g4 . prop[k1] . g2 ]]qq2 atr[SUNTrace[g3]],0]cqq[a1,d1]cqq[a2,b2];

(*-----------------------------*)
(* fourth diagram *)
d2d=If[b1===d2&&c2===a1&&a2===c1&&d1===b2,atr[SUNTrace[g1 . prop[k1+q] . g4]]qq2 atr[SUNTrace[g3 . prop[k1] . g2]],0]cqq[a1,c2]cqq[d1,b2]+If[b1===d2&&c2===b2&&a2===c1&&d1===a1,-atr[SUNTrace[g1 . prop[k1+q] . g4 . g3 . prop[k1] . g2 qq2]],0]cqq[a1,d1]cqq[b2,c2];

*)
(*-----------------------------*)
dia2=IntegrateP[d2a+d2b+d2c+d2d,k1]//SUNSimplify;
QEvaluate[DD2 I dia2 ScaleMu^(4-D),q,Parallelized->False]/.q->qq
]

(*-------------------------------------------------*)
(*-------------------------------------------------*)
*)


(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_3 *)
type3[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},atr_]:=Block[{tmp0,x,q,lora,lorb,cola,colb,k1,k2,k3,l,d3a,d3b,d3c,d3d,dia3},

(* first diagram *)
d3a=(If[a2===a1&&b1===b2&&d1===d2&&c2===c1,atr[SUNTrace[g1.g3.prop[k1]]]qq2 atr[SUNTrace[g2.prop[k1+q].g4]],0]cqq[b1,b2]cqq[c1,c2]
+If[a1===a2&&b1===c1&&d1===d2&&c2===b2,-atr[SUNTrace[g3.prop[k1].g1.g2.prop[k1+q].g4 qq2]],0]cqq[b1,c1]cqq[c2,b2]);

(* second diagram *)
d3b=(If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[g3.g2.prop[k1+q].g4.prop[k1].g1 qq2]],0]cqq[b1,b2]cqq[c1,a2]
+If[b2===a2&&d1===d2&&c2===a1&&b1===c1,atr[SUNTrace[g3]]qq2 atr[SUNTrace[g1.g2.prop[k1+q].g4.prop[k1]]],0]cqq[b1,c1]cqq[a2,b2]);

(* third diagram *)
d3c=(If[a2===a1&&b1===d2&&c2===c1&&d1===b2,-atr[SUNTrace[g1.g4.g2.prop[k1+q].g3.prop[k1]qq2]],0]cqq[c1,c2]cqq[b1,d2]
+If[b1===c1&&d1===b2&&a2===a1&&d2===c2,atr[SUNTrace[g1.g2.prop[k1+q].g3.prop[k1] ]]qq2 atr[SUNTrace[g4]] ,0]cqq[b1,c1]cqq[c2,d2]);

(* fourth diagram *)
d3d=(If[{a1,b1}==={c2,d2}&&{c1,d1}==={a2,b2},atr[SUNTrace[g1.g4.prop[k1] ]]qq2 atr[SUNTrace[g2.prop[k1+q].g3]],0]cqq[b1,d2]cqq[c1,a2]
+If[b1===c1&&d1===b2&&a2===d2&&c2===a1,-atr[SUNTrace[g2.prop[k1+q].g3.g4.prop[k1].g1 qq2]],0]cqq[b1,c1]cqq[a2,d2]);



(*(*-----------------------------*)
(* first diagram *)
d3a=(If[a2===a1&&b1===b2&&d1===d2&&c2===c1,atr[SUNTrace[g1 . g3 . prop[k1]]]qq2 atr[SUNTrace[g2 . prop[k1+q] . g4]],0]cqq[b1,b2]cqq[c1,c2]+If[a1===a2&&b1===c1&&d1===d2&&c2===b2,-atr[SUNTrace[g3 . prop[k1] . g1 . g2 . prop[k1+q] . g4 qq2]],0]cqq[b1,c1]cqq[c2,b2]);


(*-----------------------------*)
(* second diagram *)
d3b=(If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[g3 . g2 . prop[k1+q] . g4 . prop[k1] . g1 qq2]],0]cqq[b1,b2]cqq[c1,a2]+If[b2===a2&&d1===d2&&c2===a1&&b1===c1,atr[SUNTrace[g3]]qq2 atr[SUNTrace[g1 . g2 . prop[k1+q] . g4 . prop[k1]]],0]cqq[b1,c1]cqq[a2,b2]);


(*-----------------------------*)
(* third diagram *)
d3c=(If[a2===a1&&b1===d2&&c2===c1&&d1===b2,-atr[SUNTrace[g1 . g4 . g2 . prop[k1+q] . g3 . prop[k1]qq2]],0]cqq[c1,c2]cqq[b1,d2]+If[b1===c1&&d1===b2&&a2===a1&&d2===c2,atr[SUNTrace[g1 . g2 . prop[k1+q] . g3 . prop[k1]]]qq2 atr[SUNTrace[g4]] ,0]cqq[b1,c1]cqq[c2,d2]);


(*-----------------------------*)
(* fourth diagram *)
d3d=(If[{a1,b1}==={c2,d2}&&{c1,d1}==={a2,b2},atr[SUNTrace[g1 . g4 . prop[k1]]]qq2 atr[SUNTrace[g2 . prop[k1+q] . g3]],0]cqq[b1,d2]cqq[c1,a2]+If[b1===c1&&d1===b2&&a2===d2&&c2===a1,-atr[SUNTrace[g2 . prop[k1+q] . g3 . g4 . prop[k1] . g1 qq2]],0]cqq[b1,c1]cqq[a2,d2]);
*)
(*-----------------------------*)
dia3=IntegrateP[d3a+d3b+d3c+d3d,k1]//SUNSimplify;
QEvaluate[DD3 I dia3 ScaleMu^(4-D),q,Parallelized->False]/.q->qq

]

(*-------------------------------------------------*)
(*-------------------------------------------------*)*)


(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_4 *)
type4[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},atr_]:=Block[{tmp0,x,q,lora,lorb,cola,colb,k1,k2,k3,l,d4a,d4b,d4c,d4d,dia4},


(* first diagram *)
d4a=If[{c1,d1}==={c2,d2},-atr[SUNTrace[g2.prop[k1+q].g4.prop[k1]]],0]*(
If[{a1,b1}==={a2,b2},-atr[SUNTrace[g1.g3 qq2]],0]cqq[b1,b2]cqq[a1,a2]
+If[a1===b1&&a2===b2,atr[SUNTrace[g1]]qq2 atr[SUNTrace[g3]],0]cqq[a1,b1]cqq[a2,b2]);

(* second diagram *)
d4b=If[d1===b2&&a2===a1&&b1===d2&&c2===c1,-atr[SUNTrace[g2.prop[k1+q].g3.g1.g4.prop[k1] qq2]],0]cqq[a1,a2]cqq[b1,d2]
+If[d1===b2&&a2===d2&&c2===c1&&a1===b1,atr[SUNTrace[g3.g4.prop[k1].g2.prop[k1+q]]]qq2 atr[SUNTrace[g1]],0]cqq[a1,b1]cqq[a2,d2];

(* third diagram *)
d4c=If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[g1.g3.prop[k1].g2.prop[k1+q].g4 qq2]],0]cqq[b1,b2]cqq[a1,c2]
+If[a2===c1&&d1===d2&&c2===b2&&a1===b1,atr[SUNTrace[g3.prop[k1].g2.prop[k1+q].g4]]qq2 atr[SUNTrace[g1]],0]cqq[a1,b1]cqq[b2,c2];

(* fourth diagram *)
d4d=If[{c1,d1}==={a2,b2},-atr[SUNTrace[g2.prop[k1+q].g3.prop[k1]]],0] *(
If[b1===d2&&a1===c2,-atr[SUNTrace[g1.g4 qq2]],0]cqq[b1,d2]cqq[a1,c2]
+If[a1===b1&&c2===d2,atr[SUNTrace[g1]]qq2 atr[SUNTrace[g4]]]cqq[a1,b1]cqq[c2,d2]);


(*(*-----------------------------*)
(* first diagram *)
d4a=If[{c1,d1}==={c2,d2},-atr[SUNTrace[g2 . prop[k1+q] . g4 . prop[k1]]],0](If[{a1,b1}==={a2,b2},-atr[SUNTrace[g1 . g3 qq2]],0]cqq[b1,b2]cqq[a1,a2]+If[a1===b1&&a2===b2,atr[SUNTrace[g1]]qq2 atr[SUNTrace[g3]],0]cqq[a1,b1]cqq[a2,b2]);

(*-----------------------------*)
(* second diagram *)
d4b=If[d1===b2&&a2===a1&&b1===d2&&c2===c1,-atr[SUNTrace[g2 . prop[k1+q] . g3 . g1 . g4 . prop[k1] qq2]],0]cqq[a1,a2]cqq[b1,d2]+If[d1===b2&&a2===d2&&c2===c1&&a1===b1,atr[SUNTrace[g3 . g4 . prop[k1] . g2 . prop[k1+q]]]qq2 atr[SUNTrace[g1]],0]cqq[a1,b1]cqq[a2,d2];

(*-----------------------------*)
(* third diagram *)
d4c=If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[g1 . g3 . prop[k1] . g2 . prop[k1+q] . g4 qq2]],0]cqq[b1,b2]cqq[a1,c2]+If[a2===c1&&d1===d2&&c2===b2&&a1===b1,atr[SUNTrace[g3 . prop[k1] . g2 . prop[k1+q] . g4]]qq2 atr[SUNTrace[g1]],0]cqq[a1,b1]cqq[b2,c2];

(*-----------------------------*)
(* fourth diagram *)
d4d=If[{c1,d1}==={a2,b2},-atr[SUNTrace[g2 . prop[k1+q] . g3 . prop[k1]]],0] (If[b1===d2&&a1===c2,-atr[SUNTrace[g1 . g4 qq2]],0]cqq[b1,d2]cqq[a1,c2]+If[a1===b1&&c2===d2,atr[SUNTrace[g1]]qq2 atr[SUNTrace[g4]]]cqq[a1,b1]cqq[c2,d2]);
*)
(*-----------------------------*)
dia4=IntegrateP[d4a+d4b+d4c+d4d,k1]//SUNSimplify;
QEvaluate[DD4 I dia4 ScaleMu^(4-D),q,Parallelized->False]/.q->qq
]
(*-------------------------------------------------*)
(*-------------------------------------------------*)*)


(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_5 *)
type5[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},atr_]:=Block[{tmp0,x,q,lora,lorb,cola,colb,k1,k2,k3,l,d5a,d5b,d5c,d5d,dia5},

(* first diagram *)
d5a=If[{a1,b1}==={a2,b2}&&{c1,d1}==={c2,d2},atr[SUNTrace[g1.prop[k1+q].g3]]qq2 atr[SUNTrace[g2.prop[-k1].g4]],0]cqq[a1,a2]cqq[c1,c2]
+If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[g1.prop[k1+q].g3.g2.prop[-k1].g4 qq2]],0]cqq[a1,c2]cqq[c1,a2];

(* second diagram *)
d5b=If[{a1,b1}==={c2,d2}&&{c1,d1}==={a2,b2},atr[SUNTrace[g1.prop[-k1].g4]]qq2 atr[SUNTrace[g2.prop[k1+q].g3]],0]cqq[a1,c2]cqq[c1,a2]
+If[b1===d2&&c2===c1&&d1===b2&&a2===a1,-atr[SUNTrace[g1.prop[-k1].g4.g2.prop[k1+q].g3 qq2]],0]cqq[a1,a2]cqq[c1,c2];

(* third diagram *)
d5c=If[{a1,b1}==={a2,b2}&&{c1,d1}==={c2,d2},atr[SUNTrace[g3.prop[k1-q].g1]]qq2 atr[SUNTrace[g4.prop[-k1].g2]],0]cqq[b1,b2]cqq[d1,d2]
+If[{a1,b1}==={a2,d2}&&{c1,d1}==={c2,b2},-atr[SUNTrace[g3.prop[k1-q].g1.g4.prop[-k1].g2 qq2]],0]cqq[d1,b2]cqq[b1,d2];

(* fourth diagram *)
d5d=If[{a1,b1}==={c2,d2}&&{c1,d1}==={a2,b2},atr[SUNTrace[g4.prop[-k1].g1]]qq2 atr[SUNTrace[g3.prop[k1-q].g2]],0]cqq[b1,d2]cqq[d1,b2]
+If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[g1.g3.prop[k1-q].g2.g4.prop[-k1]qq2]],0]cqq[b1,b2]cqq[d1,d2];


(*(*-----------------------------*)
(* first diagram *)
d5a=If[{a1,b1}==={a2,b2}&&{c1,d1}==={c2,d2},atr[SUNTrace[g1 . prop[k1+q] . g3]]qq2 atr[SUNTrace[g2 . prop[-k1] . g4]],0]cqq[a1,a2]cqq[c1,c2]+If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[g1 . prop[k1+q] . g3 . g2 . prop[-k1] . g4 qq2]],0]cqq[a1,c2]cqq[c1,a2];


(*-----------------------------*)
(* second diagram *)
d5b=If[{a1,b1}==={c2,d2}&&{c1,d1}==={a2,b2},atr[SUNTrace[g1 . prop[-k1] . g4]]qq2 atr[SUNTrace[g2 . prop[k1+q] . g3]],0]cqq[a1,c2]cqq[c1,a2]+If[b1===d2&&c2===c1&&d1===b2&&a2===a1,-atr[SUNTrace[g1 . prop[-k1] . g4 . g2 . prop[k1+q] . g3 qq2]],0]cqq[a1,a2]cqq[c1,c2];


(*-----------------------------*)
(* third diagram *)
d5c=If[{a1,b1}==={a2,b2}&&{c1,d1}==={c2,d2},atr[SUNTrace[g3 . prop[k1-q] . g1]]qq2 atr[SUNTrace[g4 . prop[-k1] . g2]],0]cqq[b1,b2]cqq[d1,d2]+If[{a1,b1}==={a2,d2}&&{c1,d1}==={c2,b2},-atr[SUNTrace[g3 . prop[k1-q] . g1 . g4 . prop[-k1] . g2 qq2]],0]cqq[d1,b2]cqq[b1,d2];


(*-----------------------------*)
(* fourth diagram *)
d5d=If[{a1,b1}==={c2,d2}&&{c1,d1}==={a2,b2},atr[SUNTrace[g4 . prop[-k1] . g1]]qq2 atr[SUNTrace[g3 . prop[k1-q] . g2]],0]cqq[b1,d2]cqq[d1,b2]+If[b1===b2&&a2===c1&&d1===d2&&c2===a1,-atr[SUNTrace[g1 . g3 . prop[k1-q] . g2 . g4 . prop[-k1]qq2]],0]cqq[b1,b2]cqq[d1,d2];

*)
(*-----------------------------*)
dia5=IntegrateP[d5a+d5b+d5c+d5d,k1]//SUNSimplify;
QEvaluate[DD5 I dia5 ScaleMu^(4-D),q,Parallelized->False]/.q->qq

]
(*-------------------------------------------------*)
(*-------------------------------------------------*)*)


End[]
(*EndPackage[]*)
