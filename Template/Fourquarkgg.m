(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)



(* ::Code::Initialization::Plain:: *)
Fourquarkgg::usage="Fourquarkgg[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}}] give the \[LeftAngleBracket]GG\[RightAngleBracket] contribution of OPE result about <J1 J2>, with J1=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a1\)]\)\!\(\*SubscriptBox[\(v1\[CapitalPsi]\), \(b1\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c1\)]\)\!\(\*SubscriptBox[\(v2\[CapitalPsi]\), \(d1\)]\) and J2=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(b2\)]\)\!\(\*SubscriptBox[\(v3\[CapitalPsi]\), \(a2\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(d2\)]\)\!\(\*SubscriptBox[\(v4\[CapitalPsi]\), \(c2\)]\) "

Fourquarkgg::inderr="Dummy indices conflict!"

Begin["`Private`Fourquarkgg`"]



(* ::Code::Initialization::Plain:: *)
Options[Fourquarkgg]={
	Parallelized->False,
	AutoNDR->True,
	HoldFlavor->False,
	ToD4->True
}



(* ::Code::Initialization::Plain:: *)
atr[gs_,ndr_,setd_]:=If[ndr,TR5[gs,ToD4->setd],TR[gs]](* How to evaluate the trace *)


(* ::Code::Initialization::Plain:: *)
Fourquarkgg[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,dot,tmp,tmp1,tmp2,hv3,hv4,holdf=OptionValue[HoldFlavor],ndr=OptionValue[AutoNDR],diagrams,files,fdir,pall,setd=OptionValue[ToD4]},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[Fourquarkgg::inderr];
	Abort[]
];

(*--------------------*)


(* B A^+ B *)
hv3=v3//ComplexConjugate;
hv4=null v4//ComplexConjugate//FCI;
(* ComplexConjugate[sigma] will expand DiracSigma, recover it to DiracSigma *)
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


diagrams={type1,type2,type3,type4,type5,type6};


If[pall===True,

	DistributeDefinitions[v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
	tmp= Plus@@WaitAll[ParallelSubmit[{hv3,hv4,holdf,ndr,setd},
									#[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,ndr,setd]]&/@diagrams];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,

	If[pall==="External",
	
		DistributeDefinitions[qq,v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
		
		If[fdir==="None",
		
			ParallelSubmit[{hv3,hv4,holdf,ndr,setd},#[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,ndr,setd]]&/@diagrams
		,
		
			files=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files=("Fourquarkgg_"<>#)&/@files;
			
			ImExport[fdir,
					files,
					{{hv3,hv4,holdf,ndr,setd},
					{qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,ndr,setd},
					diagrams}
					]
		]
		(* !!! do not touch tmp before WaitAll[tmp], any function may modifiy the things in EvaluationObject[...] should not be used here. !!! *)		
	,
	
		tmp= Plus@@(#[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,ndr,setd]&/@diagrams);			
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	
		
	]
]


]


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(* propagators *)

(*xprop[x_]=FourierPX[prop[q],{q,x}];*)
xprop[x_]=1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];



(* ::Code::Initialization::Plain:: *)
(* x ---<--- 0 *)
(*xprog[x_,lora_,lorb_,cola_]=FourierPX[prog[q,lora,lorb,cola],{q,x}];*)
xprog[x_,cola_,lora_,lorb_]=(GAD[lora] . GAD[lorb] . GSD[x] qfact1[1/32 I^-D \[Pi]^(-D/2) qGamma[-1+D/2]] SPD[x,x]^(1-D/2) SUNT[cola]
							+GSD[x] . GAD[lorb] . GAD[lora] qfact1[qfact2[-(1/32) I^-D \[Pi]^(-D/2)] qGamma[-1+D/2]] SPD[x,x]^(1-D/2) SUNT[cola]
							+2 FVD[x,lora] FVD[x,lorb] GSD[x] qfact1[-(1/16) I^-D \[Pi]^(-D/2) qGamma[D/2]] SPD[x,x]^(-D/2) SUNT[cola]);
							
							
(* two gluons in one fermion propagator is vanished after vacuum averaged *)
(* progg[q_,lora_,lorb_,lorc_,lord_,cola_,colb_]=-I/4FAD[{q}]GSD[q].GAD[lorb].SUNT[cola].SUNT[colb].FourDivergence[FAD[{q}]GSD[q].GAD[lord].FourDivergence[GSD[q]FAD[{q}],FVD[q,lorc]],FVD[q,lora]];
(progg[q,lora,lorb,lorc,lord,cola,colb]SUNDelta[cola,colb](MTD[lora,lorc]MTD[lorb,lord]-MTD[lora,lord]MTD[lorb,lorc])//Contract//SUNSimplify//DiracSimplify)//FCI;

%//.{FeynAmpDenominator[PropagatorDenominator[Momentum[q,D],0]]->1/SPD[q],FeynAmpDenominator[PropagatorDenominator[Momentum[q,D],0],p__]->1/SPD[q]FeynAmpDenominator[p]};
%//FCI//Simplify
=0 *)


(* ::Code::Initialization::Plain:: *)
(* \[LeftAngleBracket]\!\(
\*SubsuperscriptBox[\(G\), \(\[Mu]\[Nu]\), \(a\)]
\*SubsuperscriptBox[\(G\), \(\[Alpha]\[Beta]\), \(b\)]\)\[RightAngleBracket] *)
congg[lora_,lorb_,lorc_,lord_]=1/(2 CA CF D(D-1))(MTD[lora,lorc]MTD[lorb,lord]-MTD[lora,lord]MTD[lorb,lorc]);
(* the appear Condensate["gg"] = g^2<G^n_uv G^n uv> by the definition of the diagrams below, times 4Pi so that Condensate["gg"] = \alpha_s <G^n_uv G^n uv> *)



(* ::Code::Initialization::Plain:: *)
(* -------------The following are generated by algorithem----------------- *)


(* ::Input::Initialization::Plain:: *)
type1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,setd_]:=Block[{x,q,dia1,sun,lora,lorb,lorc,lord,str,tr,contract},


dia1=Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv3 . xprog[-x, sun, lorc, lord] . v2 . xprog[x, sun, lora, lorb]]]*tr[str[hv4 . xprop[-x] . v1 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv3 . xprog[-x, sun, lorc, lord] . v1 . xprop[x] . hv4 . xprog[-x, sun, lora, lorb] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv4 . xprop[-x] . v1 . xprog[x, sun, lorc, lord] . hv3 . xprop[-x] . v2 . xprog[x, sun, lora, lorb]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["gg"]4Pi*contract[congg[lora, lorb, lorc, lord]*tr[str[hv3 . xprog[-x, sun, lorc, lord] . v1 . xprop[x]]]*tr[str[hv4 . xprog[-x, sun, lora, lorb] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia1=dia1/.FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]/.str->SUNTrace/.{tr[gs_]:>atr[gs,ndr,setd],contract->Contract};


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,setd_]:=Block[{x,q,dia2,sun,lora,lorb,lorc,lord,str,tr,contract},


dia2=Condensate["gg"]4Pi*contract[congg[lora, lorb, lorc, lord]*tr[str[hv3 . xprog[-x, sun, lora, lorb] . v2 . xprop[x]]]*tr[str[hv4 . xprog[-x, sun, lorc, lord] . v1 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["gg"]4Pi*contract[congg[lora, lorb, lorc, lord]*tr[str[hv3 . xprog[-x, sun, lora, lorb] . v2 . xprop[x]]]*tr[str[hv4 . xprop[-x] . v1 . xprog[x, sun, lorc, lord]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv3 . xprog[-x, sun, lorc, lord] . v1 . xprog[x, sun, lora, lorb] . hv4 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv3 . xprog[-x, sun, lorc, lord] . v1 . xprop[x] . hv4 . xprop[-x] . v2 . xprog[x, sun, lora, lorb]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia2=dia2/.FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]/.str->SUNTrace/.{tr[gs_]:>atr[gs,ndr,setd],contract->Contract};


dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,setd_]:=Block[{x,q,dia3,sun,lora,lorb,lorc,lord,str,tr,contract},


dia3=Condensate["gg"]4Pi*contract[tr[str[hv3 . xprop[-x] . v2 . xprop[x]]]*tr[congg[lora, lorb, lorc, lord]*str[hv4 . xprog[-x, sun, lorc, lord] . v1 . xprog[x, sun, lora, lorb]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["gg"]4Pi*contract[congg[lora, lorb, lorc, lord]*tr[str[hv3 . xprop[-x] . v2 . xprog[x, sun, lora, lorb]]]*tr[str[hv4 . xprop[-x] . v1 . xprog[x, sun, lorc, lord]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv4 . xprog[-x, sun, lorc, lord] . v1 . xprop[x] . hv3 . xprop[-x] . v2 . xprog[x, sun, lora, lorb]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["gg"]4Pi*contract[congg[lora, lorb, lorc, lord]*tr[str[hv3 . xprop[-x] . v1 . xprog[x, sun, lorc, lord]]]*tr[str[hv4 . xprog[-x, sun, lora, lorb] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia3=dia3/.FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]/.str->SUNTrace/.{tr[gs_]:>atr[gs,ndr,setd],contract->Contract};


dia3=FourierXP[dia3,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=QEvaluate[I ScaleMu^(3(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,setd_]:=Block[{x,q,dia4,sun,lora,lorb,lorc,lord,str,tr,contract},


dia4=Condensate["gg"]4Pi*contract[congg[lora, lorb, lorc, lord]*tr[str[hv3 . xprop[-x] . v2 . xprog[x, sun, lora, lorb]]]*tr[str[hv4 . xprog[-x, sun, lorc, lord] . v1 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . xprog[-x, sun, lorc, lord] . v2 . xprog[x, sun, lora, lorb]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["gg"]4Pi*contract[tr[str[hv3 . xprop[-x] . v1 . xprop[x]]]*tr[congg[lora, lorb, lorc, lord]*str[hv4 . xprog[-x, sun, lorc, lord] . v2 . xprog[x, sun, lora, lorb]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["gg"]4Pi*contract[congg[lora, lorb, lorc, lord]*tr[str[hv3 . xprog[-x, sun, lorc, lord] . v1 . xprop[x]]]*tr[str[hv4 . xprop[-x] . v2 . xprog[x, sun, lora, lorb]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia4=dia4/.FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]/.str->SUNTrace/.{tr[gs_]:>atr[gs,ndr,setd],contract->Contract};


dia4=FourierXP[dia4,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=QEvaluate[I ScaleMu^(3(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,setd_]:=Block[{x,q,dia5,sun,lora,lorb,lorc,lord,str,tr,contract},


dia5=-(Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv3 . xprop[-x] . v1 . xprog[x, sun, lorc, lord] . hv4 . xprop[-x] . v2 . xprog[x, sun, lora, lorb]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv4 . xprog[-x, sun, lorc, lord] . v1 . xprop[x] . hv3 . xprog[-x, sun, lora, lorb] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["gg"]4Pi*contract[congg[lora, lorb, lorc, lord]*tr[str[hv3 . xprop[-x] . v1 . xprog[x, sun, lorc, lord]]]*tr[str[hv4 . xprop[-x] . v2 . xprog[x, sun, lora, lorb]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv3 . xprog[-x, sun, lorc, lord] . v1 . xprog[x, sun, lora, lorb]]]*tr[str[hv4 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia5=dia5/.FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]/.str->SUNTrace/.{tr[gs_]:>atr[gs,ndr,setd],contract->Contract};


dia5=FourierXP[dia5,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=QEvaluate[I ScaleMu^(3(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,setd_]:=Block[{x,q,dia6,sun,lora,lorb,lorc,lord,str,tr,contract},


dia6=-(Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv3 . xprop[-x] . v1 . xprog[x, sun, lorc, lord] . hv4 . xprog[-x, sun, lora, lorb] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . xprog[-x, sun, lorc, lord] . v2 . xprog[x, sun, lora, lorb]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv4 . xprog[-x, sun, lorc, lord] . v1 . xprog[x, sun, lora, lorb] . hv3 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] - Condensate["gg"]4Pi*contract[tr[congg[lora, lorb, lorc, lord]*str[hv4 . xprop[-x] . v1 . xprog[x, sun, lorc, lord] . hv3 . xprog[-x, sun, lora, lorb] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia6=dia6/.FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]/.str->SUNTrace/.{tr[gs_]:>atr[gs,ndr,setd],contract->Contract};


dia6=FourierXP[dia6,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia6=QEvaluate[I ScaleMu^(3(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_1 *)
type1[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,setd_]:=Block[{tmp0,x,q,lora,lorb,lorc,lord,cola,colb,k1,k2,k3,l,tra1,tra2,dia1},

(*-----------------------------*)
(* first diagram *)
tra1=(FlavorDelta[a1,a2]FlavorDelta[b1,b2]FlavorDelta[c1,c2]FlavorDelta[d1,d2]atr[SUNTrace[g1 . xprog[x,lora,lorb,cola] . g3 . xprog[-x,lorc,lord,colb]GG[lora,lorb,lorc,lord,cola,colb]]]atr[SUNTrace[g2 . xprop[x] . g4 . xprop[-x]]]
+
(* \[CapitalGamma]_3\[TwoWayRule]\[CapitalGamma]_4 *)
FlavorDelta[a1,c2]FlavorDelta[b1,d2]FlavorDelta[c1,a2]FlavorDelta[d1,b2]atr[SUNTrace[g1 . xprog[x,lora,lorb,cola] . g4 . xprog[-x,lorc,lord,colb]GG[lora,lorb,lorc,lord,cola,colb]]]atr[SUNTrace[g2 . xprop[x] . g3 . xprop[-x]]]);

(*-----------------------------*)
(* second diagram *)
tra2=(FlavorDelta[a1,a2]FlavorDelta[b1,b2]FlavorDelta[c1,c2]FlavorDelta[d1,d2]atr[SUNTrace[g1 . xprop[x] . g3 . xprop[-x]]]atr[SUNTrace[g2 . xprog[x,lora,lorb,cola] . g4 . xprog[-x,lorc,lord,colb]GG[lora,lorb,lorc,lord,cola,colb]]]
+
(* \[CapitalGamma]_3\[TwoWayRule]\[CapitalGamma]_4 *)
FlavorDelta[a1,c2]FlavorDelta[b1,d2]FlavorDelta[c1,a2]FlavorDelta[d1,b2]atr[SUNTrace[g1 . xprop[x] . g4 . xprop[-x]]]atr[SUNTrace[g2 . xprog[x,lora,lorb,cola] . g3 . xprog[-x,lorc,lord,colb]GG[lora,lorb,lorc,lord,cola,colb]]]);


(*-----------------------------*)
dia1=Condensate["gg"]4PiFourierXP[tra1+tra2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};
QEvaluate[I ScaleMu^(3(4-D))dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq

]
(*-------------------------------------------------*)
(*-------------------------------------------------*)


*)


(* ::Code::Initialization::Plain:: *)
(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_2 *)
type2[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,setd_]:=Block[{tmp0,x,q,lora,lorb,lorc,lord,cola,colb,k1,k2,k3,l,tra21,tra22,tra23,tra24,dia2},

(*-----------------------------*)
(* first diagram *)
tra21=(FlavorDelta[a1,a2]FlavorDelta[b1,b2]FlavorDelta[c1,c2]FlavorDelta[d1,d2]atr[SUNTrace[g1 . xprog[x,lora,lorb,cola] . g3 . xprop[-x]GG[lora,lorb,lorc,lord,cola,colb]]]atr[SUNTrace[g2 . xprog[x,lorc,lord,colb] . g4 . xprop[-x]]]
+
(* \[CapitalGamma]_3\[TwoWayRule]\[CapitalGamma]_4 *)
FlavorDelta[a1,c2]FlavorDelta[b1,d2]FlavorDelta[c1,a2]FlavorDelta[d1,b2]atr[SUNTrace[g1 . xprog[x,lora,lorb,cola] . g4 . xprop[-x]GG[lora,lorb,lorc,lord,cola,colb]]]atr[SUNTrace[g2 . xprog[x,lorc,lord,colb] . g3 . xprop[-x]]]);


(*-----------------------------*)
(* second diagram *)
tra22=(FlavorDelta[a1,a2]FlavorDelta[b1,b2]FlavorDelta[c1,c2]FlavorDelta[d1,d2]atr[SUNTrace[g1 . xprog[x,lora,lorb,cola] . g3 . xprop[-x]GG[lora,lorb,lorc,lord,cola,colb]]]atr[SUNTrace[g2 . xprop[x] . g4 . xprog[-x,lorc,lord,colb]]]
+
(* \[CapitalGamma]_3\[TwoWayRule]\[CapitalGamma]_4 *)
FlavorDelta[a1,c2]FlavorDelta[b1,d2]FlavorDelta[c1,a2]FlavorDelta[d1,b2]atr[SUNTrace[g1 . xprog[x,lora,lorb,cola] . g4 . xprop[-x]GG[lora,lorb,lorc,lord,cola,colb]]]atr[SUNTrace[g2 . xprop[x] . g3 . xprog[-x,lorc,lord,colb]]]);



(*-----------------------------*)
(* third diagram *)
tra23=(FlavorDelta[a1,a2]FlavorDelta[b1,b2]FlavorDelta[c1,c2]FlavorDelta[d1,d2]atr[SUNTrace[g1 . xprop[x] . g3 . xprog[-x,lora,lorb,cola]GG[lora,lorb,lorc,lord,cola,colb]]]atr[SUNTrace[g2 . xprog[x,lorc,lord,colb] . g4 . xprop[-x]]]
+
(* \[CapitalGamma]_3\[TwoWayRule]\[CapitalGamma]_4 *)
FlavorDelta[a1,c2]FlavorDelta[b1,d2]FlavorDelta[c1,a2]FlavorDelta[d1,b2]atr[SUNTrace[g1 . xprop[x] . g4 . xprog[-x,lora,lorb,cola]GG[lora,lorb,lorc,lord,cola,colb]]]atr[SUNTrace[g2 . xprog[x,lorc,lord,colb] . g3 . xprop[-x]]]);


(*-----------------------------*)
(* fourth diagram *)
tra24=(FlavorDelta[a1,c2]FlavorDelta[b1,d2]FlavorDelta[c1,a2]FlavorDelta[d1,b2]atr[SUNTrace[g1 . xprop[x] . g3 . xprog[-x,lora,lorb,cola]GG[lora,lorb,lorc,lord,cola,colb]]]atr[SUNTrace[g2 . xprop[x] . g4 . xprog[-x,lorc,lord,colb]]]
+
(* \[CapitalGamma]_3\[TwoWayRule]\[CapitalGamma]_4 *)
FlavorDelta[a1,c2]FlavorDelta[b1,d2]FlavorDelta[c1,a2]FlavorDelta[d1,b2]atr[SUNTrace[g1 . xprop[x] . g4 . xprog[-x,lora,lorb,cola]GG[lora,lorb,lorc,lord,cola,colb]]]atr[SUNTrace[g2 . xprop[x] . g3 . xprog[-x,lorc,lord,colb]]]);





dia2=Condensate["gg"]4PiFourierXP[tra21+tra22+tra23+tra24,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};
QEvaluate[I ScaleMu^(3(4-D))dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]
(*-------------------------------------------------*)
(*-------------------------------------------------*)

*)


(* ::Code::Initialization::Plain:: *)
(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_3_1 *)
type31[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,setd_]:=Block[{tmp0,x,q,lora,lorb,lorc,lord,cola,colb,k1,k2,k3,l,trb1,trb2,trb3,trb4,trb5,trb6,dib1},


(*-----------------------------*)
(*  first diagram *)
trb1=(FlavorDelta[b1,b2]FlavorDelta[a2,c1]FlavorDelta[d1,d2]FlavorDelta[c2,a1]atr[SUNTrace[g1 . xprog[x,lora,lorb,cola] . g3 . xprop[-x] . g2 . xprop[x] . g4 . xprog[-x,lorc,lord,colb]GG[lora,lorb,lorc,lord,cola,colb]]]
+
(* \[CapitalGamma]_3\[TwoWayRule]\[CapitalGamma]_4 *)
FlavorDelta[b1,d2]FlavorDelta[c2,c1]FlavorDelta[d1,b2]FlavorDelta[a2,a1]atr[SUNTrace[g1 . xprog[x,lora,lorb,cola] . g4 . xprop[-x] . g2 . xprop[x] . g3 . xprog[-x,lorc,lord,colb]GG[lora,lorb,lorc,lord,cola,colb]]]
);


(*-----------------------------*)
(*  second diagram *)
trb2=(FlavorDelta[b1,b2]FlavorDelta[a2,c1]FlavorDelta[d1,d2]FlavorDelta[c2,a1]atr[SUNTrace[g1 . xprog[x,lora,lorb,cola] . g3 . xprog[-x,lorc,lord,colb] . g2 . xprop[x] . g4 . xprop[-x]GG[lora,lorb,lorc,lord,cola,colb]]]
+
(* \[CapitalGamma]_3\[TwoWayRule]\[CapitalGamma]_4 *)
FlavorDelta[b1,d2]FlavorDelta[c2,c1]FlavorDelta[d1,b2]FlavorDelta[a2,a1]atr[SUNTrace[g1 . xprog[x,lora,lorb,cola] . g4 . xprog[-x,lorc,lord,colb] . g2 . xprop[x] . g3 . xprop[-x]GG[lora,lorb,lorc,lord,cola,colb]]]
);



(*-----------------------------*)
dib1=Condensate["gg"]4PiFourierXP[trb1+trb2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};
QEvaluate[-I ScaleMu^(3(4-D))dib1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
(* odd number of fermion loops *)



]
(*-------------------------------------------------*)
(*-------------------------------------------------*)


*)


(* ::Code::Initialization::Plain:: *)
(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_3_2 *)
type32[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,setd_]:=Block[{tmp0,x,q,lora,lorb,lorc,lord,cola,colb,k1,k2,k3,l,trb1,trb2,trb3,trb4,trb5,trb6,dib1},

(*-----------------------------*)
(*  third diagram *)
trb3=(FlavorDelta[b1,b2]FlavorDelta[a2,c1]FlavorDelta[d1,d2]FlavorDelta[c2,a1]atr[SUNTrace[g1 . xprog[x,lora,lorb,cola] . g3 . xprop[-x] . g2 . xprog[x,lorc,lord,colb] . g4 . xprop[-x]GG[lora,lorb,lorc,lord,cola,colb]]]
+
(* \[CapitalGamma]_3\[TwoWayRule]\[CapitalGamma]_4 *)
FlavorDelta[b1,d2]FlavorDelta[c2,c1]FlavorDelta[d1,b2]FlavorDelta[a2,a1]atr[SUNTrace[g1 . xprog[x,lora,lorb,cola] . g4 . xprop[-x] . g2 . xprog[x,lorc,lord,colb] . g3 . xprop[-x]GG[lora,lorb,lorc,lord,cola,colb]]]
);


(*-----------------------------*)
(*  4th diagram *)
trb4=(FlavorDelta[b1,b2]FlavorDelta[a2,c1]FlavorDelta[d1,d2]FlavorDelta[c2,a1]atr[SUNTrace[g1 . xprop[x] . g3 . xprog[-x,lora,lorb,cola] . g2 . xprop[x] . g4 . xprog[-x,lorc,lord,colb]GG[lora,lorb,lorc,lord,cola,colb]]]
+
(* \[CapitalGamma]_3\[TwoWayRule]\[CapitalGamma]_4 *)
FlavorDelta[b1,d2]FlavorDelta[c2,c1]FlavorDelta[d1,b2]FlavorDelta[a2,a1]atr[SUNTrace[g1 . xprop[x] . g4 . xprog[-x,lora,lorb,cola] . g2 . xprop[x] . g3 . xprog[-x,lorc,lord,colb]GG[lora,lorb,lorc,lord,cola,colb]]]
);



(*-----------------------------*)
dib1=Condensate["gg"]4PiFourierXP[trb3+trb4,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};
QEvaluate[-I ScaleMu^(3(4-D))dib1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
(* odd number of fermion loops *)

]
(*-------------------------------------------------*)
(*-------------------------------------------------*)


*)


(* ::Code::Initialization::Plain:: *)
(*(*-------------------------------------------------*)
(*-------------------------------------------------*)
(* type_3_3 *)
type33[qq_,{g1_,g2_,{a1_,b1_,c1_,d1_}},{g3_,g4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,setd_]:=Block[{tmp0,x,q,lora,lorb,lorc,lord,cola,colb,k1,k2,k3,l,trb1,trb2,trb3,trb4,trb5,trb6,dib1},


(*-----------------------------*)
(*  5th diagram *)
trb5=(FlavorDelta[b1,b2]FlavorDelta[a2,c1]FlavorDelta[d1,d2]FlavorDelta[c2,a1]atr[SUNTrace[g1 . xprop[x] . g3 . xprop[-x] . g2 . xprog[x,lorc,lord,colb] . g4 . xprog[-x,lora,lorb,cola]GG[lora,lorb,lorc,lord,cola,colb]]]
+
(* \[CapitalGamma]_3\[TwoWayRule]\[CapitalGamma]_4 *)
FlavorDelta[b1,d2]FlavorDelta[c2,c1]FlavorDelta[d1,b2]FlavorDelta[a2,a1]atr[SUNTrace[g1 . xprop[x] . g4 . xprop[-x] . g2 . xprog[x,lorc,lord,colb] . g3 . xprog[-x,lora,lorb,cola]GG[lora,lorb,lorc,lord,cola,colb]]]
);


(*-----------------------------*)
(*  6th diagram *)
trb6=(FlavorDelta[b1,b2]FlavorDelta[a2,c1]FlavorDelta[d1,d2]FlavorDelta[c2,a1]atr[SUNTrace[g1 . xprop[x] . g3 . xprog[-x,lorc,lord,colb] . g2 . xprog[x,lora,lorb,cola] . g4 . xprop[-x]GG[lora,lorb,lorc,lord,cola,colb]]]
+
(* \[CapitalGamma]_3\[TwoWayRule]\[CapitalGamma]_4 *)
FlavorDelta[b1,d2]FlavorDelta[c2,c1]FlavorDelta[d1,b2]FlavorDelta[a2,a1]atr[SUNTrace[g1 . xprop[x] . g4 . xprog[-x,lorc,lord,colb] . g2 . xprog[x,lora,lorb,cola] . g3 . xprop[-x]GG[lora,lorb,lorc,lord,cola,colb]]]
);


(*-----------------------------*)
dib1=Condensate["gg"]4PiFourierXP[trb5+trb6,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};
QEvaluate[-I ScaleMu^(3(4-D))dib1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
(* odd number of fermion loops *)



]
(*-------------------------------------------------*)
(*-------------------------------------------------*)


*)


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
