(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
FourquarkAlphasType1::usage="FourquarkAlphasType1[q_,j1_,j2_] give the next leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents.
Here the gluon propagator connect with two quark propagators, the corresponding renormalization involve tetraquark currents with same flavor structure."

FourquarkAlphasType1::inderr="Dummy indices conflict!"

FourquarkAlphasType1::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`FourquarkAlphasType1`"]


(* ::Code::Initialization::Plain:: *)
Options[FourquarkAlphasType1] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->False,
	Renormalization->True,
	Strategy->"Fourier"
}


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

prop[q_]=I GSD[q]FAD[q];

(*xprop[x_]=FourierPX[prop[q],{q,x}];*)
xprop[x_]=1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];
			


FourquarkAlphasType1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,atr,hv21,hv22,listpole1,listpole2,holdf=OptionValue[HoldFlavor],
strategy=OptionValue[Strategy],diagrams,ndr=OptionValue[AutoNDR],ren=OptionValue[Renormalization],fdir,files1,pall1,files2,pall2,pall,dot},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[FourquarkAlphasType1::inderr];
	Abort[]
];


(*-------------------------------*)
(* B A^+ B *)
hv21=v3//ComplexConjugate//FCI;
hv22=v4//ComplexConjugate//FCI;
(* ComplexConjugate[sigma] will expand the DiracSigma to gamma matrices, recover it to DiracSigma *)
hv22=(hv22/.Dot->dot)/.f1_ dot[aa___,DiracGamma[LorentzIndex[lor1_,dim___],dim___],DiracGamma[LorentzIndex[lor2_,dim___],dim___],bb___]+f2_ dot[aa___,DiracGamma[LorentzIndex[lor2_,dim___],dim___],DiracGamma[LorentzIndex[lor1_,dim___],dim___],bb___]/;(f1+f2===0):>-2I f1 dot[aa,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],bb];
hv22=hv22/.{null->1,dot->Dot};

(*-------------------------------*)

If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];


(*(*---------------------------------------------------*)

If[OptionValue[Renormalization]===True,
	
	listpole1={#[[1]], {#[[2,1]],#[[2,2]],{a1,b1,c1,d1}}}&/@FourquarkPole[v1,v2];
	listpole2={ComplexConjugate[#[[1]]],{#[[2,1]],#[[2,2]],{a2,b2,c2,d2}}}&/@FourquarkPole[v3,v4];

	reno=(Plus@@(Contract[#[[1]]TetraquarkLeading[qq,#[[2]],{v3,v4,{a2,b2,c2,d2}},EpsOrder->1,AutoNDR->OptionValue[AutoNDR],HoldFlavor->holdf]]&/@listpole1))+
			Plus@@(Contract[#[[1]]TetraquarkLeading[qq,{v1,v2,{a1,b1,c1,d1}},#[[2]],EpsOrder->1,AutoNDR->OptionValue[AutoNDR],HoldFlavor->holdf]]&/@listpole2);			
];
*)
(*---------------------------------------------------*)


If[ToLowerCase[ToString[strategy]]=="fourier",
	diagrams={xtype1,xtype2,xtype3,xtype4,xtype5,xtype6,xtype7,xtype8,xtype9,xtype10,xtype11,xtype12}
,
	If[ToLowerCase[ToString[strategy]]=="integratefirst",
		diagrams={iptype1,iptype2,iptype3,iptype4,iptype5,iptype6,iptype7,iptype8,iptype9,iptype10,iptype11,iptype12}
	,
		diagrams={ptype1,ptype2,ptype3,ptype4,ptype5,ptype6,ptype7,ptype8}
	]
];

(*---------------------------------------------------*)

If[OptionValue[AutoNDR]===True,atr=TR5,atr=TR];

(*---------------------------------------------------*)
If[pall===True,

	DistributeDefinitions[v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
	
	tmp= Plus@@WaitAll[Join[ParallelSubmit[{hv21,hv22,holdf,atr},
									#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr]]&/@diagrams,
				ParallelSubmit[{ndr,holdf,ren},-#[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},ndr,holdf,ren]]&/@{reno1,reno2}
				]];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,
	If[pall==="External",

		DistributeDefinitions[v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];

		If[fdir=="None",
		(* evaluation, no import and export *)
			Join[ParallelSubmit[{hv21,hv22,holdf,atr},
									#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr]]&/@diagrams,
				ParallelSubmit[{ndr,holdf,ren},-#[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},ndr,holdf,ren]]&/@{reno1,reno2}
				]

		,
		(* evaluation, import and export *)
			files1=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files1=("FourquarkAlphasType1_"<>#)&/@files1;
			
			
			pall1=ImExport[fdir,
						files1,
						{{hv21,hv22,holdf,atr},
						{qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr},
						diagrams}
						];
		
		
			files2=(StringSplit[ToString[#],"`"][[-1]])&/@{reno1,reno2};
			files2=("FourquarkAlphasType1_"<>#)&/@files2;
			
			If[ren===True,
				pall2=ImExport[fdir,
						files2,
						{{ndr,holdf,ren},
						{qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},ndr,holdf,ren},
						{reno1,reno2}}
						]
			,
				pall2={{2,2},{0,0}}			
			];


			{pall1[[1]]+pall2[[1]],Join[pall1[[2]],-pall2[[2]]]}		
		](* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
	
	,
	
		
		tmp= Plus@@(-#[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},ndr,holdf,ren]&/@{reno1,reno2})+ Plus@@Map[#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr]&,diagrams];
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]

]

(*------------------------------------------------------------------*)


(* ::Code::Initialization::Plain:: *)
(*reno[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ndr_,holdf_,True]:=Block[{tmp,listpole1,listpole2},
	
	listpole1={#[[1]], {#[[2]],#[[3]],{a1,b1,c1,d1}}}&/@FourquarkPole[v1,v2];
	listpole2={ComplexConjugate[#[[1]]],{#[[2]],#[[3]],{a2,b2,c2,d2}}}&/@FourquarkPole[v3,v4];
	
	
	QGather[(Plus@@(Contract[FourquarkLeading[qq,#[[2]],{v3,v4,{a2,b2,c2,d2}},Pole->#[[1]],EpsOrder->0,HoldFlavor->holdf,AutoNDR->ndr]]&/@listpole1))+
			Plus@@(Contract[FourquarkLeading[qq,{v1,v2,{a1,b1,c1,d1}},#[[2]],Pole->#[[1]],EpsOrder->0,HoldFlavor->holdf,AutoNDR->ndr]]&/@listpole2),qq,ShowasTable->False]/.{CA-2CF->1/CA,2CA CF-CA^2+1->0}/.CF->(CA^2-1)/(2CA)
					
]*)


(* ::Code::Initialization::Plain:: *)
reno1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ndr_,holdf_,False]=0
reno2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ndr_,holdf_,False]=0


(* ::Code::Initialization::Plain:: *)
reno1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ndr_,holdf_,True]:=Block[{tmp,listpole1},
	tmp=FourquarkCurrent[Current[a1,v1,b1],Current[c1,v2,d1]];
	
	tmp=FourquarkPole[tmp,ShowasTable->True,VertexOnly->False,ReduceSUN->True];
	listpole1={#[[1]], {#[[2,2]],#[[3,2]],{#[[2,1]],#[[2,3]],#[[3,1]],#[[3,3]]}}}&/@tmp;
		(*FourquarkPole[v1,v2]*)
		(*listpole1={#[[1]], {#[[2]],#[[3]],{a1,b1,c1,d1}}}&/@tmp*)
	
	Plus@@(FourquarkLeading[qq,#[[2]],{v3,v4,{a2,b2,c2,d2}},Pole->#[[1]],EpsOrder->0,HoldFlavor->holdf,AutoNDR->ndr]&/@listpole1)/.{CA-2CF->1/CA,2CA CF-CA^2+1->0}/.CF->(CA^2-1)/(2CA)					
]


(* ::Code::Initialization::Plain:: *)
reno2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ndr_,holdf_,True]:=Block[{tmp,listpole2},
	tmp=FourquarkCurrent[Current[a2,v3,b2],Current[c2,v4,d2]];	

	tmp=FourquarkPole[tmp,ShowasTable->True,VertexOnly->False,ReduceSUN->True];
	listpole2={ComplexConjugate[#[[1]]], {#[[2,2]],#[[3,2]],{#[[2,1]],#[[2,3]],#[[3,1]],#[[3,3]]}}}&/@tmp;
	(*FourquarkPole[v1,v2]*)

	Plus@@(FourquarkLeading[qq,{v1,v2,{a1,b1,c1,d1}},#[[2]],Pole->#[[1]],EpsOrder->0,HoldFlavor->holdf,AutoNDR->ndr]&/@listpole2)/.{CA-2CF->1/CA,2CA CF-CA^2+1->0}/.CF->(CA^2-1)/(2CA)				
]


(* ::Code::Initialization::Plain:: *)
(* The following are generated by algorithem *)


(* ::Code::Initialization::Plain:: *)
(*---------------------------------------------------------------------*)


(* ::Input::Initialization:: *)
xtype1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia1=-(FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x], hv4, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]) + contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia2=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia3= contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x]]]]*tr[str[dot[hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2], hv3, xprop[-x], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia3=FourierXP[dia3,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=QEvaluate[I ScaleMu^(4(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia4=contract[tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x], hv3, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia4=FourierXP[dia4,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=QEvaluate[I ScaleMu^(4(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia5=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]*tr[str[dot[hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv3, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia5=FourierXP[dia5,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=QEvaluate[I ScaleMu^(4(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia6=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia6=FourierXP[dia6,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia6=QEvaluate[I ScaleMu^(4(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia7=contract[tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv4, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia7=FourierXP[dia7,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia7=QEvaluate[I ScaleMu^(4(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia8=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]*tr[str[dot[hv4, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]  - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x], hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia8=FourierXP[dia8,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia8=QEvaluate[I ScaleMu^(4(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia9=-(FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x], hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]])+ contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia9=FourierXP[dia9,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia9=QEvaluate[I ScaleMu^(4(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia10=-(FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]])+ contract[tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia10=FourierXP[dia10,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia10=QEvaluate[I ScaleMu^(4(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia11=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia11=FourierXP[dia11,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia11=QEvaluate[I ScaleMu^(4(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia12=-(FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2], hv4, xprop[-x], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]) + contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]*tr[str[dot[hv4, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia12=FourierXP[dia12,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia12=QEvaluate[I ScaleMu^(4(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*---------------------------------------------------------------------*)


(* ::Input::Initialization:: *)
(*xtype1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia1=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia2=-(FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x], hv4, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]) - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia3=-(FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x], hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]) - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2], hv3, xprop[-x], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia3=FourierXP[dia3,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=QEvaluate[I ScaleMu^(4(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia4=contract[tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x], hv3, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia4=FourierXP[dia4,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=QEvaluate[I ScaleMu^(4(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia5=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]*tr[str[dot[hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv3, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia5=FourierXP[dia5,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=QEvaluate[I ScaleMu^(4(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia6=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia6=FourierXP[dia6,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia6=QEvaluate[I ScaleMu^(4(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia7=contract[tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv4, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia7=FourierXP[dia7,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia7=QEvaluate[I ScaleMu^(4(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia8=-(FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]) - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x], hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia8=FourierXP[dia8,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia8=QEvaluate[I ScaleMu^(4(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia9=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x]]]]*tr[str[dot[hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia9=FourierXP[dia9,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia9=QEvaluate[I ScaleMu^(4(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia10=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]*tr[str[dot[hv4, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + contract[tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia10=FourierXP[dia10,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia10=QEvaluate[I ScaleMu^(4(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia11=-(FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2], hv4, xprop[-x], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]) - FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia11=FourierXP[dia11,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia11=QEvaluate[I ScaleMu^(4(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia12=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]*tr[str[dot[hv4, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia12=FourierXP[dia12,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia12=QEvaluate[I ScaleMu^(4(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Code::Initialization::Plain:: *)
(*---------------------------------------------------------------------*)


(* ::Input::Initialization:: *)
(*xtype1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia1=-(contract[tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x], hv4, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) + contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]*tr[str[dot[hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia2=-(contract[tr[str[dot[hv4, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv3, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) + contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia3=-(contract[tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) + contract[tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia3=FourierXP[dia3,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=QEvaluate[I ScaleMu^(4(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia4=-(contract[tr[str[dot[hv3, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) + contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x]]]]*tr[str[dot[hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia4=FourierXP[dia4,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=QEvaluate[I ScaleMu^(4(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia5=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - contract[tr[str[dot[hv3, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv4, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia5=FourierXP[dia5,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=QEvaluate[I ScaleMu^(4(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia6=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]*tr[str[dot[hv4, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + contract[tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia6=FourierXP[dia6,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia6=QEvaluate[I ScaleMu^(4(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia7=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - contract[tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x], hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia7=FourierXP[dia7,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia7=QEvaluate[I ScaleMu^(4(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia8=-(contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) + contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia8=FourierXP[dia8,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia8=QEvaluate[I ScaleMu^(4(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia9=-(contract[tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x], hv3, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) - contract[tr[str[dot[hv4, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia9=FourierXP[dia9,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia9=QEvaluate[I ScaleMu^(4(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia10=contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - contract[tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, xprop[x], hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia10=FourierXP[dia10,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia10=QEvaluate[I ScaleMu^(4(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia11=-(contract[tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2], hv4, xprop[-x], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) + contract[TLPropagatorX[x, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, xprop[-x], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]*tr[str[dot[hv4, xprop[-x], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia11=FourierXP[dia11,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia11=QEvaluate[I ScaleMu^(4(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,sun,lor1,lor2,lor3,lor4,lorz,str,tr,dot,contract},


dia12=contract[tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - contract[tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2], hv3, xprop[-x], v2, xprop[x]]]*TLPropagatorX[x, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia12=FourierXP[dia12,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia12=QEvaluate[I ScaleMu^(4(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Code::Initialization::Plain:: *)
(*---------------------------------------*)


(* ::Input::Initialization:: *)
(*xtype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,lorz,sun,k1,k2,l,k},


dia3=Contract[atr[SUNTrace[hv4 . GAD[lor4] . GAD[lorz] . SUNT[sun] . GAD[lor3] . v1 . xprop[x]]]*atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv3 . GAD[lor2] . GAD[lorz] . SUNT[sun] . GAD[lor1] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Contract[atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . xprop[x] . hv3 . GAD[lor4] . GAD[lorz] . SUNT[sun] . GAD[lor3] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv4 . GAD[lor2] . GAD[lorz] . SUNT[sun] . GAD[lor1] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia3=FourierXP[dia3,{x,q}]/.SUNN->CA;

dia3=QEvaluate[I ScaleMu^(4(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,lorz,sun,k1,k2,l,k},


dia4=-(Contract[atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4] . hv4 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Contract[atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv4 . GAD[lor4] . GAD[lorz] . SUNT[sun] . GAD[lor3] . v1 . xprop[x] . hv3 . GAD[lor2] . GAD[lorz] . SUNT[sun] . GAD[lor1] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Contract[atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv3 . GAD[lor4] . GAD[lorz] . SUNT[sun] . GAD[lor3] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia4=FourierXP[dia4,{x,q}]/.SUNN->CA;

dia4=QEvaluate[I ScaleMu^(4(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,lorz,sun,k1,k2,l,k},


dia5=Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv3 . GAD[lor2] . GAD[lorz] . SUNT[sun] . GAD[lor1] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Contract[atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . xprop[x] . hv4 . GAD[lor4] . GAD[lorz] . SUNT[sun] . GAD[lor3] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] + Contract[atr[SUNTrace[hv3 . GAD[lor4] . GAD[lorz] . SUNT[sun] . GAD[lor3] . v1 . xprop[x]]]*atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv4 . GAD[lor2] . GAD[lorz] . SUNT[sun] . GAD[lor1] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia5=FourierXP[dia5,{x,q}]/.SUNN->CA;

dia5=QEvaluate[I ScaleMu^(4(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,lorz,sun,k1,k2,l,k},


dia6=Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv3 . GAD[lor4] . GAD[lorz] . SUNT[sun] . GAD[lor3] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Contract[atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv3 . GAD[lor4] . GAD[lorz] . SUNT[sun] . GAD[lor3] . v1 . xprop[x] . hv4 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Contract[atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4] . hv4 . GAD[lor2] . GAD[lorz] . SUNT[sun] . GAD[lor1] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia6=FourierXP[dia6,{x,q}]/.SUNN->CA;

dia6=QEvaluate[I ScaleMu^(4(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,lorz,sun,k1,k2,l,k},


dia7=Contract[atr[SUNTrace[hv4 . GAD[lor4] . GAD[lorz] . SUNT[sun] . GAD[lor3] . v1 . xprop[x]]]*atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Contract[atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv4 . GAD[lor4] . GAD[lorz] . SUNT[sun] . GAD[lor3] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . hv3 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia7=FourierXP[dia7,{x,q}]/.SUNN->CA;

dia7=QEvaluate[I ScaleMu^(4(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,lorz,sun,k1,k2,l,k},


dia8=-(Contract[atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4] . hv3 . GAD[lor2] . GAD[lorz] . SUNT[sun] . GAD[lor1] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) - Contract[atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4] . hv3 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Contract[atr[SUNTrace[hv3 . GAD[lor4] . GAD[lorz] . SUNT[sun] . GAD[lor3] . v1 . xprop[x]]]*atr[SUNTrace[TLPropagatorX[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia8=FourierXP[dia8,{x,q}]/.SUNN->CA;

dia8=QEvaluate[I ScaleMu^(4(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Code::Initialization::Plain:: *)
(*-------------------------------------------*)


(* ::Input::Initialization:: *)
ptype1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,lorz,sun,k1,k2,l,k},


dia1=(-I)*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[k - q] . v1 . prop[k1] . GAD[lorz] . SUNT[sun] . prop[k2] . hv4 . prop[-k + l] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[k1 - k2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - I*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v1 . prop[-k + q] . hv3 . prop[-k + l] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[k1 - k2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia1=IntegrateP[dia1,{k1,k2,k,l}]/.SUNN->CA;

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,lorz,sun,k1,k2,l,k},


dia2=I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[-k + l] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]*atr[SUNTrace[hv4 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v1 . prop[-k + q]]]]*FAD[k1 - k2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[k - q] . v1 . prop[k - l]]]*atr[SUNTrace[hv4 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[k1 - k2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia2=IntegrateP[dia2,{k1,k2,k,l}]/.SUNN->CA;

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,lorz,sun,k1,k2,l,k},


dia3=(-I)*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v1 . prop[-k + q] . hv4 . prop[-k + l] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[k1 - k2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - I*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[k - q] . v1 . prop[k1] . GAD[lorz] . SUNT[sun] . prop[k2] . hv3 . prop[-k + l] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[k1 - k2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia3=IntegrateP[dia3,{k1,k2,k,l}]/.SUNN->CA;

dia3=QEvaluate[I ScaleMu^(4(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,lorz,sun,k1,k2,l,k},


dia4=I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[-k2 - l] . GAD[lorz] . SUNT[sun] . prop[-k1 - l] . v2 . prop[k - l]]]*atr[SUNTrace[hv4 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v1 . prop[-k + q]]]]*FAD[k1 - k2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - I*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[k - q] . v1 . prop[k1] . GAD[lorz] . SUNT[sun] . prop[k2] . hv3 . prop[-k2 - l] . GAD[lorz] . SUNT[sun] . prop[-k1 - l] . v2 . prop[k - l]]]]*FAD[k1 - k2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia4=IntegrateP[dia4,{k1,k2,k,l}]/.SUNN->CA;

dia4=QEvaluate[I ScaleMu^(4(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,lorz,sun,k1,k2,l,k},


dia5=I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[-k2 - l] . GAD[lorz] . SUNT[sun] . prop[-k1 - l] . v2 . prop[k - l]]]*atr[SUNTrace[hv4 . prop[k - q] . v1 . prop[k1] . GAD[lorz] . SUNT[sun] . prop[k2]]]]*FAD[k1 - k2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[k - q] . v1 . prop[k1] . GAD[lorz] . SUNT[sun] . prop[k2]]]*atr[SUNTrace[hv4 . prop[-k2 - l] . GAD[lorz] . SUNT[sun] . prop[-k1 - l] . v2 . prop[k - l]]]]*FAD[k1 - k2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia5=IntegrateP[dia5,{k1,k2,k,l}]/.SUNN->CA;

dia5=QEvaluate[I ScaleMu^(4(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,lorz,sun,k1,k2,l,k},


dia6=I*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[k - q] . v1 . prop[k - l]]]*atr[SUNTrace[hv3 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[k1 - k2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + I*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[k - q] . v2 . prop[k - l]]]*atr[SUNTrace[hv3 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[k1 - k2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia6=IntegrateP[dia6,{k1,k2,k,l}]/.SUNN->CA;

dia6=QEvaluate[I ScaleMu^(4(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,lorz,sun,k1,k2,l,k},


dia7=I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[k - q] . v2 . prop[k - l]]]*atr[SUNTrace[hv4 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[k1 - k2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - I*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v1 . prop[-k + q] . hv3 . prop[-k2 - l] . GAD[lorz] . SUNT[sun] . prop[-k1 - l] . v2 . prop[k - l]]]]*FAD[k1 - k2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia7=IntegrateP[dia7,{k1,k2,k,l}]/.SUNN->CA;

dia7=QEvaluate[I ScaleMu^(4(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,lorz,sun,k1,k2,l,k},


dia8=(-I)*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v1 . prop[-k + q] . hv4 . prop[-k2 - l] . GAD[lorz] . SUNT[sun] . prop[-k1 - l] . v2 . prop[k - l]]]]*FAD[k1 - k2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[k - q] . v1 . prop[k1] . GAD[lorz] . SUNT[sun] . prop[k2] . hv4 . prop[-k2 - l] . GAD[lorz] . SUNT[sun] . prop[-k1 - l] . v2 . prop[k - l]]]]*FAD[k1 - k2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia8=IntegrateP[dia8,{k1,k2,k,l}]/.SUNN->CA;

dia8=QEvaluate[I ScaleMu^(4(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,lorz,sun,k1,k2,l,k},


dia9=(-I)*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[k - q] . v1 . prop[k - l] . hv4 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[k1 - k2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] + I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v1 . prop[-k + q]]]*atr[SUNTrace[hv4 . prop[-k + l] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[k1 - k2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia9=IntegrateP[dia9,{k1,k2,k,l}]/.SUNN->CA;

dia9=QEvaluate[I ScaleMu^(4(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,lorz,sun,k1,k2,l,k},


dia10=I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[-k + l] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]*atr[SUNTrace[hv4 . prop[k - q] . v1 . prop[k1] . GAD[lorz] . SUNT[sun] . prop[k2]]]]*FAD[k1 - k2]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[k - q] . v1 . prop[k1] . GAD[lorz] . SUNT[sun] . prop[k2]]]*atr[SUNTrace[hv4 . prop[-k + l] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[k1 - k2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia10=IntegrateP[dia10,{k1,k2,k,l}]/.SUNN->CA;

dia10=QEvaluate[I ScaleMu^(4(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,lorz,sun,k1,k2,l,k},


dia11=(-I)*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l] . hv4 . prop[k - q] . v2 . prop[k - l]]]]*FAD[k1 - k2]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - I*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[k - q] . v1 . prop[k - l] . hv3 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[k1 - k2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia11=IntegrateP[dia11,{k1,k2,k,l}]/.SUNN->CA;

dia11=QEvaluate[I ScaleMu^(4(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,lorz,sun,k1,k2,l,k},


dia12=(-I)*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k2 + l] . hv3 . prop[k - q] . v2 . prop[k - l]]]]*FAD[k1 - k2]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[-k2] . GAD[lorz] . SUNT[sun] . prop[-k1] . v1 . prop[-k + q]]]*atr[SUNTrace[hv4 . prop[-k2 - l] . GAD[lorz] . SUNT[sun] . prop[-k1 - l] . v2 . prop[k - l]]]]*FAD[k1 - k2]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia12=IntegrateP[dia12,{k1,k2,k,l}]/.SUNN->CA;

dia12=QEvaluate[I ScaleMu^(4(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*-------------------------------------------*)


(* ::Input::Initialization:: *)
iptype1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{k,l,q,dia1,sun,lor1,lor2,lor3,lor4,lorz,dot,str,tr,contract},


dia1=contract[TLPropagatorP[k, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, prop[-l - q], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, prop[-k - l]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + contract[TLPropagatorP[k, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, prop[k + l], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]*tr[str[dot[hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, prop[l + q]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia1=IntegrateP[dia1,{k,l}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{k,l,q,dia2,sun,lor1,lor2,lor3,lor4,lorz,dot,str,tr,contract},


dia2=contract[tr[str[dot[hv4, prop[k + l], v1, prop[l + q]]]]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - contract[tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2], hv3, prop[k + l], v2, prop[l + q]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia2=IntegrateP[dia2,{k,l}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{k,l,q,dia3,sun,lor1,lor2,lor3,lor4,lorz,dot,str,tr,contract},


dia3=contract[TLPropagatorP[k, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, prop[l + q]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, prop[-k - l]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + contract[TLPropagatorP[k, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, prop[k + l], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]*tr[str[dot[hv4, prop[-l - q], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia3=IntegrateP[dia3,{k,l}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=QEvaluate[I ScaleMu^(4(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{k,l,q,dia4,sun,lor1,lor2,lor3,lor4,lorz,dot,str,tr,contract},


dia4=-(contract[tr[str[dot[hv4, prop[k + l], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, prop[l + q]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) - contract[tr[str[dot[hv4, prop[k + l], v1, prop[l + q], hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia4=IntegrateP[dia4,{k,l}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=QEvaluate[I ScaleMu^(4(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{k,l,q,dia5,sun,lor1,lor2,lor3,lor4,lorz,dot,str,tr,contract},


dia5=-(contract[tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, prop[-k - l], hv4, prop[-l - q], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - contract[tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, prop[-k - l], hv3, prop[-l - q], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia5=IntegrateP[dia5,{k,l}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=QEvaluate[I ScaleMu^(4(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{k,l,q,dia6,sun,lor1,lor2,lor3,lor4,lorz,dot,str,tr,contract},


dia6=-(contract[tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2], hv4, prop[k + l], v2, prop[l + q]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) + contract[tr[str[dot[hv3, prop[k + l], v1, prop[l + q]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia6=IntegrateP[dia6,{k,l}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia6=QEvaluate[I ScaleMu^(4(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{k,l,q,dia7,sun,lor1,lor2,lor3,lor4,lorz,dot,str,tr,contract},


dia7=contract[TLPropagatorP[k, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, prop[l + q]]]]*tr[str[dot[hv4, prop[k + l], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + contract[TLPropagatorP[k, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, prop[-l - q], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]*tr[str[dot[hv4, prop[k + l], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia7=IntegrateP[dia7,{k,l}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia7=QEvaluate[I ScaleMu^(4(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{k,l,q,dia8,sun,lor1,lor2,lor3,lor4,lorz,dot,str,tr,contract},


dia8=contract[tr[str[dot[hv3, prop[k + l], v2, prop[l + q]]]]*tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - contract[tr[str[dot[hv4, prop[k + l], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv3, prop[-l - q], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia8=IntegrateP[dia8,{k,l}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia8=QEvaluate[I ScaleMu^(4(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{k,l,q,dia9,sun,lor1,lor2,lor3,lor4,lorz,dot,str,tr,contract},


dia9=contract[TLPropagatorP[k, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, prop[-k - l]]]]*tr[str[dot[hv4, prop[-l - q], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + contract[tr[str[dot[hv4, prop[k + l], v2, prop[l + q]]]]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia9=IntegrateP[dia9,{k,l}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia9=QEvaluate[I ScaleMu^(4(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{k,l,q,dia10,sun,lor1,lor2,lor3,lor4,lorz,dot,str,tr,contract},


dia10=-(contract[tr[str[dot[hv3, prop[k + l], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, prop[l + q]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - contract[tr[str[dot[hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, prop[-k - l], hv3, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, prop[l + q]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia10=IntegrateP[dia10,{k,l}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia10=QEvaluate[I ScaleMu^(4(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{k,l,q,dia11,sun,lor1,lor2,lor3,lor4,lorz,dot,str,tr,contract},


dia11=-(contract[tr[str[dot[hv3, prop[k + l], v1, prop[l + q], hv4, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) + contract[TLPropagatorP[k, lor1, lor2, lor3, lor4]*tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, prop[-k - l]]]]*tr[str[dot[hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, prop[l + q]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia11=IntegrateP[dia11,{k,l}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia11=QEvaluate[I ScaleMu^(4(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{k,l,q,dia12,sun,lor1,lor2,lor3,lor4,lorz,dot,str,tr,contract},


dia12=-(contract[tr[str[dot[hv3, GAD[lor4], GAD[lorz], SUNT[sun], GAD[lor3], v1, prop[-k - l], hv4, GAD[lor2], GAD[lorz], SUNT[sun], GAD[lor1], v2, prop[l + q]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - contract[tr[str[dot[hv3, prop[k + l], v1, GAD[lor3], GAD[lorz], SUNT[sun], GAD[lor4], hv4, prop[-l - q], v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]*TLPropagatorP[k, lor1, lor2, lor3, lor4]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia12=IntegrateP[dia12,{k,l}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia12=QEvaluate[I ScaleMu^(4(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*-------------------------------------------*)


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
