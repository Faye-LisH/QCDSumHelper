(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
Fourquarkd10type2::usage="Fourquarkd10type2[q_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}}] gives the \[InvisibleComma]second type \[LeftAngleBracket]\!\(\*OverscriptBox[\(\[Psi]\), \(_\)]\)G\[Psi]\!\(\*SuperscriptBox[\(\[RightAngleBracket]\), \(2\)]\)\[InvisibleComma] contribution of OPE result about <J1 J2>, with J1=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a1\)]\)\!\(\*SubscriptBox[\(v1\[CapitalPsi]\), \(b1\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c1\)]\)\!\(\*SubscriptBox[\(v2\[CapitalPsi]\), \(d1\)]\) and J2=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(b2\)]\)\!\(\*SubscriptBox[\(v3\[CapitalPsi]\), \(a2\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(d2\)]\)\!\(\*SubscriptBox[\(v4\[CapitalPsi]\), \(c2\)]\), here the background gluon comes from quark propatator. "

Fourquarkd10type2::inderr="Dummy indices conflict!"

Fourquarkd10type2::curerr="Unknow current structure!"



(* ::Code::Initialization::Plain:: *)
Begin["`Private`Fourquarkd10type2`"]


(* ::Code::Initialization::Plain:: *)
Options[Fourquarkd10type2] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True,
	ToD4->"Auto"
}


(* ::Code::Initialization::Plain:: *)
(* VS-first *)


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(* propagators *)

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
Fourquarkd10type2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmp,tmp1,tmp2,hv3,hv4,holdf=OptionValue[HoldFlavor],atr,diagrams,fdir,pall,files,dot,setd},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[Fourquarkd10type2::inderr];
	Abort[]
];


If[OptionValue[AutoNDR]===True,
	atr=TR5
	,
	atr=TR
];


If[OptionValue[ToD4]==="Auto"||OptionValue[ToD4]===4,
	setd=4
,
	setd=D
];

(* B A^+ B *)
hv3=v3//ComplexConjugate;
hv4=v4//ComplexConjugate//FCI;
(* ComplexConjugate[sigma] will expand the DiracSigma to gamma matrices, recover it to DiracSigma *)
hv4=(hv4/.Dot->dot)/.(f1_ dot[aa___,DiracGamma[LorentzIndex[lor1_,dim___],dim___],DiracGamma[LorentzIndex[lor2_,dim___],dim___],bb___]+f2_ dot[aa___,DiracGamma[LorentzIndex[lor2_,dim___],dim___],DiracGamma[LorentzIndex[lor1_,dim___],dim___],bb___])/;(f1+f2===0):>-2I f1 dot[aa,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],bb];
hv4=hv4/. {null->1,dot->Dot};


diagrams={xtype1,xtype2,xtype3,xtype4,xtype5,xtype6,xtype7,xtype8,xtype9,xtype10,xtype11,xtype12,xtype13,xtype14,xtype15,xtype16,xtype17,xtype18,xtype19,xtype20,xtype21,xtype22,xtype23,xtype24};

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
	
	tmp=Plus@@WaitAll[ParallelSubmit[{hv3,hv4,holdf,atr,setd},#[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,setd]]&/@diagrams];
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	
,

	If[pall==="External",
	
		DistributeDefinitions[qq,v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
		
		If[fdir==="None",
		
			ParallelSubmit[{hv3,hv4,holdf,atr,setd},#[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,setd]]&/@diagrams
		,
		
			files=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files=("Fourquarkd10type2_"<>#)&/@files;
			
			ImExport[fdir,
					files,
					{{hv3,hv4,holdf,atr,setd},
					{qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,setd},
					diagrams}
					]
		]
		(* !!! do not touch tmp before WaitAll[tmp], any function may modifiy the things in EvaluationObject[...] should not be used here. !!! *)
	,
	
	
		tmp=Plus@@Map[#[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,setd]&,diagrams];
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	]
]

]



(* ::Code::Initialization::Plain:: *)
(*--- The following are generated by algorithm ---*)


(* ::Input::Initialization::Plain:: *)
xtype1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia1,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia1=(Condensate[{c1, "G", a2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v1]]]*tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", a2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], hv4, xprog[-x, lor1, lor2, sun], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", d1}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], hv4, xprop[-x], v2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv3, v1, xprop[x]]]]*tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia1=FCI[dia1]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia1=FourierXP[dia1,{x,q}];


 dia1=QEvaluate[I ScaleMu^(1(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia2,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia2=(Condensate[{a1, "G", c2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv4, v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]*tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", a2}]*Condensate[{b2, "G", d1}]*contract[tr[str[dot[hv3, v1, xprog[x, lor1, lor2, sun], hv4, xprop[-x], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{b2, "G", b1}]*contract[tr[str[dot[hv4, v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], hv3, xprog[-x, lor1, lor2, sun], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv4, v2, xprop[x]]]]*tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia2=FCI[dia2]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia2=FourierXP[dia2,{x,q}];


 dia2=QEvaluate[I ScaleMu^(1(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia3,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia3=(Condensate[{a1, "G", c2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv4, v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]*tr[str[dot[hv3, xprop[-x], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", d1}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2]]]*tr[str[dot[hv4, xprop[-x], v1, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{c1, "G", c2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprop[-x], v1, hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", b1}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v1]]]*tr[str[dot[hv4, xprop[-x], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia3=FCI[dia3]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia3=FourierXP[dia3,{x,q}];


 dia3=QEvaluate[I ScaleMu^(1(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia4,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia4=(Condensate[{a1, "G", c2}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv4, v1, xprog[x, lor1, lor2, sun]]]]*tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{c1, "G", a2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, v2, xprog[x, lor1, lor2, sun]]]]*tr[str[dot[hv4, xprop[-x], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", a2}]*Condensate[{b2, "G", d1}]*contract[tr[str[dot[hv3, v1, xprop[x], hv4, xprog[-x, lor1, lor2, sun], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", b1}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v1]]]*tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia4=FCI[dia4]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia4=FourierXP[dia4,{x,q}];


 dia4=QEvaluate[I ScaleMu^(1(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia5,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia5=(Condensate[{a1, "G", c2}]*Condensate[{b2, "G", d1}]*contract[tr[str[dot[hv4, v1, xprog[x, lor1, lor2, sun]]]]*tr[str[dot[hv3, xprop[-x], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", d1}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv3, v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]*tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", d1}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog[x, lor1, lor2, sun], hv4, v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", b1}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv4, xprop[-x], v1, hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia5=FCI[dia5]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia5=FourierXP[dia5,{x,q}];


 dia5=QEvaluate[I ScaleMu^(1(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia6,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia6=(Condensate[{a1, "G", c2}]*Condensate[{b2, "G", d1}]*contract[tr[str[dot[hv3, xprop[-x], v2]]]*tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", d1}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v2]]]*tr[str[dot[hv4, xprop[-x], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{c1, "G", c2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], hv4, v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprop[x], hv3, xprog[-x, lor1, lor2, sun], v2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia6=FCI[dia6]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia6=FourierXP[dia6,{x,q}];


 dia6=QEvaluate[I ScaleMu^(1(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia7,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia7=(Condensate[{b2, "G", b1}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprop[-x], v2]]]*tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", b1}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v2]]]*tr[str[dot[hv3, xprop[-x], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{c1, "G", c2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]*tr[str[dot[hv3, xprop[-x], v1, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", b1}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, xprop[-x], v1]]]*tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia7=FCI[dia7]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia7=FourierXP[dia7,{x,q}];


 dia7=QEvaluate[I ScaleMu^(1(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia8,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia8=-1/64*(Condensate[{a1, "G", a2}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprog[x, lor1, lor2, sun], hv4, v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", d1}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v1, xprop[x], hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{c1, "G", a2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v1, xprop[x], hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v2]]]*tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia8=FCI[dia8]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia8=FourierXP[dia8,{x,q}];


 dia8=QEvaluate[I ScaleMu^(1(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia9,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia9=(Condensate[{a1, "G", c2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1]]]*tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv4, v1, xprop[x], hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{c1, "G", c2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2]]]*tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v1, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{b2, "G", b1}]*contract[tr[str[dot[hv3, v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]*tr[str[dot[hv4, xprop[-x], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia9=FCI[dia9]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia9=FourierXP[dia9,{x,q}];


 dia9=QEvaluate[I ScaleMu^(1(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia10,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia10=(Condensate[{a1, "G", c2}]*Condensate[{b2, "G", d1}]*contract[tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v2]]]*tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", d1}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprop[-x], v2]]]*tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{b2, "G", b1}]*contract[tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, hv3, xprog[-x, lor1, lor2, sun], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{b2, "G", b1}]*contract[tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1]]]*tr[str[dot[hv4, xprop[-x], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia10=FCI[dia10]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia10=FourierXP[dia10,{x,q}];


 dia10=QEvaluate[I ScaleMu^(1(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia11,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia11=(Condensate[{b2, "G", d1}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv4, xprop[-x], v1]]]*tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", d1}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprop[-x], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], hv4, xprog[-x, lor1, lor2, sun], v2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", d1}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog[x, lor1, lor2, sun], hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprop[-x], v2]]]*tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia11=FCI[dia11]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia11=FourierXP[dia11,{x,q}];


 dia11=QEvaluate[I ScaleMu^(1(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia12,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia12=-1/64*(Condensate[{c1, "G", c2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprop[-x], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], hv4, v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprog[x, lor1, lor2, sun], hv3, v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{c1, "G", a2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprop[-x], v1, xprog[x, lor1, lor2, sun], hv3, v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, v1, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia12=FCI[dia12]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia12=FourierXP[dia12,{x,q}];


 dia12=QEvaluate[I ScaleMu^(1(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia13,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia13=(Condensate[{c1, "G", a2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv4, xprop[-x], v1]]]*tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", a2}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprop[x], hv4, v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", b1}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], hv3, xprop[-x], v2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{c1, "G", a2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprop[-x], v1, xprog[x, lor1, lor2, sun], hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia13=FCI[dia13]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia13=FourierXP[dia13,{x,q}];


 dia13=QEvaluate[I ScaleMu^(1(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia14,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia14=(Condensate[{a1, "G", c2}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv4, v1, xprop[x]]]]*tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", a2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, hv4, xprop[-x], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{c1, "G", c2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v1, hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprop[x], hv3, v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia14=FCI[dia14]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia14=FourierXP[dia14,{x,q}];


 dia14=QEvaluate[I ScaleMu^(1(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia15,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia15=(Condensate[{a1, "G", c2}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv3, v2, xprog[x, lor1, lor2, sun]]]]*tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", b1}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv4, xprop[-x], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], hv3, v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, v1, xprog[x, lor1, lor2, sun]]]]*tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{b2, "G", b1}]*contract[tr[str[dot[hv3, v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]*tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia15=FCI[dia15]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia15=FourierXP[dia15,{x,q}];


 dia15=QEvaluate[I ScaleMu^(1(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia16,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia16=-1/64*(Condensate[{a1, "G", a2}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, v1, xprog[x, lor1, lor2, sun], hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", d1}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v1, hv4, xprop[-x], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", b1}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v1, hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", b1}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv3, xprop[-x], v1]]]*tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia16=FCI[dia16]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia16=FourierXP[dia16,{x,q}];


 dia16=QEvaluate[I ScaleMu^(1(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype17[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia17,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia17=(Condensate[{c1, "G", a2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, v2, xprop[x]]]]*tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", b1}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprop[-x], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], hv3, xprog[-x, lor1, lor2, sun], v2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv4, v2, xprog[x, lor1, lor2, sun]]]]*tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{b2, "G", b1}]*contract[tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1]]]*tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia17=FCI[dia17]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia17=dia17/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia17=FourierXP[dia17,{x,q}];


 dia17=QEvaluate[I ScaleMu^(1(4-D)) dia17,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype18[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia18,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia18=(Condensate[{b2, "G", d1}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v1]]]*tr[str[dot[hv3, xprop[-x], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", a2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], hv4, xprop[-x], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{b2, "G", b1}]*contract[tr[str[dot[hv4, v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], hv3, xprop[-x], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", b1}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], hv3, v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia18=FCI[dia18]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia18=dia18/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia18=FourierXP[dia18,{x,q}];


 dia18=QEvaluate[I ScaleMu^(1(4-D)) dia18,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype19[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia19,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia19=-1/64*(Condensate[{b2, "G", d1}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprop[-x], v1, hv4, xprog[-x, lor1, lor2, sun], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv4, v1, xprog[x, lor1, lor2, sun], hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, v1, xprop[x], hv3, xprog[-x, lor1, lor2, sun], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{b2, "G", b1}]*contract[tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, hv3, xprop[-x], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia19=FCI[dia19]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia19=dia19/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia19=FourierXP[dia19,{x,q}];


 dia19=QEvaluate[I ScaleMu^(1(4-D)) dia19,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype20[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia20,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia20=(Condensate[{a1, "G", c2}]*Condensate[{b2, "G", d1}]*contract[tr[str[dot[hv4, v1, xprop[x]]]]*tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", d1}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2]]]*tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", a2}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, v1, xprop[x], hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", a2}]*Condensate[{b2, "G", d1}]*contract[tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprop[x], hv4, xprog[-x, lor1, lor2, sun], v2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia20=FCI[dia20]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia20=dia20/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia20=FourierXP[dia20,{x,q}];


 dia20=QEvaluate[I ScaleMu^(1(4-D)) dia20,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype21[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia21,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia21=(Condensate[{a1, "G", c2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1]]]*tr[str[dot[hv3, xprop[-x], v2, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", a2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, hv4, xprog[-x, lor1, lor2, sun], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprog[x, lor1, lor2, sun], hv3, xprop[-x], v2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", b1}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv4, v2, xprop[x]]]]*tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia21=FCI[dia21]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia21=dia21/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia21=FourierXP[dia21,{x,q}];


 dia21=QEvaluate[I ScaleMu^(1(4-D)) dia21,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype22[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia22,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia22=(Condensate[{a1, "G", c2}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv3, v2, xprop[x]]]]*tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", d1}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv3, v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]*tr[str[dot[hv4, xprop[-x], v1, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", a2}]*Condensate[{b2, "G", d1}]*contract[tr[str[dot[hv3, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v1, xprog[x, lor1, lor2, sun], hv4, xprop[-x], v2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{c1, "G", c2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, DiracSigma[GAD[lor1, lor2]] . SUNT[sun], v2]]]*tr[str[dot[hv3, xprop[-x], v1, xprog[x, lor1, lor2, sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia22=FCI[dia22]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia22=dia22/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia22=FourierXP[dia22,{x,q}];


 dia22=QEvaluate[I ScaleMu^(1(4-D)) dia22,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype23[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia23,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia23=-1/64*(Condensate[{b2, "G", d1}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v1, xprop[x], hv4, v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x])/(CA^2*CF*(-1 + D)*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, v1, xprog[x, lor1, lor2, sun], hv3, xprop[-x], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) - (Condensate[{b2, "G", b1}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprop[-x], v1, hv3, xprog[-x, lor1, lor2, sun], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{b2, "G", b1}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv4, v2, xprog[x, lor1, lor2, sun]]]]*tr[str[dot[hv3, xprop[-x], v1, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia23=FCI[dia23]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia23=dia23/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia23=FourierXP[dia23,{x,q}];


 dia23=QEvaluate[I ScaleMu^(1(4-D)) dia23,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype24[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia24,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia24=-1/64*(Condensate[{b2, "G", b1}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v1, hv3, xprop[-x], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(CA^2*CF*(-1 + D)*D^2) - (Condensate[{c1, "G", a2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprog[-x, lor1, lor2, sun], v1, xprop[x], hv3, v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{c1, "G", c2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]*tr[str[dot[hv3, xprog[-x, lor1, lor2, sun], v1, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv3, v1, xprog[x, lor1, lor2, sun]]]]*tr[str[dot[hv4, xprop[-x], v2, DiracSigma[GAD[lor1, lor2]] . SUNT[sun]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x])/(64*CA^2*CF*(-1 + D)*D^2);


 If[setd===4,dia24=FCI[dia24]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia24=dia24/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia24=FourierXP[dia24,{x,q}];


 dia24=QEvaluate[I ScaleMu^(1(4-D)) dia24,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
