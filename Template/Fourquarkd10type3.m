(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
Fourquarkd10type3::usage="Fourquarkd10type3[q_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}}] gives the third type \[LeftAngleBracket]\!\(\*OverscriptBox[\(\[Psi]\), \(_\)]\)G\[Psi]\!\(\*SuperscriptBox[\(\[RightAngleBracket]\), \(2\)]\)\[InvisibleComma] contribution of OPE result about <J1 J2>, with J1=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a1\)]\)\!\(\*SubscriptBox[\(v1\[CapitalPsi]\), \(b1\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c1\)]\)\!\(\*SubscriptBox[\(v2\[CapitalPsi]\), \(d1\)]\) and J2=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(b2\)]\)\!\(\*SubscriptBox[\(v3\[CapitalPsi]\), \(a2\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(d2\)]\)\!\(\*SubscriptBox[\(v4\[CapitalPsi]\), \(c2\)]\), here the background gluon comes from quark propatator. "

Fourquarkd10type3::inderr="Dummy indices conflict!"

Fourquarkd10type3::curerr="Unknow current structure!"



(* ::Code::Initialization::Plain:: *)
Begin["`Private`Fourquarkd10type3`"]


(* ::Code::Initialization::Plain:: *)
Options[Fourquarkd10type3] = {
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
Fourquarkd10type3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmp,tmp1,tmp2,hv3,hv4,holdf=OptionValue[HoldFlavor],atr,diagrams,fdir,pall,files,dot,setd},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[Fourquarkd10type3::inderr];
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


diagrams={xtype1,xtype2,xtype3,xtype4,xtype5,xtype6};

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
			files=("Fourquarkd10type3_"<>#)&/@files;
			
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


 dia1=-1/256*(Condensate[{a1, "G", a2}]*Condensate[{b2, "G", d1}]*contract[tr[str[dot[hv3, v1, xprop[x], hv4, xprop[-x], v2, 1]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x]^2)/(CA^2*D^2) - (Condensate[{b2, "G", b1}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv4, xprop[-x], v1, hv3, 1, v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x]^2)/(256*CA^2*D^2) + (Condensate[{c1, "G", c2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, v2, 1]]]*tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x]^2)/(256*CA^2*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{b2, "G", b1}]*contract[tr[str[dot[hv3, v1, 1]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x]^2)/(256*CA^2*D^2);


 If[setd===4,dia1=FCI[dia1]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia1=FourierXP[dia1,{x,q}];


 dia1=QEvaluate[I ScaleMu^(1(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia2,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia2=(Condensate[{a1, "G", c2}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv3, v2, xprop[x]]]]*tr[str[dot[hv4, v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x]^2)/(256*CA^2*D^2) - (Condensate[{c1, "G", c2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprop[-x], v1, hv4, 1, v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x]^2)/(256*CA^2*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{b2, "G", b1}]*contract[tr[str[dot[hv4, v1, 1, hv3, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x]^2)/(256*CA^2*D^2) + (Condensate[{b2, "G", b1}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, xprop[-x], v1]]]*tr[str[dot[hv4, v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x]^2)/(256*CA^2*D^2);


 If[setd===4,dia2=FCI[dia2]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia2=FourierXP[dia2,{x,q}];


 dia2=QEvaluate[I ScaleMu^(1(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia3,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia3=(Condensate[{b2, "G", d1}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv3, v2, 1]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x]^2)/(256*CA^2*D^2) - (Condensate[{a1, "G", a2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, v1, 1, hv4, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x]^2)/(256*CA^2*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, v1, xprop[x], hv3, xprop[-x], v2, 1]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x]^2)/(256*CA^2*D^2) - (Condensate[{c1, "G", a2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, v2, 1]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x]^2)/(256*CA^2*D^2);


 If[setd===4,dia3=FCI[dia3]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia3=FourierXP[dia3,{x,q}];


 dia3=QEvaluate[I ScaleMu^(1(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia4,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia4=(Condensate[{c1, "G", a2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x]^2)/(256*CA^2*D^2) - (Condensate[{a1, "G", a2}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, v1, xprop[x], hv4, 1, v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x]^2)/(256*CA^2*D^2) - (Condensate[{b2, "G", b1}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv4, xprop[-x], v1, hv3, xprop[-x], v2, 1]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x]^2)/(256*CA^2*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv3, v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x]^2)/(256*CA^2*D^2);


 If[setd===4,dia4=FCI[dia4]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia4=FourierXP[dia4,{x,q}];


 dia4=QEvaluate[I ScaleMu^(1(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia5,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia5=(Condensate[{a1, "G", c2}]*Condensate[{b2, "G", d1}]*contract[tr[str[dot[hv3, xprop[-x], v2]]]*tr[str[dot[hv4, v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x]^2)/(256*CA^2*D^2) + (Condensate[{a1, "G", c2}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv4, v1, 1]]]*tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x]^2)/(256*CA^2*D^2) - (Condensate[{a1, "G", c2}]*Condensate[{c1, "G", a2}]*contract[tr[str[dot[hv4, v1, xprop[x], hv3, 1, v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]*SPD[x]^2)/(256*CA^2*D^2) + (Condensate[{b2, "G", b1}]*Condensate[{d2, "G", d1}]*contract[tr[str[dot[hv3, xprop[-x], v1]]]*tr[str[dot[hv4, xprop[-x], v2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x]^2)/(256*CA^2*D^2);


 If[setd===4,dia5=FCI[dia5]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia5=FourierXP[dia5,{x,q}];


 dia5=QEvaluate[I ScaleMu^(1(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Input::Initialization::Plain:: *)
xtype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia6,lor1,lor2,lora,lorb,sun,sun1,sun2,tr,str,contract,dot},


 dia6=(Condensate[{b2, "G", d1}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprop[-x], v2]]]*tr[str[dot[hv4, xprop[-x], v1]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]*SPD[x]^2)/(256*CA^2*D^2) - (Condensate[{b2, "G", d1}]*Condensate[{d2, "G", b1}]*contract[tr[str[dot[hv3, xprop[-x], v1, hv4, xprop[-x], v2, 1]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x]^2)/(256*CA^2*D^2) - (Condensate[{b2, "G", d1}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, v2, 1]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]*SPD[x]^2)/(256*CA^2*D^2) + (Condensate[{a1, "G", a2}]*Condensate[{c1, "G", c2}]*contract[tr[str[dot[hv3, v1, xprop[x]]]]*tr[str[dot[hv4, v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]*SPD[x]^2)/(256*CA^2*D^2);


 If[setd===4,dia6=FCI[dia6]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


 dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia6=FourierXP[dia6,{x,q}];


 dia6=QEvaluate[I ScaleMu^(1(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)


]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
