(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
Tetraquarkd10type0::usage="Tetraquarkd10type0[q_,j1_,j2_] gives \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)q\!\(\*SuperscriptBox[\(\[RightAngleBracket]\), \(2\)]\)\[LeftAngleBracket]GG\[RightAngleBracket] contribution of \[LeftAngleBracket]j1 j2\[RightAngleBracket], with j1, j2 are tetraquark currents. "

Tetraquarkd10type0::inderr="Dummy indices conflict!"

Tetraquarkd10type0::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`Tetraquarkd10type0`"]


(* ::Code::Initialization::Plain:: *)
Options[Tetraquarkd10type0] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True,
	ToD4->"Auto"
}


(* ::Code::Initialization::Plain:: *)
(* VS-first *)


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

prop[q_]=I GSD[q]FAD[q];


(* ::Code::Initialization::Plain:: *)
xprop[x_]=(I^(1 - D)*DiracGamma[Momentum[x, D], D]*qGamma[D/2])/(2*Pi^(D/2)*Pair[Momentum[x, D], Momentum[x, D]]^(D/2));(*FourierPX[prop[q],{q,x}]*)


(* ::Code::Initialization::Plain:: *)
(*--- propagator with one background gluon ---*)

xprog[x_,lora_,lorb_]=(GAD[lora] . GAD[lorb] . GSD[x] qfact1[1/32 I^-D \[Pi]^(-D/2) qGamma[-1+D/2]] SPD[x,x]^(1-D/2) 
							+GSD[x] . GAD[lorb] . GAD[lora] qfact1[qfact2[-(1/32) I^-D \[Pi]^(-D/2)] qGamma[-1+D/2]] SPD[x,x]^(1-D/2) 
							+2 FVD[x,lora] FVD[x,lorb] GSD[x] qfact1[-(1/16) I^-D \[Pi]^(-D/2) qGamma[D/2]] SPD[x,x]^(-D/2) );


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
Tetraquarkd10type0[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 Tetraquarkd10type0[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])


(* ::Code::Initialization::Plain:: *)
Tetraquarkd10type0[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,holdf=OptionValue[HoldFlavor],atr,fdir,files,pall,diagrams,setd},



(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[Tetraquarkd10type0::curerr];
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

(*--- check dummy indices in input ---*)
tmp2=Length[#]&/@Gather[({c11,c12,c13,c14,c21,c22,c23,c24}//FCI)/.SUNFIndex[aa_]:>aa];

If[DummyCheck[v11,v12,v21,v22]==0||!FreeQ[tmp2,nn_Integer/;nn>2],
	Message[Tetraquarkd10type0::inderr];
	Abort[]
];



(*-------------------------------*)
(* B A^+ B *)
hv21=v21//ComplexConjugate;
hv22=v22//ComplexConjugate;


(*---------------------------------------------------*)
If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];


diagrams={xtype1,xtype2,xtype3,xtype4,xtype5,xtype6,xtype7,xtype8,xtype9,xtype10};


(*---------------------------------------------------*)
If[pall===True,

	
	tmp=Plus@@WaitAll[ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd]]&/@diagrams];
	
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

,
	If[pall==="External",

	DistributeDefinitions[qq,current1,current2];
	
		If[fdir==="None",
			(* evaluation, no import and export *)
			 ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd]]&/@diagrams
		
		,
		(* evaluation, import and export *)
			files=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files=("Tetraquarkd10type0_"<>#)&/@files;
			
			
			ImExport[fdir,
						files,
						{{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd},
						{qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd},
						diagrams}
						]
				
		](* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
	
	,
	
		
		tmp=Plus@@(#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd]&/@diagrams);
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]


]

(*------------------------------------------------------------------*)



(* ::Code::Initialization::Plain:: *)
(*--- The following are generated by algorithm ---*)


(* ::Input::Initialization::Plain:: *)
xtype1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia1,cc,tr,dot,lor1,lor2,lora,lorb,sun,sunSimplify,contract},


dia1=(Pi*Condensate["gg"]*Condensate[{f13, f24}]*Condensate[{f14, f23}]*contract[tr[dot[cc[dot[1, v12]], hv22, 1]]*tr[dot[cc[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]], v11, xprog[x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f12}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[1, hv21]], hv22, xprog[-x, lor1, lor2], cc[dot[v11, 1, v12]], (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c12, c13]*ColorDelta[c24, c22]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f21, f12}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[1, hv21]], v11, 1]]*tr[dot[cc[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[-x, lora, lorb], v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f23}]*Condensate[{f21, f11}]*contract[tr[dot[1, cc[dot[v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]]]]*tr[dot[cc[dot[hv22, 1, v12]], xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14]])/(8*CA^3*CF*(-1 + D)*D);

If[setd===4,dia1=FCI[dia1]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia1=FourierXP[dia1,{x,q}];

dia1=QEvaluate[I ScaleMu^(4-D) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia2,cc,tr,dot,lor1,lor2,lora,lorb,sun,sunSimplify,contract},


dia2=(Pi*Condensate["gg"]*Condensate[{f14, f24}]*Condensate[{f21, f11}]*contract[tr[dot[1, cc[dot[v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]]]]*tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], 1]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f14, f23}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[1, hv21]], v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb]]]*tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, 1]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f14, f11}]*Condensate[{f22, f23}]*contract[tr[dot[cc[dot[1]], v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], cc[dot[hv22, 1, hv21]], xprog[-x, lor1, lor2], v12]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c23, c22]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f23}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[hv22, 1, v12]], xprog[-x, lor1, lor2]]]*tr[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], cc[dot[v11, 1, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14]])/(8*CA^3*CF*(-1 + D)*D);

If[setd===4,dia2=FCI[dia2]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia2=FourierXP[dia2,{x,q}];

dia2=QEvaluate[I ScaleMu^(4-D) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia3,cc,tr,dot,lor1,lor2,lora,lorb,sun,sunSimplify,contract},


dia3=(Pi*Condensate["gg"]*Condensate[{f14, f24}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], 1]]*tr[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], cc[dot[v11, 1, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f14, f12}]*Condensate[{f21, f24}]*contract[tr[dot[cc[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]], v11, 1, cc[dot[hv22, xprog[-x, lor1, lor2], v12]], 1]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c12, c14]*ColorDelta[c24, c21]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f14, f24}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], 1]]*tr[dot[cc[dot[1, hv21]], v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f24}]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[1, v12]], hv22, xprog[-x, lor1, lor2]]]*tr[dot[cc[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]], v11, 1]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14]])/(8*CA^3*CF*(-1 + D)*D);

If[setd===4,dia3=FCI[dia3]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia3=FourierXP[dia3,{x,q}];

dia3=QEvaluate[I ScaleMu^(4-D) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia4,cc,tr,dot,lor1,lor2,lora,lorb,sun,sunSimplify,contract},


dia4=(Pi*Condensate["gg"]*Condensate[{f13, f11}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[1, hv21]], hv22, xprog[-x, lor1, lor2], cc[dot[1, v12]], v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c24, c22]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f14, f11}]*Condensate[{f21, f23}]*contract[tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, 1, cc[dot[v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]], 1]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c23, c21]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f11}]*Condensate[{f22, f23}]*contract[tr[dot[cc[dot[1, v12]], v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], cc[dot[hv22, 1, hv21]], xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c23, c22]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f21, f12}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv22, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[-x, lora, lorb], v12]], xprog[-x, lor1, lor2]]]*tr[dot[cc[dot[1, hv21]], v11, 1]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]])/(8*CA^3*CF*(-1 + D)*D);

If[setd===4,dia4=FCI[dia4]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia4=FourierXP[dia4,{x,q}];

dia4=QEvaluate[I ScaleMu^(4-D) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia5,cc,tr,dot,lor1,lor2,lora,lorb,sun,sunSimplify,contract},


dia5=(Pi*Condensate["gg"]*Condensate[{f14, f11}]*Condensate[{f22, f24}]*contract[tr[dot[1, cc[dot[hv22, xprog[-x, lor1, lor2], v12]], 1, hv21, cc[dot[v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb]]]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c24, c22]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f12}]*Condensate[{f21, f24}]*contract[tr[dot[cc[dot[1]], hv22, xprog[-x, lor1, lor2], cc[dot[v11, 1, v12]], (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c12, c13]*ColorDelta[c24, c21]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f12}]*Condensate[{f22, f23}]*contract[tr[dot[cc[dot[hv22, 1, hv21]], xprog[-x, lor1, lor2], cc[dot[v11, 1, v12]], (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c12, c13]*ColorDelta[c23, c22]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f21, f11}]*Condensate[{f22, f12}]*contract[tr[dot[1, cc[dot[v11, 1, hv21]]]]*tr[dot[cc[dot[hv22, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[-x, lora, lorb], v12]], xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]])/(8*CA^3*CF*(-1 + D)*D);

If[setd===4,dia5=FCI[dia5]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia5=FourierXP[dia5,{x,q}];

dia5=QEvaluate[I ScaleMu^(4-D) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia6,cc,tr,dot,lor1,lor2,lora,lorb,sun,sunSimplify,contract},


dia6=(Pi*Condensate["gg"]*Condensate[{f14, f12}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb]]], v11, 1, cc[dot[hv22, xprog[-x, lor1, lor2], v12]], 1, hv21]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c12, c14]*ColorDelta[c24, c22]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f24}]*Condensate[{f22, f12}]*contract[tr[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], cc[dot[v11, 1, hv21]]]]*tr[dot[cc[dot[1, v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f11}]*Condensate[{f21, f24}]*contract[tr[dot[1, cc[dot[v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]], 1, v12, cc[dot[hv22, xprog[-x, lor1, lor2]]]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c24, c21]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f14, f23}]*Condensate[{f22, f12}]*contract[tr[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], cc[dot[v11, 1, hv21]]]]*tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, 1]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13]])/(8*CA^3*CF*(-1 + D)*D);

If[setd===4,dia6=FCI[dia6]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia6=FourierXP[dia6,{x,q}];

dia6=QEvaluate[I ScaleMu^(4-D) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia7,cc,tr,dot,lor1,lor2,lora,lorb,sun,sunSimplify,contract},


dia7=(Pi*Condensate["gg"]*Condensate[{f13, f23}]*Condensate[{f14, f24}]*contract[tr[dot[cc[dot[hv22, 1, v12]], 1]]*tr[dot[cc[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]], v11, xprog[x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f24}]*Condensate[{f21, f11}]*contract[tr[dot[1, cc[dot[v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]]]]*tr[dot[cc[dot[1, v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f14, f23}]*Condensate[{f21, f11}]*contract[tr[dot[1, cc[dot[v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]]]]*tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, 1]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f12}]*Condensate[{f21, f23}]*contract[tr[dot[cc[dot[hv22, 1]], xprog[-x, lor1, lor2], cc[dot[v11, 1, v12]], (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c12, c13]*ColorDelta[c23, c21]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14])])/(8*CA^3*CF*(-1 + D)*D);

If[setd===4,dia7=FCI[dia7]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia7=FourierXP[dia7,{x,q}];

dia7=QEvaluate[I ScaleMu^(4-D) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia8,cc,tr,dot,lor1,lor2,lora,lorb,sun,sunSimplify,contract},


dia8=(Pi*Condensate["gg"]*Condensate[{f13, f23}]*Condensate[{f14, f24}]*contract[tr[dot[cc[dot[hv22, 1, v12]], 1]]*tr[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f24}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[1, hv21]], v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb]]]*tr[dot[cc[dot[1, v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f14, f23}]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, 1]]*tr[dot[cc[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]], v11, 1]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f23}]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[hv22, 1, v12]], xprog[-x, lor1, lor2]]]*tr[dot[cc[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]], v11, 1]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14]])/(8*CA^3*CF*(-1 + D)*D);

If[setd===4,dia8=FCI[dia8]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia8=FourierXP[dia8,{x,q}];

dia8=QEvaluate[I ScaleMu^(4-D) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype9[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia9,cc,tr,dot,lor1,lor2,lora,lorb,sun,sunSimplify,contract},


dia9=(Pi*Condensate["gg"]*Condensate[{f13, f24}]*Condensate[{f14, f23}]*contract[tr[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]*tr[dot[cc[dot[1, v12]], hv22, 1]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f14, f11}]*Condensate[{f21, f24}]*contract[tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], 1, cc[dot[v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]], 1]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c24, c21]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f21, f11}]*Condensate[{f22, f12}]*contract[tr[dot[1, cc[dot[v11, 1, hv21]]]]*tr[dot[cc[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[-x, lora, lorb], v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f23}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv22, 1, v12]], xprog[-x, lor1, lor2]]]*tr[dot[cc[dot[1, hv21]], v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14]])/(8*CA^3*CF*(-1 + D)*D);

If[setd===4,dia9=FCI[dia9]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia9=FourierXP[dia9,{x,q}];

dia9=QEvaluate[I ScaleMu^(4-D) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype10[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia10,cc,tr,dot,lor1,lor2,lora,lorb,sun,sunSimplify,contract},


dia10=(Pi*Condensate["gg"]*Condensate[{f14, f24}]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], 1]]*tr[dot[cc[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]], v11, 1]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13]])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f14, f12}]*Condensate[{f22, f23}]*contract[tr[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], cc[dot[hv22, 1, hv21]], xprog[-x, lor1, lor2], v12, cc[dot[v11, 1]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c12, c14]*ColorDelta[c23, c22]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f14, f12}]*Condensate[{f21, f23}]*contract[tr[dot[cc[dot[(-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]], v11, 1, cc[dot[xprog[-x, lor1, lor2], v12]], hv22, 1]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c12, c14]*ColorDelta[c23, c21]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13])])/(8*CA^3*CF*(-1 + D)*D) + (Pi*Condensate["gg"]*Condensate[{f13, f11}]*Condensate[{f21, f23}]*contract[tr[dot[cc[dot[xprog[-x, lor1, lor2]]], hv22, 1, cc[dot[v11, (-(MTD[lor1, lorb]*MTD[lor2, lora]) + MTD[lor1, lora]*MTD[lor2, lorb])*xprog[x, lora, lorb], hv21]], 1, v12]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c23, c21]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14])])/(8*CA^3*CF*(-1 + D)*D);

If[setd===4,dia10=FCI[dia10]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia10=FourierXP[dia10,{x,q}];

dia10=QEvaluate[I ScaleMu^(4-D) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
