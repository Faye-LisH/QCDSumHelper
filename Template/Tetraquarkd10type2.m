(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
Tetraquarkd10type2::usage="Tetraquarkd10type2[q_,j1_,j2_] gives the second type of \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)Gq\!\(\*SuperscriptBox[\(\[RightAngleBracket]\), \(2\)]\) contribution of \[LeftAngleBracket]j1 j2\[RightAngleBracket], with j1, j2 are tetraquark currents. "

Tetraquarkd10type2::inderr="Dummy indices conflict!"

Tetraquarkd10type2::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`Tetraquarkd10type2`"]


(* ::Code::Initialization::Plain:: *)
Options[Tetraquarkd10type2] = {
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
Tetraquarkd10type2[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 Tetraquarkd10type2[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])


(* ::Code::Initialization::Plain:: *)
Tetraquarkd10type2[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,holdf=OptionValue[HoldFlavor],atr,fdir,files,pall,diagrams,setd},



(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[Tetraquarkd10type2::curerr];
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
	Message[Tetraquarkd10type2::inderr];
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


diagrams={xtype1,xtype2,xtype3,xtype4,xtype5,xtype6,xtype7,xtype8,xtype9,xtype10,xtype11,xtype12,xtype13,xtype14,xtype15,xtype16,xtype17,xtype18,xtype19,xtype20,xtype21,xtype22,xtype23,xtype24};


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
			files=("Tetraquarkd10type2_"<>#)&/@files;
			
			
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
xtype1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia1,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia1=(Condensate[{f13, "G", f23}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], hv21]], v11, xprog[x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[cc[dot[xprog[x, lor1, lor2], hv21]], v11, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprop[x], hv21]], v11]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f23}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv21]], v11, xprop[x]]]*tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia1=FCI[dia1]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia1=FourierXP[dia1,{x,q}];

dia1=QEvaluate[I ScaleMu^(4-D) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia2,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia2=(Condensate[{f14, "G", f24}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f12}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv21]], v11, DiracSigma[GAD[lor1, lor2]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], v12]], hv22]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog[x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprop[x], hv21]], v11]]*tr[dot[cc[dot[hv22, DiracSigma[GAD[lor1, lor2]], v12]], xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia2=FCI[dia2]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia2=FourierXP[dia2,{x,q}];

dia2=QEvaluate[I ScaleMu^(4-D) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia3,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia3=(Condensate[{f14, "G", f23}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprog[x, lor1, lor2], hv21]], v11]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f22, "G", f12}]*contract[tr[dot[xprog[x, lor1, lor2], cc[dot[v11, hv21]]]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f11}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog[-x, lor1, lor2]]]*tr[dot[DiracSigma[GAD[lor1, lor2]], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f12}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], hv21]], v11]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia3=FCI[dia3]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia3=FourierXP[dia3,{x,q}];

dia3=QEvaluate[I ScaleMu^(4-D) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia4,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia4=(Condensate[{f14, "G", f24}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], hv21]], v11, xprog[x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]]]]*tr[dot[DiracSigma[GAD[lor1, lor2]], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lor1, lor2]], v12]]]]*tr[dot[xprog[x, lor1, lor2], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f23}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprop[x], hv21]], v11]]*tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia4=FCI[dia4]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia4=FourierXP[dia4,{x,q}];

dia4=QEvaluate[I ScaleMu^(4-D) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia5,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia5=(Condensate[{f21, "G", f11}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[v11, DiracSigma[GAD[lor1, lor2]], hv21]]]]*tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]]]]*tr[dot[xprop[x], cc[dot[v11, DiracSigma[GAD[lor1, lor2]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f23}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22]]*tr[dot[cc[dot[xprop[x], hv21]], v11, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia5=FCI[dia5]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia5=FourierXP[dia5,{x,q}];

dia5=QEvaluate[I ScaleMu^(4-D) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia6,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia6=(Condensate[{f13, "G", f23}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprog[x, lor1, lor2], hv21]], v11]]*tr[dot[cc[dot[hv22, DiracSigma[GAD[lor1, lor2]], v12]], xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f14, "G", f23}]*contract[tr[dot[xprog[x, lor1, lor2], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[v12]], hv22, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[hv22, v12]], xprog[-x, lor1, lor2]]]*tr[dot[DiracSigma[GAD[lor1, lor2]], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lor1, lor2]], v12]], xprog[-x, lor1, lor2]]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia6=FCI[dia6]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia6=FourierXP[dia6,{x,q}];

dia6=QEvaluate[I ScaleMu^(4-D) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia7,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia7=(Condensate[{f13, "G", f23}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[DiracSigma[GAD[lor1, lor2]], cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[v12]], hv22, DiracSigma[GAD[lor1, lor2]]]]*tr[dot[cc[dot[xprog[x, lor1, lor2], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f23}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f11}]*Condensate[{f22, "G", f12}]*contract[tr[dot[DiracSigma[GAD[lor1, lor2]], cc[dot[v11, hv21]]]]*tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia7=FCI[dia7]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia7=FourierXP[dia7,{x,q}];

dia7=QEvaluate[I ScaleMu^(4-D) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia8,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia8=(Condensate[{f14, "G", f24}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[xprog[x, lor1, lor2], hv21]], v11, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f12}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], hv21]], v11]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f23}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22]]*tr[dot[xprop[x], cc[dot[v11, DiracSigma[GAD[lor1, lor2]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f11}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[v11, DiracSigma[GAD[lor1, lor2]], hv21]]]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia8=FCI[dia8]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia8=FourierXP[dia8,{x,q}];

dia8=QEvaluate[I ScaleMu^(4-D) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype9[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia9,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia9=(Condensate[{f13, "G", f24}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], hv21]], v11, xprog[x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[xprog[x, lor1, lor2], cc[dot[v11, DiracSigma[GAD[lor1, lor2]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f12}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], xprop[-x]]]*tr[dot[cc[dot[hv21]], v11, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, v12]], DiracSigma[GAD[lor1, lor2]]]]*tr[dot[cc[dot[xprog[x, lor1, lor2], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia9=FCI[dia9]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia9=FourierXP[dia9,{x,q}];

dia9=QEvaluate[I ScaleMu^(4-D) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype10[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia10,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia10=(Condensate[{f13, "G", f23}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lor1, lor2]], v12]], xprop[-x]]]*tr[dot[cc[dot[hv21]], v11, xprog[x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprog[x, lor1, lor2], hv21]], v11]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprop[x], hv21]], v11]]*tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], DiracSigma[GAD[lor1, lor2]]]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia10=FCI[dia10]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia10=FourierXP[dia10,{x,q}];

dia10=QEvaluate[I ScaleMu^(4-D) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype11[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia11,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia11=(Condensate[{f13, "G", f23}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lor1, lor2]], v12]]]]*tr[dot[xprop[x], cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], DiracSigma[GAD[lor1, lor2]]]]*tr[dot[xprog[x, lor1, lor2], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], DiracSigma[GAD[lor1, lor2]]]]*tr[dot[cc[dot[hv21]], v11, xprog[x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, v12]], DiracSigma[GAD[lor1, lor2]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog[x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia11=FCI[dia11]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia11=FourierXP[dia11,{x,q}];

dia11=QEvaluate[I ScaleMu^(4-D) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype12[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia12,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia12=(Condensate[{f13, "G", f24}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f11}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[v11, DiracSigma[GAD[lor1, lor2]], hv21]]]]*tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lor1, lor2]], v12]], xprog[-x, lor1, lor2]]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], DiracSigma[GAD[lor1, lor2]]]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia12=FCI[dia12]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia12=FourierXP[dia12,{x,q}];

dia12=QEvaluate[I ScaleMu^(4-D) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype13[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia13,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia13=(Condensate[{f13, "G", f23}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]*tr[dot[cc[dot[hv22, DiracSigma[GAD[lor1, lor2]], v12]], xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f23}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv21]], v11, xprog[x, lor1, lor2]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f12}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], hv21]], v11]]*tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, v12]], xprog[-x, lor1, lor2]]]*tr[dot[xprop[x], cc[dot[v11, DiracSigma[GAD[lor1, lor2]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia13=FCI[dia13]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia13=FourierXP[dia13,{x,q}];

dia13=QEvaluate[I ScaleMu^(4-D) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype14[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia14,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia14=(Condensate[{f14, "G", f23}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[xprog[x, lor1, lor2], cc[dot[v11, DiracSigma[GAD[lor1, lor2]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f21, "G", f11}]*contract[tr[dot[DiracSigma[GAD[lor1, lor2]], cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv21]], v11, xprop[x]]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, v12]], DiracSigma[GAD[lor1, lor2]]]]*tr[dot[xprop[x], cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia14=FCI[dia14]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia14=FourierXP[dia14,{x,q}];

dia14=QEvaluate[I ScaleMu^(4-D) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype15[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia15,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia15=(Condensate[{f13, "G", f23}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lor1, lor2]], v12]]]]*tr[dot[cc[dot[xprog[x, lor1, lor2], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f23}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprog[x, lor1, lor2], hv21]], v11]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[hv22, v12]], xprog[-x, lor1, lor2]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia15=FCI[dia15]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia15=FourierXP[dia15,{x,q}];

dia15=QEvaluate[I ScaleMu^(4-D) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype16[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia16,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia16=(Condensate[{f13, "G", f23}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lor1, lor2]], v12]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog[x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[v12]], hv22, xprog[-x, lor1, lor2]]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f23}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, v12]], DiracSigma[GAD[lor1, lor2]]]]*tr[dot[xprog[x, lor1, lor2], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia16=FCI[dia16]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia16=FourierXP[dia16,{x,q}];

dia16=QEvaluate[I ScaleMu^(4-D) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype17[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia17,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia17=(Condensate[{f21, "G", f12}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], hv21]], v11]]*tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f23}]*Condensate[{f22, "G", f12}]*contract[tr[dot[xprog[x, lor1, lor2], cc[dot[v11, hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f11}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[v11, DiracSigma[GAD[lor1, lor2]], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia17=FCI[dia17]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia17=dia17/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia17=FourierXP[dia17,{x,q}];

dia17=QEvaluate[I ScaleMu^(4-D) dia17,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype18[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia18,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia18=(Condensate[{f13, "G", f24}]*Condensate[{f22, "G", f12}]*contract[tr[dot[xprog[x, lor1, lor2], cc[dot[v11, DiracSigma[GAD[lor1, lor2]], hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[xprog[x, lor1, lor2], cc[dot[v11, DiracSigma[GAD[lor1, lor2]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f11}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]], xprop[-x]]]*tr[dot[DiracSigma[GAD[lor1, lor2]], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[hv22, DiracSigma[GAD[lor1, lor2]], v12]], xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia18=FCI[dia18]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia18=dia18/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia18=FourierXP[dia18,{x,q}];

dia18=QEvaluate[I ScaleMu^(4-D) dia18,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype19[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia19,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia19=(Condensate[{f13, "G", f24}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprog[x, lor1, lor2], hv21]], v11, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f11}]*Condensate[{f22, "G", f12}]*contract[tr[dot[DiracSigma[GAD[lor1, lor2]], cc[dot[v11, hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f22, "G", f12}]*contract[tr[dot[xprop[x], cc[dot[v11, DiracSigma[GAD[lor1, lor2]], hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv21]], v11, xprog[x, lor1, lor2]]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia19=FCI[dia19]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia19=dia19/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia19=FourierXP[dia19,{x,q}];

dia19=QEvaluate[I ScaleMu^(4-D) dia19,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype20[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia20,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia20=(Condensate[{f13, "G", f23}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lor1, lor2]], v12]], xprop[-x]]]*tr[dot[xprog[x, lor1, lor2], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f14, "G", f23}]*contract[tr[dot[xprop[x], cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]*tr[dot[cc[dot[v12]], hv22, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], v12]], hv22]]*tr[dot[cc[dot[xprog[x, lor1, lor2], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia20=FCI[dia20]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia20=dia20/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia20=FourierXP[dia20,{x,q}];

dia20=QEvaluate[I ScaleMu^(4-D) dia20,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype21[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia21,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia21=(Condensate[{f14, "G", f23}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], hv21]], v11, xprog[x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[v12]], hv22, xprog[-x, lor1, lor2]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f23}]*Condensate[{f22, "G", f12}]*contract[tr[dot[xprop[x], cc[dot[v11, hv21]]]]*tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f23}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, v12]], xprog[-x, lor1, lor2]]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia21=FCI[dia21]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia21=dia21/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia21=FourierXP[dia21,{x,q}];

dia21=QEvaluate[I ScaleMu^(4-D) dia21,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype22[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia22,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia22=(Condensate[{f14, "G", f23}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[cc[dot[xprog[x, lor1, lor2], hv21]], v11, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, xprog[-x, lor1, lor2], v12]]]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f21, "G", f11}]*contract[tr[dot[DiracSigma[GAD[lor1, lor2]], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f12}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv21]], v11, DiracSigma[GAD[lor1, lor2]]]]*tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia22=FCI[dia22]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia22=dia22/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia22=FourierXP[dia22,{x,q}];

dia22=QEvaluate[I ScaleMu^(4-D) dia22,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype23[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia23,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia23=(Condensate[{f14, "G", f23}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[DiracSigma[GAD[lor1, lor2]], cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f14, "G", f24}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[DiracSigma[GAD[lor1, lor2]], cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[v12]], hv22, DiracSigma[GAD[lor1, lor2]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog[x, lor1, lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], v12]], hv22]]*tr[dot[xprog[x, lor1, lor2], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia23=FCI[dia23]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia23=dia23/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia23=FourierXP[dia23,{x,q}];

dia23=QEvaluate[I ScaleMu^(4-D) dia23,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype24[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia24,cc,tr,dot,lor1,lor2,sun,sunSimplify,contract},


dia24=(Condensate[{f14, "G", f23}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[xprog[-x, lor1, lor2], v12]], hv22]]*tr[dot[DiracSigma[GAD[lor1, lor2]], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], v12]], hv22]]*tr[dot[xprop[x], cc[dot[v11, xprog[x, lor1, lor2], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f13, "G", f24}]*Condensate[{f22, "G", f12}]*contract[tr[dot[xprop[x], cc[dot[v11, hv21]]]]*tr[dot[cc[dot[DiracSigma[GAD[lor1, lor2]], v12]], hv22, xprog[-x, lor1, lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SPD[x]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]])/(64*CA^2*CF*(-1 + D)*D^2) + (Condensate[{f21, "G", f12}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog[-x, lor1, lor2]]]*tr[dot[cc[dot[hv21]], v11, DiracSigma[GAD[lor1, lor2]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SPD[x]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14]])/(64*CA^2*CF*(-1 + D)*D^2);

If[setd===4,dia24=FCI[dia24]/.dg_DiracGamma:>(dg/.D->4)/.ld_LorentzIndex:>(ld/.D->4)];


dia24=dia24/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia24=FourierXP[dia24,{x,q}];

dia24=QEvaluate[I ScaleMu^(4-D) dia24,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
