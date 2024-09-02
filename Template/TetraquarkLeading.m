(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


TetraquarkLeading::usage="TetraquarkLeading[q_,j1_,j2_] gives leading Order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents. "

TetraquarkLeading::inderr="Dummy indices conflict!"

TetraquarkLeading::curerr="Unknow current structure!"


Begin["`Private`TetraquarkLeading`"]

Options[TetraquarkLeading] = {
	EpsOrder->0,
	HoldFlavor->False,
	AutoNDR->True,
	Parallelized->False,
	Pole->1
}




(*------------------------------------------------------------------*)
(* propagators *)

prop[q_] = I GSD[q]FAD[q];


(* xprop[x_]=FourierPX[prop[q],{q,x}]; *)
xprop[x_] = 1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];



(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
TetraquarkLeading[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 TetraquarkLeading[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])
TetraquarkLeading[qq_,current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor2 TetraquarkLeading[qq,current1,current2,ops]/;FreeQ[factor2,Current]
TetraquarkLeading[qq_,factor1_ current1_FourquarkCurrent,current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 TetraquarkLeading[qq,current1,current2,ops]/;FreeQ[factor1,Current]


TetraquarkLeading[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,trs,dia,x,q,v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,
v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,ord=OptionValue[EpsOrder],holdf=OptionValue[HoldFlavor],pole=OptionValue[Pole],atr,pall,fdir,files},



(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[TetraquarkLeading::curerr];
	Abort[]
];


(*--- check dummy indices in input ---*)
tmp2=Length[#]&/@Gather[({c11,c12,c13,c14,c21,c22,c23,c24}//FCI)/.SUNFIndex[aa_]:>aa];

If[DummyCheck[v11,v12,v21,v22]==0||!FreeQ[tmp2,nn_Integer/;nn>2],
	Message[TetraquarkAlphasType1::inderr];
	Abort[]
];



(*-------------------------------*)
(* B A^+ B *)
hv21=v21//ComplexConjugate;
hv22=v22//ComplexConjugate;

(*-------------------------------*)

If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];

(*-------------------------------*)

If[OptionValue[AutoNDR]===True,atr=TR5,atr=TR];

(*-------------------------------*)
If[pall===True,

	DistributeDefinitions[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,ord,pole];
	tmp= WaitAll[ParallelSubmit[leading[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,ord,pole]]];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,
	If[pall==="External",

	DistributeDefinitions[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,ord,pole];
	
		If[fdir==="None",
			(* evaluation, no import and export *)
			 ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,ord,pole},
									leading[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,ord,pole]]

		
		,
		(* evaluation, import and export *)
			files=(StringSplit[ToString[#],"`"][[-1]])&/@{leading};
			files=("TetraquarkLeading_"<>#)&/@files;
			
			
			ImExport[fdir,
						files,
						{{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,ord,pole},
						{qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,ord,pole},
						{leading}}
						]
				
		](* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
	,
	
		
		tmp=leading[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,ord,pole];
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]

]


(*------------------------------------------------------------------*)

fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]
(*------------------------------------------------------------------*)


leading[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,ord_,pole_]:=Block[{x,q,dia1,cc,tr,dot,sunSimplify,contract},


dia1=contract[tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + contract[tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract,sunSimplify->SUNSimplify};

dia1=FourierXP[dia1,{x,q}];

dia1=QEvaluate[I ScaleMu^(3(4-D)) pole dia1, q, HoldFlavor->holdf, EpsOrder->ord, Parallelized->False]/.q->qq/.CA-2CF->1/CA
]


(*leading[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,ord_,pole_]:=Block[{x,q,trs,dia},

(*------------------------------------------------------------------*)
(* figure *)
trs=(FlavorDelta[f11,f21]FlavorDelta[f12,f22]ColorDelta[c11,c21]ColorDelta[c12,c22]atr[v11 . xprop[x] . hv21 . FCChargeConjugateTransposed[xprop[x]]]
	-FlavorDelta[f11,f22]FlavorDelta[f12,f21]ColorDelta[c11,c22]ColorDelta[c12,c21]atr[v11 . xprop[x] . FCChargeConjugateTransposed[xprop[x] . hv21]])
	(FlavorDelta[f13,f23]FlavorDelta[f14,f24]ColorDelta[c13,c23]ColorDelta[c14,c24]atr[hv22 . xprop[-x] . v12 . FCChargeConjugateTransposed[xprop[-x]]]
	-FlavorDelta[f13,f24]FlavorDelta[f14,f23]ColorDelta[c13,c24]ColorDelta[c14,c23]atr[xprop[-x] . v12 . FCChargeConjugateTransposed[hv22 . xprop[-x]]]);
	
	
dia=FourierXP[trs,{x,q}];
dia=QEvaluate[I ScaleMu^(3(4-D)) dia pole, q, HoldFlavor->holdf,EpsOrder->ord,Parallelized->False];


(*------------------------------------------------------------------*)
QGather[dia//SUNSimplify,q,ShowasTable->False]/.q->qq
]*)


End[]
(*EndPackage[]*)
