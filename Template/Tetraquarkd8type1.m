(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
Tetraquarkd8type1::usage="Tetraquarkd8type1[q_,j1_,j2_] gives the first type of \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)q\[RightAngleBracket]\[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)Gq\[RightAngleBracket] contribution of \[LeftAngleBracket]j1 j2\[RightAngleBracket], with j1, j2 are tetraquark currents. "

Tetraquarkd8type1::inderr="Dummy indices conflict!"

Tetraquarkd8type1::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`Tetraquarkd8type1`"]


(* ::Code::Initialization::Plain:: *)
Options[Tetraquarkd8type1] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True
}


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

prop[q_]=I GSD[q]FAD[q];


(* ::Code::Initialization::Plain:: *)
xprop[x_]= 1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];(*FourierPX[prop[q],{q,x}]*)

(*--- propagator with one background gluon ---*)
(*prog[q_,lora_,lorb_,cola_]=-1/2 FAD[{q}]GSD[q].GAD[lorb].SUNT[cola].FourDivergence[GSD[q]FAD[{q}],FVD[q,lora]];*)
prog[q_,lora_,lorb_]= -1/2 FAD[{q}]GSD[q] . GAD[lorb] . (FAD[q] GAD[lora]-2FAD[q,q] FVD[q,lora] GSD[q]);

(*xprog[x_,lora_,lorb_,cola_]=FourierPX[prog[q,lora,lorb,cola],{q,x}];*)
xprog[x_,lora_,lorb_]=(GAD[lora] . GAD[lorb]. GSD[x] qfact1[1/32 I^-D \[Pi]^(-D/2) qGamma[-1+D/2]] SPD[x,x]^(1-D/2) 
							+GSD[x] . GAD[lorb] . GAD[lora] qfact1[qfact2[-(1/32) I^-D \[Pi]^(-D/2)] qGamma[-1+D/2]] SPD[x,x]^(1-D/2) 
							+2 FVD[x,lora] FVD[x,lorb] GSD[x] qfact1[-(1/16) I^-D \[Pi]^(-D/2) qGamma[D/2]] SPD[x,x]^(-D/2) );

(*--- Condensate ---*)

qgq[f1_,f2_]:=-1/(4CA CF (D-1)D) Condensate[{f1,"G",f2}];
cqq[f1_,f2_]:=-1/(4CA)Condensate[{f1,f2}];		


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]



(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
Tetraquarkd8type1[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 Tetraquarkd8type1[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])




Tetraquarkd8type1[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,tr1a,tr1b,tr2a1,tr2a2,tr2b1,tr2b2,tr3a1,tr3a2,tr3b1,tr3b2,dia1,dia2,dia3,x,q,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,holdf=OptionValue[HoldFlavor],atr,fdir,files,pall,diagrams},


(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[Tetraquarkd8type1::curerr];
	Abort[]
];


If[OptionValue[AutoNDR]===True,
	atr=TR5
,
	atr=TR
];


(*--- check dummy indices in input ---*)
tmp2=Length[#]&/@Gather[({c11,c12,c13,c14,c21,c22,c23,c24}//FCI)/.SUNFIndex[aa_]:>aa];

If[DummyCheck[v11,v12,v21,v22]==0||!FreeQ[tmp2,nn_Integer/;nn>2],
	Message[TetraquarkAlphasType2::inderr];
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


diagrams={type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16};


(*---------------------------------------------------*)
If[pall===True,

	
	tmp=Plus@@WaitAll[ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr]]&/@diagrams];
	
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

,
	If[pall==="External",

	DistributeDefinitions[qq,current1,current2];
	
		If[fdir==="None",
			(* evaluation, no import and export *)
			 ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr]]&/@diagrams
		
		,
		(* evaluation, import and export *)
			files=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files=("Tetraquarkd8type1_"<>#)&/@files;
			
			
			ImExport[fdir,
						files,
						{{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
						{qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
						diagrams}
						]
				
		](* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
	
	,
	
		
		tmp=Plus@@(#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr]&/@diagrams);
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]


]
(*------------------------------------------------------------------*)



(* ::Code::Initialization::Plain:: *)
(*--- The following are generated by algorithm ---*)


(* ::Input::Initialization:: *)
type1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,sun,lora,lorb},


dia1=-(atr[fcc[hv21] . hv22 . xprop[-x] . fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . v11 . xprog[x, lora, lorb]]*ColorDelta[c14, c23]*ColorDelta[c22, c24]*cqq[f24, f22]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f14, f23]*FlavorDelta[f22, f24]*qgq[f11, f13]*SUNTF[{sun}, c11, c13]*SUNTF[{sun}, c12, c21]) - atr[fcc[xprop[-x]] . hv22 . fcc[v11 . xprog[x, lora, lorb] . hv21] . DiracSigma[GAD[lora], GAD[lorb]] . v12]*ColorDelta[c14, c24]*ColorDelta[c21, c23]*cqq[f23, f21]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f14, f24]*FlavorDelta[f21, f23]*qgq[f11, f13]*SUNTF[{sun}, c11, c13]*SUNTF[{sun}, c12, c22] + atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . xprop[x]]*atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . v12]]*ColorDelta[c12, c21]*ColorDelta[c14, c24]*cqq[f24, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f23, f13]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c23, c13] + atr[fcc[xprop[x] . hv21] . v11 . xprog[x, lora, lorb]]*atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . v12]]*ColorDelta[c11, c22]*ColorDelta[c14, c24]*cqq[f24, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f23, f13]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c23, c13] + atr[xprop[x] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . hv21]]*atr[fcc[v12] . hv22 . xprog[-x, lora, lorb]]*ColorDelta[c11, c21]*ColorDelta[c13, c24]*cqq[f24, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f12, f22]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c23, c14] + atr[xprog[x, lora, lorb] . fcc[v11 . hv21]]*atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . hv22 . xprop[-x]]*ColorDelta[c12, c22]*ColorDelta[c14, c23]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f24, f13]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c24, c13] - atr[fcc[xprop[x] . hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]] . fcc[xprog[-x, lora, lorb] . v12] . hv22]*ColorDelta[c11, c22]*ColorDelta[c21, c23]*cqq[f23, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f13, f24]*FlavorDelta[f21, f23]*qgq[f12, f14]*SUNTF[{sun}, c12, c14]*SUNTF[{sun}, c24, c13] + atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[v11 . xprop[x] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c12, c22]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f23, f14]*SUNTF[{sun}, c23, c14]*SUNTF[{sun}, c24, c13] - atr[fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora], GAD[lorb]] . hv21 . fcc[v11 . xprog[x, lora, lorb]]]*ColorDelta[c11, c14]*ColorDelta[c13, c23]*cqq[f11, f14]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f22, f24]*qgq[f24, f22]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c24, c22] - atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . hv22 . xprog[-x, lora, lorb] . fcc[v12] . v11 . xprop[x]]*ColorDelta[c11, c13]*ColorDelta[c12, c21]*cqq[f11, f13]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f14, f23]*FlavorDelta[f22, f24]*qgq[f24, f22]*SUNTF[{sun}, c23, c14]*SUNTF[{sun}, c24, c22];


dia1=(FourierXP[dia1,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia1=QEvaluate[I ScaleMu^(4-D) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,sun,lora,lorb},


dia2=-(atr[1 . hv22 . xprop[-x] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprog[x, lora, lorb] . hv21]*ColorDelta[c14, c23]*ColorDelta[c21, c24]*cqq[f24, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f14, f23]*FlavorDelta[f21, f24]*qgq[f12, f13]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c12, c13]) + atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . v11]*atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . xprop[-x]]*ColorDelta[c12, c21]*ColorDelta[c14, c24]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f11, f22]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c23, c13] + atr[xprog[x, lora, lorb] . fcc[v11 . xprop[x] . hv21]]*atr[fcc[v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c12, c22]*ColorDelta[c13, c24]*cqq[f24, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f23, f14]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c23, c14] + atr[fcc[hv21] . v11 . xprog[x, lora, lorb]]*atr[fcc[xprop[-x] . v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c22]*ColorDelta[c13, c24]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f23, f14]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c23, c14] + atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprop[x] . hv21]]*atr[fcc[xprog[-x, lora, lorb] . v12] . hv22]*ColorDelta[c12, c22]*ColorDelta[c14, c23]*cqq[f23, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f11, f21]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c24, c13] + atr[fcc[xprog[-x, lora, lorb] . v12] . hv22]*atr[fcc[xprop[x] . hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c22]*ColorDelta[c14, c23]*cqq[f23, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f12, f21]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c24, c13] + atr[fcc[hv21] . v11 . xprop[x]]*atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c22]*ColorDelta[c12, c21]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f23, f14]*SUNTF[{sun}, c23, c14]*SUNTF[{sun}, c24, c13] + atr[fcc[xprop[x] . hv21] . v11]*atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . hv22 . xprog[-x, lora, lorb]]*ColorDelta[c11, c22]*ColorDelta[c12, c21]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f24, f13]*SUNTF[{sun}, c23, c14]*SUNTF[{sun}, c24, c13] - atr[fcc[hv22] . xprog[-x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprop[x] . hv21]*ColorDelta[c11, c22]*ColorDelta[c21, c23]*cqq[f23, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f14, f24]*FlavorDelta[f21, f23]*qgq[f12, f13]*SUNTF[{sun}, c12, c13]*SUNTF[{sun}, c24, c14] - atr[fcc[hv22 . hv21] . xprog[-x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprop[x]]*ColorDelta[c11, c21]*ColorDelta[c22, c23]*cqq[f23, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f14, f24]*FlavorDelta[f22, f23]*qgq[f12, f13]*SUNTF[{sun}, c12, c13]*SUNTF[{sun}, c24, c14];


dia2=(FourierXP[dia2,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia2=QEvaluate[I ScaleMu^(4-D) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,sun,lora,lorb},


dia3=-(atr[fcc[xprop[-x] . v12] . hv22 . fcc[v11 . xprog[x, lora, lorb] . hv21] . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c13, c24]*ColorDelta[c21, c23]*cqq[f23, f21]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f21, f23]*qgq[f11, f14]*SUNTF[{sun}, c11, c14]*SUNTF[{sun}, c12, c22]) + atr[fcc[xprop[-x] . v12] . hv22]*atr[xprog[x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . hv21]]*ColorDelta[c13, c24]*ColorDelta[c14, c23]*cqq[f23, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f12, f22]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c12, c22] + atr[xprog[x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . hv21]]*atr[fcc[v12] . hv22 . xprop[-x]]*ColorDelta[c13, c24]*ColorDelta[c14, c23]*cqq[f24, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f12, f22]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c12, c22] + atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . xprop[-x]]*atr[fcc[hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c22]*ColorDelta[c14, c24]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f12, f21]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c23, c13] + atr[fcc[v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[xprop[x] . hv21] . v11 . xprog[x, lora, lorb]]*ColorDelta[c11, c22]*ColorDelta[c13, c24]*cqq[f24, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f23, f14]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c23, c14] + atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . hv21]]*atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . xprop[-x]]*ColorDelta[c12, c22]*ColorDelta[c14, c23]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f11, f21]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c24, c13] + atr[fcc[xprog[-x, lora, lorb] . v12] . hv22]*atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . v11 . xprop[x]]*ColorDelta[c12, c21]*ColorDelta[c14, c23]*cqq[f23, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f11, f22]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c24, c13] + atr[fcc[xprog[x, lora, lorb] . hv21] . v11]*atr[fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c12, c21]*ColorDelta[c13, c23]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f24, f14]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . v12] . xprog[-x, lora, lorb]]*atr[fcc[xprop[x] . hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c22]*ColorDelta[c13, c23]*cqq[f23, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f12, f21]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c24, c14] - atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21] . v12 . fcc[hv22 . xprop[-x]]]*ColorDelta[c11, c13]*ColorDelta[c14, c23]*cqq[f11, f13]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f14, f23]*FlavorDelta[f21, f24]*qgq[f24, f21]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c24, c21];


dia3=(FourierXP[dia3,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia3=QEvaluate[I ScaleMu^(4-D) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,sun,lora,lorb},


dia4=atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprop[x] . hv21]]*atr[fcc[hv22 . xprog[-x, lora, lorb] . v12]]*ColorDelta[c12, c22]*ColorDelta[c14, c24]*cqq[f24, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f11, f21]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c23, c13] + atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprop[-x]]*atr[fcc[hv21] . v11 . xprog[x, lora, lorb]]*ColorDelta[c11, c22]*ColorDelta[c14, c24]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f23, f13]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c23, c13] + atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . hv21]]*atr[fcc[xprop[-x] . v12] . hv22 . xprog[-x, lora, lorb]]*ColorDelta[c12, c22]*ColorDelta[c13, c24]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f11, f21]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c23, c14] + atr[fcc[xprog[x, lora, lorb] . hv21] . v11]*atr[fcc[xprop[-x] . v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c12, c21]*ColorDelta[c13, c24]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f23, f14]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c23, c14] + atr[fcc[v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . xprop[x]]*ColorDelta[c12, c21]*ColorDelta[c13, c24]*cqq[f24, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f23, f14]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c23, c14] + atr[fcc[hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[xprop[-x] . v12] . hv22 . xprog[-x, lora, lorb]]*ColorDelta[c11, c22]*ColorDelta[c13, c24]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f12, f21]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c23, c14] + atr[xprop[x] . fcc[v11 . hv21]]*atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c21]*ColorDelta[c12, c22]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f23, f14]*SUNTF[{sun}, c23, c14]*SUNTF[{sun}, c24, c13] - atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprop[x] . hv21]]*ColorDelta[c11, c14]*ColorDelta[c12, c22]*cqq[f11, f14]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f21, f23]*qgq[f23, f21]*SUNTF[{sun}, c23, c21]*SUNTF[{sun}, c24, c13] - atr[1 . v11 . xprop[x] . fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . hv21] . xprog[-x, lora, lorb] . v12]*ColorDelta[c11, c14]*ColorDelta[c12, c21]*cqq[f11, f14]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f22, f23]*qgq[f23, f22]*SUNTF[{sun}, c23, c22]*SUNTF[{sun}, c24, c13] - atr[fcc[xprop[x] . hv21] . v11 . fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c22]*ColorDelta[c12, c14]*cqq[f12, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f13, f23]*FlavorDelta[f21, f24]*qgq[f24, f21]*SUNTF[{sun}, c23, c13]*SUNTF[{sun}, c24, c21];


dia4=(FourierXP[dia4,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia4=QEvaluate[I ScaleMu^(4-D) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia5,sun,lora,lorb},


dia5=-(atr[fcc[v11 . xprog[x, lora, lorb] . hv21] . DiracSigma[GAD[lora], GAD[lorb]] . v12 . fcc[hv22 . xprop[-x]]]*ColorDelta[c14, c23]*ColorDelta[c21, c24]*cqq[f24, f21]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f14, f23]*FlavorDelta[f21, f24]*qgq[f11, f13]*SUNTF[{sun}, c11, c13]*SUNTF[{sun}, c12, c22]) + atr[xprog[x, lora, lorb] . fcc[v11 . xprop[x] . hv21]]*atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . v12]]*ColorDelta[c12, c22]*ColorDelta[c14, c24]*cqq[f24, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f23, f13]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c23, c13] - atr[fcc[xprop[x] . hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]] . fcc[hv22 . xprog[-x, lora, lorb] . v12]]*ColorDelta[c11, c22]*ColorDelta[c21, c24]*cqq[f24, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f13, f23]*FlavorDelta[f21, f24]*qgq[f12, f14]*SUNTF[{sun}, c12, c14]*SUNTF[{sun}, c23, c13] + atr[fcc[xprop[x] . hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[hv22 . xprog[-x, lora, lorb] . v12]]*ColorDelta[c11, c22]*ColorDelta[c14, c24]*cqq[f24, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f12, f21]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c23, c13] + atr[fcc[v12] . hv22 . xprog[-x, lora, lorb]]*atr[fcc[xprop[x] . hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c22]*ColorDelta[c13, c24]*cqq[f24, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f12, f21]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c23, c14] - atr[fcc[xprop[-x] . v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21]]*ColorDelta[c11, c14]*ColorDelta[c13, c24]*cqq[f11, f14]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f21, f23]*qgq[f23, f21]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c23, c21] + atr[fcc[hv22 . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*atr[xprog[x, lora, lorb] . fcc[v11 . xprop[x] . hv21]]*ColorDelta[c12, c22]*ColorDelta[c13, c23]*cqq[f23, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f24, f14]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c24, c14] + atr[fcc[xprop[x] . hv21] . v11]*atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprog[-x, lora, lorb]]*ColorDelta[c11, c22]*ColorDelta[c12, c21]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f23, f13]*SUNTF[{sun}, c23, c13]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprog[-x, lora, lorb]]*atr[fcc[v11 . xprop[x] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c12, c22]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f23, f13]*SUNTF[{sun}, c23, c13]*SUNTF[{sun}, c24, c14] - atr[fcc[xprop[x]] . v11 . fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora], GAD[lorb]] . hv21]*ColorDelta[c11, c21]*ColorDelta[c12, c14]*cqq[f12, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f13, f23]*FlavorDelta[f22, f24]*qgq[f24, f22]*SUNTF[{sun}, c23, c13]*SUNTF[{sun}, c24, c22];


dia5=(FourierXP[dia5,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia5=QEvaluate[I ScaleMu^(4-D) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia6,sun,lora,lorb},


dia6=-(atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]] . fcc[hv22 . xprop[-x] . v12]]*ColorDelta[c13, c23]*ColorDelta[c21, c24]*cqq[f24, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f13, f23]*FlavorDelta[f21, f24]*qgq[f12, f14]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c12, c14]) - atr[fcc[DiracSigma[GAD[lora], GAD[lorb]]] . v11 . xprog[x, lora, lorb] . fcc[hv22 . hv21] . xprop[-x] . v12]*ColorDelta[c13, c24]*ColorDelta[c22, c23]*cqq[f23, f22]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f22, f23]*qgq[f11, f14]*SUNTF[{sun}, c11, c14]*SUNTF[{sun}, c12, c21] + atr[fcc[xprop[-x] . v12] . hv22]*atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . v11 . xprog[x, lora, lorb]]*ColorDelta[c13, c24]*ColorDelta[c14, c23]*cqq[f23, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f11, f22]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c12, c21] + atr[fcc[v12] . hv22 . xprop[-x]]*atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . v11 . xprog[x, lora, lorb]]*ColorDelta[c13, c24]*ColorDelta[c14, c23]*cqq[f24, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f11, f22]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c12, c21] + atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21]]*atr[fcc[xprop[-x] . v12] . hv22]*ColorDelta[c13, c24]*ColorDelta[c14, c23]*cqq[f23, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f11, f21]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c12, c22] + atr[xprop[x] . fcc[v11 . xprog[x, lora, lorb] . hv21]]*atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . v12]]*ColorDelta[c11, c21]*ColorDelta[c14, c24]*cqq[f24, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f23, f13]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c23, c13] - atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . fcc[xprop[-x] . v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c12, c14]*ColorDelta[c13, c24]*cqq[f12, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f13, f24]*FlavorDelta[f21, f23]*qgq[f23, f21]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c23, c21] + atr[fcc[xprog[-x, lora, lorb] . v12] . hv22]*atr[xprop[x] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c14, c23]*cqq[f23, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f12, f22]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c24, c13] + atr[fcc[hv22 . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*atr[xprop[x] . fcc[v11 . xprog[x, lora, lorb] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c13, c23]*cqq[f23, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f24, f14]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c24, c14] - atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprop[x] . hv21]]*ColorDelta[c11, c14]*ColorDelta[c12, c22]*cqq[f11, f14]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f21, f24]*qgq[f24, f21]*SUNTF[{sun}, c23, c13]*SUNTF[{sun}, c24, c21];


dia6=(FourierXP[dia6,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia6=QEvaluate[I ScaleMu^(4-D) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia7,sun,lora,lorb},


dia7=-(atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]] . fcc[xprop[-x] . v12] . hv22]*ColorDelta[c13, c24]*ColorDelta[c21, c23]*cqq[f23, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f13, f24]*FlavorDelta[f21, f23]*qgq[f12, f14]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c12, c14]) - atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . v11 . xprog[x, lora, lorb] . fcc[hv22 . hv21] . xprop[-x]]*ColorDelta[c14, c24]*ColorDelta[c22, c23]*cqq[f23, f22]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f14, f24]*FlavorDelta[f22, f23]*qgq[f11, f13]*SUNTF[{sun}, c11, c13]*SUNTF[{sun}, c12, c21] + atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21]]*atr[fcc[hv22 . xprop[-x] . v12]]*ColorDelta[c13, c23]*ColorDelta[c14, c24]*cqq[f24, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f11, f21]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c12, c22] + atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . hv21]]*atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . xprop[-x]]*ColorDelta[c12, c22]*ColorDelta[c14, c24]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f11, f21]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c23, c13] + atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprop[-x]]*atr[xprog[x, lora, lorb] . fcc[v11 . hv21]]*ColorDelta[c12, c22]*ColorDelta[c14, c24]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f23, f13]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c23, c13] - atr[fcc[v11 . xprop[x] . hv21] . DiracSigma[GAD[lora], GAD[lorb]] . v12 . fcc[hv22 . xprog[-x, lora, lorb]]]*ColorDelta[c12, c22]*ColorDelta[c21, c24]*cqq[f24, f21]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f14, f23]*FlavorDelta[f21, f24]*qgq[f11, f13]*SUNTF[{sun}, c11, c13]*SUNTF[{sun}, c23, c14] + atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprop[x] . hv21]]*atr[fcc[v12] . hv22 . xprog[-x, lora, lorb]]*ColorDelta[c12, c22]*ColorDelta[c13, c24]*cqq[f24, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f11, f21]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c23, c14] - atr[fcc[hv21] . hv22 . xprog[-x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprop[x]]*ColorDelta[c11, c21]*ColorDelta[c22, c24]*cqq[f24, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f14, f23]*FlavorDelta[f22, f24]*qgq[f12, f13]*SUNTF[{sun}, c12, c13]*SUNTF[{sun}, c23, c14] - atr[fcc[v12] . v11 . xprog[x, lora, lorb] . fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . hv21] . xprop[-x]]*ColorDelta[c11, c13]*ColorDelta[c14, c24]*cqq[f11, f13]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f14, f24]*FlavorDelta[f22, f23]*qgq[f23, f22]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c23, c22] + atr[fcc[hv21] . v11 . xprop[x]]*atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . hv22 . xprog[-x, lora, lorb]]*ColorDelta[c11, c22]*ColorDelta[c12, c21]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f24, f13]*SUNTF[{sun}, c23, c14]*SUNTF[{sun}, c24, c13];


dia7=(FourierXP[dia7,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia7=QEvaluate[I ScaleMu^(4-D) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia8,sun,lora,lorb},


dia8=-(atr[fcc[hv22] . xprop[-x] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprog[x, lora, lorb] . hv21]*ColorDelta[c14, c24]*ColorDelta[c21, c23]*cqq[f23, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f14, f24]*FlavorDelta[f21, f23]*qgq[f12, f13]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c12, c13]) - atr[fcc[xprog[x, lora, lorb]] . v11 . DiracSigma[GAD[lora], GAD[lorb]] . fcc[hv22 . xprop[-x] . v12] . hv21]*ColorDelta[c13, c23]*ColorDelta[c22, c24]*cqq[f24, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f13, f23]*FlavorDelta[f22, f24]*qgq[f12, f14]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c12, c14] - atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . fcc[v11 . xprop[x] . hv21] . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c12, c22]*ColorDelta[c21, c24]*cqq[f24, f21]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f21, f24]*qgq[f11, f14]*SUNTF[{sun}, c11, c14]*SUNTF[{sun}, c23, c13] - atr[fcc[xprop[x]] . v11 . DiracSigma[GAD[lora], GAD[lorb]] . fcc[hv22 . xprog[-x, lora, lorb] . v12] . hv21]*ColorDelta[c11, c21]*ColorDelta[c22, c24]*cqq[f24, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f13, f23]*FlavorDelta[f22, f24]*qgq[f12, f14]*SUNTF[{sun}, c12, c14]*SUNTF[{sun}, c23, c13] + atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprop[-x]]*atr[fcc[v11 . xprog[x, lora, lorb] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c14, c24]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f23, f13]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c23, c13] + atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . v11]*atr[fcc[xprop[-x] . v12] . hv22 . xprog[-x, lora, lorb]]*ColorDelta[c12, c21]*ColorDelta[c13, c24]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f11, f22]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c23, c14] - atr[xprog[x, lora, lorb] . fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . hv21] . xprop[-x] . v12 . fcc[v11]]*ColorDelta[c12, c14]*ColorDelta[c13, c24]*cqq[f12, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f13, f24]*FlavorDelta[f22, f23]*qgq[f23, f22]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c23, c22] + atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . v11]*atr[fcc[hv22 . xprop[-x] . v12] . xprog[-x, lora, lorb]]*ColorDelta[c12, c21]*ColorDelta[c13, c23]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f11, f22]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c24, c14] - atr[fcc[xprog[x, lora, lorb]] . v11 . fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora], GAD[lorb]] . hv21]*ColorDelta[c12, c14]*ColorDelta[c13, c23]*cqq[f12, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f13, f23]*FlavorDelta[f22, f24]*qgq[f24, f22]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c24, c22] - atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . hv22 . xprop[-x] . fcc[v12] . v11 . xprog[x, lora, lorb]]*ColorDelta[c11, c13]*ColorDelta[c14, c23]*cqq[f11, f13]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f14, f23]*FlavorDelta[f22, f24]*qgq[f24, f22]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c24, c22];


dia8=(FourierXP[dia8,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia8=QEvaluate[I ScaleMu^(4-D) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type9[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia9,sun,lora,lorb},


dia9=-(atr[fcc[hv22 . hv21] . xprop[-x] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprog[x, lora, lorb]]*ColorDelta[c14, c24]*ColorDelta[c22, c23]*cqq[f23, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f14, f24]*FlavorDelta[f22, f23]*qgq[f12, f13]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c12, c13]) - atr[xprog[x, lora, lorb] . fcc[hv22 . hv21] . xprop[-x] . v12 . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]]]]*ColorDelta[c13, c24]*ColorDelta[c22, c23]*cqq[f23, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f13, f24]*FlavorDelta[f22, f23]*qgq[f12, f14]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c12, c14] - atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[hv22 . xprop[-x] . v12] . hv21 . fcc[v11 . xprog[x, lora, lorb]]]*ColorDelta[c13, c23]*ColorDelta[c22, c24]*cqq[f24, f22]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f22, f24]*qgq[f11, f14]*SUNTF[{sun}, c11, c14]*SUNTF[{sun}, c12, c21] + atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[hv22 . xprop[-x] . v12]]*ColorDelta[c13, c23]*ColorDelta[c14, c24]*cqq[f24, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f12, f21]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c12, c21] + atr[fcc[xprog[x, lora, lorb] . hv21] . v11]*atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprop[-x]]*ColorDelta[c12, c21]*ColorDelta[c14, c24]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f23, f13]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c23, c13] + atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . hv22]*atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . xprop[x]]*ColorDelta[c12, c21]*ColorDelta[c14, c23]*cqq[f23, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f24, f13]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c24, c13] + atr[fcc[hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . xprop[-x]]*ColorDelta[c11, c22]*ColorDelta[c14, c23]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f12, f21]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c24, c13] + atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . hv22]*atr[xprop[x] . fcc[v11 . xprog[x, lora, lorb] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c14, c23]*cqq[f23, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f24, f13]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c24, c13] + atr[xprop[x] . fcc[v11 . hv21]]*atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . hv22 . xprog[-x, lora, lorb]]*ColorDelta[c11, c21]*ColorDelta[c12, c22]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f24, f13]*SUNTF[{sun}, c23, c14]*SUNTF[{sun}, c24, c13] + atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[hv21] . v11 . xprop[x]]*ColorDelta[c11, c22]*ColorDelta[c12, c21]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f24, f14]*SUNTF[{sun}, c23, c13]*SUNTF[{sun}, c24, c14];


dia9=(FourierXP[dia9,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia9=QEvaluate[I ScaleMu^(4-D) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type10[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia10,sun,lora,lorb},


dia10=atr[xprog[x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . hv21]]*atr[fcc[hv22 . xprop[-x] . v12]]*ColorDelta[c13, c23]*ColorDelta[c14, c24]*cqq[f24, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f12, f22]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c12, c22] - atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[hv22 . xprog[-x, lora, lorb] . v12] . hv21 . fcc[v11 . xprop[x]]]*ColorDelta[c12, c21]*ColorDelta[c22, c24]*cqq[f24, f22]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f22, f24]*qgq[f11, f14]*SUNTF[{sun}, c11, c14]*SUNTF[{sun}, c23, c13] + atr[fcc[xprop[x] . hv21] . v11]*atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c22]*ColorDelta[c12, c21]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f23, f14]*SUNTF[{sun}, c23, c14]*SUNTF[{sun}, c24, c13] + atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . hv21]]*atr[fcc[hv22 . xprop[-x] . v12] . xprog[-x, lora, lorb]]*ColorDelta[c12, c22]*ColorDelta[c13, c23]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f11, f21]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . v12] . xprog[-x, lora, lorb]]*atr[xprop[x] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c13, c23]*cqq[f23, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f12, f22]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprog[-x, lora, lorb]]*atr[xprop[x] . fcc[v11 . hv21]]*ColorDelta[c11, c21]*ColorDelta[c12, c22]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f23, f13]*SUNTF[{sun}, c23, c13]*SUNTF[{sun}, c24, c14] + atr[fcc[xprop[x] . hv21] . v11]*atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c22]*ColorDelta[c12, c21]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f24, f14]*SUNTF[{sun}, c23, c13]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[v11 . xprop[x] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c12, c22]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f24, f14]*SUNTF[{sun}, c23, c13]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*atr[xprop[x] . fcc[v11 . hv21]]*ColorDelta[c11, c21]*ColorDelta[c12, c22]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f24, f14]*SUNTF[{sun}, c23, c13]*SUNTF[{sun}, c24, c14] - atr[fcc[DiracSigma[GAD[lora], GAD[lorb]]] . hv22 . xprog[-x, lora, lorb] . fcc[v11 . v12] . xprop[x] . hv21]*ColorDelta[c11, c22]*ColorDelta[c12, c13]*cqq[f12, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f14, f23]*FlavorDelta[f21, f24]*qgq[f24, f21]*SUNTF[{sun}, c23, c14]*SUNTF[{sun}, c24, c21];


dia10=(FourierXP[dia10,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia10=QEvaluate[I ScaleMu^(4-D) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type11[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia11,sun,lora,lorb},


dia11=atr[fcc[hv22 . v12] . xprop[-x]]*atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . v11 . xprog[x, lora, lorb]]*ColorDelta[c13, c23]*ColorDelta[c14, c24]*cqq[f23, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f11, f22]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c12, c21] + atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . v11 . xprop[x]]*atr[fcc[hv22 . xprog[-x, lora, lorb] . v12]]*ColorDelta[c12, c21]*ColorDelta[c14, c24]*cqq[f24, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f11, f22]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c23, c13] + atr[xprog[x, lora, lorb] . fcc[v11 . hv21]]*atr[fcc[xprop[-x] . v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c12, c22]*ColorDelta[c13, c24]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f23, f14]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c23, c14] + atr[fcc[v12] . hv22 . xprog[-x, lora, lorb]]*atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . v11 . xprop[x]]*ColorDelta[c12, c21]*ColorDelta[c13, c24]*cqq[f24, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f11, f22]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c23, c14] - atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]]] . xprop[-x] . fcc[v11 . v12] . xprog[x, lora, lorb] . hv21]*ColorDelta[c12, c13]*ColorDelta[c14, c24]*cqq[f12, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f14, f24]*FlavorDelta[f21, f23]*qgq[f23, f21]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c23, c21] - atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . fcc[v11 . xprop[x] . hv21] . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c12, c22]*ColorDelta[c21, c23]*cqq[f23, f21]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f21, f23]*qgq[f11, f14]*SUNTF[{sun}, c11, c14]*SUNTF[{sun}, c24, c13] + atr[fcc[xprog[x, lora, lorb] . hv21] . v11]*atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . hv22 . xprop[-x]]*ColorDelta[c12, c21]*ColorDelta[c14, c23]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f24, f13]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c24, c13] + atr[fcc[hv21] . v11 . xprog[x, lora, lorb]]*atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . hv22 . xprop[-x]]*ColorDelta[c11, c22]*ColorDelta[c14, c23]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f24, f13]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c24, c13] - atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . hv22 . xprop[-x] . fcc[v11 . v12] . xprog[x, lora, lorb]]*ColorDelta[c12, c13]*ColorDelta[c14, c23]*cqq[f12, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f14, f23]*FlavorDelta[f22, f24]*qgq[f24, f22]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c24, c22] - atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora], GAD[lorb]] . hv21 . fcc[v11 . xprop[x]]]*ColorDelta[c11, c14]*ColorDelta[c12, c21]*cqq[f11, f14]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f22, f24]*qgq[f24, f22]*SUNTF[{sun}, c23, c13]*SUNTF[{sun}, c24, c22];


dia11=(FourierXP[dia11,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia11=QEvaluate[I ScaleMu^(4-D) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type12[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia12,sun,lora,lorb},


dia12=atr[fcc[hv22 . v12] . xprop[-x]]*atr[xprog[x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . hv21]]*ColorDelta[c13, c23]*ColorDelta[c14, c24]*cqq[f23, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f12, f22]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c12, c22] + atr[xprop[x] . fcc[v11 . xprog[x, lora, lorb] . hv21]]*atr[fcc[v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c21]*ColorDelta[c13, c24]*cqq[f24, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f23, f14]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c23, c14] - atr[fcc[xprop[-x]] . hv22 . DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21] . v12]*ColorDelta[c11, c13]*ColorDelta[c14, c24]*cqq[f11, f13]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f14, f24]*FlavorDelta[f21, f23]*qgq[f23, f21]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c23, c21] - atr[xprop[x] . fcc[hv22 . hv21] . xprog[-x, lora, lorb] . v12 . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]]]]*ColorDelta[c11, c21]*ColorDelta[c22, c23]*cqq[f23, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f13, f24]*FlavorDelta[f22, f23]*qgq[f12, f14]*SUNTF[{sun}, c12, c14]*SUNTF[{sun}, c24, c13] - atr[fcc[xprop[x] . hv21] . v11 . fcc[xprog[-x, lora, lorb] . v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c22]*ColorDelta[c12, c14]*cqq[f12, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f13, f24]*FlavorDelta[f21, f23]*qgq[f23, f21]*SUNTF[{sun}, c23, c21]*SUNTF[{sun}, c24, c13] - atr[fcc[xprog[-x, lora, lorb]] . hv22 . fcc[v11 . xprop[x] . hv21] . DiracSigma[GAD[lora], GAD[lorb]] . v12]*ColorDelta[c12, c22]*ColorDelta[c21, c23]*cqq[f23, f21]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f14, f24]*FlavorDelta[f21, f23]*qgq[f11, f13]*SUNTF[{sun}, c11, c13]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*atr[xprog[x, lora, lorb] . fcc[v11 . hv21]]*ColorDelta[c12, c22]*ColorDelta[c13, c23]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f24, f14]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[xprop[x] . hv21] . v11 . xprog[x, lora, lorb]]*ColorDelta[c11, c22]*ColorDelta[c13, c23]*cqq[f23, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f24, f14]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c24, c14] - atr[fcc[xprog[-x, lora, lorb]] . hv22 . DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprop[x] . hv21] . v12]*ColorDelta[c11, c13]*ColorDelta[c12, c22]*cqq[f11, f13]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f14, f24]*FlavorDelta[f21, f23]*qgq[f23, f21]*SUNTF[{sun}, c23, c21]*SUNTF[{sun}, c24, c14] - atr[fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21]]*ColorDelta[c11, c14]*ColorDelta[c13, c23]*cqq[f11, f14]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f21, f24]*qgq[f24, f21]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c24, c21];


dia12=(FourierXP[dia12,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia12=QEvaluate[I ScaleMu^(4-D) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type13[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia13,sun,lora,lorb},


dia13=atr[fcc[xprop[-x] . v12] . hv22]*atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c13, c24]*ColorDelta[c14, c23]*cqq[f23, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f12, f21]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c12, c21] + atr[xprop[x] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . hv21]]*atr[fcc[hv22 . xprog[-x, lora, lorb] . v12]]*ColorDelta[c11, c21]*ColorDelta[c14, c24]*cqq[f24, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f12, f22]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c23, c13] - atr[1 . v11 . xprog[x, lora, lorb] . fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . hv21] . xprop[-x] . v12]*ColorDelta[c11, c14]*ColorDelta[c13, c24]*cqq[f11, f14]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f22, f23]*qgq[f23, f22]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c23, c22] + atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . hv22]*atr[xprog[x, lora, lorb] . fcc[v11 . xprop[x] . hv21]]*ColorDelta[c12, c22]*ColorDelta[c14, c23]*cqq[f23, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f24, f13]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c24, c13] + atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . v11]*atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . xprop[-x]]*ColorDelta[c12, c21]*ColorDelta[c14, c23]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f11, f22]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c24, c13] - atr[xprop[x] . fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . hv21] . xprog[-x, lora, lorb] . v12 . fcc[v11]]*ColorDelta[c11, c21]*ColorDelta[c12, c14]*cqq[f12, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f13, f24]*FlavorDelta[f22, f23]*qgq[f23, f22]*SUNTF[{sun}, c23, c22]*SUNTF[{sun}, c24, c13] + atr[fcc[hv22 . xprop[-x] . v12] . xprog[-x, lora, lorb]]*atr[fcc[hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c11, c22]*ColorDelta[c13, c23]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f12, f21]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[v11 . xprog[x, lora, lorb] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c13, c23]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f24, f14]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c24, c14] - atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]]] . xprog[-x, lora, lorb] . fcc[v11 . v12] . xprop[x] . hv21]*ColorDelta[c11, c22]*ColorDelta[c12, c13]*cqq[f12, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f14, f24]*FlavorDelta[f21, f23]*qgq[f23, f21]*SUNTF[{sun}, c23, c21]*SUNTF[{sun}, c24, c14] - atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . hv21] . xprog[-x, lora, lorb] . fcc[v11 . v12] . xprop[x]]*ColorDelta[c11, c21]*ColorDelta[c12, c13]*cqq[f12, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f14, f24]*FlavorDelta[f22, f23]*qgq[f23, f22]*SUNTF[{sun}, c23, c22]*SUNTF[{sun}, c24, c14];


dia13=(FourierXP[dia13,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia13=QEvaluate[I ScaleMu^(4-D) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type14[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia14,sun,lora,lorb},


dia14=atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21]]*atr[fcc[v12] . hv22 . xprop[-x]]*ColorDelta[c13, c24]*ColorDelta[c14, c23]*cqq[f24, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f11, f21]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c12, c22] + atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21]]*atr[fcc[hv22 . v12] . xprop[-x]]*ColorDelta[c13, c23]*ColorDelta[c14, c24]*cqq[f23, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f11, f21]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c12, c22] + atr[fcc[xprop[-x] . v12] . hv22 . xprog[-x, lora, lorb]]*atr[fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c13, c24]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f12, f22]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c23, c14] + atr[fcc[xprop[-x] . v12] . hv22 . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[v11 . xprog[x, lora, lorb] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c13, c24]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f23, f14]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c23, c14] + atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . hv22]*atr[fcc[xprop[x] . hv21] . v11 . xprog[x, lora, lorb]]*ColorDelta[c11, c22]*ColorDelta[c14, c23]*cqq[f23, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f24, f13]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c24, c13] + atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . hv22 . xprog[-x, lora, lorb]]*atr[fcc[v11 . xprop[x] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c12, c22]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f24, f13]*SUNTF[{sun}, c23, c14]*SUNTF[{sun}, c24, c13] + atr[fcc[hv22 . v12] . xprog[-x, lora, lorb]]*atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . v11 . xprop[x]]*ColorDelta[c12, c21]*ColorDelta[c13, c23]*cqq[f23, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f11, f22]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[hv21] . v11 . xprog[x, lora, lorb]]*ColorDelta[c11, c22]*ColorDelta[c13, c23]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f24, f14]*SUNTF[{sun}, c12, c21]*SUNTF[{sun}, c24, c14] - atr[fcc[v12] . v11 . xprop[x] . fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . hv21] . xprog[-x, lora, lorb]]*ColorDelta[c11, c13]*ColorDelta[c12, c21]*cqq[f11, f13]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f14, f24]*FlavorDelta[f22, f23]*qgq[f23, f22]*SUNTF[{sun}, c23, c22]*SUNTF[{sun}, c24, c14] - atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . hv22 . xprog[-x, lora, lorb] . fcc[v11 . v12] . xprop[x]]*ColorDelta[c11, c21]*ColorDelta[c12, c13]*cqq[f12, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f14, f23]*FlavorDelta[f22, f24]*qgq[f24, f22]*SUNTF[{sun}, c23, c14]*SUNTF[{sun}, c24, c22];


dia14=(FourierXP[dia14,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia14=QEvaluate[I ScaleMu^(4-D) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type15[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia15,sun,lora,lorb},


dia15=-(atr[fcc[hv21] . hv22 . xprop[-x] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprog[x, lora, lorb]]*ColorDelta[c14, c23]*ColorDelta[c22, c24]*cqq[f24, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f14, f23]*FlavorDelta[f22, f24]*qgq[f12, f13]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c12, c13]) + atr[fcc[hv22 . v12] . xprop[-x]]*atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c13, c23]*ColorDelta[c14, c24]*cqq[f23, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f12, f21]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c12, c21] - atr[fcc[hv21] . hv22 . xprog[-x, lora, lorb] . fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . v11 . xprop[x]]*ColorDelta[c12, c21]*ColorDelta[c22, c24]*cqq[f24, f22]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f14, f23]*FlavorDelta[f22, f24]*qgq[f11, f13]*SUNTF[{sun}, c11, c13]*SUNTF[{sun}, c23, c14] - atr[1 . hv22 . xprog[-x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprop[x] . hv21]*ColorDelta[c11, c22]*ColorDelta[c21, c24]*cqq[f24, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f14, f23]*FlavorDelta[f21, f24]*qgq[f12, f13]*SUNTF[{sun}, c12, c13]*SUNTF[{sun}, c23, c14] - atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . hv21] . xprop[-x] . fcc[v11 . v12] . xprog[x, lora, lorb]]*ColorDelta[c12, c13]*ColorDelta[c14, c24]*cqq[f12, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f14, f24]*FlavorDelta[f22, f23]*qgq[f23, f22]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c23, c22] - atr[fcc[DiracSigma[GAD[lora], GAD[lorb]]] . v11 . xprop[x] . fcc[hv22 . hv21] . xprog[-x, lora, lorb] . v12]*ColorDelta[c12, c21]*ColorDelta[c22, c23]*cqq[f23, f22]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f22, f23]*qgq[f11, f14]*SUNTF[{sun}, c11, c14]*SUNTF[{sun}, c24, c13] + atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . xprop[-x]]*atr[fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c14, c23]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f12, f22]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c24, c13] - atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c12, c14]*ColorDelta[c13, c23]*cqq[f12, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f13, f23]*FlavorDelta[f21, f24]*qgq[f24, f21]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c24, c21] - atr[fcc[DiracSigma[GAD[lora], GAD[lorb]]] . hv22 . xprop[-x] . fcc[v11 . v12] . xprog[x, lora, lorb] . hv21]*ColorDelta[c12, c13]*ColorDelta[c14, c23]*cqq[f12, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f14, f23]*FlavorDelta[f21, f24]*qgq[f24, f21]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c24, c21] - atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprop[x] . hv21] . v12 . fcc[hv22 . xprog[-x, lora, lorb]]]*ColorDelta[c11, c13]*ColorDelta[c12, c22]*cqq[f11, f13]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f14, f23]*FlavorDelta[f21, f24]*qgq[f24, f21]*SUNTF[{sun}, c23, c14]*SUNTF[{sun}, c24, c21];


dia15=(FourierXP[dia15,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia15=QEvaluate[I ScaleMu^(4-D) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type16[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia16,sun,lora,lorb},


dia16=atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . hv21] . v11 . xprog[x, lora, lorb]]*atr[fcc[hv22 . xprop[-x] . v12]]*ColorDelta[c13, c23]*ColorDelta[c14, c24]*cqq[f24, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f11, f22]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c12, c21] + atr[fcc[v12] . hv22 . xprop[-x]]*atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c13, c24]*ColorDelta[c14, c23]*cqq[f24, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f12, f21]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c12, c21] - atr[fcc[hv22 . xprop[-x] . v12] . fcc[v11 . xprog[x, lora, lorb] . hv21] . DiracSigma[GAD[lora], GAD[lorb]]]*ColorDelta[c13, c23]*ColorDelta[c21, c24]*cqq[f24, f21]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f21, f24]*qgq[f11, f14]*SUNTF[{sun}, c11, c14]*SUNTF[{sun}, c12, c22] + atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . xprop[-x]]*atr[fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c14, c24]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f12, f22]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c23, c13] + atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . hv22 . xprop[-x]]*atr[fcc[v11 . xprog[x, lora, lorb] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c14, c23]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*qgq[f24, f13]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c24, c13] - atr[fcc[DiracSigma[GAD[lora], GAD[lorb]] . v12] . v11 . xprop[x] . fcc[hv22 . hv21] . xprog[-x, lora, lorb]]*ColorDelta[c12, c21]*ColorDelta[c22, c23]*cqq[f23, f22]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f14, f24]*FlavorDelta[f22, f23]*qgq[f11, f13]*SUNTF[{sun}, c11, c13]*SUNTF[{sun}, c24, c14] + atr[DiracSigma[GAD[lora], GAD[lorb]] . fcc[v11 . xprop[x] . hv21]]*atr[fcc[hv22 . v12] . xprog[-x, lora, lorb]]*ColorDelta[c12, c22]*ColorDelta[c13, c23]*cqq[f23, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f11, f21]*SUNTF[{sun}, c11, c21]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . v12] . DiracSigma[GAD[lora], GAD[lorb]]]*atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . xprop[x]]*ColorDelta[c12, c21]*ColorDelta[c13, c23]*cqq[f23, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f24, f14]*SUNTF[{sun}, c11, c22]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . xprop[-x] . v12] . xprog[-x, lora, lorb]]*atr[fcc[v11 . DiracSigma[GAD[lora], GAD[lorb]] . hv21]]*ColorDelta[c11, c21]*ColorDelta[c13, c23]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f12, f22]*SUNTF[{sun}, c12, c22]*SUNTF[{sun}, c24, c14] + atr[fcc[hv22 . DiracSigma[GAD[lora], GAD[lorb]] . v12] . xprog[-x, lora, lorb]]*atr[fcc[hv21] . v11 . xprop[x]]*ColorDelta[c11, c22]*ColorDelta[c12, c21]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*qgq[f23, f13]*SUNTF[{sun}, c23, c13]*SUNTF[{sun}, c24, c14];


dia16=(FourierXP[dia16,{x,q}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia16=QEvaluate[I ScaleMu^(4-D) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*-----------------------------------------------*)


(* ::Input::Initialization:: *)
(*type1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,sun,lora,lorb},


dia1=ColorDelta[c23, c13] ColorDelta[c24, c14] cqq[f23, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f11, f21] SUNTF[sun, c11, c21] SUNTF[sun, c12, c22] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21]] atr[fcc[hv22 . v12] . xprop[-x]] + ColorDelta[c11, c22] ColorDelta[c24, c14] cqq[f11, f22] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f23, f13] SUNTF[sun, c12, c21] SUNTF[sun, c23, c13] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . v12] . xprop[-x]] atr[fcc[hv21] . v11 . xprog[x, lora, lorb]] + ColorDelta[c12, c21] ColorDelta[c23, c14] cqq[f12, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f24, f13] SUNTF[sun, c11, c22] SUNTF[sun, c24, c13] atr[fcc[xprog[x, lora, lorb] . hv21] . v11] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . hv22 . xprop[-x]] + ColorDelta[c12, c21] ColorDelta[c23, c14] cqq[f12, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f11, f22] SUNTF[sun, c11, c22] SUNTF[sun, c24, c13] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . v11] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . xprop[-x]] + ColorDelta[c12, c21] ColorDelta[c24, c13] cqq[f24, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f23, f14] SUNTF[sun, c11, c22] SUNTF[sun, c23, c14] atr[fcc[v12] . hv22 . DiracSigma[GAD[lora, lorb]]] atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . xprop[x]] - ColorDelta[c11, c14] ColorDelta[c12, c22] cqq[f11, f14] FlavorDelta[f11, f14] FlavorDelta[f12, f22] FlavorDelta[f23, f21] FlavorDelta[f24, f13] qgq[f23, f21] SUNTF[sun, c23, c21] SUNTF[sun, c24, c13] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprop[x] . hv21]] - ColorDelta[c11, c21] ColorDelta[c23, c22] cqq[f23, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f14] FlavorDelta[f23, f22] FlavorDelta[f24, f13] qgq[f12, f14] SUNTF[sun, c12, c14] SUNTF[sun, c24, c13] atr[xprop[x] . fcc[hv22 . hv21] . xprog[-x, lora, lorb] . v12 . fcc[v11 . DiracSigma[GAD[lora, lorb]]]] - ColorDelta[c12, c13] ColorDelta[c23, c14] cqq[f12, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f13] FlavorDelta[f23, f14] FlavorDelta[f24, f21] qgq[f24, f21] SUNTF[sun, c11, c22] SUNTF[sun, c24, c21] atr[fcc[DiracSigma[GAD[lora, lorb]]] . hv22 . xprop[-x] . fcc[v11 . v12] . xprog[x, lora, lorb] . hv21] - ColorDelta[c12, c21] ColorDelta[c23, c22] cqq[f23, f22] FlavorDelta[f11, f14] FlavorDelta[f12, f21] FlavorDelta[f23, f22] FlavorDelta[f24, f13] qgq[f11, f14] SUNTF[sun, c11, c14] SUNTF[sun, c24, c13] atr[fcc[DiracSigma[GAD[lora, lorb]]] . v11 . xprop[x] . fcc[hv22 . hv21] . xprog[-x, lora, lorb] . v12] + ColorDelta[c11, c21] ColorDelta[c12, c22] cqq[f11, f21] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f23, f14] SUNTF[sun, c23, c14] SUNTF[sun, c24, c13] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . DiracSigma[GAD[lora, lorb]]] atr[fcc[v11 . xprop[x] . hv21]];


dia1=FourierXP[dia1,{x,q}]//SUNSimplify;

dia1=QEvaluate[I ScaleMu^(4-D) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,sun,lora,lorb},


dia2=ColorDelta[c23, c14] ColorDelta[c24, c13] cqq[f23, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f11, f21] SUNTF[sun, c11, c21] SUNTF[sun, c12, c22] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21]] atr[fcc[xprop[-x] . v12] . hv22] + ColorDelta[c12, c22] ColorDelta[c24, c14] cqq[f12, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f11, f21] SUNTF[sun, c11, c21] SUNTF[sun, c23, c13] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . hv21]] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . xprop[-x]] + ColorDelta[c23, c13] ColorDelta[c24, c14] cqq[f23, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f12, f22] SUNTF[sun, c11, c21] SUNTF[sun, c12, c22] atr[fcc[hv22 . v12] . xprop[-x]] atr[xprog[x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . hv21]] + ColorDelta[c11, c21] ColorDelta[c23, c14] cqq[f23, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f24, f13] SUNTF[sun, c12, c22] SUNTF[sun, c24, c13] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . hv22] atr[xprop[x] . fcc[v11 . xprog[x, lora, lorb] . hv21]] - ColorDelta[c11, c14] ColorDelta[c24, c13] cqq[f11, f14] FlavorDelta[f11, f14] FlavorDelta[f12, f22] FlavorDelta[f23, f21] FlavorDelta[f24, f13] qgq[f23, f21] SUNTF[sun, c12, c22] SUNTF[sun, c23, c21] atr[fcc[xprop[-x] . v12] . hv22 . DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21]] - ColorDelta[c11, c13] ColorDelta[c24, c14] cqq[f11, f13] FlavorDelta[f11, f13] FlavorDelta[f12, f21] FlavorDelta[f23, f22] FlavorDelta[f24, f14] qgq[f23, f22] SUNTF[sun, c12, c21] SUNTF[sun, c23, c22] atr[fcc[v12] . v11 . xprog[x, lora, lorb] . fcc[hv22 . DiracSigma[GAD[lora, lorb]] . hv21] . xprop[-x]] - ColorDelta[c11, c21] ColorDelta[c12, c13] cqq[f12, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f13] FlavorDelta[f23, f14] FlavorDelta[f24, f22] qgq[f24, f22] SUNTF[sun, c23, c14] SUNTF[sun, c24, c22] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . hv22 . xprog[-x, lora, lorb] . fcc[v11 . v12] . xprop[x]] - ColorDelta[c23, c22] ColorDelta[c24, c13] cqq[f23, f22] FlavorDelta[f11, f14] FlavorDelta[f12, f21] FlavorDelta[f23, f22] FlavorDelta[f24, f13] qgq[f11, f14] SUNTF[sun, c11, c14] SUNTF[sun, c12, c21] atr[fcc[DiracSigma[GAD[lora, lorb]]] . v11 . xprog[x, lora, lorb] . fcc[hv22 . hv21] . xprop[-x] . v12] + ColorDelta[c23, c13] ColorDelta[c24, c14] cqq[f24, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f12, f21] SUNTF[sun, c11, c22] SUNTF[sun, c12, c21] atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . DiracSigma[GAD[lora, lorb]]] atr[fcc[hv22 . xprop[-x] . v12]] + ColorDelta[c11, c21] ColorDelta[c23, c13] cqq[f11, f21] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f12, f22] SUNTF[sun, c12, c22] SUNTF[sun, c24, c14] atr[fcc[hv22 . xprop[-x] . v12] . xprog[-x, lora, lorb]] atr[fcc[v11 . DiracSigma[GAD[lora, lorb]] . hv21]];


dia2=FourierXP[dia2,{x,q}]//SUNSimplify;

dia2=QEvaluate[I ScaleMu^(4-D) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,sun,lora,lorb},


dia3=ColorDelta[c12, c22] ColorDelta[c23, c13] cqq[f12, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f11, f21] SUNTF[sun, c11, c21] SUNTF[sun, c24, c14] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . hv21]] atr[fcc[hv22 . xprop[-x] . v12] . xprog[-x, lora, lorb]] + ColorDelta[c12, c21] ColorDelta[c23, c13] cqq[f12, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f11, f22] SUNTF[sun, c11, c22] SUNTF[sun, c24, c14] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . v11] atr[fcc[hv22 . xprop[-x] . v12] . xprog[-x, lora, lorb]] + ColorDelta[c11, c21] ColorDelta[c23, c13] cqq[f23, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f12, f22] SUNTF[sun, c12, c22] SUNTF[sun, c24, c14] atr[fcc[hv22 . v12] . xprog[-x, lora, lorb]] atr[xprop[x] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . hv21]] + ColorDelta[c11, c22] ColorDelta[c24, c14] cqq[f11, f22] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f12, f21] SUNTF[sun, c12, c21] SUNTF[sun, c23, c13] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . xprop[-x]] atr[fcc[hv21] . v11 . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c12, c21] ColorDelta[c23, c14] cqq[f23, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f24, f13] SUNTF[sun, c11, c22] SUNTF[sun, c24, c13] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . hv22] atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . xprop[x]] + ColorDelta[c11, c22] ColorDelta[c24, c13] cqq[f24, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f23, f14] SUNTF[sun, c12, c21] SUNTF[sun, c23, c14] atr[fcc[v12] . hv22 . DiracSigma[GAD[lora, lorb]]] atr[fcc[xprop[x] . hv21] . v11 . xprog[x, lora, lorb]] - ColorDelta[c11, c21] ColorDelta[c23, c22] cqq[f23, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f13] FlavorDelta[f23, f22] FlavorDelta[f24, f14] qgq[f12, f13] SUNTF[sun, c12, c13] SUNTF[sun, c24, c14] atr[fcc[hv22 . hv21] . xprog[-x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . v12] . xprop[x]] - ColorDelta[c11, c22] ColorDelta[c23, c21] cqq[f23, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f14] FlavorDelta[f23, f21] FlavorDelta[f24, f13] qgq[f12, f14] SUNTF[sun, c12, c14] SUNTF[sun, c24, c13] atr[fcc[xprop[x] . hv21] . v11 . DiracSigma[GAD[lora, lorb]] . fcc[xprog[-x, lora, lorb] . v12] . hv22] - ColorDelta[c12, c14] ColorDelta[c23, c13] cqq[f12, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f14] FlavorDelta[f23, f13] FlavorDelta[f24, f22] qgq[f24, f22] SUNTF[sun, c11, c21] SUNTF[sun, c24, c22] atr[fcc[xprog[x, lora, lorb]] . v11 . fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora, lorb]] . hv21] + ColorDelta[c11, c22] ColorDelta[c24, c14] cqq[f24, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f12, f21] SUNTF[sun, c12, c21] SUNTF[sun, c23, c13] atr[fcc[xprop[x] . hv21] . v11 . DiracSigma[GAD[lora, lorb]]] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12]];


dia3=FourierXP[dia3,{x,q}]//SUNSimplify;

dia3=QEvaluate[I ScaleMu^(4-D) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,sun,lora,lorb},


dia4=ColorDelta[c12, c21] ColorDelta[c23, c13] cqq[f12, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f24, f14] SUNTF[sun, c11, c22] SUNTF[sun, c24, c14] atr[fcc[xprog[x, lora, lorb] . hv21] . v11] atr[fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c12, c22] ColorDelta[c23, c13] cqq[f12, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f24, f14] SUNTF[sun, c11, c21] SUNTF[sun, c24, c14] atr[fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora, lorb]]] atr[xprog[x, lora, lorb] . fcc[v11 . hv21]] + ColorDelta[c11, c22] ColorDelta[c23, c13] cqq[f11, f22] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f24, f14] SUNTF[sun, c12, c21] SUNTF[sun, c24, c14] atr[fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora, lorb]]] atr[fcc[hv21] . v11 . xprog[x, lora, lorb]] + ColorDelta[c11, c22] ColorDelta[c24, c13] cqq[f11, f22] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f12, f21] SUNTF[sun, c12, c21] SUNTF[sun, c23, c14] atr[fcc[hv21] . v11 . DiracSigma[GAD[lora, lorb]]] atr[fcc[xprop[-x] . v12] . hv22 . xprog[-x, lora, lorb]] - ColorDelta[c12, c14] ColorDelta[c23, c13] cqq[f12, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f14] FlavorDelta[f23, f13] FlavorDelta[f24, f21] qgq[f24, f21] SUNTF[sun, c11, c22] SUNTF[sun, c24, c21] atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c11, c22] ColorDelta[c23, c21] cqq[f23, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f13] FlavorDelta[f23, f21] FlavorDelta[f24, f14] qgq[f12, f13] SUNTF[sun, c12, c13] SUNTF[sun, c24, c14] atr[fcc[hv22] . xprog[-x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . v12] . xprop[x] . hv21] - ColorDelta[c11, c13] ColorDelta[c12, c21] cqq[f11, f13] FlavorDelta[f11, f13] FlavorDelta[f12, f21] FlavorDelta[f23, f22] FlavorDelta[f24, f14] qgq[f23, f22] SUNTF[sun, c23, c22] SUNTF[sun, c24, c14] atr[fcc[v12] . v11 . xprop[x] . fcc[hv22 . DiracSigma[GAD[lora, lorb]] . hv21] . xprog[-x, lora, lorb]] - ColorDelta[c12, c21] ColorDelta[c23, c22] cqq[f23, f22] FlavorDelta[f11, f13] FlavorDelta[f12, f21] FlavorDelta[f23, f22] FlavorDelta[f24, f14] qgq[f11, f13] SUNTF[sun, c11, c13] SUNTF[sun, c24, c14] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . v11 . xprop[x] . fcc[hv22 . hv21] . xprog[-x, lora, lorb]] - ColorDelta[c11, c21] ColorDelta[c24, c22] cqq[f24, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f14] FlavorDelta[f23, f13] FlavorDelta[f24, f22] qgq[f12, f14] SUNTF[sun, c12, c14] SUNTF[sun, c23, c13] atr[fcc[xprop[x]] . v11 . DiracSigma[GAD[lora, lorb]] . fcc[hv22 . xprog[-x, lora, lorb] . v12] . hv21] - ColorDelta[c23, c22] ColorDelta[c24, c13] cqq[f23, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f14] FlavorDelta[f23, f22] FlavorDelta[f24, f13] qgq[f12, f14] SUNTF[sun, c11, c21] SUNTF[sun, c12, c14] atr[xprog[x, lora, lorb] . fcc[hv22 . hv21] . xprop[-x] . v12 . fcc[v11 . DiracSigma[GAD[lora, lorb]]]];


dia4=FourierXP[dia4,{x,q}]//SUNSimplify;

dia4=QEvaluate[I ScaleMu^(4-D) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia5,sun,lora,lorb},


dia5=ColorDelta[c23, c14] ColorDelta[c24, c13] cqq[f23, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f12, f22] SUNTF[sun, c11, c21] SUNTF[sun, c12, c22] atr[fcc[xprop[-x] . v12] . hv22] atr[xprog[x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . hv21]] + ColorDelta[c11, c21] ColorDelta[c23, c14] cqq[f23, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f12, f22] SUNTF[sun, c12, c22] SUNTF[sun, c24, c13] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22] atr[xprop[x] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . hv21]] + ColorDelta[c11, c21] ColorDelta[c24, c13] cqq[f24, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f12, f22] SUNTF[sun, c12, c22] SUNTF[sun, c23, c14] atr[xprop[x] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . hv21]] atr[fcc[v12] . hv22 . xprog[-x, lora, lorb]] + ColorDelta[c23, c14] ColorDelta[c24, c13] cqq[f24, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f12, f22] SUNTF[sun, c11, c21] SUNTF[sun, c12, c22] atr[xprog[x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . hv21]] atr[fcc[v12] . hv22 . xprop[-x]] + ColorDelta[c23, c14] ColorDelta[c24, c13] cqq[f24, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f11, f22] SUNTF[sun, c11, c22] SUNTF[sun, c12, c21] atr[fcc[v12] . hv22 . xprop[-x]] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . v11 . xprog[x, lora, lorb]] + ColorDelta[c11, c22] ColorDelta[c12, c21] cqq[f12, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f23, f14] SUNTF[sun, c23, c14] SUNTF[sun, c24, c13] atr[fcc[xprop[x] . hv21] . v11] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c12, c22] ColorDelta[c24, c21] cqq[f24, f21] FlavorDelta[f11, f14] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f21] qgq[f11, f14] SUNTF[sun, c11, c14] SUNTF[sun, c23, c13] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . fcc[v11 . xprop[x] . hv21] . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c11, c13] ColorDelta[c23, c14] cqq[f11, f13] FlavorDelta[f11, f13] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f21] qgq[f24, f21] SUNTF[sun, c12, c22] SUNTF[sun, c24, c21] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21] . v12 . fcc[hv22 . xprop[-x]]] - ColorDelta[c12, c22] ColorDelta[c23, c21] cqq[f23, f21] FlavorDelta[f11, f14] FlavorDelta[f12, f22] FlavorDelta[f23, f21] FlavorDelta[f24, f13] qgq[f11, f14] SUNTF[sun, c11, c14] SUNTF[sun, c24, c13] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . fcc[v11 . xprop[x] . hv21] . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c23, c13] ColorDelta[c24, c22] cqq[f24, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f14] FlavorDelta[f23, f13] FlavorDelta[f24, f22] qgq[f12, f14] SUNTF[sun, c11, c21] SUNTF[sun, c12, c14] atr[fcc[xprog[x, lora, lorb]] . v11 . DiracSigma[GAD[lora, lorb]] . fcc[hv22 . xprop[-x] . v12] . hv21];


dia5=FourierXP[dia5,{x,q}]//SUNSimplify;

dia5=QEvaluate[I ScaleMu^(4-D) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia6,sun,lora,lorb},


dia6=ColorDelta[c11, c21] ColorDelta[c24, c13] cqq[f24, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f23, f14] SUNTF[sun, c12, c22] SUNTF[sun, c23, c14] atr[xprop[x] . fcc[v11 . xprog[x, lora, lorb] . hv21]] atr[fcc[v12] . hv22 . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c11, c22] ColorDelta[c12, c21] cqq[f11, f22] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f24, f13] SUNTF[sun, c23, c14] SUNTF[sun, c24, c13] atr[fcc[hv21] . v11 . xprop[x]] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . hv22 . xprog[-x, lora, lorb]] + ColorDelta[c11, c22] ColorDelta[c23, c14] cqq[f11, f22] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f24, f13] SUNTF[sun, c12, c21] SUNTF[sun, c24, c13] atr[fcc[hv21] . v11 . xprog[x, lora, lorb]] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . hv22 . xprop[-x]] - ColorDelta[c23, c13] ColorDelta[c24, c21] cqq[f24, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f14] FlavorDelta[f23, f13] FlavorDelta[f24, f21] qgq[f12, f14] SUNTF[sun, c11, c22] SUNTF[sun, c12, c14] atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . DiracSigma[GAD[lora, lorb]] . fcc[hv22 . xprop[-x] . v12]] - ColorDelta[c11, c14] ColorDelta[c23, c13] cqq[f11, f14] FlavorDelta[f11, f14] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f22] qgq[f24, f22] SUNTF[sun, c12, c21] SUNTF[sun, c24, c22] atr[fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora, lorb]] . hv21 . fcc[v11 . xprog[x, lora, lorb]]] - ColorDelta[c12, c22] ColorDelta[c24, c21] cqq[f24, f21] FlavorDelta[f11, f13] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f21] qgq[f11, f13] SUNTF[sun, c11, c13] SUNTF[sun, c23, c14] atr[fcc[v11 . xprop[x] . hv21] . DiracSigma[GAD[lora, lorb]] . v12 . fcc[hv22 . xprog[-x, lora, lorb]]] - ColorDelta[c12, c14] ColorDelta[c24, c13] cqq[f12, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f14] FlavorDelta[f23, f21] FlavorDelta[f24, f13] qgq[f23, f21] SUNTF[sun, c11, c22] SUNTF[sun, c23, c21] atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . fcc[xprop[-x] . v12] . hv22 . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c11, c22] ColorDelta[c12, c14] cqq[f12, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f14] FlavorDelta[f23, f21] FlavorDelta[f24, f13] qgq[f23, f21] SUNTF[sun, c23, c21] SUNTF[sun, c24, c13] atr[fcc[xprop[x] . hv21] . v11 . fcc[xprog[-x, lora, lorb] . v12] . hv22 . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c23, c14] ColorDelta[c24, c21] cqq[f24, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f13] FlavorDelta[f23, f14] FlavorDelta[f24, f21] qgq[f12, f13] SUNTF[sun, c11, c22] SUNTF[sun, c12, c13] atr[1 . hv22 . xprop[-x] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . v12] . xprog[x, lora, lorb] . hv21] + ColorDelta[c23, c13] ColorDelta[c24, c14] cqq[f24, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f11, f22] SUNTF[sun, c11, c22] SUNTF[sun, c12, c21] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . v11 . xprog[x, lora, lorb]] atr[fcc[hv22 . xprop[-x] . v12]];


dia6=FourierXP[dia6,{x,q}]//SUNSimplify;

dia6=QEvaluate[I ScaleMu^(4-D) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia7,sun,lora,lorb},


dia7=ColorDelta[c12, c21] ColorDelta[c24, c14] cqq[f12, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f23, f13] SUNTF[sun, c11, c22] SUNTF[sun, c23, c13] atr[fcc[xprog[x, lora, lorb] . hv21] . v11] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . v12] . xprop[-x]] + ColorDelta[c12, c21] ColorDelta[c24, c13] cqq[f24, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f11, f22] SUNTF[sun, c11, c22] SUNTF[sun, c23, c14] atr[fcc[v12] . hv22 . xprog[-x, lora, lorb]] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . v11 . xprop[x]] + ColorDelta[c11, c22] ColorDelta[c12, c21] cqq[f11, f22] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f23, f14] SUNTF[sun, c23, c14] SUNTF[sun, c24, c13] atr[fcc[hv21] . v11 . xprop[x]] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c12, c21] ColorDelta[c24, c13] cqq[f12, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f23, f14] SUNTF[sun, c11, c22] SUNTF[sun, c23, c14] atr[fcc[xprog[x, lora, lorb] . hv21] . v11] atr[fcc[xprop[-x] . v12] . hv22 . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c12, c22] ColorDelta[c24, c13] cqq[f12, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f23, f14] SUNTF[sun, c11, c21] SUNTF[sun, c23, c14] atr[xprog[x, lora, lorb] . fcc[v11 . hv21]] atr[fcc[xprop[-x] . v12] . hv22 . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c11, c22] ColorDelta[c24, c13] cqq[f11, f22] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f23, f14] SUNTF[sun, c12, c21] SUNTF[sun, c23, c14] atr[fcc[hv21] . v11 . xprog[x, lora, lorb]] atr[fcc[xprop[-x] . v12] . hv22 . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c23, c14] ColorDelta[c24, c22] cqq[f24, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f13] FlavorDelta[f23, f14] FlavorDelta[f24, f22] qgq[f12, f13] SUNTF[sun, c11, c21] SUNTF[sun, c12, c13] atr[fcc[hv21] . hv22 . xprop[-x] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . v12] . xprog[x, lora, lorb]] - ColorDelta[c23, c21] ColorDelta[c24, c14] cqq[f23, f21] FlavorDelta[f11, f13] FlavorDelta[f12, f22] FlavorDelta[f23, f21] FlavorDelta[f24, f14] qgq[f11, f13] SUNTF[sun, c11, c13] SUNTF[sun, c12, c22] atr[fcc[xprop[-x]] . hv22 . fcc[v11 . xprog[x, lora, lorb] . hv21] . DiracSigma[GAD[lora, lorb]] . v12] + ColorDelta[c11, c22] ColorDelta[c24, c14] cqq[f24, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f23, f13] SUNTF[sun, c12, c21] SUNTF[sun, c23, c13] atr[fcc[xprop[x] . hv21] . v11 . xprog[x, lora, lorb]] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . v12]] + ColorDelta[c11, c21] ColorDelta[c12, c22] cqq[f11, f21] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f23, f13] SUNTF[sun, c23, c13] SUNTF[sun, c24, c14] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . v12] . xprog[-x, lora, lorb]] atr[fcc[v11 . xprop[x] . hv21]];


dia7=FourierXP[dia7,{x,q}]//SUNSimplify;

dia7=QEvaluate[I ScaleMu^(4-D) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia8,sun,lora,lorb},


dia8=ColorDelta[c12, c22] ColorDelta[c23, c13] cqq[f23, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f24, f14] SUNTF[sun, c11, c21] SUNTF[sun, c24, c14] atr[fcc[hv22 . v12] . DiracSigma[GAD[lora, lorb]]] atr[xprog[x, lora, lorb] . fcc[v11 . xprop[x] . hv21]] + ColorDelta[c11, c21] ColorDelta[c23, c13] cqq[f23, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f24, f14] SUNTF[sun, c12, c22] SUNTF[sun, c24, c14] atr[fcc[hv22 . v12] . DiracSigma[GAD[lora, lorb]]] atr[xprop[x] . fcc[v11 . xprog[x, lora, lorb] . hv21]] + ColorDelta[c11, c22] ColorDelta[c23, c14] cqq[f11, f22] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f12, f21] SUNTF[sun, c12, c21] SUNTF[sun, c24, c13] atr[fcc[hv21] . v11 . DiracSigma[GAD[lora, lorb]]] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . xprop[-x]] + ColorDelta[c23, c13] ColorDelta[c24, c14] cqq[f23, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f12, f21] SUNTF[sun, c11, c22] SUNTF[sun, c12, c21] atr[fcc[hv22 . v12] . xprop[-x]] atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c11, c22] ColorDelta[c23, c14] cqq[f23, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f12, f21] SUNTF[sun, c12, c21] SUNTF[sun, c24, c13] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22] atr[fcc[xprop[x] . hv21] . v11 . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c23, c13] ColorDelta[c24, c21] cqq[f24, f21] FlavorDelta[f11, f14] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f21] qgq[f11, f14] SUNTF[sun, c11, c14] SUNTF[sun, c12, c22] atr[fcc[hv22 . xprop[-x] . v12] . fcc[v11 . xprog[x, lora, lorb] . hv21] . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c12, c22] ColorDelta[c23, c21] cqq[f23, f21] FlavorDelta[f11, f13] FlavorDelta[f12, f22] FlavorDelta[f23, f21] FlavorDelta[f24, f14] qgq[f11, f13] SUNTF[sun, c11, c13] SUNTF[sun, c24, c14] atr[fcc[xprog[-x, lora, lorb]] . hv22 . fcc[v11 . xprop[x] . hv21] . DiracSigma[GAD[lora, lorb]] . v12] - ColorDelta[c23, c14] ColorDelta[c24, c22] cqq[f24, f22] FlavorDelta[f11, f13] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f22] qgq[f11, f13] SUNTF[sun, c11, c13] SUNTF[sun, c12, c21] atr[fcc[hv21] . hv22 . xprop[-x] . fcc[DiracSigma[GAD[lora, lorb]] . v12] . v11 . xprog[x, lora, lorb]] + ColorDelta[c12, c22] ColorDelta[c24, c14] cqq[f24, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f11, f21] SUNTF[sun, c11, c21] SUNTF[sun, c23, c13] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprop[x] . hv21]] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12]] + ColorDelta[c23, c13] ColorDelta[c24, c14] cqq[f24, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f12, f22] SUNTF[sun, c11, c21] SUNTF[sun, c12, c22] atr[xprog[x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . hv21]] atr[fcc[hv22 . xprop[-x] . v12]];


dia8=FourierXP[dia8,{x,q}]//SUNSimplify;

dia8=QEvaluate[I ScaleMu^(4-D) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type9[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia9,sun,lora,lorb},


dia9=ColorDelta[c12, c22] ColorDelta[c23, c13] cqq[f23, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f11, f21] SUNTF[sun, c11, c21] SUNTF[sun, c24, c14] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprop[x] . hv21]] atr[fcc[hv22 . v12] . xprog[-x, lora, lorb]] + ColorDelta[c11, c22] ColorDelta[c12, c21] cqq[f12, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f24, f14] SUNTF[sun, c23, c13] SUNTF[sun, c24, c14] atr[fcc[xprop[x] . hv21] . v11] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c11, c22] ColorDelta[c12, c21] cqq[f11, f22] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f23, f13] SUNTF[sun, c23, c13] SUNTF[sun, c24, c14] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . v12] . xprog[-x, lora, lorb]] atr[fcc[hv21] . v11 . xprop[x]] + ColorDelta[c11, c22] ColorDelta[c23, c13] cqq[f23, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f12, f21] SUNTF[sun, c12, c21] SUNTF[sun, c24, c14] atr[fcc[hv22 . v12] . xprog[-x, lora, lorb]] atr[fcc[xprop[x] . hv21] . v11 . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c11, c22] ColorDelta[c23, c14] cqq[f23, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f24, f13] SUNTF[sun, c12, c21] SUNTF[sun, c24, c13] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . hv22] atr[fcc[xprop[x] . hv21] . v11 . xprog[x, lora, lorb]] - ColorDelta[c11, c14] ColorDelta[c12, c22] cqq[f11, f14] FlavorDelta[f11, f14] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f21] qgq[f24, f21] SUNTF[sun, c23, c13] SUNTF[sun, c24, c21] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprop[x] . hv21]] - ColorDelta[c23, c13] ColorDelta[c24, c22] cqq[f24, f22] FlavorDelta[f11, f14] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f22] qgq[f11, f14] SUNTF[sun, c11, c14] SUNTF[sun, c12, c21] atr[DiracSigma[GAD[lora, lorb]] . fcc[hv22 . xprop[-x] . v12] . hv21 . fcc[v11 . xprog[x, lora, lorb]]] - ColorDelta[c23, c22] ColorDelta[c24, c14] cqq[f23, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f13] FlavorDelta[f23, f22] FlavorDelta[f24, f14] qgq[f12, f13] SUNTF[sun, c11, c21] SUNTF[sun, c12, c13] atr[fcc[hv22 . hv21] . xprop[-x] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . v12] . xprog[x, lora, lorb]] + ColorDelta[c23, c13] ColorDelta[c24, c14] cqq[f24, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f11, f21] SUNTF[sun, c11, c21] SUNTF[sun, c12, c22] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21]] atr[fcc[hv22 . xprop[-x] . v12]] + ColorDelta[c11, c21] ColorDelta[c12, c22] cqq[f11, f21] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f24, f14] SUNTF[sun, c23, c13] SUNTF[sun, c24, c14] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora, lorb]]] atr[fcc[v11 . xprop[x] . hv21]];


dia9=FourierXP[dia9,{x,q}]//SUNSimplify;

dia9=QEvaluate[I ScaleMu^(4-D) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type10[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia10,sun,lora,lorb},


dia10=ColorDelta[c11, c22] ColorDelta[c12, c21] cqq[f12, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f23, f13] SUNTF[sun, c23, c13] SUNTF[sun, c24, c14] atr[fcc[xprop[x] . hv21] . v11] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . v12] . xprog[-x, lora, lorb]] + ColorDelta[c11, c22] ColorDelta[c12, c21] cqq[f12, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f24, f13] SUNTF[sun, c23, c14] SUNTF[sun, c24, c13] atr[fcc[xprop[x] . hv21] . v11] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . hv22 . xprog[-x, lora, lorb]] + ColorDelta[c11, c21] ColorDelta[c12, c22] cqq[f12, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f24, f13] SUNTF[sun, c23, c14] SUNTF[sun, c24, c13] atr[xprop[x] . fcc[v11 . hv21]] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . hv22 . xprog[-x, lora, lorb]] - ColorDelta[c11, c14] ColorDelta[c23, c13] cqq[f11, f14] FlavorDelta[f11, f14] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f21] qgq[f24, f21] SUNTF[sun, c12, c22] SUNTF[sun, c24, c21] atr[fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21]] - ColorDelta[c11, c22] ColorDelta[c12, c13] cqq[f12, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f13] FlavorDelta[f23, f21] FlavorDelta[f24, f14] qgq[f23, f21] SUNTF[sun, c23, c21] SUNTF[sun, c24, c14] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]]] . xprog[-x, lora, lorb] . fcc[v11 . v12] . xprop[x] . hv21] - ColorDelta[c12, c13] ColorDelta[c23, c14] cqq[f12, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f13] FlavorDelta[f23, f14] FlavorDelta[f24, f22] qgq[f24, f22] SUNTF[sun, c11, c21] SUNTF[sun, c24, c22] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . hv22 . xprop[-x] . fcc[v11 . v12] . xprog[x, lora, lorb]] - ColorDelta[c11, c22] ColorDelta[c24, c21] cqq[f24, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f13] FlavorDelta[f23, f14] FlavorDelta[f24, f21] qgq[f12, f13] SUNTF[sun, c12, c13] SUNTF[sun, c23, c14] atr[1 . hv22 . xprog[-x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . v12] . xprop[x] . hv21] - ColorDelta[c11, c14] ColorDelta[c12, c21] cqq[f11, f14] FlavorDelta[f11, f14] FlavorDelta[f12, f21] FlavorDelta[f23, f22] FlavorDelta[f24, f13] qgq[f23, f22] SUNTF[sun, c23, c22] SUNTF[sun, c24, c13] atr[1 . v11 . xprop[x] . fcc[hv22 . DiracSigma[GAD[lora, lorb]] . hv21] . xprog[-x, lora, lorb] . v12] - ColorDelta[c11, c13] ColorDelta[c12, c21] cqq[f11, f13] FlavorDelta[f11, f13] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f22] qgq[f24, f22] SUNTF[sun, c23, c14] SUNTF[sun, c24, c22] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . hv22 . xprog[-x, lora, lorb] . fcc[v12] . v11 . xprop[x]] + ColorDelta[c11, c21] ColorDelta[c24, c14] cqq[f24, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f12, f22] SUNTF[sun, c12, c22] SUNTF[sun, c23, c13] atr[xprop[x] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . hv21]] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12]];


dia10=FourierXP[dia10,{x,q}]//SUNSimplify;

dia10=QEvaluate[I ScaleMu^(4-D) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type11[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia11,sun,lora,lorb},


dia11=ColorDelta[c12, c22] ColorDelta[c23, c14] cqq[f23, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f24, f13] SUNTF[sun, c11, c21] SUNTF[sun, c24, c13] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . hv22] atr[xprog[x, lora, lorb] . fcc[v11 . xprop[x] . hv21]] + ColorDelta[c11, c22] ColorDelta[c12, c21] cqq[f11, f22] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f24, f14] SUNTF[sun, c23, c13] SUNTF[sun, c24, c14] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora, lorb]]] atr[fcc[hv21] . v11 . xprop[x]] + ColorDelta[c23, c14] ColorDelta[c24, c13] cqq[f24, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f11, f21] SUNTF[sun, c11, c21] SUNTF[sun, c12, c22] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21]] atr[fcc[v12] . hv22 . xprop[-x]] + ColorDelta[c12, c22] ColorDelta[c24, c13] cqq[f12, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f11, f21] SUNTF[sun, c11, c21] SUNTF[sun, c23, c14] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . hv21]] atr[fcc[xprop[-x] . v12] . hv22 . xprog[-x, lora, lorb]] - ColorDelta[c11, c13] ColorDelta[c24, c14] cqq[f11, f13] FlavorDelta[f11, f13] FlavorDelta[f12, f22] FlavorDelta[f23, f21] FlavorDelta[f24, f14] qgq[f23, f21] SUNTF[sun, c12, c22] SUNTF[sun, c23, c21] atr[fcc[xprop[-x]] . hv22 . DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprog[x, lora, lorb] . hv21] . v12] + ColorDelta[c12, c21] ColorDelta[c24, c14] cqq[f24, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f23, f13] SUNTF[sun, c11, c22] SUNTF[sun, c23, c13] atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . xprop[x]] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . v12]] + ColorDelta[c12, c21] ColorDelta[c24, c14] cqq[f24, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f11, f22] SUNTF[sun, c11, c22] SUNTF[sun, c23, c13] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . v11 . xprop[x]] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12]] + ColorDelta[c11, c21] ColorDelta[c24, c14] cqq[f11, f21] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f12, f22] SUNTF[sun, c12, c22] SUNTF[sun, c23, c13] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . xprop[-x]] atr[fcc[v11 . DiracSigma[GAD[lora, lorb]] . hv21]] + ColorDelta[c11, c21] ColorDelta[c24, c13] cqq[f11, f21] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f23, f14] SUNTF[sun, c12, c22] SUNTF[sun, c23, c14] atr[fcc[xprop[-x] . v12] . hv22 . DiracSigma[GAD[lora, lorb]]] atr[fcc[v11 . xprog[x, lora, lorb] . hv21]] + ColorDelta[c11, c21] ColorDelta[c12, c22] cqq[f11, f21] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f24, f13] SUNTF[sun, c23, c14] SUNTF[sun, c24, c13] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . hv22 . xprog[-x, lora, lorb]] atr[fcc[v11 . xprop[x] . hv21]];


dia11=FourierXP[dia11,{x,q}]//SUNSimplify;

dia11=QEvaluate[I ScaleMu^(4-D) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type12[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia12,sun,lora,lorb},


dia12=ColorDelta[c12, c21] ColorDelta[c23, c13] cqq[f23, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f11, f22] SUNTF[sun, c11, c22] SUNTF[sun, c24, c14] atr[fcc[hv22 . v12] . xprog[-x, lora, lorb]] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . v11 . xprop[x]] + ColorDelta[c23, c14] ColorDelta[c24, c13] cqq[f23, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f12, f21] SUNTF[sun, c11, c22] SUNTF[sun, c12, c21] atr[fcc[xprop[-x] . v12] . hv22] atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c11, c22] ColorDelta[c24, c13] cqq[f24, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f12, f21] SUNTF[sun, c12, c21] SUNTF[sun, c23, c14] atr[fcc[v12] . hv22 . xprog[-x, lora, lorb]] atr[fcc[xprop[x] . hv21] . v11 . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c11, c22] ColorDelta[c12, c14] cqq[f12, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f14] FlavorDelta[f23, f13] FlavorDelta[f24, f21] qgq[f24, f21] SUNTF[sun, c23, c13] SUNTF[sun, c24, c21] atr[fcc[xprop[x] . hv21] . v11 . fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c12, c13] ColorDelta[c24, c14] cqq[f12, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f13] FlavorDelta[f23, f22] FlavorDelta[f24, f14] qgq[f23, f22] SUNTF[sun, c11, c21] SUNTF[sun, c23, c22] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . hv21] . xprop[-x] . fcc[v11 . v12] . xprog[x, lora, lorb]] - ColorDelta[c11, c21] ColorDelta[c24, c22] cqq[f24, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f13] FlavorDelta[f23, f14] FlavorDelta[f24, f22] qgq[f12, f13] SUNTF[sun, c12, c13] SUNTF[sun, c23, c14] atr[fcc[hv21] . hv22 . xprog[-x, lora, lorb] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . v12] . xprop[x]] - ColorDelta[c12, c13] ColorDelta[c24, c14] cqq[f12, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f13] FlavorDelta[f23, f21] FlavorDelta[f24, f14] qgq[f23, f21] SUNTF[sun, c11, c22] SUNTF[sun, c23, c21] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]]] . xprop[-x] . fcc[v11 . v12] . xprog[x, lora, lorb] . hv21] - ColorDelta[c11, c21] ColorDelta[c12, c14] cqq[f12, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f14] FlavorDelta[f23, f13] FlavorDelta[f24, f22] qgq[f24, f22] SUNTF[sun, c23, c13] SUNTF[sun, c24, c22] atr[fcc[xprop[x]] . v11 . fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora, lorb]] . hv21] + ColorDelta[c11, c21] ColorDelta[c24, c13] cqq[f11, f21] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f12, f22] SUNTF[sun, c12, c22] SUNTF[sun, c23, c14] atr[fcc[xprop[-x] . v12] . hv22 . xprog[-x, lora, lorb]] atr[fcc[v11 . DiracSigma[GAD[lora, lorb]] . hv21]] + ColorDelta[c11, c21] ColorDelta[c24, c14] cqq[f11, f21] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f23, f13] SUNTF[sun, c12, c22] SUNTF[sun, c23, c13] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . v12] . xprop[-x]] atr[fcc[v11 . xprog[x, lora, lorb] . hv21]];


dia12=FourierXP[dia12,{x,q}]//SUNSimplify;

dia12=QEvaluate[I ScaleMu^(4-D) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type13[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia13,sun,lora,lorb},


dia13=ColorDelta[c12, c21] ColorDelta[c24, c14] cqq[f12, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f11, f22] SUNTF[sun, c11, c22] SUNTF[sun, c23, c13] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . v11] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . xprop[-x]] + ColorDelta[c11, c22] ColorDelta[c23, c13] cqq[f11, f22] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f12, f21] SUNTF[sun, c12, c21] SUNTF[sun, c24, c14] atr[fcc[hv22 . xprop[-x] . v12] . xprog[-x, lora, lorb]] atr[fcc[hv21] . v11 . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c23, c14] ColorDelta[c24, c13] cqq[f23, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f11, f22] SUNTF[sun, c11, c22] SUNTF[sun, c12, c21] atr[fcc[xprop[-x] . v12] . hv22] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . v11 . xprog[x, lora, lorb]] + ColorDelta[c12, c22] ColorDelta[c23, c14] cqq[f12, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f24, f13] SUNTF[sun, c11, c21] SUNTF[sun, c24, c13] atr[xprog[x, lora, lorb] . fcc[v11 . hv21]] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . hv22 . xprop[-x]] + ColorDelta[c11, c22] ColorDelta[c23, c13] cqq[f23, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f24, f14] SUNTF[sun, c12, c21] SUNTF[sun, c24, c14] atr[fcc[hv22 . v12] . DiracSigma[GAD[lora, lorb]]] atr[fcc[xprop[x] . hv21] . v11 . xprog[x, lora, lorb]] - ColorDelta[c12, c21] ColorDelta[c24, c22] cqq[f24, f22] FlavorDelta[f11, f14] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f22] qgq[f11, f14] SUNTF[sun, c11, c14] SUNTF[sun, c23, c13] atr[DiracSigma[GAD[lora, lorb]] . fcc[hv22 . xprog[-x, lora, lorb] . v12] . hv21 . fcc[v11 . xprop[x]]] - ColorDelta[c11, c13] ColorDelta[c12, c22] cqq[f11, f13] FlavorDelta[f11, f13] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f21] qgq[f24, f21] SUNTF[sun, c23, c14] SUNTF[sun, c24, c21] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprop[x] . hv21] . v12 . fcc[hv22 . xprog[-x, lora, lorb]]] - ColorDelta[c11, c22] ColorDelta[c24, c21] cqq[f24, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f14] FlavorDelta[f23, f13] FlavorDelta[f24, f21] qgq[f12, f14] SUNTF[sun, c12, c14] SUNTF[sun, c23, c13] atr[fcc[xprop[x] . hv21] . v11 . DiracSigma[GAD[lora, lorb]] . fcc[hv22 . xprog[-x, lora, lorb] . v12]] - ColorDelta[c12, c21] ColorDelta[c24, c22] cqq[f24, f22] FlavorDelta[f11, f13] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f22] qgq[f11, f13] SUNTF[sun, c11, c13] SUNTF[sun, c23, c14] atr[fcc[hv21] . hv22 . xprog[-x, lora, lorb] . fcc[DiracSigma[GAD[lora, lorb]] . v12] . v11 . xprop[x]] + ColorDelta[c11, c21] ColorDelta[c24, c14] cqq[f24, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f23, f13] SUNTF[sun, c12, c22] SUNTF[sun, c23, c13] atr[xprop[x] . fcc[v11 . xprog[x, lora, lorb] . hv21]] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . v12]];


dia13=FourierXP[dia13,{x,q}]//SUNSimplify;

dia13=QEvaluate[I ScaleMu^(4-D) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type14[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia14,sun,lora,lorb},


dia14=ColorDelta[c11, c21] ColorDelta[c12, c22] cqq[f12, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f24, f14] SUNTF[sun, c23, c13] SUNTF[sun, c24, c14] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora, lorb]]] atr[xprop[x] . fcc[v11 . hv21]] + ColorDelta[c12, c21] ColorDelta[c23, c14] cqq[f23, f14] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f11, f22] SUNTF[sun, c11, c22] SUNTF[sun, c24, c13] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . v11 . xprop[x]] + ColorDelta[c12, c21] ColorDelta[c24, c13] cqq[f12, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f11, f22] SUNTF[sun, c11, c22] SUNTF[sun, c23, c14] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . v11] atr[fcc[xprop[-x] . v12] . hv22 . xprog[-x, lora, lorb]] - ColorDelta[c23, c21] ColorDelta[c24, c13] cqq[f23, f21] FlavorDelta[f11, f14] FlavorDelta[f12, f22] FlavorDelta[f23, f21] FlavorDelta[f24, f13] qgq[f11, f14] SUNTF[sun, c11, c14] SUNTF[sun, c12, c22] atr[fcc[xprop[-x] . v12] . hv22 . fcc[v11 . xprog[x, lora, lorb] . hv21] . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c23, c21] ColorDelta[c24, c13] cqq[f23, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f14] FlavorDelta[f23, f21] FlavorDelta[f24, f13] qgq[f12, f14] SUNTF[sun, c11, c22] SUNTF[sun, c12, c14] atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . DiracSigma[GAD[lora, lorb]] . fcc[xprop[-x] . v12] . hv22] - ColorDelta[c12, c14] ColorDelta[c24, c13] cqq[f12, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f14] FlavorDelta[f23, f22] FlavorDelta[f24, f13] qgq[f23, f22] SUNTF[sun, c11, c21] SUNTF[sun, c23, c22] atr[xprog[x, lora, lorb] . fcc[hv22 . DiracSigma[GAD[lora, lorb]] . hv21] . xprop[-x] . v12 . fcc[v11]] - ColorDelta[c11, c21] ColorDelta[c12, c14] cqq[f12, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f14] FlavorDelta[f23, f22] FlavorDelta[f24, f13] qgq[f23, f22] SUNTF[sun, c23, c22] SUNTF[sun, c24, c13] atr[xprop[x] . fcc[hv22 . DiracSigma[GAD[lora, lorb]] . hv21] . xprog[-x, lora, lorb] . v12 . fcc[v11]] - ColorDelta[c11, c14] ColorDelta[c24, c13] cqq[f11, f14] FlavorDelta[f11, f14] FlavorDelta[f12, f21] FlavorDelta[f23, f22] FlavorDelta[f24, f13] qgq[f23, f22] SUNTF[sun, c12, c21] SUNTF[sun, c23, c22] atr[1 . v11 . xprog[x, lora, lorb] . fcc[hv22 . DiracSigma[GAD[lora, lorb]] . hv21] . xprop[-x] . v12] + ColorDelta[c12, c22] ColorDelta[c24, c14] cqq[f24, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f23, f13] SUNTF[sun, c11, c21] SUNTF[sun, c23, c13] atr[xprog[x, lora, lorb] . fcc[v11 . xprop[x] . hv21]] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . v12]] + ColorDelta[c11, c21] ColorDelta[c23, c14] cqq[f11, f21] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f24, f13] SUNTF[sun, c12, c22] SUNTF[sun, c24, c13] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . hv22 . xprop[-x]] atr[fcc[v11 . xprog[x, lora, lorb] . hv21]];


dia14=FourierXP[dia14,{x,q}]//SUNSimplify;

dia14=QEvaluate[I ScaleMu^(4-D) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type15[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia15,sun,lora,lorb},


dia15=ColorDelta[c12, c22] ColorDelta[c24, c14] cqq[f12, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f23, f13] SUNTF[sun, c11, c21] SUNTF[sun, c23, c13] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . v12] . xprop[-x]] atr[xprog[x, lora, lorb] . fcc[v11 . hv21]] + ColorDelta[c11, c21] ColorDelta[c12, c22] cqq[f12, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f23, f13] SUNTF[sun, c23, c13] SUNTF[sun, c24, c14] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . v12] . xprog[-x, lora, lorb]] atr[xprop[x] . fcc[v11 . hv21]] + ColorDelta[c12, c22] ColorDelta[c24, c13] cqq[f24, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f11, f21] SUNTF[sun, c11, c21] SUNTF[sun, c23, c14] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprop[x] . hv21]] atr[fcc[v12] . hv22 . xprog[-x, lora, lorb]] + ColorDelta[c12, c22] ColorDelta[c23, c14] cqq[f12, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f11, f21] SUNTF[sun, c11, c21] SUNTF[sun, c24, c13] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . hv21]] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . xprop[-x]] + ColorDelta[c23, c14] ColorDelta[c24, c13] cqq[f24, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f12, f21] SUNTF[sun, c11, c22] SUNTF[sun, c12, c21] atr[fcc[v12] . hv22 . xprop[-x]] atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . DiracSigma[GAD[lora, lorb]]] - ColorDelta[c23, c14] ColorDelta[c24, c21] cqq[f24, f21] FlavorDelta[f11, f13] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f21] qgq[f11, f13] SUNTF[sun, c11, c13] SUNTF[sun, c12, c22] atr[fcc[v11 . xprog[x, lora, lorb] . hv21] . DiracSigma[GAD[lora, lorb]] . v12 . fcc[hv22 . xprop[-x]]] - ColorDelta[c23, c21] ColorDelta[c24, c14] cqq[f23, f21] FlavorDelta[f11, f22] FlavorDelta[f12, f13] FlavorDelta[f23, f21] FlavorDelta[f24, f14] qgq[f12, f13] SUNTF[sun, c11, c22] SUNTF[sun, c12, c13] atr[fcc[hv22] . xprop[-x] . fcc[v11 . DiracSigma[GAD[lora, lorb]] . v12] . xprog[x, lora, lorb] . hv21] - ColorDelta[c23, c22] ColorDelta[c24, c14] cqq[f23, f22] FlavorDelta[f11, f13] FlavorDelta[f12, f21] FlavorDelta[f23, f22] FlavorDelta[f24, f14] qgq[f11, f13] SUNTF[sun, c11, c13] SUNTF[sun, c12, c21] atr[fcc[DiracSigma[GAD[lora, lorb]] . v12] . v11 . xprog[x, lora, lorb] . fcc[hv22 . hv21] . xprop[-x]] - ColorDelta[c11, c22] ColorDelta[c12, c13] cqq[f12, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f13] FlavorDelta[f23, f14] FlavorDelta[f24, f21] qgq[f24, f21] SUNTF[sun, c23, c14] SUNTF[sun, c24, c21] atr[fcc[DiracSigma[GAD[lora, lorb]]] . hv22 . xprog[-x, lora, lorb] . fcc[v11 . v12] . xprop[x] . hv21] - ColorDelta[c11, c13] ColorDelta[c23, c14] cqq[f11, f13] FlavorDelta[f11, f13] FlavorDelta[f12, f21] FlavorDelta[f23, f14] FlavorDelta[f24, f22] qgq[f24, f22] SUNTF[sun, c12, c21] SUNTF[sun, c24, c22] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . hv22 . xprop[-x] . fcc[v12] . v11 . xprog[x, lora, lorb]];


dia15=FourierXP[dia15,{x,q}]//SUNSimplify;

dia15=QEvaluate[I ScaleMu^(4-D) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type16[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia16,sun,lora,lorb},


dia16=ColorDelta[c12, c22] ColorDelta[c23, c14] cqq[f23, f14] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f11, f21] SUNTF[sun, c11, c21] SUNTF[sun, c24, c13] atr[DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprop[x] . hv21]] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22] + ColorDelta[c12, c22] ColorDelta[c24, c13] cqq[f24, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f23, f14] SUNTF[sun, c11, c21] SUNTF[sun, c23, c14] atr[xprog[x, lora, lorb] . fcc[v11 . xprop[x] . hv21]] atr[fcc[v12] . hv22 . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c23, c13] ColorDelta[c24, c14] cqq[f23, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f11, f22] SUNTF[sun, c11, c22] SUNTF[sun, c12, c21] atr[fcc[hv22 . v12] . xprop[-x]] atr[fcc[DiracSigma[GAD[lora, lorb]] . hv21] . v11 . xprog[x, lora, lorb]] + ColorDelta[c11, c21] ColorDelta[c12, c22] cqq[f12, f22] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f23, f14] SUNTF[sun, c23, c14] SUNTF[sun, c24, c13] atr[xprop[x] . fcc[v11 . hv21]] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . DiracSigma[GAD[lora, lorb]]] + ColorDelta[c12, c21] ColorDelta[c23, c13] cqq[f23, f13] FlavorDelta[f11, f22] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f24, f14] SUNTF[sun, c11, c22] SUNTF[sun, c24, c14] atr[fcc[hv22 . v12] . DiracSigma[GAD[lora, lorb]]] atr[fcc[xprog[x, lora, lorb] . hv21] . v11 . xprop[x]] - ColorDelta[c11, c21] ColorDelta[c12, c13] cqq[f12, f13] FlavorDelta[f11, f21] FlavorDelta[f12, f13] FlavorDelta[f23, f22] FlavorDelta[f24, f14] qgq[f23, f22] SUNTF[sun, c23, c22] SUNTF[sun, c24, c14] atr[fcc[hv22 . DiracSigma[GAD[lora, lorb]] . hv21] . xprog[-x, lora, lorb] . fcc[v11 . v12] . xprop[x]] - ColorDelta[c11, c14] ColorDelta[c12, c21] cqq[f11, f14] FlavorDelta[f11, f14] FlavorDelta[f12, f21] FlavorDelta[f23, f13] FlavorDelta[f24, f22] qgq[f24, f22] SUNTF[sun, c23, c13] SUNTF[sun, c24, c22] atr[fcc[hv22 . xprog[-x, lora, lorb] . v12] . DiracSigma[GAD[lora, lorb]] . hv21 . fcc[v11 . xprop[x]]] - ColorDelta[c11, c13] ColorDelta[c12, c22] cqq[f11, f13] FlavorDelta[f11, f13] FlavorDelta[f12, f22] FlavorDelta[f23, f21] FlavorDelta[f24, f14] qgq[f23, f21] SUNTF[sun, c23, c21] SUNTF[sun, c24, c14] atr[fcc[xprog[-x, lora, lorb]] . hv22 . DiracSigma[GAD[lora, lorb]] . fcc[v11 . xprop[x] . hv21] . v12] + ColorDelta[c11, c21] ColorDelta[c23, c14] cqq[f11, f21] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f14] FlavorDelta[f24, f13] qgq[f12, f22] SUNTF[sun, c12, c22] SUNTF[sun, c24, c13] atr[fcc[xprog[-x, lora, lorb] . v12] . hv22 . xprop[-x]] atr[fcc[v11 . DiracSigma[GAD[lora, lorb]] . hv21]] + ColorDelta[c11, c21] ColorDelta[c23, c13] cqq[f11, f21] FlavorDelta[f11, f21] FlavorDelta[f12, f22] FlavorDelta[f23, f13] FlavorDelta[f24, f14] qgq[f24, f14] SUNTF[sun, c12, c22] SUNTF[sun, c24, c14] atr[fcc[hv22 . xprop[-x] . v12] . DiracSigma[GAD[lora, lorb]]] atr[fcc[v11 . xprog[x, lora, lorb] . hv21]];


dia16=FourierXP[dia16,{x,q}]//SUNSimplify;

dia16=QEvaluate[I ScaleMu^(4-D) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
