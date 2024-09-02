(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)



(* ::Code::Initialization::Plain:: *)
Tetraquarkgg::usage="Tetraquarkgg[q_,j1_,j2_] gives the \[LeftAngleBracket]GG\[RightAngleBracket] contribution of \[LeftAngleBracket]j1 j2\[RightAngleBracket], with j1, j2 are tetraquark currents. "

Tetraquarkgg::inderr="Dummy indices conflict!"

Tetraquarkgg::curerr="Unknow current structure!"

Begin["`Private`Tetraquarkgg`"]



(* ::Code::Initialization::Plain:: *)
Options[Tetraquarkgg]={
	Parallelized->False,
	HoldFlavor->False,
	AutoNDR->True,
	ToD4->True
}



(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]

atr[gs_,ndr_,setd_]:=If[ndr,TR5[gs,ToD4->setd],TR[gs]]


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
Tetraquarkgg[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 Tetraquarkgg[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])


(* ::Code::Initialization::Plain:: *)
Tetraquarkgg[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,holdf=OptionValue[HoldFlavor],ndr=OptionValue[AutoNDR],fdir,files,pall,diagrams,setd=OptionValue[ToD4]},



(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[Tetraquarkgg::curerr];
	Abort[]
];


(*If[OptionValue[AutoNDR]===True,
	atr[gs_,op___]:=TR5[gs,op]
	,
	atr[gs_,op___]:=TR[gs]
];*)

(*--- check dummy indices in input ---*)
tmp2=Length[#]&/@Gather[({c11,c12,c13,c14,c21,c22,c23,c24}//FCI)/.SUNFIndex[aa_]:>aa];

If[DummyCheck[v11,v12,v21,v22]==0||!FreeQ[tmp2,nn_Integer/;nn>2],
	Message[Tetraquarkgg::inderr];
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


diagrams={type1,type2,type3};


(*---------------------------------------------------*)
If[pall===True,

	
	tmp=Plus@@WaitAll[ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,ndr,setd},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,ndr,setd]]&/@diagrams];
	
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

,
	If[pall==="External",

	DistributeDefinitions[qq,current1,current2];
	
		If[fdir==="None",
			(* evaluation, no import and export *)
			 ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,ndr,setd},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,ndr,setd]]&/@diagrams
		
		,
		(* evaluation, import and export *)
			files=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files=("Tetraquarkgg_"<>#)&/@files;
			
			
			ImExport[fdir,
						files,
						{{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,ndr,setd},
						{qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,ndr,setd},
						diagrams}
						]
				
		](* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
	
	,
	
		
		tmp=Plus@@(#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,ndr,setd]&/@diagrams);
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]


]

(*------------------------------------------------------------------*)



(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

(*xprop[x_]=FourierPX[prop[q],{q,x}];*)
xprop[x_]=1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];
(* the appear Condensate["gg"] = g^2<G^n_uv G^n uv> by the definition of the diagrams below, times 4Pi so that Condensate["gg"] = \alpha_s <G^n_uv G^n uv> *)


(* ::Code::Initialization::Plain:: *)
(*--- The following are generated by algorithm ---*)


(* ::Input::Initialization::Plain:: *)
type1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,ndr_,setd_]:=Block[{x,q,dia1,sun,lora,lorb,lorc,lord},


dia1=ColorDelta[c23, c14]*ColorDelta[c24, c13]*Condensate["gg"]4Pi*Contract[atr[fcc[(((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)))/(2*CA*CF*(-1 + D)*D)) . hv21] . v11 . (GAD[lorc] . GAD[lord] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lorc]*FVD[x, lord]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd]*atr[fcc[xprop[-x] . v12] . hv22 . xprop[-x],ndr,setd]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21] + ColorDelta[c23, c13]*ColorDelta[c24, c14]*Condensate["gg"]4Pi*Contract[atr[fcc[hv22 . xprop[-x] . v12] . xprop[-x],ndr,setd]*atr[fcc[(((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)))/(2*CA*CF*(-1 + D)*D)) . hv21] . v11 . (GAD[lorc] . GAD[lord] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lorc]*FVD[x, lord]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21] + ColorDelta[c23, c14]*ColorDelta[c24, c13]*Condensate["gg"]4Pi*Contract[atr[(((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)))/(2*CA*CF*(-1 + D)*D)) . fcc[v11 . (GAD[lorc] . GAD[lord] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lorc]*FVD[x, lord]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . hv21],ndr,setd]*atr[fcc[xprop[-x] . v12] . hv22 . xprop[-x],ndr,setd]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22] + ColorDelta[c23, c13]*ColorDelta[c24, c14]*Condensate["gg"]4Pi*Contract[atr[fcc[hv22 . xprop[-x] . v12] . xprop[-x],ndr,setd]*atr[(((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)))/(2*CA*CF*(-1 + D)*D)) . fcc[v11 . (GAD[lorc] . GAD[lord] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lorc]*FVD[x, lord]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . hv21],ndr,setd]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22] + ColorDelta[c11, c21]*ColorDelta[c24, c14]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[fcc[hv22 . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . v12] . xprop[-x],ndr,setd]*atr[xprop[x] . fcc[v11 . (GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . hv21],ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13] + ColorDelta[c12, c22]*ColorDelta[c23, c14]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . fcc[v11 . xprop[x] . hv21],ndr,setd]*atr[fcc[(GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . v12] . hv22 . xprop[-x],ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13] + ColorDelta[c12, c21]*ColorDelta[c23, c13]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[fcc[hv22 . xprop[-x] . v12] . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd]*atr[fcc[(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . hv21] . v11 . xprop[x],ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*Condensate["gg"]4Pi*Contract[atr[fcc[hv22 . (((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*(GAD[lora] . GAD[lorb] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lora]*FVD[-x, lorb]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)))/(2*CA*CF*(-1 + D)*D)) . v12] . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd]*atr[fcc[xprop[x] . hv21] . v11 . xprop[x],ndr,setd]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14];


dia1=FourierXP[dia1,{x,q}]//SUNSimplify;

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,ndr_,setd_]:=Block[{x,q,dia2,sun,lora,lorb,lorc,lord},


dia2=ColorDelta[c12, c21]*ColorDelta[c24, c14]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[fcc[hv22 . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . v12] . xprop[-x],ndr,setd]*atr[fcc[(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . hv21] . v11 . xprop[x],ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13] + ColorDelta[c11, c22]*ColorDelta[c24, c14]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[fcc[hv22 . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . v12] . xprop[-x],ndr,setd]*atr[fcc[xprop[x] . hv21] . v11 . (GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13] + ColorDelta[c12, c21]*ColorDelta[c24, c13]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[fcc[(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . hv21] . v11 . xprop[x],ndr,setd]*atr[fcc[xprop[-x] . v12] . hv22 . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14] + ColorDelta[c11, c22]*ColorDelta[c24, c13]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[fcc[xprop[-x] . v12] . hv22 . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd]*atr[fcc[xprop[x] . hv21] . v11 . (GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14] + ColorDelta[c12, c21]*ColorDelta[c23, c14]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[fcc[(GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . v12] . hv22 . xprop[-x],ndr,setd]*atr[fcc[(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . hv21] . v11 . xprop[x],ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13] + ColorDelta[c11, c21]*ColorDelta[c23, c14]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[xprop[x] . fcc[v11 . (GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . hv21],ndr,setd]*atr[fcc[(GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . v12] . hv22 . xprop[-x],ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13] + ColorDelta[c11, c22]*ColorDelta[c23, c13]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[fcc[hv22 . xprop[-x] . v12] . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd]*atr[fcc[xprop[x] . hv21] . v11 . (GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14] + ColorDelta[c11, c21]*ColorDelta[c23, c13]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[fcc[hv22 . xprop[-x] . v12] . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd]*atr[xprop[x] . fcc[v11 . (GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . hv21],ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14];


dia2=FourierXP[dia2,{x,q}]//SUNSimplify;

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,ndr_,setd_]:=Block[{x,q,dia3,sun,lora,lorb,lorc,lord},


dia3=ColorDelta[c12, c22]*ColorDelta[c24, c14]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[fcc[hv22 . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . v12] . xprop[-x],ndr,setd]*atr[(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . fcc[v11 . xprop[x] . hv21],ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13] + ColorDelta[c12, c22]*ColorDelta[c24, c13]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . fcc[v11 . xprop[x] . hv21],ndr,setd]*atr[fcc[xprop[-x] . v12] . hv22 . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14] + ColorDelta[c11, c21]*ColorDelta[c24, c13]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[xprop[x] . fcc[v11 . (GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . hv21],ndr,setd]*atr[fcc[xprop[-x] . v12] . hv22 . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14] + ColorDelta[c11, c22]*ColorDelta[c23, c14]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[fcc[(GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . v12] . hv22 . xprop[-x],ndr,setd]*atr[fcc[xprop[x] . hv21] . v11 . (GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*Condensate["gg"]4Pi*Contract[atr[fcc[(((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*(GAD[lora] . GAD[lorb] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lora]*FVD[-x, lorb]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)))/(2*CA*CF*(-1 + D)*D)) . v12] . hv22 . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd]*atr[fcc[xprop[x] . hv21] . v11 . xprop[x],ndr,setd]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*Condensate["gg"]4Pi*Contract[atr[xprop[x] . fcc[v11 . xprop[x] . hv21],ndr,setd]*atr[fcc[(((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*(GAD[lora] . GAD[lorb] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lora]*FVD[-x, lorb]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)))/(2*CA*CF*(-1 + D)*D)) . v12] . hv22 . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13] + ColorDelta[c12, c22]*ColorDelta[c23, c13]*Condensate["gg"]4Pi*Contract[((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*atr[fcc[hv22 . xprop[-x] . v12] . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd]*atr[(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)) . fcc[v11 . xprop[x] . hv21],ndr,setd])/(2*CA*CF*(-1 + D)*D)]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*Condensate["gg"]4Pi*Contract[atr[fcc[hv22 . (((-(MTD[lora, lord]*MTD[lorb, lorc]) + MTD[lora, lorc]*MTD[lorb, lord])*(GAD[lora] . GAD[lorb] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lora]*FVD[-x, lorb]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)))/(2*CA*CF*(-1 + D)*D)) . v12] . (GAD[lorc] . GAD[lord] . GSD[-x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[-x] . GAD[lord] . GAD[lorc]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*SPD[x, x]^(1 - D/2) + (2*FVD[-x, lorc]*FVD[-x, lord]*GSD[-x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2)),ndr,setd]*atr[xprop[x] . fcc[v11 . xprop[x] . hv21],ndr,setd]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14];


dia3=FourierXP[dia3,{x,q}]//SUNSimplify;

dia3=QEvaluate[I ScaleMu^(3(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
