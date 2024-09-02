(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
Tetraquarkqq2::usage = "Tetraquarkqq2[q_,j1_,j2_] gives leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents. "

Tetraquarkqq2::inderr = "Dummy indices conflict!"

Tetraquarkqq2::curerr = "Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`Tetraquarkqq2`"]



(* ::Code::Initialization::Plain:: *)
Options[Tetraquarkqq2] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True
}



(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(* propagators *)

prop[q_] = I GSD[q]FAD[q];

(* xprop[x_]=FourierPX[prop[q],{q,x}]; *)
xprop[x_] = 1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];

			
(*--- Condensate ---*)
cqq2[f1_,f2_,f3_,f4_]:=1/(4CA)^2 Condensate[{f1,f2}]Condensate[{f3,f4}];
cqq[f1_,f2_]:=1/(4CA)Condensate[{f1,f2}];




(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)

Tetraquarkqq2[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 Tetraquarkqq2[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])


Tetraquarkqq2[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,tr1a,tr1b,tr2a1,tr2a2,tr2b1,tr2b2,tr3a1,tr3a2,tr3b1,tr3b2,dia1,dia2,dia3,x,q,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,holdf=OptionValue[HoldFlavor],atr,fdir,files,pall,diagrams},



(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[Tetraquarkqq2::curerr];
	Abort[]
];


(*--- check dummy indices in input ---*)
tmp2=Length[#]&/@Gather[({c11,c12,c13,c14,c21,c22,c23,c24}//FCI)/.SUNFIndex[aa_]:>aa];

If[DummyCheck[v11,v12,v21,v22]==0||!FreeQ[tmp2,nn_Integer/;nn>2],
	Message[Tetraquarkqq2::inderr];
	Abort[]
];


If[OptionValue[AutoNDR]===True,
	atr=TR5
	,
	atr=TR
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


diagrams={xtype1,xtype2,xtype3,xtype4,xtype5};


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
			files=("Tetraquarkqq2_"<>#)&/@files;
			
			
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
xtype1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,cc,tr,dot,sunSimplify,contract},


 dia1=Condensate[{f21, f12}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv21]], v11]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f13, f23}]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f14, f24}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f14, f24}]*Condensate[{f21, f11}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f21, f11}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[v11, hv21]]]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f14, f24}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f13, f23}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f13, f11}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[hv21]], hv22, xprop[-x], cc[dot[v12]], v11, xprop[x]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c22])/(16*CA^2)];


 dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia1=FourierXP[dia1,{x,q}];


 dia1=QEvaluate[I ScaleMu^(1(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
xtype2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,cc,tr,dot,sunSimplify,contract},


 dia2=Condensate[{f14, f23}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f14, f23}]*Condensate[{f21, f11}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f14, f23}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f13, f24}]*Condensate[{f22, f12}]*contract[tr[dot[xprop[x], cc[dot[v11, hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f13, f23}]*Condensate[{f21, f11}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[hv22, v12]], xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f13, f11}]*Condensate[{f22, f23}]*contract[tr[dot[cc[dot[v12]], v11, xprop[x], cc[dot[hv22, hv21]], xprop[-x]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c21]*ColorDelta[c23, c22]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f14, f12}]*Condensate[{f21, f24}]*contract[tr[dot[cc[dot[xprop[x], hv21]], v11, cc[dot[hv22, xprop[-x], v12]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c14]*ColorDelta[c23, c13]*ColorDelta[c24, c21])/(16*CA^2)] + Condensate[{f14, f11}]*Condensate[{f21, f24}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c21])/(16*CA^2)];


 dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia2=FourierXP[dia2,{x,q}];


 dia2=QEvaluate[I ScaleMu^(1(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
xtype3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,cc,tr,dot,sunSimplify,contract},


 dia3=Condensate[{f13, f24}]*Condensate[{f14, f23}]*contract[tr[dot[cc[dot[v12]], hv22]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f13, f24}]*Condensate[{f14, f23}]*contract[tr[dot[cc[dot[v12]], hv22]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f14, f12}]*Condensate[{f21, f23}]*contract[tr[dot[cc[dot[xprop[x], hv21]], v11, cc[dot[xprop[-x], v12]], hv22]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c14]*ColorDelta[c23, c21]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f13, f23}]*Condensate[{f14, f24}]*contract[tr[dot[cc[dot[hv22, v12]]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f13, f12}]*Condensate[{f21, f23}]*contract[tr[dot[cc[dot[hv22]], xprop[-x], cc[dot[v11, v12]], xprop[x], hv21]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c13]*ColorDelta[c23, c21]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f13, f12}]*Condensate[{f21, f24}]*contract[tr[dot[cc[dot[]], hv22, xprop[-x], cc[dot[v11, v12]], xprop[x], hv21]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c13]*ColorDelta[c23, c14]*ColorDelta[c24, c21])/(16*CA^2)] + Condensate[{f13, f11}]*Condensate[{f21, f24}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]], v12, cc[dot[hv22, xprop[-x]]]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c21])/(16*CA^2)] + Condensate[{f13, f12}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[hv21]], hv22, xprop[-x], cc[dot[v11, v12]], xprop[x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c13]*ColorDelta[c23, c14]*ColorDelta[c24, c22])/(16*CA^2)];


 dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia3=FourierXP[dia3,{x,q}];


 dia3=QEvaluate[I ScaleMu^(1(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
xtype4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,cc,tr,dot,sunSimplify,contract},


 dia4=Condensate[{f14, f23}]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[cc[dot[xprop[x], hv21]], v11]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f13, f24}]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[xprop[x], hv21]], v11]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f13, f24}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv21]], v11, xprop[x]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f14, f11}]*Condensate[{f21, f23}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22, cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c22]*ColorDelta[c23, c21]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f14, f11}]*Condensate[{f22, f23}]*contract[tr[dot[cc[dot[]], v11, xprop[x], cc[dot[hv22, hv21]], xprop[-x], v12]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c21]*ColorDelta[c23, c22]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f21, f12}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv21]], v11]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f13, f23}]*Condensate[{f14, f24}]*contract[tr[dot[cc[dot[hv22, v12]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f13, f12}]*Condensate[{f22, f23}]*contract[tr[dot[cc[dot[hv22, hv21]], xprop[-x], cc[dot[v11, v12]], xprop[x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c13]*ColorDelta[c23, c22]*ColorDelta[c24, c14])/(16*CA^2)];


 dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia4=FourierXP[dia4,{x,q}];


 dia4=QEvaluate[I ScaleMu^(1(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
xtype5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia5,cc,tr,dot,sunSimplify,contract},


 dia5=Condensate[{f13, f24}]*Condensate[{f21, f11}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f21, f11}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[v11, hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f14, f12}]*Condensate[{f22, f23}]*contract[tr[dot[xprop[x], cc[dot[hv22, hv21]], xprop[-x], v12, cc[dot[v11]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c14]*ColorDelta[c23, c22]*ColorDelta[c24, c13])/(16*CA^2)] + Condensate[{f14, f24}]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f13, f23}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f13, f11}]*Condensate[{f21, f23}]*contract[tr[dot[cc[dot[xprop[-x]]], hv22, cc[dot[v11, xprop[x], hv21]], v12]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c22]*ColorDelta[c23, c21]*ColorDelta[c24, c14])/(16*CA^2)] + Condensate[{f14, f12}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[xprop[x]]], v11, cc[dot[hv22, xprop[-x], v12]], hv21]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c14]*ColorDelta[c23, c13]*ColorDelta[c24, c22])/(16*CA^2)] + Condensate[{f14, f11}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], hv21, cc[dot[v11, xprop[x]]]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c22])/(16*CA^2)];


 dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia5=FourierXP[dia5,{x,q}];


 dia5=QEvaluate[I ScaleMu^(1(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
(*xtype1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,cc,tr,dot,sunSimplify,contract},


dia1=contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*cqq[f11, f22]*cqq[f23, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c13, c24]*ColorDelta[c14, c23]] + contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*cqq[f12, f22]*cqq[f23, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c13, c24]*ColorDelta[c14, c23]] + contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*cqq[f11, f22]*cqq[f23, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c13, c23]*ColorDelta[c14, c24]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11]]]*cqq[f12, f21]*cqq[f24, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c13, c23]*ColorDelta[c14, c24]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*cqq[f12, f22]*cqq[f24, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c13, c23]*ColorDelta[c14, c24]] + contract[tr[dot[cc[dot[xprop[x], hv21]], v11, cc[dot[hv22, xprop[-x], v12]]]]]*cqq[f12, f14]*cqq[f24, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f13, f23]*FlavorDelta[f21, f24]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c14]*ColorDelta[c13, c23]*ColorDelta[c21, c24])] + contract[tr[dot[v11, xprop[x], cc[dot[hv22, hv21]], xprop[-x], v12]]]*cqq[f11, f14]*cqq[f23, f22]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f22, f23]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c21]*ColorDelta[c13, c24]*ColorDelta[c22, c23])] + contract[tr[dot[cc[dot[xprop[x]]], v11, cc[dot[hv22, xprop[-x], v12]], hv21]]]*cqq[f12, f14]*cqq[f24, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f13, f23]*FlavorDelta[f22, f24]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c14]*ColorDelta[c13, c23]*ColorDelta[c22, c24])];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


dia1=FourierXP[dia1,{x,q}];

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*xtype2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,cc,tr,dot,sunSimplify,contract},


dia2=contract[tr[dot[xprop[x], cc[dot[v11, hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*cqq[f12, f22]*cqq[f24, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c13, c24]*ColorDelta[c14, c23]] + contract[tr[dot[cc[dot[v12]], hv22]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*cqq[f23, f14]*cqq[f24, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c13, c24]*ColorDelta[c14, c23]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*cqq[f11, f22]*cqq[f24, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c13, c23]*ColorDelta[c14, c24]] + contract[tr[dot[cc[dot[hv22, v12]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*cqq[f23, f13]*cqq[f24, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c13, c23]*ColorDelta[c14, c24]] + contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*cqq[f12, f22]*cqq[f23, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c13, c23]*ColorDelta[c14, c24]] + contract[tr[dot[cc[dot[xprop[-x]]], hv22, cc[dot[v11, xprop[x], hv21]], v12]]]*cqq[f11, f13]*cqq[f23, f21]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f14, f24]*FlavorDelta[f21, f23]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c22]*ColorDelta[c14, c24]*ColorDelta[c21, c23])] + contract[tr[dot[hv22, xprop[-x], cc[dot[v11, v12]], xprop[x], hv21]]]*cqq[f12, f13]*cqq[f24, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f14, f23]*FlavorDelta[f21, f24]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c13]*ColorDelta[c14, c23]*ColorDelta[c21, c24])] + contract[tr[dot[xprop[x], cc[dot[hv22, hv21]], xprop[-x], v12, cc[dot[v11]]]]]*cqq[f12, f14]*cqq[f23, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f13, f24]*FlavorDelta[f22, f23]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c14]*ColorDelta[c13, c24]*ColorDelta[c22, c23])];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


dia2=FourierXP[dia2,{x,q}];

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*xtype3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,cc,tr,dot,sunSimplify,contract},


dia3=contract[tr[dot[cc[dot[hv21]], v11]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*cqq[f11, f22]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c13, c24]*ColorDelta[c14, c23]] + contract[tr[dot[cc[dot[v12]], hv22]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*cqq[f23, f14]*cqq[f24, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c13, c24]*ColorDelta[c14, c23]] + contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22]]]*cqq[f11, f21]*cqq[f23, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c13, c24]*ColorDelta[c14, c23]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[v11, xprop[x], hv21]]]]]*cqq[f11, f21]*cqq[f24, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c13, c23]*ColorDelta[c14, c24]] + contract[tr[dot[cc[dot[hv22, v12]]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*cqq[f23, f13]*cqq[f24, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c13, c23]*ColorDelta[c14, c24]] + contract[tr[dot[cc[dot[hv22]], xprop[-x], cc[dot[v11, v12]], xprop[x], hv21]]]*cqq[f12, f13]*cqq[f23, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f14, f24]*FlavorDelta[f21, f23]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c13]*ColorDelta[c14, c24]*ColorDelta[c21, c23])] + contract[tr[dot[cc[dot[hv22, hv21]], xprop[-x], cc[dot[v11, v12]], xprop[x]]]]*cqq[f12, f13]*cqq[f23, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f14, f24]*FlavorDelta[f22, f23]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c13]*ColorDelta[c14, c24]*ColorDelta[c22, c23])] + contract[tr[dot[cc[dot[hv21]], hv22, xprop[-x], cc[dot[v11, v12]], xprop[x]]]]*cqq[f12, f13]*cqq[f24, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f14, f23]*FlavorDelta[f22, f24]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c13]*ColorDelta[c14, c23]*ColorDelta[c22, c24])];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


dia3=FourierXP[dia3,{x,q}];

dia3=QEvaluate[I ScaleMu^(3(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*xtype4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,cc,tr,dot,sunSimplify,contract},


dia4=contract[tr[dot[cc[dot[hv21]], v11, xprop[x]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*cqq[f11, f22]*cqq[f24, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c13, c24]*ColorDelta[c14, c23]] + contract[tr[dot[cc[dot[v11, hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*cqq[f11, f21]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c13, c24]*ColorDelta[c14, c23]] + contract[tr[dot[cc[dot[hv21]], v11]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]]*cqq[f11, f22]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c13, c23]*ColorDelta[c14, c24]] + contract[tr[dot[cc[dot[v11, hv21]]]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]]*cqq[f11, f21]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c13, c23]*ColorDelta[c14, c24]] + contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[hv22, v12]], xprop[-x]]]]*cqq[f11, f21]*cqq[f23, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c13, c23]*ColorDelta[c14, c24]] + contract[tr[dot[cc[dot[xprop[x], hv21]], v11, cc[dot[xprop[-x], v12]], hv22]]]*cqq[f12, f14]*cqq[f23, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f13, f24]*FlavorDelta[f21, f23]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c14]*ColorDelta[c13, c24]*ColorDelta[c21, c23])] + contract[tr[dot[cc[dot[xprop[-x], v12]], hv22, cc[dot[v11, xprop[x], hv21]]]]]*cqq[f11, f14]*cqq[f23, f21]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f21, f23]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c22]*ColorDelta[c13, c24]*ColorDelta[c21, c23])] + contract[tr[dot[cc[dot[v12]], v11, xprop[x], cc[dot[hv22, hv21]], xprop[-x]]]]*cqq[f11, f13]*cqq[f23, f22]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f14, f24]*FlavorDelta[f22, f23]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c21]*ColorDelta[c14, c24]*ColorDelta[c22, c23])];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


dia4=FourierXP[dia4,{x,q}];

dia4=QEvaluate[I ScaleMu^(3(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*xtype5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia5,cc,tr,dot,sunSimplify,contract},


dia5=contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[cc[dot[xprop[x], hv21]], v11]]]*cqq[f12, f21]*cqq[f23, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c13, c24]*ColorDelta[c14, c23]] + contract[tr[dot[cc[dot[xprop[x], hv21]], v11]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*cqq[f12, f21]*cqq[f24, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c13, c24]*ColorDelta[c14, c23]] + contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*cqq[f11, f21]*cqq[f24, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f13, f24]*FlavorDelta[f14, f23]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c13, c24]*ColorDelta[c14, c23]] + contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11]]]*cqq[f12, f21]*cqq[f23, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f14, f24]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c13, c23]*ColorDelta[c14, c24]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], cc[dot[v11, xprop[x], hv21]]]]]*cqq[f11, f14]*cqq[f24, f21]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f13, f23]*FlavorDelta[f21, f24]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c22]*ColorDelta[c13, c23]*ColorDelta[c21, c24])] + contract[tr[dot[cc[dot[v11, xprop[x], hv21]], v12, cc[dot[hv22, xprop[-x]]]]]]*cqq[f11, f13]*cqq[f24, f21]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f14, f23]*FlavorDelta[f21, f24]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c22]*ColorDelta[c14, c23]*ColorDelta[c21, c24])] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], hv21, cc[dot[v11, xprop[x]]]]]]*cqq[f11, f14]*cqq[f24, f22]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f13, f23]*FlavorDelta[f22, f24]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c21]*ColorDelta[c13, c23]*ColorDelta[c22, c24])] + contract[tr[dot[cc[dot[hv21]], hv22, xprop[-x], cc[dot[v12]], v11, xprop[x]]]]*cqq[f11, f13]*cqq[f24, f22]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f14, f23]*FlavorDelta[f22, f24]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c21]*ColorDelta[c14, c23]*ColorDelta[c22, c24])];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


dia5=FourierXP[dia5,{x,q}];

dia5=QEvaluate[I ScaleMu^(3(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Code::Initialization::Plain:: *)
End[]

(*EndPackage[]*)
