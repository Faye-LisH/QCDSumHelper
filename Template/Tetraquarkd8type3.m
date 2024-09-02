(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


Tetraquarkd8type3::usage="Tetraquarkd8type3[q_,j1_,j2_] gives the second type of \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)q\[RightAngleBracket]\[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)Gq\[RightAngleBracket] contribution of \[LeftAngleBracket]j1 j2\[RightAngleBracket], with j1, j2 are tetraquark currents. "

Tetraquarkd8type3::inderr="Dummy indices conflict!"

Tetraquarkd8type3::curerr="Unknow current structure!"


Begin["`Private`Tetraquarkd8type3`"]


Options[Tetraquarkd8type3] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True
}


(* VS-first *)


(*------------------------------------------------------------------*)
(*--- propagators ---*)


xprop[x_]=(I^(1 - D)*DiracGamma[Momentum[x, D], D]*qGamma[D/2])/(2*Pi^(D/2)*Pair[Momentum[x, D], Momentum[x, D]]^(D/2));(*FourierPX[prop[q],{q,x}]*)

(*--- propagator with one background gluon ---*)

(*xprog[x_,lora_,lorb_,cola_]=FourierPX[prog[q,lora,lorb,cola],{q,x}];*)
xprog[x_,lora_,lorb_]=(GAD[lora] . GAD[lorb] . GSD[x] qfact1[1/32 I^-D \[Pi]^(-D/2) qGamma[-1+D/2]] SPD[x,x]^(1-D/2) 
							+GSD[x] . GAD[lorb] . GAD[lora] qfact1[qfact2[-(1/32) I^-D \[Pi]^(-D/2)] qGamma[-1+D/2]] SPD[x,x]^(1-D/2) 
							+2 FVD[x,lora] FVD[x,lorb] GSD[x] qfact1[-(1/16) I^-D \[Pi]^(-D/2) qGamma[D/2]] SPD[x,x]^(-D/2) );

(*--- Condensate ---*)

cqq[f2_,f1_]=-1/(4CA)Condensate[{f2,f1}];
cqgq[x_,f2_,f1_]=-SPD[x]/(2^4 D CA)Condensate[{f2,"G",f1}];	



(*------------------------------------------------------------------*)

fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]


(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
Tetraquarkd8type3[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 Tetraquarkd8type3[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])


Tetraquarkd8type3[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,holdf=OptionValue[HoldFlavor],atr,fdir,files,pall,diagrams},



(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[Tetraquarkd8type3::curerr];
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
	Message[Tetraquarkd8type3::inderr];
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


diagrams={xtype1,xtype2,xtype3,xtype4,xtype5,xtype6,xtype7,xtype8};


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
			files=("Tetraquarkd8type3_"<>#)&/@files;
			
			
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



(*--- The following are generated by algorithm ---*)


(* ::Input::Initialization:: *)
xtype1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,cc,tr,dot,sunSimplify,contract},


 dia1=Condensate[{f22, f11}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[1, hv21]], v11]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f13, f24}]*Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[1, v12]], hv22]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f13, f24}]*Condensate[{f22, "G", f12}]*contract[tr[dot[xprop[x], cc[dot[v11, hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f21, f11}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[v11, 1, hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f13, f23}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f14, f24}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)];


 dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia1=FourierXP[dia1,{x,q}];


 dia1=QEvaluate[I ScaleMu^(1(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)/.q->qq


]


(* ::Input::Initialization:: *)
xtype2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,cc,tr,dot,sunSimplify,contract},


 dia2=Condensate[{f13, f24}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv21]], v11, xprop[x]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f14, f23}]*Condensate[{f13, "G", f24}]*contract[tr[dot[cc[dot[1, v12]], hv22]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f21, f11}]*Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f14, f23}]*Condensate[{f13, "G", f24}]*contract[tr[dot[cc[dot[1, v12]], hv22]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f13, f23}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[hv22, v12]], xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f14, f24}]*Condensate[{f13, "G", f23}]*contract[tr[dot[cc[dot[hv22, 1, v12]]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)];


 dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia2=FourierXP[dia2,{x,q}];


 dia2=QEvaluate[I ScaleMu^(1(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)/.q->qq


]


(* ::Input::Initialization:: *)
xtype3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,cc,tr,dot,sunSimplify,contract},


 dia3=Condensate[{f21, f12}]*Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[cc[dot[xprop[x], hv21]], v11]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f13, f24}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprop[x], hv21]], v11]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f21, f11}]*Condensate[{f13, "G", f24}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f13, f23}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, 1, v12]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f22, f12}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f13, f23}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)];


 dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia3=FourierXP[dia3,{x,q}];


 dia3=QEvaluate[I ScaleMu^(1(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)/.q->qq


]


(* ::Input::Initialization:: *)
xtype4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,cc,tr,dot,sunSimplify,contract},


 dia4=Condensate[{f21, f12}]*Condensate[{f13, "G", f24}]*contract[tr[dot[cc[dot[xprop[x], hv21]], v11]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f13, f24}]*Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[1, v12]], hv22]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f22, f12}]*Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f14, f23}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f22, f11}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f13, f23}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)];


 dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia4=FourierXP[dia4,{x,q}];


 dia4=QEvaluate[I ScaleMu^(1(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)/.q->qq


]


(* ::Input::Initialization:: *)
xtype5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia5,cc,tr,dot,sunSimplify,contract},


 dia5=Condensate[{f14, f23}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f21, f12}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[1, hv21]], v11]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f22, f12}]*Condensate[{f13, "G", f24}]*contract[tr[dot[xprop[x], cc[dot[v11, hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f22, f12}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, 1, hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f21, f12}]*Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[1, hv21]], v11]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f22, f11}]*Condensate[{f13, "G", f23}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)];


 dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia5=FourierXP[dia5,{x,q}];


 dia5=QEvaluate[I ScaleMu^(1(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)/.q->qq


]


(* ::Input::Initialization:: *)
xtype6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia6,cc,tr,dot,sunSimplify,contract},


 dia6=Condensate[{f14, f23}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f14, f24}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f21, f12}]*Condensate[{f13, "G", f23}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f21, f11}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f14, f24}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f22, f12}]*Condensate[{f13, "G", f23}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)];


 dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia6=FourierXP[dia6,{x,q}];


 dia6=QEvaluate[I ScaleMu^(1(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)/.q->qq


]


(* ::Input::Initialization:: *)
xtype7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia7,cc,tr,dot,sunSimplify,contract},


 dia7=Condensate[{f14, f23}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[cc[dot[xprop[x], hv21]], v11]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f22, f11}]*Condensate[{f13, "G", f24}]*contract[tr[dot[cc[dot[hv21]], v11, xprop[x]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f21, f12}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f22, f11}]*Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[1, hv21]], v11]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f14, f24}]*Condensate[{f13, "G", f23}]*contract[tr[dot[cc[dot[hv22, 1, v12]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f14, f24}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)];


 dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia7=FourierXP[dia7,{x,q}];


 dia7=QEvaluate[I ScaleMu^(1(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)/.q->qq


]


(* ::Input::Initialization:: *)
xtype8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia8,cc,tr,dot,sunSimplify,contract},


 dia8=Condensate[{f22, f11}]*Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f13, f24}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*SPD[x])/(64*CA^2*D)] + Condensate[{f21, f11}]*Condensate[{f13, "G", f23}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[hv22, v12]], xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f22, f12}]*Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[v11, 1, hv21]]]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f21, f11}]*Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[v11, 1, hv21]]]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)] + Condensate[{f13, f23}]*Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, 1, v12]]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[(ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*SPD[x])/(64*CA^2*D)];


 dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia8=FourierXP[dia8,{x,q}];


 dia8=QEvaluate[I ScaleMu^(1(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)/.q->qq


]


(* ::Input::Initialization:: *)
(*type1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,sun,lora,lorb},


dia1=ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f23, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f12, f21]*SPD[x]*atr[fcc[hv22 . v12] . xprop[-x]]*atr[fcc[xprop[x] . hv21] . v11] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f23, f13]*SPD[x]*atr[fcc[hv22 . v12] . xprop[-x]]*atr[xprop[x] . fcc[v11 . hv21]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f23, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f24, f13]*SPD[x]*atr[fcc[1 . v12] . hv22]*atr[xprop[x] . fcc[v11 . xprop[x] . hv21]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f23, f13]*SPD[x]*atr[fcc[hv22 . v12] . xprop[-x]]*atr[fcc[hv21] . v11 . xprop[x]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f23, f14]*SPD[x]*atr[fcc[xprop[-x] . v12] . hv22]*atr[fcc[hv21] . v11 . xprop[x]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f24, f14]*SPD[x]*atr[fcc[xprop[x] . hv21] . v11]*atr[fcc[hv22 . xprop[-x] . v12]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f24, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f12, f22]*SPD[x]*atr[xprop[x] . fcc[v11 . hv21]]*atr[fcc[hv22 . xprop[-x] . v12]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f24, f14]*SPD[x]*atr[fcc[hv21] . v11 . xprop[x]]*atr[fcc[hv22 . xprop[-x] . v12]];


dia1=(FourierXP[dia1,{x,q}]//SUNSimplify);

dia1=QEvaluate[I ScaleMu^(4-D) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,sun,lora,lorb},


dia2=ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f23, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f12, f21]*SPD[x]*atr[fcc[xprop[-x] . v12] . hv22]*atr[fcc[xprop[x] . hv21] . v11] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f23, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f11, f22]*SPD[x]*atr[fcc[hv22 . v12] . xprop[-x]]*atr[fcc[hv21] . v11 . xprop[x]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f24, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f12, f22]*SPD[x]*atr[xprop[x] . fcc[v11 . hv21]]*atr[fcc[v12] . hv22 . xprop[-x]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f24, f13]*SPD[x]*atr[fcc[hv21] . v11 . xprop[x]]*atr[fcc[v12] . hv22 . xprop[-x]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f11, f22]*SPD[x]*atr[fcc[1 . hv21] . v11]*atr[fcc[xprop[-x] . v12] . hv22 . xprop[-x]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f23, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f24, f13]*SPD[x]*atr[fcc[1 . v12] . hv22]*atr[fcc[xprop[x] . hv21] . v11 . xprop[x]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f23, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f24, f14]*SPD[x]*atr[fcc[xprop[x] . hv21] . v11 . xprop[x]]*atr[fcc[hv22 . 1 . v12]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f23, f13]*SPD[x]*atr[fcc[hv22 . v12] . xprop[-x]]*atr[fcc[v11 . xprop[x] . hv21]];


dia2=(FourierXP[dia2,{x,q}]//SUNSimplify);

dia2=QEvaluate[I ScaleMu^(4-D) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,sun,lora,lorb},


dia3=ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f23, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f11, f22]*SPD[x]*atr[fcc[xprop[-x] . v12] . hv22]*atr[fcc[hv21] . v11 . xprop[x]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f24, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f12, f21]*SPD[x]*atr[fcc[xprop[x] . hv21] . v11]*atr[fcc[v12] . hv22 . xprop[-x]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f12, f21]*SPD[x]*atr[fcc[1 . hv21] . v11]*atr[fcc[xprop[-x] . v12] . hv22 . xprop[-x]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f24, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f23, f14]*SPD[x]*atr[fcc[1 . v12] . hv22]*atr[fcc[xprop[x] . hv21] . v11 . xprop[x]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f24, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f11, f22]*SPD[x]*atr[fcc[hv21] . v11 . xprop[x]]*atr[fcc[hv22 . xprop[-x] . v12]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f11, f21]*SPD[x]*atr[fcc[xprop[-x] . v12] . hv22 . xprop[-x]]*atr[fcc[v11 . 1 . hv21]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f23, f14]*SPD[x]*atr[fcc[xprop[-x] . v12] . hv22]*atr[fcc[v11 . xprop[x] . hv21]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f24, f13]*SPD[x]*atr[fcc[v12] . hv22 . xprop[-x]]*atr[fcc[v11 . xprop[x] . hv21]];


dia3=(FourierXP[dia3,{x,q}]//SUNSimplify);

dia3=QEvaluate[I ScaleMu^(4-D) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,sun,lora,lorb},


dia4=ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f23, f14]*SPD[x]*atr[fcc[xprop[-x] . v12] . hv22]*atr[fcc[xprop[x] . hv21] . v11] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f11, f22]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f12, f21]*SPD[x]*atr[fcc[1 . hv21] . v11]*atr[fcc[hv22 . xprop[-x] . v12] . xprop[-x]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f23, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f24, f14]*SPD[x]*atr[xprop[x] . fcc[v11 . xprop[x] . hv21]]*atr[fcc[hv22 . 1 . v12]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f24, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f23, f13]*SPD[x]*atr[fcc[xprop[x] . hv21] . v11 . xprop[x]]*atr[fcc[hv22 . 1 . v12]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f11, f21]*SPD[x]*atr[fcc[hv22 . xprop[-x] . v12] . xprop[-x]]*atr[fcc[v11 . 1 . hv21]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f12, f22]*SPD[x]*atr[fcc[hv22 . xprop[-x] . v12] . xprop[-x]]*atr[fcc[v11 . 1 . hv21]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f23, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f11, f21]*SPD[x]*atr[fcc[hv22 . v12] . xprop[-x]]*atr[fcc[v11 . xprop[x] . hv21]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f24, f14]*SPD[x]*atr[fcc[hv22 . xprop[-x] . v12]]*atr[fcc[v11 . xprop[x] . hv21]];


dia4=(FourierXP[dia4,{x,q}]//SUNSimplify);

dia4=QEvaluate[I ScaleMu^(4-D) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia5,sun,lora,lorb},


dia5=ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f23, f13]*SPD[x]*atr[fcc[hv22 . v12] . xprop[-x]]*atr[fcc[xprop[x] . hv21] . v11] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f11, f22]*SPD[x]*atr[fcc[1 . hv21] . v11]*atr[fcc[hv22 . xprop[-x] . v12] . xprop[-x]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f23, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f12, f22]*SPD[x]*atr[fcc[xprop[-x] . v12] . hv22]*atr[xprop[x] . fcc[v11 . hv21]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f24, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f23, f14]*SPD[x]*atr[fcc[1 . v12] . hv22]*atr[xprop[x] . fcc[v11 . xprop[x] . hv21]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f12, f21]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f24, f13]*SPD[x]*atr[fcc[xprop[x] . hv21] . v11]*atr[fcc[v12] . hv22 . xprop[-x]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f24, f13]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f11, f22]*SPD[x]*atr[fcc[hv21] . v11 . xprop[x]]*atr[fcc[v12] . hv22 . xprop[-x]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f23, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f11, f21]*SPD[x]*atr[fcc[xprop[-x] . v12] . hv22]*atr[fcc[v11 . xprop[x] . hv21]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f24, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f11, f21]*SPD[x]*atr[fcc[hv22 . xprop[-x] . v12]]*atr[fcc[v11 . xprop[x] . hv21]];


dia5=(FourierXP[dia5,{x,q}]//SUNSimplify);

dia5=QEvaluate[I ScaleMu^(4-D) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*type6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia6,sun,lora,lorb},


dia6=ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f23, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f12, f22]*SPD[x]*atr[fcc[hv22 . v12] . xprop[-x]]*atr[xprop[x] . fcc[v11 . hv21]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f23, f14]*SPD[x]*atr[fcc[xprop[-x] . v12] . hv22]*atr[xprop[x] . fcc[v11 . hv21]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f24, f13]*SPD[x]*atr[xprop[x] . fcc[v11 . hv21]]*atr[fcc[v12] . hv22 . xprop[-x]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f24, f14]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f23, f13]*SPD[x]*atr[xprop[x] . fcc[v11 . xprop[x] . hv21]]*atr[fcc[hv22 . 1 . v12]] + ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f24, f14]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f12, f21]*SPD[x]*atr[fcc[xprop[x] . hv21] . v11]*atr[fcc[hv22 . xprop[-x] . v12]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]*cqq[f12, f22]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*qgq3[f24, f14]*SPD[x]*atr[xprop[x] . fcc[v11 . hv21]]*atr[fcc[hv22 . xprop[-x] . v12]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f11, f21]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f12, f22]*SPD[x]*atr[fcc[xprop[-x] . v12] . hv22 . xprop[-x]]*atr[fcc[v11 . 1 . hv21]] + ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]*cqq[f24, f13]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*qgq3[f11, f21]*SPD[x]*atr[fcc[v12] . hv22 . xprop[-x]]*atr[fcc[v11 . xprop[x] . hv21]];


dia6=(FourierXP[dia6,{x,q}]//SUNSimplify);

dia6=QEvaluate[I ScaleMu^(4-D) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


End[]
(*EndPackage[]*)
