(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
TetraquarkAlphasType3::usage="TetraquarkAlphasType3[q_,j1_,j2_] give the next leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents.
Here the gluon propagator comes from quark self-energy, the renormalization involve \!\(\*SubscriptBox[\(Z\), \(2\)]\) factor."

TetraquarkAlphasType3::inderr="Dummy indices conflict!"

TetraquarkAlphasType3::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`TetraquarkAlphasType3`"]


(* ::Code::Initialization::Plain:: *)
Options[TetraquarkAlphasType3] = {
	Parallelized->True,
	HoldFlavor->False,
	Renormalization->True,
	AutoNDR->True
}


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(* propagators *)

prop[q_]=I GSD[q]FAD[q];

(*xprop[x_]=FourierPX[prop[q],{q,x}];*)
xprop[x_]=1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];

			
(* 1-loop fermion propagator *)
(*aprop[q_]=(I gStrong)^2 IntegrateP[prop[q].SUNT[colzz].GAD[loraa].prop[q+ll].GAD[lorbb].SUNT[colzz].prop[q]propg[ll,loraa,lorbb],ll];*)
(*aprop[q_] = FAD[q]^2 GSD[q] qfact1[(qfact2[-I I^D 2^(-1-D) CF (-2+D)^2 gStrong^2 \[Pi]^(-D/2)] qGamma[2-D/2] qGamma[-1+D/2]^2)/qGamma[-1+D]] SPD[q,q]^(-1+D/2);*)
(*xaprop[x_]=FourierPX[aprop[q],{q,x}];*)
xaprop[x_] = GSD[x] qfact1[(qGamma[2-D/2] qGamma[-1+D/2]^2 (qfact2[-(1/128) I CF (-2+D)^2 (2+D) gStrong^2 (-\[Pi])^-D] qGamma[-2+D]
						+qfact2[1/64 I CF (-2+D)^2 gStrong^2 (-\[Pi])^-D] qGamma[-1+D]))/(qGamma[4-D/2] qGamma[-1+D])] SPD[x,x]^(2-D);


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
TetraquarkAlphasType3[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 TetraquarkAlphasType3[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])
TetraquarkAlphasType3[qq_,factor1_ current1_FourquarkCurrent,current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 TetraquarkAlphasType3[qq,current1,current2,ops]/;FreeQ[factor1,Current]
TetraquarkAlphasType3[qq_,current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor2 TetraquarkAlphasType3[qq,current1,current2,ops]/;FreeQ[factor2,Current]


TetraquarkAlphasType3[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,atr,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,holdf=OptionValue[HoldFlavor],diagrams,ndr=OptionValue[AutoNDR],ren=OptionValue[Renormalization],fdir,files1,pall1,files2,pall2,pall},



(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[TetraquarkAlphasType3::curerr];
	Abort[]
];



(*--- check dummy indices in input ---*)
tmp2=Length[#]&/@Gather[({c11,c12,c13,c14,c21,c22,c23,c24}//FCI)/.SUNFIndex[aa_]:>aa];

If[DummyCheck[v11,v12,v21,v22]==0||!FreeQ[tmp2,nn_Integer/;nn>2],
	Message[TetraquarkAlphasType3::inderr];
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


diagrams={type1,type2};


If[ToLowerCase[ToString[OptionValue[AutoNDR]]]=="true",atr=TR5,atr=TR];

(*---------------------------------------------------*)
If[pall===True,

	DistributeDefinitions[qq,current1,current2];
	tmp= Plus@@WaitAll[Join[ParallelSubmit[{ndr,holdf,ren},#[qq,current1,current2,ndr,holdf,ren]]&/@{reno}, ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr]]&/@diagrams]];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,
	If[pall==="External",

	DistributeDefinitions[qq,current1,current2];
	
		If[fdir==="None",
			(* evaluation, no import and export *)
			Join[ParallelSubmit[{ndr,holdf,ren},#[qq,current1,current2,ndr,holdf,ren]]&/@{reno}, ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr]]&/@diagrams]
		
		,
		(* evaluation, import and export *)
			files1=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files1=("TetraquarkAlphasType3_"<>#)&/@files1;
			
			
			pall1=ImExport[fdir,
						files1,
						{{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
						{qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
						diagrams}
						];
		
		
			files2=(StringSplit[ToString[#],"`"][[-1]])&/@{reno};
			files2=("TetraquarkAlphasType3_"<>#)&/@files2;
			
			If[ren===True,
				pall2=ImExport[fdir,
						files2,
						{{ndr,holdf,ren},
						{qq,current1,current2,ndr,holdf,ren},
						{reno}}
						]
			,
				pall2={{1,1},{0}}			
			];


			{pall1[[1]]+pall2[[1]],Join[pall1[[2]],pall2[[2]]]}		
				
		](* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
(*					
		tmp= - reno + Plus@@(ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr]]&/@diagrams)*)
		
	
	,
	
		
		tmp= Plus@@Join[reno[qq,current1,current2,ndr,holdf,ren],#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr]&/@diagrams];
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]

]
(*------------------------------------------------------------------*)



(* ::Code::Initialization::Plain:: *)
reno[qq_,current1_,current2_,ndr_,holdf_,False]=0


(* ::Code::Initialization::Plain:: *)
reno[qq_,current1_,current2_,ndr_,holdf_,True]:=QGather[(CA^2-1)/(2CA) gStrong^2/(Epsilon 4 Pi^2)TetraquarkLeading[qq,current1,current2,EpsOrder->1,HoldFlavor->holdf,AutoNDR->ndr],qq,ShowasTable->False]


(* ::Code::Initialization::Plain:: *)
(*--- The following are generated by algorithm ---*)


(* ::Input::Initialization:: *)
type1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,cc,tr,dot,sunSimplify,contract},


dia1=contract[tr[dot[cc[dot[xaprop[-x], v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + contract[tr[dot[cc[dot[xprop[-x], v12]], hv22, xaprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + contract[tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xaprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + contract[tr[dot[xaprop[x], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[cc[dot[xaprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xaprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[xaprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, xaprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia1=FourierXP[dia1,{x,q}];

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
type2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,cc,tr,dot,sunSimplify,contract},


dia2=contract[tr[dot[cc[dot[xaprop[x], hv21]], v11, xprop[x]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + contract[tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xaprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + contract[tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xaprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + contract[tr[dot[xprop[x], cc[dot[v11, xaprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xaprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + contract[tr[dot[cc[dot[hv22, xaprop[-x], v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + contract[tr[dot[cc[dot[hv22, xaprop[-x], v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xaprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia2=FourierXP[dia2,{x,q}];

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
