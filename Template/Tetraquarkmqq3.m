(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
Tetraquarkmqq3::usage="Tetraquarkmqq3[q_,j1_,j2_] gives the third type of \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)Gq\!\(\*SuperscriptBox[\(\[RightAngleBracket]\), \(2\)]\) contribution of \[LeftAngleBracket]j1 j2\[RightAngleBracket], with j1, j2 are tetraquark currents. "

Tetraquarkmqq3::inderr="Dummy indices conflict!"

Tetraquarkmqq3::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`Tetraquarkmqq3`"]


(* ::Code::Initialization::Plain:: *)
Options[Tetraquarkmqq3] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True,
	ToD4->"Auto"
}


(* ::Code::Initialization::Plain:: *)
(* VS-first *)


(* ::Code::Initialization::Plain:: *)
(*xprop[x_]=(I^(1 - D)*DiracGamma[Momentum[x, D], D]*qGamma[D/2])/(2*Pi^(D/2)*Pair[Momentum[x, D], Momentum[x, D]]^(D/2));(*FourierPX[prop[q],{q,x}]*)*)

cqqm[x_,f_]=-1/4*1/CA + ((I/4)*DiracGamma[Momentum[x, D], D]*quarkMass[f])/(CA*D);

(*-----------------------------------*)
(* propagators and condensates *)

xprop[x_,f_]=(DiracGamma[Momentum[x, D], D]*qfact1[(I^(1 - D)*qGamma[D/2])/(2*Pi^(D/2))])/Pair[Momentum[x, D], Momentum[x, D]]^(D/2) + Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[(I^(2 - D)*qGamma[-1 + D/2]*quarkMass[f])/(4*Pi^(D/2))];


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
Tetraquarkmqq3[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 Tetraquarkmqq3[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])


(* ::Code::Initialization::Plain:: *)
Tetraquarkmqq3[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,holdf=OptionValue[HoldFlavor],atr,fdir,files,pall,diagrams,setd},



(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[Tetraquarkmqq3::curerr];
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
	Message[Tetraquarkmqq3::inderr];
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


diagrams={xtype1,xtype2,xtype3,xtype4,xtype5,xtype6};


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
			files=("Tetraquarkmqq3_"<>#)&/@files;
			
			
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



(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(*--- The following are generated by algorithm ---*)
(*------------------------------------------------------------------*)


(* ::Input::Initialization::Plain:: *)
xtype1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,cc,tr,dot,lor,sun,sunSimplify,contract},


dia1=Condensate[{f13, f24}]*Condensate[{f14, f23}]*Condensate[{f21, f11}]*contract[tr[dot[cqqm[x, f11], cc[dot[v11, xprop[x, f12], hv21]]]]*tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, cqqm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f13, f24}]*Condensate[{f14, f12}]*Condensate[{f21, f23}]*contract[tr[dot[cc[dot[xprop[x, f11], hv21]], v11, cqqm[0, f12], cc[dot[cqqm[-x, f24], v12]], hv22, cqqm[0, f23]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c14]*ColorDelta[c23, c21]*ColorDelta[c24, c13])] + Condensate[{f14, f12}]*Condensate[{f21, f11}]*Condensate[{f22, f23}]*contract[tr[dot[cqqm[x, f11], cc[dot[hv22, cqqm[0, f23], hv21]], xprop[-x, f24], v12, cc[dot[v11, cqqm[0, f12]]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c14]*ColorDelta[c23, c22]*ColorDelta[c24, c13])] + Condensate[{f13, f24}]*Condensate[{f14, f12}]*Condensate[{f22, f23}]*contract[tr[dot[xprop[x, f11], cc[dot[hv22, cqqm[0, f23], hv21]], cqqm[-x, f24], v12, cc[dot[v11, cqqm[0, f12]]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c14]*ColorDelta[c23, c22]*ColorDelta[c24, c13])] + Condensate[{f13, f12}]*Condensate[{f14, f23}]*Condensate[{f21, f24}]*contract[tr[dot[cc[dot[cqqm[0, f24]]], hv22, cqqm[-x, f23], cc[dot[v11, cqqm[0, f12], v12]], xprop[x, f11], hv21]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c13]*ColorDelta[c23, c14]*ColorDelta[c24, c21])] + Condensate[{f13, f11}]*Condensate[{f21, f24}]*Condensate[{f22, f12}]*contract[tr[dot[cqqm[0, f24], cc[dot[v11, cqqm[x, f12], hv21]], cqqm[0, f11], v12, cc[dot[hv22, xprop[-x, f23]]]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c21])] + Condensate[{f13, f23}]*Condensate[{f14, f11}]*Condensate[{f22, f24}]*contract[tr[dot[cqqm[0, f11], cc[dot[hv22, cqqm[-x, f23], v12]], cqqm[0, f24], hv21, cc[dot[v11, xprop[x, f12]]]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c22])] + Condensate[{f14, f11}]*Condensate[{f21, f12}]*Condensate[{f22, f24}]*contract[tr[dot[cqqm[0, f11], cc[dot[hv22, xprop[-x, f23], v12]], cqqm[0, f24], hv21, cc[dot[v11, cqqm[x, f12]]]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c22])];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia1=FourierXP[Expand[dia1]/.{Power[_quarkMass,n_]->0, mm1_quarkMass mm2_quarkMass->0},{x,q}];

dia1=QEvaluate[I ScaleMu^(4-D) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,cc,tr,dot,lor,sun,sunSimplify,contract},


dia2=Condensate[{f13, f24}]*Condensate[{f14, f23}]*Condensate[{f22, f12}]*contract[tr[dot[xprop[x, f11], cc[dot[v11, cqqm[x, f12], hv21]]]]*tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, cqqm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f13, f11}]*Condensate[{f21, f12}]*Condensate[{f22, f23}]*contract[tr[dot[cc[dot[cqqm[0, f11], v12]], v11, cqqm[x, f12], cc[dot[hv22, cqqm[0, f23], hv21]], xprop[-x, f24]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c21]*ColorDelta[c23, c22]*ColorDelta[c24, c14])] + Condensate[{f14, f12}]*Condensate[{f21, f24}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, cqqm[0, f12], cc[dot[hv22, xprop[-x, f23], v12]], cqqm[0, f24]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c14]*ColorDelta[c23, c13]*ColorDelta[c24, c21])] + Condensate[{f14, f11}]*Condensate[{f21, f24}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[hv22, xprop[-x, f23], v12]], cqqm[0, f24], cc[dot[v11, cqqm[x, f12], hv21]], cqqm[0, f11]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c21])] + Condensate[{f14, f12}]*Condensate[{f21, f11}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[cqqm[x, f11]]], v11, cqqm[0, f12], cc[dot[hv22, xprop[-x, f23], v12]], cqqm[0, f24], hv21]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c14]*ColorDelta[c23, c13]*ColorDelta[c24, c22])] + Condensate[{f13, f23}]*Condensate[{f14, f12}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[xprop[x, f11]]], v11, cqqm[0, f12], cc[dot[hv22, cqqm[-x, f23], v12]], cqqm[0, f24], hv21]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c14]*ColorDelta[c23, c13]*ColorDelta[c24, c22])] + Condensate[{f13, f11}]*Condensate[{f14, f23}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[cqqm[0, f24], hv21]], hv22, cqqm[-x, f23], cc[dot[cqqm[0, f11], v12]], v11, xprop[x, f12]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c22])] + Condensate[{f13, f11}]*Condensate[{f21, f12}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[cqqm[0, f24], hv21]], hv22, xprop[-x, f23], cc[dot[cqqm[0, f11], v12]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c22])];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia2=FourierXP[Expand[dia2]/.{Power[_quarkMass,n_]->0, mm1_quarkMass mm2_quarkMass->0},{x,q}];

dia2=QEvaluate[I ScaleMu^(4-D) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,cc,tr,dot,lor,sun,sunSimplify,contract},


dia3=Condensate[{f13, f24}]*Condensate[{f14, f23}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, cqqm[-x, f23]]]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, xprop[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f14, f23}]*Condensate[{f21, f12}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, cqqm[x, f12]]]*tr[dot[cc[dot[xprop[-x, f24], v12]], hv22, cqqm[-x, f23]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f13, f24}]*Condensate[{f21, f11}]*Condensate[{f22, f12}]*contract[tr[dot[cqqm[x, f11], cc[dot[v11, cqqm[x, f12], hv21]]]]*tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, xprop[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f13, f23}]*Condensate[{f14, f24}]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], cqqm[-x, f24]]]*tr[dot[cc[dot[xprop[x, f11], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate[{f13, f23}]*Condensate[{f14, f24}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], cqqm[-x, f24]]]*tr[dot[xprop[x, f11], cc[dot[v11, cqqm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate[{f13, f11}]*Condensate[{f14, f24}]*Condensate[{f21, f23}]*contract[tr[dot[cc[dot[cqqm[-x, f24]]], hv22, cqqm[0, f23], cc[dot[v11, xprop[x, f12], hv21]], cqqm[0, f11], v12]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c22]*ColorDelta[c23, c21]*ColorDelta[c24, c14])] + Condensate[{f13, f11}]*Condensate[{f21, f23}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[xprop[-x, f24]]], hv22, cqqm[0, f23], cc[dot[v11, cqqm[x, f12], hv21]], cqqm[0, f11], v12]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c22]*ColorDelta[c23, c21]*ColorDelta[c24, c14])] + Condensate[{f13, f12}]*Condensate[{f21, f11}]*Condensate[{f22, f23}]*contract[tr[dot[cc[dot[hv22, cqqm[0, f23], hv21]], xprop[-x, f24], cc[dot[v11, cqqm[0, f12], v12]], cqqm[x, f11]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c13]*ColorDelta[c23, c22]*ColorDelta[c24, c14])];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia3=FourierXP[Expand[dia3]/.{Power[_quarkMass,n_]->0, mm1_quarkMass mm2_quarkMass->0},{x,q}];

dia3=QEvaluate[I ScaleMu^(4-D) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,cc,tr,dot,lor,sun,sunSimplify,contract},


dia4=Condensate[{f13, f24}]*Condensate[{f21, f12}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, xprop[-x, f23]]]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f13, f23}]*Condensate[{f21, f11}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], xprop[-x, f24]]]*tr[dot[cqqm[x, f11], cc[dot[v11, cqqm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate[{f14, f24}]*Condensate[{f21, f11}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[hv22, xprop[-x, f23], v12]], cqqm[-x, f24]]]*tr[dot[cqqm[x, f11], cc[dot[v11, cqqm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate[{f13, f23}]*Condensate[{f14, f24}]*Condensate[{f21, f11}]*contract[tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], cqqm[-x, f24]]]*tr[dot[cqqm[x, f11], cc[dot[v11, xprop[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate[{f13, f12}]*Condensate[{f21, f23}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv22, cqqm[0, f23]]], xprop[-x, f24], cc[dot[v11, cqqm[0, f12], v12]], cqqm[x, f11], hv21]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c13]*ColorDelta[c23, c21]*ColorDelta[c24, c14])] + Condensate[{f13, f23}]*Condensate[{f14, f12}]*Condensate[{f21, f24}]*contract[tr[dot[cc[dot[xprop[x, f11], hv21]], v11, cqqm[0, f12], cc[dot[hv22, cqqm[-x, f23], v12]], cqqm[0, f24]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c14]*ColorDelta[c23, c13]*ColorDelta[c24, c21])] + Condensate[{f13, f12}]*Condensate[{f21, f24}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[cqqm[0, f24]]], hv22, xprop[-x, f23], cc[dot[v11, cqqm[0, f12], v12]], cqqm[x, f11], hv21]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c13]*ColorDelta[c23, c14]*ColorDelta[c24, c21])] + Condensate[{f13, f11}]*Condensate[{f14, f23}]*Condensate[{f21, f24}]*contract[tr[dot[cqqm[0, f24], cc[dot[v11, xprop[x, f12], hv21]], cqqm[0, f11], v12, cc[dot[hv22, cqqm[-x, f23]]]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c21])];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia4=FourierXP[Expand[dia4]/.{Power[_quarkMass,n_]->0, mm1_quarkMass mm2_quarkMass->0},{x,q}];

dia4=QEvaluate[I ScaleMu^(4-D) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia5,cc,tr,dot,lor,sun,sunSimplify,contract},


dia5=Condensate[{f14, f12}]*Condensate[{f21, f23}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, cqqm[0, f12], cc[dot[xprop[-x, f24], v12]], hv22, cqqm[0, f23]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c14]*ColorDelta[c23, c21]*ColorDelta[c24, c13])] + Condensate[{f13, f24}]*Condensate[{f14, f11}]*Condensate[{f21, f23}]*contract[tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, cqqm[0, f23], cc[dot[v11, xprop[x, f12], hv21]], cqqm[0, f11]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c22]*ColorDelta[c23, c21]*ColorDelta[c24, c13])] + Condensate[{f13, f24}]*Condensate[{f14, f11}]*Condensate[{f22, f23}]*contract[tr[dot[cc[dot[cqqm[0, f11]]], v11, xprop[x, f12], cc[dot[hv22, cqqm[0, f23], hv21]], cqqm[-x, f24], v12]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c21]*ColorDelta[c23, c22]*ColorDelta[c24, c13])] + Condensate[{f13, f12}]*Condensate[{f14, f24}]*Condensate[{f21, f23}]*contract[tr[dot[cc[dot[hv22, cqqm[0, f23]]], cqqm[-x, f24], cc[dot[v11, cqqm[0, f12], v12]], xprop[x, f11], hv21]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c12, c13]*ColorDelta[c23, c21]*ColorDelta[c24, c14])] + Condensate[{f13, f12}]*Condensate[{f14, f24}]*Condensate[{f22, f23}]*contract[tr[dot[cc[dot[hv22, cqqm[0, f23], hv21]], cqqm[-x, f24], cc[dot[v11, cqqm[0, f12], v12]], xprop[x, f11]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c13]*ColorDelta[c23, c22]*ColorDelta[c24, c14])] + Condensate[{f13, f11}]*Condensate[{f14, f24}]*Condensate[{f22, f23}]*contract[tr[dot[cc[dot[cqqm[0, f11], v12]], v11, xprop[x, f12], cc[dot[hv22, cqqm[0, f23], hv21]], cqqm[-x, f24]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c13]*ColorDelta[c12, c21]*ColorDelta[c23, c22]*ColorDelta[c24, c14])] + Condensate[{f13, f23}]*Condensate[{f14, f11}]*Condensate[{f21, f24}]*contract[tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], cqqm[0, f24], cc[dot[v11, xprop[x, f12], hv21]], cqqm[0, f11]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c21])] + Condensate[{f13, f12}]*Condensate[{f14, f23}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[cqqm[0, f24], hv21]], hv22, cqqm[-x, f23], cc[dot[v11, cqqm[0, f12], v12]], xprop[x, f11]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c13]*ColorDelta[c23, c14]*ColorDelta[c24, c22])];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia5=FourierXP[Expand[dia5]/.{Power[_quarkMass,n_]->0, mm1_quarkMass mm2_quarkMass->0},{x,q}];

dia5=QEvaluate[I ScaleMu^(4-D) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization::Plain:: *)
xtype6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia6,cc,tr,dot,lor,sun,sunSimplify,contract},


dia6=Condensate[{f13, f24}]*Condensate[{f14, f23}]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, cqqm[-x, f23]]]*tr[dot[cc[dot[xprop[x, f11], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f14, f23}]*Condensate[{f21, f11}]*Condensate[{f22, f12}]*contract[tr[dot[cqqm[x, f11], cc[dot[v11, cqqm[x, f12], hv21]]]]*tr[dot[cc[dot[xprop[-x, f24], v12]], hv22, cqqm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f14, f11}]*Condensate[{f21, f23}]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[xprop[-x, f24], v12]], hv22, cqqm[0, f23], cc[dot[v11, cqqm[x, f12], hv21]], cqqm[0, f11]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c22]*ColorDelta[c23, c21]*ColorDelta[c24, c13])] + Condensate[{f14, f11}]*Condensate[{f21, f12}]*Condensate[{f22, f23}]*contract[tr[dot[cc[dot[cqqm[0, f11]]], v11, cqqm[x, f12], cc[dot[hv22, cqqm[0, f23], hv21]], xprop[-x, f24], v12]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c14]*ColorDelta[c12, c21]*ColorDelta[c23, c22]*ColorDelta[c24, c13])] + Condensate[{f13, f23}]*Condensate[{f21, f12}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], xprop[-x, f24]]]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate[{f14, f24}]*Condensate[{f21, f12}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv22, xprop[-x, f23], v12]], cqqm[-x, f24]]]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate[{f13, f23}]*Condensate[{f14, f24}]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], cqqm[-x, f24]]]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, xprop[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate[{f13, f12}]*Condensate[{f21, f11}]*Condensate[{f22, f24}]*contract[tr[dot[cc[dot[cqqm[0, f24], hv21]], hv22, xprop[-x, f23], cc[dot[v11, cqqm[0, f12], v12]], cqqm[x, f11]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c12, c13]*ColorDelta[c23, c14]*ColorDelta[c24, c22])];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia6=FourierXP[Expand[dia6]/.{Power[_quarkMass,n_]->0, mm1_quarkMass mm2_quarkMass->0},{x,q}];

dia6=QEvaluate[I ScaleMu^(4-D) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
