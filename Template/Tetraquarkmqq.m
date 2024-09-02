(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


Tetraquarkmqq::usage="Tetraquarkmqq[q_,j1_,j2_] gives m\[InvisibleComma]\[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)q\[RightAngleBracket] contribution of <j1 j2>, with j1, j2 are tetraquark currents. "

Tetraquarkmqq::inderr="Dummy indices conflict!"

Tetraquarkmqq::curerr="Unknow structure of the current!"


Begin["`Private`Tetraquarkmqq`"]

Options[Tetraquarkmqq] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True
}




(*------------------------------------------------------------------*)

xprop[x_] = 1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];


xprop[x_,a_] = (Pair[Momentum[x,D],Momentum[x,D]]^(1-D/2) quarkMass[a]qfact1[1/4 I^(2-D) \[Pi]^(-D/2) qGamma[-1+D/2]] 
				+ DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qfact1[1/2 I^(1-D) \[Pi]^(-D/2) qGamma[D/2]]);
			
			
(*--- Condensate ---*)

mqq[f1_,f2_,x_]:=Condensate[{f1,f2}](-1/(4CA)+I/(4D CA) quarkMass[f1]GSD[x]);
(* \[LeftAngleBracket]q(x)Overscript[q, _](0)\[RightAngleBracket] = \[LeftAngleBracket]Overscript[q, _]q\[RightAngleBracket](-1/(4CA) + I/(4 D CA)m GSD[x]) *)
(* \[LeftAngleBracket]q(0)Overscript[q, _](x)\[RightAngleBracket] = \[LeftAngleBracket]Overscript[q, _]q\[RightAngleBracket](-1/(4CA) + I/(4 D CA)m GSD[-x]) *)			



cqqm[x_,f_]=-1/4*1/CA + ((I/4)*DiracGamma[Momentum[x, D], D]*quarkMass[f])/(CA*D);

(*-----------------------------------*)
(* propagators and condensates *)

xpropm[x_,f_]=(DiracGamma[Momentum[x, D], D]*qfact1[(I^(1 - D)*qGamma[D/2])/(2*Pi^(D/2))])/Pair[Momentum[x, D], Momentum[x, D]]^(D/2) + 
 Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[(I^(2 - D)*qGamma[-1 + D/2]*quarkMass[f])/(4*Pi^(D/2))];


fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]


(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
Tetraquarkmqq[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 Tetraquarkmqq[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])


Tetraquarkmqq[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Module[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,trs1,trs2,trs3,trs4,dia,null0,x,q,holdf=OptionValue[HoldFlavor],ndr=OptionValue[AutoNDR],
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,atr,fdir,files,pall,diagrams},



(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[Tetraquarkmqq::curerr];
	Abort[]
];


(*--- check dummy indices in input ---*)
tmp2=Length[#]&/@Gather[({c11,c12,c13,c14,c21,c22,c23,c24}//FCI)/.SUNFIndex[aa_]:>aa];

If[DummyCheck[v11,v12,v21,v22]==0||!FreeQ[tmp2,nn_Integer/;nn>2],
	Message[Tetraquarkmqq::inderr];
	Abort[]
];



(*-------------------------------*)
(* B A^+ B *)
hv21=v21//ComplexConjugate;
hv22=v22//ComplexConjugate;


If[OptionValue[AutoNDR]===True,atr=TR5,atr=TR];


(*---------------------------------------------------*)
If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];


diagrams={xtype1,xtype2,xtype3,xtype4};


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
			files=("Tetraquarkmqq_"<>#)&/@files;
			
			
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
(*------------------------------------------------------------------*)




(* ::Code::Initialization::Plain:: *)
(*-------------The following are generated by algorithem----------------*)


(* ::Input::Initialization:: *)
xtype1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,cc,tr,dot,sunSimplify,contract},


 dia1=Condensate[{f13, f24}]*contract[tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, xpropm[-x, f23]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f14, f23}]*contract[tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, cqqm[-x, f23]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f21, f11}]*contract[tr[dot[cqqm[x, f11], cc[dot[v11, xpropm[x, f12], hv21]]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f21, f12}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xpropm[-x, f24]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]];


 dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia1=FourierXP[dia1,{x,q}];


 dia1=dia1/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia1=QEvaluate[I ScaleMu^(2(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia1=QGather[Expand[dia1]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,cc,tr,dot,sunSimplify,contract},


 dia2=Condensate[{f22, f12}]*contract[tr[dot[xpropm[x, f11], cc[dot[v11, cqqm[x, f12], hv21]]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f13, f23}]*contract[tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], xpropm[-x, f24]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate[{f21, f11}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xpropm[-x, f24]]]*tr[dot[cqqm[x, f11], cc[dot[v11, xpropm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate[{f14, f24}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], cqqm[-x, f24]]]*tr[dot[xpropm[x, f11], cc[dot[v11, xpropm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]];


 dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia2=FourierXP[dia2,{x,q}];


 dia2=dia2/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia2=QEvaluate[I ScaleMu^(2(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia2=QGather[Expand[dia2]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,cc,tr,dot,sunSimplify,contract},


 dia3=Condensate[{f14, f23}]*contract[tr[dot[xpropm[x, f11], cc[dot[v11, xpropm[x, f12], hv21]]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, cqqm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f14, f24}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], cqqm[-x, f24]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate[{f22, f12}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xpropm[-x, f24]]]*tr[dot[xpropm[x, f11], cc[dot[v11, cqqm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate[{f13, f23}]*contract[tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], xpropm[-x, f24]]]*tr[dot[xpropm[x, f11], cc[dot[v11, xpropm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]];


 dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia3=FourierXP[dia3,{x,q}];


 dia3=dia3/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia3=QEvaluate[I ScaleMu^(2(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia3=QGather[Expand[dia3]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,cc,tr,dot,sunSimplify,contract},


 dia4=Condensate[{f22, f11}]*contract[tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, xpropm[x, f12]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f21, f12}]*contract[tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xpropm[-x, f23]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f13, f24}]*contract[tr[dot[xpropm[x, f11], cc[dot[v11, xpropm[x, f12], hv21]]]]*tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xpropm[-x, f24]]]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]];


 dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia4=FourierXP[dia4,{x,q}];


 dia4=dia4/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia4=QEvaluate[I ScaleMu^(2(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia4=QGather[Expand[dia4]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(*type1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia,null},

(*------------------------------------------------------------------*)
(* figure_1 *)
dia=(FlavorDelta[f11,f21]FlavorDelta[f12,f22]ColorDelta[c11,c21]ColorDelta[c12,c22]atr[v11.xprop[x,f12].hv21.fcc[mqq[f11,f21,x]]]
	-FlavorDelta[f11,f22]FlavorDelta[f12,f21]ColorDelta[c11,c22]ColorDelta[c12,c21]atr[v11.xprop[x,f12].fcc[mqq[f11,f22,x].hv21]])
	(FlavorDelta[f13,f23]FlavorDelta[f14,f24]ColorDelta[c13,c23]ColorDelta[c14,c24]atr[hv22.xprop[-x,f13].v12.fcc[xprop[-x,f14]]]
	-FlavorDelta[f13,f24]FlavorDelta[f14,f23]ColorDelta[c13,c24]ColorDelta[c14,c23]atr[xprop[-x,f13].v12.fcc[hv22.xprop[-x,f14]]]);

	
(*--------------------*)	
dia=Expand[null dia]/.quarkMass[aa_]:>null quarkMass[aa];
dia=dia/.Power[null,pow_/;pow>2]->0;(* only keep the term linear with quarkmass *)

dia=FourierXP[dia/.null->1,{x,q}];
dia=QEvaluate[I ScaleMu^(2(4-D)) dia, q, HoldFlavor->holdf,Parallelized->False]/.{q->qq,SUNN->CA}
]
*)


(*type2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia,null},

(*------------------------------------------------------------------*)
(* figure_2 *)
dia=(FlavorDelta[f11,f21]FlavorDelta[f12,f22]ColorDelta[c11,c21]ColorDelta[c12,c22]atr[v11.mqq[f12,f22,x].hv21.fcc[xprop[x,f11]]]
	-FlavorDelta[f11,f22]FlavorDelta[f12,f21]ColorDelta[c11,c22]ColorDelta[c12,c21]atr[v11.mqq[f12,f21,x].fcc[xprop[x,f11].hv21]])
	(FlavorDelta[f13,f23]FlavorDelta[f14,f24]ColorDelta[c13,c23]ColorDelta[c14,c24]atr[hv22.xprop[-x,f13].v12.fcc[xprop[-x,f14]]]
	-FlavorDelta[f13,f24]FlavorDelta[f14,f23]ColorDelta[c13,c24]ColorDelta[c14,c23]atr[xprop[-x,f13].v12.fcc[hv22.xprop[-x,f14]]]);


	
(*--------------------*)	
dia=Expand[null dia]/.quarkMass[aa_]:>null quarkMass[aa];
dia=dia/.Power[null,pow_/;pow>2]->0;(* only keep the term linear with quarkmass *)

dia=FourierXP[dia/.null->1,{x,q}];
dia=QEvaluate[I ScaleMu^(2(4-D)) dia, q, HoldFlavor->holdf,Parallelized->False]/.{q->qq,SUNN->CA}
]
*)


(*type3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia,null},

(*------------------------------------------------------------------*)
(* figure_3 *)
dia=(FlavorDelta[f11,f21]FlavorDelta[f12,f22]ColorDelta[c11,c21]ColorDelta[c12,c22]atr[v11.xprop[x,f12].hv21.fcc[xprop[x,f11]]]
	-FlavorDelta[f11,f22]FlavorDelta[f12,f21]ColorDelta[c11,c22]ColorDelta[c12,c21]atr[v11.xprop[x,f12].fcc[xprop[x,f11].hv21]])
	(FlavorDelta[f13,f23]FlavorDelta[f14,f24]ColorDelta[c13,c23]ColorDelta[c14,c24]atr[hv22.mqq[f13,f23,-x].v12.fcc[xprop[-x,f14]]]
	-FlavorDelta[f13,f24]FlavorDelta[f14,f23]ColorDelta[c13,c24]ColorDelta[c14,c23]atr[mqq[f13,f24,-x].v12.fcc[hv22.xprop[-x,f14]]]);


	
(*--------------------*)	
dia=Expand[null dia]/.quarkMass[aa_]:>null quarkMass[aa];
dia=dia/.Power[null,pow_/;pow>2]->0;(* only keep the term linear with quarkmass *)

dia=FourierXP[dia/.null->1,{x,q}];
dia=QEvaluate[I ScaleMu^(2(4-D)) dia, q, HoldFlavor->holdf,Parallelized->False]/.{q->qq,SUNN->CA}
]

*)


(*type4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia,null},

(*------------------------------------------------------------------*)
(* figure_4 *)
dia=(FlavorDelta[f11,f21]FlavorDelta[f12,f22]ColorDelta[c11,c21]ColorDelta[c12,c22]atr[v11.xprop[x,f12].hv21.fcc[xprop[x,f11]]]
	-FlavorDelta[f11,f22]FlavorDelta[f12,f21]ColorDelta[c11,c22]ColorDelta[c12,c21]atr[v11.xprop[x,f12].fcc[xprop[x,f11].hv21]])
	(FlavorDelta[f13,f23]FlavorDelta[f14,f24]ColorDelta[c13,c23]ColorDelta[c14,c24]atr[hv22.xprop[-x,f13].v12.fcc[mqq[f14,f24,-x]]]
	-FlavorDelta[f13,f24]FlavorDelta[f14,f23]ColorDelta[c13,c24]ColorDelta[c14,c23]atr[xprop[-x,f13].v12.fcc[hv22.mqq[f14,f23,-x]]]);
	

	
(*--------------------*)	
dia=Expand[null dia]/.quarkMass[aa_]:>null quarkMass[aa];
dia=dia/.Power[null,pow_/;pow>2]->0;(* only keep the term linear with quarkmass *)

dia=FourierXP[dia/.null->1,{x,q}];
dia=QEvaluate[I ScaleMu^(2(4-D)) dia, q, HoldFlavor->holdf,Parallelized->False]/.{q->qq,SUNN->CA}
]*)



(*------------------*)


(*
(*------------------------------------------------------------------*)
(* figure *)
trs1=(FlavorDelta[f11,f21]FlavorDelta[f12,f22]ColorDelta[c11,c21]ColorDelta[c12,c22]atr[v11.xprop[x,f12].hv21.fcc[mqq[f11,f21,x]]]
	-FlavorDelta[f11,f22]FlavorDelta[f12,f21]ColorDelta[c11,c22]ColorDelta[c12,c21]atr[v11.xprop[x,f12].fcc[mqq[f11,f22,x].hv21]])
	(FlavorDelta[f13,f23]FlavorDelta[f14,f24]ColorDelta[c13,c23]ColorDelta[c14,c24]atr[hv22.xprop[-x,f13].v12.fcc[xprop[-x,f14]]]
	-FlavorDelta[f13,f24]FlavorDelta[f14,f23]ColorDelta[c13,c24]ColorDelta[c14,c23]atr[xprop[-x,f13].v12.fcc[hv22.xprop[-x,f14]]]);

trs2=(FlavorDelta[f11,f21]FlavorDelta[f12,f22]ColorDelta[c11,c21]ColorDelta[c12,c22]atr[v11.mqq[f12,f22,x].hv21.fcc[xprop[x,f11]]]
	-FlavorDelta[f11,f22]FlavorDelta[f12,f21]ColorDelta[c11,c22]ColorDelta[c12,c21]atr[v11.mqq[f12,f21,x].fcc[xprop[x,f11].hv21]])
	(FlavorDelta[f13,f23]FlavorDelta[f14,f24]ColorDelta[c13,c23]ColorDelta[c14,c24]atr[hv22.xprop[-x,f13].v12.fcc[xprop[-x,f14]]]
	-FlavorDelta[f13,f24]FlavorDelta[f14,f23]ColorDelta[c13,c24]ColorDelta[c14,c23]atr[xprop[-x,f13].v12.fcc[hv22.xprop[-x,f14]]]);

trs3=(FlavorDelta[f11,f21]FlavorDelta[f12,f22]ColorDelta[c11,c21]ColorDelta[c12,c22]atr[v11.xprop[x,f12].hv21.fcc[xprop[x,f11]]]
	-FlavorDelta[f11,f22]FlavorDelta[f12,f21]ColorDelta[c11,c22]ColorDelta[c12,c21]atr[v11.xprop[x,f12].fcc[xprop[x,f11].hv21]])
	(FlavorDelta[f13,f23]FlavorDelta[f14,f24]ColorDelta[c13,c23]ColorDelta[c14,c24]atr[hv22.mqq[f13,f23,-x].v12.fcc[xprop[-x,f14]]]
	-FlavorDelta[f13,f24]FlavorDelta[f14,f23]ColorDelta[c13,c24]ColorDelta[c14,c23]atr[mqq[f13,f24,-x].v12.fcc[hv22.xprop[-x,f14]]]);

trs4=(FlavorDelta[f11,f21]FlavorDelta[f12,f22]ColorDelta[c11,c21]ColorDelta[c12,c22]atr[v11.xprop[x,f12].hv21.fcc[xprop[x,f11]]]
	-FlavorDelta[f11,f22]FlavorDelta[f12,f21]ColorDelta[c11,c22]ColorDelta[c12,c21]atr[v11.xprop[x,f12].fcc[xprop[x,f11].hv21]])
	(FlavorDelta[f13,f23]FlavorDelta[f14,f24]ColorDelta[c13,c23]ColorDelta[c14,c24]atr[hv22.xprop[-x,f13].v12.fcc[mqq[f14,f24,-x]]]
	-FlavorDelta[f13,f24]FlavorDelta[f14,f23]ColorDelta[c13,c24]ColorDelta[c14,c23]atr[xprop[-x,f13].v12.fcc[hv22.mqq[f14,f23,-x]]]);
	
	
(*--------------------*)	
dia=Expand[null0(trs1+trs2+trs3+trs4)]/.quarkMass[aa_]:>null0 quarkMass[aa];
dia=dia/.Power[null0,pow_/;pow>2]->0;(* only keep the term linear with quarkmass *)

dia=FourierXP[dia/.null0->1,{x,q}];
dia=QEvaluate[I ScaleMu^(2(4-D)) dia, q, HoldFlavor->holdf,Parallelized->False];


(*------------------------------------------------------------------*)
QGather[dia//SUNSimplify,q,ShowasTable->False]/.{q->qq,SUNN->CA}*)


End[]
(*EndPackage[]*)
