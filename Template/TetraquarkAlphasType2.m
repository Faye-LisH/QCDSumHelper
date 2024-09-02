(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
TetraquarkAlphasType2::usage="TetraquarkAlphasType2[q_,j1_,j2_] give the next leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents.
 Here the gluon propagator connect with j1 & j2, the corresponding renormalization invole hybird-like current."

TetraquarkAlphasType2::inderr="Dummy indices conflict!"

TetraquarkAlphasType2::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`TetraquarkAlphasType2`"]


(* ::Code::Initialization::Plain:: *)
Options[TetraquarkAlphasType2] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True,
	Renormalization->True,
	Strategy->"Fourier"
}


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(* propagators *)

prop[q_]=I GSD[q]FAD[q];

(*xprop[x_]=FourierPX[prop[q],{q,x}];*)
xprop[x_]=1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];
		
			
xprop2[x_,lor1_,lor2_,lor3_,lor4_]=Pair[LorentzIndex[lor1, D], LorentzIndex[lor2, D]]*Pair[LorentzIndex[lor3, D], LorentzIndex[lor4, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(3 - (3*D)/2)*
  qfact1[(qfact2[((20 - 8*D + D^2)*gStrong^2)/(1024*E^(((3*I)/2)*D*Pi)*Pi^((3*D)/2))]*qGamma[1 - D/2]^2*qGamma[D/2]^4*qGamma[-3 + (3*D)/2])/(qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lor4, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lor3, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(3 - (3*D)/2)*
  qfact1[(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-3 + (3*D)/2])/(256*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lor3, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lor4, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(3 - (3*D)/2)*
  qfact1[(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-3 + (3*D)/2])/(256*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lor3, D], LorentzIndex[lor4, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[(qfact2[((12 - 8*D + D^2)*gStrong^2)/(512*E^(((3*I)/2)*D*Pi)*Pi^((3*D)/2))]*qGamma[1 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lor2, D]]*Pair[LorentzIndex[lor3, D], Momentum[x, D]]*Pair[LorentzIndex[lor4, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[(qfact2[((12 - 8*D + D^2)*gStrong^2)/(512*E^(((3*I)/2)*D*Pi)*Pi^((3*D)/2))]*qGamma[1 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lor4, D]]*Pair[LorentzIndex[lor3, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[-(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(128*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lor4, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lor3, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[-(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(128*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lor3, D]]*Pair[LorentzIndex[lor4, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[-(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(128*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lor3, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lor4, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[-(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(128*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lor3, D], Momentum[x, D]]*Pair[LorentzIndex[lor4, D], Momentum[x, D]]*
  Pair[Momentum[x, D], Momentum[x, D]]^(1 - (3*D)/2)*qfact1[(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-1 + (3*D)/2])/(64*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)];


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
TetraquarkAlphasType2[qq_,factor1_ current1_FourquarkCurrent, current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 TetraquarkAlphasType2[qq,current1,current2,ops]/;FreeQ[factor1,Current]
TetraquarkAlphasType2[qq_,current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor2 TetraquarkAlphasType2[qq,current1,current2,ops]/;FreeQ[factor2,Current]
TetraquarkAlphasType2[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 TetraquarkAlphasType2[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])


TetraquarkAlphasType2[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,atr,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,listpole1,listpole2,holdf=OptionValue[HoldFlavor],strategy=OptionValue[Strategy],diagrams,
ren=OptionValue[Renormalization],fdir,files1,pall1,files2,pall2,pall,ndr=OptionValue[AutoNDR]},


(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[TetraquarkAlphasType2::curerr];
	Abort[]
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




If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];

(*---------------------------------------------------*)


If[ToLowerCase[ToString[strategy]]=="fourier",
	diagrams={xtype1,xtype2,xtype3,xtype4,xtype5,xtype6,xtype7,xtype8}
,
	If[ToLowerCase[ToString[strategy]]=="integratefirst",
		diagrams={iptype1,iptype2,iptype3,iptype4,iptype5,iptype6,iptype7,iptype8}
	,
		diagrams={ptype1,ptype2,ptype3,ptype4,ptype5,ptype6,ptype7,ptype8}
	]
];


(*----------------------------------------------*)

If[OptionValue[AutoNDR]===True,atr=TR5,atr=TR];


(*----------------------------------------------*)
If[pall===True,

	DistributeDefinitions[qq,current1,current2];
	tmp= Plus@@WaitAll[Join[ParallelSubmit[{f11,f12,f13,f14,f21,f22,f23,f24,ndr,holdf,ren},-#[qq,current1,current2,ndr,holdf,ren]]&/@{reno1,reno2}, ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr]]&/@diagrams]];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,
	If[pall==="External",

	DistributeDefinitions[qq,current1,current2];
	
		If[fdir==="None",
			(* evaluation, no import and export *)
			Join[ParallelSubmit[{f11,f12,f13,f14,f21,f22,f23,f24,ndr,holdf,ren},-#[qq,current1,current2,ndr,holdf,ren]]&/@{reno1,reno2}, ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr]]&/@diagrams]
		
		,
		(* evaluation, import and export *)
			files1=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files1=("TetraquarkAlphasType2_"<>#)&/@files1;
			
			
			pall1=ImExport[fdir,
						files1,
						{{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
						{qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr},
						diagrams}
						];
		
		
			files2=(StringSplit[ToString[#],"`"][[-1]])&/@{reno1,reno2};
			files2=("TetraquarkAlphasType2_"<>#)&/@files2;
			
			If[ren===True,
				pall2=ImExport[fdir,
						files2,
						{{ndr,holdf,ren},
						{qq,current1,current2,ndr,holdf,ren},
						{reno1,reno2}}
						]
			,
				pall2={{2,2},{0,0}}			
			];


			{pall1[[1]]+pall2[[1]],Join[pall1[[2]],-pall2[[2]]]}		
				
		](* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
		
	
	,
	
		
		tmp= Plus@@Join[-#[qq,current1,current2,ndr,holdf,ren]&/@{reno1,reno2},#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr]&/@diagrams];
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]

]

(*------------------------------------------------------------------*)



(* ::Code::Initialization::Plain:: *)
(*	
reno[qq_,cur1_,cur2_,ndr_,holdf_,ndr_,holdf_,False]=0

		
reno[qq_,cur1_,cur2_,ndr_,holdf_,ndr_,holdf_,True]:=Block[{listpole1,listpole2},

listpole1=Plus@@(TetraquarkHybridlike1[qq,#[[2]],cur2,Pole->#[[1]],EpsOrder->0,AutoNDR->ndr,HoldFlavor->holdf,Parallelized->False]&/@FourquarkHPole[cur1]);
listpole2=Plus@@(TetraquarkHybridlike2[qq,cur1,#[[2]],Pole->#[[1]],EpsOrder->0,AutoNDR->ndr,HoldFlavor->holdf,Parallelized->False]&/@FourquarkHPole[cur2]);
	
QGather[listpole1+listpole2,qq,ShowasTable->False]

]
*)


(* ::Code::Initialization::Plain:: *)
reno1[qq_,cur1_,cur2_,ndr_,holdf_,False]=0

reno1[qq_,cur1_,cur2_,ndr_,holdf_,True]:=Block[{listpole1,listpole2},

Plus@@(TetraquarkHybridlike1[qq,#[[2]],cur2,Pole->#[[1]],EpsOrder->0,AutoNDR->ndr,HoldFlavor->holdf,Parallelized->False]&/@FourquarkHPole[cur1])

]


reno2[qq_,cur1_,cur2_,ndr_,holdf_,False]=0

reno2[qq_,cur1_,cur2_,ndr_,holdf_,True]:=Block[{listpole1,listpole2},

Plus@@(TetraquarkHybridlike2[qq,cur1,#[[2]],Pole->#[[1]],EpsOrder->0,AutoNDR->ndr,HoldFlavor->holdf,Parallelized->False]&/@FourquarkHPole[cur2])

]


(* ::Code::Initialization::Plain:: *)
(* The following are generated by algorithem *)


(* ::Input::Initialization:: *)
xtype1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,lor1,lor2,lor3,lor4,cc,tr,dot,sunSimplify,contract},


dia1=contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[GAD[lor3], GAD[lorz], GAD[lor4], hv21]], hv22, xprop[-x], cc[dot[GAD[lor1], GAD[lorz], GAD[lor2], v12]], v11, xprop[x]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c13]*SUNTF[sun, c24, c22])] + contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[GAD[lor3], GAD[lorz], GAD[lor4], hv21]], hv22, xprop[-x], cc[dot[v11, GAD[lor1], GAD[lorz], GAD[lor2], v12]], xprop[x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c13]*SUNTF[sun, c24, c22])];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia1=FourierXP[dia1,{x,q}];

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
xtype2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,lor1,lor2,lor3,lor4,cc,tr,dot,sunSimplify,contract},


dia2=contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[hv22, xprop[-x], v12]], GAD[lor3], GAD[lorz], GAD[lor4], cc[dot[v11, xprop[x], hv21]], GAD[lor1], GAD[lorz], GAD[lor2]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c14]*SUNTF[sun, c24, c21])] + contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[GAD[lor3], GAD[lorz], GAD[lor4]]], hv22, xprop[-x], cc[dot[v11, GAD[lor1], GAD[lorz], GAD[lor2], v12]], xprop[x], hv21]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c13]*SUNTF[sun, c24, c21])];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia2=FourierXP[dia2,{x,q}];

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
xtype3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,lor1,lor2,lor3,lor4,cc,tr,dot,sunSimplify,contract},


dia3=contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], GAD[lor3], GAD[lorz], GAD[lor4], cc[dot[v11, xprop[x], hv21]], GAD[lor1], GAD[lorz], GAD[lor2], v12, cc[dot[hv22, xprop[-x]]]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c13]*SUNTF[sun, c24, c21])] + contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[xprop[x]]], v11, GAD[lor1], GAD[lorz], GAD[lor2], cc[dot[hv22, xprop[-x], v12]], GAD[lor3], GAD[lorz], GAD[lor4], hv21]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c14]*SUNTF[sun, c24, c22])];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia3=FourierXP[dia3,{x,q}];

dia3=QEvaluate[I ScaleMu^(4(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
xtype4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,lor1,lor2,lor3,lor4,cc,tr,dot,sunSimplify,contract},


dia4=contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[xprop[-x]]], hv22, GAD[lor3], GAD[lorz], GAD[lor4], cc[dot[v11, xprop[x], hv21]], GAD[lor1], GAD[lorz], GAD[lor2], v12]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c13]*SUNTF[sun, c23, c21])] + contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[hv22, GAD[lor3], GAD[lorz], GAD[lor4], hv21]], xprop[-x], cc[dot[v11, GAD[lor1], GAD[lorz], GAD[lor2], v12]], xprop[x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c13]*SUNTF[sun, c23, c22])];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia4=FourierXP[dia4,{x,q}];

dia4=QEvaluate[I ScaleMu^(4(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
xtype5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia5,lor1,lor2,lor3,lor4,cc,tr,dot,sunSimplify,contract},


dia5=contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[xprop[x], hv21]], v11, GAD[lor1], GAD[lorz], GAD[lor2], cc[dot[xprop[-x], v12]], hv22, GAD[lor3], GAD[lorz], GAD[lor4]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c14]*SUNTF[sun, c23, c21])] + contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[xprop[x], hv21]], v11, GAD[lor1], GAD[lorz], GAD[lor2], cc[dot[hv22, xprop[-x], v12]], GAD[lor3], GAD[lorz], GAD[lor4]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c14]*SUNTF[sun, c24, c21])];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia5=FourierXP[dia5,{x,q}];

dia5=QEvaluate[I ScaleMu^(4(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
xtype6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia6,lor1,lor2,lor3,lor4,cc,tr,dot,sunSimplify,contract},


dia6=contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[GAD[lor1], GAD[lorz], GAD[lor2], v12]], v11, xprop[x], cc[dot[hv22, GAD[lor3], GAD[lorz], GAD[lor4], hv21]], xprop[-x]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c13]*SUNTF[sun, c23, c22])] + contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], GAD[lor1], GAD[lorz], GAD[lor2], cc[dot[hv22, xprop[-x], v12]], GAD[lor3], GAD[lorz], GAD[lor4], hv21, cc[dot[v11, xprop[x]]]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*sunSimplify[-(ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c14]*SUNTF[sun, c24, c22])];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia6=FourierXP[dia6,{x,q}];

dia6=QEvaluate[I ScaleMu^(4(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
xtype7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia7,lor1,lor2,lor3,lor4,cc,tr,dot,sunSimplify,contract},


dia7=contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[xprop[-x], v12]], hv22, GAD[lor3], GAD[lorz], GAD[lor4], cc[dot[v11, xprop[x], hv21]], GAD[lor1], GAD[lorz], GAD[lor2]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c14]*SUNTF[sun, c23, c21])] + contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], xprop[x], cc[dot[hv22, GAD[lor3], GAD[lorz], GAD[lor4], hv21]], xprop[-x], v12, cc[dot[v11, GAD[lor1], GAD[lorz], GAD[lor2]]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c14]*SUNTF[sun, c23, c22])];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia7=FourierXP[dia7,{x,q}];

dia7=QEvaluate[I ScaleMu^(4(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
xtype8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia8,lor1,lor2,lor3,lor4,cc,tr,dot,sunSimplify,contract},


dia8=contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[hv22, GAD[lor3], GAD[lorz], GAD[lor4]]], xprop[-x], cc[dot[v11, GAD[lor1], GAD[lorz], GAD[lor2], v12]], xprop[x], hv21]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*sunSimplify[-(ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c13]*SUNTF[sun, c23, c21])] + contract[tr[dot[xprop2[x, lor1, lor2, lor3, lor4], cc[dot[GAD[lor1], GAD[lorz], GAD[lor2]]], v11, xprop[x], cc[dot[hv22, GAD[lor3], GAD[lorz], GAD[lor4], hv21]], xprop[-x], v12]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*sunSimplify[-(ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c14]*SUNTF[sun, c23, c22])];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];(* SUNSimplify[SUNTrace[SUNT[a,b,c]]SUNF[a,b,c]] doesn't work *)


dia8=FourierXP[dia8,{x,q}];

dia8=QEvaluate[I ScaleMu^(4(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Code::Initialization::Plain:: *)
(*-----------------------------------------------*)


(* ::Input::Initialization:: *)
iptype1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,l,k,dia1,sun,lor1,lor2,lor3,lor4,lorz},


dia1=-(ColorDelta[c11, c22]*ColorDelta[c23, c14]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[GAD[lor3] . GAD[lorz] . GAD[lor4]] . hv22 . prop[-k + l] . fcc[v11 . GAD[lor1] . GAD[lorz] . GAD[lor2] . v12] . prop[-k + q] . hv21]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*SUNTF[sun, c12, c13]*SUNTF[sun, c24, c21]) - ColorDelta[c11, c22]*ColorDelta[c23, c13]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[prop[-k + q] . hv21] . v11 . GAD[lor1] . GAD[lorz] . GAD[lor2] . fcc[hv22 . prop[-k + l] . v12] . GAD[lor3] . GAD[lorz] . GAD[lor4]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*SUNTF[sun, c12, c14]*SUNTF[sun, c24, c21];


dia1=(IntegrateP[dia1,{l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,l,k,dia2,sun,lor1,lor2,lor3,lor4,lorz},


dia2=-(ColorDelta[c11, c22]*ColorDelta[c24, c13]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[prop[-k + q] . hv21] . v11 . GAD[lor1] . GAD[lorz] . GAD[lor2] . fcc[prop[-k + l] . v12] . hv22 . GAD[lor3] . GAD[lorz] . GAD[lor4]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*SUNTF[sun, c12, c14]*SUNTF[sun, c23, c21]) - ColorDelta[c12, c22]*ColorDelta[c23, c13]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[hv22 . prop[-k + l] . v12] . GAD[lor3] . GAD[lorz] . GAD[lor4] . fcc[v11 . prop[-k + q] . hv21] . GAD[lor1] . GAD[lorz] . GAD[lor2]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*SUNTF[sun, c11, c14]*SUNTF[sun, c24, c21];


dia2=(IntegrateP[dia2,{l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,l,k,dia3,sun,lor1,lor2,lor3,lor4,lorz},


dia3=-(ColorDelta[c12, c21]*ColorDelta[c24, c14]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[GAD[lor1] . GAD[lorz] . GAD[lor2] . v12] . v11 . prop[-k + q] . fcc[hv22 . GAD[lor3] . GAD[lorz] . GAD[lor4] . hv21] . prop[-k + l]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*SUNTF[sun, c11, c13]*SUNTF[sun, c23, c22]) - ColorDelta[c11, c21]*ColorDelta[c24, c13]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . prop[-k + q] . fcc[hv22 . GAD[lor3] . GAD[lorz] . GAD[lor4] . hv21] . prop[-k + l] . v12 . fcc[v11 . GAD[lor1] . GAD[lorz] . GAD[lor2]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*SUNTF[sun, c12, c14]*SUNTF[sun, c23, c22];


dia3=(IntegrateP[dia3,{l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia3=QEvaluate[I ScaleMu^(4(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,l,k,dia4,sun,lor1,lor2,lor3,lor4,lorz},


dia4=-(ColorDelta[c12, c21]*ColorDelta[c23, c13]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . GAD[lor1] . GAD[lorz] . GAD[lor2] . fcc[hv22 . prop[-k + l] . v12] . GAD[lor3] . GAD[lorz] . GAD[lor4] . hv21 . fcc[v11 . prop[-k + q]]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*SUNTF[sun, c11, c14]*SUNTF[sun, c24, c22]) - ColorDelta[c11, c21]*ColorDelta[c23, c14]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[GAD[lor3] . GAD[lorz] . GAD[lor4] . hv21] . hv22 . prop[-k + l] . fcc[v11 . GAD[lor1] . GAD[lorz] . GAD[lor2] . v12] . prop[-k + q]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*SUNTF[sun, c12, c13]*SUNTF[sun, c24, c22];


dia4=(IntegrateP[dia4,{l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia4=QEvaluate[I ScaleMu^(4(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,l,k,dia5,sun,lor1,lor2,lor3,lor4,lorz},


dia5=-(ColorDelta[c12, c22]*ColorDelta[c24, c13]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[prop[-k + l] . v12] . hv22 . GAD[lor3] . GAD[lorz] . GAD[lor4] . fcc[v11 . prop[-k + q] . hv21] . GAD[lor1] . GAD[lorz] . GAD[lor2]]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*SUNTF[sun, c11, c14]*SUNTF[sun, c23, c21]) - ColorDelta[c11, c21]*ColorDelta[c24, c14]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[hv22 . GAD[lor3] . GAD[lorz] . GAD[lor4] . hv21] . prop[-k + l] . fcc[v11 . GAD[lor1] . GAD[lorz] . GAD[lor2] . v12] . prop[-k + q]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*SUNTF[sun, c12, c13]*SUNTF[sun, c23, c22];


dia5=(IntegrateP[dia5,{l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia5=QEvaluate[I ScaleMu^(4(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,l,k,dia6,sun,lor1,lor2,lor3,lor4,lorz},


dia6=-(ColorDelta[c12, c21]*ColorDelta[c23, c14]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[GAD[lor3] . GAD[lorz] . GAD[lor4] . hv21] . hv22 . prop[-k + l] . fcc[GAD[lor1] . GAD[lorz] . GAD[lor2] . v12] . v11 . prop[-k + q]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*SUNTF[sun, c11, c13]*SUNTF[sun, c24, c22]) - ColorDelta[c11, c21]*ColorDelta[c23, c13]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[prop[-k + q]] . v11 . GAD[lor1] . GAD[lorz] . GAD[lor2] . fcc[hv22 . prop[-k + l] . v12] . GAD[lor3] . GAD[lorz] . GAD[lor4] . hv21]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*SUNTF[sun, c12, c14]*SUNTF[sun, c24, c22];


dia6=(IntegrateP[dia6,{l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia6=QEvaluate[I ScaleMu^(4(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,l,k,dia7,sun,lor1,lor2,lor3,lor4,lorz},


dia7=-(ColorDelta[c12, c22]*ColorDelta[c24, c14]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[prop[-k + l]] . hv22 . GAD[lor3] . GAD[lorz] . GAD[lor4] . fcc[v11 . prop[-k + q] . hv21] . GAD[lor1] . GAD[lorz] . GAD[lor2] . v12]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*SUNTF[sun, c11, c13]*SUNTF[sun, c23, c21]) - ColorDelta[c11, c22]*ColorDelta[c24, c14]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[hv22 . GAD[lor3] . GAD[lorz] . GAD[lor4]] . prop[-k + l] . fcc[v11 . GAD[lor1] . GAD[lorz] . GAD[lor2] . v12] . prop[-k + q] . hv21]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*SUNTF[sun, c12, c13]*SUNTF[sun, c23, c21];


dia7=(IntegrateP[dia7,{l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia7=QEvaluate[I ScaleMu^(4(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,l,k,dia8,sun,lor1,lor2,lor3,lor4,lorz},


dia8=-(ColorDelta[c12, c21]*ColorDelta[c24, c13]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . fcc[GAD[lor1] . GAD[lorz] . GAD[lor2]] . v11 . prop[-k + q] . fcc[hv22 . GAD[lor3] . GAD[lorz] . GAD[lor4] . hv21] . prop[-k + l] . v12]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*SUNTF[sun, c11, c14]*SUNTF[sun, c23, c22]) - ColorDelta[c12, c22]*ColorDelta[c23, c14]*Contract[atr[prop2[l, lor1, lor2, lor3, lor4] . GAD[lor3] . GAD[lorz] . GAD[lor4] . fcc[v11 . prop[-k + q] . hv21] . GAD[lor1] . GAD[lorz] . GAD[lor2] . v12 . fcc[hv22 . prop[-k + l]]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*SUNTF[sun, c11, c13]*SUNTF[sun, c24, c21];


dia8=(IntegrateP[dia8,{l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia8=QEvaluate[I ScaleMu^(4(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*-----------------------------------------------*)


(* ::Input::Initialization:: *)
ptype1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,k1,k2,l,k,dia1,sun,lor1,lor2,lor3,lor4,lorz},


dia1=(-I)*gStrong^2*ColorDelta[c12, c21]*ColorDelta[c24, c13]*Contract[atr[fcc[prop[k1 + l] . GAD[lorz] . prop[k1]] . v11 . prop[-k + q] . fcc[hv22 . prop[k2] . GAD[lorz] . prop[k2 + l] . hv21] . prop[-k + l] . v12]]*FAD[l]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*SUNTF[sun, c11, c14]*SUNTF[sun, c23, c22] - I*gStrong^2*ColorDelta[c11, c22]*ColorDelta[c23, c14]*Contract[atr[fcc[prop[k2] . GAD[lorz] . prop[k2 + l]] . hv22 . prop[-k + l] . fcc[v11 . prop[k1 + l] . GAD[lorz] . prop[k1] . v12] . prop[-k + q] . hv21]]*FAD[l]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*SUNTF[sun, c12, c13]*SUNTF[sun, c24, c21];


dia1=(IntegrateP[dia1,{k1,k2,l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,k1,k2,l,k,dia2,sun,lor1,lor2,lor3,lor4,lorz},


dia2=(-I)*gStrong^2*ColorDelta[c11, c22]*ColorDelta[c23, c13]*Contract[atr[fcc[prop[-k + q] . hv21] . v11 . prop[k1 + l] . GAD[lorz] . prop[k1] . fcc[hv22 . prop[-k + l] . v12] . prop[k2] . GAD[lorz] . prop[k2 + l]]]*FAD[l]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*SUNTF[sun, c12, c14]*SUNTF[sun, c24, c21] - I*gStrong^2*ColorDelta[c12, c21]*ColorDelta[c23, c14]*Contract[atr[fcc[prop[k2] . GAD[lorz] . prop[k2 + l] . hv21] . hv22 . prop[-k + l] . fcc[prop[k1 + l] . GAD[lorz] . prop[k1] . v12] . v11 . prop[-k + q]]]*FAD[l]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*SUNTF[sun, c11, c13]*SUNTF[sun, c24, c22];


dia2=(IntegrateP[dia2,{k1,k2,l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,k1,k2,l,k,dia3,sun,lor1,lor2,lor3,lor4,lorz},


dia3=(-I)*gStrong^2*ColorDelta[c11, c21]*ColorDelta[c24, c13]*Contract[atr[prop[-k + q] . fcc[hv22 . prop[k2] . GAD[lorz] . prop[k2 + l] . hv21] . prop[-k + l] . v12 . fcc[v11 . prop[k1 + l] . GAD[lorz] . prop[k1]]]]*FAD[l]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f22]*FlavorDelta[f24, f13]*SUNTF[sun, c12, c14]*SUNTF[sun, c23, c22] - I*gStrong^2*ColorDelta[c12, c22]*ColorDelta[c23, c14]*Contract[atr[prop[k2] . GAD[lorz] . prop[k2 + l] . fcc[v11 . prop[-k + q] . hv21] . prop[k1 + l] . GAD[lorz] . prop[k1] . v12 . fcc[hv22 . prop[-k + l]]]]*FAD[l]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f21]*SUNTF[sun, c11, c13]*SUNTF[sun, c24, c21];


dia3=(IntegrateP[dia3,{k1,k2,l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia3=QEvaluate[I ScaleMu^(4(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,k1,k2,l,k,dia4,sun,lor1,lor2,lor3,lor4,lorz},


dia4=(-I)*gStrong^2*ColorDelta[c11, c21]*ColorDelta[c24, c14]*Contract[atr[fcc[hv22 . prop[k2] . GAD[lorz] . prop[k2 + l] . hv21] . prop[-k + l] . fcc[v11 . prop[k1 + l] . GAD[lorz] . prop[k1] . v12] . prop[-k + q]]]*FAD[l]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*SUNTF[sun, c12, c13]*SUNTF[sun, c23, c22] - I*gStrong^2*ColorDelta[c12, c22]*ColorDelta[c23, c13]*Contract[atr[fcc[hv22 . prop[-k + l] . v12] . prop[k2] . GAD[lorz] . prop[k2 + l] . fcc[v11 . prop[-k + q] . hv21] . prop[k1 + l] . GAD[lorz] . prop[k1]]]*FAD[l]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f21]*SUNTF[sun, c11, c14]*SUNTF[sun, c24, c21];


dia4=(IntegrateP[dia4,{k1,k2,l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia4=QEvaluate[I ScaleMu^(4(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,k1,k2,l,k,dia5,sun,lor1,lor2,lor3,lor4,lorz},


dia5=(-I)*gStrong^2*ColorDelta[c12, c21]*ColorDelta[c23, c13]*Contract[atr[prop[k1 + l] . GAD[lorz] . prop[k1] . fcc[hv22 . prop[-k + l] . v12] . prop[k2] . GAD[lorz] . prop[k2 + l] . hv21 . fcc[v11 . prop[-k + q]]]]*FAD[l]*FlavorDelta[f11, f14]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*SUNTF[sun, c11, c14]*SUNTF[sun, c24, c22] - I*gStrong^2*ColorDelta[c11, c21]*ColorDelta[c23, c13]*Contract[atr[fcc[prop[-k + q]] . v11 . prop[k1 + l] . GAD[lorz] . prop[k1] . fcc[hv22 . prop[-k + l] . v12] . prop[k2] . GAD[lorz] . prop[k2 + l] . hv21]]*FAD[l]*FlavorDelta[f11, f21]*FlavorDelta[f12, f14]*FlavorDelta[f23, f13]*FlavorDelta[f24, f22]*SUNTF[sun, c12, c14]*SUNTF[sun, c24, c22];


dia5=(IntegrateP[dia5,{k1,k2,l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia5=QEvaluate[I ScaleMu^(4(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,k1,k2,l,k,dia6,sun,lor1,lor2,lor3,lor4,lorz},


dia6=(-I)*gStrong^2*ColorDelta[c12, c22]*ColorDelta[c24, c14]*Contract[atr[fcc[prop[-k + l]] . hv22 . prop[k2] . GAD[lorz] . prop[k2 + l] . fcc[v11 . prop[-k + q] . hv21] . prop[k1 + l] . GAD[lorz] . prop[k1] . v12]]*FAD[l]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*SUNTF[sun, c11, c13]*SUNTF[sun, c23, c21] - I*gStrong^2*ColorDelta[c12, c21]*ColorDelta[c24, c14]*Contract[atr[fcc[prop[k1 + l] . GAD[lorz] . prop[k1] . v12] . v11 . prop[-k + q] . fcc[hv22 . prop[k2] . GAD[lorz] . prop[k2 + l] . hv21] . prop[-k + l]]]*FAD[l]*FlavorDelta[f11, f13]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f14]*SUNTF[sun, c11, c13]*SUNTF[sun, c23, c22];


dia6=(IntegrateP[dia6,{k1,k2,l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia6=QEvaluate[I ScaleMu^(4(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,k1,k2,l,k,dia7,sun,lor1,lor2,lor3,lor4,lorz},


dia7=(-I)*gStrong^2*ColorDelta[c11, c22]*ColorDelta[c24, c14]*Contract[atr[fcc[hv22 . prop[k2] . GAD[lorz] . prop[k2 + l]] . prop[-k + l] . fcc[v11 . prop[k1 + l] . GAD[lorz] . prop[k1] . v12] . prop[-k + q] . hv21]]*FAD[l]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f23, f21]*FlavorDelta[f24, f14]*SUNTF[sun, c12, c13]*SUNTF[sun, c23, c21] - I*gStrong^2*ColorDelta[c11, c21]*ColorDelta[c23, c14]*Contract[atr[fcc[prop[k2] . GAD[lorz] . prop[k2 + l] . hv21] . hv22 . prop[-k + l] . fcc[v11 . prop[k1 + l] . GAD[lorz] . prop[k1] . v12] . prop[-k + q]]]*FAD[l]*FlavorDelta[f11, f21]*FlavorDelta[f12, f13]*FlavorDelta[f23, f14]*FlavorDelta[f24, f22]*SUNTF[sun, c12, c13]*SUNTF[sun, c24, c22];


dia7=(IntegrateP[dia7,{k1,k2,l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia7=QEvaluate[I ScaleMu^(4(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,k1,k2,l,k,dia8,sun,lor1,lor2,lor3,lor4,lorz},


dia8=(-I)*gStrong^2*ColorDelta[c12, c22]*ColorDelta[c24, c13]*Contract[atr[fcc[prop[-k + l] . v12] . hv22 . prop[k2] . GAD[lorz] . prop[k2 + l] . fcc[v11 . prop[-k + q] . hv21] . prop[k1 + l] . GAD[lorz] . prop[k1]]]*FAD[l]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*SUNTF[sun, c11, c14]*SUNTF[sun, c23, c21] - I*gStrong^2*ColorDelta[c11, c22]*ColorDelta[c24, c13]*Contract[atr[fcc[prop[-k + q] . hv21] . v11 . prop[k1 + l] . GAD[lorz] . prop[k1] . fcc[prop[-k + l] . v12] . hv22 . prop[k2] . GAD[lorz] . prop[k2 + l]]]*FAD[l]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f23, f21]*FlavorDelta[f24, f13]*SUNTF[sun, c12, c14]*SUNTF[sun, c23, c21];


dia8=(IntegrateP[dia8,{k1,k2,l,k}]//SUNSimplify)/.CF->(CA^2-1)/(2CA);

dia8=QEvaluate[I ScaleMu^(4(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*-----------------------------------------------*)


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
