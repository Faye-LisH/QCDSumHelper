(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
Tetraquarkmqgq::usage="Tetraquarkmqgq[q_,j1_,j2_] gives the second type of \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)q\[RightAngleBracket]\[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)Gq\[RightAngleBracket] contribution of \[LeftAngleBracket]j1 j2\[RightAngleBracket], with j1, j2 are tetraquark currents. "

Tetraquarkmqgq::inderr="Dummy indices conflict!"

Tetraquarkmqgq::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`Tetraquarkmqgq`"]


(* ::Code::Initialization::Plain:: *)
Options[Tetraquarkmqgq] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True
}



(* ::Code::Initialization::Plain:: *)
(*-----------------------------------*)
(* propagators and condensates *)

xpropm[x_,f_]=(DiracGamma[Momentum[x, D], D]*qfact1[(I^(1 - D)*qGamma[D/2])/(2*Pi^(D/2))])/Pair[Momentum[x, D], Momentum[x, D]]^(D/2) + 
 Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[(I^(2 - D)*qGamma[-1 + D/2]*quarkMass[f])/(4*Pi^(D/2))];


(* ::Code::Initialization::Plain:: *)
xprogm[x_,f_,sun_,lora_,lorb_]=DiracGamma[LorentzIndex[lora, D], D] . DiracGamma[LorentzIndex[lorb, D], D] . DiracGamma[Momentum[x, D], D]*
  Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[-1/32*qGamma[-1 + D/2]/(I^D*Pi^(D/2))]*SUNT[SUNIndex[sun]] + 
 DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorb, D], D] . DiracGamma[LorentzIndex[lora, D], D]*
  Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[1/(32*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[sun]] + 
 (2*DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lora, D], Momentum[x, D]]*Pair[LorentzIndex[lorb, D], Momentum[x, D]]*
   qfact1[qGamma[D/2]/(16*I^D*Pi^(D/2))]*SUNT[SUNIndex[sun]])/Pair[Momentum[x, D], Momentum[x, D]]^(D/2) + 
 Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*(DiracGamma[LorentzIndex[lorb, D], D] . DiracGamma[LorentzIndex[lora, D], D]*
    qfact1[qfact2[((I/64)*quarkMass[f])/(E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
   DiracGamma[LorentzIndex[lora, D], D] . DiracGamma[LorentzIndex[lorb, D], D]*qfact1[((-1/64*I)*qGamma[-2 + D/2]*quarkMass[f])/(I^D*Pi^(D/2))])*
  SUNT[SUNIndex[sun]] + DiracGamma[LorentzIndex[lorb, D], D] . DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lora, D], Momentum[x, D]]*
  Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[((I/32)*qGamma[-1 + D/2]*quarkMass[f])/(I^D*Pi^(D/2))]*SUNT[SUNIndex[sun]] + 
 DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorb, D], D]*Pair[LorentzIndex[lora, D], Momentum[x, D]]*
  Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[((I/32)*qGamma[-1 + D/2]*quarkMass[f])/(I^D*Pi^(D/2))]*SUNT[SUNIndex[sun]];



(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
Tetraquarkmqgq[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 Tetraquarkmqgq[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])


(* ::Code::Initialization::Plain:: *)
Tetraquarkmqgq[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,holdf=OptionValue[HoldFlavor],atr,fdir,files,pall,diagrams},



(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[Tetraquarkmqgq::curerr];
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
	Message[Tetraquarkmqgq::inderr];
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


diagrams={type1,type2,type3,type4,type5,type6,type7,type8};


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
			files=("Tetraquarkmqgq_"<>#)&/@files;
			
			
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


(* ::Input::Initialization::Plain:: *)
type1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,cc,tr,dot,sunSimplify,contract},


 dia1=Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprogm[x, f11, sun, lora, lorb], hv21]], v11, DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21])/(CA*CF*(-1 + D)*D)] + Condensate[{f13, "G", f23}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lora, lorb]] . SUNT[sun], v12]], xpropm[-x, f24]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, xprogm[x, f12, sun, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f13, "G", f24}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], v12]], hv22, xpropm[-x, f23]]]*tr[dot[cc[dot[xprogm[x, f11, sun, lora, lorb], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f22, "G", f12}]*contract[tr[dot[xpropm[x, f11], cc[dot[v11, DiracSigma[GAD[lora, lorb]] . SUNT[sun], hv21]]]]*tr[dot[cc[dot[xprogm[-x, f24, sun, lora, lorb], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xprogm[-x, f24, sun, lora, lorb]]]*tr[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], cc[dot[v11, xpropm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14])/(CA*CF*(-1 + D)*D)] + Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, xprogm[-x, f23, sun, lora, lorb], v12]], DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]*tr[dot[xpropm[x, f11], cc[dot[v11, xpropm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14])/(CA*CF*(-1 + D)*D)];


 dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia1=FourierXP[dia1,{x,q}];


 dia1=dia1/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia1=QEvaluate[I ScaleMu^(2(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia1=QGather[Expand[dia1]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,cc,tr,dot,sunSimplify,contract},


 dia2=Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[hv22, xprogm[-x, f23, sun, lora, lorb], v12]], xpropm[-x, f24]]]*tr[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], cc[dot[v11, xpropm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], hv21]], v11, xpropm[x, f12]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xprogm[-x, f23, sun, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14])/(CA*CF*(-1 + D)*D)] + Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xprogm[-x, f23, sun, lora, lorb]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14])/(CA*CF*(-1 + D)*D)] + Condensate[{f14, "G", f23}]*contract[tr[dot[xpropm[x, f11], cc[dot[v11, xprogm[x, f12, sun, lora, lorb], hv21]]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14])/(CA*CF*(-1 + D)*D)] + Condensate[{f13, "G", f24}]*contract[tr[dot[xprogm[x, f11, sun, lora, lorb], cc[dot[v11, xpropm[x, f12], hv21]]]]*tr[dot[cc[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f13, "G", f23}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lora, lorb]] . SUNT[sun], v12]], xprogm[-x, f24, sun, lora, lorb]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14])/(CA*CF*(-1 + D)*D)];


 dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia2=FourierXP[dia2,{x,q}];


 dia2=dia2/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia2=QEvaluate[I ScaleMu^(2(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia2=QGather[Expand[dia2]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,cc,tr,dot,sunSimplify,contract},


 dia3=Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], hv21]], v11, xprogm[x, f12, sun, lora, lorb]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21])/(CA*CF*(-1 + D)*D)] + Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xpropm[-x, f24]]]*tr[dot[xprogm[x, f11, sun, lora, lorb], cc[dot[v11, DiracSigma[GAD[lora, lorb]] . SUNT[sun], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22])/(CA*CF*(-1 + D)*D)] + Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, xprogm[-x, f23, sun, lora, lorb], v12]], xpropm[-x, f24]]]*tr[dot[cc[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f13, "G", f23}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lora, lorb]] . SUNT[sun], v12]], xpropm[-x, f24]]]*tr[dot[xpropm[x, f11], cc[dot[v11, xprogm[x, f12, sun, lora, lorb], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f14, "G", f23}]*contract[tr[dot[xpropm[x, f11], cc[dot[v11, xpropm[x, f12], hv21]]]]*tr[dot[cc[dot[xprogm[-x, f24, sun, lora, lorb], v12]], hv22, DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, xprogm[-x, f23, sun, lora, lorb], v12]], DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14])/(CA*CF*(-1 + D)*D)];


 dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia3=FourierXP[dia3,{x,q}];


 dia3=dia3/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia3=QEvaluate[I ScaleMu^(2(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia3=QGather[Expand[dia3]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,cc,tr,dot,sunSimplify,contract},


 dia4=Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xpropm[-x, f24]]]*tr[dot[cc[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], hv21]], v11, xprogm[x, f12, sun, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21])/(CA*CF*(-1 + D)*D)] + Condensate[{f21, "G", f11}]*contract[tr[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], cc[dot[v11, xpropm[x, f12], hv21]]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xprogm[-x, f23, sun, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14])/(CA*CF*(-1 + D)*D)] + Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, xprogm[x, f12, sun, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14])/(CA*CF*(-1 + D)*D)] + Condensate[{f13, "G", f24}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], v12]], hv22, xpropm[-x, f23]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, xprogm[x, f12, sun, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]*tr[dot[xprogm[x, f11, sun, lora, lorb], cc[dot[v11, xpropm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14])/(CA*CF*(-1 + D)*D)] + Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, xprogm[x, f12, sun, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14])/(CA*CF*(-1 + D)*D)];


 dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia4=FourierXP[dia4,{x,q}];


 dia4=dia4/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia4=QEvaluate[I ScaleMu^(2(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia4=QGather[Expand[dia4]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia5,cc,tr,dot,sunSimplify,contract},


 dia5=Condensate[{f13, "G", f23}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lora, lorb]] . SUNT[sun], v12]], xpropm[-x, f24]]]*tr[dot[cc[dot[xprogm[x, f11, sun, lora, lorb], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, xprogm[-x, f23, sun, lora, lorb], v12]], xpropm[-x, f24]]]*tr[dot[xpropm[x, f11], cc[dot[v11, DiracSigma[GAD[lora, lorb]] . SUNT[sun], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f14, "G", f23}]*contract[tr[dot[xprogm[x, f11, sun, lora, lorb], cc[dot[v11, xpropm[x, f12], hv21]]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14])/(CA*CF*(-1 + D)*D)] + Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[xprogm[x, f11, sun, lora, lorb], hv21]], v11, xpropm[x, f12]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14])/(CA*CF*(-1 + D)*D)] + Condensate[{f13, "G", f24}]*contract[tr[dot[xpropm[x, f11], cc[dot[v11, xprogm[x, f12, sun, lora, lorb], hv21]]]]*tr[dot[cc[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xprogm[-x, f24, sun, lora, lorb]]]*tr[dot[cc[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14])/(CA*CF*(-1 + D)*D)];


 dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia5=FourierXP[dia5,{x,q}];


 dia5=dia5/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia5=QEvaluate[I ScaleMu^(2(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia5=QGather[Expand[dia5]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia6,cc,tr,dot,sunSimplify,contract},


 dia6=Condensate[{f22, "G", f12}]*contract[tr[dot[xprogm[x, f11, sun, lora, lorb], cc[dot[v11, DiracSigma[GAD[lora, lorb]] . SUNT[sun], hv21]]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22])/(CA*CF*(-1 + D)*D)] + Condensate[{f22, "G", f11}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], hv21]], v11, xpropm[x, f12]]]*tr[dot[cc[dot[xprogm[-x, f24, sun, lora, lorb], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f13, "G", f24}]*contract[tr[dot[cc[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], v12]], hv22, xprogm[-x, f23, sun, lora, lorb]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f13, "G", f24}]*contract[tr[dot[xpropm[x, f11], cc[dot[v11, xpropm[x, f12], hv21]]]]*tr[dot[cc[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], v12]], hv22, xprogm[-x, f23, sun, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]*tr[dot[cc[dot[xprogm[x, f11, sun, lora, lorb], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14])/(CA*CF*(-1 + D)*D)] + Condensate[{f14, "G", f24}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]*tr[dot[xpropm[x, f11], cc[dot[v11, xprogm[x, f12, sun, lora, lorb], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14])/(CA*CF*(-1 + D)*D)];


 dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia6=FourierXP[dia6,{x,q}];


 dia6=dia6/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia6=QEvaluate[I ScaleMu^(2(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia6=QGather[Expand[dia6]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia7,cc,tr,dot,sunSimplify,contract},


 dia7=Condensate[{f21, "G", f11}]*contract[tr[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], cc[dot[v11, xprogm[x, f12, sun, lora, lorb], hv21]]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22])/(CA*CF*(-1 + D)*D)] + Condensate[{f21, "G", f11}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xpropm[-x, f24]]]*tr[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], cc[dot[v11, xprogm[x, f12, sun, lora, lorb], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22])/(CA*CF*(-1 + D)*D)] + Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[hv22, xprogm[-x, f23, sun, lora, lorb], v12]], xpropm[-x, f24]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f21, "G", f11}]*contract[tr[dot[DiracSigma[GAD[lora, lorb]] . SUNT[sun], cc[dot[v11, xpropm[x, f12], hv21]]]]*tr[dot[cc[dot[xprogm[-x, f24, sun, lora, lorb], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[xprogm[-x, f24, sun, lora, lorb], v12]], hv22, xpropm[-x, f23]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xprogm[-x, f24, sun, lora, lorb]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14])/(CA*CF*(-1 + D)*D)];


 dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia7=FourierXP[dia7,{x,q}];


 dia7=dia7/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia7=QEvaluate[I ScaleMu^(2(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia7=QGather[Expand[dia7]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
type8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia8,cc,tr,dot,sunSimplify,contract},


 dia8=Condensate[{f21, "G", f12}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xpropm[-x, f24]]]*tr[dot[cc[dot[xprogm[x, f11, sun, lora, lorb], hv21]], v11, DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21])/(CA*CF*(-1 + D)*D)] + Condensate[{f13, "G", f23}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lora, lorb]] . SUNT[sun], v12]], xpropm[-x, f24]]]*tr[dot[xprogm[x, f11, sun, lora, lorb], cc[dot[v11, xpropm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f22, "G", f12}]*contract[tr[dot[xpropm[x, f11], cc[dot[v11, DiracSigma[GAD[lora, lorb]] . SUNT[sun], hv21]]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xprogm[-x, f23, sun, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14])/(CA*CF*(-1 + D)*D)] + Condensate[{f14, "G", f23}]*contract[tr[dot[cc[dot[xprogm[-x, f24, sun, lora, lorb], v12]], hv22, DiracSigma[GAD[lora, lorb]] . SUNT[sun]]]*tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[-1/4*(ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13])/(CA*CF*(-1 + D)*D)] + Condensate[{f22, "G", f12}]*contract[tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xprogm[-x, f24, sun, lora, lorb]]]*tr[dot[xpropm[x, f11], cc[dot[v11, DiracSigma[GAD[lora, lorb]] . SUNT[sun], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14])/(CA*CF*(-1 + D)*D)] + Condensate[{f13, "G", f23}]*contract[tr[dot[cc[dot[hv22, DiracSigma[GAD[lora, lorb]] . SUNT[sun], v12]], xprogm[-x, f24, sun, lora, lorb]]]*tr[dot[xpropm[x, f11], cc[dot[v11, xpropm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[-1/4*(ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14])/(CA*CF*(-1 + D)*D)];


 dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia8=FourierXP[dia8,{x,q}];


 dia8=dia8/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia8=QEvaluate[I ScaleMu^(2(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia8=QGather[Expand[dia8]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
