(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
Tetraquarkmqqgg::usage="Tetraquarkmqqgg[q_,j1_,j2_] gives the second type of \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)q\[RightAngleBracket]\[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)Gq\[RightAngleBracket] contribution of \[LeftAngleBracket]j1 j2\[RightAngleBracket], with j1, j2 are tetraquark currents. "

Tetraquarkmqqgg::inderr="Dummy indices conflict!"

Tetraquarkmqqgg::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`Tetraquarkmqqgg`"]


(* ::Code::Initialization::Plain:: *)
Options[Tetraquarkmqqgg] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True
}


(* ::Code::Initialization::Plain:: *)
(*-----------------------------------*)
(* condensate *)
(* \[LeftAngleBracket]\!\(
\*SubsuperscriptBox[\(G\), \(\[Mu]\[Nu]\), \(a\)]
\*SubsuperscriptBox[\(G\), \(\[Alpha]\[Beta]\), \(b\)]\)\[RightAngleBracket] *)
congg[lora_,lorb_,lorc_,lord_]=1/(2 CA CF D(D-1))(MTD[lora,lorc]MTD[lorb,lord]-MTD[lora,lord]MTD[lorb,lorc]);
(* m\[LeftAngleBracket]Overscript[q, _]q\[RightAngleBracket] *)
cqqm[x_,f_]=-1/4*1/CA + ((I/4)*DiracGamma[Momentum[x, D], D]*quarkMass[f])/(CA*D);

(*-----------------------------------*)
(* propagators and condensates *)

xpropm[x_,f_]=(DiracGamma[Momentum[x, D], D]*qfact1[(I^(1 - D)*qGamma[D/2])/(2*Pi^(D/2))])/Pair[Momentum[x, D], Momentum[x, D]]^(D/2) + 
 Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[(I^(2 - D)*qGamma[-1 + D/2]*quarkMass[f])/(4*Pi^(D/2))];
 
 (*-------------------------------*)
 xprop[x_]=(I^(1 - D)*GSD[x]*qGamma[D/2])/(2*Pi^(D/2)*SPD[x, x]^(D/2));
 
 xpromgg[x_,f_]=(qGamma[-3 + D/2]*quarkMass[f]*SPD[x, x]^(3 - D/2))/(128*I^D*CA*D*Pi^(D/2));
 (*I/4FAD[{q}](GSD[q]+quarkMass[f]).GAD[lorb].SUNT[cola].SUNT[colb].FourDivergence[FAD[{q}](GSD[q]+quarkMass[f]).GAD[lord].FourDivergence[(GSD[q]+quarkMass[f])FAD[{q}],FVD[q,lorc]],FVD[q,lora]];

(MTD[lora,lorc]MTD[lorb,lord]-MTD[lora,lord]MTD[lorb,lorc])SUNDelta[cola,colb]%//Contract//SUNSimplify//DiracSimplify;
%/.Power[_quarkMass,_]:>0;
%//FCI;
%/.Power[FeynAmpDenominator[PropagatorDenominator[Momentum[q,D],0]],n_]:>SPD[q]^(-n)/.FeynAmpDenominator[PropagatorDenominator[Momentum[q,D],0],pp___]:>SPD[q]^(-Length[{pp}]-1)//FCI//Simplify;
%/(2 CA CF D(D-1))

 xpromgg=FourierPX[%,{q,x}]//FCE//InputForm*)


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
Tetraquarkmqqgg[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 Tetraquarkmqqgg[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])


(* ::Code::Initialization::Plain:: *)
Tetraquarkmqqgg[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,holdf=OptionValue[HoldFlavor],atr,fdir,files,pall,diagrams},



(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[Tetraquarkmqqgg::curerr];
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
	Message[Tetraquarkmqqgg::inderr];
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


diagrams={typeA1,typeA2,typeA3,typeA4,typeA5,typeA6,typeA7,typeA8,typeB1,typeB2,typeB3,typeB4,typeB5,typeB6,typeB7,typeB8};


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
			files=("Tetraquarkmqqgg_"<>#)&/@files;
			
			
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
typeA1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,cc,tr,dot,sunSimplify,contract},


 dia1=4*Pi*Condensate["gg"]*Condensate[{f22, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, xprogm[-x, f23, sun, lora, lorb], v12]], xpropm[-x, f24]]]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, xprogm[x, f12, sun, lorc, lord]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13]] + 4*Pi*Condensate["gg"]*Condensate[{f13, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], xprogm[-x, f24, sun, lora, lorb]]]*tr[dot[xprogm[x, f11, sun, lorc, lord], cc[dot[v11, xpropm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14]] + 4*Pi*Condensate["gg"]*Condensate[{f21, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xprogm[-x, f24, sun, lora, lorb]]]*tr[dot[cc[dot[xprogm[x, f11, sun, lorc, lord], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14]];


 dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia1=FourierXP[dia1,{x,q}];


 dia1=dia1/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia1=QEvaluate[I ScaleMu^(2(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia1=QGather[Expand[dia1]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeA2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,cc,tr,dot,sunSimplify,contract},


 dia2=4*Pi*Condensate["gg"]*Condensate[{f21, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, xprogm[-x, f23, sun, lora, lorb], v12]], xpropm[-x, f24]]]*tr[dot[cqqm[x, f11], cc[dot[v11, xprogm[x, f12, sun, lorc, lord], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13]] + 4*Pi*Condensate["gg"]*Condensate[{f22, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, xprogm[x, f12, sun, lorc, lord]]]*tr[dot[cc[dot[xprogm[-x, f24, sun, lora, lorb], v12]], hv22, xpropm[-x, f23]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13]] + 4*Pi*Condensate["gg"]*Condensate[{f22, f12}]*contract[tr[dot[xpropm[x, f11], cc[dot[v11, cqqm[x, f12], hv21]]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[hv22, xprogm[-x, f23, sun, lorc, lord], v12]], xprogm[-x, f24, sun, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]];


 dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia2=FourierXP[dia2,{x,q}];


 dia2=dia2/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia2=QEvaluate[I ScaleMu^(2(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia2=QGather[Expand[dia2]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeA3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,cc,tr,dot,sunSimplify,contract},


 dia3=4*Pi*Condensate["gg"]*Condensate[{f14, f23}]*contract[tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, cqqm[-x, f23]]]*tr[dot[congg[lora, lorb, lorc, lord], xprogm[x, f11, sun, lorc, lord], cc[dot[v11, xprogm[x, f12, sun, lora, lorb], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]] + 4*Pi*Condensate["gg"]*Condensate[{f21, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[xprogm[x, f11, sun, lorc, lord], hv21]], v11, cqqm[x, f12]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xprogm[-x, f23, sun, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14]] + 4*Pi*Condensate["gg"]*Condensate[{f14, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[xprogm[-x, f24, sun, lora, lorb], v12]], hv22, cqqm[-x, f23]]]*tr[dot[cc[dot[xprogm[x, f11, sun, lorc, lord], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13]];


 dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia3=FourierXP[dia3,{x,q}];


 dia3=dia3/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia3=QEvaluate[I ScaleMu^(2(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia3=QGather[Expand[dia3]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeA4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,cc,tr,dot,sunSimplify,contract},


 dia4=4*Pi*Condensate["gg"]*Condensate[{f13, f23}]*contract[tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], xpropm[-x, f24]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[xprogm[x, f11, sun, lorc, lord], hv21]], v11, xprogm[x, f12, sun, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]] + 4*Pi*Condensate["gg"]*Condensate[{f14, f24}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, xprogm[-x, f23, sun, lora, lorb], v12]], cqqm[-x, f24]]]*tr[dot[xpropm[x, f11], cc[dot[v11, xprogm[x, f12, sun, lorc, lord], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13]] + 4*Pi*Condensate["gg"]*Condensate[{f22, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[xprogm[x, f11, sun, lorc, lord], cc[dot[v11, cqqm[x, f12], hv21]]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xprogm[-x, f23, sun, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14]];


 dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia4=FourierXP[dia4,{x,q}];


 dia4=dia4/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia4=QEvaluate[I ScaleMu^(2(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia4=QGather[Expand[dia4]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeA5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia5,cc,tr,dot,sunSimplify,contract},


 dia5=4*Pi*Condensate["gg"]*Condensate[{f21, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, xprogm[-x, f23, sun, lora, lorb], v12]], xpropm[-x, f24]]]*tr[dot[cc[dot[xprogm[x, f11, sun, lorc, lord], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13]] + 4*Pi*Condensate["gg"]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[xpropm[x, f11], hv21]], v11, cqqm[x, f12]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[xprogm[-x, f24, sun, lorc, lord], v12]], hv22, xprogm[-x, f23, sun, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]] + 4*Pi*Condensate["gg"]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, xpropm[x, f12]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[hv22, xprogm[-x, f23, sun, lorc, lord], v12]], xprogm[-x, f24, sun, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]];


 dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia5=FourierXP[dia5,{x,q}];


 dia5=dia5/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia5=QEvaluate[I ScaleMu^(2(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia5=QGather[Expand[dia5]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeA6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia6,cc,tr,dot,sunSimplify,contract},


 dia6=4*Pi*Condensate["gg"]*Condensate[{f14, f24}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, xprogm[-x, f23, sun, lora, lorb], v12]], cqqm[-x, f24]]]*tr[dot[xprogm[x, f11, sun, lorc, lord], cc[dot[v11, xpropm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13]] + 4*Pi*Condensate["gg"]*Condensate[{f14, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[xpropm[x, f11], cc[dot[v11, xprogm[x, f12, sun, lorc, lord], hv21]]]]*tr[dot[cc[dot[xprogm[-x, f24, sun, lora, lorb], v12]], hv22, cqqm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13]] + 4*Pi*Condensate["gg"]*Condensate[{f21, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xprogm[-x, f24, sun, lora, lorb]]]*tr[dot[cqqm[x, f11], cc[dot[v11, xprogm[x, f12, sun, lorc, lord], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14]];


 dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia6=FourierXP[dia6,{x,q}];


 dia6=dia6/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia6=QEvaluate[I ScaleMu^(2(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia6=QGather[Expand[dia6]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeA7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia7,cc,tr,dot,sunSimplify,contract},


 dia7=4*Pi*Condensate["gg"]*Condensate[{f13, f24}]*contract[tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, xpropm[-x, f23]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[xprogm[x, f11, sun, lorc, lord], hv21]], v11, xprogm[x, f12, sun, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]] + 4*Pi*Condensate["gg"]*Condensate[{f21, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[xprogm[-x, f24, sun, lora, lorb], v12]], hv22, xpropm[-x, f23]]]*tr[dot[cc[dot[xprogm[x, f11, sun, lorc, lord], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13]] + 4*Pi*Condensate["gg"]*Condensate[{f13, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], xprogm[-x, f24, sun, lora, lorb]]]*tr[dot[cc[dot[xprogm[x, f11, sun, lorc, lord], hv21]], v11, xpropm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14]];


 dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia7=FourierXP[dia7,{x,q}];


 dia7=dia7/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia7=QEvaluate[I ScaleMu^(2(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia7=QGather[Expand[dia7]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeA8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia8,cc,tr,dot,sunSimplify,contract},


 dia8=4*Pi*Condensate["gg"]*Condensate[{f21, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cqqm[x, f11], cc[dot[v11, xprogm[x, f12, sun, lorc, lord], hv21]]]]*tr[dot[cc[dot[xpropm[-x, f24], v12]], hv22, xprogm[-x, f23, sun, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14]] + 4*Pi*Condensate["gg"]*Condensate[{f14, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[xprogm[x, f11, sun, lorc, lord], cc[dot[v11, xpropm[x, f12], hv21]]]]*tr[dot[cc[dot[xprogm[-x, f24, sun, lora, lorb], v12]], hv22, cqqm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13]] + 4*Pi*Condensate["gg"]*Condensate[{f22, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, xpropm[-x, f23], v12]], xprogm[-x, f24, sun, lora, lorb]]]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, xprogm[x, f12, sun, lorc, lord]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14]];


 dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia8=FourierXP[dia8,{x,q}];


 dia8=dia8/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia8=QEvaluate[I ScaleMu^(2(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia8=QGather[Expand[dia8]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*--------------------------------------*)


(* ::Input::Initialization::Plain:: *)
typeB1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,cc,tr,dot,sunSimplify,contract},


 dia1=-((Pi*Condensate["gg"]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[xprop[x], hv21]], v11]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xpromgg[-x, f23]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA) - (Pi*Condensate["gg"]*Condensate[{f14, f24}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]]]]*tr[dot[xpromgg[x, f11], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]])/CA - (Pi*Condensate["gg"]*Condensate[{f13, f23}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[xpromgg[x, f11], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]])/CA;


 dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia1=FourierXP[dia1,{x,q}];


 dia1=QEvaluate[I ScaleMu^(2(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia1=QGather[Expand[dia1],q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeB2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,cc,tr,dot,sunSimplify,contract},


 dia2=-((Pi*Condensate["gg"]*Condensate[{f14, f23}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[xpromgg[x, f11], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA) - (Pi*Condensate["gg"]*Condensate[{f21, f11}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xpromgg[-x, f24], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA - (Pi*Condensate["gg"]*Condensate[{f22, f12}]*contract[tr[dot[xprop[x], cc[dot[v11, hv21]]]]*tr[dot[cc[dot[xpromgg[-x, f24], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA;


 dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia2=FourierXP[dia2,{x,q}];


 dia2=QEvaluate[I ScaleMu^(2(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia2=QGather[Expand[dia2],q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeB3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,cc,tr,dot,sunSimplify,contract},


 dia3=-((Pi*Condensate["gg"]*Condensate[{f13, f24}]*contract[tr[dot[cc[dot[v12]], hv22, xpromgg[-x, f23]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA) - (Pi*Condensate["gg"]*Condensate[{f21, f11}]*contract[tr[dot[cc[dot[v11, xpromgg[x, f12], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA - (Pi*Condensate["gg"]*Condensate[{f13, f23}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xpromgg[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]])/CA;


 dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia3=FourierXP[dia3,{x,q}];


 dia3=QEvaluate[I ScaleMu^(2(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia3=QGather[Expand[dia3],q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeB4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,cc,tr,dot,sunSimplify,contract},


 dia4=-((Pi*Condensate["gg"]*Condensate[{f13, f24}]*contract[tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[v12]], hv22, xpromgg[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA) - (Pi*Condensate["gg"]*Condensate[{f13, f24}]*contract[tr[dot[xpromgg[x, f11], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA - (Pi*Condensate["gg"]*Condensate[{f13, f23}]*contract[tr[dot[cc[dot[hv22, v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, xpromgg[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]])/CA;


 dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia4=FourierXP[dia4,{x,q}];


 dia4=QEvaluate[I ScaleMu^(2(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia4=QGather[Expand[dia4],q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeB5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia5,cc,tr,dot,sunSimplify,contract},


 dia5=-((Pi*Condensate["gg"]*Condensate[{f13, f24}]*contract[tr[dot[cc[dot[v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xpromgg[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA) - (Pi*Condensate["gg"]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[hv22, xpromgg[-x, f23], v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]])/CA - (Pi*Condensate["gg"]*Condensate[{f22, f12}]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xpromgg[-x, f24]]]*tr[dot[xprop[x], cc[dot[v11, hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]])/CA;


 dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia5=FourierXP[dia5,{x,q}];


 dia5=QEvaluate[I ScaleMu^(2(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia5=QGather[Expand[dia5],q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeB6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia6,cc,tr,dot,sunSimplify,contract},


 dia6=-((Pi*Condensate["gg"]*Condensate[{f14, f23}]*contract[tr[dot[cc[dot[xpromgg[-x, f24], v12]], hv22]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA) - (Pi*Condensate["gg"]*Condensate[{f14, f23}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[xprop[x], cc[dot[v11, xpromgg[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA - (Pi*Condensate["gg"]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[xpromgg[x, f11], hv21]], v11]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]])/CA;


 dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia6=FourierXP[dia6,{x,q}];


 dia6=QEvaluate[I ScaleMu^(2(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia6=QGather[Expand[dia6],q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeB7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia7,cc,tr,dot,sunSimplify,contract},


 dia7=-((Pi*Condensate["gg"]*Condensate[{f14, f23}]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22]]*tr[dot[cc[dot[xpromgg[x, f11], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA) - (Pi*Condensate["gg"]*Condensate[{f21, f11}]*contract[tr[dot[cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xpromgg[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]])/CA - (Pi*Condensate["gg"]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[hv22, xpromgg[-x, f23], v12]], xprop[-x]]]*tr[dot[cc[dot[hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]])/CA;


 dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia7=FourierXP[dia7,{x,q}];


 dia7=QEvaluate[I ScaleMu^(2(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia7=QGather[Expand[dia7],q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeB8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia8,cc,tr,dot,sunSimplify,contract},


 dia8=-((Pi*Condensate["gg"]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[xprop[x], hv21]], v11]]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xpromgg[-x, f24]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]])/CA) - (Pi*Condensate["gg"]*Condensate[{f14, f24}]*contract[tr[dot[cc[dot[hv22, xpromgg[-x, f23], v12]]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]])/CA - (Pi*Condensate["gg"]*Condensate[{f13, f23}]*contract[tr[dot[cc[dot[hv22, v12]], xpromgg[-x, f24]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]])/CA;


 dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia8=FourierXP[dia8,{x,q}];


 dia8=QEvaluate[I ScaleMu^(2(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia8=QGather[Expand[dia8],q,ShowasTable->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
(*type1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia1,cc,tr,dot,sunSimplify,contract},


 dia1=Condensate["gg"]*Condensate[{f13, f24}]*contract[tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f23]), {q, -x}]]]*tr[dot[congg[lora, lorb, lorc, lord], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], cc[dot[v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lora]], {q, x}]]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]] + Condensate["gg"]*Condensate[{f22, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]]]]*tr[dot[cc[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f24]), {q, -x}], v12]], hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14]] + Condensate["gg"]*Condensate[{f14, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f11]), {q, x}], cc[dot[v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]], hv21]]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]], v12]], hv22, cqqm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13]];


 dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia1=FourierXP[dia1,{x,q}];


 dia1=dia1/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia1=QEvaluate[I ScaleMu^(2(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia1=QGather[Expand[dia1]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia2,cc,tr,dot,sunSimplify,contract},


 dia2=Condensate["gg"]*Condensate[{f14, f23}]*contract[tr[dot[cc[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f24]), {q, -x}], v12]], hv22, cqqm[-x, f23]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], hv21]], v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lora]], {q, x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]] + Condensate["gg"]*Condensate[{f21, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], hv21]], v11, cqqm[x, f12]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]], v12]], hv22, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f23]), {q, -x}]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13]] + Condensate["gg"]*Condensate[{f13, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], hv21]], v11, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f12]), {q, x}]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14]];


 dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia2=FourierXP[dia2,{x,q}];


 dia2=dia2/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia2=QEvaluate[I ScaleMu^(2(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia2=QGather[Expand[dia2]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia3,cc,tr,dot,sunSimplify,contract},


 dia3=Condensate["gg"]*Condensate[{f13, f24}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f11]), {q, x}], cc[dot[v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]], hv21]]]]*tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14]] + Condensate["gg"]*Condensate[{f21, f11}]*contract[tr[dot[cqqm[x, f11], cc[dot[v11, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f12]), {q, x}], hv21]]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lorc]], {q, -x}]]], v12]], hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]] + Condensate["gg"]*Condensate[{f21, f11}]*contract[tr[dot[cqqm[x, f11], cc[dot[v11, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f12]), {q, x}], hv21]]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lorc]], {q, -x}]]], v12]], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]];


 dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia3=FourierXP[dia3,{x,q}];


 dia3=dia3/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia3=QEvaluate[I ScaleMu^(2(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia3=QGather[Expand[dia3]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia4,cc,tr,dot,sunSimplify,contract},


 dia4=Condensate["gg"]*Condensate[{f14, f24}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]], v12]], cqqm[-x, f24]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], hv21]], v11, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f12]), {q, x}]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13]] + Condensate["gg"]*Condensate[{f22, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]], v12]], FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f24]), {q, -x}]]]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13]] + Condensate["gg"]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f12]), {q, x}]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lorc]], {q, -x}]]], v12]], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]];


 dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia4=FourierXP[dia4,{x,q}];


 dia4=dia4/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia4=QEvaluate[I ScaleMu^(2(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia4=QGather[Expand[dia4]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia5,cc,tr,dot,sunSimplify,contract},


 dia5=Condensate["gg"]*Condensate[{f13, f23}]*contract[tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f24]), {q, -x}]]]*tr[dot[congg[lora, lorb, lorc, lord], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], cc[dot[v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lora]], {q, x}]]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]] + Condensate["gg"]*Condensate[{f21, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f23]), {q, -x}], v12]], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14]] + Condensate["gg"]*Condensate[{f22, f12}]*contract[tr[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f11]), {q, x}], cc[dot[v11, cqqm[x, f12], hv21]]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lorc]], {q, -x}]]], v12]], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]];


 dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia5=FourierXP[dia5,{x,q}];


 dia5=dia5/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia5=QEvaluate[I ScaleMu^(2(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia5=QGather[Expand[dia5]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia6,cc,tr,dot,sunSimplify,contract},


 dia6=Condensate["gg"]*Condensate[{f14, f24}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]], v12]], cqqm[-x, f24]]]*tr[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], cc[dot[v11, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f12]), {q, x}], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13]] + Condensate["gg"]*Condensate[{f21, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]], v12]], FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f24]), {q, -x}]]]*tr[dot[cqqm[x, f11], cc[dot[v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13]] + Condensate["gg"]*Condensate[{f21, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cqqm[x, f11], cc[dot[v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]], hv21]]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]], v12]], hv22, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f23]), {q, -x}]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13]];


 dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia6=FourierXP[dia6,{x,q}];


 dia6=dia6/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia6=QEvaluate[I ScaleMu^(2(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia6=QGather[Expand[dia6]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia7,cc,tr,dot,sunSimplify,contract},


 dia7=Condensate["gg"]*Condensate[{f13, f24}]*contract[tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f23]), {q, -x}]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], hv21]], v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lora]], {q, x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]] + Condensate["gg"]*Condensate[{f14, f24}]*contract[tr[dot[cc[dot[hv22, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f23]), {q, -x}], v12]], cqqm[-x, f24]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], hv21]], v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lora]], {q, x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]] + Condensate["gg"]*Condensate[{f22, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], cc[dot[v11, cqqm[x, f12], hv21]]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]], v12]], hv22, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f23]), {q, -x}]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13]];


 dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia7=FourierXP[dia7,{x,q}];


 dia7=dia7/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia7=QEvaluate[I ScaleMu^(2(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia7=QGather[Expand[dia7]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia8,cc,tr,dot,sunSimplify,contract},


 dia8=Condensate["gg"]*Condensate[{f21, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f24]), {q, -x}], v12]], hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14]] + Condensate["gg"]*Condensate[{f14, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f11]), {q, x}], hv21]], v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]], v12]], hv22, cqqm[-x, f23]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13]] + Condensate["gg"]*Condensate[{f22, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f23]), {q, -x}], v12]], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]]]]*tr[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], cc[dot[v11, cqqm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14]];


 dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia8=FourierXP[dia8,{x,q}];


 dia8=dia8/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia8=QEvaluate[I ScaleMu^(2(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia8=QGather[Expand[dia8]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type9[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia9,cc,tr,dot,sunSimplify,contract},


 dia9=Condensate["gg"]*Condensate[{f14, f24}]*contract[tr[dot[cc[dot[hv22, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f23]), {q, -x}], v12]], cqqm[-x, f24]]]*tr[dot[congg[lora, lorb, lorc, lord], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], cc[dot[v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lora]], {q, x}]]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]] + Condensate["gg"]*Condensate[{f21, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]], v12]], FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f24]), {q, -x}]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], hv21]], v11, cqqm[x, f12]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13]] + Condensate["gg"]*Condensate[{f13, f24}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]]]]*tr[dot[cc[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f11]), {q, x}], hv21]], v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14]];


 dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia9=FourierXP[dia9,{x,q}];


 dia9=dia9/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia9=QEvaluate[I ScaleMu^(2(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia9=QGather[Expand[dia9]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type10[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia10,cc,tr,dot,sunSimplify,contract},


 dia10=Condensate["gg"]*Condensate[{f14, f23}]*contract[tr[dot[cc[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f24]), {q, -x}], v12]], hv22, cqqm[-x, f23]]]*tr[dot[congg[lora, lorb, lorc, lord], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], cc[dot[v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lora]], {q, x}]]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]] + Condensate["gg"]*Condensate[{f14, f24}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]], v12]], cqqm[-x, f24]]]*tr[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f11]), {q, x}], cc[dot[v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13]] + Condensate["gg"]*Condensate[{f14, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], hv21]], v11, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f12]), {q, x}]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]], v12]], hv22, cqqm[-x, f23]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13]];


 dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia10=FourierXP[dia10,{x,q}];


 dia10=dia10/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia10=QEvaluate[I ScaleMu^(2(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia10=QGather[Expand[dia10]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type11[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia11,cc,tr,dot,sunSimplify,contract},


 dia11=Condensate["gg"]*Condensate[{f14, f24}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]], v12]], cqqm[-x, f24]]]*tr[dot[cc[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f11]), {q, x}], hv21]], v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13]] + Condensate["gg"]*Condensate[{f21, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f23]), {q, -x}], v12]], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]]]]*tr[dot[cqqm[x, f11], cc[dot[v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14]] + Condensate["gg"]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f11]), {q, x}], hv21]], v11, cqqm[x, f12]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lorc]], {q, -x}]]], v12]], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]];


 dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia11=FourierXP[dia11,{x,q}];


 dia11=dia11/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia11=QEvaluate[I ScaleMu^(2(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia11=QGather[Expand[dia11]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type12[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia12,cc,tr,dot,sunSimplify,contract},


 dia12=Condensate["gg"]*Condensate[{f14, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], cc[dot[v11, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f12]), {q, x}], hv21]]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]], v12]], hv22, cqqm[-x, f23]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13]] + Condensate["gg"]*Condensate[{f22, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f23]), {q, -x}], v12]], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]]]]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14]] + Condensate["gg"]*Condensate[{f13, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]]]]*tr[dot[cc[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f11]), {q, x}], hv21]], v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14]];


 dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia12=FourierXP[dia12,{x,q}];


 dia12=dia12/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia12=QEvaluate[I ScaleMu^(2(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia12=QGather[Expand[dia12]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type13[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia13,cc,tr,dot,sunSimplify,contract},


 dia13=Condensate["gg"]*Condensate[{f13, f24}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], cc[dot[v11, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f12]), {q, x}], hv21]]]]*tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14]] + Condensate["gg"]*Condensate[{f22, f11}]*contract[tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f12]), {q, x}]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lorc]], {q, -x}]]], v12]], hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]] + Condensate["gg"]*Condensate[{f21, f12}]*contract[tr[dot[cc[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f11]), {q, x}], hv21]], v11, cqqm[x, f12]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lorc]], {q, -x}]]], v12]], hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]];


 dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia13=FourierXP[dia13,{x,q}];


 dia13=dia13/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia13=QEvaluate[I ScaleMu^(2(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia13=QGather[Expand[dia13]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type14[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia14,cc,tr,dot,sunSimplify,contract},


 dia14=Condensate["gg"]*Condensate[{f13, f23}]*contract[tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f24]), {q, -x}]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], hv21]], v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lora]], {q, x}]]]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]] + Condensate["gg"]*Condensate[{f13, f24}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[cqqm[-x, f24], v12]], hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], hv21]], v11, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f12]), {q, x}]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14]] + Condensate["gg"]*Condensate[{f22, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[cqqm[x, f11], hv21]], v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]]]]*tr[dot[cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]], v12]], hv22, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f23]), {q, -x}]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13]];


 dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia14=FourierXP[dia14,{x,q}];


 dia14=dia14/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia14=QEvaluate[I ScaleMu^(2(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia14=QGather[Expand[dia14]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type15[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia15,cc,tr,dot,sunSimplify,contract},


 dia15=Condensate["gg"]*Condensate[{f22, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]], v12]], FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f24]), {q, -x}]]]*tr[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], cc[dot[v11, cqqm[x, f12], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13]] + Condensate["gg"]*Condensate[{f22, f12}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], cc[dot[v11, cqqm[x, f12], hv21]]]]*tr[dot[cc[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f24]), {q, -x}], v12]], hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14]] + Condensate["gg"]*Condensate[{f13, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]]]]*tr[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f11])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f11]), FVD[q, lorc]], {q, x}]]], cc[dot[v11, FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f12]), {q, x}], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14]];


 dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia15=FourierXP[dia15,{x,q}];


 dia15=dia15/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia15=QEvaluate[I ScaleMu^(2(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia15=QGather[Expand[dia15]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*type16[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_]:=Block[{x,q,dia16,cc,tr,dot,sunSimplify,contract},


 dia16=Condensate["gg"]*Condensate[{f21, f11}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cqqm[x, f11], cc[dot[v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]], hv21]]]]*tr[dot[cc[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f24]), {q, -x}], v12]], hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14]] + Condensate["gg"]*Condensate[{f22, f12}]*contract[tr[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f11]), {q, x}], cc[dot[v11, cqqm[x, f12], hv21]]]]*tr[dot[congg[lora, lorb, lorc, lord], cc[dot[QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lorc]], {q, -x}]]], v12]], hv22, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f23])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f23]), FVD[q, lora]], {q, -x}]]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]] + Condensate["gg"]*Condensate[{f13, f23}]*contract[congg[lora, lorb, lorc, lord]*tr[dot[cc[dot[hv22, cqqm[-x, f23], v12]], QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f24])) . GAD[lorb] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f24]), FVD[q, lora]], {q, -x}]]]]]*tr[dot[FourierPX[I*FAD[q]*(GSD[q] + quarkMass[f11]), {q, x}], cc[dot[v11, QSimplify[DiracSimplify[FourierPX[-1/2*(I*FAD[q]*(GSD[q] + quarkMass[f12])) . GAD[lord] . SUNT[sun] . FourDivergence[I*FAD[q]*(GSD[q] + quarkMass[f12]), FVD[q, lorc]], {q, x}]]], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14]];


 dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.{tr->atr,contract->Contract};


 dia16=FourierXP[dia16,{x,q}];


 dia16=dia16/.Power[quarkMass[_],_]->0/.quarkMass[f1_] quarkMass[f2_]->0;


 dia16=QEvaluate[I ScaleMu^(2(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia16=QGather[Expand[dia16]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq
]*)


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
