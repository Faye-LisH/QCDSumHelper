(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
FourquarkHybridlike2::usage="FourquarkHybridlike2[q_,j1_,j2_] gives leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents. "

FourquarkHybridlike2::inderr="Dummy indices conflict!"

FourquarkHybridlike2::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`FourquarkHybridlike2`"]

Options[FourquarkHybridlike2] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->False,
	EpsOrder->0,
	Strategy->"Fourier",
	Pole->1
}


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(* propagators *)

prop[q_] = I GSD[q]FAD[q];

(* xprop[x_]=FourierPX[prop[q],{q,x}]; *)
xprop[x_] = 1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];
			
			
(* DG - quark loop *)	
			
								
xprop2[x_,lorb_,lorz_,lor1_,lor2_]=(Pair[LorentzIndex[lor1, D], LorentzIndex[lor2, D]]*Pair[LorentzIndex[lorb, D], Momentum[x, D]]*Pair[LorentzIndex[lorz, D], Momentum[x, D]]*
   qfact1[(qfact2[-gStrong/(8*E^(I*D*Pi)*Pi^D)]*qGamma[1 - D/2]*qGamma[D/2]^2)/qGamma[3 - D/2]])/Pair[Momentum[x, D], Momentum[x, D]]^D + 
 (Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lorz, D]]*Pair[LorentzIndex[lorb, D], Momentum[x, D]]*
   qfact1[-((-1)^(1 - D)*gStrong*qGamma[2 - D/2]*qGamma[D/2]^2)/(8*Pi^D*qGamma[3 - D/2])])/Pair[Momentum[x, D], Momentum[x, D]]^D + 
 (Pair[LorentzIndex[lor1, D], LorentzIndex[lorz, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lorb, D], Momentum[x, D]]*
   qfact1[-((-1)^(1 - D)*gStrong*qGamma[2 - D/2]*qGamma[D/2]^2)/(8*Pi^D*qGamma[3 - D/2])])/Pair[Momentum[x, D], Momentum[x, D]]^D + 
 (Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lorb, D]]*Pair[LorentzIndex[lorz, D], Momentum[x, D]]*
   qfact1[-((-1)^(1 - D)*gStrong*qGamma[2 - D/2]*qGamma[D/2]^2)/(8*Pi^D*qGamma[3 - D/2])])/Pair[Momentum[x, D], Momentum[x, D]]^D + 
 (Pair[LorentzIndex[lor1, D], LorentzIndex[lorb, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lorz, D], Momentum[x, D]]*
   qfact1[-((-1)^(1 - D)*gStrong*qGamma[2 - D/2]*qGamma[D/2]^2)/(8*Pi^D*qGamma[3 - D/2])])/Pair[Momentum[x, D], Momentum[x, D]]^D + 
 (Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lorb, D], LorentzIndex[lorz, D]]*
   qfact1[(I^(2 + D)*gStrong*qGamma[D/2]^2)/(4*E^(((3*I)/2)*D*Pi)*Pi^D) - ((-1)^(1 - D)*gStrong*qGamma[2 - D/2]*qGamma[D/2]^2)/(8*Pi^D*qGamma[3 - D/2])])/Pair[Momentum[x, D], Momentum[x, D]]^D + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lor2, D]]*Pair[LorentzIndex[lorb, D], LorentzIndex[lorz, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D)*
  qfact1[(qfact2[gStrong/(16*E^(I*D*Pi)*Pi^D)]*qGamma[1 - D/2]*qGamma[-1 + D]*qGamma[D/2]^2)/(qGamma[3 - D/2]*qGamma[D])] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lorz, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lorb, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D)*
  qfact1[((-1)^(1 - D)*gStrong*qGamma[2 - D/2]*qGamma[-1 + D]*qGamma[D/2]^2)/(16*Pi^D*qGamma[3 - D/2]*qGamma[D])] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lorb, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lorz, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D)*
  qfact1[((-1)^(1 - D)*gStrong*qGamma[2 - D/2]*qGamma[-1 + D]*qGamma[D/2]^2)/(16*Pi^D*qGamma[3 - D/2]*qGamma[D])] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lorb, D], Momentum[x, D]]*Pair[LorentzIndex[lorz, D], Momentum[x, D]]*
  Pair[Momentum[x, D], Momentum[x, D]]^(-1 - D)*qfact1[((-1)^(1 - D)*gStrong*qGamma[2 - D/2]*qGamma[D/2]^2*qGamma[1 + D])/(4*Pi^D*qGamma[3 - D/2]*qGamma[D])];


FourquarkHybridlike2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},0,OptionsPattern[]]=0
FourquarkHybridlike2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{a2_,sut_,gs_,0,b2_},OptionsPattern[]]=0


FourquarkHybridlike2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{a2_,sut_,gs_,gv_,b2_},OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,trs,dia,x,q,sun=sut[[1,1]],lorb,sign,v3=sut gv,
hv3,hv4,order=OptionValue[EpsOrder],holdf=OptionValue[HoldFlavor],strategy=OptionValue[Strategy],diagrams,atr,pole=OptionValue[Pole]},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,0]==0,
	Message[FourquarkAlphasType2::inderr];
	Abort[]
];


If[OptionValue[AutoNDR]===True,
	atr=TR5
	,
	atr=TR
];


(*-------------------------------*)
(* B A^+ B *)
hv3=v3//ComplexConjugate;
pole=ComplexConjugate[pole];
gs/.{GluonStrength[{lo1_,lo2_,col_},lo2_]:>(lorb=lo1;sign=-1;null),GluonStrength[{lo1_,lo2_,col_},lo1_]:>(lorb=lo2;sign=1;null)};
(*------------------------------------------------------------------*)

If[ToLowerCase[ToString[strategy]]=="fourier",
	diagrams={xtype1}
,
	diagrams={typetp1}
];



If[OptionValue[AutoNDR]===True,atr=TR5,atr=TR];


(*---------------------------------------------------*)


If[OptionValue[Parallelized]===True,

	DistributeDefinitions[qq,v1,v2,a1,b1,c1,d1,a2,b2];
	tmp= Plus@@WaitAll[ParallelSubmit[{hv3,sun,lorb,holdf,atr,order,pole},
									#[qq,{v1,v2,{a1,b1,c1,d1}},{a2,hv3,b2},sun,lorb,holdf,atr,order,pole]]&/@diagrams];

	sign QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	
	
,
	If[ToString[OptionValue[Parallelized]]==="External",
	
		DistributeDefinitions[qq,v1,v2,a1,b1,c1,d1,a2,b2];
		tmp= sign Plus@@(ParallelSubmit[{hv3,sun,lorb,holdf,atr,order,pole},
									#[qq,{v1,v2,{a1,b1,c1,d1}},{a2,hv3,b2},sun,lorb,holdf,atr,order,pole]]&/@diagrams)
		(* !!! do not touch tmp before WaitAll[tmp], any function may modifiy the things in EvaluationObject[...] should not be used here. !!! *)
	
	,
	

		tmp= Plus@@Map[#[qq,{v1,v2,{a1,b1,c1,d1}},{a2,hv3,b2},sun,lorb,holdf,atr,order,pole]&,diagrams];
		sign QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]


]


(* ::Code::Initialization::Plain:: *)
(*--- The following are generated by algorithm ---*)



xtype1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{a2_,hv3_,b2_},sun_,lorb_,holdf_,atr_,order_,pole_]:=Block[{x,q,dia1,str,tr,dot,contract},


dia1=contract[tr[str[dot[hv3, xprop[-x], v2, xprop[x], xprop2[x, lorb, lorz, lor1, lor2]]]]*tr[str[dot[v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[d1, b2] + contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x], xprop2[x, lorb, lorz, lor1, lor2]]]]*tr[str[dot[v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[d1, c1] - FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[d1, b2]*tr[str[dot[hv3, xprop[-x], v1, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2], v2, xprop[x], xprop2[x, lorb, lorz, lor1, lor2]]]] - FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[d1, a1]*tr[str[dot[v2, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2], v1, xprop[x], hv3, xprop[-x], xprop2[x, lorb, lorz, lor1, lor2]]]];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(3(4-D))pole dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
typetx1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{a2_,hv3_,b2_},sun_,lorb_,holdf_,atr_,order_,pole_]:=Block[{x,q,dia1,lorz},


dia1=-(SUNTrace[atr[xprop2[x, lorb, lorz, lor1, lor2] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . xprop[-x]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[d1, a1]) + SUNTrace[atr[hv3 . xprop[-x] . v2 . xprop[x]]]*SUNTrace[atr[xprop2[x, lorb, lorz, lor1, lor2] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[d1, b2] - SUNTrace[atr[xprop2[x, lorb, lorz, lor1, lor2] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[d1, b2] + SUNTrace[atr[hv3 . xprop[-x] . v1 . xprop[x]]]*SUNTrace[atr[xprop2[x, lorb, lorz, lor1, lor2] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[d1, c1];


dia1=FourierXP[dia1,{x,q}]//SUNSimplify;

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1 pole,q,HoldFlavor->holdf,EpsOrder->order,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
typetp1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{a2_,hv3_,b2_},sun_,lorb_,holdf_,atr_,order_,pole_]:=Block[{k1,l,k,dia1,lorz},


dia1=-(gStrong*SUNTrace[atr[(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l]) . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1] . v1 . prop[-k + q] . hv3 . prop[-k + l]]]*FAD[l]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[d1, a1]) + dia[SUNTrace[atr[(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l]) . v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1]]], SUNTrace[atr[hv3 . prop[k - q] . v2 . prop[k - l]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[d1, b2] - gStrong*SUNTrace[atr[(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l]) . hv3 . prop[k - q] . v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1] . v2 . prop[k - l]]]*FAD[l]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[d1, b2] + dia[SUNTrace[atr[hv3 . prop[k - q] . v1 . prop[k - l]]], SUNTrace[atr[(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l]) . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[d1, c1];


dia1=IntegrateP[dia1,{k1,l,k}]//SUNSimplify;

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1 pole,q,HoldFlavor->holdf,EpsOrder->order,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
