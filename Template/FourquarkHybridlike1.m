(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
FourquarkHybridlike1::usage="FourquarkHybridlike1[q_,j1_,j2_] gives leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents. "

FourquarkHybridlike1::inderr="Dummy indices conflict!"

FourquarkHybridlike1::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`FourquarkHybridlike1`"]

Options[FourquarkHybridlike1] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->False,
	EpsOrder->0,
	Strategy->"Fourier",
	Pole->1
}


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

prop[q_] = I GSD[q]FAD[q];

(* xprop[x_]=FourierPX[prop[q],{q,x}]; *)
xprop[x_] = 1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];
			
			
(* DG - quark loop : D^loru G_{loru lorb} GAD[lor1, lorz, lor2]*)			
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


FourquarkHybridlike1[qq_,0,{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]=0
FourquarkHybridlike1[qq_,{a1_,sut_,gs_,0,b1_},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]=0


FourquarkHybridlike1[qq_,{a1_,sut_,gs_,gv_,b1_},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,trs,dia,x,q,sun=sut[[1,1]],lorb,sign,v1=sut gv,
hv3,hv4,order=OptionValue[EpsOrder],holdf=OptionValue[HoldFlavor],strategy=OptionValue[Strategy],diagrams,atr,plus,waitAll,parallelSubmit,pole=OptionValue[Pole]},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,0,v3,v4]==0,
	Message[FourquarkAlphasType2::inderr];
	Abort[]
];



(*-------------------------------*)
If[OptionValue[AutoNDR]===True,
	atr=TR5
	,
	atr=TR
];


(*-------------------------------*)
(* B A^+ B *)
hv3=v3//ComplexConjugate;
hv4=v4//ComplexConjugate;

gs/.{GluonStrength[{lo1_,lo2_,col_},lo2_]:>(lorb=lo1;sign=-1;null),GluonStrength[{lo1_,lo2_,col_},lo1_]:>(lorb=lo2;sign=1;null)};


(*------------------------------------------------------------------*)

If[ToLowerCase[ToString[strategy]]=="fourier",
	diagrams={typetx1}
,
	diagrams={typetp1}
];



If[OptionValue[AutoNDR]===True,atr=TR5,atr=TR];


(*---------------------------------------------------*)


If[OptionValue[Parallelized]===True,

	DistributeDefinitions[qq,a1,b1,v3,v4,a2,b2,c2,d2];
	tmp= Plus@@WaitAll[ParallelSubmit[{v1,hv3,hv4,sun,lorb,holdf,order,atr,pole},
									#[qq,{a1,v1,b1},{hv3,hv4,{a2,b2,c2,d2}},sun,lorb,holdf,atr,order,pole]]&/@diagrams];

	sign QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	
	

,
	If[ToString[OptionValue[Parallelized]]==="External",
	
		DistributeDefinitions[qq,a1,b1,v3,v4,a2,b2,c2,d2];
		tmp= sign Plus@@(ParallelSubmit[{v1,hv3,hv4,sun,lorb,holdf,order,atr,pole},
									#[qq,{a1,v1,b1},{hv3,hv4,{a2,b2,c2,d2}},sun,lorb,holdf,atr,order,pole]]&/@diagrams)
		(* !!! do not touch tmp before WaitAll[tmp], any function may modifiy the things in EvaluationObject[...] should not be used here. !!! *)
	
	,
	

		tmp= Plus@@Map[#[qq,{a1,v1,b1},{hv3,hv4,{a2,b2,c2,d2}},sun,lorb,holdf,atr,order,pole]&,diagrams];
		sign QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]


]


(* ::Code::Initialization::Plain:: *)
(* The following are generated by algorithem *)



xtype1[qq_,{a1_,v1_,b1_},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},sun_,lorb_,holdf_,atr_,order_,pole_]:=Block[{x,q,dia1,str,tr,dot,contract},


dia1=contract[tr[str[dot[hv3, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], xprop2[x, lorb, lorz, lor1, lor2]]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1] + contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x], xprop2[x, lorb, lorz, lor1, lor2]]]]*tr[str[dot[hv4, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2] - FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2], xprop2[x, lorb, lorz, lor1, lor2]]]] - FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, GAD[lor1], GAD[lorz], SUNT[sun], GAD[lor2], xprop2[x, lorb, lorz, lor1, lor2]]]];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace,tr->atr,contract->Contract};


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(3(4-D))pole dia1,q,HoldFlavor->holdf,EpsOrder->order,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
typetx1[qq_,{a1_,v1_,b1_},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},sun_,lorb_,holdf_,atr_,order_,pole_]:=Block[{x,q,dia1,lor1,lor2,lorz},


dia1=-(SUNTrace[atr[xprop2[x, lorb, lorz, lor1, lor2] . hv4 . xprop[-x] . v1 . xprop[x] . hv3 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]) + SUNTrace[atr[hv4 . xprop[-x] . v1 . xprop[x]]]*SUNTrace[atr[xprop2[x, lorb, lorz, lor1, lor2] . hv3 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1] - SUNTrace[atr[xprop2[x, lorb, lorz, lor1, lor2] . hv3 . xprop[-x] . v1 . xprop[x] . hv4 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2] + SUNTrace[atr[hv3 . xprop[-x] . v1 . xprop[x]]]*SUNTrace[atr[xprop2[x, lorb, lorz, lor1, lor2] . hv4 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2];


dia1=FourierXP[dia1,{x,q}]//SUNSimplify;

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1 pole,q,HoldFlavor->holdf,Parallelized->False,EpsOrder->order]/.q->qq
]


(* ::Input::Initialization:: *)
typetp1[qq_,{a1_,v1_,b1_},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},sun_,lorb_,holdf_,atr_,order_,pole_]:=Block[{k1,l,k,q,dia1,lor1,lor2,lorz},


dia1=-(gStrong*SUNTrace[atr[(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l]) . hv4 . prop[k - q] . v1 . prop[k - l] . hv3 . prop[k1] . GAD[lorz] . SUNT[sun] . prop[k1 + l]]]*FAD[l]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]) + gStrong*SUNTrace[atr[hv4 . prop[k - q] . v1 . prop[k - l]]]*SUNTrace[atr[(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l]) . hv3 . prop[k1] . GAD[lorz] . SUNT[sun] . prop[k1 + l]]]*FAD[l]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1] - gStrong*SUNTrace[atr[(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l]) . hv3 . prop[k - q] . v1 . prop[k - l] . hv4 . prop[k1] . GAD[lorz] . SUNT[sun] . prop[k1 + l]]]*FAD[l]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2] + gStrong*SUNTrace[atr[hv3 . prop[k - q] . v1 . prop[k - l]]]*SUNTrace[atr[(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l]) . hv4 . prop[k1] . GAD[lorz] . SUNT[sun] . prop[k1 + l]]]*FAD[l]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2];


dia1=IntegrateP[dia1,{k1,l,k}]//SUNSimplify;

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1 pole,q,HoldFlavor->holdf,Parallelized->False,EpsOrder->order]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
