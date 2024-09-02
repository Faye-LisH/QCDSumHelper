(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
TetraquarkHybridlike2::usage="TetraquarkHybridlike2[q_,j1_,j2_] gives leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents. "

TetraquarkHybridlike2::inderr="Dummy indices conflict!"

TetraquarkHybridlike2::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`TetraquarkHybridlike2`"]

Options[TetraquarkHybridlike2] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True,
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

prop2[l_,lorb_,lorz_,lor1_,lor2_]=FeynAmpDenominator[PropagatorDenominator[Momentum[l, D], 0]]*Pair[LorentzIndex[lor1, D], LorentzIndex[lor2, D]]*Pair[LorentzIndex[lorb, D], Momentum[l, D]]*
  Pair[LorentzIndex[lorz, D], Momentum[l, D]]*Pair[Momentum[l, D], Momentum[l, D]]^(-1 + D/2)*qfact1[-((I^(1 + D)*2^(-1 - D)*gStrong*qGamma[1 - D/2]*qGamma[D/2]^2)/(Pi^(D/2)*qGamma[D]))] + 
 FeynAmpDenominator[PropagatorDenominator[Momentum[l, D], 0]]*Pair[LorentzIndex[lor1, D], LorentzIndex[lor2, D]]*Pair[LorentzIndex[lorb, D], LorentzIndex[lorz, D]]*
  Pair[Momentum[l, D], Momentum[l, D]]^(D/2)*qfact1[(I^(1 + D)*2^(-1 - D)*gStrong*qGamma[1 - D/2]*qGamma[D/2]^2)/(Pi^(D/2)*qGamma[D])] + 
 FeynAmpDenominator[PropagatorDenominator[Momentum[l, D], 0]]*Pair[LorentzIndex[lor1, D], Momentum[l, D]]*Pair[LorentzIndex[lor2, D], Momentum[l, D]]*
  Pair[LorentzIndex[lorb, D], LorentzIndex[lorz, D]]*Pair[Momentum[l, D], Momentum[l, D]]^(-1 + D/2)*qfact1[-((I^(1 + D)*gStrong*qGamma[2 - D/2]*qGamma[D/2]^2)/(2^D*Pi^(D/2)*qGamma[D]))] + 
 FeynAmpDenominator[PropagatorDenominator[Momentum[l, D], 0]]*Pair[LorentzIndex[lor1, D], Momentum[l, D]]*Pair[LorentzIndex[lor2, D], Momentum[l, D]]*Pair[LorentzIndex[lorb, D], Momentum[l, D]]*
  Pair[LorentzIndex[lorz, D], Momentum[l, D]]*Pair[Momentum[l, D], Momentum[l, D]]^(-2 + D/2)*qfact1[(I^(1 + D)*gStrong*qGamma[2 - D/2]*qGamma[D/2]^2)/(2^D*Pi^(D/2)*qGamma[D])];
				
								
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


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
TetraquarkHybridlike2[qq_,factor1_ current1_FourquarkCurrent, factor2_ current2_HybirdCurrent,ops:OptionsPattern[]]:=factor1 factor2 TetraquarkHybridlike2[qq,current1,current2,ops]/;(FreeQ[factor1,Current|HybirdCurrent]&&FreeQ[factor2,Current|HybirdCurrent])
TetraquarkHybridlike2[qq_,factor1_ current1_FourquarkCurrent, current2_HybirdCurrent,ops:OptionsPattern[]]:=factor1 TetraquarkHybridlike2[qq,current1,current2,ops]/;FreeQ[factor1,Current|HybirdCurrent]
TetraquarkHybridlike2[qq_,current1_FourquarkCurrent, factor2_ current2_HybirdCurrent,ops:OptionsPattern[]]:=factor2 TetraquarkHybridlike2[qq,current1,current2,ops]/;FreeQ[factor2,Current|HybirdCurrent]


TetraquarkHybridlike2[qq_,current1_FourquarkCurrent,current2_HybirdCurrent,ops:OptionsPattern[]]:=Block[{c21,c22,f21,f22,sut,gs,gv},
f21=current2[[1]];
c21=current2[[2]];
sut=SUNTF[current2[[3,1,3]],current2[[2]],current2[[6]]];
gs=GluonStrength@@current2[[3]];
gv=current2[[4]];
f22=current2[[5]];
c22=current2[[6]];


TetraquarkHybridlike2[qq,current1,{f21,c21,sut,gs,gv,f22,c22},ops]
]



(* ::Code::Initialization::Plain:: *)
TetraquarkHybridlike2[qq_,_,0,ops:OptionsPattern[]]=0
TetraquarkHybridlike2[qq_,_,{f3_,c3_,sut_,gs_,0,f4_,c4_},ops:OptionsPattern[]]=0
TetraquarkHybridlike2[qq_,factor1_ current1_FourquarkCurrent,{f3_,c3_,sut_,gs_,gv_,f4_,c4_},ops:OptionsPattern[]]:=factor1 TetraquarkHybridlike2[qq,current1,{f3,c3,sut,gs,gv,f4,c4},ops]/;FreeQ[factor1,Current|HybirdCurrent]


(* ParallelSubmit is trick, ParallelSubmit[{...,f21,f22,c21,c22,...}, expr] will submit f21,f22,... to expr before sumbiting it for evaluation,
 however, here if writing the function as TetraquarkHybridlike2[qq_,current1_FourquarkCurrent,{f21_,c21_,sut_,gs_,v21_,f22_,c22_},OptionsPattern[]]:=..., the f21,f22,... are not values.  *)
TetraquarkHybridlike2[qq_,current1_FourquarkCurrent,{f3_,c3_,sut_,gs_,gv_,f4_,c4_},OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,trs,dia,x,q,sun=sut[[1,1]],lorb,sign,f21=f3,c21=c3,v21=gv,f22=f4,c22=c4,hv21,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,order=OptionValue[EpsOrder],holdf=OptionValue[HoldFlavor],lorz,lor1,lor2,strategy=OptionValue[Strategy],diagrams,atr,pole=OptionValue[Pole]},


(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[TetraquarkHybridlike2::curerr];
	Abort[]
];


(*--- check dummy indices in input ---*)
tmp2=Length[#]&/@Gather[({c11,c12,c13,c14,c21,c22}//FCI)/.SUNFIndex[aa_]:>aa];

If[DummyCheck[v11,v12,v21,0]==0!FreeQ[tmp2,nn_Integer/;nn>2],
	Message[Tetraquarkggg::inderr];
	Abort[]
];


(*-------------------------------*)
(* B A^+ B *)
hv21=v21//ComplexConjugate;
pole=ComplexConjugate[pole];
gs/.{GluonStrength[{lo1_,lo2_,col_},lo2_]:>(lorb=lo1;sign=-1;null),GluonStrength[{lo1_,lo2_,col_},lo1_]:>(lorb=lo2;sign=1;null)};
(*------------------------------------------------------------------*)

If[ToLowerCase[ToString[strategy]]=="fourier",
	diagrams={typetx1}
,
	If[ToLowerCase[ToString[strategy]]=="integratefirst",
		diagrams={typep1}
	,
		diagrams={typetp1}
	]
];


If[OptionValue[AutoNDR]===True,atr=TR5,atr=TR];


(*---------------------------------------------------*)
If[OptionValue[Parallelized]===True,

	
	tmp=Plus@@WaitAll[ParallelSubmit[{hv21,sun,lorb,lorz,lor1,lor2,v11,v12,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,f21,f22,holdf,order,atr,pole},
									#[qq,hv21,sun,lorb,lorz,lor1,lor2,v11,v12,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,f21,f22,holdf,order,atr,pole]]&/@diagrams];
	
	sign QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,
	If[ToString[OptionValue[Parallelized]]==="External",
	

		
		tmp=sign Plus@@(ParallelSubmit[{hv21,sun,lorb,lorz,lor1,lor2,v11,v12,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,f21,f22,holdf,order,atr,pole},
									#[qq,hv21,sun,lorb,lorz,lor1,lor2,v11,v12,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,f21,f22,holdf,order,atr,pole]]&/@diagrams)
		(* !!! do not touch tmp before WaitAll[tmp], any function may modifiy the things in EvaluationObject[...] should not be used here. !!! *)
	
	,
	
		
		tmp=Plus@@Map[#[qq,hv21,sun,lorb,lorz,lor1,lor2,v11,v12,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,f21,f22,holdf,order,atr,pole]&,diagrams];
		
		sign QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]


]


(* ::Code::Initialization::Plain:: *)
(*--- The following are generated by algorithm ---*)


(* ::Input::Initialization:: *)
typetx1[qq_,hv21_,sun_,lorb_,lorz_,lor1_,lor2_,v11_,v12_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,f21_,f22_,holdf_,order_,atr_,pole_]:=Block[{x,q,dia1},


dia1=-(ColorDelta[c12, c22]*ColorDelta[c21, c14]*Contract[SUNTF[sun, c22, c21]*atr[fcc[GAD[lor1] . GAD[lorz] . GAD[lor2] . v12] . v11 . xprop[x] . hv21 . xprop[-x]]*xprop2[x, lorb, lorz, lor1, lor2]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f21, f14]*SUNTF[sun, c11, c13]) - ColorDelta[c12, c22]*ColorDelta[c21, c13]*Contract[SUNTF[sun, c22, c21]*atr[fcc[GAD[lor1] . GAD[lorz] . GAD[lor2]] . v11 . xprop[x] . hv21 . xprop[-x] . v12]*xprop2[x, lorb, lorz, lor1, lor2]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f21, f13]*SUNTF[sun, c11, c14] - ColorDelta[c11, c22]*ColorDelta[c21, c14]*Contract[SUNTF[sun, c22, c21]*atr[fcc[v11 . GAD[lor1] . GAD[lorz] . GAD[lor2] . v12] . xprop[x] . hv21 . xprop[-x]]*xprop2[x, lorb, lorz, lor1, lor2]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f21, f14]*SUNTF[sun, c12, c13] - ColorDelta[c11, c22]*ColorDelta[c21, c13]*Contract[SUNTF[sun, c22, c21]*atr[xprop[x] . hv21 . xprop[-x] . v12 . fcc[v11 . GAD[lor1] . GAD[lorz] . GAD[lor2]]]*xprop2[x, lorb, lorz, lor1, lor2]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f21, f13]*SUNTF[sun, c12, c14];


dia1=FourierXP[dia1,{x,q}]//SUNSimplify;

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1 pole,q,HoldFlavor->holdf,EpsOrder->order,Parallelized->False]/.q->qq
]



typep1[qq_,hv21_,sun_,lorb_,lorz_,lor1_,lor2_,v11_,v12_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,f21_,f22_,holdf_,order_,atr_,pole_]:=Block[{x,q,k1,k2,l,k,dia1},


dia1=-(ColorDelta[c12, c22]*ColorDelta[c21, c14]*Contract[prop2[l, lorb, lorz, lor1, lor2]*SUNTF[sun, c22, c21]*atr[fcc[GAD[lor1] . GAD[lorz] . GAD[lor2] . v12] . v11 . prop[-k + q] . hv21 . prop[-k + l]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f21, f14]*SUNTF[sun, c11, c13]) - ColorDelta[c12, c22]*ColorDelta[c21, c13]*Contract[prop2[l, lorb, lorz, lor1, lor2]*SUNTF[sun, c22, c21]*atr[fcc[GAD[lor1] . GAD[lorz] . GAD[lor2]] . v11 . prop[-k + q] . hv21 . prop[-k + l] . v12]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f21, f13]*SUNTF[sun, c11, c14] - ColorDelta[c11, c22]*ColorDelta[c21, c14]*Contract[prop2[l, lorb, lorz, lor1, lor2]*SUNTF[sun, c22, c21]*atr[fcc[v11 . GAD[lor1] . GAD[lorz] . GAD[lor2] . v12] . prop[-k + q] . hv21 . prop[-k + l]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f21, f14]*SUNTF[sun, c12, c13] - ColorDelta[c11, c22]*ColorDelta[c21, c13]*Contract[prop2[l, lorb, lorz, lor1, lor2]*SUNTF[sun, c22, c21]*atr[prop[-k + q] . hv21 . prop[-k + l] . v12 . fcc[v11 . GAD[lor1] . GAD[lorz] . GAD[lor2]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f21, f13]*SUNTF[sun, c12, c14];


dia1=IntegrateP[dia1,{l,k}]//SUNSimplify;

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1 pole,q,HoldFlavor->holdf,EpsOrder->order,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
typetp1[qq_,hv21_,sun_,lorb_,lorz_,lor1_,lor2_,v11_,v12_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,f21_,f22_,simp_,holdf_,order_,atr_,pole_]:=Block[{x,q,k1,k2,l,k,dia1},


dia1=-(ColorDelta[c12, c22]*ColorDelta[c21, c14]*Contract[gStrong*FAD[l]*(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l])*SUNTF[sun, c22, c21]*atr[fcc[prop[k1 + l] . GAD[lorz] . prop[k1] . v12] . v11 . prop[-k + q] . hv21 . prop[-k + l]]]*FlavorDelta[f11, f13]*FlavorDelta[f12, f22]*FlavorDelta[f21, f14]*SUNTF[sun, c11, c13]) - ColorDelta[c12, c22]*ColorDelta[c21, c13]*Contract[gStrong*FAD[l]*(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l])*SUNTF[sun, c22, c21]*atr[fcc[prop[k1 + l] . GAD[lorz] . prop[k1]] . v11 . prop[-k + q] . hv21 . prop[-k + l] . v12]]*FlavorDelta[f11, f14]*FlavorDelta[f12, f22]*FlavorDelta[f21, f13]*SUNTF[sun, c11, c14] - ColorDelta[c11, c22]*ColorDelta[c21, c14]*Contract[gStrong*FAD[l]*(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l])*SUNTF[sun, c22, c21]*atr[fcc[v11 . prop[k1 + l] . GAD[lorz] . prop[k1] . v12] . prop[-k + q] . hv21 . prop[-k + l]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f13]*FlavorDelta[f21, f14]*SUNTF[sun, c12, c13] - ColorDelta[c11, c22]*ColorDelta[c21, c13]*Contract[gStrong*FAD[l]*(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l])*SUNTF[sun, c22, c21]*atr[prop[-k + q] . hv21 . prop[-k + l] . v12 . fcc[v11 . prop[k1 + l] . GAD[lorz] . prop[k1]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f14]*FlavorDelta[f21, f13]*SUNTF[sun, c12, c14];


dia1=IntegrateP[dia1,{k1,l,k}]//SUNSimplify;

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1 pole,q,HoldFlavor->holdf,EpsOrder->order,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
