(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
TetraquarkHybridlike1::usage="TetraquarkHybridlike1[q_,j1_,j2_] gives leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents. "

TetraquarkHybridlike1::inderr="Dummy indices conflict!"

TetraquarkHybridlike1::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`TetraquarkHybridlike1`"]

Options[TetraquarkHybridlike1] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True,
	EpsOrder->0,
	Strategy->"Fourier",
	Pole->1
}


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

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


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)

fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(*------------------------------------------------------------------*)
TetraquarkHybridlike1[qq_,factor1_ current1_HybirdCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 TetraquarkHybridlike1[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current|HybirdCurrent])
TetraquarkHybridlike1[qq_,current1_HybirdCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor2 TetraquarkHybridlike1[qq,current1,current2,ops]/;FreeQ[factor2,Current|HybirdCurrent]
TetraquarkHybridlike1[qq_,factor1_ current1_HybirdCurrent,current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 TetraquarkHybridlike1[qq,current1,current2,ops]/;FreeQ[factor1,Current|HybirdCurrent]


TetraquarkHybridlike1[qq_,current1_HybirdCurrent,current2_FourquarkCurrent,ops:OptionsPattern[]]:=Block[{c11,c12,f11,f12,sut,gs,gv,lorb,tmpc1},
f11=current1[[1]];
c11=current1[[2]];
sut=SUNTF[current1[[3,1,3]],current1[[2]],current1[[6]]];
gs=GluonStrength@@current1[[3]];
gv=current1[[4]];
f12=current1[[5]];
c12=current1[[6]];


TetraquarkHybridlike1[qq,{f11,c11,sut,gs,gv,f12,c12},current2,ops]
]



(* ::Code::Initialization::Plain:: *)
TetraquarkHybridlike1[qq_,0,_,ops:OptionsPattern[]]=0
TetraquarkHybridlike1[qq_,{f1_,c1_,sut_,gs_,0,f2_,c2_},_,ops:OptionsPattern[]]=0
TetraquarkHybridlike1[qq_,{f1_,c1_,sut_,gs_,gv_,f2_,c2_},factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor2 TetraquarkHybridlike1[qq,{f1,c1,sut,gs,gv,f2,c2},current2,ops]/;FreeQ[factor2,Current|HybirdCurrent]


TetraquarkHybridlike1[qq_,{f1_,c1_,sut_,gs_,gv_,f2_,c2_},current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,trs,dia,x,q,sun=sut[[1,1]],lorb,sign,f11=f1,c11=c1,v11=gv,f12=f2,c12=c2,
v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,order=OptionValue[EpsOrder],holdf=OptionValue[HoldFlavor],lorz,lor1,lor2,strategy=OptionValue[Strategy],diagrams,atr,pole=OptionValue[Pole]},


(*--- get the falovrs, vertices, and color indices ---*)
tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[TetraquarkHybridlike1::curerr];
	Abort[]
];

(*--- check dummy indices in input ---*)
tmp2=Length[#]&/@Gather[({c11,c12,c21,c22,c23,c24}//FCI)/.SUNFIndex[aa_]:>aa];

If[DummyCheck[v11,0,v21,v22]==0||!FreeQ[tmp2,nn_Integer/;nn>2],
	Message[Tetraquarkggg::inderr];
	Abort[]
];


(*-------------------------------*)
(* B A^+ B *)
hv21=v21//ComplexConjugate;
hv22=v22//ComplexConjugate;

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


	tmp=Plus@@WaitAll[ParallelSubmit[{v11,sun,lorb,lorz,lor1,lor2,hv21,hv22,c11,c12,f11,f12,c21,c22,c23,c24,f21,f22,f23,f24,holdf,order,atr,pole},
									#[qq,v11,sun,lorb,lorz,lor1,lor2,hv21,hv22,c11,c12,f11,f12,c21,c22,c23,c24,f21,f22,f23,f24,holdf,order,atr,pole]]&/@diagrams];
	
	sign QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	
	

,
	If[ToString[OptionValue[Parallelized]]==="External",
	
		
		tmp=sign Plus@@(ParallelSubmit[{v11,sun,lorb,lorz,lor1,lor2,hv21,hv22,c11,c12,f11,f12,c21,c22,c23,c24,f21,f22,f23,f24,holdf,order,atr,pole},
								#[qq,v11,sun,lorb,lorz,lor1,lor2,hv21,hv22,c11,c12,f11,f12,c21,c22,c23,c24,f21,f22,f23,f24,holdf,order,atr,pole]]&/@diagrams)
		(* !!! do not touch tmp before WaitAll[tmp], any function may modifiy the things in EvaluationObject[...] should not be used here. !!! *)
	
	,
	
		
		tmp= Plus@@Map[#[qq,v11,sun,lorb,lorz,lor1,lor2,hv21,hv22,c11,c12,f11,f12,c21,c22,c23,c24,f21,f22,f23,f24,holdf,order,atr,pole]&,diagrams];
		
		sign QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]


]


(* ::Code::Initialization::Plain:: *)
(* The following are generated by algorithem *)


typetp1[qq_,v11_,sun_,lorb_,lorz_,lor1_,lor2_,hv21_,hv22_,c11_,c12_,f11_,f12_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,order_,atr_,pole_]:=Block[{x,q,k1,k2,l,k,dia1},


dia1=-(ColorDelta[c12, c22]*ColorDelta[c24, c11]*Contract[gStrong*FAD[l]*(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l])*SUNTF[sun, c11, c12]*atr[fcc[prop[-k + q] . v11 . prop[-k + l] . hv21] . hv22 . prop[k1] . GAD[lorz] . prop[k1 + l]]]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f11]*SUNTF[sun, c23, c21]) - ColorDelta[c12, c21]*ColorDelta[c24, c11]*Contract[gStrong*FAD[l]*(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l])*SUNTF[sun, c11, c12]*atr[prop[-k + q] . v11 . prop[-k + l] . fcc[hv22 . prop[k1] . GAD[lorz] . prop[k1 + l] . hv21]]]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f11]*SUNTF[sun, c23, c22] - ColorDelta[c12, c22]*ColorDelta[c23, c11]*Contract[gStrong*FAD[l]*(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l])*SUNTF[sun, c11, c12]*atr[fcc[hv22 . prop[-k + q] . v11 . prop[-k + l] . hv21] . prop[k1] . GAD[lorz] . prop[k1 + l]]]*FlavorDelta[f12, f22]*FlavorDelta[f23, f11]*FlavorDelta[f24, f21]*SUNTF[sun, c24, c21] - ColorDelta[c12, c21]*ColorDelta[c23, c11]*Contract[gStrong*FAD[l]*(FVD[l, lorb]*FVD[l, lorz] - MTD[lorb, lorz]*SPD[l])*SUNTF[sun, c11, c12]*atr[prop[k1] . GAD[lorz] . prop[k1 + l] . hv21 . fcc[hv22 . prop[-k + q] . v11 . prop[-k + l]]]]*FlavorDelta[f12, f21]*FlavorDelta[f23, f11]*FlavorDelta[f24, f22]*SUNTF[sun, c24, c22];


dia1=IntegrateP[dia1,{k1,k,l}]//SUNSimplify;

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1 pole,q,HoldFlavor->holdf,EpsOrder->order,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
typetx1[qq_,v11_,sun_,lorb_,lorz_,lor1_,lor2_,hv21_,hv22_,c11_,c12_,f11_,f12_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,order_,atr_,pole_]:=Block[{x,q,dia1},


dia1=-(ColorDelta[c12, c22]*ColorDelta[c24, c11]*Contract[SUNTF[sun, c11, c12]*atr[fcc[xprop[-x] . v11 . xprop[x] . hv21] . hv22 . GAD[lor1] . GAD[lorz] . GAD[lor2]]*xprop2[x, lorb, lorz, lor1, lor2]]*FlavorDelta[f12, f22]*FlavorDelta[f23, f21]*FlavorDelta[f24, f11]*SUNTF[sun, c23, c21]) - ColorDelta[c12, c21]*ColorDelta[c24, c11]*Contract[SUNTF[sun, c11, c12]*atr[xprop[-x] . v11 . xprop[x] . fcc[hv22 . GAD[lor1] . GAD[lorz] . GAD[lor2] . hv21]]*xprop2[x, lorb, lorz, lor1, lor2]]*FlavorDelta[f12, f21]*FlavorDelta[f23, f22]*FlavorDelta[f24, f11]*SUNTF[sun, c23, c22] - ColorDelta[c12, c22]*ColorDelta[c23, c11]*Contract[SUNTF[sun, c11, c12]*atr[fcc[hv22 . xprop[-x] . v11 . xprop[x] . hv21] . GAD[lor1] . GAD[lorz] . GAD[lor2]]*xprop2[x, lorb, lorz, lor1, lor2]]*FlavorDelta[f12, f22]*FlavorDelta[f23, f11]*FlavorDelta[f24, f21]*SUNTF[sun, c24, c21] - ColorDelta[c12, c21]*ColorDelta[c23, c11]*Contract[SUNTF[sun, c11, c12]*atr[GAD[lor1] . GAD[lorz] . GAD[lor2] . hv21 . fcc[hv22 . xprop[-x] . v11 . xprop[x]]]*xprop2[x, lorb, lorz, lor1, lor2]]*FlavorDelta[f12, f21]*FlavorDelta[f23, f11]*FlavorDelta[f24, f22]*SUNTF[sun, c24, c22];


dia1=FourierXP[dia1,{x,q}]//SUNSimplify;

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1 pole,q,HoldFlavor->holdf,EpsOrder->order,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
