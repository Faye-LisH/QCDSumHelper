(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)



(* ::Code::Initialization::Plain:: *)
Tetraquarkggg::usage="Tetraquarkggg[q_,j1_,j2_] gives the \[LeftAngleBracket]GGG\[RightAngleBracket] contribution of \[LeftAngleBracket]j1 j2\[RightAngleBracket], with j1, j2 are tetraquark currents. "

Tetraquarkggg::inderr="Dummy indices conflict!"

Begin["`Private`Tetraquarkggg`"]



(* ::Code::Initialization::Plain:: *)
Options[Tetraquarkggg]={
	Parallelized->True,
	AutoNDR->True,
	HoldFlavor->False,
	TypeTags->False,
	ToD4->"Auto"
}



(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
fcc[expr_]:=FCChargeConjugateTransposed[expr,Explicit->True]


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
Tetraquarkggg[qq_,factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 Tetraquarkggg[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])



Tetraquarkggg[qq_,current1_FourquarkCurrent,current2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,
v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,holdf=OptionValue[HoldFlavor],setd=D,atr,diagrams,fdir,files,pall},



(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=current1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=current2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[Tetraquarkgg::curerr];
	Abort[]
];



(*--- check dummy indices in input ---*)
tmp2=Length[#]&/@Gather[({c11,c12,c13,c14,c21,c22,c23,c24}//FCI)/.SUNFIndex[aa_]:>aa];

If[DummyCheck[v11,v12,v21,v22]==0||!FreeQ[tmp2,nn_Integer/;nn>2],
	Message[Tetraquarkggg::inderr];
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



If[OptionValue[AutoNDR]===True,
	atr=TR5
	,
	atr=TR
];

(*---------------------------------------------------*)

(* since the highest \[Epsilon]-pole is 1/\[Epsilon], it's unnecessary to evaluate trace at D-dimension *)
If[OptionValue[ToD4]===True,
	setd=4
,
	If[OptionValue[ToD4]==="Auto",
		diagrams={typeA1,typeA2,typeBS1,typeBS2,typeBS3,typeBS4,typeBS5,typeBS6,typeBS7,typeBS8,typeBS9,typeBS10,typeBS11,typeBS12,typeBS13,typeBS14,typeBS15,typeBS16,
				fdtypeC1,fdtypeC2,fdtypeC3,fdtypeC4,fdtypeC5,fdtypeC6,fdtypeC7,fdtypeC8,fdtypeC9,fdtypeC10,fdtypeC11,fdtypeC12,fdtypeC13,fdtypeC14,fdtypeC15,fdtypeC16,typeDS1,typeDS2,typeDS3,typeDS4,typeDS5,typeDS6,typeDS7,typeDS8}
	,
		diagrams={typeA1,typeA2,typeBS1,typeBS2,typeBS3,typeBS4,typeBS5,typeBS6,typeBS7,typeBS8,typeBS9,typeBS10,typeBS11,typeBS12,typeBS13,typeBS14,typeBS15,typeBS16,
				typeC1,typeC2,typeC3,typeC4,typeC5,typeC6,typeC7,typeC8,typeC9,typeC10,typeC11,typeC12,typeC13,typeC14,typeC15,typeC16,typeDS1,typeDS2,typeDS3,typeDS4,typeDS5,typeDS6,typeDS7,typeDS8}
	]
];


(*---------------------------------------------------*)
If[pall===True,

	
	tmp=Plus@@WaitAll[ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd]]&/@diagrams];
	
	tmp=QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	
,

	If[pall==="External",

	DistributeDefinitions[qq,current1,current2];
	
		If[fdir==="None",
		
			(* evaluation, no import and export *)
			 tmp=ParallelSubmit[{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd},
									#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd]]&/@diagrams
		
		,
		(* evaluation, import and export *)
			files=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files=("Tetraquarkggg_"<>#)&/@files;
			
			
			tmp=ImExport[fdir,
						files,
						{{v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd},
						{qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd},
						diagrams}
					]
				
		](* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied after the WaitAll[tmp]. !!! *)
	
	,
	
		
		tmp=Plus@@(#[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,setd]&/@diagrams);
		
		tmp=QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA;
		

	]
];


If[OptionValue[TypeTags]===True,
	tmp
,
	tmp/.Condensate[typ_/;StringTake[typ,3]=="Typ"]->1
]

]

(*------------------------------------------------------------------*)



(* ::Code::Initialization::Plain:: *)
(*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*)
(*-------------------------------------------------*)
(* propagators and condensates *)
(* see https://github.com/QSSRHelper/fourquark_generation/blob/main/GGG_condensate.nb for detials *)
(*-------------------------------------------------*)
(*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*)


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------*)
(*------------------------------------------------------------*)
(* gluon condensate *)


(* ::Code::Initialization::Plain:: *)
(* \[LeftAngleBracket]G_\[Mu]\[Nu] G_\[Alpha]\[Beta] G_\[Rho]\[Sigma]\[RightAngleBracket] ~ (((g_\[Nu]\[Alpha] g_\[Beta]\[Rho] g_\[Mu]\[Sigma] - {\[Mu]\[TwoWayRule]\[Nu]}) - {\[Alpha]\[TwoWayRule]\[Beta] }) - {\[Rho]\[TwoWayRule]\[Sigma]})\[LeftAngleBracket]GGG\[RightAngleBracket] *)

(* g_\[Nu]\[Alpha] g_\[Beta]\[Rho] g_\[Mu]\[Sigma] *)
gg3[lora_,lorb_,lorc_,lord_,lore_,lorf_]=MTD[lorb,lorc]MTD[lord,lore]MTD[lorf,lora];

(* (((g_\[Nu]\[Alpha] g_\[Beta]\[Rho] g_\[Mu]\[Sigma] - {\[Mu]\[TwoWayRule]\[Nu]}) - {\[Alpha]\[TwoWayRule]\[Beta] }) - {\[Rho]\[TwoWayRule]\[Sigma]}) *)
ga0[lora_,lorb_,lorc_,lord_,lore_,lorf_]=(-(MTD[lora, lorf]*MTD[lorb, lord]*MTD[lorc, lore]) + MTD[lora, lord]*MTD[lorb, lorf]*MTD[lorc, lore] + 
 MTD[lora, lore]*MTD[lorb, lord]*MTD[lorc, lorf] - MTD[lora, lord]*MTD[lorb, lore]*MTD[lorc, lorf] + 
 MTD[lora, lorf]*MTD[lorb, lorc]*MTD[lord, lore] - MTD[lora, lorc]*MTD[lorb, lorf]*MTD[lord, lore] - 
 MTD[lora, lore]*MTD[lorb, lorc]*MTD[lord, lorf] + MTD[lora, lorc]*MTD[lorb, lore]*MTD[lord, lorf])


(* ::Code::Initialization::Plain:: *)
(* \[LeftAngleBracket]G_\[Mu]\[Nu] D_\[Alpha] D_\[Beta] G_\[Rho]\[Sigma]\[RightAngleBracket] ~( A(((g_\[Nu]\[Alpha] g_\[Beta]\[Rho] g_\[Mu]\[Sigma] - {\[Mu]\[TwoWayRule]\[Nu]}) - {\[Alpha]\[TwoWayRule]\[Beta] }) - {\[Rho]\[TwoWayRule]\[Sigma]}) + B(((g_\[Nu]\[Alpha] g_\[Beta]\[Rho] g_\[Mu]\[Sigma] - {\[Mu]\[TwoWayRule]\[Nu]}) + {\[Alpha]\[TwoWayRule]\[Beta] }) - {\[Rho]\[TwoWayRule]\[Sigma]}) + C g_\[Alpha]\[Beta](g_\[Mu]\[Rho] g_\[Nu]\[Sigma] - g_\[Mu]\[Sigma] g_\[Nu]\[Rho]) )\[LeftAngleBracket]GGG\[RightAngleBracket] *)

(* (((g_\[Nu]\[Alpha] g_\[Beta]\[Rho] g_\[Sigma]\[Mu] - {\[Mu]\[TwoWayRule]\[Nu]}) + {\[Alpha]\[TwoWayRule]\[Beta] }) - {\[Rho]\[TwoWayRule]\[Sigma]}) *)
gs0[lora_,lorb_,lorc_,lord_,lore_,lorf_]=(MTD[lora, lorf]*MTD[lorb, lord]*MTD[lorc, lore] - MTD[lora, lord]*MTD[lorb, lorf]*MTD[lorc, lore] - 
 MTD[lora, lore]*MTD[lorb, lord]*MTD[lorc, lorf] + MTD[lora, lord]*MTD[lorb, lore]*MTD[lorc, lorf] + 
 MTD[lora, lorf]*MTD[lorb, lorc]*MTD[lord, lore] - MTD[lora, lorc]*MTD[lorb, lorf]*MTD[lord, lore] - 
 MTD[lora, lore]*MTD[lorb, lorc]*MTD[lord, lorf] + MTD[lora, lorc]*MTD[lorb, lore]*MTD[lord, lorf])


(* ::Code::Initialization::Plain:: *)
(* g_\[Alpha]\[Beta](g_\[Mu]\[Rho] g_\[Nu]\[Sigma] - g_\[Mu]\[Sigma] g_\[Nu]\[Rho] *)
g0[lora_,lorb_,lorc_,lord_,lore_,lorf_]=MTD[lorc,lord](MTD[lora,lore]MTD[lorb,lorf]-MTD[lora,lorf]MTD[lorb,lore]);


(* ::Code::Initialization::Plain:: *)
(*---------------------------------------------------*)
ggg[lora_,lorb_,lorc_,lord_,lore_,lorf_]=1/(2CA^2CF(D-2)(D-1)D)ga0[lora,lorb,lorc,lord,lore,lorf]

gddg[lora_,lorb_,lorc_,lord_,lore_,lorf_]=1/(4CA CF(D^2-4)(D-1)D)((D+2)ga0[lora,lorb,lorc,lord,lore,lorf]-(D-4)gs0[lora,lorb,lorc,lord,lore,lorf]+4(D-1)g0[lora,lorb,lorc,lord,lore,lorf])


(* ::Code::Initialization::Plain:: *)
(* \[LeftAngleBracket](D_\[Alpha] G_\[Mu]\[Nu])^a (D_\[Beta] G_\[Rho]\[Sigma])^b\[RightAngleBracket]\[Delta]^ab = -\[LeftAngleBracket] G^a_\[Mu]\[Nu] (D_\[Alpha] D_\[Beta] G_\[Rho]\[Sigma])^b\[RightAngleBracket]^ab *)
dgdg[lora_,lorb_,lorc_,lord_,lore_,lorf_]=-gddg[lorb,lorc,lora,lord,lore,lorf]


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------*)
(* basic propagators *)

(*xprop[x_]=FourierPX[prop[q],{q,x}];*)
xprop[x_]=1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];


(* ::Code::Initialization::Plain:: *)
(* quark propagator wiht background gluon *)
(*-----------------------------------------------------------------------------------------------------------*)

(* x \[Rule] -x if translate x \[Rule] -x *)
(*prog[q_,lora_,lorb_]=-1/2FAD[{q}]GSD[q].GAD[lorb].FourDivergence[GSD[q]FAD[{q}],FVD[q,lora]];*)


(* ::Code::Initialization::Plain:: *)
(*(*xprog[x_,lora_,lorb_]=FourierPX[prog[q,lora,lorb,cola],{q,x}]//DiracSimplify//QSimplify;*)
xprog[x_,lora_,lorb_]=(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*
  SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2));
  *)


(* ::Code::Initialization::Plain:: *)
(* xprodg[x_,lora_,lorb_,lorc_,cola_]=FourierPX[prodg[q,lora,lorb,lorc,cola],{q,x}]; *)
xprodg[x_,lora_,lorb_,lorc_] = (DiracGamma[LorentzIndex[lorc, D], D]*Pair[LorentzIndex[lora, D], LorentzIndex[lorb, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[-1/144*1/(E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
 DiracGamma[LorentzIndex[lora, D], D] . DiracGamma[LorentzIndex[lorc, D], D] . DiracGamma[LorentzIndex[lorb, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[1/(288*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
 DiracGamma[LorentzIndex[lorb, D], D] . DiracGamma[LorentzIndex[lorc, D], D] . DiracGamma[LorentzIndex[lora, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[1/(288*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + DiracGamma[LorentzIndex[lorb, D], D] . DiracGamma[LorentzIndex[lorc, D], D] . DiracGamma[Momentum[x, D], D]*
  Pair[LorentzIndex[lora, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qGamma[-1 + D/2]/(72*I^D*Pi^(D/2))] + 
 DiracGamma[LorentzIndex[lora, D], D] . DiracGamma[LorentzIndex[lorc, D], D] . DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lorb, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qGamma[-1 + D/2]/(72*I^D*Pi^(D/2))] + DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorc, D], D] . DiracGamma[LorentzIndex[lorb, D], D]*
  Pair[LorentzIndex[lora, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[-1/144*1/(E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]] + 
 DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorc, D], D] . DiracGamma[LorentzIndex[lora, D], D]*Pair[LorentzIndex[lorb, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qfact2[-1/144*1/(E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]] + 2*DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lora, D], LorentzIndex[lorb, D]]*
  Pair[LorentzIndex[lorc, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[-1/144*1/(E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]] + 
 (2*DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lora, D], Momentum[x, D]]*Pair[LorentzIndex[lorb, D], Momentum[x, D]]*Pair[LorentzIndex[lorc, D], Momentum[x, D]]*
   qfact1[-1/36*qGamma[D/2]/(I^D*Pi^(D/2))])/Pair[Momentum[x, D], Momentum[x, D]]^(D/2) + DiracGamma[LorentzIndex[lorc, D], D]*Pair[LorentzIndex[lora, D], Momentum[x, D]]*
  Pair[LorentzIndex[lorb, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[-1/72*(-2 + D)/(E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2] + qGamma[D/2]/(36*I^D*Pi^(D/2))]*
  SUNT[SUNIndex[cola]]);
  


(* ::Code::Initialization::Plain:: *)
(*(* prodg[q_,lora_,lorb_,lorc_,cola_]=I/3FAD[{q}]GSD[q].GAD[lorc].SUNT[cola].FourDivergence[FourDivergence[GSD[q]FAD[{q}],FVD[q,lora]],FVD[q,lorb]]; *)*)


(* ::Code::Initialization::Plain:: *)
(*(* proddg[q_,lora_,lorb_,lorc_,lord_,cola_]=1/8FAD[{q}]GSD[q].GAD[lord].SUNT[cola].FourDivergence[FourDivergence[FourDivergence[GSD[q]FAD[{q}],FVD[q,lora]],FVD[q,lorb]],FVD[q,lorc]]; *)
(* xproddg[x_,lora_,lorb_,lorc_,lord_]=FourierPX[proddg[q,lora,lorb,lorc,lord],{q,x}]; *)*)


(* ::Code::Initialization::Plain:: *)
(*(*progg[q_,lora_,lorb_,lorc_,lord_,cola_,colb_]=-I/4FAD[{q}]GSD[q].GAD[lorb].SUNT[cola].SUNT[colb].FourDivergence[FAD[{q}]GSD[q].GAD[lord].FourDivergence[GSD[q]FAD[{q}],FVD[q,lorc]],FVD[q,lora]];*)
(*xprogg[x_,lora_,lorb_,lorc_,lord_]=FourierPX[progg[q,lora,lorb,lorc,lord],{q,x}];*)*)


(* ::Code::Initialization::Plain:: *)
(*(* xprog2 SUNTF[colc] = xprogg[xx,lora,lorb,lorc,lord,cola,colb]ggg[lora,lorb,lorc,lord,lore,lorf,cola,colb,colc]+xproddg[xx,lora,lorb,lorc,lord,cola]gddg[lore,lorf,lora,lorb,lorc,lord,cola,colc] *)

xprog2[x_,lore_,lorf_]=*)


(* ::Code::Initialization::Plain:: *)
(*xprog[x_,lora_,lorb_]=FourierPX[prog[q,lora,lorb,cola],{q,x}]//DiracSimplify//QSimplify;*)
xprog[x_,lora_,lorb_]=(GAD[lora] . GAD[lorb] . GSD[x]*qfact1[qGamma[-1 + D/2]/(32*I^D*Pi^(D/2))]*SPD[x, x]^(1 - D/2) + GSD[x] . GAD[lorb] . GAD[lora]*qfact1[qfact2[-1/(32*I^D*Pi^(D/2))]*qGamma[-1 + D/2]]*
  SPD[x, x]^(1 - D/2) + (2*FVD[x, lora]*FVD[x, lorb]*GSD[x]*qfact1[-qGamma[D/2]/(16*I^D*Pi^(D/2))])/SPD[x, x]^(D/2));
  


(* ::Code::Initialization::Plain:: *)
(* xprodg1 ~ Contract[ xprodg <DGDG> ] *)
xprodg1[x_,lord_,lore_,lorf_]=(DiracGamma[LorentzIndex[lorf, D], D]*Pair[LorentzIndex[lord, D], LorentzIndex[lore, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[-((CA - 2*CF)*(-16 + 7*D))/(192*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
 DiracGamma[LorentzIndex[lore, D], D]*Pair[LorentzIndex[lord, D], LorentzIndex[lorf, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[((CA - 2*CF)*(-16 + 7*D))/(192*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
 DiracGamma[LorentzIndex[lord, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lorf, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(CA - 2*CF)/(192*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
 DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lord, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(CA - 2*CF)/(192*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
 DiracGamma[LorentzIndex[lord, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lore, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(-CA + 2*CF)/(192*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
 DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lord, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(-CA + 2*CF)/(192*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lord, D], LorentzIndex[lore, D]]*
  Pair[LorentzIndex[lorf, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[-((CA - 2*CF)*(-4 + D))/(96*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]
 + DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lord, D], LorentzIndex[lorf, D]]*Pair[LorentzIndex[lore, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qfact2[((CA - 2*CF)*(-4 + D))/(96*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]] + 
 DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lore, D], D]*Pair[LorentzIndex[lord, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qfact2[(CA - 2*CF)/(96*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]] + DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lord, D], D]*
  Pair[LorentzIndex[lore, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[(CA - 2*CF)/(96*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]] + 
 DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lord, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qfact2[(CA - 2*CF)/(48*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]] + DiracGamma[LorentzIndex[lord, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[Momentum[x, D], D]*
  Pair[LorentzIndex[lorf, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[(CA - 2*CF)/(48*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]] + 
 DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lorf, D], D]*Pair[LorentzIndex[lord, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qfact2[(-CA + 2*CF)/(96*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]] + DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lord, D], D]*
  Pair[LorentzIndex[lorf, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[(-CA + 2*CF)/(96*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]] + 
 DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lord, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qfact2[(-CA + 2*CF)/(48*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]] + DiracGamma[LorentzIndex[lord, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[Momentum[x, D], D]*
  Pair[LorentzIndex[lore, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[(-CA + 2*CF)/(48*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]])



(* ::Code::Initialization::Plain:: *)
(*  xprog2[x_,lore_,lorf_,colc_] = xprogg[x_,lora_,lorb_,lorc_,lord_,cola_,colb_] ggg[lora_,lorb_,lorc_,lord_,lore_,lorf_,cola_,colb_,colc_] + xproddg[x_,lora_,lorb_,lorc_,lord_,cola_] gddg[lore_,lorf_,lora_,lorb_,lorc_,lord_,cola_,colc_]  *)
xprog2[x_,lore_,lorf_]=(DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lore, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(9 - 2*D)/(6144*CA*CF*(-2 + D)*(-1 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
 DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorf, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(-9 + 2*D)/(6144*CA*CF*(-2 + D)*(-1 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
 DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lore, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(-152 + 134*D + D^2 - 4*D^3)/(12288*CA*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
 DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lorf, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(152 - 134*D - D^2 + 4*D^3)/(12288*CA*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
 DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[Momentum[x, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(280 - 202*D + 27*D^2 + 3*D^3)/(12288*CA*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2] + 
    qfact2[-1/1024*1/(CA*CF*(-2 + D)*(-1 + D)*D*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[D/2]] + DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[Momentum[x, D], D]*
  Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*qfact1[qfact2[-1/12288*(280 - 202*D + 27*D^2 + 3*D^3)/(CA*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2] + 
    qfact2[1/(1024*CA*CF*(-2 + D)*(-1 + D)*D*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[D/2]] + DiracGamma[LorentzIndex[lore, D], D]*Pair[LorentzIndex[lorf, D], Momentum[x, D]]*
  Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*qfact1[qfact2[(-32 + 140*D - 124*D^2 + 19*D^3)/(6144*CA*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2] + 
    qfact2[(14 - 5*D)/(512*CA*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[D/2]] + DiracGamma[LorentzIndex[lorf, D], D]*Pair[LorentzIndex[lore, D], Momentum[x, D]]*
  Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*qfact1[qfact2[(32 - 140*D + 124*D^2 - 19*D^3)/(6144*CA*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2] + 
    qfact2[(-14 + 5*D)/(512*CA*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[D/2]])
    


(* ::Code::Initialization::Plain:: *)
(* prog3 = proggg \[LeftAngleBracket]GGG\[RightAngleBracket] + progddg \[LeftAngleBracket]GDDG\[RightAngleBracket] + proddgg \[LeftAngleBracket]DDGG\[RightAngleBracket] + prodgdg \[LeftAngleBracket]DGDG\[RightAngleBracket] *)
xprog3[x_]= GSD[x] qfact1[qfact2[(I I^(-D) (CA-2 CF) (24-10 D+D^2) \[Pi]^(-D/2))/(9216(2+D)D)] qGamma[-3+D/2]] SPD[x,x]^(3-D/2);



(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)
(* ------------- The following are generated by algorithem -----------------*)
(*----------------------------------------------------------------------*)


(* ::Input::Initialization:: *)
typeA1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia1,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia1=Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog3[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog3[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog3[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[xprog3[x], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog3[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog3[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[xprog3[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog3[-x], v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia1=FourierXP[dia1,{x,q}];

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeA2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia2,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia2=Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog3[x], hv21]], v11, xprop[x]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog3[-x], v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprog3[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprog3[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c14]*ColorDelta[c24, c13]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[cc[dot[xprog3[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog3[-x], v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, xprog3[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]] + Condensate["TypeA"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog3[-x]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*ColorDelta[c23, c13]*ColorDelta[c24, c14]];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia2=FourierXP[dia2,{x,q}];

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)


(* ::Input::Initialization:: *)
typeBS1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia1,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia1=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[xprog2[x, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog[-x, lore, lorf]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprog[x, lore, lorf], hv21]]]]*tr[dot[cc[dot[xprog2[-x, lore, lorf], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog2[-x, lore, lorf]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14]];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia1=FourierXP[dia1,{x,q}];

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia2,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia2=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog2[x, lore, lorf], hv21]], v11, xprog[x, lore, lorf]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[cc[dot[xprog2[x, lore, lorf], hv21]], v11, xprog[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprop[x]]]*tr[dot[cc[dot[xprog2[-x, lore, lorf], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13]];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia2=FourierXP[dia2,{x,q}];

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia3,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia3=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprog2[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprog2[x, lore, lorf], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprog2[x, lore, lorf], hv21]]]]*tr[dot[cc[dot[xprog[-x, lore, lorf], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13]];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia3=FourierXP[dia3,{x,q}];

dia3=QEvaluate[I ScaleMu^(3(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia4,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia4=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog2[-x, lore, lorf], v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprog2[x, lore, lorf], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog[-x, lore, lorf]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprog[-x, lore, lorf], v12]], hv22, xprog2[-x, lore, lorf]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia4=FourierXP[dia4,{x,q}];

dia4=QEvaluate[I ScaleMu^(3(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia5,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia5=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog[-x, lore, lorf], v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, xprog2[x, lore, lorf], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprop[x]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog2[-x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog[-x, lore, lorf]]]*tr[dot[xprop[x], cc[dot[v11, xprog2[x, lore, lorf], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14]];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia5=FourierXP[dia5,{x,q}];

dia5=QEvaluate[I ScaleMu^(3(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia6,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia6=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog[-x, lore, lorf], v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprog2[x, lore, lorf], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog[-x, lore, lorf]]]*tr[dot[xprog2[x, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog[-x, lore, lorf]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog2[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14]];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia6=FourierXP[dia6,{x,q}];

dia6=QEvaluate[I ScaleMu^(3(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia7,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia7=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprog2[x, lore, lorf], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog[-x, lore, lorf], v12]], xprop[-x]]]*tr[dot[cc[dot[xprog2[x, lore, lorf], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog2[-x, lore, lorf]]]*tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14]];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia7=FourierXP[dia7,{x,q}];

dia7=QEvaluate[I ScaleMu^(3(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia8,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia8=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog2[-x, lore, lorf]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog2[-x, lore, lorf], v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog[-x, lore, lorf], v12]], xprog2[-x, lore, lorf]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia8=FourierXP[dia8,{x,q}];

dia8=QEvaluate[I ScaleMu^(3(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS9[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia9,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia9=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[xprog2[x, lore, lorf], cc[dot[v11, xprog[x, lore, lorf], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[xprog2[x, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprog[-x, lore, lorf], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog2[-x, lore, lorf], v12]], xprog[-x, lore, lorf]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia9=FourierXP[dia9,{x,q}];

dia9=QEvaluate[I ScaleMu^(3(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS10[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia10,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia10=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[xprog2[x, lore, lorf], cc[dot[v11, xprog[x, lore, lorf], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog[-x, lore, lorf]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog2[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprog2[-x, lore, lorf], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13]];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia10=FourierXP[dia10,{x,q}];

dia10=QEvaluate[I ScaleMu^(3(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS11[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia11,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia11=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog[-x, lore, lorf], v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog2[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog2[-x, lore, lorf], v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, xprog[x, lore, lorf], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprog2[-x, lore, lorf], v12]], hv22, xprog[-x, lore, lorf]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia11=FourierXP[dia11,{x,q}];

dia11=QEvaluate[I ScaleMu^(3(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS12[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia12,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia12=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog2[-x, lore, lorf], v12]], xprop[-x]]]*tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog2[x, lore, lorf], hv21]], v11, xprop[x]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog[-x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog[-x, lore, lorf], v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog2[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13]];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia12=FourierXP[dia12,{x,q}];

dia12=QEvaluate[I ScaleMu^(3(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS13[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia13,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia13=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprog2[x, lore, lorf]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog2[-x, lore, lorf]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog[-x, lore, lorf], v12]], xprog2[-x, lore, lorf]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]];


dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia13=FourierXP[dia13,{x,q}];

dia13=QEvaluate[I ScaleMu^(3(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS14[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia14,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia14=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog2[-x, lore, lorf], v12]], xprop[-x]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprog[x, lore, lorf], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog2[-x, lore, lorf]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog2[-x, lore, lorf]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14]];


dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia14=FourierXP[dia14,{x,q}];

dia14=QEvaluate[I ScaleMu^(3(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS15[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia15,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia15=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog[-x, lore, lorf], v12]], hv22, xprog2[-x, lore, lorf]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog2[-x, lore, lorf], v12]], hv22, xprog[-x, lore, lorf]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog[-x, lore, lorf]]]*tr[dot[cc[dot[xprog2[x, lore, lorf], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14]];


dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia15=FourierXP[dia15,{x,q}];

dia15=QEvaluate[I ScaleMu^(3(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeBS16[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia16,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia16=Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog[-x, lore, lorf], v12]], xprop[-x]]]*tr[dot[xprog2[x, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog2[-x, lore, lorf]]]*tr[dot[xprop[x], cc[dot[v11, xprog[x, lore, lorf], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14]] + Condensate["TypeB"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog2[-x, lore, lorf], v12]], xprog[-x, lore, lorf]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]];


dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia16=FourierXP[dia16,{x,q}];

dia16=QEvaluate[I ScaleMu^(3(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)


(* ::Code::Initialization::Plain:: *)
(* set the dimension of gamma-matrices as 4 *)
fdtypeC1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC1[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC2[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC3[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC4[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC5[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC6[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC7[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC8[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC9[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC9[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC10[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC10[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC11[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC11[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC12[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC12[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC13[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC13[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC14[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC14[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC15[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC15[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]
fdtypeC16[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=typeC16[qq,v11,v12,hv21,hv22,c11,c12,c13,c14,f11,f12,f13,f14,c21,c22,c23,c24,f21,f22,f23,f24,holdf,atr,4]


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)


(* ::Input::Initialization:: *)
typeC1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia1,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia1=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprog[x, lorc, lord], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog[-x, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c24, c13]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c14]*SUNTF[sun[2], c12, c22]*SUNTF[sun[3], c11, c21]];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia1=FourierXP[dia1,{x,q}];

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia2,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia2=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprog[x, lorc, lord]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog[-x, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c24, c13]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c14]*SUNTF[sun[2], c12, c21]*SUNTF[sun[3], c11, c22]];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia2=FourierXP[dia2,{x,q}];

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia3,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia3=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[cc[dot[hv22, xprog[-x, lorc, lord], v12]], xprog[-x, lora, lorb]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c14]*SUNTF[sun[2], c23, c13]*SUNTF[sun[3], c12, c21]];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia3=FourierXP[dia3,{x,q}];

dia3=QEvaluate[I ScaleMu^(3(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia4,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia4=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[cc[dot[hv22, xprog[-x, lorc, lord], v12]], xprog[-x, lora, lorb]]]*tr[dot[xprop[x], cc[dot[v11, xprog[x, lore, lorf], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c14]*SUNTF[sun[2], c23, c13]*SUNTF[sun[3], c12, c22]];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia4=FourierXP[dia4,{x,q}];

dia4=QEvaluate[I ScaleMu^(3(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia5,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia5=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[cc[dot[hv22, xprog[-x, lora, lorb], v12]], xprop[-x]]]*tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprog[x, lorc, lord], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c24, c14]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c13]*SUNTF[sun[2], c12, c22]*SUNTF[sun[3], c11, c21]];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia5=FourierXP[dia5,{x,q}];

dia5=QEvaluate[I ScaleMu^(3(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia6,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia6=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprog[x, lorc, lord], hv21]]]]*tr[dot[cc[dot[xprog[-x, lora, lorb], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c13]*SUNTF[sun[2], c12, c22]*SUNTF[sun[3], c11, c21]];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia6=FourierXP[dia6,{x,q}];

dia6=QEvaluate[I ScaleMu^(3(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia7,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia7=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog[-x, lora, lorb]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprog[x, lorc, lord]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c14]*SUNTF[sun[2], c12, c21]*SUNTF[sun[3], c11, c22]];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia7=FourierXP[dia7,{x,q}];

dia7=QEvaluate[I ScaleMu^(3(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia8,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia8=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[cc[dot[hv22, xprog[-x, lorc, lord], v12]], xprog[-x, lora, lorb]]]*tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c14]*SUNTF[sun[2], c23, c13]*SUNTF[sun[3], c11, c21]];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia8=FourierXP[dia8,{x,q}];

dia8=QEvaluate[I ScaleMu^(3(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC9[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia9,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia9=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[cc[dot[xprog[-x, lorc, lord], v12]], hv22, xprog[-x, lora, lorb]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c14]*SUNTF[sun[2], c24, c13]*SUNTF[sun[3], c12, c21]];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia9=FourierXP[dia9,{x,q}];

dia9=QEvaluate[I ScaleMu^(3(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC10[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia10,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia10=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[xprop[x], cc[dot[v11, xprog[x, lore, lorf], hv21]]]]*tr[dot[cc[dot[xprog[-x, lorc, lord], v12]], hv22, xprog[-x, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c14]*SUNTF[sun[2], c24, c13]*SUNTF[sun[3], c12, c22]];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia10=FourierXP[dia10,{x,q}];

dia10=QEvaluate[I ScaleMu^(3(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC11[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia11,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia11=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[cc[dot[xprog[-x, lorc, lord], v12]], hv22, xprog[-x, lora, lorb]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c14]*SUNTF[sun[2], c24, c13]*SUNTF[sun[3], c11, c22]];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia11=FourierXP[dia11,{x,q}];

dia11=QEvaluate[I ScaleMu^(3(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC12[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia12,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia12=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[cc[dot[xprog[-x, lora, lorb], v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprog[x, lorc, lord]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c13]*SUNTF[sun[2], c12, c21]*SUNTF[sun[3], c11, c22]];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia12=FourierXP[dia12,{x,q}];

dia12=QEvaluate[I ScaleMu^(3(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC13[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia13,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia13=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprog[-x, lorc, lord], v12]], hv22, xprog[-x, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c14]*SUNTF[sun[2], c24, c13]*SUNTF[sun[3], c11, c21]];


dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia13=FourierXP[dia13,{x,q}];

dia13=QEvaluate[I ScaleMu^(3(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC14[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia14,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia14=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog[-x, lora, lorb]]]*tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprog[x, lorc, lord], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c14]*SUNTF[sun[2], c12, c22]*SUNTF[sun[3], c11, c21]];


dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia14=FourierXP[dia14,{x,q}];

dia14=QEvaluate[I ScaleMu^(3(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC15[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia15,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia15=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[cc[dot[hv22, xprog[-x, lora, lorb], v12]], xprop[-x]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprog[x, lorc, lord]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c24, c14]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c13]*SUNTF[sun[2], c12, c21]*SUNTF[sun[3], c11, c22]];


dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia15=FourierXP[dia15,{x,q}];

dia15=QEvaluate[I ScaleMu^(3(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeC16[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia16,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia16=Condensate["TypeC"]Condensate["ggg"]*contract[ggg[lora, lorb, lorc, lord, lore, lorf]*tr[dot[cc[dot[hv22, xprog[-x, lorc, lord], v12]], xprog[-x, lora, lorb]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c14]*SUNTF[sun[2], c23, c13]*SUNTF[sun[3], c11, c22]];


dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia16=FourierXP[dia16,{x,q}];

dia16=QEvaluate[I ScaleMu^(3(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
(*typeC1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia1,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia1=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[x, lore, lorf], cc[dot[v11, xprog[x, lorc, lord], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog[-x, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c24, c13]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c14]*SUNTF[sun[2], c12, c22]*SUNTF[sun[3], c11, c21]];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia1=FourierXP[dia1,{x,q}];

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia2,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia2=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[x, lorc, lord]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprog[-x, lora, lorb]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c24, c13]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c14]*SUNTF[sun[2], c12, c21]*SUNTF[sun[3], c11, c22]];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia2=FourierXP[dia2,{x,q}];

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia3,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia3=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog[-x, lorc, lord], v12]], ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[-x, lora, lorb]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c14]*SUNTF[sun[2], c23, c13]*SUNTF[sun[3], c12, c21]];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia3=FourierXP[dia3,{x,q}];

dia3=QEvaluate[I ScaleMu^(3(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia4,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia4=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog[-x, lorc, lord], v12]], ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[-x, lora, lorb]]]*tr[dot[xprop[x], cc[dot[v11, xprog[x, lore, lorf], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c14]*SUNTF[sun[2], c23, c13]*SUNTF[sun[3], c12, c22]];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia4=FourierXP[dia4,{x,q}];

dia4=QEvaluate[I ScaleMu^(3(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia5,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia5=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog[-x, lora, lorb], v12]], xprop[-x]]]*tr[dot[ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[x, lore, lorf], cc[dot[v11, xprog[x, lorc, lord], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c24, c14]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c13]*SUNTF[sun[2], c12, c22]*SUNTF[sun[3], c11, c21]];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia5=FourierXP[dia5,{x,q}];

dia5=QEvaluate[I ScaleMu^(3(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia6,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia6=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[x, lore, lorf], cc[dot[v11, xprog[x, lorc, lord], hv21]]]]*tr[dot[cc[dot[xprog[-x, lora, lorb], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c13]*SUNTF[sun[2], c12, c22]*SUNTF[sun[3], c11, c21]];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia6=FourierXP[dia6,{x,q}];

dia6=QEvaluate[I ScaleMu^(3(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia7,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia7=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog[-x, lora, lorb]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[x, lorc, lord]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c14]*SUNTF[sun[2], c12, c21]*SUNTF[sun[3], c11, c22]];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia7=FourierXP[dia7,{x,q}];

dia7=QEvaluate[I ScaleMu^(3(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia8,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia8=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog[-x, lorc, lord], v12]], ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[-x, lora, lorb]]]*tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c14]*SUNTF[sun[2], c23, c13]*SUNTF[sun[3], c11, c21]];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia8=FourierXP[dia8,{x,q}];

dia8=QEvaluate[I ScaleMu^(3(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC9[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia9,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia9=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog[-x, lorc, lord], v12]], hv22, ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[-x, lora, lorb]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprog[x, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c14]*SUNTF[sun[2], c24, c13]*SUNTF[sun[3], c12, c21]];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia9=FourierXP[dia9,{x,q}];

dia9=QEvaluate[I ScaleMu^(3(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC10[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia10,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia10=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprog[x, lore, lorf], hv21]]]]*tr[dot[cc[dot[xprog[-x, lorc, lord], v12]], hv22, ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[-x, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c14]*SUNTF[sun[2], c24, c13]*SUNTF[sun[3], c12, c22]];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia10=FourierXP[dia10,{x,q}];

dia10=QEvaluate[I ScaleMu^(3(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC11[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia11,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia11=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog[-x, lorc, lord], v12]], hv22, ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[-x, lora, lorb]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c14]*SUNTF[sun[2], c24, c13]*SUNTF[sun[3], c11, c22]];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia11=FourierXP[dia11,{x,q}];

dia11=QEvaluate[I ScaleMu^(3(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC12[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia12,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia12=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprog[-x, lora, lorb], v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[x, lorc, lord]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c13]*SUNTF[sun[2], c12, c21]*SUNTF[sun[3], c11, c22]];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia12=FourierXP[dia12,{x,q}];

dia12=QEvaluate[I ScaleMu^(3(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC13[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia13,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia13=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[xprog[x, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprog[-x, lorc, lord], v12]], hv22, ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[-x, lora, lorb]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c14]*SUNTF[sun[2], c24, c13]*SUNTF[sun[3], c11, c21]];


dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia13=FourierXP[dia13,{x,q}];

dia13=QEvaluate[I ScaleMu^(3(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC14[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia14,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia14=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprog[-x, lora, lorb]]]*tr[dot[ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[x, lore, lorf], cc[dot[v11, xprog[x, lorc, lord], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c14]*SUNTF[sun[2], c12, c22]*SUNTF[sun[3], c11, c21]];


dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia14=FourierXP[dia14,{x,q}];

dia14=QEvaluate[I ScaleMu^(3(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC15[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia15,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia15=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog[-x, lora, lorb], v12]], xprop[-x]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[x, lorc, lord]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c24, c14]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c23, c13]*SUNTF[sun[2], c12, c21]*SUNTF[sun[3], c11, c22]];


dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia15=FourierXP[dia15,{x,q}];

dia15=QEvaluate[I ScaleMu^(3(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Input::Initialization:: *)
(*typeC16[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia16,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia16=Condensate["TypeC"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprog[-x, lorc, lord], v12]], ggg[lora, lorb, lorc, lord, lore, lorf]*xprog[-x, lora, lorb]]]*tr[dot[cc[dot[xprog[x, lore, lorf], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*SUNF[sun[1], sun[2], sun[3]]*SUNTF[sun[1], c24, c14]*SUNTF[sun[2], c23, c13]*SUNTF[sun[3], c11, c22]];


dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.{dot->Dot,cc->fcc}/.{tr->atr,contract->Contract}/.sf_sunSimplify:>SUNSimplify[FCI[sf/.sunSimplify->SUNSimplify]/.SUNTrace[sunts_,op___Rule]ff_SUNF:>(SUNTrace[sunts ff,##]&@@({op}/.(Explicit->False)->(Explicit->True)))];


dia16=FourierXP[dia16,{x,q}];

dia16=QEvaluate[I ScaleMu^(3(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]*)


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)


(* ::Input::Initialization:: *)
typeDS1[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia1,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia1=Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprodg1[-x, lord, lore, lorf], v12]], hv22, xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprodg[x, lord, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c13]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprodg[x, lord, lore, lorf], hv21]]]]*tr[dot[cc[dot[xprodg1[-x, lord, lore, lorf], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c13]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprodg1[-x, lord, lore, lorf]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprodg[x, lord, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c23, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c24, c14]];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia1=FourierXP[dia1,{x,q}];

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeDS2[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia2,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia2=Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[xprodg[x, lord, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprodg1[-x, lord, lore, lorf]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c14]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprodg[x, lord, lore, lorf], hv21]], v11, xprop[x]]]*tr[dot[cc[dot[xprodg1[-x, lord, lore, lorf], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c13]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprodg[-x, lord, lore, lorf], v12]], xprodg1[-x, lord, lore, lorf]]]*tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia2=FourierXP[dia2,{x,q}];

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeDS3[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia3,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia3=Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[xprodg[x, lord, lore, lorf], cc[dot[v11, xprodg1[x, lord, lore, lorf], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprop[-x], v12]], hv22, xprodg1[-x, lord, lore, lorf]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprodg[x, lord, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c13]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c14]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprodg[x, lord, lore, lorf], hv21]]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprodg1[-x, lord, lore, lorf]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c14]];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia3=FourierXP[dia3,{x,q}];

dia3=QEvaluate[I ScaleMu^(3(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeDS4[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia4,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia4=Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[xprodg[x, lord, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprodg1[-x, lord, lore, lorf], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c13]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprodg1[-x, lord, lore, lorf]]]*tr[dot[xprodg[x, lord, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c23, c13]*SUNTF[sun, c11, c21]*SUNTF[sun, c24, c14]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprodg[-x, lord, lore, lorf], v12]], xprodg1[-x, lord, lore, lorf]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c13]*SUNTF[sun, c24, c14]];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia4=FourierXP[dia4,{x,q}];

dia4=QEvaluate[I ScaleMu^(3(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeDS5[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia5,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia5=Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprodg1[-x, lord, lore, lorf], v12]], xprop[-x]]]*tr[dot[cc[dot[xprodg[x, lord, lore, lorf], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c13]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprodg[x, lord, lore, lorf], hv21]], v11, xprop[x]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprodg1[-x, lord, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c23, c14]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprodg1[-x, lord, lore, lorf]]]*tr[dot[xprop[x], cc[dot[v11, xprodg[x, lord, lore, lorf], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c23, c13]*SUNTF[sun, c12, c22]*SUNTF[sun, c24, c14]];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia5=FourierXP[dia5,{x,q}];

dia5=QEvaluate[I ScaleMu^(3(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeDS6[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia6,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia6=Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprodg1[-x, lord, lore, lorf], v12]], xprop[-x]]]*tr[dot[xprodg[x, lord, lore, lorf], cc[dot[v11, xprop[x], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c22]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c23, c13]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprodg1[-x, lord, lore, lorf], v12]], xprop[-x]]]*tr[dot[xprop[x], cc[dot[v11, xprodg[x, lord, lore, lorf], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c24, c14]*SUNTF[sun, c12, c22]*SUNTF[sun, c23, c13]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprodg1[-x, lord, lore, lorf]]]*tr[dot[cc[dot[xprodg[x, lord, lore, lorf], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c12, c21]*ColorDelta[c23, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c24, c14]];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia6=FourierXP[dia6,{x,q}];

dia6=QEvaluate[I ScaleMu^(3(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeDS7[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia7,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia7=Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprodg1[-x, lord, lore, lorf], v12]], xprop[-x]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprodg[x, lord, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c24, c14]*SUNTF[sun, c12, c21]*SUNTF[sun, c23, c13]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprodg[-x, lord, lore, lorf], v12]], hv22, xprodg1[-x, lord, lore, lorf]]]*tr[dot[cc[dot[xprop[x], hv21]], v11, xprop[x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c22]*ColorDelta[c12, c21]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[xprop[x], cc[dot[v11, xprop[x], hv21]]]]*tr[dot[cc[dot[xprodg[-x, lord, lore, lorf], v12]], hv22, xprodg1[-x, lord, lore, lorf]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c11, c21]*ColorDelta[c12, c22]*SUNTF[sun, c23, c14]*SUNTF[sun, c24, c13]];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia7=FourierXP[dia7,{x,q}];

dia7=QEvaluate[I ScaleMu^(3(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Input::Initialization:: *)
typeDS8[qq_,v11_,v12_,hv21_,hv22_,c11_,c12_,c13_,c14_,f11_,f12_,f13_,f14_,c21_,c22_,c23_,c24_,f21_,f22_,f23_,f24_,holdf_,atr_,setd_]:=Block[{x,q,dia8,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,cc,tr,dot,sunSimplify,contract},


dia8=Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[xprodg[x, lord, lore, lorf], hv21]], v11, xprodg1[x, lord, lore, lorf]]]*tr[dot[cc[dot[xprop[-x], v12]], hv22, xprop[-x]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f14]*FlavorDelta[f24, f13]*sunSimplify[ColorDelta[c23, c14]*ColorDelta[c24, c13]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[cc[dot[xprodg[x, lord, lore, lorf], hv21]], v11, xprodg1[x, lord, lore, lorf]]]]*FlavorDelta[f11, f22]*FlavorDelta[f12, f21]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c22]*SUNTF[sun, c12, c21]] + Condensate["TypeD"]Condensate["ggg"]*contract[tr[dot[cc[dot[hv22, xprop[-x], v12]], xprop[-x]]]*tr[dot[xprodg[x, lord, lore, lorf], cc[dot[v11, xprodg1[x, lord, lore, lorf], hv21]]]]]*FlavorDelta[f11, f21]*FlavorDelta[f12, f22]*FlavorDelta[f23, f13]*FlavorDelta[f24, f14]*sunSimplify[ColorDelta[c23, c13]*ColorDelta[c24, c14]*SUNTF[sun, c11, c21]*SUNTF[sun, c12, c22]];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{sunSimplify->SUNSimplify,cc->fcc}/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.contract->Contract;


dia8=FourierXP[dia8,{x,q}];

dia8=QEvaluate[I ScaleMu^(3(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA)
]


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
