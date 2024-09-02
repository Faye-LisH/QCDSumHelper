(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)



(* ::Code::Initialization::Plain:: *)
Fourquarkggg::usage="Fourquarkggg[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}}] give the \[LeftAngleBracket]GGG\[RightAngleBracket] contribution of OPE result about <J1 J2>, with J1=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a1\)]\)\!\(\*SubscriptBox[\(v1\[CapitalPsi]\), \(b1\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c1\)]\)\!\(\*SubscriptBox[\(v2\[CapitalPsi]\), \(d1\)]\) and J2=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(b2\)]\)\!\(\*SubscriptBox[\(v3\[CapitalPsi]\), \(a2\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(d2\)]\)\!\(\*SubscriptBox[\(v4\[CapitalPsi]\), \(c2\)]\) "

Fourquarkggg::inderr="Dummy indices conflict!"

Begin["`Private`Fourquarkggg`"]



(* ::Code::Initialization::Plain:: *)
Options[Fourquarkggg]={
	Parallelized->True,
	AutoNDR->True,
	HoldFlavor->False,
	TraceFirst->False,
	TypeTags->False,
	ToD4->"Auto"
}


(* ::Code::Initialization::Plain:: *)
(*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*)
(*-------------------------------------------------*)
(* propagators and condensates *)
(*-------------------------------------------------*)
(*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*)


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------*)
(*------------------------------------------------------------*)
(* gluon condensate *)


(* ::Code::Initialization::Plain:: *)
(* i\[LeftAngleBracket]tr[G_\[Mu]\[Nu] G_\[Nu]\[Rho] G_\[Rho]\[Mu]]\[RightAngleBracket] = -(1/4)\[LeftAngleBracket]g^3  f^def  G^d_\[Mu]\[Nu] G^e_\[Nu]\[Rho] G^f_\[Rho]\[Mu]\[RightAngleBracket] *)
(* g_\[Nu]\[Alpha] g_\[Beta]\[Rho] g_\[Mu]\[Sigma] *)
gg3[lora_,lorb_,lorc_,lord_,lore_,lorf_]=MTD[lorb,lorc]MTD[lord,lore]MTD[lorf,lora];

(* \[LeftAngleBracket]\!\(
\*SubsuperscriptBox[\(G\), \(\[Mu]\[Nu]\), \(a\)]
\*SubsuperscriptBox[\(G\), \(\[Alpha]\[Beta]\), \(b\)]
\*SubsuperscriptBox[\(G\), \(\[Rho]\[Sigma]\), \(c\)]\)\[RightAngleBracket]=\[LeftAngleBracket]f^def G^d_\[Mu]\[Nu] G^e_\[Nu]\[Rho] G^f_\[Rho]\[Mu]\[RightAngleBracket]*(f^abc/(2 (C_A)^2  C_F  D(D-1)(D-2)))(((g_\[Nu]\[Alpha] g_\[Beta]\[Rho] g_\[Mu]\[Sigma]-{\[Mu]\[TwoWayRule]\[Nu]})-{\[Alpha]\[TwoWayRule]\[Beta] })-{\[Rho]\[TwoWayRule]\[Sigma]}) *)

ga0[lora_,lorb_,lorc_,lord_,lore_,lorf_]=(((gg3[lora,lorb,lorc,lord,lore,lorf]-gg3[lorb,lora,lorc,lord,lore,lorf])-(gg3[lora,lorb,lord,lorc,lore,lorf]-gg3[lorb,lora,lord,lorc,lore,lorf]))
										-((gg3[lora,lorb,lorc,lord,lorf,lore]-gg3[lorb,lora,lorc,lord,lorf,lore])-(gg3[lora,lorb,lord,lorc,lorf,lore]-gg3[lorb,lora,lord,lorc,lorf,lore])));


(* ::Code::Initialization::Plain:: *)
(* \[LeftAngleBracket]G^a_\[Mu]\[Nu] D_\[Alpha] D_\[Beta] G^b_\[Rho]\[Sigma]\[RightAngleBracket]=\[LeftAngleBracket]f^def  G^d_\[Mu]\[Nu] G^e_\[Nu]\[Rho] G^f_\[Rho]\[Mu]\[RightAngleBracket](\[Delta]^ab/(4 C_A] C_F D(D-1)(D^2-4)))[(D+2)(((g_\[Nu]\[Alpha] g_\[Beta]\[Rho] g_\[Mu]\[Sigma]-{\[Mu]\[TwoWayRule]\[Nu]})-{\[Alpha]\[TwoWayRule]\[Beta]})-{\[Rho]\[TwoWayRule]\[Sigma]})-(D-4)(((g_\[Nu]\[Alpha] g_\[Beta]\[Rho] g_\[Mu]\[Sigma]-{\[Mu]\[TwoWayRule]\[Nu]})+{\[Alpha]\[TwoWayRule]\[Beta]})-{\[Rho]\[TwoWayRule]\[Sigma]})+4(D-1) g_\[Alpha]\[Beta]( g_\[Mu]\[Rho] g_\[Nu]\[Sigma] - g_\[Mu]\[Sigma] g_\[Nu]\[Rho])] *)
gs0[lora_,lorb_,lorc_,lord_,lore_,lorf_]=(((gg3[lora,lorb,lorc,lord,lore,lorf]-gg3[lorb,lora,lorc,lord,lore,lorf])+(gg3[lora,lorb,lord,lorc,lore,lorf]-gg3[lorb,lora,lord,lorc,lore,lorf]))
										-((gg3[lora,lorb,lorc,lord,lorf,lore]-gg3[lorb,lora,lorc,lord,lorf,lore])+(gg3[lora,lorb,lord,lorc,lorf,lore]-gg3[lorb,lora,lord,lorc,lorf,lore])));


(* ::Code::Initialization::Plain:: *)
g0[lora_,lorb_,lorc_,lord_,lore_,lorf_]=MTD[lorc,lord](MTD[lora,lore]MTD[lorb,lorf]-MTD[lora,lorf]MTD[lorb,lore]);



(* ::Code::Initialization::Plain:: *)
(*---------------------------------------------------*)
ggg[lora_,lorb_,lorc_,lord_,lore_,lorf_,cola_,colb_,colc_]=SUNF[cola,colb,colc]/(2CA^2CF(D-2)(D-1)D)ga0[lora,lorb,lorc,lord,lore,lorf];

gddg[lora_,lorb_,lorc_,lord_,lore_,lorf_,cola_,colb_]=SUNDelta[cola,colb]/(4CA CF(D^2-4)(D-1)D)((D+2)ga0[lora,lorb,lorc,lord,lore,lorf]-(D-4)gs0[lora,lorb,lorc,lord,lore,lorf]+4(D-1)g0[lora,lorb,lorc,lord,lore,lorf]);


(* ::Code::Initialization::Plain:: *)
(* \[LeftAngleBracket](D_\[Alpha] G_\[Mu]\[Nu])^a (D_\[Beta] G_\[Rho]\[Sigma])^b\[RightAngleBracket]=-\[LeftAngleBracket] G^a_\[Mu]\[Nu] (D_\[Alpha] D_\[Beta] G_\[Rho]\[Sigma])^b\[RightAngleBracket] *)
dgdg[lora_,lorb_,lorc_,lord_,lore_,lorf_,cola_,colb_]=-gddg[lorb,lorc,lora,lord,lore,lorf,cola,colb];


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------*)
(*------------------------------------------------------------*)
(* basic propagators *)

(*xprop[x_]=FourierPX[prop[q],{q,x}];*)
xprop[x_]=1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];



(* ::Code::Initialization::Plain:: *)
(* quark propagator wiht background gluon *)
(*-----------------------------------------------------------------------------------------------------------*)

(*(* x \[Rule] -x if translate x \[Rule] -x *)
(*prog[q_,lora_,lorb_,cola_]=-1/2FAD[{q}]GSD[q].GAD[lorb].SUNT[cola].FourDivergence[GSD[q]FAD[{q}],FVD[q,lora]];*)
prog[q_,lora_,lorb_,cola_]=-1/2FAD[{q}]GSD[q] . GAD[lorb] . SUNT[cola] . (FAD[q] GAD[lora]-2 FAD[q,q] FVD[q,lora] GSD[q]);*)


(* ::Code::Initialization::Plain:: *)
(*xprog[x_,lora_,lorb_,cola_]=FourierPX[prog[q,lora,lorb,cola],{q,x}];*)
xprog[x_,lora_,lorb_,cola_]=(GAD[lora] . GAD[lorb] . GSD[x] qfact1[1/32 I^(-D) \[Pi]^(-D/2) qGamma[-1+D/2]] SPD[x,x]^(1-D/2) SUNT[cola]
							+GSD[x] . GAD[lorb] . GAD[lora] qfact1[qfact2[-(1/32) I^(-D) \[Pi]^(-D/2)] qGamma[-1+D/2]] SPD[x,x]^(1-D/2) SUNT[cola]
							+2 FVD[x,lora] FVD[x,lorb] GSD[x] qfact1[-(1/16) I^(-D) \[Pi]^(-D/2) qGamma[D/2]] SPD[x,x]^(-D/2) SUNT[cola]);


(* ::Code::Initialization::Plain:: *)
(* Note, e.g. FourDivergence[SUNT[colb].(FAD[{q}]+FVD[q,lora]GSD[q]FAD[{q}]),FVD[q,lorb]] = (T^colb g^(loralorb) \[Gamma]\[CenterDot]q+T^colb q^lora \[Gamma]^lorb)/q^2+(-2 T^colb q^lora q^lorb \[Gamma]\[CenterDot]q-2 T^colb q^lorb)/(q^2)^2 , The Color matrix not in the Dot product, which cause problem when later tr[] *)


(* ::Code::Initialization::Plain:: *)
(*(*-----------------------------------------------------------------------------------------------------------*)
(* prodg[q_,lora_,lorb_,lorc_,cola_]=I/3FAD[{q}]GSD[q].GAD[lorc].SUNT[cola].FourDivergence[FourDivergence[GSD[q]FAD[{q}],FVD[q,lora]],FVD[q,lorb]]; *)
*)


(* ::Code::Initialization::Plain:: *)
(* xprodg[x_,lora_,lorb_,lorc_,cola_]=FourierPX[prodg[q,lora,lorb,lorc,cola],{q,x}]; *)
xprodg[x_,lora_,lorb_,lorc_,cola_] = (DiracGamma[LorentzIndex[lorc, D], D]*Pair[LorentzIndex[lora, D], LorentzIndex[lorb, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[-1/144*1/(E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]]*SUNT[SUNIndex[cola]] + 
 DiracGamma[LorentzIndex[lora, D], D] . DiracGamma[LorentzIndex[lorc, D], D] . DiracGamma[LorentzIndex[lorb, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[1/(288*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]]*SUNT[SUNIndex[cola]] + 
 DiracGamma[LorentzIndex[lorb, D], D] . DiracGamma[LorentzIndex[lorc, D], D] . DiracGamma[LorentzIndex[lora, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[1/(288*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]]*SUNT[SUNIndex[cola]] + DiracGamma[LorentzIndex[lorb, D], D] . DiracGamma[LorentzIndex[lorc, D], D] . DiracGamma[Momentum[x, D], D]*
  Pair[LorentzIndex[lora, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qGamma[-1 + D/2]/(72*I^D*Pi^(D/2))]*SUNT[SUNIndex[cola]] + 
 DiracGamma[LorentzIndex[lora, D], D] . DiracGamma[LorentzIndex[lorc, D], D] . DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lorb, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qGamma[-1 + D/2]/(72*I^D*Pi^(D/2))]*SUNT[SUNIndex[cola]] + DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorc, D], D] . DiracGamma[LorentzIndex[lorb, D], D]*
  Pair[LorentzIndex[lora, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[-1/144*1/(E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[cola]] + 
 DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorc, D], D] . DiracGamma[LorentzIndex[lora, D], D]*Pair[LorentzIndex[lorb, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qfact2[-1/144*1/(E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[cola]] + 2*DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lora, D], LorentzIndex[lorb, D]]*
  Pair[LorentzIndex[lorc, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[-1/144*1/(E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[cola]] + 
 (2*DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lora, D], Momentum[x, D]]*Pair[LorentzIndex[lorb, D], Momentum[x, D]]*Pair[LorentzIndex[lorc, D], Momentum[x, D]]*
   qfact1[-1/36*qGamma[D/2]/(I^D*Pi^(D/2))]*SUNT[SUNIndex[cola]])/Pair[Momentum[x, D], Momentum[x, D]]^(D/2) + DiracGamma[LorentzIndex[lorc, D], D]*Pair[LorentzIndex[lora, D], Momentum[x, D]]*
  Pair[LorentzIndex[lorb, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[-1/72*(-2 + D)/(E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2] + qGamma[D/2]/(36*I^D*Pi^(D/2))]*
  SUNT[SUNIndex[cola]]);
  


(* ::Code::Initialization::Plain:: *)
(* xprodg1 ~ Contract[ xprodg <DGDG> ] *)
xprodg1[x_,lord_,lore_,lorf_,colc_]=(DiracGamma[LorentzIndex[lorf, D], D]*Pair[LorentzIndex[lord, D], LorentzIndex[lore, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[-((CA - 2*CF)*(-16 + 7*D))/(192*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]]*SUNT[SUNIndex[colc]] + 
 DiracGamma[LorentzIndex[lore, D], D]*Pair[LorentzIndex[lord, D], LorentzIndex[lorf, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[((CA - 2*CF)*(-16 + 7*D))/(192*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]]*SUNT[SUNIndex[colc]] + 
 DiracGamma[LorentzIndex[lord, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lorf, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(CA - 2*CF)/(192*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]]*SUNT[SUNIndex[colc]] + 
 DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lord, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(CA - 2*CF)/(192*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]]*SUNT[SUNIndex[colc]] + 
 DiracGamma[LorentzIndex[lord, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lore, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(-CA + 2*CF)/(192*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]]*SUNT[SUNIndex[colc]] + 
 DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lord, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
  qfact1[qfact2[(-CA + 2*CF)/(192*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]]*SUNT[SUNIndex[colc]] + DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lord, D], LorentzIndex[lore, D]]*
  Pair[LorentzIndex[lorf, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[-((CA - 2*CF)*(-4 + D))/(96*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*
  SUNT[SUNIndex[colc]] + DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lord, D], LorentzIndex[lorf, D]]*Pair[LorentzIndex[lore, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qfact2[((CA - 2*CF)*(-4 + D))/(96*CF*(-2 + D)*(-1 + D)*D*(2 + D)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[colc]] + 
 DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lore, D], D]*Pair[LorentzIndex[lord, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qfact2[(CA - 2*CF)/(96*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[colc]] + DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lord, D], D]*
  Pair[LorentzIndex[lore, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[(CA - 2*CF)/(96*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[colc]] + 
 DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lord, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qfact2[(CA - 2*CF)/(48*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[colc]] + DiracGamma[LorentzIndex[lord, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[Momentum[x, D], D]*
  Pair[LorentzIndex[lorf, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[(CA - 2*CF)/(48*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[colc]] + 
 DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lorf, D], D]*Pair[LorentzIndex[lord, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qfact2[(-CA + 2*CF)/(96*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[colc]] + DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lord, D], D]*
  Pair[LorentzIndex[lorf, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[(-CA + 2*CF)/(96*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[colc]] + 
 DiracGamma[LorentzIndex[lore, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lord, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*
  qfact1[qfact2[(-CA + 2*CF)/(48*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[colc]] + DiracGamma[LorentzIndex[lord, D], D] . DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[Momentum[x, D], D]*
  Pair[LorentzIndex[lore, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[(-CA + 2*CF)/(48*CF*D*(-4 + D^2)*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[colc]])



(* ::Code::Initialization::Plain:: *)
(*
(* proddg[q_,lora_,lorb_,lorc_,lord_,cola_]=1/8FAD[{q}]GSD[q].GAD[lord].SUNT[cola].FourDivergence[FourDivergence[FourDivergence[GSD[q]FAD[{q}],FVD[q,lora]],FVD[q,lorb]],FVD[q,lorc]]; *)
*)


(* ::Code::Initialization::Plain:: *)
(* xproddg[x_,lora_,lorb_,lorc_,lord_,cola_]=FourierPX[proddg[q,lora,lorb,lorc,lord,cola],{q,x}]; *)



(* ::Code::Initialization::Plain:: *)
(*
(*proddga[q_,lora_,lorb_,lorc_,lord_,cola_]=-1/8FourDivergence[FourDivergence[FourDivergence[FAD[{q}]GSD[q],FVD[q,lora]],FVD[q,lorb]],FVD[q,lorc]].GAD[lord].SUNT[cola].GSD[q]FAD[{q}];*)
*)
																																						
(*xproddga[x_,lora_,lorb_,lorc_,lord_,cola_]=FourierPX[proddga[q,lora,lorb,lorc,lord,cola],{q,x}];*)



(* ::Code::Initialization::Plain:: *)
(*-----------------------------------------------------------------------------------------------------------*)
(*
(*progg[q_,lora_,lorb_,lorc_,lord_,cola_,colb_]=-I/4FAD[{q}]GSD[q].GAD[lorb].SUNT[cola].SUNT[colb].FourDivergence[FAD[{q}]GSD[q].GAD[lord].FourDivergence[GSD[q]FAD[{q}],FVD[q,lorc]],FVD[q,lora]];*)
*)


(* ::Code::Initialization::Plain:: *)
(*xprogg[x_,lora_,lorb_,lorc_,lord_,cola_,colb_]=FourierPX[progg[q,lora,lorb,lorc,lord,cola,colb],{q,x}];*)


(* ::Code::Initialization::Plain:: *)
(*-----------------------------------------------------------------------------------------------------------*)
(*prodgdg[q_,lora_,lorb_,lorc_,lord_,lore_,lorf_]=I/9FAD[{q}]GSD[q].GAD[lorc].FourDivergence[FourDivergence[FAD[{q}]GSD[q].GAD[lorf].FourDivergence[FourDivergence[FAD[{q}]GSD[q],FVD[q,lore]],FVD[q,lord]],FVD[q,lora]],FVD[q,lorb]];
xprodgdg[x_,lora_,lorb_,lorc_,lord_,lore_,lorf_]=FourierPX[prodgdg[q,lora,lorb,lorc,lord,lore,lorf],{q,x}];*)

(*
proggg[q_,lora_,lorb_,lorc_,lord_,lore_,lorf_]=1/8FAD[{q}]GSD[q].GAD[lorb].FourDivergence[FAD[{q}]GSD[q].GAD[lord].FourDivergence[FAD[{q}]GSD[q].GAD[lorf].FourDivergence[FAD[{q}]GSD[q],FVD[q,lore]],FVD[q,lorc]],FVD[q,lora]];
xproggg[x_,lora_,lorb_,lorc_,lord_,lore_,lorf_]=FourierPX[proggg[q,lora,lorb,lorc,lord,lore,lorf],{q,x}];
*)



(* ::Code::Initialization::Plain:: *)
(*-----------------------------------------------------------------------------------------------------------*)
(* quark propagator with dimension=6 background glouns  *)
(*xprog3[x_]=FourierPX[prog3[q],{q,x}]//DiracSimplify//QSimplify;*)


(* ::Code::Initialization::Plain:: *)
(*  xprog2[x_,lore_,lorf_,colc_] = xprogg[x_,lora_,lorb_,lorc_,lord_,cola_,colb_] ggg[lora_,lorb_,lorc_,lord_,lore_,lorf_,cola_,colb_,colc_] + xproddg[x_,lora_,lorb_,lorc_,lord_,cola_] gddg[lore_,lorf_,lora_,lorb_,lorc_,lord_,cola_,colc_]  *)
xprog2[x_,lore_,lorf_,colc_]=SUNT[colc](DiracGamma[LorentzIndex[lorf, D], D] . DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lore, D], D]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*
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
(* --------------------- *)


(* ::Code::Initialization::Plain:: *)
Fourquarkggg[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,atr,hv21,hv22,listpole1,listpole2,
holdf=OptionValue[HoldFlavor],diagrams,ndr=OptionValue[AutoNDR],fdir,files,pall,dot,trfirst=OptionValue[TraceFirst],setd=D},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[Fourquarkggg::inderr];
	Abort[]
];



(*-------------------------------*)
(* B A^+ B *)
hv21=v3//ComplexConjugate;
hv22=v4//ComplexConjugate//FCI;
(* ComplexConjugate[sigma] will expand the DiracSigma to gamma matrices, recover it to DiracSigma *)
hv22=(hv22/.Dot->dot)/.f1_ dot[aa___,DiracGamma[LorentzIndex[lor1_,dim___],dim___],DiracGamma[LorentzIndex[lor2_,dim___],dim___],bb___]+f2_ dot[aa___,DiracGamma[LorentzIndex[lor2_,dim___],dim___],DiracGamma[LorentzIndex[lor1_,dim___],dim___],bb___]/;(f1+f2===0):>-2I f1 dot[aa,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],bb];
hv22=hv22/.{null->1,dot->Dot};
(*--------------------------------------------*)		

					
If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];


(*---------------------------------------------------*)


(*diagrams={typeA1,typeA2,typeBS1,typeBS2,typeBS3,typeBS4,typeBS5,typeBS6,typeBS7,typeBS8,typeBS9,typeBS10,typeBS11,typeBS12,typeBS13,typeBS14,typeBS15,typeBS16,
typeC1,typeC2,typeC3,typeC4,typeC5,typeC6,typeC7,typeC8,typeC9,typeC10,typeC11,typeC12,typeC13,typeC14,typeC15,typeC16,typeDS1,typeDS2,typeDS3,typeDS4,typeDS5,typeDS6,typeDS7,typeDS8};
*)

(*---------------------------------------------------*)

If[OptionValue[AutoNDR]===True,atr=TR5,atr=TR];

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

	DistributeDefinitions[v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
	tmp= Plus@@WaitAll[ParallelSubmit[{hv21,hv22,holdf,atr,setd},
									#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr,setd]]&/@diagrams];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,

	If[pall==="External",
	
		DistributeDefinitions[qq,v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
		
		If[fdir==="None",
		
			ParallelSubmit[{hv21,hv22,holdf,atr,setd},#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr,setd]]&/@diagrams
		,
		
			files=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files=("Fourquarkggg_"<>#)&/@files;
			
			ImExport[fdir,
						files,
						{{hv21,hv22,holdf,atr,setd},
						{qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr,setd},
						diagrams}
						]
		]
		(* !!! do not touch tmp before WaitAll[tmp], any function may modifiy the things in EvaluationObject[...] should not be used here. !!! *)
	
	,
	
		
		tmp= Plus@@Map[#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr,setd]&,diagrams];
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]

]


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)
(*-------------The following are generated by algorithem----------------*)
(*----------------------------------------------------------------------*)


(* ::Input::Initialization::Plain:: *)
typeA1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia1,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia1=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprog3[x]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog3[x], hv4, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, xprop[-x], v2, xprog3[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprog3[-x], v1, xprop[x], hv3, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv4, xprog3[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, xprog3[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprog3[-x], v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog3[x]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeA2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia2,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia2=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv4, xprog3[-x], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, xprog3[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprog3[-x], v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprog3[-x], v1, xprop[x], hv4, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, xprog3[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, xprog3[x], hv3, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, xprog3[-x], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, xprop[-x], v2, xprog3[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)


(* ::Input::Initialization::Plain:: *)
typeBS1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia1,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia1=Condensate["ggg"]*contract[tr[str[dot[hv3, xprog[-x, lore, lorf, colc], v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, xprog2[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprog2[-x, lore, lorf, colc], v1, xprop[x], hv4, xprog[-x, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprog2[-x, lore, lorf, colc], v1, xprop[x], hv3, xprog[-x, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia2,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia2=-(Condensate["ggg"]*contract[tr[str[dot[hv4, xprog[-x, lore, lorf, colc], v1, xprop[x], hv3, xprop[-x], v2, xprog2[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv4, xprog2[-x, lore, lorf, colc], v2, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprog[-x, lore, lorf, colc], v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, xprog2[x, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia3,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia3=Condensate["ggg"]*contract[tr[str[dot[hv3, xprog[-x, lore, lorf, colc], v2, xprog2[x, lore, lorf, colc]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprog[-x, lore, lorf, colc], v1, xprop[x], hv4, xprop[-x], v2, xprog2[x, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprog[-x, lore, lorf, colc], v1, xprop[x]]]]*tr[str[dot[hv4, xprog2[-x, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia3=FourierXP[dia3,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=QEvaluate[I ScaleMu^(3(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia4,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia4=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprog[x, lore, lorf, colc]]]]*tr[str[dot[hv4, xprog2[-x, lore, lorf, colc], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprog2[-x, lore, lorf, colc], v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog2[x, lore, lorf, colc], hv4, xprog[-x, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia4=FourierXP[dia4,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=QEvaluate[I ScaleMu^(3(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia5,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia5=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprog2[x, lore, lorf, colc]]]]*tr[str[dot[hv4, xprog[-x, lore, lorf, colc], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog[x, lore, lorf, colc], hv4, xprog2[-x, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprog2[-x, lore, lorf, colc], v1, xprog[x, lore, lorf, colc]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia5=FourierXP[dia5,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=QEvaluate[I ScaleMu^(3(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia6,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia6=Condensate["ggg"]*contract[tr[str[dot[hv3, xprog2[-x, lore, lorf, colc], v2, xprog[x, lore, lorf, colc]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprog[-x, lore, lorf, colc], v1, xprog2[x, lore, lorf, colc], hv4, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprog2[-x, lore, lorf, colc], v1, xprog[x, lore, lorf, colc], hv4, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia6=FourierXP[dia6,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia6=QEvaluate[I ScaleMu^(3(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia7,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia7=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv4, xprog[-x, lore, lorf, colc], v1, xprog2[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprog2[-x, lore, lorf, colc], v1, xprog[x, lore, lorf, colc], hv3, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv4, xprog[-x, lore, lorf, colc], v2, xprog2[x, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia7=FourierXP[dia7,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia7=QEvaluate[I ScaleMu^(3(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia8,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia8=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv4, xprog2[-x, lore, lorf, colc], v1, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprog[-x, lore, lorf, colc], v2, xprop[x]]]]*tr[str[dot[hv4, xprog2[-x, lore, lorf, colc], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, xprog[x, lore, lorf, colc], hv3, xprop[-x], v2, xprog2[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia8=FourierXP[dia8,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia8=QEvaluate[I ScaleMu^(3(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia9,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia9=-(Condensate["ggg"]*contract[tr[str[dot[hv4, xprog[-x, lore, lorf, colc], v1, xprog2[x, lore, lorf, colc], hv3, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) - Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, xprog2[x, lore, lorf, colc], hv3, xprog[-x, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog2[x, lore, lorf, colc]]]]*tr[str[dot[hv4, xprog[-x, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia9=FourierXP[dia9,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia9=QEvaluate[I ScaleMu^(3(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia10,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia10=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprog[x, lore, lorf, colc]]]]*tr[str[dot[hv4, xprop[-x], v1, xprog2[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog[x, lore, lorf, colc], hv4, xprop[-x], v2, xprog2[x, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprog2[-x, lore, lorf, colc], v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia10=FourierXP[dia10,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia10=QEvaluate[I ScaleMu^(3(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia11,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia11=-(Condensate["ggg"]*contract[tr[str[dot[hv4, xprog[-x, lore, lorf, colc], v1, xprop[x], hv3, xprog2[-x, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog[x, lore, lorf, colc]]]]*tr[str[dot[hv4, xprog2[-x, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprog[-x, lore, lorf, colc], v1, xprog2[x, lore, lorf, colc]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia11=FourierXP[dia11,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia11=QEvaluate[I ScaleMu^(3(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia12,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia12=Condensate["ggg"]*contract[tr[str[dot[hv3, xprog2[-x, lore, lorf, colc], v2, xprop[x]]]]*tr[str[dot[hv4, xprog[-x, lore, lorf, colc], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, xprog[-x, lore, lorf, colc], v2, xprog2[x, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprog2[-x, lore, lorf, colc], v1, xprop[x], hv3, xprop[-x], v2, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia12=FourierXP[dia12,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia12=QEvaluate[I ScaleMu^(3(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia13,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia13=-(Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, xprog2[-x, lore, lorf, colc], v2, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, xprog2[x, lore, lorf, colc], hv3, xprop[-x], v2, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, xprog2[-x, lore, lorf, colc], v2, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia13=FourierXP[dia13,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia13=QEvaluate[I ScaleMu^(3(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia14,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia14=-(Condensate["ggg"]*contract[tr[str[dot[hv3, xprog2[-x, lore, lorf, colc], v1, xprop[x], hv4, xprop[-x], v2, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog2[x, lore, lorf, colc], hv4, xprop[-x], v2, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, xprog[x, lore, lorf, colc], hv3, xprog2[-x, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia14=FourierXP[dia14,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia14=QEvaluate[I ScaleMu^(3(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia15,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia15=-(Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, xprog[-x, lore, lorf, colc], v2, xprog2[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog2[x, lore, lorf, colc]]]]*tr[str[dot[hv4, xprop[-x], v2, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog[x, lore, lorf, colc]]]]*tr[str[dot[hv4, xprop[-x], v2, xprog2[x, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia15=FourierXP[dia15,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia15=QEvaluate[I ScaleMu^(3(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeBS16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia16,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia16=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprog2[x, lore, lorf, colc]]]]*tr[str[dot[hv4, xprop[-x], v1, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprog[-x, lore, lorf, colc], v1, xprop[x], hv4, xprog2[-x, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprog2[-x, lore, lorf, colc], v1, xprop[x]]]]*tr[str[dot[hv4, xprog[-x, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia16=FourierXP[dia16,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia16=QEvaluate[I ScaleMu^(3(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)


(* ::Input::Initialization::Plain:: *)
(*typeC1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia1,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia1=-(Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[x, lore, lorf, colc], hv3, xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]);


If[trfirst===True,

    dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia2,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia2=Condensate["ggg"]*contract[tr[str[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]*tr[str[dot[hv4, xprop[-x], v1, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


If[trfirst===True,

    dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia3,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia3=-(Condensate["ggg"]*contract[tr[str[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprop[x], hv3, xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]);


If[trfirst===True,

    dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia3=FourierXP[dia3,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia3=FourierXP[dia3,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia3=QEvaluate[I ScaleMu^(3(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia4,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia4=-(Condensate["ggg"]*contract[tr[str[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprop[x], hv4, xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]);


If[trfirst===True,

    dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia4=FourierXP[dia4,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia4=FourierXP[dia4,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia4=QEvaluate[I ScaleMu^(3(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia5,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia5=Condensate["ggg"]*contract[tr[str[dot[hv3, xprog[-x, lora, lorb, cola], v2, xprop[x]]]]*tr[str[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


If[trfirst===True,

    dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia5=FourierXP[dia5,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia5=FourierXP[dia5,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia5=QEvaluate[I ScaleMu^(3(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia6,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia6=-(Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[x, lore, lorf, colc], hv4, xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]);


If[trfirst===True,

    dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia6=FourierXP[dia6,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia6=FourierXP[dia6,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia6=QEvaluate[I ScaleMu^(3(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia7,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia7=Condensate["ggg"]*contract[tr[str[dot[hv3, xprog[-x, lore, lorf, colc], v1, xprop[x]]]]*tr[str[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


If[trfirst===True,

    dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia7=FourierXP[dia7,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia7=FourierXP[dia7,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia7=QEvaluate[I ScaleMu^(3(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia8,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia8=-(Condensate["ggg"]*contract[tr[str[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb], hv4, xprop[-x], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]);


If[trfirst===True,

    dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia8=FourierXP[dia8,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia8=FourierXP[dia8,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia8=QEvaluate[I ScaleMu^(3(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia9,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia9=Condensate["ggg"]*contract[tr[str[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]*tr[str[dot[hv4, xprog[-x, lore, lorf, colc], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


If[trfirst===True,

    dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia9=FourierXP[dia9,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia9=FourierXP[dia9,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia9=QEvaluate[I ScaleMu^(3(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia10,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia10=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprog[x, lore, lorf, colc]]]]*tr[str[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


If[trfirst===True,

    dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia10=FourierXP[dia10,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia10=FourierXP[dia10,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia10=QEvaluate[I ScaleMu^(3(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia11,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia11=Condensate["ggg"]*contract[tr[str[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb]]]]*tr[str[dot[hv4, xprop[-x], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


If[trfirst===True,

    dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia11=FourierXP[dia11,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia11=FourierXP[dia11,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia11=QEvaluate[I ScaleMu^(3(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia12,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia12=Condensate["ggg"]*contract[tr[str[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb]]]]*tr[str[dot[hv4, xprog[-x, lora, lorb, cola], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


If[trfirst===True,

    dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia12=FourierXP[dia12,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia12=FourierXP[dia12,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia12=QEvaluate[I ScaleMu^(3(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia13,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia13=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprog[x, lora, lorb, cola]]]]*tr[str[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


If[trfirst===True,

    dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia13=FourierXP[dia13,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia13=FourierXP[dia13,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia13=QEvaluate[I ScaleMu^(3(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia14,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia14=-(Condensate["ggg"]*contract[tr[str[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb], hv4, xprog[-x, lora, lorb, cola], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]);


If[trfirst===True,

    dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia14=FourierXP[dia14,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia14=FourierXP[dia14,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia14=QEvaluate[I ScaleMu^(3(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia15,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia15=-(Condensate["ggg"]*contract[tr[str[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb], hv3, xprog[-x, lora, lorb, cola], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]);


If[trfirst===True,

    dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia15=FourierXP[dia15,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia15=FourierXP[dia15,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia15=QEvaluate[I ScaleMu^(3(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization::Plain:: *)
(*typeC16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,trfirst_]:=Block[{x,q,dia16,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia16=-(Condensate["ggg"]*contract[tr[str[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb], hv3, xprop[-x], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]);


If[trfirst===True,

    dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]\[RuleDelayed]If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


    dia16=FourierXP[dia16,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}
,

    dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.{str->SUNTrace}/.tr[0]->0/.contract[0]->0;

    dia16=FourierXP[dia16,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)}/.{tr->atr,contract->Contract}
];

dia16=QEvaluate[I ScaleMu^(3(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)


(* ::Code::Initialization::Plain:: *)
(* set the dimension of gamma-matrices as 4 *)
fdtypeC1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC1[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC2[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC3[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC4[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC5[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC6[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC7[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC8[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC9[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC10[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC11[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC12[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC13[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC14[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC15[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]
fdtypeC16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=typeC16[qq,{v1,v2,{a1,b1,c1,d1}},{hv3,hv4,{a2,b2,c2,d2}},holdf,atr,4]


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)


(* ::Input::Initialization::Plain:: *)
typeC1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia1,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia1=Condensate["ggg"]*contract[str[tr[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb]]]]*str[tr[dot[hv4, xprog[-x, lora, lorb, cola], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia2,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia2=-(Condensate["ggg"]*contract[str[tr[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb], hv4, xprop[-x], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]);


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia3,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia3=Condensate["ggg"]*contract[str[tr[dot[hv3, xprop[-x], v1, xprog[x, lore, lorf, colc]]]]*str[tr[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia3=FourierXP[dia3,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=QEvaluate[I ScaleMu^(3(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia4,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia4=-(Condensate["ggg"]*contract[str[tr[dot[hv4, xprop[-x], v1, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[x, lore, lorf, colc], hv3, xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]);


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia4=FourierXP[dia4,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=QEvaluate[I ScaleMu^(3(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia5,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia5=Condensate["ggg"]*contract[str[tr[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]*str[tr[dot[hv4, xprop[-x], v1, xprog[x, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia5=FourierXP[dia5,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=QEvaluate[I ScaleMu^(3(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia6,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia6=-(Condensate["ggg"]*contract[str[tr[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprop[x], hv4, xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]);


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia6=FourierXP[dia6,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia6=QEvaluate[I ScaleMu^(3(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia7,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia7=-(Condensate["ggg"]*contract[str[tr[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb], hv4, xprog[-x, lora, lorb, cola], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]);


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia7=FourierXP[dia7,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia7=QEvaluate[I ScaleMu^(3(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia8,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia8=Condensate["ggg"]*contract[str[tr[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]*str[tr[dot[hv4, xprog[-x, lore, lorf, colc], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia8=FourierXP[dia8,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia8=QEvaluate[I ScaleMu^(3(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia9,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia9=-(Condensate["ggg"]*contract[str[tr[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprop[x], hv3, xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]);


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia9=FourierXP[dia9,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia9=QEvaluate[I ScaleMu^(3(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia10,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia10=Condensate["ggg"]*contract[str[tr[dot[hv3, xprog[-x, lore, lorf, colc], v1, xprop[x]]]]*str[tr[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia10=FourierXP[dia10,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia10=QEvaluate[I ScaleMu^(3(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia11,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia11=Condensate["ggg"]*contract[str[tr[dot[hv3, xprog[-x, lora, lorb, cola], v2, xprop[x]]]]*str[tr[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia11=FourierXP[dia11,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia11=QEvaluate[I ScaleMu^(3(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia12,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia12=Condensate["ggg"]*contract[str[tr[dot[hv3, xprop[-x], v2, xprog[x, lora, lorb, cola]]]]*str[tr[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia12=FourierXP[dia12,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia12=QEvaluate[I ScaleMu^(3(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia13,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia13=-(Condensate["ggg"]*contract[str[tr[dot[hv3, xprop[-x], v1, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[x, lore, lorf, colc], hv4, xprog[-x, lorc, lord, colb], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]);


dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia13=FourierXP[dia13,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia13=QEvaluate[I ScaleMu^(3(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia14,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia14=-(Condensate["ggg"]*contract[str[tr[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb], hv3, xprog[-x, lora, lorb, cola], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]);


dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia14=FourierXP[dia14,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia14=QEvaluate[I ScaleMu^(3(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia15,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia15=-(Condensate["ggg"]*contract[str[tr[dot[hv4, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb], hv3, xprop[-x], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]);


dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia15=FourierXP[dia15,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia15=QEvaluate[I ScaleMu^(3(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeC16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia16,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia16=Condensate["ggg"]*contract[str[tr[dot[hv3, ggg[lora, lorb, lorc, lord, lore, lorf, cola, colb, colc]*xprog[-x, lore, lorf, colc], v1, xprog[x, lorc, lord, colb]]]]*str[tr[dot[hv4, xprop[-x], v2, xprog[x, lora, lorb, cola]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.tr[gs_]:>If[setd===D,atr[gs],Tr[gs/.{DiracGamma[LorentzIndex[lo_,D],D]:>DiracGamma[LorentzIndex[lo]],DiracGamma[Momentum[mo_,D],D]:>DiracGamma[Momentum[mo]],LorentzIndex[lo_,D]:>LorentzIndex[lo]}]]/.{str->SUNTrace,contract->Contract};


dia16=FourierXP[dia16,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia16=QEvaluate[I ScaleMu^(3(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)


(* ::Input::Initialization::Plain:: *)
typeDS1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia1,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia1=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprop[x]]]]*tr[str[dot[hv4, xprodg[-x, lord, lore, lorf, colc], v1, xprodg1[x, lord, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprodg[-x, lord, lore, lorf, colc], v1, xprodg1[x, lord, lore, lorf, colc], hv3, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprodg[x, lord, lore, lorf, colc]]]]*tr[str[dot[hv4, xprop[-x], v2, xprodg1[x, lord, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(3(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeDS2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia2,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia2=-(Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprodg[x, lord, lore, lorf, colc], hv4, xprop[-x], v2, xprodg1[x, lord, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x], hv4, xprodg[-x, lord, lore, lorf, colc], v2, xprodg1[x, lord, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprodg[x, lord, lore, lorf, colc]]]]*tr[str[dot[hv4, xprodg1[-x, lord, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(3(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeDS3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia3,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia3=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprodg1[x, lord, lore, lorf, colc]]]]*tr[str[dot[hv4, xprodg[-x, lord, lore, lorf, colc], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprodg[-x, lord, lore, lorf, colc], v2, xprodg1[x, lord, lore, lorf, colc]]]]*tr[str[dot[hv4, xprop[-x], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprodg[-x, lord, lore, lorf, colc], v1, xprop[x], hv3, xprop[-x], v2, xprodg1[x, lord, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia3=FourierXP[dia3,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=QEvaluate[I ScaleMu^(3(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeDS4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia4,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia4=Condensate["ggg"]*contract[tr[str[dot[hv3, xprodg1[-x, lord, lore, lorf, colc], v2, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v1, xprodg[x, lord, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, xprodg[x, lord, lore, lorf, colc], hv3, xprop[-x], v2, xprodg1[x, lord, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprodg[-x, lord, lore, lorf, colc], v1, xprop[x]]]]*tr[str[dot[hv4, xprop[-x], v2, xprodg1[x, lord, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia4=FourierXP[dia4,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=QEvaluate[I ScaleMu^(3(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeDS5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia5,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia5=Condensate["ggg"]*contract[tr[str[dot[hv3, xprodg1[-x, lord, lore, lorf, colc], v2, xprop[x]]]]*tr[str[dot[hv4, xprodg[-x, lord, lore, lorf, colc], v1, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv3, xprodg[-x, lord, lore, lorf, colc], v1, xprodg1[x, lord, lore, lorf, colc], hv4, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprodg[-x, lord, lore, lorf, colc], v1, xprop[x], hv3, xprodg1[-x, lord, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia5=FourierXP[dia5,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=QEvaluate[I ScaleMu^(3(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeDS6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia6,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia6=-(Condensate["ggg"]*contract[tr[str[dot[hv3, xprodg[-x, lord, lore, lorf, colc], v1, xprop[x], hv4, xprodg1[-x, lord, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprodg[x, lord, lore, lorf, colc], hv4, xprodg1[-x, lord, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprodg[-x, lord, lore, lorf, colc], v1, xprop[x]]]]*tr[str[dot[hv4, xprodg1[-x, lord, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia6=FourierXP[dia6,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia6=QEvaluate[I ScaleMu^(3(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeDS7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia7,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia7=Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v2, xprodg1[x, lord, lore, lorf, colc]]]]*tr[str[dot[hv4, xprop[-x], v1, xprodg[x, lord, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, xprop[x], hv3, xprodg[-x, lord, lore, lorf, colc], v2, xprodg1[x, lord, lore, lorf, colc]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprodg[-x, lord, lore, lorf, colc], v1, xprodg1[x, lord, lore, lorf, colc]]]]*tr[str[dot[hv4, xprop[-x], v2, xprop[x]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia7=FourierXP[dia7,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia7=QEvaluate[I ScaleMu^(3(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization::Plain:: *)
typeDS8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_,setd_]:=Block[{x,q,dia8,lora,lorb,lorc,lord,lore,lorf,cola,colb,colc,tr,str,contract,dot},


dia8=-(Condensate["ggg"]*contract[tr[str[dot[hv3, xprodg[-x, lord, lore, lorf, colc], v1, xprop[x], hv4, xprop[-x], v2, xprodg1[x, lord, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["ggg"]*contract[tr[str[dot[hv4, xprop[-x], v1, xprodg[x, lord, lore, lorf, colc], hv3, xprodg1[-x, lord, lore, lorf, colc], v2, xprop[x]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["ggg"]*contract[tr[str[dot[hv3, xprop[-x], v1, xprop[x]]]]*tr[str[dot[hv4, xprodg[-x, lord, lore, lorf, colc], v2, xprodg1[x, lord, lore, lorf, colc]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia8=FourierXP[dia8,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia8=QEvaluate[I ScaleMu^(3(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
