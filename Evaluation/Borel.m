(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)




Borel::usage = 
	"Borel[expr,{momentum,parameter}] is Borel transformation"


Borel::evawarn = "The input seems hasn't been evaluated. Evaluate or refine the expression."
Borel::polewarn="Non-local pole ~log/\[Epsilon] involved, which will be discard."
Borel::imwarn="Imaginary comes other than s^n log^m involved, the result may not desired."
Borel::argerr="No `1` involved in expression."


Begin["`Private`Borel`"]



Options[Borel] ={
	Renormalization->"Auto"}
(* ,deriv=OptionValue[Derivate],sub=OptionValue[Subtraction]
Derivate->deriv,Subtraction->sub, *)



Borel[expr_,{p_Symbol,t_,s0_},OptionsPattern[]]/;FreeQ[expr,p]:=Message[Borel::argerr,p]
Borel[expr_,{p_Symbol,t_},OptionsPattern[]]/;FreeQ[expr,p]:=Message[Borel::argerr,p]
Borel[expr_,{p_Symbol,{t_,nn_}s0_},OptionsPattern[]]/;FreeQ[expr,p]:=Message[Borel::argerr,p]
Borel[expr_,{p_Symbol,{t_,nn_}},OptionsPattern[]]/;FreeQ[expr,p]:=Message[Borel::argerr,p]


(* Integrate from 0 to s0 *)
Borel[expr:Except[_?ListQ],{p_Symbol,t_,s0_},OptionsPattern[]]/;MatchQ[t,_Symbol|_Positive]:=Block[
{tmp,re=OptionValue[Renormalization]},

tmp=expr//FCI;

If[!FreeQ[tmp,Momentum[aa_Plus]/;!FreeQ[aa,p]],Message[Borel::evawarn];Abort[]];

tmp=Borel[QGather[tmp,p,ShowasTable->True],{p,t,s0},Renormalization->re]
]


(*-----------*)
Borel[expr:Except[_?ListQ],{p_Symbol,{t_,n_Integer},s0_},OptionsPattern[]]/;MatchQ[t,_Symbol|_Positive]&&Positive[n]:=
Borel[expr Pair[Momentum[p,D],Momentum[p,D]]^n,{p,t,s0},Renormalization->OptionValue[Renormalization]]




(* s0 -> Infinity *)
Borel[expr:Except[_?ListQ],{p_Symbol,t_},OptionsPattern[]]/;MatchQ[t,_Symbol|_Positive]:=Block[{tmp,re=OptionValue[Renormalization]},

tmp=expr//FCI;

If[!FreeQ[tmp,Momentum[aa_Plus]/;!FreeQ[aa,p]],Message[Borel::evawarn];Abort[]];

tmp=Borel[QGather[tmp,p,ShowasTable->True],{p,t},Renormalization->re]


]

(*-----------*)
Borel[expr:Except[_?ListQ],{p_Symbol,{t_,n_Integer}},OptionsPattern[]]/;MatchQ[t,_Symbol|_Positive]&&Positive[n]:=
Borel[expr Pair[Momentum[p,D],Momentum[p,D]]^n ,{p,t},Renormalization->OptionValue[Renormalization]]




(* s0 -> Infinity *)
(* for input match the output of QGather *)
Borel[expr_List,{p_Symbol,{t_,n_Integer}},OptionsPattern[]]/;(And@@(MatchQ[#,{_,_List}]&/@expr))&&MatchQ[t,_Symbol|_Positive]&&Positive[n]:=
Borel[({1,Pair[Momentum[p,D],Momentum[p,D]]^n}#)&/@expr ,{p,t},Renormalization->OptionValue[Renormalization]]

(*-----------*)
Borel[expr_List,{p_Symbol,t_},OptionsPattern[]]/;(And@@(MatchQ[#,{_,_List}]&/@expr))&&MatchQ[t,_Symbol|_Positive]:=
Borel[expr,{p,t,Infinity},Renormalization->OptionValue[Renormalization]]




(* Integrate from 0 to s0 *)
(* for input match the output of QGather *)

Borel[expr_List,{p_Symbol,{t_,n_Integer},s0_},OptionsPattern[]]/;(And@@(MatchQ[#,{_,_List}]&/@expr))&&MatchQ[t,_Symbol|_Positive]&&Positive[n]:=
Borel[({1,Pair[Momentum[p,D],Momentum[p,D]]^n}#)&/@expr,{p,t,s0},Renormalization->OptionValue[Renormalization]]

(*-----------*)
Borel[expr_List,{p_Symbol,t_,s0_},OptionsPattern[]]/;(And@@(MatchQ[#,{_,_List}]&/@expr))&&MatchQ[t,_Symbol|_Positive]:=Block[
{tmp,tmp1,re=OptionValue[Renormalization],s,null,null1,null0,r,v2,tmp2},

tmp=expr//FCI;



If[!FreeQ[tmp,Momentum[aa_Plus]/;!FreeQ[aa,p]],Message[Borel::evawarn];Abort[]];
If[!FreeQ[tmp,Log[_]Power[Epsilon,_?Negative]|Power[Log[_],_]Power[Epsilon,_?Negative]],Message[Borel::polewarn]];
(*If[FreeQ[tmp,p],Message[Borel::argerr];Abort[]];*)


(* remove 1/epsiln pole *)
tmp=Replace[tmp/.Power[Epsilon,_?Negative]->0,aa_List:>DeleteCases[aa,0],{2}];


(*--- p^2 -> s; extract imaginary part; Borel transformation ---*)
tmp=Replace[tmp,aa_List:>Expand[null1(Plus@@aa)/.Pair[Momentum[p,___],Momentum[p,___]]->s]+null0+null0^2,{2}];

(* treat 4Pi mu^2 ~ mu^2 , since the result differ by rescale mu (or tau) *)
tmp={#[[1]],#[[2]]/.{Log[-s/(4Pi ScaleMu^2)]:>Log[-s/v2],Log[s/(4Pi ScaleMu^2)]:>Log[s/v2],Log[-s/(ScaleMu^2)]:>Log[-s/v2],Log[s/(ScaleMu^2)]:>Log[s/v2]}}&/@tmp;

(*--- label different terms: s^n log^m <-> null[0,_] ; log^m/s^n <-> null[1,_] ; 1/s^n <-> null[2,_] ; others <-> null[3,_] ---*)
tmp={#[[1]],Replace[#[[2]]/.{aa_ Power[s,nn_?Positive] Log[-s/v2]^mm_:>null[0,aa s^nn Log[-s/v2]^mm],
							aa_ Power[s,nn_?Positive] Log[-s/v2]:>null[0,aa s^nn Log[-s/v2]],
							aa_ Power[s,nn_?Negative] Log[-s/v2]^mm_:>null[1,aa s^nn Log[-s/v2]^mm],
							aa_ Power[s,nn_?Negative] Log[-s/v2]:>null[1,aa s^nn Log[-s/v2]]},
							
								{aa_ Power[Log[-s/v2],mm_]:>null[0,aa Log[-s/v2]^mm],
								aa_ Log[-s/v2]:>null[0,aa Log[-s/v2]],
								aa_ Power[s,nn_?Negative]:>null[2,aa s^nn]},{1}]}&/@tmp;


tmp={#[[1]], #[[2]]-(#[[2]]/.null[_,_]->0) + null[3,#[[2]]/.null[_,_]->0]}&/@tmp;


If[!FreeQ[tmp,null[3,aa_/;!FreeQ[aa,Complex]]],Message[Borel::imwarn]];


(* combine same type of term *)
tmp=tmp//.null[ii_,aa_]+null[ii_,bb_]:>null[ii,aa+bb];


(* extrace imaginary part for null[0,_] and null[3,_] before evaluation to save the time; for null[1,_] and null[2,_], no significant improvement *)
tmp=tmp/.null[0,aa_]:>null[0,(aa/.Log[-s/v2]:>Log[s/v2]-Pi I)//Expand];
tmp=tmp/.{null[0,aa_]:>null[0,(aa-(aa/.Complex[_,_]->0))/.{bb_ Complex[rr_,ii_]:>bb ii,Complex[rr_,ii_]:>ii}],
		null[3,aa_]:>null[3,(aa-(aa/.Complex[_,_]->0))/.{bb_ Complex[rr_,ii_]:>bb ii,Complex[rr_,ii_]:>ii}]};
	
	
(*--- Borel transform each type of term ---*)
tmp=tmp/.{null[0,aa_]:>(tmp1=IntegrateLog[Exp[-t s]aa,s];null[0,Limit[tmp1,s->s0,Direction->1,Assumptions->{v2>0,t>0}]-Limit[tmp1,s->0,Direction->-1,Assumptions->{v2>0,t>0}]]),
																	
		null[2,aa_]:>null[2,Expand[aa/.Power[s,nn_Integer]/;Negative[nn] :> -I Pi/((-nn-1)!) (-t)^(-nn-1)]],(* integrate 1/s^n e^(-st) from 0- to 0+ *)
	
		null[3,aa_]:>(tmp1=IntegrateLog[Exp[-t s]aa,s];null[3,Limit[tmp1,s->s0,Direction->1,Assumptions->{v2>0,t>0}]-Limit[tmp1,s->0,Direction->-1,Assumptions->{v2>0,t>0}]])};


(* for log^m/s^n, Im[log^m/s^n] divergence at s=0, to get the result, Borel transform it (times Pi I to recover imaginary part) then -{Integrate from s0 to Infinity} *)
If[s0===Infinity,
	tmp=tmp/.null[1,aa_]:>null[1,aa/.{s^nn_ Log[-s/v2] :>Pi I (-1)^(-nn)(t^(-nn-1)/Gamma[-nn](-Log[t v2] +PolyGamma[-nn])),
	
					s^nn_ Log[-s/v2]^2 :>Pi I (-1)^(-nn)(t^(-nn-1)/Gamma[-nn](Log[t v2]^2-2PolyGamma[-nn]Log[t v2] +PolyGamma[-nn]^2-PolyGamma[1,-nn])),
										
															(* arg and sign for dQ^2 <-> -ds; times 1/2 since other integrations only account the part from 0- to s0 above the x-axis *)
					s^nn_ Log[-s/v2]^mm_Integer/;mm>2 :>1/2 Limit[-Integrate[(Log[r/v2]+I theta)^mm (-1)^(-nn)(r Exp[I theta])^(1+nn) Normal[Series[Exp[r Exp[I theta]t],{r,0,-nn-1}]] I ,{theta,Pi,-Pi}]+
																	IntegrateLog[(Log[s/v2]-I Pi)^mm s^nn Exp[-s t]-(Log[s/v2]+I Pi)^mm s^nn Exp[-s t]//Expand,{s,r,s0},Assumptions->{s0>r,r>0,t>0}],
																	r->0,Direction->-1,Assumptions->{v2>0,t>0}]}
					]
		
(* add assumptions {v2>0,s0>r,r>0,t>0} make the Integrate slow, but t>0 is necessary when s0 = Infinity *)		
,
	tmp=tmp/.null[1,aa_]:>null[1,aa/.{s^nn_ Log[-s/v2] :>Pi I (-1)^(-nn)(t^(-nn-1)/Gamma[-nn](-Log[t v2] +PolyGamma[-nn]))+Pi I t^(-nn-1)Gamma[1+nn,s0 t],
	
					s^nn_ Log[-s/v2]^2 :>Pi I (-1)^(-nn)(t^(-nn-1)/Gamma[-nn](Log[t v2]^2-2PolyGamma[-nn]Log[t v2] +PolyGamma[-nn]^2-PolyGamma[1,-nn]))-
											IntegrateLog[Exp[-s t](-Pi I)2 Log[s/v2]s^nn,{s,s0,Infinity},Assumptions->t>0],
										
															(* arg and sign for dQ^2 <-> -ds; times 1/2 since other integrations only account the part from 0- to s0 above the x-axis *)
					s^nn_ Log[-s/v2]^mm_Integer/;mm>2 :>1/2 Limit[-Integrate[(Log[r/v2]+I theta)^mm (-1)^(-nn)(r Exp[I theta])^(1+nn) Normal[Series[Exp[r Exp[I theta]t],{r,0,-nn-1}]] I ,{theta,Pi,-Pi}]+
																	IntegrateLog[(Log[s/v2]-I Pi)^mm s^nn Exp[-s t]-(Log[s/v2]+I Pi)^mm s^nn Exp[-s t]//Expand,{s,r,s0},Assumptions->{s0>r,r>0}],
																	r->0,Direction->-1,Assumptions->{v2>0,t>0}]}
					]
];


tmp=tmp/.{null[0,aa_]:>aa,null[1,aa_]:>((aa-(aa/.Complex[_,_]->0))/.{bb_ Complex[rr_,ii_]:>bb ii,Complex[rr_,ii_]:>ii}),
		null[2,aa_]:>((aa-(aa/.Complex[_,_]->0))/.{bb_ Complex[rr_,ii_]:>bb ii,Complex[rr_,ii_]:>ii}),null[3,aa_]:>aa};


(* renormalization improvement v2 -> 1/t; or set v2 to specifical value *)
If[ToLowerCase[ToString[re]]==="auto",

	tmp=tmp/.{Log[1/v2]->Log[t],Log[v2]->-Log[t]}/.v2->1/t
,
	If[MatchQ[re,Rule[_,_]],tmp=tmp/.v2->re[[2]]]
];


(*---            ---*)
tmp={#[[1]],Expand[1/Pi #[[2]]/.{null1->1,null0->0}]//Simplify}&/@tmp;


If[!FreeQ[tmp,Condensate],
	tmp={#[[1]],Collect[Collect[#[[2]],t],Condensate[_]]}&/@tmp,
	tmp={#[[1]],Collect[#[[2]],t]}&/@tmp
];


If[Length[tmp]==1&&tmp[[1,1]]==1,tmp[[1,2]],tmp]



]



End[]
