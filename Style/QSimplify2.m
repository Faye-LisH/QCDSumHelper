(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



(* to do:
the Current version can't deal with the terms like Gamma[n - qdelta + D/2] properly,
 *)


QSimplify2::usage= "QSimplify2[expr_,ruless___] is a function to simplify the experssion by sorting the Lorentz structure and simplify the Gamma functions.";


Begin["`Private`QSimplify2`"]


Options[QSimplify2] = {
	Separate->{"null"},
	SimplifyGamma->True}



(*---------------------------------------------------------------*)

QSimplify2[expr_,OptionsPattern[]]:=Block[{tmp,identity,ruless=OptionValue[Separate],rules0,rules,null0,null,nogamm,dot,gamm,sgamma=OptionValue[SimplifyGamma]},

tmp=expr//FCI//Expand;
tmp=tmp/.Dot->dot;


If[Head[ruless]=!=List,ruless=List[ruless]];


rules0=If[And@@(MatchQ[#,_Rule]&/@ruless),#[[1]]&/@ruless, ruless];
rules=Join[{Pair,FeynAmpDenominator,DiracGamma,DiracSigma,Eps,SUNT,SUNTF,SUND,SUNF,SUNDelta,dot,ColorDelta,FlavorDelta},rules0]//DeleteDuplicates;(* the terms with these heads should not be including in qfact1[] *)
rules=DeleteCases[rules,_String];(* delete the "null" in the rules *)

(*----------------------------------------------------------------*)

If[Head[tmp]===Plus,
	tmp=List@@tmp;

	(*--- seprate commutative and non-commutative parts ---*)
	tmp={Replace[#,aa_/;!FreeQ[aa,Alternatives@@rules]->1,{1}],Replace[#,aa_/;FreeQ[aa,Alternatives@@rules]->1,{1}]}&/@tmp;
	tmp=Gather[tmp,Last[#1]==Last[#2]&];(* treat different non-commutative parts as different terms *)


	tmp=(Transpose[#]&/@tmp)/.{qfact1->Identity,dot->Dot};
	(**--- to a list of {{a list of commutative parts},non-commutative part} ---**)
	tmp=qfact1[(Plus@@#[[1]])//Expand] #[[2,1]]&/@tmp;
];


(*--------------------------------------------------------------------------------*)
If[Depth[tmp]<=2,
	(* doesn't need to simplify *)
	Replace[tmp/.dot->Dot, List->Plus,{1},Heads->True]
	
,

	tmp=(tmp/.dot->Dot)/.qfact1[aa_]:>qfact1[aa/.qfact2[bb_]qfact2[cc_]:>qfact2[Simplify[bb cc/.qfact2->Identity]]];
	
(*---    ------------------------------------------------ Simplify Gamma functions ---------------------------------------------    ---*)
	
	(* isolate the terms without Gamma functions *)
	tmp=tmp/.qfact1[aa_Plus]:>(nogamm=(aa/.Power[qGamma[___],___]->0)/.qGamma[___]->0; qfact1[qfact2[nogamm] + aa - nogamm]);
	
	
	tmp=tmp/.qfact1[aa_Plus]:>qfact1[List@@aa];
	tmp=tmp/.qGamma[aa_]:>qGamma[Expand[aa]]/.qGamma[aa_Integer]/;Positive[aa]:>Gamma[aa];
	
	
(*---     -----------------------------------  prepareing for simplify -------------------------------------------     ---*)

		
(**--- for each term in qfact1, write it to {Gamma functions without prefactors, the whole term } ---**)
(* inside the gamma function, seprate the numerical terms and symbolic terms, i.e. qGamma[term] \[Rule] gamm[symbolic, numerical, exponent] *)		
	tmp=tmp/.qfact1[aa_List]:>qfact1[({#/(#/.qGamma[___]->1), #}/.Power[qGamma[cc_],pow_]:>gamm[cc+null0,pow]/.qGamma[cc_]:>gamm[cc+null0,1]/.{
																	gamm[cc_,pow_]:>gamm[Replace[cc,{_Rational->0,_Integer->0},{1}],cc-Replace[cc,{_Rational->0,_Integer->0},{1}],pow]
																	}/.null0->0
	
(*	/.{Power[qGamma[cc_Plus],dd_]:>gamm[Replace[cc,{_Rational->0,_Integer->0,_Complex},{1}],cc-Replace[cc,{_Rational->0,_Integer->0,_Complex},{1}],dd],
																	Power[qGamma[cc:Except[_Plus|_Rational|_Integer|_Complex]],dd_]:>gamm[cc,0,dd]
																	}/.{
																		qGamma[cc_Plus]:>gamm[Replace[cc,{_Rational->0,_Integer->0,_Complex},{1}],cc-Replace[cc,{_Rational->0,_Integer->0,_Complex},{1}],1],
																		qGamma[cc:Except[_Plus|_Rational|_Integer|_Complex]]:>gamm[cc,0,1]
																		}*)
										)&/@aa];
										
										(*{Power[qGamma[cc_Plus],dd_]:>gamm[Replace[cc,{_Rational->0,_Integer->0},{1}],cc-Replace[cc,{_Rational->0,_Integer->0},{1}],dd],
																	Power[qGamma[cc:Except[_Plus]],dd_]:>gamm[cc,0,dd]
																	}/.{
																		qGamma[cc_Plus]:>gamm[Replace[cc,{_Rational->0,_Integer->0},{1}],cc-Replace[cc,{_Rational->0,_Integer->0},{1}],1],
																		qGamma[cc:Except[_Plus]]:>gamm[cc,0,1]
																		}*)
						
						
(**---   seprate the gamma functions i.e. {{the gamma functions in the numerator, the gamma functions in the denominator}, the whole term} for each term in qfactor1    ---**)
	tmp=tmp/.qfact1[aa_List]:>qfact1[{{#[[1]]/.gamm[_,_,_?Negative]->1,#[[1]]/.gamm[_,_,_?Positive]->1},#[[2]]}&/@aa];
	
	
(* recall that gamm have the form: gamm[symbolic term, numerical term, exponent];
 for each term in qfactor1, modify it to the form:
 {{postive power of gamma functions but omit the numeric, negative power of gamma functions but omit the numeric},
 {postive power of gamma functions but only keep the gamma function with minimum numeric and omit the exponent,
  negative power of gamma functions but only keep the gamma function with maximum numeric and omit the exponent },
  the whole term}
  , 
  e.g.
  {{gamm[D/2-qdelta,0,3],gamm[3D/2,0,-2]},
  {gamm[D/2-qdelta,1/2,_],gamm[3D/2,-1,_]},
  factor1 gamm[D/2-qdelta,1/2,2]gamm[D/2-qdelta,3/2,1]/(gamm[3D/2,-1,-1]gamm[3D/2,-2,-1])};
  
  if two terms with same {postive power of gamma function but omit the numerical term, negative power of gamma function but omit the numerical term},
   it's possible to combine them to a single term by recursive property of Gamma function *)
  
	tmp=tmp/.qfact1[aa_List]:>qfact1[{(#[[1]]//.gamm[symbol_,_,exponent1_]gamm[symbol_,_,exponent2_]:>gamm[symbol,0,exponent1+exponent2])/.gamm[symbol_,_,exponent_]:>gamm[symbol,0,exponent],
									{#[[1,1]]//.gamm[symbol_,numeric1_,_]gamm[symbol_,numeric2_,_]/;numeric1>=numeric2:>gamm[symbol,numeric2,0],(* to generate a list of rules for replacement later *)
									#[[1,2]]//.gamm[symbol_,numeric1_,_]gamm[symbol_,numeric2_,_]/;numeric1>=numeric2:>gamm[symbol,numeric1,0]},
									#[[2]]
									}&/@aa];
									
							
(**--- gather the gamma functions by {postive power of gamma functions but omit the numeric, negative power of gamma functions but omit the numeric} ---**)							
	tmp=tmp/.qfact1[aa_List]:>qfact1[Gather[aa,First[#1]===First[#2]&]];
(*------------------------ Combine the Gamma functions ---------------------------*)	

	If[Head[tmp]===List,	
		tmp=Total[tmp/.qfact1[aa_List]:>qfact1[(qSimplify3[#,sgamma]&/@aa)//Total]]
		,	
		tmp=tmp/.qfact1[aa_List]:>qfact1[(qSimplify3[#,sgamma]&/@aa)//Total]
	]
]

]


(*---------------------------------------------------------------*)


qSimplify3[expr_List,sgamma_]:=Block[{tmp,list,nlist1,nlist2,dlist1,dlist2,result,rule1,rule2,gamm,e},

If[Length[expr]===1,
	
	tmp=Expand[expr[[1,3]]]/.qfact2[aa_]qfact2[bb_]:>qfact2[Simplify[aa bb]];
	tmp=tmp/.gamm[aa_,bb_,cc_]:>qGamma[aa+bb]^cc

,
(*--- gather the gamma functions by recursive property of gamma function ---*)
(* see lines 107 - 132 *)

	tmp=expr//Transpose;
(**--- get the symbolic terms in gamm[] ---**)
	nlist1=tmp[[1,1,1]];
	nlist1=If[Head[nlist1]===Times,(List@@nlist1)/.gamm[aa_,_,_]:>aa,{nlist1/.gamm[aa_,_,_]:>aa}];(* symbolic term in numerator *)
	dlist1=tmp[[1,1,2]];
	dlist1=If[Head[dlist1]===Times,(List@@dlist1)/.gamm[aa_,_,_]:>aa,{dlist1/.gamm[aa_,_,_]:>aa}];(* symbolic term in denominator *)

(**---  find the maximun, minimun of numerical terms in gamm[]  ---**)
	list=tmp[[2]]//Transpose;

	nlist2=(Cases[list[[1]],gamm[#,aa_,_]:>aa,Infinity]//Min)&/@nlist1; (* minimum of numerical terms for corresponding symbolic term *)
	dlist2=(Cases[list[[2]],gamm[#,aa_,_]:>aa,Infinity]//Max)&/@dlist1; (* maximun of numerical terms for corresponding symbolic term *)

(*---   replace gamm[sym_,num_,exp_] by recursive property of gamma function  ---*)
	
	rule1=(gamm[#[[1]],bb_,exp_?Positive]/;IntegerQ[bb-#[[2]]]:>Product[#[[1]]+i,{i,#[[2]],bb-1}]^exp gamm[#[[1]],#[[2]],exp])&/@Transpose[{nlist1,nlist2}];
	rule2=(gamm[#[[1]],bb_,exp_?Negative]/;IntegerQ[bb-#[[2]]]:>gamm[#[[1]],#[[2]],exp]Product[#[[1]]+i,{i,bb,#[[2]]-1}]^(-exp))&/@Transpose[{dlist1,dlist2}];

	result=tmp[[3]]/.rule1/.rule2;

(*---  recover the expression  ---*)

	result={(#/(#/._gamm->1))/.gamm[aa_,bb_,exp_]:>qGamma[aa+bb]^exp, Simplify[#/._gamm->1]}&/@result;
	result=Transpose[result]/.{Power[Complex[0,im_],exp_]:>im^exp Exp[Pi/2 I exp],Power[-1,exp_]:>Exp[Pi I exp]}/.{
							factor_ Power[E,exp_] Complex[0,im_]/;!FreeQ[exp,_Complex]:>factor im Exp[exp+Pi I/2 //Simplify],
							factor_ Power[Complex[0,im1_],exp_] Complex[0,im2_]:>factor im1^exp im2 Exp[(exp+1)Pi I/2 //Simplify]};(* unify the (-1)^n or I^n to ~ e^(Pi I n) *)
	
(*	result=Transpose[result]/.{Power[Complex[0,im_],exp_]:>im e[Pi/2 I exp],Power[-1,exp_]:>e[Pi I exp]}/.{
							factor_ Power[E,exp_] Complex[0,im_]/;!FreeQ[exp,_Complex]:>factor im e[exp+Pi I/2 //Simplify],
							factor_ Power[e,exp_] Complex[0,im_]/;!FreeQ[exp,_Complex]:>factor im e[exp+Pi I/2 //Simplify]}/.{
							factor_ Power[Complex[0,im1_],exp_] Complex[0,im2_]:>factor im1 im2 e[(exp+1)Pi I/2 //Simplify]};(* unify the (-1)^n or I^n to ~ e^(Pi I n) *)
*)		
(*--- simplify the generated rational function of D ---*)
	If[sgamma===True,
		result[[1,1]]qfact2[Simplify[Plus@@result[[2]]/.qfact2[factors_]:>factors] ](* avoid nesting qfact2 *)
	,
		result[[1,1]](Plus@@result[[2]])/.qfact2[factors_]:>factors
	]
]
]




(*----------------------------------------------------------------------------*)


End[]
