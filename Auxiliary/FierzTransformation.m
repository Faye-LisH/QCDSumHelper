(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



FierzTransformation::usage = "FierzTransformation[{a_,gA_,b_},{c_,gB_,d_}] give the Fierz transformation of \!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(A\)]\)\!\(\*SubscriptBox[\(\[CapitalPsi]\), \(b\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(B\)]\)\!\(\*SubscriptBox[\(\[CapitalPsi]\), \(d\)]\), the output is a list of terms ~{{a,gA',d},{c,gB',b}};
FierzTransformation[{{a_,n_},gA_,{b_,m_}},{{c_,m_},gB_,{d_,n_}}] give the Fierz transformation of (\!\(\*SubsuperscriptBox[\(\[CapitalPsi]\), \(a\), \(n\)]\)\!\(\*SuperscriptBox[\()\), \(T\)]\)C \!\(\*SubscriptBox[\(\[CapitalGamma]\), \(A\)]\) \!\(\*SubsuperscriptBox[\(\[CapitalPsi]\), \(b\), \(m\)]\) \!\(\*SubsuperscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c\), \(m\)]\) \!\(\*SubscriptBox[\(\[CapitalGamma]\), \(B\)]\)C(\!\(\*SubsuperscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(d\), \(n\)]\)\!\(\*SuperscriptBox[\()\), \(T\)]\),
  the output is a list of fourquark molecule currents ~{{d,gA',a},{c,gB',b}}"

FierzTransformation::dimerr="Dimension inconsistent!"

FierzTransformation::err="Unknown expression!"


Begin["`Private`FierzTransformation`"]


Options[FierzTransformation]={
	Ga5toEps->"Auto"}(* replace GA[5] by definition -i/4! GA[a,b,c,d]Eps[a,b,c,d] *)




FierzTransformation[expr_/;!FreeQ[expr,FourquarkCurrent],ops:OptionsPattern[]]:=Block[{tmp,tmp1,null,aa,bb,current,psi=Identity,dot,dott},

tmp= expr/.fc_FourquarkCurrent:>((*--- a very long rule ---*)


null=fc/.FourquarkCurrent[ca_,cb_]:>(aa=current@@ca; bb=current@@cb; 1);


If[And@@(MatchQ[#,current[{_},{_}]|current[Except[_List],Except[_List]]|current[{_},_,{_}]|current[Except[_List],_,Except[_List]]]&/@{aa,bb}),
(* meson-meson current *)

(* unify the form of currents *)
	{aa,bb}={aa,bb}/.current[a1_,b1_]:>current[a1,1,b1]/.current[List[a1_],g_,List[b1_]]:>(psi=List;current[a1,g,b1]);(* use psi to label the style of Current *)

	tmp1=FierzTransformation[{aa[[1]],aa[[2]],aa[[3]]},{bb[[1]],bb[[2]],bb[[3]]},ops];
	
,

(* diquark-antidiquark current *)
 	

	If[MatchQ[aa,current[{_,_},{_,_}]|current[_,_,_,_]|current[{_,_},_,{_,_}]|current[_,_,_,_,_]]&&MatchQ[bb,current[{-1,_,_},{-1,_,_}]|current[-1,_,_,-1,_,_]|current[{-1,_,_},_,{-1,_,_}]|current[-1,_,_,_,-1,_,_]],
		
		aa=aa/.current[a1_List,b1_List]:>(psi=List;current[a1,1,b1])/.{current[a1_,b1_,c1_,d1_]:>current[{a1,b1},1,{c1,d1}], current[a1_List,g_,b1_List]:>(psi=List;current[a1,g,b1]), current[a1_,b1_,g_,c1_,d1_]:>current[{a1,b1},g,{c1,d1}]};
		bb=bb/.current[a1_List,b1_List]:>(psi=List;current[a1[[2;;3]],1,b1[[2;;3]]])/.{
					current[-1,a1_,b1_,-1,c1_,d1_]:>current[{a1,b1},1,{c1,d1}], current[{-1,a1_,b1_},g_,{-1,c1_,d1_}]:>(psi=List;current[{a1,b1},g,{c1,d1}]), current[-1,a1_,b1_,g_,-1,c1_,d1_]:>current[{a1,b1},g,{c1,d1}]};
		
		tmp1=FierzTransformation[List@@aa,List@@bb,ops];
		
	,

		Message[FierzTransformation::err];
		(*Hold[FierzTransformation[fc]]*)
		Abort[]
	]

	
];


(*--- recover the results to currents form ---*)
tmp1=(current[#[[1,1]],#[[1,2]],#[[1,3]]]current[#[[2,1]],#[[2,2]]//Expand,#[[2,3]]]&/@tmp1)//FCI;
tmp1=(SUNSimplify[Contract[#]]&/@tmp1)/.current[a1_,b1_,c1_]:>current[a1,b1//Expand,c1]/.current[a1_,b1_Plus,c1_]:>(current[a1,#,c1]&/@b1);
(* e.g. q^bar (g_1+ g_2 + ...) q -> q^bar g_1 q + q^bar g_2 q + ... *)

(* isolate non-communicate part and Eps  *)
tmp1=(tmp1//.current[a_,b1_ b2_,c_]/;And@@(FreeQ[b1,#]&/@{DiracGamma,DiracSigma,SUNT,Eps}):>b1 current[a,b2,c])//.current[a_,b_,c_]/;And@@(FreeQ[b,#]&/@{DiracGamma,DiracSigma,SUNT,Eps}):>b current[a,1,c];	


tmp1=FCI[tmp1//Expand]/.Dot->dot;

(* protect the term like GAD[a,b,c]LCD[u,a,b,c] *)
tmp1=tmp1/.gg_ ee_Eps/;(Or@@(!FreeQ[Cases[gg,LorentzIndex[ll_,___]:>ll,Infinity],#]&/@Cases[ee,LorentzIndex[ll_,___]:>ll,Infinity])):>dott[gg,ee];

tmp1=tmp1/.{current[aa1_,ds_DiracSigma,cc1_]:>current[aa1,dot[ds],cc1],current[aa1_,ds_DiracSigma ss_SUNT,cc1_]:>current[aa1,dot[ds,ss],cc1]}/.{current[aa1_,dot[do__] ss_SUNT,cc1_]:>current[aa1,dot[do,ss],cc1]}/.{
				current[aa1_,dot[ad___,DiracSigma[DiracGamma[LorentzIndex[lora_,di___],di___],DiracGamma[LorentzIndex[lorb_,di___],di___]],bd___],cc1_]current[aa2_,dot[aad___,DiracGamma[LorentzIndex[lora_,di___],di___],DiracGamma[LorentzIndex[lorb_,di___],di___],bbd___],cc2_]:>
					-I current[aa1,dot[ad,DiracSigma[DiracGamma[LorentzIndex[lora,di],di],DiracGamma[LorentzIndex[lorb,di],di]],bd],cc1]current[aa2,dot[aad,DiracSigma[DiracGamma[LorentzIndex[lora,di],di],DiracGamma[LorentzIndex[lorb,di],di]],bbd],cc2]
				,
				current[aa1_,dot[ad___,DiracSigma[DiracGamma[LorentzIndex[lora_,di___],di___],DiracGamma[LorentzIndex[lorb_,di___],di___]],bd___],cc1_]current[aa2_,dot[aad___,DiracGamma[LorentzIndex[lorb_,di___],di___],DiracGamma[LorentzIndex[lora_,di___],di___],bbd___],cc2_]:>
					I current[aa1,dot[ad,DiracSigma[DiracGamma[LorentzIndex[lora,di],di],DiracGamma[LorentzIndex[lorb,di],di]],bd],cc1]current[aa2,dot[aad,DiracSigma[DiracGamma[LorentzIndex[lora,di],di],DiracGamma[LorentzIndex[lorb,di],di]],bbd],cc2]
				}/.dot->Dot;(* e.g. (q^bar sigma^ab q)(q^bar g^a g^b q) -> -i (q^bar sigma^ab q)(q^bar sigma^ab q) *)
				
					
					
(* It's possible that -sigma^ab generated after use Qsimplify, act it twice *)

tmp1=QSimplify[QSimplify[Plus@@tmp1,Separate->current]//.current[a_,b1_ b2_,c_]/;And@@(FreeQ[b1,#]&/@{DiracGamma,DiracSigma,SUNT}):>b1 current[a,b2,c],Separate->current]//Expand//QNormal;

tmp1=tmp1/.dott->Dot;(* recover the dott to Dot *)

tmp1=(tmp1//FCI)//.{Pair[LorentzIndex[lor1_,dim1_:4],LorentzIndex[lor2_,dim1_:4]]current[aa_,bb_,cc_]/;!FreeQ[bb,LorentzIndex[lor1,___]]:>current[aa,bb/.LorentzIndex[lor1,dim2_:4]:>LorentzIndex[lor2,dimcontract[dim1,dim2]],cc],
			Pair[LorentzIndex[lor1_,dim1_:4],LorentzIndex[lor2_,dim1_:4]]current[aa_,bb_,cc_]/;!FreeQ[bb,LorentzIndex[lor2,___]]:>current[aa,bb/.LorentzIndex[lor2,dim2_:4]:>LorentzIndex[lor1,dimcontract[dim1,dim2]],cc]};


(* restore the style of Current *)
tmp1/.current[aa1_,bb1_,cc1_]current[aa2_,bb2_,cc2_]:>FourquarkCurrent[Current[psi[aa1],bb1,psi[cc1]],Current[psi[aa2],bb2,psi[cc2]]]
)


]




(* ::Code::Initialization::Plain:: *)
dimcontract[d1_,d2_]:=(Times@@({d1,d2}/.D-4->{0,1}/.{D->{1,1},4->{1,0}}))/.{{1,1}->D,{1,0}->4,{0,1}->D-4}
(* get the overlap of two dimensions *)


(* ::Code::Initialization::Plain:: *)
vsimp[expr_]:=Simplify[DiracSimplify[SUNSimplify[expr]]](* vertex simplify *)



(* ::Code::Initialization::Plain:: *)
FierzTransformation[{flavor1_,gammaA_,flavor2_},{flavor3_,gammaB_,flavor4_},OptionsPattern[]]:=Block[{tmp,ga,dim,dim1=D,dim2=D,lor1=$AL[Unique[]],lor2=$AL[Unique[]],cola=$AC[Unique[]],g5e=OptionValue[Ga5toEps]},

(* get the dimension *)
FCI[gammaA]/.DiracGamma[LorentzIndex[__],di_:4]:>(dim1=di;1);
FCI[gammaB]/.DiracGamma[LorentzIndex[__],di_:4]:>(dim2=di;1);


If[FreeQ[gammaA//FCI,DiracGamma],dim1=dim2];
If[FreeQ[gammaB//FCI,DiracGamma],dim2=dim1];

If[dim1=!=dim2,
	Message[FierzTransformation::dimerr];
	dim=D
,
	dim=dim1
];



(*---------------------------------------------*)

(* dim = 4 or D, depending on the input *)
	tmp={{{flavor1,1,flavor4},{flavor3,-vsimp[gammaB . gammaA]/(4CA),flavor2}}(*TR[SUNTrace[1]]= 4 CA *)
,
		{{flavor1,DiracGamma[LorentzIndex[lor1,dim],dim],flavor4},{flavor3,-vsimp[gammaB . DiracGamma[LorentzIndex[lor1,dim],dim] . gammaA]/(4 CA),flavor2}}(*TR[SUNTrace[ga[a].ga[b]]]= 4 CA g^ab *)
,
		{{flavor1,DiracGamma[5],flavor4},{flavor3,-vsimp[gammaB . DiracGamma[5] . gammaA]/(4CA),flavor2}}(*TR[SUNTrace[DiracGamma[5].DiracGamma[5]]]= 4 CA *)
,
		{{flavor1,DiracGamma[LorentzIndex[lor1,dim],dim] . DiracGamma[5],flavor4},{flavor3,-vsimp[gammaB . DiracGamma[LorentzIndex[lor1,dim],dim] . DiracGamma[5] . gammaA]/(-4 CA),flavor2}}(*TR[SUNTrace[ga[a].DiracGamma[5].ga[b].DiracGamma[5]]]= -4 CA g^ab *)
,
		{{flavor1,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],flavor4},{flavor3,-vsimp[gammaB . DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]] . gammaA]/(8 CA),flavor2}}(*TR[SUNTrace[DiracSigma[ga[a],ga[b]].DiracSigma[ga[c],ga[d]]]]= 8 CA (g^ac g^bd - g^ad g^bc)/2 *)

		,(*------------------------------------------------*)

		{{flavor1,SUNT[cola],flavor4},{flavor3,-vsimp[gammaB . SUNT[cola] . gammaA]/(2),flavor2}}(*TR[SUNTrace[SUNT[n].SUNT[m]]] = 2 delta^nm *)
,
		{{flavor1,DiracGamma[LorentzIndex[lor1,dim],dim] . SUNT[cola],flavor4},{flavor3,-vsimp[gammaB . DiracGamma[LorentzIndex[lor1,dim],dim] . SUNT[cola] . gammaA]/(2),flavor2}}(*TR[SUNTrace[ga[a].SUNT[n].ga[b].SUNT[m]]] = 2 g^ab delta^nm *)
,
		{{flavor1,DiracGamma[5] . SUNT[cola],flavor4},{flavor3,-vsimp[gammaB . DiracGamma[5] . SUNT[cola] . gammaA]/(2),flavor2}}(*TR[SUNTrace[DiracGamma[5].SUNT[cola].DiracGamma[5].SUNT[cola]]] = 2 g^ab delta^nm *)
,
		{{flavor1,DiracGamma[LorentzIndex[lor1,dim],dim] . DiracGamma[5] . SUNT[cola],flavor4},{flavor3,-vsimp[gammaB . DiracGamma[LorentzIndex[lor1,dim],dim] . DiracGamma[5] . SUNT[cola] . gammaA]/(-2),flavor2}}(*TR[SUNTrace[ga[lor1].DiracGamma[5].SUNT[cola].ga[lor1].DiracGamma[5].SUNT[cola]]] = -2 g^ab delta^nm *)
,
		{{flavor1,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]] . SUNT[cola],flavor4},{flavor3,-vsimp[gammaB . DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]] . SUNT[cola] . gammaA]/(4),flavor2}}(*TR[SUNTrace[DiracSigma[ga[lor1],ga[lor2]].SUNT[cola].DiracSigma[ga[lor1],ga[lor2].SUNT[cola]]]] = 4 (g^ac g^bd - g^ad g^bc)/2 delta^nm *)
		}/.{_,{_,0,_}}->0;





(* Contract dummy indices *)
(*/.{{{f1_,vgammaA_,f4_},{f3_,SUNDelta[cc1_,cc2_]vgammaB_,f2_}}/;!FreeQ[vgammaA,cc1]:>{{f1,vgammaA/.cc1->cc2,f4},{f3,vgammaB,f2}},
								{{f1_,vgammaA_,f4_},{f3_,SUNDelta[cc1_,cc2_]vgammaB_,f2_}}/;!FreeQ[vgammaA,cc2]:>{{f1,vgammaA/.cc2->cc1,f4},{f3,vgammaB,f2}}}//.
								{{{f1_,vgammaA_,f4_},{f3_,Pair[LorentzIndex[ll1_],LorentzIndex[ll2_]]vgammaB_,f2_}}/;!FreeQ[vgammaA,ll1]:>{{f1,vgammaA/.ll1->ll2,f4},{f3,vgammaB,f2}},
								{{f1_,vgammaA_,f4_},{f3_,Pair[LorentzIndex[ll1_],LorentzIndex[ll2_]]vgammaB_,f2_}}/;!FreeQ[vgammaA,ll2]:>{{f1,vgammaA/.ll2->ll1,f4},{f3,vgammaB,f2}}}*)
						
tmp=replacega5[tmp,g5e,dim];
		
tmp=DeleteCases[tmp,0]/.{CA-2CF->1/CA, 2CF-CA->-1/CA}								

]






(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------------------------------------------*)




(* ::Code::Initialization::Plain:: *)
(* for tetraquark current *)



(* ::Code::Initialization::Plain:: *)
FierzTransformation[{{flavor1_,col1_},gammaA_,{flavor2_,col2_}},{{flavor3_,col2_},gammaB_,{flavor4_,col1_}},OptionsPattern[]]:=Block[{tmp,ga,dim,dim1=D,dim2=D,lor1=$AL[Unique[]],lor2=$AL[Unique[]],cola=$AC[Unique[]],g5e=OptionValue[Ga5toEps]},

(* get the dimension *)
FCI[gammaA]/.DiracGamma[LorentzIndex[__],di_:4]:>(dim1=di;1);
FCI[gammaB]/.DiracGamma[LorentzIndex[__],di_:4]:>(dim2=di;1);


If[FreeQ[gammaA//FCI,DiracGamma],dim1=dim2];
If[FreeQ[gammaB//FCI,DiracGamma],dim2=dim1];


If[dim1=!=dim2,
	Message[FierzTransformation::dimerr];
	dim=D
,
	dim=dim1
];



(*---------------------------------------------*)


(* C^-1 = -C ; C gamma[u]^T C = gamma[u] ; C gamma5^T C = -gamma5 ;  C gamma[u]^T gamma5^T C = gamma5 gamma[u] ; C sigma^ab C = sigma^ab *)
tmp={
{{flavor4,1,flavor1},{flavor3,-vsimp[gammaB . gammaA]/(4),flavor2}}(* TR[1] = 4 *)
,
	{{flavor4,DiracGamma[LorentzIndex[lor1,dim],dim],flavor1},{flavor3,vsimp[gammaB . DiracGamma[LorentzIndex[lor1,dim],dim] . gammaA]/(4),flavor2}}(* TR[ga[a].ga[b]] = 4 g^ab *)
,
	{{flavor4,DiracGamma[5],flavor1},{flavor3,-vsimp[gammaB . DiracGamma[5] . gammaA]/(4),flavor2}}(* TR[DiracGamma[5].DiracGamma[5]] = 4 *)
,
	{{flavor4,DiracGamma[LorentzIndex[lor1,dim],dim] . DiracGamma[5],flavor1}
,{flavor3,vsimp[gammaB . DiracGamma[5] . DiracGamma[LorentzIndex[lor1,dim],dim] . gammaA]/(-4),flavor2}}(* TR[ga[a].DiracGamma[5].ga[b].DiracGamma[5]] = -4 g^ab *)
,
	{{flavor4,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],flavor1},{flavor3,vsimp[gammaB . DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]] . gammaA]/(8),flavor2}}(*TR[DiracSigma[ga[a],ga[b]].DiracSigma[ga[c],ga[d]]]= 8 (g^ac g^bd - g^ad g^bc)/2 *)
}/.{_,{_,0,_}}->0;





tmp=replacega5[tmp,g5e,dim];
								
tmp=DeleteCases[tmp,0]/.{ CA-2CF->1/CA, 2CF-CA->-1/CA}							

]




(* ::Code::Initialization::Plain:: *)
(* transpose the second, the color structure will be same as previous *)

FierzTransformation[{{flavor1_,col1_},gammaA_,{flavor2_,col2_}},{{flavor3_,col1_},gammaB_,{flavor4_,col2_}},OptionsPattern[]]:=Block[{tmp,ga,dim,dim1=D,dim2=D,lor1=$AL[Unique[]],lor2=$AL[Unique[]],cola=$AC[Unique[]],g5e=OptionValue[Ga5toEps]},

(* get the dimension *)
FCI[gammaA]/.DiracGamma[LorentzIndex[__],di_:4]:>(dim1=di;1);
FCI[gammaB]/.DiracGamma[LorentzIndex[__],di_:4]:>(dim2=di;1);


If[FreeQ[gammaA//FCI,DiracGamma],dim1=dim2];
If[FreeQ[gammaB//FCI,DiracGamma],dim2=dim1];

If[dim1=!=dim2,
	Message[FierzTransformation::dimerr];
	dim=D
,
	dim=dim1
];



(*---------------------------------------------*)
(* extra minus sign for communtate Fermions *)
(* C^T (g^B)^T g^T C g^A  *)
tmp={{{flavor3,1,flavor1},{flavor4,-vsimp[FCChargeConjugateTransposed[gammaB,Explicit->True] . gammaA]/(4),flavor2}}(*TR[1] = 4 *)
,
	{{flavor3,DiracGamma[LorentzIndex[lor1,dim],dim],flavor1},{flavor4,-vsimp[FCChargeConjugateTransposed[DiracGamma[LorentzIndex[lor1,dim],dim].gammaB,Explicit->True] . gammaA]/(4),flavor2}}(*TR[ga[a].ga[b]] = 4 g^ab *)
,
	{{flavor3,DiracGamma[5],flavor1},{flavor4,-vsimp[ FCChargeConjugateTransposed[DiracGamma[5] .gammaB ,Explicit->True]. gammaA]/(4),flavor2}}(*TR[DiracGamma[5].DiracGamma[5]] = 4 *)
,
	{{flavor3,DiracGamma[LorentzIndex[lor1,dim],dim] . DiracGamma[5],flavor1},{flavor4,-vsimp[FCChargeConjugateTransposed[DiracGamma[LorentzIndex[lor1,dim],dim] . DiracGamma[5] .gammaB,Explicit->True].  gammaA]/(-4),flavor2}}(*TR[ga[a].DiracGamma[5].ga[b].DiracGamma[5]] = -4 g^ab *)
,
	{{flavor3,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],flavor1},{flavor4,-vsimp[FCChargeConjugateTransposed[DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]] .gammaB ,Explicit->True].  gammaA]/(8),flavor2}}(*TR[DiracSigma[ga[a],ga[b]].DiracSigma[ga[c],ga[d]]]= 4 (g^ac g^bd - g^ad g^bc)*)
}/.{_,{_,0,_}}->0;



tmp=replacega5[tmp,g5e,dim];
								
tmp=DeleteCases[tmp,0]/.{CA-2CF->1/CA, 2CF-CA->-1/CA}		
]






(* ::Code::Initialization::Plain:: *)
replacega5[expr_,g5e_,dim_]:=Block[{tmp,dot,lore1,lore2,lore3,lore4},

tmp=expr//FCI;

(* recover the expanded DiracSigma[GAD[a,b]] *)
tmp=tmp/.Dot->dot/.dot[aa___,DiracGamma[LorentzIndex[lora_,di___],di___],DiracGamma[LorentzIndex[lorb_,di___],di___],bb___]-dot[aa___,DiracGamma[LorentzIndex[lorb_,di___],di___],DiracGamma[LorentzIndex[lora_,di___],di___],bb___]
								:>-2I dot[aa,DiracSigma[DiracGamma[LorentzIndex[lora,di],di],DiracGamma[LorentzIndex[lorb,di],di]],bb];


(* remove the GA[5] in DiracSigma[GAD[a,b]].GA[5] *)
(* the dott should be recover to Dot after all simplification are done; here it should not be touched *)
If[g5e===True||g5e==="All"||g5e==="Auto",

tmp=tmp//.
dot[aa___,DiracSigma[DiracGamma[LorentzIndex[lora_,di___],di___],DiracGamma[LorentzIndex[lorb_,di___],di___]],DiracGamma[5],bb___]:>(lore1=$AL[Unique[]];lore2=$AL[Unique[]];-1/2 dot[aa,DiracGamma[LorentzIndex[lore1,di],di].DiracGamma[LorentzIndex[lore2,di],di]Eps[LorentzIndex[lora,di],LorentzIndex[lorb,di],LorentzIndex[lore1,di],LorentzIndex[lore2,di]],bb]);

tmp=tmp//.
dot[aa___,DiracGamma[5],DiracSigma[DiracGamma[LorentzIndex[lora_,di___],di___],DiracGamma[LorentzIndex[lorb_,di___],di___]],bb___]:>(lore1=$AL[Unique[]];lore2=$AL[Unique[]];-1/2 dot[aa,DiracGamma[LorentzIndex[lore1,di],di].DiracGamma[LorentzIndex[lore2,di],di]Eps[LorentzIndex[lora,di],LorentzIndex[lorb,di],LorentzIndex[lore1,di],LorentzIndex[lore2,di]],bb]);
];


(* remove all GA[5] *)
If[g5e===True||g5e==="All",

tmp=tmp//.dot[aa___,DiracGamma[LorentzIndex[lora_,di___],di___],DiracGamma[5],bb___]:>(lore1=$AL[Unique[]];lore2=$AL[Unique[]];lore3=$AL[Unique[]];-I/6 dot[aa,DiracGamma[LorentzIndex[lore1,di],di].DiracGamma[LorentzIndex[lore2,di],di].DiracGamma[LorentzIndex[lore3,di],di]Eps[LorentzIndex[lora,di],LorentzIndex[lore1,di],LorentzIndex[lore2,di],LorentzIndex[lore3,di]],bb]);

tmp=tmp//.dot[aa___,DiracGamma[5],DiracGamma[LorentzIndex[lora_,di___],di___],bb___]:>(lore1=$AL[Unique[]];lore2=$AL[Unique[]];lore3=$AL[Unique[]];-I/6 dot[aa,DiracGamma[LorentzIndex[lore1,di],di].DiracGamma[LorentzIndex[lore2,di],di].DiracGamma[LorentzIndex[lore3,di],di]Eps[LorentzIndex[lora,di],LorentzIndex[lore1,di],LorentzIndex[lore2,di],LorentzIndex[lore3,di]],bb]);

tmp=tmp//.DiracGamma[5]->(lore1=$AL[Unique[]];lore2=$AL[Unique[]];lore3=$AL[Unique[]];lore4=$AL[Unique[]];DiracGamma[LorentzIndex[lore1,dim],dim].DiracGamma[LorentzIndex[lore2,dim],dim].DiracGamma[LorentzIndex[lore3,dim],dim].DiracGamma[LorentzIndex[lore4,dim],dim]Eps[LorentzIndex[lore1,dim],LorentzIndex[lore2,dim],LorentzIndex[lore3,dim],LorentzIndex[lore4,dim]](-I/24))

];

tmp/.dot->Dot
]


 End[]
 (* End the package *)
