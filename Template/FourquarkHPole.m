(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)




FourquarkHPole::usage = "FourquarkHPole[{v1_,{f1_,f2_}},{v2_,{f3_,f4_}}] give the \[Epsilon]-pole at 1-loop level (at O(g)) about the Green function \[LeftAngleBracket]\!\(\*SubscriptBox[\(J\), \(4  q\)]\) \!\(\*SubscriptBox[OverscriptBox[\(\[Psi]\), \(_\)], \(a\)]\)\!\(\*SubsuperscriptBox[\(A\), \(\[Mu]\), \(n\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(b\)]\)\[RightAngleBracket] with \!\(\*SubscriptBox[\(J\), \(4  q\)]\)= \!\(\*SubscriptBox[OverscriptBox[\(\[Psi]\), \(_\)], \(1\)]\)\!\(\*SubscriptBox[\(v\), \(1\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(2\)]\) \!\(\*SubscriptBox[OverscriptBox[\(\[Psi]\), \(_\)], \(3\)]\)\!\(\*SubscriptBox[\(v\), \(2\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(d\)]\). } "

FourquarkHPole::dimerr="Dimension inconsistent!"
FourquarkHPole::dimwarn="Can't get the dimension, it will be set to D by default."


Begin["`Private`FourquarkHPole`"]

Options[FourquarkHPole]={
	HoldFlavor->True,
	Massless->True,
	ShowasTable->True,
	KeepGamma->False,
	Style2->True,
	LessDummy->False,
	SetDimension->D}



FourquarkHPole[factor_ cur_FourquarkCurrent,OptionsPattern[]]:=Block[{tb=OptionValue[ShowasTable],mless=OptionValue[Massless],hf=OptionValue[HoldFlavor],sdim=OptionValue[SetDimension]},
If[ToLowerCase[ToString[tb]]=="true",
	{#[[1]]factor,#[[2]]}&/@FourquarkHPole[cur,ShowasTable->True,Massless->mless,HoldFlavor->hf,SetDimension->sdim]
,
	factor FourquarkHPole[cur,ShowasTable->False,Massless->mless,HoldFlavor->hf,SetDimension->sdim]
]
]


FourquarkHPole[cur_FourquarkCurrent,OptionsPattern[]]:=Module[{q,mless=OptionValue[Massless],mass,v1,v2,pole,pole1,polem1,pole2,polem2,pole3,polem3,pole4,polem4,f1,f2,f3,f4,
cola,colb,colc,cold,cdelta,lora=$AL[Unique[]],lorb=$AL[Unique[]],sun=$AC[Unique[]],astb=OptionValue[ShowasTable],hf=OptionValue[HoldFlavor],sdim=OptionValue[SetDimension]},


(* get the gamma matices *)
v1=If[Length[cur[[1]]]==2,1,cur[[1,2]] ];
v2=If[Length[cur[[2]]]==2,1,cur[[2,2]] ];

(* get the flavors *)
{f1,f2,f3,f4}={cur[[1,1]],cur[[1,-1]],cur[[2,1]],cur[[2,-1]]};


pole=FourquarkHPole[{v1,{f1,f2}},{v2,{f3,f4}},Massless->mless,ShowasTable->True,HoldFlavor->hf,SetDimension->sdim];


If[ToLowerCase[ToString[astb]]=="true"||pole===0,
	pole
,
	Plus@@(#[[1]]HybirdCurrent[#[[2,1]],(List@@(#[[2,2]]/(#[[2,2]]/._GluonStrength->1)))/.{{a_,b_,_},c___}:>{a,b,c},#[[2,2]]/.{_GluonStrength->1,_SUNT->1},#[[2,3]]]&/@pole)
]


]/;Length[cur[[2]]]<=3




FourquarkHPole[{va_,{{f1_},{f2_}}},{vb_,{{f3_},{f4_}}},ops:OptionsPattern[]]:=FourquarkHPole[{va,{f1,f2}},{vb,{f3,f4}},ops]

FourquarkHPole[{va_,{f1_,f2_}},{vb_,{f3_,f4_}},OptionsPattern[]]:=Module[{tmp,lora=$AL[Unique[]],lorb=$AL[Unique[]],cola=$AC[Unique[]],null1,current,mless=OptionValue[Massless],hf=OptionValue[HoldFlavor],tr,v1,v2,dim,sdim=OptionValue[SetDimension]},

(*--- get the indices, dimension ---*)
{v1,v2}=ChangeDimension[#,4]&/@{va,vb};
dim=Cases[{va,vb}//FCI,LorentzIndex[_,dd_:4]:>dd,Infinity];
dim=DeleteDuplicates[dim];

If[sdim===D,
	dim=D
,
	If[ToLowerCase[sdim]==="asis"&&Length[dim]>1,
		Message[FourquarkHPole::dimerr];
		Abort[]
	,
		If[Length[dim]==0,(* if the vertex is gamma[5] or 1, the dimension is unclear *)
			dim=D;
			Message[FourquarkHPole::dimwarn]
		,
			dim=dim[[1]]
		]
	]
];

(*--- diagrams ---*)
tmp={
	{-gStrong FlavorDelta[f2,f3]/(Epsilon 48 Pi^2),{f1,(v1 . GA[lorb] . SUNT[cola] . v2//DiracSimplify//SUNSimplify) GluonStrength[{lora,lorb,cola},lora],f4}},

	{-gStrong FlavorDelta[f1,f4]/(Epsilon 48 Pi^2),{f3,(v2 . GA[lorb] . SUNT[cola] . v1//DiracSimplify//SUNSimplify) GluonStrength[{lora,lorb,cola},lora],f2}},

	{gStrong FlavorDelta[f1,f2]/(Epsilon 48 Pi^2),{f3,(tr[SUNTrace[v1 . GA[lorb] . SUNT[cola]]]v2//Contract//DiracSimplify//SUNSimplify) GluonStrength[{lora,lorb,cola},lora],f4}},

	{gStrong FlavorDelta[f3,f4]/(Epsilon 48 Pi^2),{f1,(tr[SUNTrace[v2 . GA[lorb] . SUNT[cola]]]v1//Contract//DiracSimplify//SUNSimplify) GluonStrength[{lora,lorb,cola},lora],f2}},


(**--- non-zero quark mass ---**)
	{quarkMass[f2]FlavorDelta[f2,f3]gStrong/(Epsilon 32 Pi^2),{f1,(v1 . DiracSigma[GA[lora,lorb]] . SUNT[cola] . v2//DiracSimplify//SUNSimplify) GluonStrength[{lora,lorb,cola}],f4}},

	{quarkMass[f1]FlavorDelta[f1,f4]gStrong/(Epsilon 32 Pi^2),{f3,(v2 . DiracSigma[GA[lora,lorb]] . SUNT[cola] . v1//DiracSimplify//SUNSimplify) GluonStrength[{lora,lorb,cola}],f2}},

	{-quarkMass[f1]FlavorDelta[f1,f2]gStrong/(Epsilon 32 Pi^2),{f3,(tr[SUNTrace[v1 . DiracSigma[GA[lora,lorb]] . SUNT[cola]]]v2//Contract//DiracSimplify//SUNSimplify) GluonStrength[{lora,lorb,cola}],f4}},

	{-quarkMass[f3]FlavorDelta[f3,f4]gStrong/(Epsilon 32 Pi^2),{f1,(tr[SUNTrace[v2 . DiracSigma[GA[lora,lorb]] . SUNT[cola]]]v1//Contract//DiracSimplify//SUNSimplify) GluonStrength[{lora,lorb,cola}],f2}}
};


 
(*--- simplify above---*)


If[OptionValue[KeepGamma]===True,
	tmp=tmp/.tr[aa_]:>DiracTrace[aa,DiracTraceEvaluate->False]
,
	tmp=tmp/.tr->TR
];

(*-----------------------------------------*)

If[ToLowerCase[ToString[mless]]==="true",
	tmp=tmp/.quarkMass[_]->0
];


tmp=tmp/.FlavorDelta[a_,b_,___]:>FlavorDelta[a,b,HoldFlavor->hf];


(*--- Contract lorentzindices ---*)
tmp=(tmp//FCI)//.{Pair[LorentzIndex[lo1_,di___],LorentzIndex[lo2_,di___]]GluonStrength[{aa_,lo1_,bb_},other___]:>GluonStrength[{aa,lo2,bb},other],
				SUNDelta[SUNIndex[coa_],SUNIndex[cob_]]GluonStrength[{aa_,bb_,coa_},other___]:>GluonStrength[{aa,bb,cob},other]};
				
tmp=DeleteCases[Replace[tmp,{0,{_,_,_}}|{_,{_,0,_}}->0,{1}],0,{1}];


tmp=Expand[Plus@@(#[[1]] current[#[[2,1]],#[[2,2]],#[[2,3]]]&/@tmp)];

tmp=tmp/.current[aa_,bb_,cc_]/;!FreeQ[bb,Plus]:>current[aa,Expand[bb],cc];

tmp=Expand[tmp/.current[aa_,bb_Plus,cc_]:>(current[aa,#,cc]&/@bb)]//.current[aa_,bb_ cc_,dd_]/;FreeQ[bb,DiracGamma|DiracSigma|SUNT|SUNTF|GluonStrength]:>bb current[aa,cc,dd];

(*
If[OptionValue[KeepGamma]=!=True,
(* dicard D-4 terms *)
	tmp=tmp/.current[aa_,cc_,dd_]/;!FreeQ[cc,LorentzIndex[_,D-4]]->0;
	tmp=tmp/.current[aa_,cc_,dd_]/;!FreeQ[cc,LorentzIndex[_,D-4]]->0/.current[aa_,cc_,dd_]:>current[aa,ChangeDimension[cc,D],dd]
];
*)

(*--- simplify dummy indices ---*)
tmp=tmp//.{Pair[LorentzIndex[lor1_],LorentzIndex[lor2_]]current[aa_,bb_,cc_]/;!FreeQ[bb,lor2]:>current[aa,bb/.lor2->lor1,cc],
		Pair[LorentzIndex[lor1_],LorentzIndex[lor2_]]current[aa_,bb_,cc_]/;!FreeQ[bb,lor1]:>current[aa,bb/.lor1->lor2,cc]};


If[OptionValue[LessDummy]===True,
	tmp=tmp/.GluonStrength[{la_,lb_,sun_},lc___]:>(GluonStrength[{LorentzIndex[la],LorentzIndex[lb],SUNIndex[sun]},LorentzIndex[##]]&@@{lc});
	tmp=tmp//QSimplify;
	tmp=tmp/.GluonStrength[{la_,lb_,sun_},lc___]:>(GluonStrength[{la/.LorentzIndex[lo_,___]:>lo,lb/.LorentzIndex[lo_,___]:>lo,sun/.SUNIndex[su_]:>su},##]&@@({lc}/.LorentzIndex[lo_,___]:>lo));

	tmp=tmp/.current[aa_,bb_ cc_,dd_]/;FreeQ[bb,DiracGamma|DiracSigma|SUNT|SUNTF|GluonStrength]:>bb current[aa,cc,dd]
(* the sign may changes, extract it outside of the current *)
];


(*--- collect the same currents ---*)
tmp=Collect[tmp,_current]//QNormal;

If[tmp=!=0,
	FCI[If[OptionValue[Style2]===True,
		(#/.aa_ bb_current:>{aa//Simplify,{bb[[1]],bb[[2]]/(bb[[2]]/._SUNT->1), bb[[2]]/(bb[[2]]/._GluonStrength->1),bb[[2]]/.{_SUNT->1, _GluonStrength->1}, bb[[3]]}}&/@If[Head[tmp]===Plus,List@@tmp,{tmp}])
	,
		(#/.aa_ bb_current:>{aa//Simplify,{bb[[1]],bb[[2]],bb[[3]]}}&/@If[Head[tmp]===Plus,List@@tmp,{tmp}])

	]]/.DiracGamma[LorentzIndex[lo_,___],___]:>DiracGamma[LorentzIndex[lo,dim],dim]/.{vv1_,{fa_,t_,g_,-1 ,fb_}}:>{-vv1,{fa,t,g,1,fb}}
,
	0
]
]



(*----------------------------------------------------*)


(* to do change dimension D *)
FourquarkHPole[cur_FourquarkCurrent,OptionsPattern[]]:=Module[{q,mless=OptionValue[Massless],mass,v1,v2,pole,pole1,polem1,pole2,polem2,pole3,polem3,pole4,polem4,f1,f2,f3,f4,
cola,colb,colc,cold,cdelta,lora=$AL[Unique[]],lorb=$AL[Unique[]],sun=$AC[Unique[]],astb=OptionValue[ShowasTable],hf=OptionValue[HoldFlavor],dim,sdim=OptionValue[SetDimension]},


(* get the color indices *)
{cola,colb,colc,cold}={cur[[1,2]],cur[[1,-1]],cur[[2,3]],cur[[2,-1]]};


(* get the gamma matices *)
v1=If[Length[cur[[1]]]==4,1,cur[[1,3]] ];
v2=If[Length[cur[[2]]]==6,1,cur[[2,4]] ];


(*--- get the indices, dimension ---*)
dim=Cases[{v1,v2}//FCI,LorentzIndex[_,dd_:4]:>dd,Infinity];
dim=DeleteDuplicates[dim];


If[sdim===D,
	dim=D
,
	If[ToLowerCase[sdim]==="asis"&&Length[dim]>1,
		Message[FourquarkHPole::dimerr];
		Abort[]
	,
		If[Length[dim]==0,(* if the vertex is gamma[5] or 1, the dimension is unclear *)
			dim=D;
			Message[FourquarkHPole::dimwarn]
		,
			dim=dim[[1]]
		]
	]
];

{v1,v2}=ChangeDimension[#,4]&/@{v1,v2};



(* get the flavors *)
{f1,f2,f3,f4}={cur[[1,1]],cur[[1,-2]],cur[[2,2]],cur[[2,-2]]};


(**)
pole1={-FlavorDelta[f2,f3]gStrong/(Epsilon 48 Pi^2),{f4,cold,SUNTF[sun,colb,colc],GluonStrength[{lora,lorb,sun},lora],FCChargeConjugateTransposed[v1 . GA[lorb] . v2,Explicit->True]//DiracSimplify//SUNSimplify,f1,cola}};

(**)
pole2={FlavorDelta[f1,f4]gStrong/(Epsilon 48Pi^2),{f3,colc,SUNTF[sun,cola,cold],GluonStrength[{lora,lorb,sun},lora],v2 . GA[lorb] . v1//DiracSimplify//SUNSimplify,f2,colb}};

(**)
pole3={-FlavorDelta[f1,f3]gStrong/(Epsilon 48Pi^2),{f4,cold,SUNTF[sun,cola,colc],GluonStrength[{lora,lorb,sun},lora],FCChargeConjugateTransposed[GA[lorb] . v2,Explicit->True] . v1//DiracSimplify//SUNSimplify,f2,colb}};

(**)
pole4={-FlavorDelta[f2,f4]gStrong/(Epsilon 48Pi^2),{f3,colc,SUNTF[sun,colb,cold],GluonStrength[{lora,lorb,sun},lora],v2 . FCChargeConjugateTransposed[v1 . GA[lorb],Explicit->True]//DiracSimplify//SUNSimplify,f1,cola}};


pole=If[ToLowerCase[ToString[mless]]=="true",

{pole1,pole2,pole3,pole4}

,
polem1={-quarkMass[f2]FlavorDelta[f2,f3]gStrong/(Epsilon 32Pi^2),{f4,cold,SUNTF[sun,colb,colc],GluonStrength[{lora,lorb,sun}],FCChargeConjugateTransposed[v1 . DiracSigma[GA[lora,lorb]] . v2,Explicit->True]//DiracSimplify//SUNSimplify,f1,cola}};

polem2={-quarkMass[f1]FlavorDelta[f1,f4]gStrong/(Epsilon 32 Pi^2),{f3,colc,SUNTF[sun,cola,cold],GluonStrength[{lora,lorb,sun}],v2 . DiracSigma[GA[lora,lorb]] . v1//DiracSimplify//SUNSimplify,f2,colb}};

polem3={quarkMass[f1]FlavorDelta[f1,f3]gStrong/(Epsilon 32 Pi^2),{f4,cold,SUNTF[sun,cola,colc],GluonStrength[{lora,lorb,sun}],FCChargeConjugateTransposed[DiracSigma[GA[lora,lorb]] . v2,Explicit->True] . v1//DiracSimplify//SUNSimplify,f2,colb}};

polem4={quarkMass[f2]FlavorDelta[f2,f4]gStrong/(Epsilon 32Pi^2),{f3,colc,SUNTF[sun,colb,cold],GluonStrength[{lora,lorb,sun}],v2 . FCChargeConjugateTransposed[v1 . DiracSigma[GA[lora,lorb]],Explicit->True]//DiracSimplify//SUNSimplify,f1,cola}};


{pole1,pole2,pole3,pole4,polem1,polem2,polem3,polem4}
];


(*-------------------------------------------*)

pole=pole/.{SUNTF[_,cc_,cc_]:>0,FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->hf]};


pole=DeleteCases[pole,{_,{_,_,_,_,0,_,_}}|{_,{_,_,0,_,_,_,_}}|{0,_}]/.
	{po_,{fa_,ca_,-sut_SUNTF,gg_,gv_,fb_,cb_}}:>{-po,{fa,ca,sut,gg,gv,fb,cb}}/.
	{po_,{fa_,ca_,sut_,- gg_GluonStrength,gv_,fb_,cb_}}:>{-po,{fa,ca,sut,gg,gv,fb,cb}}/.
	{po_,{fa_,ca_,sut_,gg_GluonStrength,fc_ gv_,fb_,cb_}}/;(And@@(FreeQ[fc,#]&/@{DiracGamma,DiracSigma,SUNT})):>{ fc po,{fa,ca,sut,gg,gv,fb,cb}};



FCI[If[ToLowerCase[ToString[astb]]=="true",
	pole
,
	Plus@@(#[[1]]HybirdCurrent[#[[2,1]],cola,List@@#[[2,4]],#[[2,5]],#[[2,6]],colb]&/@pole)
]]/.DiracGamma[LorentzIndex[lo_,___],___]:>DiracGamma[LorentzIndex[lo,dim],dim]

]/;Length[cur[[2]]]>=6


End[]
(*EndPackage[]*)
