(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



GluonStrength::usage = 
	"GluonStrength[{mu,nu,a},b___] is the gluontensor tensor symbol \!\(\*SubscriptBox[SuperscriptBox[\(G\), \(a\)], \(\[Mu]\[Nu]\)]\)"
	
	
Begin["`Private`GluonStrength`"]	
Options[GluonStrength] = {
	Order->"All"
	}
(* allow to add a label in G^uv, so that d^uA^v-d^vA^u can be denoted as GluonStrength[n,v,Order->0], and f^nabc A^buA^cv can be denoted as GluonStrength[n,v,Order->1] *)



GluonStrength[a:Except[_List],b_,c___]:=Signature[{a,b}]GluonStrength[##,c]&@@Sort[{a,b}]/;!OrderedQ[{a,b}]
(* This definition cause the pattern match GluonStrength[lis_,__] failed, unless without pattern name, like GluonStrength[_,__]; to match GluonStrength[lis_List,__], write it as GluonStrength[{a_,b_,c_},__] *)


GluonStrength[{a_,b_,sun_},c___]:=Signature[{a,b}]GluonStrength[{##,sun},c]&@@Sort[{a,b}]/;!OrderedQ[{a,b}]

(*GluonStrength[{a_,b_},c___]:=GluonStrength[a,b,c]/;Head[a]===Head[b]*)


GluonStrength[{a_,a_,sun_},c___]=0
GluonStrength[{a_,a_},c___]=0

(* f/:g[f,OptioinsPattern[]]:= ... doesn't work, OptionsPattern doesn't work inside UpValue *)
GluonStrength[a:Except[_List],b_,c_List,OptionsPattern[]]:=If[OptionValue[Order]==="All",GluonStrength[a,b,c,False],GluonStrength[a,b,c,OptionValue[Order]]]
GluonStrength[a:Except[_List],b_,OptionsPattern[]]:=If[OptionValue[Order]==="All",GluonStrength[a,b,{},False],GluonStrength[a,b,{},OptionValue[Order]]]
GluonStrength[{a_,b_,sun_},c_List,OptionsPattern[]]:=If[OptionValue[Order]==="All",GluonStrength[{a,b,sun},c,False],GluonStrength[{a,b,sun},c,OptionValue[Order]]]
GluonStrength[{a_,b_,sun_},OptionsPattern[]]:=If[OptionValue[Order]==="All",GluonStrength[{a,b,sun},{},False],GluonStrength[{a,b,sun},{},OptionValue[Order]]]


(* SUNSimplify don't know what GluonStrength is. *)
GluonStrength/:SUNSimplify[GluonStrength[expr__]]:=GluonStrength[expr]


(*-------------------------------------------------------------------------------------------*)
(*GluonStrength[lors__LorentzIndex]:=GluonStrength[##]&@@({lors}/.LorentzIndex[lo_,___]:>lo)*)


GluonStrength/:MakeBoxes[GluonStrength[OptionsPattern[]],TraditionalForm]:=ToBoxes["G"]


(*GluonStrength/:MakeBoxes[GluonStrength[mu:Except[_List],nu_,de_List,Leading_],TraditionalForm]:=Block[{dlist,covd=True},

If[Length[de]===0,

	SubscriptBox[UnderscriptBox["G","_"],RowBox[{ToBoxes[mu,TraditionalForm],ToBoxes[nu,TraditionalForm]}]]
,
	If[And@@(MatchQ[#,_List]&/@de),covd=False];(* if no covariant derivative involved and only partial derivative is involved *)
	dlist=If[MatchQ[#,_List],SubscriptBox["\[PartialD]",ToBoxes[#[[1]],TraditionalForm]],SubscriptBox["D",ToBoxes[#,TraditionalForm]]]&/@de;

	RowBox[{##,SubscriptBox[UnderscriptBox["G","_"],RowBox[{ToBoxes[mu,TraditionalForm],ToBoxes[nu,TraditionalForm]}]]}]&@@dlist
]

]*)

GluonStrength/:MakeBoxes[GluonStrength[mu:Except[_List],nu_,de_List,order_],TraditionalForm]:=Block[{guv,dlist},
Which[order===0,
	guv=UnderscriptBox["G","_"]
,
	order===1,
	guv=UnderscriptBox[UnderscriptBox["G","_"],"_"]
,
	True,
	guv="G"
];

If[Length[de]===0,

	SubscriptBox[guv,RowBox[{ToBoxes[mu,TraditionalForm],ToBoxes[nu,TraditionalForm]}]]
,
	dlist=If[MatchQ[#,_List],SubscriptBox["\[PartialD]",ToBoxes[#[[1]],TraditionalForm]],SubscriptBox["D",ToBoxes[#,TraditionalForm]]]&/@de;

	RowBox[{##,SubscriptBox[guv,RowBox[{ToBoxes[mu,TraditionalForm],ToBoxes[nu,TraditionalForm]}]]}]&@@dlist
]

]



(*GluonStrength[{mu_LorentzIndex,nu_LorentzIndex,sun_},lors___LorentzIndex]:=GluonStrength[{mu/.LorentzIndex[lo_,___]:>lo,nu/.LorentzIndex[lo_,___]:>lo,sun},##]&@@({lors}/.LorentzIndex[lo_,___]:>lo)

*)
GluonStrength/:MakeBoxes[GluonStrength[{mu_,nu_,sun_},de_List,order_],TraditionalForm]:=Block[{guv,dlist,covd=True},

Which[order===0,
	guv=UnderscriptBox["G","_"]
,
	order===1,
	guv=UnderscriptBox[UnderscriptBox["G","_"],"_"]
, 
	True,
	guv="G"
];

If[Length[de]===0,

	SubsuperscriptBox[guv,RowBox[{ToBoxes[mu,TraditionalForm],ToBoxes[nu,TraditionalForm]}],ToBoxes[sun,TraditionalForm]]
,
	If[And@@(MatchQ[#,_List]&/@de),covd=False];(* if no covariant derivative involved and only partial derivative is involved *)
	dlist=If[MatchQ[#,_List],SubscriptBox["\[PartialD]",ToBoxes[#[[1]],TraditionalForm]],SubscriptBox["D",ToBoxes[#,TraditionalForm]]]&/@de;
	
	If[covd,
		SuperscriptBox[RowBox[{"(",##,SubscriptBox[guv,RowBox[{ToBoxes[mu,TraditionalForm],ToBoxes[nu,TraditionalForm]}]],")"}],ToBoxes[sun,TraditionalForm]]&@@dlist
	,
		RowBox[Join[dlist,{SubsuperscriptBox[guv,RowBox[{ToBoxes[mu,TraditionalForm],ToBoxes[nu,TraditionalForm]}],ToBoxes[sun,TraditionalForm]]}]]
	]
]

]


End[]
