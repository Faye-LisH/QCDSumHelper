(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



GluonStrength::usage = 
	"GluonStrength[{mu,nu,a},b___] is the gluontensor tensor symbol \!\(\*SubscriptBox[SuperscriptBox[\(G\), \(a\)], \(\[Mu]\[Nu]\)]\)"
	
	
Begin["`Private`GluonStrength`"]	
Options[GluonStrength] = {
	LeadingOnly->False
	}
(* allow to add a label in G^uv, so that d^uA^v-d^vA^u can be denoted as GluonStrength[n,v,LeadingOnly->True] *)



GluonStrength[a:Except[_List],b_,c___]:=Signature[{a,b}]GluonStrength[##,c]&@@Sort[{a,b}]/;!OrderedQ[{a,b}]
GluonStrength[{a_,b_,sun_},c___]:=Signature[{a,b}]GluonStrength[{##,sun},c]&@@Sort[{a,b}]/;!OrderedQ[{a,b}]

(*GluonStrength[{a_,b_},c___]:=GluonStrength[a,b,c]/;Head[a]===Head[b]*)


GluonStrength[{a_,a_,sun_},c___]=0
GluonStrength[{a_,a_},c___]=0

(* f/:g[f,OptioinsPattern[]]:= ... doesn't work, OptionsPattern doesn't work inside UpValue *)
GluonStrength[a:Except[_List],b_,c_List,OptionsPattern[]]:=If[OptionValue[LeadingOnly]===True,GluonStrength[a,b,c,True],GluonStrength[a,b,c,False]]
GluonStrength[a:Except[_List],b_,OptionsPattern[]]:=If[OptionValue[LeadingOnly]===True,GluonStrength[a,b,{},True],GluonStrength[a,b,{},False]]
GluonStrength[{a_,b_,sun_},c_List,OptionsPattern[]]:=If[OptionValue[LeadingOnly]===True,GluonStrength[{a,b,sun},c,True],GluonStrength[{a,b,sun},c,False]]
GluonStrength[{a_,b_,sun_},OptionsPattern[]]:=If[OptionValue[LeadingOnly]===True,GluonStrength[{a,b,sun},{},True],GluonStrength[{a,b,sun},{},False]]


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

GluonStrength/:MakeBoxes[GluonStrength[mu:Except[_List],nu_,de_List,leading_],TraditionalForm]:=Block[{guv,dlist},
If[leading===True,
	guv=UnderscriptBox["G","_"]
,
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
GluonStrength/:MakeBoxes[GluonStrength[{mu_,nu_,sun_},de_List,leading_],TraditionalForm]:=Block[{guv,dlist,covd=True},

If[leading===True,
	guv=UnderscriptBox["G","_"]
,
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
