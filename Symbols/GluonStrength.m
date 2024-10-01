(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



GluonStrength::usage = 
	"GluonStrength[{mu,nu,a},b___] is the gluontensor tensor symbol \!\(\*SubscriptBox[SuperscriptBox[\(G\), \(a\)], \(\[Mu]\[Nu]\)]\)"
	
	
Begin["`Private`GluonStrength`"]	
(*Options[GluonStrength] = {}*)



(*GluonStrength[a:Except[_List|_LorentzIndex],LorentzIndex[b_,dim___],c___]:=GluonStrength[LorentzIndex[a,dim],LorentzIndex[b,dim],c]
GluonStrength[LorentzIndex[a_,dim___],b:Except[_List|_LorentzIndex],c___]:=GluonStrength[LorentzIndex[a,dim],LorentzIndex[b,dim],c]


GluonStrength[LorentzIndex[a_,dim1_:4],LorentzIndex[b_,dim2_:4],c___]:=Block[{dim},
dim=Times@@({dim1,dim2}/.D-4\[Rule]{0,1}/.{4\[Rule]{1,0},D\[Rule]{1,1}});
dim=dim/.{{1,1}\[Rule]D,{1,0}\[Rule]4,{0,1}\[Rule]D-4,{0,0}\[Rule]0};

If[dim===0,
	0
,
	GluonStrength[LorentzIndex[a,dim],LorentzIndex[b,dim],c]
]
]



(*------------------------------------------------------*)

GluonStrength[{a:Except[_List|_LorentzIndex],LorentzIndex[b_,dim___],sun_}c___]:=GluonStrength[{LorentzIndex[a,dim],LorentzIndex[b,dim],sun},c]
GluonStrength[{LorentzIndex[a_,dim___],b:Except[_List|_LorentzIndex],sun_},c___]:=GluonStrength[{LorentzIndex[a,dim],LorentzIndex[b,dim],sun},c]


GluonStrength[{LorentzIndex[a_,dim1_:4],LorentzIndex[b_,dim2_:4],sun_},c___]:=Block[{dim},
dim=Times@@({dim1,dim2}/.D-4\[Rule]{0,1}/.{4\[Rule]{1,0},D\[Rule]{1,1}});
dim=dim/.{{1,1}\[Rule]D,{1,0}\[Rule]4,{0,1}\[Rule]D-4,{0,0}\[Rule]0};

If[dim===0,
	0
,
	GluonStrength[{LorentzIndex[a,dim],LorentzIndex[b,dim],sun},c]
]
]*)


GluonStrength[a:Except[_List],b_,c___]:=Signature[{a,b}]GluonStrength[##,c]&@@Sort[{a,b}]/;!OrderedQ[{a,b}]

GluonStrength[{a_,b_,sun_},c___]:=Signature[{a,b}]GluonStrength[{##,sun},c]&@@Sort[{a,b}]/;!OrderedQ[{a,b}]

GluonStrength[{a_,b_},c___]:=GluonStrength[a,b,c]/;Head[a]===Head[b]


GluonStrength[{a_,a_,sun_},c___]=0
GluonStrength[{a_,a_},c___]=0

(* SUNSimplify don't know what GluonStrength is. *)
GluonStrength/:SUNSimplify[GluonStrength[expr_]]:=GluonStrength[expr]


(*-------------------------------------------------------------------------------------------*)
(*GluonStrength[lors__LorentzIndex]:=GluonStrength[##]&@@({lors}/.LorentzIndex[lo_,___]:>lo)*)


GluonStrength/:MakeBoxes[GluonStrength[mu:Except[_List],nu_,a___],TraditionalForm]:=Block[{dlist},
If[Length[{a}]===0,

	SubscriptBox["G",RowBox[{ToBoxes[mu,TraditionalForm],ToBoxes[nu,TraditionalForm]}]]
,
	dlist=SubscriptBox["D",ToBoxes[#,TraditionalForm]]&/@{a};

RowBox[{##,SubscriptBox["G",RowBox[{ToBoxes[mu,TraditionalForm],ToBoxes[nu,TraditionalForm]}]]}]&@@dlist
]

]




(*GluonStrength[{mu_LorentzIndex,nu_LorentzIndex,sun_},lors___LorentzIndex]:=GluonStrength[{mu/.LorentzIndex[lo_,___]:>lo,nu/.LorentzIndex[lo_,___]:>lo,sun},##]&@@({lors}/.LorentzIndex[lo_,___]:>lo)

*)
GluonStrength/:MakeBoxes[GluonStrength[{mu_,nu_,sun_},a___],TraditionalForm]:=Block[{dlist},
If[Length[{a}]===0,

	SubsuperscriptBox["G",RowBox[{ToBoxes[mu,TraditionalForm],ToBoxes[nu,TraditionalForm]}],ToBoxes[sun,TraditionalForm]]
,
	dlist=SubscriptBox["D",ToBoxes[#,TraditionalForm]]&/@{a};

SuperscriptBox[RowBox[{"(",##,SubscriptBox["G",RowBox[{ToBoxes[mu,TraditionalForm],ToBoxes[nu,TraditionalForm]}]],")"}],ToBoxes[sun,TraditionalForm]]&@@dlist
]

]


End[]
