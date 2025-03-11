(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



Fermion::usage="Fermion field symbol"
Fermion::derr="Unkonw derivative structure"
	
Begin["`Private`Fermion`"]	


 Fermion/:MakeBoxes[ Fermion[],TraditionalForm]:= ToBoxes["\[CapitalPsi]"]
(* Fermion/:MakeBoxes[ Fermion[x:Except[_List]],TraditionalForm]:= RowBox[{"\[CapitalPsi]","(",ToBoxes[x],")"}]
 Fermion/:MakeBoxes[ Fermion[f_,x_],TraditionalForm]:= RowBox[{ToBoxes[f],"(",ToBoxes[x],")"}]*)
 Fermion/:MakeBoxes[ Fermion[f:Except[_List]],TraditionalForm]:= ToBoxes[f]
 Fermion/:MakeBoxes[ Fermion[{f_}],TraditionalForm]:= SubscriptBox["\[CapitalPsi]",ToBoxes[f]]
(* Fermion/:MakeBoxes[ Fermion[{f_},x_],TraditionalForm]:= RowBox[{SubscriptBox["\[CapitalPsi]",ToBoxes[f]],"(",ToBoxes[x],")"}]*)

Fermion/:MakeBoxes[ Fermion[flavor_,dirac_,color_,ds_List:{}],TraditionalForm]:=Block[{field,fla,tmp,dlist={},covd=True},
If[MatchQ[ds,{___,{___,_List,___},___}],
	Message[ Fermion::derr];
,

(*cdlist=ds//.{ad0___,d_List,ad1___}:>(dlist=Join[dlist,d];{ad0,ad1});(* seperate the derivaives *)*)

	If[MatchQ[flavor,{_}],
		field=SubscriptBox["\[CapitalPsi]",ToBoxes[flavor[[1]]]]
	,
		field=ToBoxes[flavor]
	];


	If[Length[ds]>0||MatchQ[flavor,{_}],
		If[And@@(MatchQ[#,_List]&/@ds),covd=False];(* if no covariant derivative involved and only partial derivative is involved *)

		(* derivative symbols *)
		tmp=ds//.{{d0___,d1:Except[_List]}:>(dlist=Prepend[dlist,SuperscriptBox["\[Del]",ToBoxes[d1]]];{d0}),{d0___,{dd0__,d1_}}:>(dlist=Prepend[dlist,SuperscriptBox["\[PartialD]",ToBoxes[d1]]];{d0,{dd0}}),{d0___,{d1_}}:>(dlist=Prepend[dlist,SuperscriptBox["\[PartialD]",ToBoxes[d1]]];{d0})};
			If[covd,
				RowBox[{"(",##,SubsuperscriptBox[")",ToBoxes[dirac],ToBoxes[color]]}]&@@Append[dlist,field]
			,
				If[MatchQ[flavor,{_}],
					RowBox[Join[dlist,{"(",field,SubsuperscriptBox[")",ToBoxes[dirac],ToBoxes[color]]}]](* avoid the flavor and dirac index appears as same subscripts *)
				,
					RowBox[Join[dlist,{SubsuperscriptBox[field,ToBoxes[dirac],ToBoxes[color]]}]]
				]
			]

	,
		SubsuperscriptBox[field,ToBoxes[dirac],ToBoxes[color]]

	]
]
]


End[]
