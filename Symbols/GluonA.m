(* ::Package:: *)

(* Wolfram Language package *)
(* Author: ShungHong Li *)



GluonA::usage="GluonA field symbol"
GluonA::derr="Unkonw derivative structure"

Begin["`Private`GluonA`"]	


GluonA/:MakeBoxes[GluonA[],TraditionalForm]:=ToBoxes["A"]

GluonA/:MakeBoxes[GluonA[lor:Except[_List],color:Except[_List],ds_List:{}],TraditionalForm]:=Block[{tmp,dlist={},covd=True},
If[MatchQ[ds,{___,{___,_List,___},___}],
    Message[GluonA::derr];
,
    If[Length[ds]>0,
        If[And@@(MatchQ[#,_List]&/@ds),covd=False];(* if no covariant derivative involved and only partial derivative is involved *)
        (* derivative symbols *)
        tmp=ds//.{{d0___,d1:Except[_List]}:>(dlist=Prepend[dlist,SuperscriptBox["D",ToBoxes[d1]]];{d0}),{d0___,{dd0__,d1_}}:>(dlist=Prepend[dlist,SuperscriptBox["\[PartialD]",ToBoxes[d1]]];{d0,{dd0}}),{d0___,{d1_}}:>(dlist=Prepend[dlist,SuperscriptBox["\[PartialD]",ToBoxes[d1]]];{d0})};

        If[covd,
            RowBox[{"(",##,SuperscriptBox[")",ToBoxes[color]]}]&@@Append[dlist,SubscriptBox["A",ToBoxes[lor]]]
        ,
            RowBox[Join[dlist,{SubsuperscriptBox["A",ToBoxes[lor],ToBoxes[color]]}]]
        ]

    ,
        SubsuperscriptBox["A",ToBoxes[lor],ToBoxes[color]]
    ]
]
]


End[]
