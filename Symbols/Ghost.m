(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



Ghost::usage="Faddeev\[Dash]Popov ghost field symbol"
Ghost::derr="Unkonw derivative structure"
	
Begin["`Private`Ghost`"]	


Ghost/:MakeBoxes[ Ghost[],TraditionalForm]:=SubscriptBox[ToBoxes["c"],"F"]

Ghost/:MakeBoxes[ Ghost[color:Except[_List],ds_List:{}],TraditionalForm]:=Block[{tmp,dlist={},covd=True},
If[MatchQ[ds,{___,{___,_List,___},___}],
    Message[ Ghost::derr];
,
    If[Length[ds]>0,
        If[And@@(MatchQ[#,_List]&/@ds),covd=False];(* if no covariant derivative involved and only partial derivative is involved *)
        (* derivative symbols *)
        tmp=ds//.{{d0___,d1:Except[_List]}:>(dlist=Prepend[dlist,SuperscriptBox["D",ToBoxes[d1]]];{d0}),{d0___,{dd0__,d1_}}:>(dlist=Prepend[dlist,SuperscriptBox["\[PartialD]",ToBoxes[d1]]];{d0,{dd0}}),{d0___,{d1_}}:>(dlist=Prepend[dlist,SuperscriptBox["\[PartialD]",ToBoxes[d1]]];{d0})};
        If[covd,
            RowBox[{"(",##,SuperscriptBox[")",ToBoxes[color]]}]&@@Append[dlist,SubscriptBox[ToBoxes["c"],"F"]]
        ,
            RowBox[Join[dlist,{SubsuperscriptBox[ToBoxes["c"],"F",ToBoxes[color]]}]]
        ]


    ,
        SuperscriptBox[SubscriptBox[ToBoxes["c"],"F"],ToBoxes[color]]
    ]
]
]


End[]
