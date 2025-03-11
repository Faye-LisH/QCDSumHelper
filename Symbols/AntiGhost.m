(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



AntiGhost::usage="Anti-Faddeev\[Dash]Popov ghost field symbol"
AntiGhost::derr="Unkonw derivative structure"
	
Begin["`Private`AntiGhost`"]	


AntiGhost/:MakeBoxes[ AntiGhost[],TraditionalForm]:=OverscriptBox[SubscriptBox[ToBoxes["c"],"F"],"_"]
AntiGhost/:MakeBoxes[ AntiGhost[color:Except[_List],ds_List:{}],TraditionalForm]:=Block[{tmp,dlist={},covd=True},
If[MatchQ[ds,{___,{___,_List,___},___}],
    Message[ AntiGhost::derr];
,
    If[Length[ds]>0,
        If[And@@(MatchQ[#,_List]&/@ds),covd=False];(* if no covariant derivative involved and only partial derivative is involved *)
        (* derivative symbols *)
        tmp=ds//.{{d0___,d1:Except[_List]}:>(dlist=Prepend[dlist,SuperscriptBox["\!\(\*OverscriptBox[\(D\), \(\[LeftArrow]\)]\)",ToBoxes[d1]]];{d0}),{d0___,{dd0__,d1_}}:>(dlist=Prepend[dlist,SuperscriptBox["\!\(\*OverscriptBox[\(\[PartialD]\), \(\[LeftArrow]\)]\)",ToBoxes[d1]]];{d0,{dd0}}),{d0___,{d1_}}:>(dlist=Prepend[dlist,SuperscriptBox["\!\(\*OverscriptBox[\(\[PartialD]\), \(\[LeftArrow]\)]\)",ToBoxes[d1]]];{d0})};
        If[covd,
            RowBox[{"(",##,SuperscriptBox[")",ToBoxes[color]]}]&@@Prepend[dlist,OverscriptBox[SubscriptBox[ToBoxes["c"],"F"],"_"]]
        ,
            RowBox[Join[{SuperscriptBox[OverscriptBox[SubscriptBox[ToBoxes["c"],"F"],"_"],ToBoxes[color]]},dlist]]
        ]


    ,
        SuperscriptBox[OverscriptBox[SubscriptBox[ToBoxes["c"],"F"],"_"],ToBoxes[color]]
    ]
]
]


End[]
