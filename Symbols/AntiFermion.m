(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



AntiFermion::usage="AntiFermion field symbol"
AntiFermion::derr="Unkonw derivative structure"
	
Begin["`Private`AntiFermion`"]	


AntiFermion/:MakeBoxes[AntiFermion[],TraditionalForm]:= ToBoxes["\!\(\*OverscriptBox[\(\[CapitalPsi]\), \(_\)]\)"]
(*AntiFermion/:MakeBoxes[AntiFermion[x:Except[_List]],TraditionalForm]:= RowBox[{"Overscript[\[CapitalPsi], _]","(",ToBoxes[x],")"}]
AntiFermion/:MakeBoxes[AntiFermion[f_,x_],TraditionalForm]:= RowBox[{OverscriptBox[ToBoxes[f],"_"],"(",ToBoxes[x],")"}]*)
AntiFermion/:MakeBoxes[ AntiFermion[f:Except[_List]],TraditionalForm]:= OverscriptBox[ToBoxes[f],"_"]
AntiFermion/:MakeBoxes[AntiFermion[{f_}],TraditionalForm]:= SubscriptBox["\!\(\*OverscriptBox[\(\[CapitalPsi]\), \(_\)]\)",ToBoxes[f]]
(*AntiFermion/:MakeBoxes[AntiFermion[{f_},x_],TraditionalForm]:= RowBox[{SubscriptBox["Overscript[\[CapitalPsi], _]",ToBoxes[f]],"(",ToBoxes[x],")"}]*)

AntiFermion/:MakeBoxes[ AntiFermion[flavor_,dirac_,color_,ds_List:{}],TraditionalForm]:=Block[{field,fla,tmp,dlist={},covd=True},
If[MatchQ[ds,{___,{___,_List,___},___}],
    Message[ AntiFermion::derr];
,

    If[MatchQ[flavor,{_}],
        field=SubscriptBox["\!\(\*OverscriptBox[\(\[CapitalPsi]\), \(_\)]\)",ToBoxes[flavor[[1]]]]
    ,
        field=OverscriptBox[ToBoxes[flavor],"_"]
    ];


    If[Length[ds]>0||MatchQ[flavor,{_}],
        If[And@@(MatchQ[#,_List]&/@ds),covd=False];(* if no covariant derivative involved and only partial derivative is involved *)
        tmp=ds//.{{d0___,d1:Except[_List]}:>(dlist=Prepend[dlist,SuperscriptBox["\!\(\*OverscriptBox[\(\[Del]\), \(\[LeftArrow]\)]\)",ToBoxes[d1]]];{d0}),{d0___,{dd0__,d1_}}:>(dlist=Prepend[dlist,SuperscriptBox["\!\(\*OverscriptBox[\(\[PartialD]\), \(\[LeftArrow]\)]\)",ToBoxes[d1]]];{d0,{dd0}}),{d0___,{d1_}}:>(dlist=Prepend[dlist,SuperscriptBox["\!\(\*OverscriptBox[\(\[PartialD]\), \(\[LeftArrow]\)]\)",ToBoxes[d1]]];{d0})};
        If[covd,
            RowBox[{"(",##,SubsuperscriptBox[")",ToBoxes[dirac],ToBoxes[color]]}]&@@Prepend[dlist,field]
        ,
            If[MatchQ[flavor,{_}],
                RowBox[Join[{"(",field,SubsuperscriptBox[")",ToBoxes[dirac],ToBoxes[color]]},dlist]](* avoid the flavor and dirac index appears as same subscripts *)
            ,
                RowBox[Join[{SubsuperscriptBox[field,ToBoxes[dirac],ToBoxes[color]]},dlist]]
            ]
        ]

    ,
        SubsuperscriptBox[field,ToBoxes[dirac],ToBoxes[color]]
    ]

]
]


End[]
