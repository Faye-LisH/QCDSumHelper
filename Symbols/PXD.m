(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



PXD::usage="Partial derivative of propagator in coordinate space"
PXD::derr="Unrecognized derivative structures"
PXD::perr="Not a propagator in coordinate space"


Begin["`Private`PXD`"]	


PXD[pp_,ds:Except[_List]]:=PXD[pp,{ds}]
PXD/:MakeBoxes[PXD[pp_,ds_List],TraditionalForm]:=Block[{tmp,null,tmpds,x1,x2,propagatorTrue=False,other,lors,sign,dlist={}},

pp/.{FermionPropagatorX[x:Except[_List],___]|GluonPropagatorX[x:Except[_List],___]|GhostPropagatorX[x:Except[_List],___]:>(x1=x;propagatorTrue=True;null),
FermionPropagatorX[{x_,y_},___]|GluonPropagatorX[{x_,y_},___]|GhostPropagatorX[{x_,y_},___]:>(x1=x;x2=y;propagatorTrue=True;null)
};(* the positions in propagator *)

If[!propagatorTrue,(* not a propagator in coordinate space *)
    Message[PXD::perr];
    0
,

    other=Cases[FCI[ds],fvd_Pair/;!MatchQ[fvd,Pair[LorentzIndex[__],Momentum[x1|x2|x1-x2|x2-x1,___]]]:>fvd,Infinity];(* irrelevant derivatives *)
    lors=DeleteCases[FCI[ds],Alternatives@@other]/.{Pair[LorentzIndex[lo_,___],Momentum[x1|x1-x2,___]]:>lo,Pair[LorentzIndex[lo_,___],Momentum[x2|x2-x1,___]]:>{-1,lo}};(* the LorentzIndices *)
    (* treat the indices v as d^v_x1 by defualt *)


    (* for free propagator, propagator[x1,x2]=propagator[x1-x2], d_x2 = -d_x1 *)
    sign=Times@@Cases[lors,{-1,_}:>-1,Infinity];
    lors=lors/.{-1,lo_}:>lo;(* keep the indices after extract the overall sign *)

    (*----------------------*)
    If[MatchQ[ds,{___,_List,___}]||Length[other]>0,
        Message[ PXD::derr];
        0
    ,
        (* derivative symbols *)
        tmp=SuperscriptBox["\[PartialD]",ToBoxes[#]]&/@lors;
        If[sign==1,
            RowBox[{##,ToBoxes[pp,TraditionalForm]}]&@@tmp
        ,
            RowBox[{"-",##,ToBoxes[pp,TraditionalForm]}]&@@tmp
        ]
    ]
]
]/;Length[ds]>0


End[]
