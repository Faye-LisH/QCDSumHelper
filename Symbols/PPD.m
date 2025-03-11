(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



PPD::usage="Partial derivative of propagator in momentum space"
PPD::derr="Unrecognized derivative structures"
PPD::perr="Not a propagator in momentum space"


Begin["`Private`PPD`"]	


PPD[pp_,ds:Except[_List]]:=PPD[pp,{ds}]
PPD/:MakeBoxes[PPD[pp_,ds_List],TraditionalForm]:=Block[{tmp,null,tmpds,p,propagatorTrue=False,other,lors,dlist={}},

pp/.FermionPropagatorP[q_,___]|GluonPropagatorP[q_,___]|GhostPropagatorP[q_,___]:>(p=q;propagatorTrue=True;null);(* the momentum in propagator *)

If[!propagatorTrue,(* not a propagator in momentum space *)
    Message[PPD::perr];
    0
,

    other=Cases[FCI[ds],fvd_Pair/;!MatchQ[fvd,Pair[LorentzIndex[__],Momentum[p,___]]]:>fvd,Infinity];(* irrelevant derivatives *)
    lors=DeleteCases[FCI[ds],Alternatives@@other]/.Pair[LorentzIndex[lo_,___],Momentum[__]]:>lo;(* the LorentzIndices *)

    (*----------------------*)
    If[MatchQ[ds,{___,_List,___}]||Length[other]>0,
        Message[ PPD::derr];
        0
    ,
        (* derivative symbols *)
        tmp=SuperscriptBox["\[PartialD]",ToBoxes[#]]&/@lors;
        RowBox[{##,ToBoxes[pp,TraditionalForm]}]&@@tmp
    ]
]
]/;Length[ds]>0


End[]
