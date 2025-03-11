(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



Operator::usage="Combine fields to an hadron operator"
Operator::serr="The indices for symmetrizing cannot be recoginzed."
Begin["`Private`Operator`"]	


Operator[x_,vertx:Except[_List],{ops1___,factor_ field_,ops2___},lorlist___]:=factor Operator[x,vertx,{ops1,field,ops2},lorlist]/;NumericQ[factor]


(*(*Operator[x_,vertx:Except[_List],fields:Except[_List]]:=Operator[x,vertx,{fields}]*)(* this part makes the pattern match Operator[xx_,vv_,fields_List] failed *)
Operator/:MakeBoxes[Operator[x_,vertx:Except[_List],fields_List],TraditionalForm]:=RowBox[{"(",MakeBoxes[vertx,TraditionalForm],##,")","(",ToBoxes[x],")"}]&@@(MakeBoxes[#,TraditionalForm]&/@fields)*)

(*Operator/:MakeBoxes[Operator[x_,vertx:Except[_List],fields_List,{{lor1_},{lor2_}}],TraditionalForm]:=RowBox[{"(","(",MakeBoxes[vertx,TraditionalForm],##,")","(",ToBoxes[x],")","-","{",ToBoxes[lor1],"\[LeftRightArrow]",ToBoxes[lor2],"}",")"}]&@@(MakeBoxes[#,TraditionalForm]&/@fields)
Operator/:MakeBoxes[Operator[x_,vertx:Except[_List],fields_List,{lor1_,lor2_}],TraditionalForm]:=RowBox[{"(","(",MakeBoxes[vertx,TraditionalForm],##,")","(",ToBoxes[x],")","+","{",ToBoxes[lor1],"\[LeftRightArrow]",ToBoxes[lor2],"}",")"}]&@@(MakeBoxes[#,TraditionalForm]&/@fields)*)


Operator/:MakeBoxes[Operator[x_,vertx:Except[_List],fields_List,lors_List:{}],TraditionalForm]:=Module[{tmp,lens,lorlist},
(*(* Unify the form: aslist = {lor1_List, lor2_List, ...}, slist = {lora_List, lorb_List, ...} *)
(* antisymmetrizing lor1, lor2, ..., and symmetrizing lora, lorb, ... *)
tmp=Boole[MatchQ[#,_List]]&/@aslist;
aslist = DeleteCases[Join[{ DeleteCases[(1-tmp) aslist,0]},DeleteCases[tmp aslist,0]],{}];

tmp=Boole[MatchQ[#,_List]]&/@slist;
slist = DeleteCases[Join[{ DeleteCases[(1-tmp) slist,0]},DeleteCases[tmp slist,0]],{}];*)


If[MatchQ[lors,{Except[_List],Except[_List]}|{{Except[_List]},{Except[_List]}}],
    tmp={lors}
,
    tmp=lors
];
lens=Length[tmp];

If[!And@@(MatchQ[#,{{Except[_List]},{Except[_List]}}|{Except[_List],Except[_List]}]&/@tmp),
    Message[Operator::serr];
    Abort[]
,

    lorlist=If[MatchQ[#,{{Except[_List]},{Except[_List]}}],{"-","{",ToBoxes[#[[1,1]]],"\[LeftRightArrow]",ToBoxes[#[[2,1]]],"}",")"},{"+","{",ToBoxes[#[[1]]],"\[LeftRightArrow]",ToBoxes[#[[2]]],"}",")"}]&/@tmp;
    lorlist=lorlist//Flatten;

    RowBox[Join[Table["(",lens],{" ","(",MakeBoxes[vertx,TraditionalForm],##,")","(",ToBoxes[x],")"},lorlist]]&@@(MakeBoxes[#,TraditionalForm]&/@fields)
]
]


End[]
