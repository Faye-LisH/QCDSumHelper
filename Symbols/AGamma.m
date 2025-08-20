(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



AGamma::usage =
	"AGamma[expr] generate normalized gamma matrices with antisymmetry indices"
	

Begin["`Private`AGamma`"]

Options[AGamma] = {
	Explicit->True}


(*-------------------------------------------------------------------------------------------*)


AGamma[expr___,ops___Rule]:=Signature[{expr}](AGamma[##,ops]&@@Sort[{expr}])/;!OrderedQ[{expr}]


AGamma[expr___,Explicit->True]:=Block[
{list,tmp,map,sign,resu,dim=OptionValue[Dimension]},
tmp={expr};

list=Table[i,{i,1,Length[tmp]}];
map=Thread[Rule[list,tmp]];

tmp=Permutations[list];
sign=Signature[#]&/@tmp;

1/Length[tmp]Total[sign(GA@@@(tmp/.map))]

]/;OrderedQ[{expr}]


AGamma/:MakeBoxes[AGamma[expr___],TraditionalForm]:=SuperscriptBox["\!\(\*OverscriptBox[\(\[Gamma]\), \(_\)]\)",RowBox[ToBoxes[#]&/@{expr}]]
AGamma/:MakeBoxes[AGamma[expr___,Explicit->False],TraditionalForm]:=SuperscriptBox["\!\(\*OverscriptBox[\(\[Gamma]\), \(_\)]\)",RowBox[ToBoxes[#]&/@{expr}]]


End[]
