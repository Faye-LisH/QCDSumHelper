(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



AGammaD::usage =
	"AGammaD[expr] generate normalized D-dimensional gamma matrices with antisymmetry indices"


Begin["`Private`AGammaD`"]

Options[AGammaD] = {
	Explicit->True}


(*-------------------------------------------------------------------------------------------*)


AGammaD[expr___,ops___Rule]:=Signature[{expr}](AGammaD[##,ops]&@@Sort[{expr}])/;!OrderedQ[{expr}]


AGammaD[expr___,Explicit->True]:=Block[
{list,tmp,map,sign,resu,dim=OptionValue[Dimension]},
tmp={expr};

list=Table[i,{i,1,Length[tmp]}];
map=Thread[Rule[list,tmp]];

tmp=Permutations[list];
sign=Signature[#]&/@tmp;

1/Length[tmp]Total[sign(GAD@@@(tmp/.map))]
]/;OrderedQ[{expr}]


AGammaD/:MakeBoxes[AGammaD[expr___],TraditionalForm]:=SuperscriptBox["\[Gamma]",RowBox[ToBoxes[#]&/@{expr}]]
AGammaD/:MakeBoxes[AGammaD[expr___,Explicit->False],TraditionalForm]:=SuperscriptBox["\[Gamma]",RowBox[ToBoxes[#]&/@{expr}]]



End[]
