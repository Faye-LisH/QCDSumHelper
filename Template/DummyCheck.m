(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


DummyCheck::usage="DummyCheck[v1_,v2_,v3_,v4_] checks whether the dummy indices are valid, if valid(invalid), it gives 1(0)."
DummyCheck::lorerr="Different sets of Lorentz indices involved in a single vertex!"
DummyCheck::sunerr="Different sets of SU(N) indices involved in a single vertex!"

(*ToFile::usage=="Give a directory to Export and Import the results."*)


Begin["`Private`DummyCheck`"]


DummyCheck[v1_,v2_,v3_,v4_]:=Block[{indices,indlist,vlist,lorlist,sunlist,consist=1},
vlist=DotSimplify[#]&/@({v1,v2,v3,v4}//FCI//Expand);
vlist=If[Head[#]===Plus,List@@#,{#}]&/@vlist;(* separate each term in vertices *)

(* the lorentz indices and sun indices for each term in each vertex *)
lorlist=Cases[#,LorentzIndex[lor_,___]:>lor,Infinity]&/@#&/@vlist;
sunlist=Cases[#,SUNIndex[sun_,___]:>sun,Infinity]&/@#&/@vlist;

(* discard terms with same indices in each vertex *)
lorlist=DeleteDuplicates[#]&/@lorlist;
sunlist=DeleteDuplicates[#]&/@sunlist;

(* each vertex should has only 1 set of indices *)
If[Or@@(Length[#]>1&/@lorlist),
Message[DummyCheck::lorerr];
consist=0
];
If[Or@@(Length[#]>1&/@sunlist),
Message[DummyCheck::sunerr];
consist=0
];

(* the number of dummy indices *)
lorlist=Gather[Flatten[lorlist]];
sunlist=Gather[Flatten[sunlist]];

lorlist=Length[#]&/@lorlist;
sunlist=Length[#]&/@sunlist;

If[FreeQ[lorlist,nn_Integer/;nn>2]&&FreeQ[sunlist,nn_Integer/;nn>2],
	1
,
	0
]

]


End[]
(*EndPackage[]*)
