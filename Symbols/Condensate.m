(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



Condensate::usage = "Condensate[condensate] generate the symbol of coorsponding condensate"



Begin["`Private`Condensate`"]	
(*Options[Condensate] = {}*)



(*-------------------------------------------------------------------------------------------*)

Condensate[xx1_,xx2__]:=Condensate[xx1]Condensate[xx2]

Condensate/: MakeBoxes[Condensate[x_],TraditionalForm]:=Block[{xx,list,list1,list2,tmp},

list={{"qq",RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["q", "_"],"q","\[ThinSpace]","\[RightAngleBracket]"}]},
	{"ss",RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["s", "_"],"s","\[ThinSpace]","\[RightAngleBracket]"}]},
	{"mqq",RowBox[{"m","\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["q", "_"],"q","\[ThinSpace]","\[RightAngleBracket]"}]},
	{"msqq",RowBox[{SubscriptBox["m","s"],"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["q", "_"],"q","\[ThinSpace]","\[RightAngleBracket]"}]},
	{"mss",RowBox[{SubscriptBox["m","s"],"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["s", "_"],"s","\[ThinSpace]","\[RightAngleBracket]"}]},
	{"qgq",RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["q", "_"],"G","q","\[ThinSpace]","\[RightAngleBracket]"}]},
	{"sgs",RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["s", "_"],"G","s","\[ThinSpace]","\[RightAngleBracket]"}]},
	{"qq2",RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["q", "_"],"q",SuperscriptBox["\[RightAngleBracket]", "2"]}]},
	{"ss2",RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["s", "_"],"s",SuperscriptBox["\[RightAngleBracket]", "2"]}]},
	{"gg",RowBox[{"\[LeftAngleBracket]","\[ThinSpace]","G","G","\[ThinSpace]","\[RightAngleBracket]"}]},
	{"ggg",RowBox[{"\[LeftAngleBracket]","",SuperscriptBox["G", "3"],"","\[RightAngleBracket]"}]},
	{"d8",RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["q", "_"],"q","","\[RightAngleBracket]","\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["q", "_"],"G","q","\[ThinSpace]","\[RightAngleBracket]"}]},
	{"qq3",RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["q", "_"],"q",SuperscriptBox["\[RightAngleBracket]", "3"]}]}};
	
	
list1=Transpose[list][[1]];
list2=Transpose[list][[2]];


(* -------------------- *)
If[!FreeQ[list1,ToString[x]],
	xx=Position[list1,ToString[x]][[1,1]];
	list2[[xx]]
	,
	If[Head[x]===List,
		tmp=x;
		If[MatchQ[tmp,{_,_}|{_,"G",_}],	
			tmp[[1]]=OverBar[tmp[[1]]]
		];
		(* however, for unrecognized condensate, only bracket it *)
		
		RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",##,"","\[RightAngleBracket]"}]&@@(ToBoxes[#]&/@tmp)
	,
		RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",ToString[x],"","\[RightAngleBracket]"}]
	]
]
	
	
]


End[]
