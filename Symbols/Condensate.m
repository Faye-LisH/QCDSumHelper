(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



Condensate::usage = "Condensate[condensate] generate the symbol of coorsponding condensate"
Condensate::err="Unrecongized condensate structure, the structure should be specified as {nq_NumberofQuarks, ng_NumberofGluons, nm_NumberofQuarkMasses}"


Begin["`Private`Condensate`"]	
(*Options[Condensate] = {
	Factorization->False
}*)



(*Condensate[{n_Integer,m_Integer}]:=Condensate[n]/;(n-3/2m==0)*)


(*-------------------------------------------------------------------------------------------*)

Condensate[xx1_,xx2__]:=Condensate[xx1]Condensate[xx2]

Condensate/: MakeBoxes[Condensate[x_],TraditionalForm]:=Block[{xx,list,matched,list1,list2,tmp},

list={
	{"qq",RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["q", "_"],"q","\[ThinSpace]","\[RightAngleBracket]"}]},
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
	{"qq3",RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["q", "_"],"q",SuperscriptBox["\[RightAngleBracket]", "3"]}]}
	,
(* bases of quark-condensates Subsuperscript[Q, i, n] in A.GROZIN, "METHODS OF CALCULATION OF HIGHER POWER CORRECTIONS IN QCD",International Journal of Modern Physics A 10 (1995) 3497. https://doi.org/10.1142/S0217751X95001674 *)
	{{"Q3",_},SuperscriptBox["Q","3"]},
	{{"Q5",_},SuperscriptBox["Q","5"]},
	{{"Q6",_},SuperscriptBox["Q","6"]},
	{{"Q71",_},SubsuperscriptBox["Q","1","7"]},
	{{"Q72",_},SubsuperscriptBox["Q","2","7"]},
	{{"Q73",_},SubsuperscriptBox["Q","3","7"]},
	{{"Q74",_},SubsuperscriptBox["Q","4","7"]},
	{{"A",_}, ToBoxes["A"]},
	{{"Q81",_},SubsuperscriptBox["Q","1","8"]},
	{{"Q82",_},SubsuperscriptBox["Q","2","8"]},
	{{"Q83",_},SubsuperscriptBox["Q","3","8"]},
	{{"Q84",_},SubsuperscriptBox["Q","4","8"]},
	{{"Q85",_},SubsuperscriptBox["Q","5","8"]},
	{{"Q86",_},SubsuperscriptBox["Q","6","8"]}
	};
	
	
list1=Transpose[list][[1]];
list2=Transpose[list][[2]];


(* -------------------- *)
matched=Position[list1,aa_/;MatchQ[x,aa],1];
If[Length[matched]>0,
	xx=matched[[1,1]];
	list2[[xx]]
	,

	If[Head[x]===List,
		tmp=x;
		If[MatchQ[tmp,{Except[_Integer],Except[_Integer]}|{Except[_Integer],"G",Except[_Integer]}],	
			tmp[[1]]=OverBar[tmp[[1]]];
			RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",##,"","\[RightAngleBracket]"}]&@@(ToBoxes[#]&/@tmp)
			,
			If[And@@(IntegerQ[#]&/@x),
			
			Which[
			MatchQ[x,{0,nn_Integer}|{0,nn_Integer,0}/;nn>1],
			RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",SuperscriptBox["G",x[[2]]],"\[ThinSpace]","\[RightAngleBracket]"}]
			,
			MatchQ[x,{2}|{2,0}|{2,0,0}],
			RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["\[CapitalPsi]", "_"],"\[CapitalPsi]","\[ThinSpace]","\[RightAngleBracket]"}]
			,
			MatchQ[x,{nn_?EvenQ}|{nn_?EvenQ,0}|{nn_?EvenQ,0,0}/;nn>3],
			RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",SuperscriptBox[RowBox[{"(",OverscriptBox["\[CapitalPsi]", "_"],"\[CapitalPsi]",")"}],x[[1]]/2],"\[ThinSpace]","\[RightAngleBracket]"}]
			,
			MatchQ[x,{2,1}|{2,1,0}],
			RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["\[CapitalPsi]", "_"],"G","\[CapitalPsi]","\[ThinSpace]","\[RightAngleBracket]"}]
			,
			MatchQ[x,{2,n2_Integer}|{2,n2_Integer,0}/;n2>1],
			RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["\[CapitalPsi]", "_"],"\[ThinSpace]",SuperscriptBox["G",x[[2]]],"\[CapitalPsi]","\[ThinSpace]","\[RightAngleBracket]"}]
			,
			MatchQ[x,{n1_?EvenQ,1}|{n1_?EvenQ,1,0}/;n1>2],
			RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",SuperscriptBox[RowBox[{"(",OverscriptBox["\[CapitalPsi]", "_"],"\[CapitalPsi]",")"}],x[[1]]/2],"\[ThinSpace]","G","\[ThinSpace]","\[RightAngleBracket]"}]
			,
			MatchQ[x,{n1_?EvenQ,n2_Integer}|{n1_?EvenQ,n2_Integer,0}/;n1>2&&n2>1],
			RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",SuperscriptBox[RowBox[{"(",OverscriptBox["\[CapitalPsi]", "_"],"\[CapitalPsi]",")"}],x[[1]]/2],"\[ThinSpace]",SuperscriptBox["G",x[[2]]],"\[ThinSpace]","\[RightAngleBracket]"}]
			
			,(*------------------------------*)
			
			MatchQ[x,{0,nn_Integer,m_Integer}/;nn>1&&m>0],
			RowBox[{SuperscriptBox["m",x[[3]]]/.{SuperscriptBox["m",1]->"m",SuperscriptBox["m",0]->""},"\[LeftAngleBracket]","\[ThinSpace]",SuperscriptBox["G",x[[2]]],"\[ThinSpace]","\[RightAngleBracket]"}]
			,
			MatchQ[x,{2,0,m_Integer}/;m>0],
			RowBox[{SuperscriptBox["m",x[[3]]]/.{SuperscriptBox["m",1]->"m",SuperscriptBox["m",0]->""},"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["\[CapitalPsi]", "_"],"\[CapitalPsi]","\[ThinSpace]","\[RightAngleBracket]"}]
			,
			MatchQ[x,{nn_?EvenQ,0,m_Integer}/;nn>3&&m>0],
			RowBox[{SuperscriptBox["m",x[[3]]]/.{SuperscriptBox["m",1]->"m",SuperscriptBox["m",0]->""},"\[LeftAngleBracket]","\[ThinSpace]",SuperscriptBox[RowBox[{"(",OverscriptBox["\[CapitalPsi]", "_"],"\[CapitalPsi]",")"}],x[[1]]/2],"\[ThinSpace]","\[RightAngleBracket]"}]
			,
			MatchQ[x,{2,1,m_Integer}/;m>0],
			RowBox[{SuperscriptBox["m",x[[3]]]/.{SuperscriptBox["m",1]->"m",SuperscriptBox["m",0]->""},"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["\[CapitalPsi]", "_"],"G","\[CapitalPsi]","\[ThinSpace]","\[RightAngleBracket]"}]
			,
			MatchQ[x,{2,n2_Integer,m_Integer}/;n2>1&&m>0],
			RowBox[{SuperscriptBox["m",x[[3]]]/.{SuperscriptBox["m",1]->"m",SuperscriptBox["m",0]->""},"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["\[CapitalPsi]", "_"],"\[ThinSpace]",SuperscriptBox["G",x[[2]]],"\[CapitalPsi]","\[ThinSpace]","\[RightAngleBracket]"}]
			,
			MatchQ[x,{n1_?EvenQ,1,m_Integer}/;n1>2&&m>0],
			RowBox[{SuperscriptBox["m",x[[3]]]/.{SuperscriptBox["m",1]->"m",SuperscriptBox["m",0]->""},"\[LeftAngleBracket]","\[ThinSpace]",SuperscriptBox[RowBox[{"(",OverscriptBox["\[CapitalPsi]", "_"],"\[CapitalPsi]",")"}],x[[1]]/2],"\[ThinSpace]","G","\[ThinSpace]","\[RightAngleBracket]"}]
			,
			MatchQ[x,{n1_?EvenQ,n2_Integer,m_Integer}/;n1>2&&n2>1&&m>0],
			RowBox[{SuperscriptBox["m",x[[3]]]/.{SuperscriptBox["m",1]->"m",SuperscriptBox["m",0]->""},"\[LeftAngleBracket]","\[ThinSpace]",SuperscriptBox[RowBox[{"(",OverscriptBox["\[CapitalPsi]", "_"],"\[CapitalPsi]",")"}],x[[1]]/2],"\[ThinSpace]",SuperscriptBox["G",x[[2]]],"\[ThinSpace]","\[RightAngleBracket]"}]
			
			,(* unrecognized input *)
			True,
			Message[Condensate::err];
			RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",SubscriptBox["\[ScriptCapitalO]",RowBox[Flatten[{ToBoxes[#],","}&/@x]]],"\[ThinSpace]","\[RightAngleBracket]"}]
			]
			,
			Message[Condensate::err];
			RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",SubscriptBox["\[ScriptCapitalO]",RowBox[Flatten[{ToBoxes[#],","}&/@x]]],"\[ThinSpace]","\[RightAngleBracket]"}]
			]
		]
		
	,
	If[EvenQ[x]&&x>1,
			If[x==2,
				RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",OverscriptBox["\[CapitalPsi]", "_"],"\[CapitalPsi]","\[ThinSpace]","\[RightAngleBracket]"}]
			,
				RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",SuperscriptBox[RowBox[{"(",OverscriptBox["\[CapitalPsi]", "_"],"\[CapitalPsi]",")"}],x/2],"\[ThinSpace]","\[RightAngleBracket]"}]
			]
		,
			Message[Condensate::err];
			RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",SubscriptBox["\[ScriptCapitalO]",ToBoxes[x]],"\[ThinSpace]","\[RightAngleBracket]"}]
		]
	]
]
	
	
]


End[]
