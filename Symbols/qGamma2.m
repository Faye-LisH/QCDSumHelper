(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



qGamma2::usage = 
"qGamma2[expr] is a special qGamma function, which show \[Epsilon]-pole explicitly";

Options[qGamma2]={
	D->4-2Epsilon,
	Expand->False};
	
Begin["`Private`qGamma2`"]


qGamma2[y_,OptionsPattern[]]:=Block[
{t,n,tmp,dim=OptionValue[D],eps=0,order=OptionValue[Expand],result},

(*If[FreeQ[y,D]&&FreeQ[y,Epsilon]&&FreeQ[y,_Integer],Print["check input!"];Abort[]];*)

t = y/.D->dim//Expand;

n = t/. Epsilon->0;
eps = t - n ;

If[n<=0 ,
	tmp=(-1)^n /eps Gamma[1-eps]Gamma[1+ eps]/Gamma[-n+1 -eps];
	
	If[IntegerQ[order]&&order>=0, 
		tmp=Series[tmp,{Epsilon,0,order}]//Normal
		(*result=Sum[tmp[[3]][[i]] tmp[[1]]^(tmp[[4]]+i-1),{i,1,tmp[[5]]-tmp[[4]]}],
		result = tmp*)
	];
	result=tmp
,
		
	result=Gamma[t]
];

result	
]



End[]
