(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)



(* ::Code::Initialization::Plain:: *)
QGather::usage =
"QGather[expr,p,ShowasTable\[Rule]True,
	SeparateTransverse\[Rule]True] is a function gather the expr to a more symmetry form about Lorentz structure"


(* ::Code::Initialization::Plain:: *)
QGather::trwarn="The DiracTrace hasn't been evaluated!"

QGather::lowarn="Lorentz structure inconsistent! make sure the input is correct."


(* ::Code::Initialization::Plain:: *)
Begin["`Private`QGather`"]


(* ::Code::Initialization::Plain:: *)
Options[QGather] = {
	ShowasTable->"Auto",
	Transverse->True,
	Subtract->"None",
	SelectD->False}	


(* ::Code::Initialization::Plain:: *)
QGather[expr_List, p:Except[_Rule],rules___Rule]/;And@@(MatchQ[#,{_,_List}]&/@expr):=Plus@@((#[[1]](Plus@@#[[2]]))&/@expr)


(* ::Code::Initialization::Plain:: *)
(*------------------------------*)
QGather[expr_,rules___Rule]:=Block[{tmp},
tmp=Cases[expr//FCI,Pair[Momentum[q_,___],Momentum[q_,___]]|Pair[Momentum[q_,___],LorentzIndex[___]]:>q,Infinity];
tmp=DeleteDuplicates[tmp];

If[Length[tmp]==1,
	QGather[expr,tmp[[1]],rules]
,
	expr
]

]


(* ::Code::Initialization::Plain:: *)
QGather[expr:Except[_List],p:Except[_Rule],OptionsPattern[]]/;!FreeQ[expr,DiracTrace]:=(Message[QGather::trwarn];expr)


QGather[expr:Except[_List],p:Except[_Rule],OptionsPattern[]]/;FreeQ[expr,DiracTrace]:=Block[
{tmp,tmp2,tmplor,tmpcoe,tmprem,tmp0,list={},list2,count,dot,log,llog,null,null1,null0,null2,nulllor,i,n={},l,nn,m,k,ldim,lorentzIndex,
tf=OptionValue[ShowasTable],separate=OptionValue[Transverse],subtract=OptionValue[Subtract],nullist,fact,spdfad=False,ord,nullcondensate},


tmp=expr//FCI//Expand;

(*--- Contract if dummy lorentzindex involved ---*)
tmp2=Cases[tmp,LorentzIndex[lo_,___]:>lo,Infinity];
If[Length[tmp2]>Length[tmp2//DeleteDuplicates],tmp=tmp//Contract];


If[FreeQ[tmp,Momentum[p,___]],
	tmp
,
(*--- dicsard the D-4 dimension terms, write the 4 dimension terms as D dimension ---*)
	If[OptionValue[SelectD],
		tmp=tmp//Expand;
		tmp=tmp/.{Pair[Momentum[p],LorentzIndex[lo_]]:>Pair[Momentum[p,D],LorentzIndex[lo,D]],Pair[Momentum[p],Momentum[zz_]]:>Pair[Momentum[p,D],Momentum[zz,D]],
					Pair[Momentum[p,D-4],LorentzIndex[lo_,D-4]]->0,Pair[Momentum[p,D-4],Momentum[zz_,D-4]]->0};
		tmp=ChangeDimension[tmp,D]	
	];
	
	
(*		(* FAD[x] \[Rule] 1/SPD[x] *)
	If[!FreeQ[tmp,FeynAmpDenominator],
		spdfad=True;
		tmp=FeynAmpDenominatorSplit[tmp]/.{FeynAmpDenominator[PropagatorDenominator[Momentum[p,di___],mass_]]:>
												1/Pair[Momentum[p,di],Momentum[p,di]],
											FeynAmpDenominator[PropagatorDenominator[aa_ Momentum[p,di___],mass_]]:>
												1/(aa^2 Pair[Momentum[p,di],Momentum[p,di]]),
											FeynAmpDenominator[PropagatorDenominator[ Momentum[bb_ p,di___],mass_]]:>
												1/(bb^2 Pair[Momentum[p,di],Momentum[p,di]]),
											FeynAmpDenominator[PropagatorDenominator[aa_ Momentum[bb_ p,di___],mass_]]:>
												1/((aa bb)^2 Pair[Momentum[p,di],Momentum[p,di]])}
	];
*)
	
(*--- check the number of free Lorentz indeices ---*)
	tmp=tmp//Expand;
	If[Head[tmp]===Plus,
		tmp=(# count[Cases[#,LorentzIndex[lo_,___]:>lo,Infinity]])&/@tmp;
		tmp=tmp/.count[aa_]:>count^(2 Length[DeleteDuplicates[aa]]-Length[aa])
	];
	
	n=Exponent[tmp,count,List];
	If[(Length[n]=!=1),Message[QGather::lowarn]];
	n=n[[1]];
	
	tmp=tmp/.count->1;

	
(*--- refine the log term make it more compact as usuall form ---*)
	tmp=tmp/.{Log[Power[Pair[Momentum[p,di___],Momentum[p,di___]],po_Integer]]/;Negative[po]:> -Log[Power[Pair[Momentum[p,D],Momentum[p,D]],-po]],
						Log[aa_ Power[Pair[Momentum[p,di___],Momentum[p,di___]],po_Integer]]/;Negative[po]:> -Log[Power[Pair[Momentum[p,D],Momentum[p,D]],-po]/aa]};
						
(*--- isolate the log term if it already show as intended ---*)					
	tmp=tmp/.{Log[-Pair[Momentum[p,dim___],Momentum[p,dim___]]Power[4Pi ScaleMu^2,-1]]:>llog[-Pair[Momentum[p,dim],Momentum[p,dim]]Power[4Pi ScaleMu^2,-1]],
				Log[-Pair[Momentum[p,dim___],Momentum[p,dim___]]Power[ScaleMu^2,-1]]:>llog[-Pair[Momentum[p,dim],Momentum[p,dim]]Power[ScaleMu^2,-1]]};				
	
											
	tmp=tmp/.Log[x_]/;FreeQ[x,Pair]:>PowerExpand[Log[x]];
	tmp=tmp/.Log[Pair[Momentum[p,dim___],Momentum[p,dim___]]]:>log[-Pair[Momentum[p,dim],Momentum[p,dim]]/(4\[Pi] ScaleMu^2)]+2 Log[2]+Log[Pi]+2 Log[ScaleMu]- Pi I;
	
	tmp=tmp/.Log[a_ Pair[Momentum[p,dim___],Momentum[p,dim___]]]:>log[-Pair[Momentum[p,dim],Momentum[p,dim]]/(4\[Pi] ScaleMu^2)]+2 Log[2]+Log[Pi]+2 Log[ScaleMu]+ PowerExpand[Log[-a]];
	
	tmp=tmp/.Log[a_ Power[Pair[Momentum[p,dim___],Momentum[p,dim___]],-1]]:>-log[-Pair[Momentum[p,dim],Momentum[p,dim]]/(4\[Pi] ScaleMu^2)]-2 Log[2]-Log[Pi]-2 Log[ScaleMu]+ PowerExpand[Log[-a]];


	If[ToLowerCase[ToString[subtract]]=="msbar",
		tmp=tmp/.{log[a_ Pair[Momentum[p,dim___],Momentum[p,dim___]]]:>
						log[-Pair[Momentum[p,dim],Momentum[p,dim]]/(ScaleMu^2)]+2 Log[2]+Log[Pi]+2 Log[ScaleMu]+ PowerExpand[Log[-a]]-EulerGamma,
				log[a_ Power[Pair[Momentum[p,dim___],Momentum[p,dim___]],-1]]:>
						-log[-Pair[Momentum[p,dim],Momentum[p,dim]]/(ScaleMu^2)]-2 Log[2]-Log[Pi]-2 Log[ScaleMu]+ PowerExpand[Log[-a]]+EulerGamma,
				llog[-Pair[Momentum[p,dim___],Momentum[p,dim___]]/(4Pi ScaleMu^2)]:>llog[-Pair[Momentum[p,dim],Momentum[p,dim]]/(ScaleMu^2)]-EulerGamma}
		];

	tmp=tmp/.{llog->log,Log[a_]:>PowerExpand[Log[a]]};
(*--------------------------------------------------------------------------------*)

	(* to separate transverse and longitudinal part,
	first MTD[a,b]-> MTD[a,b] + FVD[p,a]FVD[p,b]/SPD[p], then MTD[a,b]-> MTD[a,b]- FVD[p,a]FVD[p,b]/SPD[p] in the end *)
	
	(* times a vector so that below can also deal with scalar input *)
	tmp=tmp Pair[Momentum[p,D],LorentzIndex[nulllor,D]];
	
	If[separate==True,tmp=tmp/.Pair[LorentzIndex[a_,dim___],LorentzIndex[b_,dim___]]:>
				Pair[LorentzIndex[a,dim],LorentzIndex[b,dim]]+Pair[Momentum[p,dim],LorentzIndex[a,dim]]Pair[Momentum[p,dim],LorentzIndex[b,dim]]/Pair[Momentum[p,dim],Momentum[p,dim]]];



	(* gather to a list *)
	If[n>0,
		tmp=tmp//Expand;
		tmp=tmp/.Dot->dot;
		
	(* //////// get the terms have same Lorentz structure with first term //////////*)
		While[Head[tmp]===Plus,
			tmpcoe=tmp[[1]]/.{Eps[__]:>1,Pair[LorentzIndex[__],_]:>1,DiracGamma[LorentzIndex[__],___]:>1,DiracGamma[Momentum[_,___],___]:>1}/._dot->1;
			tmplor=tmp[[1]]/tmpcoe;

			tmprem=tmp/.tmplor:>0;
			tmp0=tmp-tmprem/.{Eps[__]:>1,Pair[LorentzIndex[__],_]:>1,DiracGamma[LorentzIndex[__],___]:>1,DiracGamma[Momentum[_,___],___]:>1}/._dot->1;

			tmp=tmprem;
			list=Append[list,{tmplor,tmp0}]
			];

		If[tmp=!=0,
			tmpcoe=tmp/.{Eps[__]:>1,Pair[LorentzIndex[__],_]:>1,DiracGamma[LorentzIndex[__],___]:>1,DiracGamma[Momentum[_,___],___]:>1}/._dot->1;
			tmplor=tmp/tmpcoe;
			list=Append[list,{tmplor,tmpcoe}]
			];
	,

		list={{1,tmp}}
	];
	
	

(* sorting by power of Epsilon  *)
	ord=Min[Cases[list,Power[Epsilon,pe_]:>pe,Infinity]];
	If[ord<0, ord=-ord , ord=0 ];
	
	For[l=1,l<Length[list]+1,l++,
		list[[l,2]]= CoefficientList[Epsilon^(ord+1)list[[l,2]],Epsilon];
		nn=Length[list[[l,2]]];
		list[[l,2]]=DeleteCases[list[[l,2]]Table[Epsilon^(k-(ord+1)),{k,0,nn-1}],0];
		];


	(* gather longitude part and transverse part*)
	If[separate==True,list=list/.Pair[LorentzIndex[a_,dim___],LorentzIndex[b_,dim___]]:>
				Pair[LorentzIndex[a,dim],LorentzIndex[b,dim]]-Pair[Momentum[p,dim],LorentzIndex[a,dim]]Pair[Momentum[p,dim],LorentzIndex[b,dim]]/Pair[Momentum[p,dim],Momentum[p,dim]]];


	(*----------- subtract Epsilon-pole  ------------*)
	If[((ToLowerCase[ToString[subtract]]=="ms")||(ToLowerCase[ToString[subtract]]=="msbar")) && 
		FreeQ[list//Expand,log[___]Power[Epsilon,_?Negative]]&&FreeQ[list//Expand,Power[log[___],___]Power[Epsilon,_?Negative]],
		list=list/.Power[Epsilon,_?Negative]->0;
		list={#[[1]],Replace[#[[2]],{0..,aa__}:>{aa},{0}]}&/@list
	];
	
	

(**--- Gather condensate ---**)
	list=list/.{gStrong->Condensate[gStrong],lg_log:>Condensate[Log]lg};
	(* write g_s as <g_s> and Log as <Log> Log, so that LO and NLO can be separated by these code *)
	list={#[[1]],Condensate[nullcondensate]#[[2]]}&/@list;
	list={#[[1]],((Plus@@#)&/@Gather[List@@(Expand[#]+null0+null0^2),#1/(#1/.Condensate[_]->1)===#2/(#2/.Condensate[_]->1)&])&/@#[[2]]}&/@list;
	(* e.g. {{q^u, {<Condensate>/Epsilon, ... }} , ...} ->  {{q^u, {{<GG>/Epsilon, <qgq>/Epsilon}, ... }}, ...} *)



(*-------------------------Simplify PolyGamma--------------------------------*)
	If[!FreeQ[list,PolyGamma],list=list//FullSimplify];
			

	(* simplify each term *)		
    For[i=1,i<Length[list]+1,i++,
        list[[i,2]]=Simplify[list[[i,2]]]
    ];			
			
	
	(***--- recover the structure of list to the form before (**--- Gather condensate ---**)  ---***)
	list={#[[1]],Plus@@#&/@#[[2]]}&/@list;
	
	
	
	(*-----------------------------------------------------------------------*)
	tmp=list/.{null0->0,null1->Identity,log->Log,Pair[Momentum[p,D],LorentzIndex[nulllor,D]]->1,Condensate[nullcondensate]->1};
	

	If[!(Length[tmp]==1&&NumberQ[tmp[[1,1]]]),
		(* if it is not a scalar, collect the terms only differ by the tensor structure *)
		tmp=Flatten[Gather[tmp,(Expand[Plus@@(#1[[2]])+Plus@@(#2[[2]])]===0)||(Expand[Plus@@(#1[[2]])-Plus@@(#2[[2]])]===0)&],1]
	];
	
	

	(*------ whether show as a table ------*)
(*	If[tf==="ForcetoTable",
		tf=True;

	,
		If[tf==="Auto"&&Length[tmp]==1&&NumberQ[tmp[[1,1]]],tf=False](* for scalar, not show it as table by default *)
	];*)
	
	If[tf==="Auto",
		If[Length[tmp]==1&&NumberQ[tmp[[1,1]]],
			tf=False
		,
			tf=True
		]
	];(* for scalar, not show it as a table by default *)
	
	
	If[tf,
 	   tmp
	,
		Plus@@(#[[1]]Plus@@(#[[2]])&/@tmp)
	]/.dot->Dot
	
]/.{Condensate[gStrong]->gStrong,Condensate[Log]->1}

]



(* ::Code::Initialization::Plain:: *)
End[]
