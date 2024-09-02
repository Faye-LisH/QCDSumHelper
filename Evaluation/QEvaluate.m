(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)




QEvaluate::usage = 
"QEvaluate[expr,p,EpsOrder->0] expand expr to \[Epsilon]-Series, for large input expr, QEvaluate automatically use Parallelize evaluation "

QEvaluate::onebyone = "Re-evaluate the experssion by expanding it totally and act Series[] on each term of it."
QEvaluate::drwarn = "Unknow Dirac structutre or the DiracTrace hasn't been evaluated!"


Begin["`Private`QEvaluate`"]




Options[QEvaluate] = {
	D->4-2Epsilon,
	EpsOrder->0,
	OnebyOne->False,
	Subtract->"None",
	Parallelized->"Auto",
	HoldFlavor->False(*,
	ShowasTable\[Rule]False*)
	}



QEvaluate[expr_,p_,symm_List:{0},OptionsPattern[]]/;!FreeQ[expr,DiracTrace]:=(Message[QEvaluate::drwarn];expr)


QEvaluate[expr_,p_,symm_List:{0},OptionsPattern[]]/;FreeQ[expr,DiracTrace]:=Block[
{tmp,tmp2,tmp3,tmperr,null,null1,null2,list,list2,nullist,tmplist,log,plus,nlocal=True,ieps,ilog,dimm=OptionValue[D],ord,parall=OptionValue[Parallelized],i,j,tmpexpr},


tmp=If[OptionValue[HoldFlavor]===True,
	expr//FCI//QSimplify2
,
	(expr/.FlavorDelta[fa_,fb_,___]:>If[fa===fb,1,0])//FCI//QSimplify2
];



tmp2=Cases[tmp,LorentzIndex[lo_,___]:>lo,Infinity];
If[Length[tmp2]>Length[tmp2//DeleteDuplicates],tmp=tmp//QContract];(* Contract if have dummy indices *)

If[symm=!={0},tmp=QSymmetry[tmp,symm]];

If[FreeQ[tmp,Momentum[p,___]],
	tmp
,
		
	(* FAD[x] \[Rule] 1/SPD[x] *)
	If[!FreeQ[tmp,FeynAmpDenominator],
		tmp=FeynAmpDenominatorSplit[tmp]/.{FeynAmpDenominator[PropagatorDenominator[Momentum[p,di___],mass_]]:>
												1/Pair[Momentum[p,di],Momentum[p,di]],
											FeynAmpDenominator[PropagatorDenominator[aa_ Momentum[p,di___],mass_]]:>
												1/(aa^2 Pair[Momentum[p,di],Momentum[p,di]]),
											FeynAmpDenominator[PropagatorDenominator[ Momentum[bb_ p,di___],mass_]]:>
												1/(bb^2 Pair[Momentum[p,di],Momentum[p,di]]),
											FeynAmpDenominator[PropagatorDenominator[aa_ Momentum[bb_ p,di___],mass_]]:>
												1/((aa bb)^2 Pair[Momentum[p,di],Momentum[p,di]])}
	];
	
	
	
	tmp=Expand[tmp];
	tmp=tmp/.{Pair[Momentum[p],LorentzIndex[lo_]]:>Pair[Momentum[p,D],LorentzIndex[lo,D]],Pair[Momentum[p],Momentum[zz_]]:>Pair[Momentum[p,D],Momentum[zz,D]],
				Pair[Momentum[p,D-4],LorentzIndex[lo_,D-4]]->0,Pair[Momentum[p,D-4],Momentum[zz_,D-4]]->0,Pair[LorentzIndex[_,D-4],LorentzIndex[_,D-4]]->0};
					
	tmp=ChangeDimension[tmp,D];

	tmp=QDimension[tmp,D->dimm];

	tmp2=Cases[tmp,Epsilon,Infinity];
	
	(*------------------------------------------------------------*)
	If[FreeQ[tmp,Epsilon],

		tmp=QGather[tmp,p,ShowasTable->False]
	,
		
		


(**--- Expand the experssion otherwise it's too complicate for Seies[]; avoid expand experssion too sparse, which make the Evaluation slow.
		Here the experssion is separated upto qfact2[...] by default.
		If qdelta involved, one must take Limit qdelta->0 before Series expansion, so separate experssion all.
		In some sitution, the Series gives warning Series::ztest1:  "Unable to decide whether numeric quantity ... is equal to zero. Assuming it is. "
		Expand the experssion totally can avoid problem.
		 ---**)
	If[FreeQ[tmp,qdelta],

			tmp=tmp/.{qGamma->Gamma,qfact1->Identity};
			tmp=tmp//Expand;
			tmp=tmp+null+null^2;
		
		
		(* seems mahtematica distribuate 'Series[#,{Epsilon,0,ord}]&' before assign the value of ord to 'ord'. *)	
		(* some times there's warning "Unable to decide whether numeric quantity ... is equal to zero. Assuming it is. ". Expand tmp can prevent it but it make the evaluation slow  *)
			If[OptionValue[OnebyOne]===True,
				tmp=List@@Expand[tmp/.{qfact1,qfact2->Identity}]
			,
				tmp=(List@@tmp)/.qfact2->Identity;			
			];
		
		,
		
		(* if involve qdelta, take the limit before Series[] *)
			tmp=tmp//QNormal;
			tmp=tmp//Expand;
			tmp=tmp+null+null^2;
		
			If[OptionValue[OnebyOne]===True,
				tmp=List@@Expand[tmp]
			,
				tmp={tmp}		
			];
		];
		
		
		ord=OptionValue[EpsOrder];
	
		(*-------weather or not use parallelize evaluation--------*)
		If[((Length[tmp2]>500)&&(ToString[Head[tmp]]==="Plus")&&(ToLowerCase[ToString[parall]]==="auto"))||(ToLowerCase[ToString[parall]]==="true"),
		
			ParallelEvaluate[ord=OptionValue[EpsOrder]];
			tmp=Check[ParallelMap[Series[#,{Epsilon,0,ord}]&,tmp]//Total//Normal,
				Message[QEvaluate::onebyone];
				 tmp3=List@@Plus@@Expand[tmp];
				 ParallelMap[Series[#,{Epsilon,0,ord}]&,tmp3]//Total//Normal,
				 {Series::ztest1}
				 ]
				 (* Expand the experssion totally if the Series gives Series::ztest1 warning*)
			,
			tmp=Check[Map[Series[#,{Epsilon,0,ord}]&,tmp]//Total//Normal,
				Message[QEvaluate::onebyone];
				 tmp3=List@@Plus@@Expand[tmp];
				 Map[Series[#,{Epsilon,0,ord}]&,tmp3]//Total//Normal,
				 {Series::ztest1}
				 ]
		];



		(*(**--- Expand the experssion otherwise it's too complicate for Seies[]; avoid expand experssion too sparse, which make the Evaluation slow.
		Here the experssion is separated upto qfact2[...] by default.
		If qdelta involved, one must take Limit qdelta->0 before Series expansion, so separate experssion all.
		In some sitution, the Series gives warning Series::ztest1:  "Unable to decide whether numeric quantity ... is equal to zero. Assuming it is. "
		Expand the experssion totally can avoid problem.
		 ---**)
		If[FreeQ[tmp,qdelta],

			tmp=tmp/.{qGamma->Gamma,qfact1->Identity};
			tmp=tmp//Expand;
			tmp=tmp+null+null^2;
		
		
		(* seems mahtematica distribuate 'Series[#,{Epsilon,0,ord}]&' before assign the value of ord to 'ord'. *)	
		(* some times there's warning "Unable to decide whether numeric quantity ... is equal to zero. Assuming it is. ". Expand tmp can prevent it but it make the evaluation slow  *)
			If[OptionValue[OnebyOne]===True,
				tmp={ord,#}&/@(List@@Expand[tmp/.{qfact1,qfact2->Identity}])
			,
				tmp=({ord,#}&/@(List@@tmp))/.qfact2->Identity;			
			];
		
		,
		
		(* if involve qdelta, take the limit before Series[] *)
			tmp=tmp//QNormal;
			tmp=tmp//Expand;
			tmp=tmp+null+null^2;
		
			If[OptionValue[OnebyOne]===True,
				tmp={ord,#}&/@(List@@Expand[tmp])
			,
				tmp={{ord,tmp}}		
			];
		];
		
		
		
		
		(*-------weather or not use parallelize evaluation--------*)
		If[((Length[tmp2]>500)&&(ToString[Head[tmp]]==="Plus")&&(ToLowerCase[ToString[parall]]==="auto"))||(ToLowerCase[ToString[parall]]==="true"),
		
			
			tmp=Check[ParallelMap[Series[#[[2]],{Epsilon,0,#[[1]]}]&,tmp]//Total//Normal,
				Message[QEvaluate::onebyone];
				 tmp3=List@@Plus@@Expand[tmp];
				 ParallelMap[Series[#[[2]],{Epsilon,0,#[[1]]}]&,tmp3]//Total//Normal,
				 {Series::ztest1}
				 ]
				 (* Expand the experssion totally if the Series gives Series::ztest1 warning*)
			,
			tmp=Check[Map[Series[#[[2]],{Epsilon,0,#[[1]]}]&,tmp]//Total//Normal,
				Message[QEvaluate::onebyone];
				 tmp3=List@@Plus@@Expand[tmp];
				 Map[Series[#[[2]],{Epsilon,0,#[[1]]}]&,tmp3]//Total//Normal,
				 {Series::ztest1}
				 ]
		];*)	
		
		
		tmp=(tmp//Expand)/.null->0;


		If[FreeQ[tmp,p],
		(* it's possible that dim[tmp]=0 and no momentum involved *)
			tmp//Simplify
		
		,
		(* make the result looks more clear *)(*;
			
			Plus@@((#[[1]]Plus@@#[[2]])&/@tmp)*)
			tmp=QGather[tmp,p,Subtract->OptionValue[Subtract],SelectD->True,ShowasTable->False]
		
		]
	]

]
]





End[]
