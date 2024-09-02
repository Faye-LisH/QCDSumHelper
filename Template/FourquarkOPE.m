(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
FourquarkOPE::usage="FourquarkOPE[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}}] give the Operator Product Expansion (OPE) of <J1 J2>, with J1=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a1\)]\)\!\(\*SubscriptBox[\(v1\[CapitalPsi]\), \(b1\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c1\)]\)\!\(\*SubscriptBox[\(v2\[CapitalPsi]\), \(d1\)]\) and J2=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(b2\)]\)\!\(\*SubscriptBox[\(v3\[CapitalPsi]\), \(a2\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(d2\)]\)\!\(\*SubscriptBox[\(v4\[CapitalPsi]\), \(c2\)]\) "
FourquarkOPE::inderr="Dummy indices conflict!"
FourquarkOPE::ipterr="Unknow expression!"
FourquarkOPE::cerr="Can't recongize the assigned condensates"
FourquarkOPE::mwarn="Momentum not match, momentum `1` in the files will be replaced to `2`"
FourquarkOPE::merr="Fail to specify the momentum"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`FourquarkOPE`"]

Options[FourquarkOPE]={
	NLO->True,
	Renormalization->True,
	Parallelized->True,
	Strategy->"Fourier",
	ShowProgress->True,
	AutoNDR->"Auto",
	Condensates->All,
	HoldFlavor->True,
	EpsOrder->0
}



FourquarkOPE[qq_:"auto",cur1:Except[List],cur2:Except[List],ops:OptionsPattern[]]/;!FreeQ[cur1,FourquarkCurrent]&&!FreeQ[cur2,FourquarkCurrent]:=Block[{null,c1,tmp1,c2,tmp2,current},


(* unify the input; turn tetraquark current to fourquark molecule current *)
tmp1=Expand[null cur1/.{FourquarkCurrent[aa_Current,bb_Current]/;
		MatchQ[aa,Current[{_,_},{_,_}]|Current[_,_,_,_]|Current[{_,_},_,{_,_}]|Current[_,_,_,_,_]]&&MatchQ[bb,Current[{-1,_,_},{-1,_,_}]|Current[-1,_,_,-1,_,_]|Current[{-1,_,_},_,{-1,_,_}]|Current[-1,_,_,_,-1,_,_]]:>
		FierzTransformation[FourquarkCurrent[aa,bb]]
	}];

		
tmp2=Expand[null cur2/.{FourquarkCurrent[aa_Current,bb_Current]/;
		MatchQ[aa,Current[{_,_},{_,_}]|Current[_,_,_,_]|Current[{_,_},_,{_,_}]|Current[_,_,_,_,_]]&&MatchQ[bb,Current[{-1,_,_},{-1,_,_}]|Current[-1,_,_,-1,_,_]|Current[{-1,_,_},_,{-1,_,_}]|Current[-1,_,_,_,-1,_,_]]:>
		FierzTransformation[FourquarkCurrent[aa,bb]]
	}];

(*------------------------*)
(* to a list of {{prefactor,{fa,v1,fb},{fc,v2,fd}},...} *)

tmp1=topair[tmp1];
tmp2=topair[tmp2];

(*------------------------------------------------------*)

Outer[#1[[1]]#2[[1]]FourquarkOPE[qq,{#1[[2,1,2]],#1[[2,2,2]],{#1[[2,1,1]],#1[[2,1,3]],#1[[2,2,1]],#1[[2,2,3]]}},{#2[[2,1,2]],#2[[2,2,2]],{#2[[2,1,1]],#2[[2,1,3]],#2[[2,2,1]],#2[[2,2,3]]}},ops]&,tmp1,tmp2,1]//Total//Total

]



topair[expr_]:=Block[{c1,null,tmp},

If[MatchQ[expr,aa_ cc_FourquarkCurrent],
	c1={expr}/._FourquarkCurrent->1;(* a list of prefactors *)
	tmp={expr}/c1(* a list of currents *)
,
	If[And@@(MatchQ[#,aa_ cc_FourquarkCurrent]&/@(List@@expr)),
	
		c1=(List@@expr)/._FourquarkCurrent->1;
		tmp=(List@@expr)/c1;
	,
	
		Message[FourquarkOPE::ipterr];
		Abort[]
	]

];


({c1/.null->1,tmp}//Transpose)/.FourquarkCurrent[aa_Current,bb_Current]:>{List@@aa,List@@bb}
]



(*-----------------------------------------------*)

FourquarkOPE[qq_:"auto",{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmp,tmp1,tmp2,tmp0,poles1,poles2,strategy,nlo=True,lofact=1,nlofact=1,renorm=OptionValue[Renormalization],evlist,evlistlength,
condensates=OptionValue[Condensates],con1=0,conmqq=0,congg=0,conqq2=0,conmqgq=0,conggg=0,cond81=0,cond82=0,cond83=0,cond84=0,cond85=0,mqqgg=0,cond100=0,cond101=0,cond102=0,cond103=0,condmqq3=0,ndr1,ndr,ndr0,pall=OptionValue[Parallelized],fdir,parallelized,holdf=OptionValue[HoldFlavor],tmpq,conds,types,d01,d02,d03},

If[qq==="auto"&&!DirectoryQ[OptionValue[Parallelized]//ToString],Message[FourquarkOPE::merr];Abort[]];
(* momentum must be specified, except when load the date from a directory *)


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[FourquarkOPE::inderr];
	Abort[]
];

(*---------------------------------------------------*)

(* The calculation about next FourquarkLeading order diagrams is time-consuming, so one may need different strategies for different cases. 
In most cases the fourier transformation can give a best performance, since the involved vector is only the x.
 In "fouriertransformation" strategy, the evaluation is: coordinate space propagator \[Rule] Trace \[Rule] fouriertransform to momentumspace;
 In "integrationfirst" strategy, the integration is evaluate before evaluation of TR if possible, 
 use this abnormal way is because if the TR involving many \[Gamma]^5 matrices, in BMHV scheme, the results involve terms in D, 4, and D-4 dimensional objects, 
 compare with other \[Gamma]^5 scheme, for a diagram have n propagators, this 3 possible situations make the complexity of the result of TR is roughly 3^(n/2) times than usual;
 integrate before evaluate TR may save the time in some situations, since after integrate, only one momentum is involved.
 The "tracefirst" strategy is the usually way: Trace \[Rule] Integrate *)
Which[
	ToLowerCase[ToString[OptionValue[Strategy]]]==="fourier"||ToLowerCase[ToString[OptionValue[Strategy]]]==="fouriertransformation",
	strategy="fourier",
	
	ToLowerCase[ToString[OptionValue[Strategy]]]==="integrate"||ToLowerCase[ToString[OptionValue[Strategy]]]==="integratefirst"||ToLowerCase[ToString[OptionValue[Strategy]]]==="integrationfirst",
	strategy="integrationfirst",
	
	ToLowerCase[ToString[OptionValue[Strategy]]]==="trace"||ToLowerCase[ToString[OptionValue[Strategy]]]==="tracefirst",
	strategy="tracefirst"
];



(*---------------- options & setting -------------------*)
Which[
	OptionValue[NLO]===False, nlofact=0,
	OptionValue[NLO]==="NLOonly", lofact=0
];


If[StringQ[condensates],condensates={condensates}];(* allow string input for one type of condensate *)
(* specify the condensates *)

conds={"d0_1","d0_2","d0_3","mqq","gg","qq2","mqgq","ggg","d8_1","d8_2","d8_3","d8_4","mqqgg","d10_0","d10_1","d10_2","d10_3","mqq3"};
If[Head[condensates]===List,
	{d01,d02,d03,conmqq,congg,conqq2,conmqgq,conggg,cond81,cond82,cond83,cond84,mqqgg,cond100,cond101,cond102,cond103,condmqq3}=Boole[!FreeQ[condensates,#]]&/@conds;	
	If[(Plus@@{d01,d02,d03})=!=0, con1=1;types={"d0_1","d0_2","d0_3"}{d01,d02,d03}];
	If[!FreeQ[condensates,"d0",{1}],con1=1;types="All"];
	If[!FreeQ[condensates,"d4",{1}],conmqq=1;congg=1];
	If[!FreeQ[condensates,"d6",{1}],conmqq=1;congg=1;conmqgq=1];
	If[!FreeQ[condensates,"d8",{1}],cond81=1;cond83=1;mqqgg=1];
	If[!FreeQ[condensates,"d10",{1}],cond100=1;cond101=1;cond102=1;cond103=1;condmqq3=1]
	
,
	If[condensates===All,
		con1=1;conmqq=1;congg=1;conqq2=1;conmqgq=1;conggg=1;cond81=1;cond83=1;mqqgg=1;cond100=1;cond101=1;cond102=1;cond103=1;condmqq3=1;types="All"
	,
		Message[FourquarkOPE::cerr];
		Abort[]
	]
];(* select the condensates *)


(* the gammma^5 scheme is relevant for evaluate perturbative contribution *)
If[OptionValue[AutoNDR]==="Auto",
	ndr1=False;
	If[OptionValue[EpsOrder]===0,ndr0=True,ndr0=False];
	ndr=True
,
	ndr=OptionValue[AutoNDR];
	ndr1=ndr;
	ndr0=ndr
];	


If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];

(*----------------------------------------*)


tmp={con1 lofact Hold[FourquarkLeading[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},AutoNDR->ndr0,Parallelized->parallelized,HoldFlavor->holdf,EpsOrder->OptionValue[EpsOrder]]],
	con1 nlofact Hold[FourquarkAlphas[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},AutoNDR->ndr1,Renormalization->renorm,Parallelized->parallelized,Strategy->strategy,HoldFlavor->holdf,AlphasTypes->types]],
	congg Hold[Fourquarkgg[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	conmqq Hold[Fourquarkmqq[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	conggg Hold[Fourquarkggg[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	conqq2 Hold[Fourquarkqq2[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond81 Hold[Fourquarkd8type1[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond82 Hold[Fourquarkd8type2[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond83 Hold[Fourquarkd8type3[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond84 Hold[Fourquarkd8type4[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond85 Hold[Fourquarkd8type5[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	mqqgg Hold[Fourquarkmqqgg[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	conmqgq Hold[Fourquarkmqgq[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond100 Hold[Fourquarkd10type0[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond101 Hold[Fourquarkd10type1[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond102 Hold[Fourquarkd10type2[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond103 Hold[Fourquarkd10type3[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	condmqq3 Hold[Fourquarkmqq3[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]]
	};


(*----------------------------------------*)

If[pall===True,

	tmp=Plus@@#&/@WaitAll[ReleaseHold[tmp/.parallelized->"External"]];
	Plus@@((QGather[(#/.{CA->SUNN,CF->(SUNN^2-1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0})&/@tmp)
	

,
	If[pall==="External",
	
		If[fdir==="None",
			tmp=ReleaseHold[DeleteCases[tmp/.parallelized->"External",0]];
			If[Length[tmp]>0,{Plus@@(#[[1]]&/@tmp),#[[2]]&/@tmp},0]	
																							
		,
			tmp=ReleaseHold[DeleteCases[tmp/.parallelized->fdir,0]];
			(* gather the results *)
			tmp=If[Length[tmp]>0,{Plus@@(#[[1]]&/@tmp),#[[2]]&/@tmp},0];
			(*  the tmp has the form: { {the number of completed diagrams, the number of totall diagrams}, { {list of results for a certain type condensate}, {....} } } *)
			
			
			(*-----------------------------------*)
			(* if the evelation has already been done, the specified momentum may differ with the momentum in data *)
			If[MatchQ[tmp[[1]],{n_Integer,n_Integer} ],
			
			(* detect the momentum in the files *)
				tmpq=DeleteDuplicates[Cases[tmp,Pair[Momentum[mm_,___],Momentum[mm_,___]]:>mm,Infinity]];
				If[Length[tmpq]>1,Message[TetraquarkOPE::merr];Abort[]];
			
				If[MatchQ[tmp[[2]],{{0..}..}],
					0(* if the result vanish *)
				,
					tmpq=tmpq[[1]];
			(* if the momentum not match, replace it *)
					If[tmpq=!=qq&&qq=!="auto",Message[TetraquarkOPE::mwarn,tmpq,qq];tmp=tmp/.tmpq->qq]
				]
			];
			
			tmp
		]	
	
	
	,
	(* If not parallelized *)
	
		tmp=ReleaseHold[tmp/.parallelized->False];		
		Plus@@((QGather[(#/.{CA->SUNN,CF->(SUNN^2-1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0})&/@tmp)
	]

]
]


(* ::Code::Initialization::Plain:: *)
(*(*----------------------------------------------------------*)
	If[OptionValue[ShowProgress]===True,
	
		Needs["Parallel`Developer`"];
		PrintTemporary[Cases[tmp0,_EvaluationObject,Infinity]//StandardForm];
		(*evlist=Cases[tmp0,_EvaluationObject,Infinity];
		evlistlength=Length[evlist];
		
		
		progress:=Plus@@Boole[MatchQ[#,_dequeued]&/@(ProcessState[#]&/@evlist)];
		Dynamic[progress];(* dynamic update the number of completed processions *)
		
		
		(*--- Show pregress bar during the evaluation ---*)
		Monitor[tmp=WaitAll[tmp0],Row[{ProgressIndicator[progress,{0,evlistlength}],ToString[progress]<>"/"<>ToString[evlistlength]}," "]];*)
		
		tmp=WaitAll[tmp0];
		tmp=Plus@@#&/@tmp;
		(* separate the leading term to make the perturbative contribution looks clear *)
		
		con1 FourquarkLeading[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}}] + Plus@@((QGather[(#/.{CA->SUNN,CF->(SUNN^2-1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.{CA-2CF->1/CA, 2CA^3CF-CA^2(4CF^2+1)+1->0})&/@tmp)
	
	,
	
	
*)


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
