(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
TetraquarkHybridlike::usage="TetraquarkHybridlike[q_,j1_,j2_] gives leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents. "

TetraquarkHybridlike::inderr="Dummy indices conflict!"

TetraquarkHybridlike::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`TetraquarkHybridlike`"]

Options[TetraquarkHybridlike] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True,
	EpsOrder->0,
	Strategy->"Fourier",
	Pole->1
}


TetraquarkHybridlike[qq_,0,_,OptionsPattern[]]=0
TetraquarkHybridlike[qq_,_,0,OptionsPattern[]]=0


TetraquarkHybridlike[qq_,jh_List,cur2_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmp,files,evaluations,fdir,pall,holdf=OptionValue[HoldFlavor],ndr=OptionValue[AutoNDR],eps=OptionValue[EpsOrder]},



(*-------------------------------*)
If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];


(*-------------------------------*)

If[pall===True,

	DistributeDefinitions[jh,cur2];
	
	tmp= Plus@@WaitAll[ParallelSubmit[{holdf,ndr,eps},TetraquarkHybridlike1[qq,#[[2]],cur2,
															Pole->#[[1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps
															]
									]&/@jh];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

,
	If[pall==="External",

		
		DistributeDefinitions[jh,cur2];
		
(*		files=Table["TetraquarkHybridlike1_"<>ToString[i],{i,1,Length[jh]}];			
		evaluations=Transpose[{files,jh}];
		DistributeDefinitions[evaluations];*)

		If[fdir==="None",

			ParallelSubmit[{holdf,ndr,eps},
				TetraquarkHybridlike1[qq,#[[2]],cur2,Pole->#[[1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps]
			]&/@jh

		,
		(* export the results for each type *)
			files=Table["TetraquarkHybridlike1_"<>ToString[i],{i,1,Length[jh]}];
			
			evaluations=Transpose[{files,jh}];
			
			pall=ImExport[fdir,
						{#[[1]]},
						{{holdf,ndr,eps},
						{qq,#[[2,2]],cur2,Pole->#[[2,1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps},
						{TetraquarkHybridlike1}}
						]&/@evaluations;
						
			pall=Transpose[pall];
			pall={Plus@@pall[[1]],Plus@@pall[[2]]}

				
		]
		(* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
	
		,
	
		
		tmp= Plus@@(TetraquarkHybridlike1[qq,#[[2]],cur2,Pole->#[[1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps
															]&/@jh);
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]

]/;And@@(MatchQ[#,{_,{_,_,_,_,_,_,_}}]&/@jh)


(*-------------------------------------------*)


TetraquarkHybridlike[qq_,cur1_FourquarkCurrent,jh_List,OptionsPattern[]]:=Block[{null,tmp,files,evaluations,fdir,pall,holdf=OptionValue[HoldFlavor],ndr=OptionValue[AutoNDR],eps=OptionValue[EpsOrder]},


(*-------------------------------*)
If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];


(*-------------------------------*)

If[pall===True,

	DistributeDefinitions[jh,cur1];
	
	tmp= Plus@@WaitAll[ParallelSubmit[{holdf,ndr,eps},TetraquarkHybridlike2[qq,cur1,#[[2]],
															Pole->#[[1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps
															]
									]&/@jh];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,
	If[pall==="External",

		
		DistributeDefinitions[jh,cur1];
		
(*		files=Table["TetraquarkHybridlike1_"<>ToString[i],{i,1,Length[jh]}];			
		evaluations=Transpose[{files,jh}];
		DistributeDefinitions[evaluations];*)

		If[fdir==="None",

			ParallelSubmit[{holdf,ndr,eps},
				TetraquarkHybridlike2[qq,cur1,#[[2]],Pole->#[[1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps]
			]&/@jh

		,
		(* export the results for each type *)
			files=Table["TetraquarkHybridlike2_"<>ToString[i],{i,1,Length[jh]}];
			
			evaluations=Transpose[{files,jh}];
			
			pall=ImExport[fdir,
						{#[[1]]},
						{{holdf,ndr,eps},
						{qq,cur1,#[[2,2]],Pole->#[[2,1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps},
						{TetraquarkHybridlike2}}
						]&/@evaluations;
						
			pall=Transpose[pall];
			pall={Plus@@pall[[1]],Plus@@pall[[2]]}

				
		]
		(* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
	
		,
	
		
		tmp= Plus@@(TetraquarkHybridlike2[qq,cur1,#[[2]],Pole->#[[1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps
															]&/@jh);
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]

]/;And@@(MatchQ[#,{_,{_,_,_,_,_,_,_}}]&/@jh)


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
