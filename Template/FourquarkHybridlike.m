(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


FourquarkHybridlike::usage="FourquarkHybridlike[qq_,J1_hybirdlikeList,J2_FourquarkCurrent] gives the perturbative contribution of <J1 J2>, with one fourquark current and one Hybirdlike current: J1 and J2 =\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a1\)]\)\!\(\*SubscriptBox[\(v1\[CapitalPsi]\), \(b1\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c1\)]\)\!\(\*SubscriptBox[\(v2\[CapitalPsi]\), \(d1\)]\) "

(*FourquarkHybridlike::inderr="Dummy indices conflict!"

FourquarkHybridlike::gerr="No gluon has found!"*)


Begin["`Private`FourquarkHybridlike`"]

Options[FourquarkHybridlike] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->False,
	EpsOrder->0,
	Pole->1
}



(* The fourquark current mix with Hybird-like current under renormalization, as long as the symmetry permit. *)



FourquarkHybridlike[qq_,_,0,OptionsPattern[]]=0
FourquarkHybridlike[qq_,0,_,OptionsPattern[]]=0


(*--- < hybird_like  fourquark >  ---*)
FourquarkHybridlike[qq_,jh_List,jf_FourquarkCurrent,OptionsPattern[]]:=Block[{null,tmp,v3,v4,a2,b2,c2,d2,files,evaluations,fdir,pall,holdf=OptionValue[HoldFlavor],ndr=OptionValue[AutoNDR],eps=OptionValue[EpsOrder]},


(*-------------------------------*)
If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];

(*-------------------------------*)

{a2,v3,b2}={jf[[1,1]],jf[[1,2]],jf[[1,3]]};
{c2,v4,d2}={jf[[2,1]],jf[[2,2]],jf[[2,3]]};


(*-------------------------------*)

If[pall===True,

	DistributeDefinitions[jh];
	
	tmp= Plus@@WaitAll[ParallelSubmit[{v3,v4,a2,b2,c2,d2,holdf,ndr,eps},FourquarkHybridlike1[qq,#[[2]],{v3,v4,{a2,b2,c2,d2}},
															Pole->#[[1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps
															]
									]&/@jh];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,
	If[pall==="External",

		
		DistributeDefinitions[jh];
		
(*		files=Table["FourquarkHybridlike1_"<>ToString[i],{i,1,Length[jh]}];			
		evaluations=Transpose[{files,jh}];
		DistributeDefinitions[evaluations];*)

		If[fdir==="None",

			ParallelSubmit[{v3,v4,a2,b2,c2,d2,holdf,ndr,eps},
				FourquarkHybridlike1[qq,#[[2]],{v3,v4,{a2,b2,c2,d2}},Pole->#[[1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps]
			]&/@jh

		,
		(* export the results for each type *)
			files=Table["FourquarkHybridlike1_"<>ToString[i],{i,1,Length[jh]}];
			
			evaluations=Transpose[{files,jh}];
			
			pall=ImExport[fdir,
						{#[[1]]},
						{{v3,v4,a2,b2,c2,d2,holdf,ndr,eps},
						{qq,#[[2,2]],{v3,v4,{a2,b2,c2,d2}},Pole->#[[2,1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps},
						{FourquarkHybridlike1}}
						]&/@evaluations;
						
			pall=Transpose[pall];
			pall={Plus@@pall[[1]],Plus@@pall[[2]]}

				
		]
		(* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
	
		,
	
		
		tmp= Plus@@(FourquarkHybridlike1[qq,#[[2]],{v3,v4,{a2,b2,c2,d2}},
															Pole->#[[1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps
															]&/@jh);
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]

(*Plus@@(FourquarkHybridlike1[qq,#[[2]],{v3,v4,{a2,b2,c2,d2}},Pole->#[[1]],Parallelized\[Rule]OptionValue[Parallelized],HoldFlavor\[Rule]OptionValue[HoldFlavor],AutoNDR\[Rule]OptionValue[AutoNDR],EpsOrder\[Rule]OptionValue[EpsOrder]]&/@jh)

*)

]/;And@@(MatchQ[#,{_,{_,_,_,_,_}}]&/@jh)




(*--- < fourquark   hybird_like >  ---*)
FourquarkHybridlike[qq_,jf_FourquarkCurrent,jh_List,OptionsPattern[]]:=Block[{null,tmp,v1,v2,a1,b1,c1,d1,files,evaluations,fdir,pall,holdf=OptionValue[HoldFlavor],ndr=OptionValue[AutoNDR],eps=OptionValue[EpsOrder]},


(*-------------------------------*)
If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];
(*-------------------------------*)



{a1,v1,b1}={jf[[1,1]],jf[[1,2]],jf[[1,3]]};
{c1,v2,d1}={jf[[2,1]],jf[[2,2]],jf[[2,3]]};


(*-------------------------------*)

If[pall===True,

	DistributeDefinitions[jh];
	
	tmp= Plus@@WaitAll[ParallelSubmit[{v1,v2,a1,b1,c1,d1,holdf,ndr,eps},FourquarkHybridlike2[qq,{v1,v2,{a1,b1,c1,d1}},#[[2]],
															Pole->#[[1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps
															]
									]&/@jh];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,
	If[pall==="External",

		
		DistributeDefinitions[jh];

		If[fdir==="None",

			ParallelSubmit[{v1,v2,a1,b1,c1,d1,holdf,ndr,eps},
				FourquarkHybridlike2[qq,{v1,v2,{a1,b1,c1,d1}},#[[2]],Pole->#[[1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps]
			]&/@jh

		,
		(* export the results for each type *)
			files=Table["FourquarkHybridlike2_"<>ToString[i],{i,1,Length[jh]}];
			
			evaluations=Transpose[{files,jh}];
			
			pall=ImExport[fdir,
						{#[[1]]},
						{{v1,v2,a1,b1,c1,d1,holdf,ndr,eps},
						{qq,{v1,v2,{a1,b1,c1,d1}},#[[2,2]],Pole->#[[2,1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps},
						{FourquarkHybridlike2}}
						]&/@evaluations;
						
			pall=Transpose[pall];
			pall={Plus@@pall[[1]],Plus@@pall[[2]]}			

				
		]
		(* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
	
		,
	
		
		tmp= Plus@@(FourquarkHybridlike2[qq,{v1,v2,{a1,b1,c1,d1}},#[[2]],
															Pole->#[[1]],Parallelized->False,HoldFlavor->holdf,AutoNDR->ndr,EpsOrder->eps
															]&/@jh);
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]

(*Plus@@(FourquarkHybridlike1[qq,{v1,v2,{a1,b1,c1,d1}},#[[2]],Pole->#[[1]],Parallelized\[Rule]OptionValue[Parallelized],HoldFlavor\[Rule]OptionValue[HoldFlavor],AutoNDR\[Rule]OptionValue[AutoNDR],EpsOrder\[Rule]OptionValue[EpsOrder]]&/@jh)
*)
]/;And@@(MatchQ[#,{_,{_,_,_,_,_}}]&/@jh)


End[]
(*EndPackage[]*)
