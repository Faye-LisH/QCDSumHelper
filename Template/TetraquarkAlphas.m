(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


TetraquarkAlphas::usage="TetraquarkAlphas[q_,j1_,j2_] give the next leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents."

TetraquarkAlphas::inderr="Dummy indices conflict!"

TetraquarkAlphas::curerr="Unknow current structure!"

TetraquarkAlphas::terr="NLO diagrams types unrecoginzed."


Begin["`Private`TetraquarkAlphas`"]


Options[TetraquarkAlphas] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->True,
	Renormalization->True,
	Strategy->"Fourier",
	AlphasTypes->"All"
}


TetraquarkAlphas[qq_,current1_,current2_,ops:OptionsPattern[]]:=Block[{tmp1,tmp2,tmp3,leading,holdf=OptionValue[HoldFlavor],
ndr=OptionValue[AutoNDR],ren=OptionValue[Renormalization],pall,fdir,strategy=OptionValue[Strategy],pall0,alphas,type1,type2,type3},

alphas=OptionValue[AlphasTypes];
If[StringQ[alphas],alphas={alphas}];

If[Head[alphas]===List,
	If[!FreeQ[alphas,"All"],
		{type1,type2,type3}={1,1,1}
	,
		If[And@@(FreeQ[{"d0_1","d0_2","d0_3"}, #]&/@alphas), Message[TetraquarkAlphas::terr];Abort[]];
		{type1,type2,type3}=Boole[!FreeQ[alphas,#]]&/@{"d0_1","d0_2","d0_3"}
	]
,
	Message[TetraquarkAlphas::terr];
	Abort[]
];

(*--------------------------*)
tmp1=If[type1==1,
		TetraquarkAlphasType1[qq,current1,current2,Parallelized->OptionValue[Parallelized],HoldFlavor->holdf,AutoNDR->ndr,Renormalization->ren,Strategy->strategy]
	,
		{0,{0}}
	];
					

tmp2=If[type2==1,
		TetraquarkAlphasType2[qq,current1,current2,Parallelized->OptionValue[Parallelized],HoldFlavor->holdf,AutoNDR->ndr,Renormalization->ren,Strategy->strategy]
	,
		{0,{0}}
	];

					
tmp3=If[type3==1,
		TetraquarkAlphasType3[qq,current1,current2,Parallelized->OptionValue[Parallelized],HoldFlavor->holdf,AutoNDR->ndr,Renormalization->ren]
	,
		{0,{0}}
	];

(*--------------------------*)
										
If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];	

							
(*--------------------------------------------*)													
																												
If[pall==="External",
	
	If[fdir=="None",
	
		Join[tmp1,tmp2,tmp3]
	,
		
		{tmp1[[1]]+tmp2[[1]]+tmp3[[1]],Join[tmp1[[2]], tmp2[[2]], tmp3[[2]]]}(* avoid apply any function here, since the expression haven't been evaluated in EvaluationObject[...] *)
	]
	
,(*--------------------------------------------*)

	
	QGather[tmp1+tmp2+tmp3,qq,ShowasTable->False]
		

]
]


End[]
(*EndPackage[]*)
