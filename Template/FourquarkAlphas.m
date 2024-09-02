(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


FourquarkAlphas::usage="FourquarkAlphas[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}}] give the next leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents.
Here the gluon propagator connect with two quark propagators, the corresponding renormalization involve tetraquark currents with same flavor structure."

FourquarkAlphas::inderr="Dummy indices conflict!"

FourquarkAlphas::curerr="Unknow current structure!"

FourquarkAlphas::terr="NLO diagrams types unrecoginzed."


Begin["`Private`FourquarkAlphas`"]


Options[FourquarkAlphas] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->False,
	Renormalization->True,
	Strategy->"Fourier",
	AlphasTypes->"All"
	(*,
	ToFile->"None"*)
}


FourquarkAlphas[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ops:OptionsPattern[]]:=Block[{tmp1,tmp2,tmp3,holdf=OptionValue[HoldFlavor],
ndr=OptionValue[AutoNDR],ren=OptionValue[Renormalization],renh,renz,pall,fdir,strategy=OptionValue[Strategy],pall0,alphas,type1,type2,type3},

Which[
	ren==="ExceptFourquark"
	,
	ren=False;
	renh=True;
	renz=True
,
	ren==="ExceptHybirdLike"
	,
	ren=True;
	renh=False;
	renz=True
,
	ren==="FermionOnly"
	,
	ren=False;
	renh=False;
	renz=True	
,
	True
	,
	renh=ren;
	renz=ren	
];


(*--------------------------*)
alphas=OptionValue[AlphasTypes];

If[StringQ[alphas],alphas={alphas}];


If[Head[alphas]===List,
	If[!FreeQ[alphas,"All"],
		{type1,type2,type3}={1,1,1}
	,
		If[And@@(FreeQ[{"d0_1","d0_2","d0_3"}, #]&/@alphas), Message[FourquarkAlphas::terr];Abort[]];
		{type1,type2,type3}=Boole[!FreeQ[alphas,#]]&/@{"d0_1","d0_2","d0_3"}
	]
,
	Message[TetraquarkAlphas::terr];
	Abort[]
];

(*--------------------------*)

tmp1=If[type1==1,
		FourquarkAlphasType1[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->OptionValue[Parallelized],HoldFlavor->holdf,AutoNDR->ndr,Renormalization->ren,Strategy->strategy]
	,
		{0,{0}}
	];
					

tmp2=If[type2==1,
		FourquarkAlphasType2[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->OptionValue[Parallelized],HoldFlavor->holdf,AutoNDR->ndr,Renormalization->renh,Strategy->strategy]
	,
		{0,{0}}
	];

					
tmp3=If[type3==1,
		FourquarkAlphasType3[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},Parallelized->OptionValue[Parallelized],HoldFlavor->holdf,AutoNDR->ndr,Renormalization->renz]
	,
		{0,{0}}
	];

										
(*--------------------------------------------*)		

					
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
		
		{tmp1[[1]]+tmp2[[1]]+tmp3[[1]],Join[ tmp1[[2]], tmp2[[2]], tmp3[[2]]]}(* avoid apply any function here, since the expression haven't been evaluated in EvaluationObject[...] *)
	]
	
,(*--------------------------------------------*)

	
	QGather[tmp1+tmp2+tmp3,qq,ShowasTable->False]
		

]

]


End[]
(*EndPackage[]*)
