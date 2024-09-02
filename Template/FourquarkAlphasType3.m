(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
FourquarkAlphasType3::usage="FourquarkAlphasType3[q_,j1_,j2_] give the next leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents.
Here the gluon propagator comes from quark self-energy, the renormalization involve \!\(\*SubscriptBox[\(Z\), \(2\)]\) factor."

FourquarkAlphasType3::inderr="Dummy indices conflict!"

FourquarkAlphasType3::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`FourquarkAlphasType3`"]


(* ::Code::Initialization::Plain:: *)
Options[FourquarkAlphasType3] = {
	Parallelized->True,
	HoldFlavor->False,
	Renormalization->True,
	AutoNDR->False
}


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(* propagators *)

prop[q_]=I GSD[q]FAD[q];

(*xprop[x_]=FourierPX[prop[q],{q,x}];*)
xprop[x_]=1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];

			
(* 1-loop fermion propagator *)
(*aprop[q_]=(I gStrong)^2 IntegrateP[prop[q].SUNT[colzz].GAD[loraa].prop[q+ll].GAD[lorbb].SUNT[colzz].prop[q]propg[ll,loraa,lorbb],ll];*)
(*aprop[q_] = FAD[q]^2 GSD[q] qfact1[(qfact2[-I I^D 2^(-1-D) CF (-2+D)^2 gStrong^2 \[Pi]^(-D/2)] qGamma[2-D/2] qGamma[-1+D/2]^2)/qGamma[-1+D]] SPD[q,q]^(-1+D/2);*)
(*xaprop[x_]=FourierPX[aprop[q],{q,x}];*)
xaprop[x_] = GSD[x] qfact1[(qGamma[2-D/2] qGamma[-1+D/2]^2 (qfact2[-(1/128) I CF (-2+D)^2 (2+D) gStrong^2 (-\[Pi])^-D] qGamma[-2+D]
						+qfact2[1/64 I CF (-2+D)^2 gStrong^2 (-\[Pi])^-D] qGamma[-1+D]))/(qGamma[4-D/2] qGamma[-1+D])] SPD[x,x]^(2-D);



FourquarkAlphasType3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,atr,
hv21,hv22,listpole1,listpole2,holdf=OptionValue[HoldFlavor],diagrams,ndr=OptionValue[AutoNDR],ren=OptionValue[Renormalization],fdir,pall,pall1,pall2,dot},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[FourquarkAlphasType3::inderr];
	Abort[]
];


(*-------------------------------*)
(* B A^+ B *)
hv21=v3//ComplexConjugate;
hv22=v4//ComplexConjugate//FCI;
(* ComplexConjugate[sigma] will expand the DiracSigma to gamma matrices, recover it to DiracSigma *)
hv22=(hv22/.Dot->dot)/.f1_ dot[aa___,DiracGamma[LorentzIndex[lor1_,dim___],dim___],DiracGamma[LorentzIndex[lor2_,dim___],dim___],bb___]+f2_ dot[aa___,DiracGamma[LorentzIndex[lor2_,dim___],dim___],DiracGamma[LorentzIndex[lor1_,dim___],dim___],bb___]/;(f1+f2===0):>-2I f1 dot[aa,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],bb];
hv22=hv22/.{null->1,dot->Dot};
(*-------------------------------*)

If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];
(*---------------------------------------------------*)


If[OptionValue[AutoNDR]===True,atr=TR5,atr=TR];

(*---------------------------------------------------*)

If[pall===True,

	DistributeDefinitions[v1,v2,v3,v4,a1,b1,c1,d1,a2,b2,c2,d2];
	tmp= Plus@@WaitAll[Join[ ParallelSubmit[{hv21,hv22,holdf,atr},
									#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr]]&/@{type1,type2},{ParallelSubmit[{holdf,ren,ndr},reno[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},holdf,ndr,ren]]} ]];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,
	If[pall==="External",

	DistributeDefinitions[v1,v2,v3,v4,a1,b1,c1,d1,a2,b2,c2,d2];

	If[fdir=="None",

		Join[ ParallelSubmit[{hv21,hv22,holdf,atr},
									#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr]]&/@{type1,type2},{ParallelSubmit[{holdf,ren,ndr},reno[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},holdf,ndr,ren]]}]

	,



		pall1=ImExport[fdir,
						{"FourquarkAlphasType3_type1","FourquarkAlphasType3_type2"},
						{{hv21,hv22,holdf,atr},
						{qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr},
						{type1,type2}}
						];

		If[ren===True,
			pall2=ImExport[fdir,
						{"FourquarkAlphasType3_reno1"},
						{{holdf,ren,ndr},
						{qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},holdf,ndr,ren},
						{reno}}
						]
		,
			pall2={{1,1},{0}}				
		];
						
						
		{pall1[[1]]+pall2[[1]],Join[pall1[[2]],pall2[[2]]]}	


(*		imexport[fdir,
				{"fourquarkAlphasType3_type1","fourquarkAlphasType3_type2","fourquarkAlphasType3_reno1"},
				{{hv21,hv22,holdf,atr},
				{qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr},
				{type1,type2}}
				,
				{{holdf,ren},
				{qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},holdf,ren},
				{reno}}
				]*)
				
	](* !!! do not touch the output here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[...]. !!! *)
	
	,
	
		
		tmp=reno[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},holdf,ndr,ren]+ Plus@@Map[#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr]&,{type1,type2}];
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]

]

(*------------------------------------------------------------------*)



(*(*--- import & export ---*)

imexport[dir_,names_,{val_List,arg_List,types_List},{renval_,renarg_List,rens_List}]:=Block[{queue1,queue2,tmpnames,files,results,filt,tmp,fdir,vals,renvals,tmp2,tmp3,dirs1,dirs2},	

	
files=FileNames["*.txt",dir];
tmpnames=StringSplit[StringSplit[#,"\\"|"/"][[-1]],"."][[1]]&/@files;(* get the file names *)
tmpnames={tmpnames,files}//Transpose;

(* import the already exist results *)
filt=Boole[!FileExistsQ[FileNameJoin[{dir,#<>".txt"}]]]&/@names;
tmp=DeleteCases[(1-filt) names,0];
tmp=(Plus@@Cases[tmpnames,{#,file_}:>ToExpression[Import[file]]])&/@tmp;


(*----------------------------*)
If[(Plus@@filt)==0,
(* if all results are already exist *)
       tmp
,

(*--- the experssions should be evaluated ---*)


(* the directory and file names for the expersssion which should be evaulated *)
(* names ~ Join[types,rens] *)
	dirs1=FileNameJoin[{dir,#<>".txt"}]&/@names[[;;Length[types]]];
	dirs2=FileNameJoin[{dir,#<>".txt"}]&/@names[[Length[types]+1;;]];
	
	dirs1=DeleteCases[filt[[;;Length[types]]]dirs1,0];
	dirs2=DeleteCases[filt[[Length[types]+1;;]]dirs2,0];
	
(*(* The ParallelSubmit is trick; here defines a directory function so that the directory can be distributed before submiting *)
	Do[dirs1[types[[i]]]=queue1[[i]],{i,1,Length[queue1]}];
	Do[dirs2[rens[[i]]]=queue2[[i]],{i,1,Length[queue2]}];*)
	vals=Append[val,fdir];
	renvals=Append[renval,fdir];


(* the experssion should be evaluated *)
	queue1=DeleteCases[filt[[;;Length[types]]]types,0];
	queue2=DeleteCases[filt[[Length[types]+1;;]]rens,0];


(* export the nonexist results(ParallelSubmit the experssions) *)
	DistributeDefinitions[#]&/@Join[val,renval];



	tmp2=Table[ParallelSubmit[{dirs1,i,queue1},Export[dirs1[[i]],FCI[queue1[[i]]@@arg]]],{i,1,Length[queue1]}];
	tmp3=Table[ParallelSubmit[{dirs2,i,queue2},Export[dirs2[[i]],FCI[queue2[[i]]@@renarg]]],{i,1,Length[queue2]}];

(* Join[ParallelSubmit[...]&/@...,ParallelSubmit[...]&/@... doesn't work, it maps ParallelSubmit[] but doesn't use the definition of queue *)
	Join[tmp2,tmp3]

]
]*)


(* ::Code::Initialization::Plain:: *)
(* renormalization *)
reno[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,False]=0
reno[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},holdf_,ndr_,True]:=(CF gStrong^2/(Epsilon 4 Pi^2)FourquarkLeading[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},EpsOrder->1,HoldFlavor->holdf,AutoNDR->ndr])/.{CA-2CF->1/CA,2CA CF-CA^2+1->0,-2CA CF+CA^2-1->0}/.CF->(CA^2-1)/(2CA)


(* ::Code::Initialization::Plain:: *)
(*--- The following are generated by algorithm ---*)


(* ::Input::Initialization:: *)
type1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1},


dia1=Contract[atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . xaprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . xaprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . xprop[-x] . v2 . xaprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . xaprop[x] . hv3 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] - Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . xaprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] - Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . xprop[-x] . v2 . xaprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv4 . xaprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv4 . xprop[-x] . v2 . xaprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia1=FourierXP[dia1,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
type2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2},


dia2=Contract[atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv4 . xaprop[-x] . v1 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Contract[atr[SUNTrace[hv3 . xaprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Contract[atr[SUNTrace[hv3 . xprop[-x] . v2 . xaprop[x]]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Contract[atr[SUNTrace[hv3 . xaprop[-x] . v1 . xprop[x] . hv4 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xaprop[x] . hv4 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Contract[atr[SUNTrace[hv4 . xaprop[-x] . v1 . xprop[x] . hv3 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Contract[atr[SUNTrace[hv3 . xaprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xaprop[x]]]*atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia2=FourierXP[dia2,{x,q}]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
