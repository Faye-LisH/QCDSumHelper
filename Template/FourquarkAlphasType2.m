(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
FourquarkAlphasType2::usage="FourquarkAlphasType2[q_,j1_,j2_] give the next leading order perturbative contribution of <j1 j2>, with j1, j2 are tetraquark currents.
 Here the gluon propagator connect with j1 & j2, the corresponding renormalization invole hybird-like current."

FourquarkAlphasType2::inderr="Dummy indices conflict!"

FourquarkAlphasType2::curerr="Unknow current structure!"


(* ::Code::Initialization::Plain:: *)
Begin["`Private`FourquarkAlphasType2`"]


(* ::Code::Initialization::Plain:: *)
Options[FourquarkAlphasType2] = {
	Parallelized->True,
	HoldFlavor->False,
	AutoNDR->False,
	Renormalization->True,
	Strategy->"Fourier"
}


(* ::Code::Initialization::Plain:: *)
(*------------------------------------------------------------------*)
(* propagators *)

prop[q_]=I GSD[q]FAD[q];

(*xprop[x_]=FourierPX[prop[q],{q,x}];*)
xprop[x_]=1/2 I^(1-D) \[Pi]^(-D/2) DiracGamma[Momentum[x,D],D] Pair[Momentum[x,D],Momentum[x,D]]^(-D/2) qGamma[D/2];
		
			
xprop2[x_,lor1_,lor2_,lor3_,lor4_]=Pair[LorentzIndex[lor1, D], LorentzIndex[lor2, D]]*Pair[LorentzIndex[lor3, D], LorentzIndex[lor4, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(3 - (3*D)/2)*
  qfact1[(qfact2[((20 - 8*D + D^2)*gStrong^2)/(1024*E^(((3*I)/2)*D*Pi)*Pi^((3*D)/2))]*qGamma[1 - D/2]^2*qGamma[D/2]^4*qGamma[-3 + (3*D)/2])/(qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lor4, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lor3, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(3 - (3*D)/2)*
  qfact1[(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-3 + (3*D)/2])/(256*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lor3, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lor4, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(3 - (3*D)/2)*
  qfact1[(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-3 + (3*D)/2])/(256*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lor3, D], LorentzIndex[lor4, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[(qfact2[((12 - 8*D + D^2)*gStrong^2)/(512*E^(((3*I)/2)*D*Pi)*Pi^((3*D)/2))]*qGamma[1 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lor2, D]]*Pair[LorentzIndex[lor3, D], Momentum[x, D]]*Pair[LorentzIndex[lor4, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[(qfact2[((12 - 8*D + D^2)*gStrong^2)/(512*E^(((3*I)/2)*D*Pi)*Pi^((3*D)/2))]*qGamma[1 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lor4, D]]*Pair[LorentzIndex[lor3, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[-(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(128*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lor4, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lor3, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[-(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(128*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lor3, D]]*Pair[LorentzIndex[lor4, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[-(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(128*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lor3, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lor4, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[-(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(128*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lor3, D], Momentum[x, D]]*Pair[LorentzIndex[lor4, D], Momentum[x, D]]*
  Pair[Momentum[x, D], Momentum[x, D]]^(1 - (3*D)/2)*qfact1[(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-1 + (3*D)/2])/(64*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)];


(* ::Code::Initialization::Plain:: *)
xprop2a[x_,lor1_,lor2_,lor3_,lor4_]=Pair[LorentzIndex[lor1, D], LorentzIndex[lor2, D]]*Pair[LorentzIndex[lor3, D], LorentzIndex[lor4, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(3 - (3*D)/2)*
 qfact1[(qfact2[((20 - 8*D + D^2)*gStrong^2)/(1024*E^(((3*I)/2)*D*Pi)*Pi^((3*D)/2))]*qGamma[1 - D/2]^2*qGamma[D/2]^4*qGamma[-3 + (3*D)/2])/(qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], LorentzIndex[lor4, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lor3, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(3 - (3*D)/2)*
  qfact1[(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-3 + (3*D)/2])/(256*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] ;


(* ::Code::Initialization::Plain:: *)
xprop2b[x_,lor1_,lor2_,lor3_,lor4_]=Pair[LorentzIndex[lor1, D], LorentzIndex[lor3, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lor4, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(3 - (3*D)/2)*
  qfact1[(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-3 + (3*D)/2])/(256*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lor3, D], LorentzIndex[lor4, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[(qfact2[((12 - 8*D + D^2)*gStrong^2)/(512*E^(((3*I)/2)*D*Pi)*Pi^((3*D)/2))]*qGamma[1 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(qGamma[5 - D]*qGamma[D]^2)] ;


(* ::Code::Initialization::Plain:: *)
xprop2c[x_,lor1_,lor2_,lor3_,lor4_]= Pair[LorentzIndex[lor1, D], LorentzIndex[lor2, D]]*Pair[LorentzIndex[lor3, D], Momentum[x, D]]*Pair[LorentzIndex[lor4, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[(qfact2[((12 - 8*D + D^2)*gStrong^2)/(512*E^(((3*I)/2)*D*Pi)*Pi^((3*D)/2))]*qGamma[1 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lor4, D]]*Pair[LorentzIndex[lor3, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[-(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(128*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)];


(* ::Code::Initialization::Plain:: *)
xprop2d[x_,lor1_,lor2_,lor3_,lor4_]=Pair[LorentzIndex[lor1, D], LorentzIndex[lor4, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lor3, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[-(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(128*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], LorentzIndex[lor3, D]]*Pair[LorentzIndex[lor4, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[-(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(128*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)];


(* ::Code::Initialization::Plain:: *)
xprop2e[x_,lor1_,lor2_,lor3_,lor4_]=Pair[LorentzIndex[lor1, D], LorentzIndex[lor3, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lor4, D], Momentum[x, D]]*Pair[Momentum[x, D], Momentum[x, D]]^(2 - (3*D)/2)*
  qfact1[-(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-2 + (3*D)/2])/(128*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)] + 
 Pair[LorentzIndex[lor1, D], Momentum[x, D]]*Pair[LorentzIndex[lor2, D], Momentum[x, D]]*Pair[LorentzIndex[lor3, D], Momentum[x, D]]*Pair[LorentzIndex[lor4, D], Momentum[x, D]]*
  Pair[Momentum[x, D], Momentum[x, D]]^(1 - (3*D)/2)*qfact1[(I^(2*D)*gStrong^2*qGamma[2 - D/2]^2*qGamma[D/2]^4*qGamma[-1 + (3*D)/2])/(64*E^(((5*I)/2)*D*Pi)*Pi^((3*D)/2)*qGamma[5 - D]*qGamma[D]^2)];


FourquarkAlphasType2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,atr,
f11,f12,f13,f14,f21,f22,f23,f24,hv21,hv22,listpole1,listpole2,holdf=OptionValue[HoldFlavor],strategy=OptionValue[Strategy],diagrams,ndr=OptionValue[AutoNDR],ren=OptionValue[Renormalization],fdir,pall,files1,pall1,files2,pall2,dot},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[FourquarkAlphasType2::inderr];
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


If[ToLowerCase[ToString[strategy]]=="fourier",
	diagrams={xtype1a,xtype1b,xtype1c,xtype1d,xtype1e,xtype2a,xtype2b,xtype2c,xtype2d,xtype2e,xtype3a,xtype3b,xtype3c,xtype3d,xtype3e,xtype4a,xtype4b,xtype4c,xtype4d,xtype4e,xtype5a,xtype5b,xtype5c,xtype5d,xtype5e,xtype6a,xtype6b,xtype6c,xtype6d,xtype6e,xtype7a,xtype7b,xtype7c,xtype7d,xtype7e,xtype8a,xtype8b,xtype8c,xtype8d,xtype8e,xtype9a,xtype9b,xtype9c,xtype9d,xtype9e,xtype10a,xtype10b,xtype10c,xtype10d,xtype10e,xtype11a,xtype11b,xtype11c,xtype11d,xtype11e,xtype12a,xtype12b,xtype12c,xtype12d,xtype12e,xtype13a,xtype13b,xtype13c,xtype13d,xtype13e,xtype14a,xtype14b,xtype14c,xtype14d,xtype14e,xtype15a,xtype15b,xtype15c,xtype15d,xtype15e,xtype16a,xtype16b,xtype16c,xtype16d,xtype16e}
	(*diagrams={xtype1,xtype2,xtype3,xtype4,xtype5,xtype6,xtype7,xtype8,xtype9,xtype10,xtype11,xtype12,xtype13,xtype14,xtype15,xtype16}*)
,
	If[ToLowerCase[ToString[strategy]]=="integratefirst",
		diagrams={iptype1,iptype2,iptype3,iptype4,iptype5,iptype6,iptype7,iptype8,iptype9,iptype10,iptype11,iptype12,iptype13,iptype14,iptype15,iptype16}
	,
		diagrams={ptype1,ptype2,ptype3,ptype4,ptype5,ptype6,ptype7,ptype8,ptype9,ptype10,ptype11,ptype12,ptype13,ptype14,ptype15,ptype16}
	]
];


(*---------------------------------------------------*)


If[OptionValue[AutoNDR]===True,atr=TR5,atr=TR];

(*---------------------------------------------------*)
If[pall===True,

	DistributeDefinitions[v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
	
	tmp= Plus@@WaitAll[Join[ParallelSubmit[{hv21,hv22,holdf,atr},
									#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr]]&/@diagrams,ParallelSubmit[{ndr,holdf,ren},-#[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},ndr,holdf,ren]]&/@{reno1,reno2}]];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,
	If[pall==="External",


		DistributeDefinitions[v1,v2,v3,v4,a1,b1,c1,d1,a2,b2,c2,d2];

		If[fdir==="None",

			Join[ParallelSubmit[{hv21,hv22,holdf,atr},
									#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr]]&/@diagrams,
				ParallelSubmit[{ndr,holdf,ren},-#[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},ndr,holdf,ren]]&/@{reno1,reno2}
				]

		,
		(* export the results for each type *)
			files1=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files1=("FourquarkAlphasType2_"<>#)&/@files1;
			
			
			pall1=ImExport[fdir,
						files1,
						{{hv21,hv22,holdf,atr},
						{qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr},
						diagrams}
						];
		
		
			files2=(StringSplit[ToString[#],"`"][[-1]])&/@{reno1,reno2};
			files2=("FourquarkAlphasType2_"<>#)&/@files2;
			
			If[ren===True,
				pall2=ImExport[fdir,
						files2,
						{{ndr,holdf,ren},
						{qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},ndr,holdf,ren},
						{reno1,reno2}}
						]
			,
				pall2={{2,2},{0,0}}			
			];


			{pall1[[1]]+pall2[[1]],Join[pall1[[2]],-pall2[[2]]]}	
		
		
		
		
		
(*		files=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files=("fourquarkAlphasType2_"<>#)&/@files;
			files=Append[files,"fourquarkAlphasType2_reno1"];
			
			
			imexport[fdir,
				files,
				{{hv21,hv22,holdf,atr},
				{qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr},
				diagrams}
				,
				{{ndr,holdf,ren},
				{qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},ndr,holdf,ren},
				{reno}}
				]*)
				
		]
		(* !!! do not touch the tmp here, any function may modifiy the things in EvaluationObject[...] should be applied outer the WaitAll[tmp]. !!! *)
	
		,
	
		
		tmp=Plus@@( -#[qq,{v1,v2,{a1,b1,c1,d1}},{v3,v4,{a2,b2,c2,d2}},ndr,holdf,ren]&/@{reno1,reno2})+ Plus@@Map[#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr]&,diagrams];
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]

]

(*------------------------------------------------------------------*)


(*	
reno[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ndr_,holdf_,False]:=0

		
reno[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ndr_,holdf_,True]:=Block[{listpole1,listpole2},
(* it's possible to generate terms like Eps[u,v,a,b]/(eps Pi^2)..., which give factor D after contarction *)
listpole1=Plus@@(FourquarkHybridlike1[qq,#[[2]],{v3,v4,{a2,b2,c2,d2}},Pole->#[[1]],EpsOrder->0,AutoNDR->ndr,HoldFlavor->holdf,Parallelized->False]&/@FourquarkHPole[{v1,{a1,b1}},{v2,{c1,d1}}]);
listpole2=Plus@@(FourquarkHybridlike2[qq,{v1,v2,{a1,b1,c1,d1}},#[[2]],Pole->#[[1]],EpsOrder->0,AutoNDR->ndr,HoldFlavor->holdf,Parallelized->False]&/@FourquarkHPole[{v3,{a2,b2}},{v4,{c2,d2}}]);
	
QGather[listpole1+listpole2/.{CA-2CF->1/CA,2CA CF-CA^2+1->0}/.CF->(CA^2-1)/(2CA),qq,ShowasTable->False]


]
*)


(* ::Code::Initialization::Plain:: *)
reno1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ndr_,holdf_,False]=0
reno1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ndr_,holdf_,True]:=Block[{listpole1,listpole2,poles},

poles = FourquarkHPole[{v1,{a1,b1}},{v2,{c1,d1}}];
If[poles===0,
	0
,
	(Plus@@(FourquarkHybridlike1[qq,#[[2]],{v3,v4,{a2,b2,c2,d2}},Pole->#[[1]],EpsOrder->0,AutoNDR->ndr,HoldFlavor->holdf,Parallelized->False]&/@poles))/.{CA-2CF->1/CA,2CA CF-CA^2+1->0}/.CF->(CA^2-1)/(2CA)
]
]

reno2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ndr_,holdf_,False]=0
reno2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ndr_,holdf_,True]:=Block[{listpole1,listpole2,poles},

poles = FourquarkHPole[{v3,{a2,b2}},{v4,{c2,d2}}];
If[poles===0,
	0
,
	(Plus@@(FourquarkHybridlike2[qq,{v1,v2,{a1,b1,c1,d1}},#[[2]],Pole->#[[1]],EpsOrder->0,AutoNDR->ndr,HoldFlavor->holdf,Parallelized->False]&/@poles))/.{CA-2CF->1/CA,2CA CF-CA^2+1->0}/.CF->(CA^2-1)/(2CA)
]
]


(* ::Code::Initialization::Plain:: *)
(* need not to be parallelized, since the evaluation of 3-loop diagram is much faster than 4-loop diagram *)
(*reno[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},ndr_,holdf_,"External",True]:=Block[{listpole1,listpole2},

listpole1=(#[[1]]/.FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf])FourquarkHybridlike1[qq,#[[2]],{v3,v4,{a2,b2,c2,d2}},EpsOrder->0,AutoNDR->ndr,HoldFlavor->holdf,Parallelized->"External"]&/@FourquarkHPole[{v1,{a1,b1}},{v2,{c1,d1}}];
listpole2=(#[[1]]/.FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf])FourquarkHybridlike2[qq,{v1,v2,{a1,b1,c1,d1}},#[[2]],EpsOrder->0,AutoNDR->ndr,HoldFlavor->holdf,Parallelized->"External"]&/@FourquarkHPole[{v3,{a2,b2}},{v4,{c2,d2}}];

		
Plus@@Join[listpole1,listpole2]

]
*)


(* ::Code::Initialization::Plain:: *)
(* The following are generated by algorithem *)


(* ::Input::Initialization:: *)
xtype1a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia1a=-(Contract[atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1]);


dia1a=FourierXP[dia1a,{x,q}]/.SUNN->CA;

dia1a=QEvaluate[I ScaleMu^(4(4-D)) dia1a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype1b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia1b=-(Contract[atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1]);


dia1b=FourierXP[dia1b,{x,q}]/.SUNN->CA;

dia1b=QEvaluate[I ScaleMu^(4(4-D)) dia1b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype1c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia1c=-(Contract[atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1]);


dia1c=FourierXP[dia1c,{x,q}]/.SUNN->CA;

dia1c=QEvaluate[I ScaleMu^(4(4-D)) dia1c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype1d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia1d=-(Contract[atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1]);


dia1d=FourierXP[dia1d,{x,q}]/.SUNN->CA;

dia1d=QEvaluate[I ScaleMu^(4(4-D)) dia1d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype1e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia1e=-(Contract[atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1]);


dia1e=FourierXP[dia1e,{x,q}]/.SUNN->CA;

dia1e=QEvaluate[I ScaleMu^(4(4-D)) dia1e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype2a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia2a=-(Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, c1]);


dia2a=FourierXP[dia2a,{x,q}]/.SUNN->CA;

dia2a=QEvaluate[I ScaleMu^(4(4-D)) dia2a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype2b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia2b=-(Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, c1]);


dia2b=FourierXP[dia2b,{x,q}]/.SUNN->CA;

dia2b=QEvaluate[I ScaleMu^(4(4-D)) dia2b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype2c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia2c=-(Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, c1]);


dia2c=FourierXP[dia2c,{x,q}]/.SUNN->CA;

dia2c=QEvaluate[I ScaleMu^(4(4-D)) dia2c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype2d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia2d=-(Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, c1]);


dia2d=FourierXP[dia2d,{x,q}]/.SUNN->CA;

dia2d=QEvaluate[I ScaleMu^(4(4-D)) dia2d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype2e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia2e=-(Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, c1]);


dia2e=FourierXP[dia2e,{x,q}]/.SUNN->CA;

dia2e=QEvaluate[I ScaleMu^(4(4-D)) dia2e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype3a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia3a=-(Contract[atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]);


dia3a=FourierXP[dia3a,{x,q}]/.SUNN->CA;

dia3a=QEvaluate[I ScaleMu^(4(4-D)) dia3a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype3b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia3b=-(Contract[atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]);


dia3b=FourierXP[dia3b,{x,q}]/.SUNN->CA;

dia3b=QEvaluate[I ScaleMu^(4(4-D)) dia3b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype3c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia3c=-(Contract[atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]);


dia3c=FourierXP[dia3c,{x,q}]/.SUNN->CA;

dia3c=QEvaluate[I ScaleMu^(4(4-D)) dia3c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype3d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia3d=-(Contract[atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]);


dia3d=FourierXP[dia3d,{x,q}]/.SUNN->CA;

dia3d=QEvaluate[I ScaleMu^(4(4-D)) dia3d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype3e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia3e=-(Contract[atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]);


dia3e=FourierXP[dia3e,{x,q}]/.SUNN->CA;

dia3e=QEvaluate[I ScaleMu^(4(4-D)) dia3e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype4a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia4a=Contract[atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv4 . xprop[-x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1];


dia4a=FourierXP[dia4a,{x,q}]/.SUNN->CA;

dia4a=QEvaluate[I ScaleMu^(4(4-D)) dia4a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype4b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia4b=Contract[atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv4 . xprop[-x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1];


dia4b=FourierXP[dia4b,{x,q}]/.SUNN->CA;

dia4b=QEvaluate[I ScaleMu^(4(4-D)) dia4b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype4c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia4c=Contract[atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv4 . xprop[-x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1];


dia4c=FourierXP[dia4c,{x,q}]/.SUNN->CA;

dia4c=QEvaluate[I ScaleMu^(4(4-D)) dia4c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype4d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia4d=Contract[atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv4 . xprop[-x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1];


dia4d=FourierXP[dia4d,{x,q}]/.SUNN->CA;

dia4d=QEvaluate[I ScaleMu^(4(4-D)) dia4d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype4e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia4e=Contract[atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv4 . xprop[-x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1];


dia4e=FourierXP[dia4e,{x,q}]/.SUNN->CA;

dia4e=QEvaluate[I ScaleMu^(4(4-D)) dia4e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype5a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia5a=-(Contract[atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, a1]);


dia5a=FourierXP[dia5a,{x,q}]/.SUNN->CA;

dia5a=QEvaluate[I ScaleMu^(4(4-D)) dia5a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype5b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia5b=-(Contract[atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, a1]);


dia5b=FourierXP[dia5b,{x,q}]/.SUNN->CA;

dia5b=QEvaluate[I ScaleMu^(4(4-D)) dia5b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype5c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia5c=-(Contract[atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, a1]);


dia5c=FourierXP[dia5c,{x,q}]/.SUNN->CA;

dia5c=QEvaluate[I ScaleMu^(4(4-D)) dia5c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype5d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia5d=-(Contract[atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, a1]);


dia5d=FourierXP[dia5d,{x,q}]/.SUNN->CA;

dia5d=QEvaluate[I ScaleMu^(4(4-D)) dia5d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype5e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia5e=-(Contract[atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, a1]);


dia5e=FourierXP[dia5e,{x,q}]/.SUNN->CA;

dia5e=QEvaluate[I ScaleMu^(4(4-D)) dia5e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype6a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia6a=-(Contract[atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2]);


dia6a=FourierXP[dia6a,{x,q}]/.SUNN->CA;

dia6a=QEvaluate[I ScaleMu^(4(4-D)) dia6a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype6b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia6b=-(Contract[atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2]);


dia6b=FourierXP[dia6b,{x,q}]/.SUNN->CA;

dia6b=QEvaluate[I ScaleMu^(4(4-D)) dia6b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype6c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia6c=-(Contract[atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2]);


dia6c=FourierXP[dia6c,{x,q}]/.SUNN->CA;

dia6c=QEvaluate[I ScaleMu^(4(4-D)) dia6c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype6d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia6d=-(Contract[atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2]);


dia6d=FourierXP[dia6d,{x,q}]/.SUNN->CA;

dia6d=QEvaluate[I ScaleMu^(4(4-D)) dia6d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype6e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia6e=-(Contract[atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2]);


dia6e=FourierXP[dia6e,{x,q}]/.SUNN->CA;

dia6e=QEvaluate[I ScaleMu^(4(4-D)) dia6e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype7a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia7a=-(Contract[atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2]);


dia7a=FourierXP[dia7a,{x,q}]/.SUNN->CA;

dia7a=QEvaluate[I ScaleMu^(4(4-D)) dia7a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype7b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia7b=-(Contract[atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2]);


dia7b=FourierXP[dia7b,{x,q}]/.SUNN->CA;

dia7b=QEvaluate[I ScaleMu^(4(4-D)) dia7b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype7c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia7c=-(Contract[atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2]);


dia7c=FourierXP[dia7c,{x,q}]/.SUNN->CA;

dia7c=QEvaluate[I ScaleMu^(4(4-D)) dia7c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype7d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia7d=-(Contract[atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2]);


dia7d=FourierXP[dia7d,{x,q}]/.SUNN->CA;

dia7d=QEvaluate[I ScaleMu^(4(4-D)) dia7d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype7e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia7e=-(Contract[atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2]);


dia7e=FourierXP[dia7e,{x,q}]/.SUNN->CA;

dia7e=QEvaluate[I ScaleMu^(4(4-D)) dia7e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype8a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia8a=Contract[atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia8a=FourierXP[dia8a,{x,q}]/.SUNN->CA;

dia8a=QEvaluate[I ScaleMu^(4(4-D)) dia8a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype8b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia8b=Contract[atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia8b=FourierXP[dia8b,{x,q}]/.SUNN->CA;

dia8b=QEvaluate[I ScaleMu^(4(4-D)) dia8b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype8c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia8c=Contract[atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia8c=FourierXP[dia8c,{x,q}]/.SUNN->CA;

dia8c=QEvaluate[I ScaleMu^(4(4-D)) dia8c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype8d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia8d=Contract[atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia8d=FourierXP[dia8d,{x,q}]/.SUNN->CA;

dia8d=QEvaluate[I ScaleMu^(4(4-D)) dia8d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype8e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia8e=Contract[atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia8e=FourierXP[dia8e,{x,q}]/.SUNN->CA;

dia8e=QEvaluate[I ScaleMu^(4(4-D)) dia8e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype9a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia9a=-(Contract[atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]);


dia9a=FourierXP[dia9a,{x,q}]/.SUNN->CA;

dia9a=QEvaluate[I ScaleMu^(4(4-D)) dia9a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype9b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia9b=-(Contract[atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]);


dia9b=FourierXP[dia9b,{x,q}]/.SUNN->CA;

dia9b=QEvaluate[I ScaleMu^(4(4-D)) dia9b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype9c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia9c=-(Contract[atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]);


dia9c=FourierXP[dia9c,{x,q}]/.SUNN->CA;

dia9c=QEvaluate[I ScaleMu^(4(4-D)) dia9c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype9d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia9d=-(Contract[atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]);


dia9d=FourierXP[dia9d,{x,q}]/.SUNN->CA;

dia9d=QEvaluate[I ScaleMu^(4(4-D)) dia9d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype9e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia9e=-(Contract[atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]);


dia9e=FourierXP[dia9e,{x,q}]/.SUNN->CA;

dia9e=QEvaluate[I ScaleMu^(4(4-D)) dia9e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype10a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia10a=Contract[atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1];


dia10a=FourierXP[dia10a,{x,q}]/.SUNN->CA;

dia10a=QEvaluate[I ScaleMu^(4(4-D)) dia10a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype10b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia10b=Contract[atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1];


dia10b=FourierXP[dia10b,{x,q}]/.SUNN->CA;

dia10b=QEvaluate[I ScaleMu^(4(4-D)) dia10b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype10c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia10c=Contract[atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1];


dia10c=FourierXP[dia10c,{x,q}]/.SUNN->CA;

dia10c=QEvaluate[I ScaleMu^(4(4-D)) dia10c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype10d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia10d=Contract[atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1];


dia10d=FourierXP[dia10d,{x,q}]/.SUNN->CA;

dia10d=QEvaluate[I ScaleMu^(4(4-D)) dia10d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype10e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia10e=Contract[atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1];


dia10e=FourierXP[dia10e,{x,q}]/.SUNN->CA;

dia10e=QEvaluate[I ScaleMu^(4(4-D)) dia10e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype11a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia11a=Contract[atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia11a=FourierXP[dia11a,{x,q}]/.SUNN->CA;

dia11a=QEvaluate[I ScaleMu^(4(4-D)) dia11a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype11b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia11b=Contract[atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia11b=FourierXP[dia11b,{x,q}]/.SUNN->CA;

dia11b=QEvaluate[I ScaleMu^(4(4-D)) dia11b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype11c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia11c=Contract[atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia11c=FourierXP[dia11c,{x,q}]/.SUNN->CA;

dia11c=QEvaluate[I ScaleMu^(4(4-D)) dia11c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype11d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia11d=Contract[atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia11d=FourierXP[dia11d,{x,q}]/.SUNN->CA;

dia11d=QEvaluate[I ScaleMu^(4(4-D)) dia11d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype11e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia11e=Contract[atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia11e=FourierXP[dia11e,{x,q}]/.SUNN->CA;

dia11e=QEvaluate[I ScaleMu^(4(4-D)) dia11e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype12a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia12a=Contract[atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2];


dia12a=FourierXP[dia12a,{x,q}]/.SUNN->CA;

dia12a=QEvaluate[I ScaleMu^(4(4-D)) dia12a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype12b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia12b=Contract[atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2];


dia12b=FourierXP[dia12b,{x,q}]/.SUNN->CA;

dia12b=QEvaluate[I ScaleMu^(4(4-D)) dia12b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype12c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia12c=Contract[atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2];


dia12c=FourierXP[dia12c,{x,q}]/.SUNN->CA;

dia12c=QEvaluate[I ScaleMu^(4(4-D)) dia12c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype12d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia12d=Contract[atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2];


dia12d=FourierXP[dia12d,{x,q}]/.SUNN->CA;

dia12d=QEvaluate[I ScaleMu^(4(4-D)) dia12d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype12e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia12e=Contract[atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2];


dia12e=FourierXP[dia12e,{x,q}]/.SUNN->CA;

dia12e=QEvaluate[I ScaleMu^(4(4-D)) dia12e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype13a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia13a=Contract[atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2];


dia13a=FourierXP[dia13a,{x,q}]/.SUNN->CA;

dia13a=QEvaluate[I ScaleMu^(4(4-D)) dia13a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype13b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia13b=Contract[atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2];


dia13b=FourierXP[dia13b,{x,q}]/.SUNN->CA;

dia13b=QEvaluate[I ScaleMu^(4(4-D)) dia13b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype13c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia13c=Contract[atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2];


dia13c=FourierXP[dia13c,{x,q}]/.SUNN->CA;

dia13c=QEvaluate[I ScaleMu^(4(4-D)) dia13c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype13d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia13d=Contract[atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2];


dia13d=FourierXP[dia13d,{x,q}]/.SUNN->CA;

dia13d=QEvaluate[I ScaleMu^(4(4-D)) dia13d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype13e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia13e=Contract[atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2];


dia13e=FourierXP[dia13e,{x,q}]/.SUNN->CA;

dia13e=QEvaluate[I ScaleMu^(4(4-D)) dia13e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype14a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia14a=Contract[atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . xprop[-x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, a1];


dia14a=FourierXP[dia14a,{x,q}]/.SUNN->CA;

dia14a=QEvaluate[I ScaleMu^(4(4-D)) dia14a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype14b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia14b=Contract[atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . xprop[-x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, a1];


dia14b=FourierXP[dia14b,{x,q}]/.SUNN->CA;

dia14b=QEvaluate[I ScaleMu^(4(4-D)) dia14b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype14c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia14c=Contract[atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . xprop[-x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, a1];


dia14c=FourierXP[dia14c,{x,q}]/.SUNN->CA;

dia14c=QEvaluate[I ScaleMu^(4(4-D)) dia14c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype14d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia14d=Contract[atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . xprop[-x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, a1];


dia14d=FourierXP[dia14d,{x,q}]/.SUNN->CA;

dia14d=QEvaluate[I ScaleMu^(4(4-D)) dia14d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype14e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia14e=Contract[atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . xprop[-x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, a1];


dia14e=FourierXP[dia14e,{x,q}]/.SUNN->CA;

dia14e=QEvaluate[I ScaleMu^(4(4-D)) dia14e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype15a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia15a=-(Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1]);


dia15a=FourierXP[dia15a,{x,q}]/.SUNN->CA;

dia15a=QEvaluate[I ScaleMu^(4(4-D)) dia15a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype15b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia15b=-(Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1]);


dia15b=FourierXP[dia15b,{x,q}]/.SUNN->CA;

dia15b=QEvaluate[I ScaleMu^(4(4-D)) dia15b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype15c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia15c=-(Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1]);


dia15c=FourierXP[dia15c,{x,q}]/.SUNN->CA;

dia15c=QEvaluate[I ScaleMu^(4(4-D)) dia15c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype15d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia15d=-(Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1]);


dia15d=FourierXP[dia15d,{x,q}]/.SUNN->CA;

dia15d=QEvaluate[I ScaleMu^(4(4-D)) dia15d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype15e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia15e=-(Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1]);


dia15e=FourierXP[dia15e,{x,q}]/.SUNN->CA;

dia15e=QEvaluate[I ScaleMu^(4(4-D)) dia15e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype16a[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16a,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia16a=Contract[atr[SUNTrace[xprop2a[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, c1];


dia16a=FourierXP[dia16a,{x,q}]/.SUNN->CA;

dia16a=QEvaluate[I ScaleMu^(4(4-D)) dia16a,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype16b[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16b,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia16b=Contract[atr[SUNTrace[xprop2b[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, c1];


dia16b=FourierXP[dia16b,{x,q}]/.SUNN->CA;

dia16b=QEvaluate[I ScaleMu^(4(4-D)) dia16b,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype16c[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16c,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia16c=Contract[atr[SUNTrace[xprop2c[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, c1];


dia16c=FourierXP[dia16c,{x,q}]/.SUNN->CA;

dia16c=QEvaluate[I ScaleMu^(4(4-D)) dia16c,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype16d[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16d,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia16d=Contract[atr[SUNTrace[xprop2d[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, c1];


dia16d=FourierXP[dia16d,{x,q}]/.SUNN->CA;

dia16d=QEvaluate[I ScaleMu^(4(4-D)) dia16d,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
xtype16e[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16e,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia16e=Contract[atr[SUNTrace[xprop2e[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, c1];


dia16e=FourierXP[dia16e,{x,q}]/.SUNN->CA;

dia16e=QEvaluate[I ScaleMu^(4(4-D)) dia16e,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
(*xtype1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia1=-(Contract[atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]);


dia1=FourierXP[dia1,{x,q}]/.SUNN->CA;

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia2=Contract[atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv4 . xprop[-x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1];


dia2=FourierXP[dia2,{x,q}]/.SUNN->CA;

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia3=-(Contract[atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2]);


dia3=FourierXP[dia3,{x,q}]/.SUNN->CA;

dia3=QEvaluate[I ScaleMu^(4(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia4=Contract[atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia4=FourierXP[dia4,{x,q}]/.SUNN->CA;

dia4=QEvaluate[I ScaleMu^(4(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia5=Contract[atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . xprop[-x]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, a1];


dia5=FourierXP[dia5,{x,q}]/.SUNN->CA;

dia5=QEvaluate[I ScaleMu^(4(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia6=Contract[atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1];


dia6=FourierXP[dia6,{x,q}]/.SUNN->CA;

dia6=QEvaluate[I ScaleMu^(4(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia7=Contract[atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia7=FourierXP[dia7,{x,q}]/.SUNN->CA;

dia7=QEvaluate[I ScaleMu^(4(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia8=-(Contract[atr[SUNTrace[hv4 . xprop[-x] . v2 . xprop[x]]]*atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]);


dia8=FourierXP[dia8,{x,q}]/.SUNN->CA;

dia8=QEvaluate[I ScaleMu^(4(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia9=-(Contract[atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, a1]);


dia9=FourierXP[dia9,{x,q}]/.SUNN->CA;

dia9=QEvaluate[I ScaleMu^(4(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia10=-(Contract[atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . hv4 . xprop[-x] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . xprop[x] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1]);


dia10=FourierXP[dia10,{x,q}]/.SUNN->CA;

dia10=QEvaluate[I ScaleMu^(4(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia11=Contract[atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, c1];


dia11=FourierXP[dia11,{x,q}]/.SUNN->CA;

dia11=QEvaluate[I ScaleMu^(4(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia12=-(Contract[atr[SUNTrace[hv4 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1]);


dia12=FourierXP[dia12,{x,q}]/.SUNN->CA;

dia12=QEvaluate[I ScaleMu^(4(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia13=Contract[atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2];


dia13=FourierXP[dia13,{x,q}]/.SUNN->CA;

dia13=QEvaluate[I ScaleMu^(4(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia14=-(Contract[atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . hv3 . xprop[-x] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2]);


dia14=FourierXP[dia14,{x,q}]/.SUNN->CA;

dia14=QEvaluate[I ScaleMu^(4(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia15=-(Contract[atr[SUNTrace[hv3 . xprop[-x] . v1 . xprop[x]]]*atr[SUNTrace[hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, c1]);


dia15=FourierXP[dia15,{x,q}]/.SUNN->CA;

dia15=QEvaluate[I ScaleMu^(4(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Input::Initialization:: *)
(*xtype16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia16=Contract[atr[SUNTrace[xprop2[x, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[hv3 . xprop[-x] . v2 . xprop[x] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2];


dia16=FourierXP[dia16,{x,q}]/.SUNN->CA;

dia16=QEvaluate[I ScaleMu^(4(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]*)


(* ::Code::Initialization::Plain:: *)
(*--------------------------------------*)


(* ::Input::Initialization:: *)
ptype1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia1=(-I)*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[-k + l] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1] . v1 . prop[-k + q] . hv4 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[l]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, a1];


dia1=IntegrateP[dia1,{k1,k2,k,l}]/.SUNN->CA;

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia2=(-I)*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[k - q] . v2 . prop[k - l]]]*atr[SUNTrace[hv4 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]*atr[SUNTrace[v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1]]]]*FAD[l]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2];


dia2=IntegrateP[dia2,{k1,k2,k,l}]/.SUNN->CA;

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia3=I*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]*atr[SUNTrace[v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1] . v1 . prop[-k + q] . hv3 . prop[-k + l]]]]*FAD[l]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, a1];


dia3=IntegrateP[dia3,{k1,k2,k,l}]/.SUNN->CA;

dia3=QEvaluate[I ScaleMu^(4(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia4=I*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]*atr[SUNTrace[hv3 . prop[k - q] . v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1] . v2 . prop[k - l]]]]*FAD[l]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2];


dia4=IntegrateP[dia4,{k1,k2,k,l}]/.SUNN->CA;

dia4=QEvaluate[I ScaleMu^(4(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia5=I*gStrong^2*Contract[atr[SUNTrace[v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1]]]*atr[SUNTrace[hv4 . prop[k - q] . v1 . prop[k - l] . hv3 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[l]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1];


dia5=IntegrateP[dia5,{k1,k2,k,l}]/.SUNN->CA;

dia5=QEvaluate[I ScaleMu^(4(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia6=I*gStrong^2*Contract[atr[SUNTrace[v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1]]]*atr[SUNTrace[hv4 . prop[k - q] . v2 . prop[k - l] . hv3 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[l]*FlavorDelta[a2, d2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia6=IntegrateP[dia6,{k1,k2,k,l}]/.SUNN->CA;

dia6=QEvaluate[I ScaleMu^(4(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia7=I*gStrong^2*Contract[atr[SUNTrace[v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1]]]*atr[SUNTrace[hv3 . prop[k - q] . v2 . prop[k - l] . hv4 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[l]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2];


dia7=IntegrateP[dia7,{k1,k2,k,l}]/.SUNN->CA;

dia7=QEvaluate[I ScaleMu^(4(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia8=I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]*atr[SUNTrace[hv4 . prop[k - q] . v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1] . v2 . prop[k - l]]]]*FAD[l]*FlavorDelta[a2, b2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia8=IntegrateP[dia8,{k1,k2,k,l}]/.SUNN->CA;

dia8=QEvaluate[I ScaleMu^(4(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia9=I*gStrong^2*Contract[atr[SUNTrace[v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1]]]*atr[SUNTrace[hv3 . prop[k - q] . v1 . prop[k - l] . hv4 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[l]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, c1];


dia9=IntegrateP[dia9,{k1,k2,k,l}]/.SUNN->CA;

dia9=QEvaluate[I ScaleMu^(4(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia10=(-I)*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[k - q] . v1 . prop[k - l]]]*atr[SUNTrace[hv3 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]*atr[SUNTrace[v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1]]]]*FAD[l]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1];


dia10=IntegrateP[dia10,{k1,k2,k,l}]/.SUNN->CA;

dia10=QEvaluate[I ScaleMu^(4(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia11=(-I)*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[k - q] . v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1] . v2 . prop[k - l] . hv3 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[l]*FlavorDelta[a2, d2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2];


dia11=IntegrateP[dia11,{k1,k2,k,l}]/.SUNN->CA;

dia11=QEvaluate[I ScaleMu^(4(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia12=(-I)*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[k - q] . v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1] . v2 . prop[k - l] . hv4 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[l]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2];


dia12=IntegrateP[dia12,{k1,k2,k,l}]/.SUNN->CA;

dia12=QEvaluate[I ScaleMu^(4(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia13=(-I)*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[k - q] . v1 . prop[k - l]]]*atr[SUNTrace[hv4 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]*atr[SUNTrace[v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1]]]]*FAD[l]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, c1];


dia13=IntegrateP[dia13,{k1,k2,k,l}]/.SUNN->CA;

dia13=QEvaluate[I ScaleMu^(4(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia14=I*gStrong^2*Contract[atr[SUNTrace[hv3 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]*atr[SUNTrace[v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1] . v1 . prop[-k + q] . hv4 . prop[-k + l]]]]*FAD[l]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1];


dia14=IntegrateP[dia14,{k1,k2,k,l}]/.SUNN->CA;

dia14=QEvaluate[I ScaleMu^(4(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia15=(-I)*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[-k + l] . v2 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1] . v1 . prop[-k + q] . hv3 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]]*FAD[l]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1];


dia15=IntegrateP[dia15,{k1,k2,k,l}]/.SUNN->CA;

dia15=QEvaluate[I ScaleMu^(4(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
ptype16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16,lorz,sun,k1,k2,l,k,lor1,lor2,lor3,lor4},


dia16=(-I)*gStrong^2*Contract[atr[SUNTrace[hv4 . prop[k - q] . v2 . prop[k - l]]]*atr[SUNTrace[hv3 . prop[k2] . GAD[lorz] . SUNT[sun] . prop[k2 + l]]]*atr[SUNTrace[v1 . prop[k1 + l] . GAD[lorz] . SUNT[sun] . prop[k1]]]]*FAD[l]*FlavorDelta[a2, b2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia16=IntegrateP[dia16,{k1,k2,k,l}]/.SUNN->CA;

dia16=QEvaluate[I ScaleMu^(4(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*--------------------------------------*)


(* ::Input::Initialization:: *)
iptype1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia1=-(Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv4 . prop[k - q] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . prop[k - l] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2]);


dia1=IntegrateP[dia1,{k,l}]/.SUNN->CA;

dia1=QEvaluate[I ScaleMu^(4(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia2=-(Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv3 . prop[k - q] . v2 . prop[k - l]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2]);


dia2=IntegrateP[dia2,{k,l}]/.SUNN->CA;

dia2=QEvaluate[I ScaleMu^(4(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia3=Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv3 . prop[k - q] . v1 . prop[k - l] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, c1];


dia3=IntegrateP[dia3,{k,l}]/.SUNN->CA;

dia3=QEvaluate[I ScaleMu^(4(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia4=-(Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv3 . prop[k - q] . v1 . prop[k - l]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, c1]);


dia4=IntegrateP[dia4,{k,l}]/.SUNN->CA;

dia4=QEvaluate[I ScaleMu^(4(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia5=-(Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv3 . prop[-k + l] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . prop[-k + q] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, b2]*FlavorDelta[d1, a1]);


dia5=IntegrateP[dia5,{k,l}]/.SUNN->CA;

dia5=QEvaluate[I ScaleMu^(4(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia6=Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv4 . prop[k - q] . v2 . prop[k - l] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia6=IntegrateP[dia6,{k,l}]/.SUNN->CA;

dia6=QEvaluate[I ScaleMu^(4(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia7=Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv3 . prop[k - q] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . prop[k - l]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, d2]*FlavorDelta[d1, b2];


dia7=IntegrateP[dia7,{k,l}]/.SUNN->CA;

dia7=QEvaluate[I ScaleMu^(4(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia8=-(Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv4 . prop[k - q] . v2 . prop[k - l]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, a1]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2]);


dia8=IntegrateP[dia8,{k,l}]/.SUNN->CA;

dia8=QEvaluate[I ScaleMu^(4(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia9=Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv4 . prop[k - q] . v1 . prop[k - l] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1];


dia9=IntegrateP[dia9,{k,l}]/.SUNN->CA;

dia9=QEvaluate[I ScaleMu^(4(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia10=-(Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv4 . prop[k - q] . v1 . prop[k - l]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, c1]);


dia10=IntegrateP[dia10,{k,l}]/.SUNN->CA;

dia10=QEvaluate[I ScaleMu^(4(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia11=-(Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv4 . prop[-k + l] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . prop[-k + q] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, d2]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1]);


dia11=IntegrateP[dia11,{k,l}]/.SUNN->CA;

dia11=QEvaluate[I ScaleMu^(4(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia12=Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv4 . prop[k - q] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . prop[k - l]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, c1]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia12=IntegrateP[dia12,{k,l}]/.SUNN->CA;

dia12=QEvaluate[I ScaleMu^(4(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia13=Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv3 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . prop[-k + q] . hv4 . prop[-k + l]]]]*FlavorDelta[a2, b2]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, a1];


dia13=IntegrateP[dia13,{k,l}]/.SUNN->CA;

dia13=QEvaluate[I ScaleMu^(4(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia14=Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . v2 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v1 . prop[-k + q] . hv3 . prop[-k + l]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, d2]*FlavorDelta[d1, a1];


dia14=IntegrateP[dia14,{k,l}]/.SUNN->CA;

dia14=QEvaluate[I ScaleMu^(4(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia15=Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2]]]*atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv3 . prop[k - q] . v2 . prop[k - l] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, a1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2];


dia15=IntegrateP[dia15,{k,l}]/.SUNN->CA;

dia15=QEvaluate[I ScaleMu^(4(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Input::Initialization:: *)
iptype16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16,lorz,sun,l,k,lor1,lor2,lor3,lor4},


dia16=-(Contract[atr[SUNTrace[prop2[l, lor1, lor2, lor3, lor4] . hv3 . prop[k - q] . v1 . GAD[lor1] . GAD[lorz] . SUNT[sun] . GAD[lor2] . v2 . prop[k - l] . hv4 . GAD[lor3] . GAD[lorz] . SUNT[sun] . GAD[lor4]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, c1]*FlavorDelta[c2, b2]*FlavorDelta[d1, d2]);


dia16=IntegrateP[dia16,{k,l}]/.SUNN->CA;

dia16=QEvaluate[I ScaleMu^(4(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
]


(* ::Code::Initialization::Plain:: *)
(*--------------------------------------*)


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
