(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)



(* ::Code::Initialization::Plain:: *)
Fourquarkmqqgg::usage="Fourquarkmqqgg[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}}] give the \[LeftAngleBracket]GGG\[RightAngleBracket] contribution of OPE result about <J1 J2>, with J1=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a1\)]\)\!\(\*SubscriptBox[\(v1\[CapitalPsi]\), \(b1\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(c1\)]\)\!\(\*SubscriptBox[\(v2\[CapitalPsi]\), \(d1\)]\) and J2=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(b2\)]\)\!\(\*SubscriptBox[\(v3\[CapitalPsi]\), \(a2\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(d2\)]\)\!\(\*SubscriptBox[\(v4\[CapitalPsi]\), \(c2\)]\) "

Fourquarkmqqgg::inderr="Dummy indices conflict!"

Begin["`Private`Fourquarkmqqgg`"]



(* ::Code::Initialization::Plain:: *)
Options[Fourquarkmqqgg]={
	Parallelized->True,
	AutoNDR->True,
	HoldFlavor->False,
	TraceFirst->False
}



(* ::Code::Initialization::Plain:: *)
(*-----------------------------------*)
(* condensate *)
(* \[LeftAngleBracket]\!\(
\*SubsuperscriptBox[\(G\), \(\[Mu]\[Nu]\), \(a\)]
\*SubsuperscriptBox[\(G\), \(\[Alpha]\[Beta]\), \(b\)]\)\[RightAngleBracket] *)
congg[lora_,lorb_,lorc_,lord_]=1/(2 CA CF D(D-1))(MTD[lora,lorc]MTD[lorb,lord]-MTD[lora,lord]MTD[lorb,lorc]);
(* m\[LeftAngleBracket]Overscript[q, _]q\[RightAngleBracket] *)
cqqm[x_,f_]=-1/4*1/CA + ((I/4)*DiracGamma[Momentum[x, D], D]*quarkMass[f])/(CA*D);

(*-----------------------------------*)
(* propagators and condensates *)

xpropm[x_,f_]=(DiracGamma[Momentum[x, D], D]*qfact1[(I^(1 - D)*qGamma[D/2])/(2*Pi^(D/2))])/Pair[Momentum[x, D], Momentum[x, D]]^(D/2) + 
 Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[(I^(2 - D)*qGamma[-1 + D/2]*quarkMass[f])/(4*Pi^(D/2))];


(* ::Code::Initialization::Plain:: *)
xprogm[x_,f_,sun_,lora_,lorb_]=DiracGamma[LorentzIndex[lora, D], D] . DiracGamma[LorentzIndex[lorb, D], D] . DiracGamma[Momentum[x, D], D]*
  Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[-1/32*qGamma[-1 + D/2]/(I^D*Pi^(D/2))]*SUNT[SUNIndex[sun]] + 
 DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorb, D], D] . DiracGamma[LorentzIndex[lora, D], D]*
  Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[qfact2[1/(32*E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-1 + D/2]]*SUNT[SUNIndex[sun]] + 
 (2*DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lora, D], Momentum[x, D]]*Pair[LorentzIndex[lorb, D], Momentum[x, D]]*
   qfact1[qGamma[D/2]/(16*I^D*Pi^(D/2))]*SUNT[SUNIndex[sun]])/Pair[Momentum[x, D], Momentum[x, D]]^(D/2) + 
 Pair[Momentum[x, D], Momentum[x, D]]^(2 - D/2)*(DiracGamma[LorentzIndex[lorb, D], D] . DiracGamma[LorentzIndex[lora, D], D]*
    qfact1[qfact2[((I/64)*quarkMass[f])/(E^((I/2)*D*Pi)*Pi^(D/2))]*qGamma[-2 + D/2]] + 
   DiracGamma[LorentzIndex[lora, D], D] . DiracGamma[LorentzIndex[lorb, D], D]*qfact1[((-1/64*I)*qGamma[-2 + D/2]*quarkMass[f])/(I^D*Pi^(D/2))])*
  SUNT[SUNIndex[sun]] + DiracGamma[LorentzIndex[lorb, D], D] . DiracGamma[Momentum[x, D], D]*Pair[LorentzIndex[lora, D], Momentum[x, D]]*
  Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[((I/32)*qGamma[-1 + D/2]*quarkMass[f])/(I^D*Pi^(D/2))]*SUNT[SUNIndex[sun]] + 
 DiracGamma[Momentum[x, D], D] . DiracGamma[LorentzIndex[lorb, D], D]*Pair[LorentzIndex[lora, D], Momentum[x, D]]*
  Pair[Momentum[x, D], Momentum[x, D]]^(1 - D/2)*qfact1[((I/32)*qGamma[-1 + D/2]*quarkMass[f])/(I^D*Pi^(D/2))]*SUNT[SUNIndex[sun]];



(* ::Code::Initialization::Plain:: *)
Fourquarkmqqgg[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{v3_,v4_,{a2_,b2_,c2_,d2_}},OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,x,q,atr,hv21,hv22,listpole1,listpole2,
holdf=OptionValue[HoldFlavor],diagrams,ndr=OptionValue[AutoNDR],fdir,files,pall,dot,trfirst=OptionValue[TraceFirst]},


(*--- check dummy indices in input ---*)
If[DummyCheck[v1,v2,v3,v4]==0,
	Message[Fourquarkmqqgg::inderr];
	Abort[]
];
;


(*-------------------------------*)
(* B A^+ B *)
hv21=v3//ComplexConjugate;
hv22=v4//ComplexConjugate//FCI;
(* ComplexConjugate[sigma] will expand the DiracSigma to gamma matrices, recover it to DiracSigma *)
hv22=(hv22/.Dot->dot)/.f1_ dot[aa___,DiracGamma[LorentzIndex[lor1_,dim___],dim___],DiracGamma[LorentzIndex[lor2_,dim___],dim___],bb___]+f2_ dot[aa___,DiracGamma[LorentzIndex[lor2_,dim___],dim___],DiracGamma[LorentzIndex[lor1_,dim___],dim___],bb___]/;(f1+f2===0):>-2I f1 dot[aa,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],bb];
hv22=hv22/.{null->1,dot->Dot};
(*--------------------------------------------*)		

					
If[DirectoryQ[OptionValue[Parallelized]//ToString],
	fdir=OptionValue[Parallelized];
	pall="External"
,
	fdir="None";
	pall=OptionValue[Parallelized]
];


(*---------------------------------------------------*)


diagrams={type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16};

(*---------------------------------------------------*)

If[OptionValue[AutoNDR]===True,atr=TR5,atr=TR];

(*---------------------------------------------------*)

If[pall===True,

	DistributeDefinitions[v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
	tmp= Plus@@WaitAll[ParallelSubmit[{hv21,hv22,holdf,atr},
									#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr]]&/@diagrams];
									
			
	QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA
	

,

	If[pall==="External",
	
		DistributeDefinitions[qq,v1,v2,a1,b1,c1,d1,a2,b2,c2,d2];
		
		If[fdir==="None",
		
			ParallelSubmit[{hv21,hv22,holdf,atr},#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr]]&/@diagrams
		,
		
			files=(StringSplit[ToString[#],"`"][[-1]])&/@diagrams;
			files=("Fourquarkmqqgg_"<>#)&/@files;
			
			ImExport[fdir,
						files,
						{{hv21,hv22,holdf,atr},
						{qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr},
						diagrams}
						]
		]
		(* !!! do not touch tmp before WaitAll[tmp], any function may modifiy the things in EvaluationObject[...] should not be used here. !!! *)
	
	,
	
		
		tmp= Plus@@Map[#[qq,{v1,v2,{a1,b1,c1,d1}},{hv21,hv22,{a2,b2,c2,d2}},holdf,atr]&,diagrams];
		
		QGather[(tmp/.{CA->SUNN,CF->(SUNN^2 - 1)/(2SUNN)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA

	]
]

]


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)
(*-------------The following are generated by algorithem----------------*)
(*----------------------------------------------------------------------*)


(* ::Input::Initialization::Plain:: *)
type1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia1=Condensate["gg"]*Condensate[{a1, c2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xpropm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]*tr[str[dot[hv4, cqqm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{b2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xpropm[x, b1], hv4, xprogm[-x, c2, sun, lora, lorb], v2, cqqm[x, d1]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{c1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xpropm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord], hv4, cqqm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


 dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia1=FourierXP[dia1,{x,q}];


 dia1=dia1/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia1=QEvaluate[I ScaleMu^(2(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia1=QGather[Expand[dia1]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia2=Condensate["gg"]*Condensate[{b2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb]]]]]*tr[str[dot[hv3, xpropm[-x, a2], v2, cqqm[x, d1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{b2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb], hv4, xpropm[-x, c2], v2, cqqm[x, d1]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{d2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb]]]]]*tr[str[dot[hv4, xpropm[-x, c2], v2, cqqm[x, d1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


 dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia2=FourierXP[dia2,{x,q}];


 dia2=dia2/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia2=QEvaluate[I ScaleMu^(2(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia2=QGather[Expand[dia2]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia3=-(Condensate["gg"]*Condensate[{a1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, cqqm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord], hv4, xprogm[-x, c2, sun, lora, lorb], v2, xpropm[x, d1]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["gg"]*Condensate[{b2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xpropm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord], hv4, xprogm[-x, c2, sun, lora, lorb], v2, cqqm[x, d1]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{c1, c2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xpropm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord]]]]*tr[str[dot[hv4, cqqm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


 dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia3=FourierXP[dia3,{x,q}];


 dia3=dia3/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia3=QEvaluate[I ScaleMu^(2(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia3=QGather[Expand[dia3]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia4=-(Condensate["gg"]*Condensate[{d2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xpropm[-x, a2], v1, cqqm[x, b1], hv4, xprogm[-x, c2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["gg"]*Condensate[{d2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb], hv3, xpropm[-x, a2], v2, cqqm[x, d1]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] - Condensate["gg"]*Condensate[{c1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xpropm[x, b1], hv3, cqqm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


 dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia4=FourierXP[dia4,{x,q}];


 dia4=dia4/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia4=QEvaluate[I ScaleMu^(2(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia4=QGather[Expand[dia4]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia5=Condensate["gg"]*Condensate[{b2, d1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lora, lorb], v2, cqqm[x, d1]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xpropm[x, b1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{c1, a2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, cqqm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]*tr[str[dot[hv4, xpropm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{b2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]*tr[str[dot[hv3, xpropm[-x, a2], v1, cqqm[x, b1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


 dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia5=FourierXP[dia5,{x,q}];


 dia5=dia5/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia5=QEvaluate[I ScaleMu^(2(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia5=QGather[Expand[dia5]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia6=Condensate["gg"]*Condensate[{c1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb]]]]]*tr[str[dot[hv3, cqqm[-x, a2], v2, xpropm[x, d1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{d2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, cqqm[x, b1], hv4, xpropm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{d2, d1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xpropm[x, b1]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lora, lorb], v2, cqqm[x, d1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


 dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia6=FourierXP[dia6,{x,q}];


 dia6=dia6/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia6=QEvaluate[I ScaleMu^(2(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia6=QGather[Expand[dia6]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia7=Condensate["gg"]*Condensate[{d2, b1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lora, lorb], v2, xpropm[x, d1]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, cqqm[x, b1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{a1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, cqqm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord], hv4, xpropm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{b2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xpropm[-x, c2], v1, cqqm[x, b1], hv3, xprogm[-x, a2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


 dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia7=FourierXP[dia7,{x,q}];


 dia7=dia7/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia7=QEvaluate[I ScaleMu^(2(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia7=QGather[Expand[dia7]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia8=-(Condensate["gg"]*Condensate[{d2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, cqqm[x, b1], hv4, xprogm[-x, c2, sun, lora, lorb], v2, xpropm[x, d1]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) + Condensate["gg"]*Condensate[{c1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb]]]]]*tr[str[dot[hv4, cqqm[-x, c2], v2, xpropm[x, d1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{a1, a2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, cqqm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord]]]]*tr[str[dot[hv4, xpropm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


 dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia8=FourierXP[dia8,{x,q}];


 dia8=dia8/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia8=QEvaluate[I ScaleMu^(2(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia8=QGather[Expand[dia8]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia9=-(Condensate["gg"]*Condensate[{c1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb], hv4, cqqm[-x, c2], v2, xpropm[x, d1]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["gg"]*Condensate[{a1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, cqqm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord], hv3, xpropm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] - Condensate["gg"]*Condensate[{c1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb], hv3, cqqm[-x, a2], v2, xpropm[x, d1]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


 dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia9=FourierXP[dia9,{x,q}];


 dia9=dia9/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia9=QEvaluate[I ScaleMu^(2(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia9=QGather[Expand[dia9]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia10=Condensate["gg"]*Condensate[{d2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]*tr[str[dot[hv4, xpropm[-x, c2], v1, cqqm[x, b1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{b2, d1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lora, lorb], v2, cqqm[x, d1]]]]*tr[str[dot[hv4, xpropm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{c1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xpropm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord], hv3, cqqm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


 dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia10=FourierXP[dia10,{x,q}];


 dia10=dia10/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia10=QEvaluate[I ScaleMu^(2(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia10=QGather[Expand[dia10]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia11=Condensate["gg"]*Condensate[{c1, a2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, cqqm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xpropm[x, b1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{d2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xpropm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord], hv3, xprogm[-x, a2, sun, lora, lorb], v2, cqqm[x, d1]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{b2, b1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, cqqm[x, b1]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lora, lorb], v2, xpropm[x, d1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


 dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia11=FourierXP[dia11,{x,q}];


 dia11=dia11/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia11=QEvaluate[I ScaleMu^(2(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia11=QGather[Expand[dia11]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia12=-(Condensate["gg"]*Condensate[{a1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, cqqm[-x, a2], v1, xpropm[x, b1], hv4, xprogm[-x, c2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["gg"]*Condensate[{b2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, cqqm[x, b1], hv3, xpropm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{a1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]*tr[str[dot[hv3, cqqm[-x, a2], v1, xpropm[x, b1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


 dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia12=FourierXP[dia12,{x,q}];


 dia12=dia12/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia12=QEvaluate[I ScaleMu^(2(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia12=QGather[Expand[dia12]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia13=-(Condensate["gg"]*Condensate[{b2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, cqqm[x, b1], hv3, xprogm[-x, a2, sun, lora, lorb], v2, xpropm[x, d1]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) - Condensate["gg"]*Condensate[{d2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xpropm[x, b1], hv3, xprogm[-x, a2, sun, lora, lorb], v2, cqqm[x, d1]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{d2, d1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xpropm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lora, lorb], v2, cqqm[x, d1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


 dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia13=FourierXP[dia13,{x,q}];


 dia13=dia13/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia13=QEvaluate[I ScaleMu^(2(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia13=QGather[Expand[dia13]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia14=Condensate["gg"]*Condensate[{d2, b1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xpropm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, cqqm[x, b1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{c1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xpropm[x, b1], hv4, cqqm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{a1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, cqqm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord], hv3, xprogm[-x, a2, sun, lora, lorb], v2, xpropm[x, d1]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


 dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia14=FourierXP[dia14,{x,q}];


 dia14=dia14/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia14=QEvaluate[I ScaleMu^(2(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia14=QGather[Expand[dia14]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia15=Condensate["gg"]*Condensate[{a1, c2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lora, lorb], v2, xpropm[x, d1]]]]*tr[str[dot[hv4, cqqm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{a1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, cqqm[-x, c2], v1, xpropm[x, b1], hv3, xprogm[-x, a2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{c1, c2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xpropm[x, b1]]]]*tr[str[dot[hv4, cqqm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


 dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia15=FourierXP[dia15,{x,q}];


 dia15=dia15/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia15=QEvaluate[I ScaleMu^(2(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia15=QGather[Expand[dia15]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Input::Initialization::Plain:: *)
type16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


 dia16=Condensate["gg"]*Condensate[{a1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]*tr[str[dot[hv4, cqqm[-x, c2], v1, xpropm[x, b1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{a1, a2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, cqqm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lora, lorb], v2, xpropm[x, d1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{b2, b1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, cqqm[x, b1]]]]*tr[str[dot[hv4, xpropm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


 dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf],dot[aa___,1,bb___]:>dot[aa,bb]}/.dot->Dot/.str->SUNTrace/.{tr->atr,contract->Contract}/.sunSimplify->SUNSimplify;


 dia16=FourierXP[dia16,{x,q}];


 dia16=dia16/.Power[_quarkMass,_]->0/. _ quarkMass[f1_] quarkMass[f2_]->0;


 dia16=QEvaluate[I ScaleMu^(2(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0}/.CF->(CA^2-1)/(2CA);


 dia16=QGather[Expand[dia16]/.Power[_quarkMass,_]->0/. quarkMass[f1_] quarkMass[f2_]->0,q,ShowasTable->False]/.q->qq


]


(* ::Code::Initialization::Plain:: *)
(*----------------------------------------------------------------------*)


(* ::Input::Initialization::Plain:: *)
(*type1[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia1,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia1=-(Condensate["gg"]*Condensate[{a1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, cqqm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord], hv3, xpropm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) + Condensate["gg"]*Condensate[{a1, a2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, cqqm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lora, lorb], v2, xpropm[x, d1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{b2, b1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, cqqm[x, b1]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lora, lorb], v2, xpropm[x, d1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia1=dia1/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia1=SUNSimplify[FourierXP[dia1,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia1=dia1/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia1=QEvaluate[I ScaleMu^(2(4-D)) dia1,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type2[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia2,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia2=Condensate["gg"]*Condensate[{b2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb]]]]]*tr[str[dot[hv3, xpropm[-x, a2], v2, cqqm[x, d1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{a1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]*tr[str[dot[hv4, cqqm[-x, c2], v1, xpropm[x, b1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{c1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb]]]]]*tr[str[dot[hv4, cqqm[-x, c2], v2, xpropm[x, d1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia2=dia2/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia2=SUNSimplify[FourierXP[dia2,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia2=dia2/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia2=QEvaluate[I ScaleMu^(2(4-D)) dia2,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type3[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia3,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia3=Condensate["gg"]*Condensate[{c1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb]]]]]*tr[str[dot[hv3, cqqm[-x, a2], v2, xpropm[x, d1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{d2, b1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xpropm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, cqqm[x, b1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{b2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb], hv4, xpropm[-x, c2], v2, cqqm[x, d1]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia3=dia3/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia3=SUNSimplify[FourierXP[dia3,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia3=dia3/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia3=QEvaluate[I ScaleMu^(2(4-D)) dia3,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type4[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia4,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia4=Condensate["gg"]*Condensate[{a1, c2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xpropm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]*tr[str[dot[hv4, cqqm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{b2, d1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lora, lorb], v2, cqqm[x, d1]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xpropm[x, b1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{a1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, cqqm[-x, a2], v1, xpropm[x, b1], hv4, xprogm[-x, c2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2];


dia4=dia4/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia4=SUNSimplify[FourierXP[dia4,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia4=dia4/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia4=QEvaluate[I ScaleMu^(2(4-D)) dia4,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type5[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia5,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia5=-(Condensate["gg"]*Condensate[{a1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, cqqm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord], hv4, xprogm[-x, c2, sun, lora, lorb], v2, xpropm[x, d1]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["gg"]*Condensate[{d2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb], hv3, xpropm[-x, a2], v2, cqqm[x, d1]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{c1, c2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xpropm[x, b1]]]]*tr[str[dot[hv4, cqqm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia5=dia5/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia5=SUNSimplify[FourierXP[dia5,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia5=dia5/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia5=QEvaluate[I ScaleMu^(2(4-D)) dia5,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type6[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia6,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia6=-(Condensate["gg"]*Condensate[{a1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, cqqm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord], hv4, xpropm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["gg"]*Condensate[{b2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xpropm[x, b1], hv4, xprogm[-x, c2, sun, lora, lorb], v2, cqqm[x, d1]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{b2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, cqqm[x, b1], hv3, xprogm[-x, a2, sun, lora, lorb], v2, xpropm[x, d1]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia6=dia6/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia6=SUNSimplify[FourierXP[dia6,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia6=dia6/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia6=QEvaluate[I ScaleMu^(2(4-D)) dia6,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type7[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia7,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia7=Condensate["gg"]*Condensate[{c1, a2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, cqqm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xpropm[x, b1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{c1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xpropm[x, b1], hv4, cqqm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{a1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]*tr[str[dot[hv3, cqqm[-x, a2], v1, xpropm[x, b1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia7=dia7/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia7=SUNSimplify[FourierXP[dia7,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia7=dia7/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia7=QEvaluate[I ScaleMu^(2(4-D)) dia7,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type8[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia8,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia8=-(Condensate["gg"]*Condensate[{c1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xpropm[x, b1], hv3, cqqm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) + Condensate["gg"]*Condensate[{b2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]*tr[str[dot[hv3, xpropm[-x, a2], v1, cqqm[x, b1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{b2, b1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, cqqm[x, b1]]]]*tr[str[dot[hv4, xpropm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia8=dia8/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia8=SUNSimplify[FourierXP[dia8,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia8=dia8/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia8=QEvaluate[I ScaleMu^(2(4-D)) dia8,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type9[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia9,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia9=Condensate["gg"]*Condensate[{d2, b1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lora, lorb], v2, xpropm[x, d1]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, cqqm[x, b1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] + Condensate["gg"]*Condensate[{b2, d1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lora, lorb], v2, cqqm[x, d1]]]]*tr[str[dot[hv4, xpropm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{d2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xpropm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord], hv3, xprogm[-x, a2, sun, lora, lorb], v2, cqqm[x, d1]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia9=dia9/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia9=SUNSimplify[FourierXP[dia9,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia9=dia9/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia9=QEvaluate[I ScaleMu^(2(4-D)) dia9,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type10[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia10,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia10=Condensate["gg"]*Condensate[{d2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]*tr[str[dot[hv4, xpropm[-x, c2], v1, cqqm[x, b1]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{b2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xpropm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord], hv4, xprogm[-x, c2, sun, lora, lorb], v2, cqqm[x, d1]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{d2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xpropm[x, b1], hv3, xprogm[-x, a2, sun, lora, lorb], v2, cqqm[x, d1]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia10=dia10/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia10=SUNSimplify[FourierXP[dia10,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia10=dia10/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia10=QEvaluate[I ScaleMu^(2(4-D)) dia10,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type11[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia11,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia11=-(Condensate["gg"]*Condensate[{d2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, cqqm[x, b1], hv4, xprogm[-x, c2, sun, lora, lorb], v2, xpropm[x, d1]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["gg"]*Condensate[{d2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xpropm[-x, a2], v1, cqqm[x, b1], hv4, xprogm[-x, c2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{b2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xpropm[-x, c2], v1, cqqm[x, b1], hv3, xprogm[-x, a2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia11=dia11/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia11=SUNSimplify[FourierXP[dia11,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia11=dia11/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia11=QEvaluate[I ScaleMu^(2(4-D)) dia11,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type12[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia12,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia12=-(Condensate["gg"]*Condensate[{c1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xpropm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord], hv4, cqqm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) - Condensate["gg"]*Condensate[{a1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, cqqm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord], hv3, xprogm[-x, a2, sun, lora, lorb], v2, xpropm[x, d1]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{c1, c2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xpropm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord]]]]*tr[str[dot[hv4, cqqm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia12=dia12/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia12=SUNSimplify[FourierXP[dia12,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia12=dia12/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia12=QEvaluate[I ScaleMu^(2(4-D)) dia12,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type13[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia13,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia13=Condensate["gg"]*Condensate[{c1, a2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, cqqm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]*tr[str[dot[hv4, xpropm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{c1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xpropm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord], hv3, cqqm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{a1, a2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, cqqm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord]]]]*tr[str[dot[hv4, xpropm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia13=dia13/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia13=SUNSimplify[FourierXP[dia13,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia13=dia13/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia13=QEvaluate[I ScaleMu^(2(4-D)) dia13,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type14[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia14,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia14=Condensate["gg"]*Condensate[{a1, c2}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lora, lorb], v2, xpropm[x, d1]]]]*tr[str[dot[hv4, cqqm[-x, c2], v1, xprogm[x, b1, sun, lorc, lord]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, d2]*FlavorDelta[c2, a1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{d2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, cqqm[x, b1], hv4, xpropm[-x, c2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2] - Condensate["gg"]*Condensate[{a1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, cqqm[-x, c2], v1, xpropm[x, b1], hv3, xprogm[-x, a2, sun, lorc, lord], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2];


dia14=dia14/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia14=SUNSimplify[FourierXP[dia14,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia14=dia14/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia14=QEvaluate[I ScaleMu^(2(4-D)) dia14,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type15[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia15,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia15=-(Condensate["gg"]*Condensate[{c1, c2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb], hv4, cqqm[-x, c2], v2, xpropm[x, d1]]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, d2]*FlavorDelta[c2, c1]*FlavorDelta[d1, b2]) + Condensate["gg"]*Condensate[{d2, d1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xpropm[x, b1]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lora, lorb], v2, cqqm[x, d1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{d2, d1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv3, xprogm[-x, a2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb]]]]]*tr[str[dot[hv4, xpropm[-x, c2], v2, cqqm[x, d1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia15=dia15/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia15=SUNSimplify[FourierXP[dia15,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia15=dia15/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia15=QEvaluate[I ScaleMu^(2(4-D)) dia15,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Input::Initialization::Plain:: *)
(*type16[qq_,{v1_,v2_,{a1_,b1_,c1_,d1_}},{hv3_,hv4_,{a2_,b2_,c2_,d2_}},holdf_,atr_]:=Block[{x,q,dia16,lora,lorb,lorc,lord,sun,tr,str,contract,dot},


dia16=-(Condensate["gg"]*Condensate[{b2, b1}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, cqqm[x, b1], hv3, xpropm[-x, a2], v2, xprogm[x, d1, sun, lora, lorb]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2]) - Condensate["gg"]*Condensate[{c1, a2}]*contract[tr[dot[congg[lora, lorb, lorc, lord], str[dot[hv4, xprogm[-x, c2, sun, lorc, lord], v1, xprogm[x, b1, sun, lora, lorb], hv3, cqqm[-x, a2], v2, xpropm[x, d1]]]]]]*FlavorDelta[a2, c1]*FlavorDelta[b1, b2]*FlavorDelta[c2, a1]*FlavorDelta[d1, d2] + Condensate["gg"]*Condensate[{d2, d1}]*contract[congg[lora, lorb, lorc, lord]*tr[str[dot[hv3, xpropm[-x, a2], v1, xprogm[x, b1, sun, lorc, lord]]]]*tr[str[dot[hv4, xprogm[-x, c2, sun, lora, lorb], v2, cqqm[x, d1]]]]]*FlavorDelta[a2, a1]*FlavorDelta[b1, b2]*FlavorDelta[c2, c1]*FlavorDelta[d1, d2];


dia16=dia16/.{FlavorDelta[fa_,fb_,___]:>FlavorDelta[fa,fb,HoldFlavor->holdf]}/.dot->Dot/.str->SUNTrace/.tr->atr/.contract->Contract;


dia16=SUNSimplify[FourierXP[dia16,{x,q}]]/.{SUNN->CA,CF->(CA^2-1)/(2CA)};

dia16=dia16/.Power[quarkMass[_],_]->0/.quarkMass[f1_]quarkMass[f2_]->0;

dia16=QEvaluate[I ScaleMu^(2(4-D)) dia16,q,HoldFlavor->holdf,Parallelized->False]/.q->qq
];*)


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
