(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)




FourquarkPole::usage = "FourquarkPole[gamma1_,gamma_2] give the \[Epsilon]-pole at 1-loop level about the Green function \[LeftAngleBracket]\!\(\*SubscriptBox[\(J\), \(4  q\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(a\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[Psi]\), \(_\)], \(b\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(c\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[Psi]\), \(_\)], \(d\)]\)\[RightAngleBracket] with \!\(\*SubscriptBox[\(J\), \(4  q\)]\)=\!\(\*SubscriptBox[OverscriptBox[\(\[CapitalPsi]\), \(_\)], \(a\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(1\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(b\)]\)\!\(\*SubscriptBox[OverscriptBox[\(\[Psi]\), \(_\)], \(c\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(2\)]\)\!\(\*SubscriptBox[\(\[Psi]\), \(d\)]\).
The output is a list of terms ~{1/\[Epsilon],\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(2\)]\)'} "

FourquarkPole::sunindexerr = "Uknow sun structure: too much dummy indices!"
FourquarkPole::sunerr = "Reduce su(n) matrices failed: unknow su(n) structure."


Begin["`Private`FourquarkPole`"]

Options[FourquarkPole]={
	DimensionFour->True,
	KeepGamma->False,
	ShowasTable->True,
	VertexOnly->True,
	ReduceSUN->False}
	


FourquarkPole[nn_Integer cur_FourquarkCurrent,ops:OptionsPattern[]]:=Block[{},
If[ToLowerCase[ToString[OptionValue[ShowasTable]]]=="true",
	{#[[1]]nn, #[[2]], #[[3]]}&/@FourquarkPole[cur,ops]
,
	nn FourquarkPole[cur,ops]
]
]


FourquarkPole[cur_FourquarkCurrent,OptionsPattern[]]:=Block[{q,f11,f12,f21,f22,pole,pole11,pole12,pole21,pole22,pole31,pole32,lor1,lor2,lora=$AL[Unique[]],lorb=$AL[Unique[]],cola=$AC[Unique[]],colb=$AC[Unique[]],lorz=$AL[Unique[]],vtx,V1,V2,gaA,gaB,dim=D(* it's possible only 1, gamma[5] involved so can't get dimension from vertex *),
asd4=OptionValue[DimensionFour],astb=OptionValue[ShowasTable],pole0,polet,null,null0},


pole = -(I gStrong)^2 MTD[lora,lorb]/(64Pi^2 Epsilon);(* hold also for massive quark by dimensional reason *)

{f11,f12,f21,f22}={cur[[1,1]],cur[[1,-1]],cur[[2,1]],cur[[2,-1]]};


(* get the gamma matices *)
V1=If[Length[cur[[1]]]==2,1,cur[[1,2]] ];
V2=If[Length[cur[[2]]]==2,1,cur[[2,2]] ];

V1/.DiracGamma[LorentzIndex[__],di_:4]:>(dim1=di;1);

(*-----------*)

(* left diagram *)
pole11=pole vtx[1,f11,SUNT[cola] . GAD[lorz,lora] . V1 . GAD[lorb,lorz] . SUNT[cola]//DiracSimplify//SUNSimplify,f12]vtx[2,f21,V2,f22];
(* right diagram *)
pole12=pole vtx[1,f11,V1,f12]vtx[2,f21,SUNT[cola] . GAD[lorz,lora] . V2 . GAD[lorb,lorz] . SUNT[cola]//DiracSimplify//SUNSimplify,f22];


(* ---------------- *)
(* left diagram *)
pole21=-pole vtx[1,f11,SUNT[cola] . GAD[lorz,lora] . V1//SUNSimplify,f12]vtx[2,f21,SUNT[cola] . GAD[lorz,lorb] . V2//SUNSimplify,f22];
(* -1 for the direction of momentum *)
(* right diagram *)
pole22=-pole vtx[1,f11,V1 . SUNT[cola] . GAD[lora,lorz]//SUNSimplify,f12]vtx[2,f21,V2 . SUNT[cola] . GAD[lorb,lorz]//SUNSimplify,f22];


(* ---------------- *)
(* left diagram *)
pole31=pole vtx[1,f11,V1 . SUNT[cola] . GAD[lora,lorz]//SUNSimplify,f12]vtx[2,f21,SUNT[cola] . GAD[lorz,lorb] . V2//SUNSimplify,f22];
(* right diagram *)
pole32=pole vtx[1,f11,SUNT[cola] . GAD[lorz,lora] . V1//SUNSimplify,f12]vtx[2,f21,V2 . SUNT[cola] . GAD[lorb,lorz]//SUNSimplify,f22];


pole={pole11 , pole12 , pole21, pole22, pole31, pole32}/.(CA-2CF)->1/CA;


pole=simplifyv[pole,OptionValue[KeepGamma]];

(*------ reduce the SU(N) T matrix -----*)
(* (T^n)^a_b (T^n)^c_d = 1/2 delta[a,d]delta[b,c] - 1/(2Nc) delta[a,b]delta[c,d] *)
If[OptionValue[ReduceSUN]===True,	
	lor1=$AL[Unique[]];
	lor2=$AL[Unique[]];
	
	pole/.Momentum[_,di_:4]|LorentzIndex[_,di_:4]:>(dim=di;null);(* get the dimension *)
	
	pole0=Boole[FreeQ[#,_SUNT]]&/@pole;
	polet=pole0/.{1->0,0->1};

	pole0=DeleteCases[pole0 pole,{0,_,_}];
	polet=DeleteCases[polet pole,{0,_,_}];

(* for the cases that each vertex involve a same SU(N) generate *)
	If[Length[polet]>0,

		polet=(#[[1]]vtx[1, #[[2,1]], #[[2,2]]//SUNSimplify, #[[2,3]]]vtx[2, #[[3,1]], #[[3,2]]//SUNSimplify, #[[3,3]]])&/@polet;
		(* Cases[SUNT[a],sun_SUNT:>sun,Infinity] = {} *)
		polet=polet/.vtx[1,f11,vA_,f12]vtx[2,f21,vB_,f22]/;MatchQ[{Cases[vA,sun_SUNT:>sun,Infinity],Cases[vB,sun_SUNT:>sun,Infinity]},{{sun_},{sun_}}]||MatchQ[{vA,vB},{sun_SUNT,sun_SUNT}]:>(
			{gaA,gaB}={vA/._SUNT->1,vB/._SUNT->1};
			{-1/(2CA)vtx[1,f11,gaA,f12]vtx[2,f21,gaB,f22]
			,(* fierz rearrangement *)
			-1/2{vtx[1,f11,1,f22]vtx[2,f21, gaB . gaA/4 //DiracSimplify,f12](* TR[1]= 4 *)
				,
				vtx[1,f11,DiracGamma[LorentzIndex[lor1,dim],dim],f22]vtx[2,f21, gaB . DiracGamma[LorentzIndex[lor1,dim],dim] . gaA/4 //DiracSimplify,f12](* TR[ga[a].ga[b]]= 4 g^ab *)
				,
				vtx[1,f11,DiracGamma[5],f22]vtx[2,f21, gaB . DiracGamma[5] . gaA/4 //DiracSimplify,f12](* TR[DiracGamma[5].DiracGamma[5]]= 4 *)
				,
				vtx[1,f11,DiracGamma[LorentzIndex[lor1,dim],dim] . DiracGamma[5],f22]vtx[2,f21, gaB . DiracGamma[LorentzIndex[lor1,dim],dim] . DiracGamma[5] . gaA/(-4) //DiracSimplify,f12](* TR[ga[a].DiracGamma[5].ga[b].DiracGamma[5]]= -4 g^ab *)
				,
				vtx[1,f11,DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]],f22]vtx[2,f21, gaB . DiracSigma[DiracGamma[LorentzIndex[lor1,dim],dim],DiracGamma[LorentzIndex[lor2,dim],dim]] . gaA/8 //DiracSimplify,f12](* TR[DiracSigma[ga[a],ga[b]].DiracSigma[ga[c],ga[d]]]= 8 (g^ac g^bd - g^ad g^bc)/2 *)
				}
			});

		If[!FreeQ[polet,vv_vtx/;!FreeQ[vv,_SUNT]],Message[FourquarkPole::sunerr]];(* unknow su(n) structure *)
	
		polet=simplifyv[Flatten[polet],OptionValue[KeepGamma]]
	];

	pole=Join[pole0,polet];
];

(*---------------------*)

If[ToLowerCase[ToString[astb]]=="true",
(* when the SUN(N) generate has not been reduced, the flavor keep same as the original so can be ommit here *)
	If[OptionValue[VertexOnly]===True&&OptionValue[ReduceSUN]=!=True,
		Replace[pole,{aa_,vv1_,vv2_}:>{aa,vv1[[2]],vv2[[2]]},{1}]
	,
		pole
	]
,
	Plus@@(#[[1]]FourquarkCurrent[Current[#[[2,1]],#[[2,2]],#[[2,3]]],Current[#[[3,1]],#[[3,2]],#[[3,3]]]]&/@pole)
]

]/;Length[cur[[2]]]<=3


simplifyv[polist_,keepga_]:=Block[{pole,vtx,null,rule1,rule2,rule3,rule4,rule5,rule6,vlist,coll,verx},

pole=(polist/.vtx[vn_,ffa_,aa_,ffb_]:>vtx[vn,ffa,aa+null+null^2//Expand,ffb])//Expand//FCI;

If[keepga===True,
	pole=pole// Contract
,

	pole=pole/.vtx[vn_,ffa_,aa_,ffb_]:>(vlist=Expand[aa];
		vlist=List@@vlist;
		vlist=vlist(Boole[FreeQ[#,Momentum[_,D-4]|LorentzIndex[_,D-4]]]&/@vlist);(* a bit of risk here ... *)
		vlist=Plus@@vlist;
		vtx[vn,ffa,FCReplaceD[DiracSimplify[ChangeDimension[vlist,4]],D->4]/.null->0// Contract,ffb]
	)
	(*pole=pole/.vtx[vn_,ffa_,aa_,ffb_]:>vtx[vn,ffa,FCReplaceD[DiracSimplify[ChangeDimension[Replace[Expand[aa],(* a bit of risk... *)bb_/;!FreeQ[bb,Momentum[_,D-4]|LorentzIndex[_,D-4]]->0,{1}],4]],D->4]/.null->0// Contract,ffb]*)
];


(* ---------------- *)
(* rules *)
rule1=vtx[vn_,ffa_,aa_,ffb_]:>vtx[vn,ffa,DiracReduce[aa]//SUNSimplify//Expand,ffb];
rule2=SUNT[aa_] . SUNT[bb_]:>(coll=$AC[Unique[]];SUNDelta[aa,bb]/(2SUNN)+SUNT[coll]/2 (SUND[aa,bb,coll]+I SUNF[aa,bb,coll]));
rule3=vtx[vn_,ffa_,aa_,ffb_]:>vtx[vn,ffa,aa//Expand,ffb];
rule4=vtx[vn_,ffa_,aa_Plus,ffb_]:>(vtx[vn,ffa,#,ffb]&/@aa);
rule5=vtx[vn_,ffa_,aa_ bb_,ffb_]/;(And@@(FreeQ[bb,#]&/@{DiracGamma,DiracSigma,SUNT})):>vtx[vn,ffa,aa,ffb]bb;
rule6=vtx[vn_,ffa_,Pair[LorentzIndex[lor1_,di___],LorentzIndex[lor2_,di___]] ,ffb_]:>vtx[vn,ffa,1,ffb]Pair[LorentzIndex[lor1,di],LorentzIndex[lor2,di]];
(*rule6=vtx[1,aa_]vtx[2,bb_]/;!FreeQ[aa,Pair[LorentzIndex[lor1_,di___],LorentzIndex[lor2_,di___]]]&&!FreeQ[bb,lor1|lor2]:>};*)


(* ---------------- *)
(* pole/.rule1/.rule4/.rule2/.rule3/.rule4/.rule5//Expand//Total//Total//Contract//SUNSimplify//QSimplify *)
pole=pole/.rule1/.rule2/.rule3/.rule4/.rule5/.rule6//Contract//SUNSimplify;
pole=pole//Total//Expand;


(* ---------------- *)
pole=pole/.vtx[vn_,ffa_,n_,ffb_]/;NumberQ[n]:>n vtx[vn,ffa,1,ffb];

pole=pole/.vtx[vn_,ffa_,Pair[LorentzIndex[lla_,dim___],LorentzIndex[llb_,dim___]],ffb_]:>Pair[LorentzIndex[lla,dim],LorentzIndex[llb,dim]]vtx[ffa,vn,1,ffb];
(* extract the vertex doesn't involving gamma-matirx of sun-matrix *)


pole=pole/.vtx[vn_,ffa_,-v1_,ffb_]:>-vtx[vn,ffa,v1,ffb];

(*-------------------------------------*)

pole=Expand[pole]/.aa_ vtx[1,vv1__] vtx[2,vv2__]:>Contract[aa] verx[vtx[1,vv1],vtx[2,vv2]];(* isolate the vertex part *)

pole=(pole//QSimplify)/.rule3/.rule5/.CA-2CF->1/CA;
pole=(pole/.verx[-v1_vtx,v2_]:>-verx[v1,v2]/.verx[v1_,-v2_vtx]:>-verx[v1,v2])//QSimplify;
(* it's possible generate terms like vtx[1,_,-a,_]vtx[2,_,-b_] and vtx[1,_,a_]vtx[2,_,b,_] after first acting QSimplify; act QSimplify twice *)
pole=pole/.aa_ vvs_Plus/;!FreeQ[vvs,_vtx]:>Plus@@(aa List@@vvs);(* Expanding terms like (a verx[_,_] + b verx[_,_]) *)

(* the terms in list will be have different vertices *)
pole=If[Head[pole]===Plus,List@@pole,{pole}]/.verx->Times;


(* show as a table *)
pole=({FCReplaceD[#/._vtx->1,D->4],#/.aa_ vtx[1,bb__]:>{bb},#/.aa_ vtx[2,bb__]:>{bb}}&/@pole)//QNormal;
pole=ChangeDimension[pole,D](* show the gamma-matrices as D-dimensional *)
]


(*--------------------------------------------------------------*)
(*--------------------------------------------------------------*)
(* tetraquark *)


FourquarkPole[cur_FourquarkCurrent,OptionsPattern[]]:=Block[{q,p1,p2,pole,pole11,pole12,pole21,pole22,pole31,pole32,q1,q2,l,lorp=$AL[Unique[]],lorz=$AL[Unique[]],cola,colb,colc,cold,cole,colf,f11,f12,f21,f22,
vl,vr,V1,V2,coll,null,simp,sun,expand,seprate,factor,transc,dummy,ruledummy,asd4=OptionValue[DimensionFour],astb=OptionValue[ShowasTable],verx},

(* get the color indices *)
{cola,colb,colc,cold}={cur[[1,2]],cur[[1,-1]],cur[[2,3]],cur[[2,-1]]};
{f11,f12,f21,f22}={cur[[1,1]],cur[[1,-2]],cur[[2,2]],cur[[2,-2]]};

dummy=Cases[Times[cola,colb,colc,cold],Power[aa_,2]:>aa,Infinity];
If[!FreeQ[dummy,_Power],Message[FourquarkPole::sunindexerr];Abort[]];


(* get the gamma matices *)
V1=If[Length[cur[[1]]]==4,1,cur[[1,3]] ];
V2=If[Length[cur[[2]]]==6,1,cur[[2,4]] ];


pole= -1/(64Pi^2 Epsilon);(* -MTD[lora,lorb]/(64Pi^2 Epsilon), also hold for massive quark by dimensional reason *)


(*--------------------------------------------------------------------------------*)

(* -1 for the direction of momentum for some diagrams *)

(*--------------------------------------------------------------------------------*)
(* (t^n)^a_b (t^n)^c_d = 1/2 delta[a,d]delta[b,c] - 1/(2Nc) delta[a,b]delta[c,d] *)

(* left diagram *)
pole11=-pole (ColorDelta[cola,colf]ColorDelta[colb,cole]/2-ColorDelta[cola,cole]ColorDelta[colb,colf]/(2CA)) vl[cole,FCChargeConjugateTransposed[GAD[lorp,lorz],Explicit->True] . V1 . GAD[lorp,lorz],colf]vr[colc,V2,cold]//Contract//DiracSimplify;


(* right diagram *)
pole12=-pole vl[cola,V1,colb]vr[cole, GAD[lorz,lorp] . V2 . FCChargeConjugateTransposed[GAD[lorz,lorp] ,Explicit->True],colf](ColorDelta[cole,cold]ColorDelta[colf,colc]/2-ColorDelta[cole,colc]ColorDelta[cold,colf]/(2CA))//Contract//DiracSimplify;



(* ---------------- *)

(* left diagram *)
pole21=pole(ColorDelta[cola,colc]ColorDelta[cole,colf]/2-ColorDelta[cole,cola]ColorDelta[colf,colc]/(2CA)) vl[cole,FCChargeConjugateTransposed[GAD[lorp,lorz],Explicit->True] . V1,colb] vr[colf,GAD[lorz,lorp] . V2,cold];

(* right diagram *)
pole22=pole (ColorDelta[colb,cold]ColorDelta[cole,colf]/2-ColorDelta[colb,cole]ColorDelta[cold,colf]/(2CA))vl[cola,V1 . GAD[lorp,lorz],cole] vr[colc,V2 . FCChargeConjugateTransposed[GAD[lorz,lorp],Explicit->True],colf];


(* ---------------- *)

(* left diagram *)
pole31=pole (ColorDelta[colb,colc]ColorDelta[cole,colf]/2-ColorDelta[colb,cole]ColorDelta[colf,colc]/(2CA))vl[cola,V1 . GAD[lorp,lorz],cole] vr[colf,GAD[lorz,lorp] . V2,cold];

(* right diagram *)
pole32=pole(ColorDelta[cola,cold]ColorDelta[cole,colf]/2-ColorDelta[cola,cole]ColorDelta[colf,cold]/(2CA)) vl[cole, FCChargeConjugateTransposed[GAD[lorp,lorz],Explicit->True] . V1,colb] vr[colc,V2 . FCChargeConjugateTransposed[GAD[lorz,lorp],Explicit->True],colf];


(*--------------------------------------------------------------------------------*)


pole=(I gStrong)^2 {pole11 , pole12 , pole21, pole22, pole31, pole32} /.{vl[aa_]:>vl[aa+null+null^2//Expand],vr[aa_]:>vr[aa+null+null^2//Expand]}//Expand//FCI;

pole=pole//.{ColorDelta[aa_,bb_]vl[aa_,vv_,cc_]:>vl[bb,vv,cc],ColorDelta[aa_,bb_]vl[cc_,vv_,aa_]:>vl[cc,vv,bb]}//.{ColorDelta[aa_,bb_]vr[bb_,vv_,cc_]:>vr[aa,vv,cc],ColorDelta[aa_,bb_]vr[cc_,vv_,aa_]:>vr[cc,vv,bb]};


(* replace the dummy indices generated by FeynCalc, use the set of dummy indices in the input expression . *)
pole=(#/.vl[cca_,_,ccb_]vr[ccc_,_,ccd_]:>(ruledummy=(Times@@DeleteDuplicates[{cca,ccb,ccc,ccd}])/(Times@@dummy);1);(* generate replace rule Rule[a,b] by a/b *)
If[ruledummy===1,#,ruledummy=ruledummy/.aa_ bb_ Power[cc_,-1]Power[dd_,-1]:>{aa->cc,bb->dd}/.aa_ Power[bb_,-1]:>{aa->bb};#/.ruledummy])&/@pole;


(*If[ToLowerCase[ToString[asd4]]=="true",pole=FCReplaceD[ChangeDimension[pole,4],D->4]];


pole=pole/.{vl[ca_,aa_,cb_]:>vl[ca,DiracSimplify[aa]/.null->0,cb], vr[ca_,aa_,cb_]:>vr[ca,DiracSimplify[aa]/.null->0,cb]} // Contract;*)
(* delete D-4 objects, since later will set all dimension to D *)

If[OptionValue[KeepGamma]===True,
	pole=pole//Contract
,
	pole=pole/.{vl[ca_,aa_,cb_]:>vl[ca,FCReplaceD[DiracSimplify[ChangeDimension[Replace[Expand[aa],bb_/;!FreeQ[bb,Momentum[_,D-4]|LorentzIndex[_,D-4]]->0,{1}],4]],D->4]/.null->0,cb],
			 vr[ca_,aa_,cb_]:>vr[ca,FCReplaceD[DiracSimplify[ChangeDimension[Replace[Expand[aa],bb_/;!FreeQ[bb,Momentum[_,D-4]|LorentzIndex[_,D-4]]->0,{1}],4]],D->4]/.null->0,cb]} // Contract
];

(* ---------------- *)
(* rules *)
simp={vl[ca_,aa_,cb_]:>vl[ca,DiracReduce[aa]//SUNSimplify//Expand,cb],vr[ca_,aa_,cb_]:>vr[ca,DiracReduce[aa]//SUNSimplify//Expand,cb]};
sun=SUNT[aa_] . SUNT[bb_]:>(coll=$AC[Unique[]];SUNDelta[aa,bb]/(2SUNN)+SUNT[coll]/2(SUND[aa,bb,coll]+I SUNF[aa,bb,coll]));
expand={vl[ca_,aa_,cb_]:>vl[ca,aa//Expand,cb],vr[ca_,aa_,cb_]:>vr[ca,aa//Expand,cb]};
seprate={vl[ca_,aa_Plus,cb_]:>(vl[ca,#,cb]&/@aa),vr[ca_,aa_Plus,cb_]:>(vr[ca,#,cb]&/@aa)};
factor={vl[ca_,aa_ bb_,cb_]/;(And@@(FreeQ[bb,#]&/@{DiracGamma,DiracSigma,SUNT})):>vl[ca,aa,cb]bb,vr[ca_,aa_ bb_,cb_]/;(And@@(FreeQ[bb,#]&/@{DiracGamma,DiracSigma,SUNT})):>vr[ca,aa,cb]bb};




(*(* unify the color structure *)
transc={vl[colb,aa_,cola]:>-vl[cola,DiracSimplify[FCChargeConjugateTransposed[aa,Explicit->True]],colb],vr[cold,aa_,colc]:>-vr[colc,DiracSimplify[FCChargeConjugateTransposed[aa,Explicit->True]],cold]};
(* use DiracSimplify to move the GA[5] to the right. *)
*)
(* ---------------- *)
(* pole/.rule1/.rule4/.rule2/.rule3/.rule4/.rule5//Expand//Total//Total//Contract//SUNSimplify//QSimplify *)
pole=pole/.simp/.sun/.expand/.seprate/.factor//Contract//SUNSimplify;
pole=pole//Total//Expand;


(* ---------------- *)
pole=pole/.{vl[ca_,n_Integer|n_Complex,cb_]:>n vl[ca,1,cb],vr[ca_,n_Integer|n_Complex,cb_]:>n vr[ca,1,cb]};
pole=pole/.{vl[ca_,Pair[LorentzIndex[lla_,dim___],LorentzIndex[llb_,dim___]],cb_]:>Pair[LorentzIndex[lla,dim],LorentzIndex[llb,dim]]vl[ca,1,cb],
			vr[ca_,Pair[LorentzIndex[lla_,dim___],LorentzIndex[llb_,dim___]],cb_]:>Pair[LorentzIndex[lla,dim],LorentzIndex[llb,dim]]vr[ca,1,cb]};
(* extract the vertex doesn't involving gamma-matirx of sun-matrix *)

pole=pole/.{vl[-v1_]:>-vl[v1],vr[-v2_]:>-vr[v2]};
pole=Expand[pole]/.aa_ v1_vl v2_vr:>Contract[aa] verx[v1,v2];(* isolate the vertex part *)

pole=(pole//QSimplify)/.expand/.factor/.CA-2CF->1/CA;
pole=(pole/.verx[-v1_vl,v2_]:>-verx[v1,v2]/.verx[v1_,-v2_vr]:>-verx[v1,v2])//QSimplify;
(* it's possible generate terms like vl[-a]vr[b] and vl[a]vr[b] after first acting QSimplify; act QSimplify twice *)
pole=pole/.aa_ vvs_Plus/;!FreeQ[vvs,_vl|_vr]:>Plus@@(aa List@@vvs);(* Expanding terms like (a verx[_,_] + b verx[_,_]) *)

(* the terms in list will be have different vertices *)
pole=If[Head[pole]===Plus,List@@pole,{pole}]/.verx->Times;


(*(* replace the dummy indices added by algorithem to the set of dummy indices in input. *)
pole=(#/.vl[cca_,_,ccb_]vr[ccc_,_,ccd_]\[RuleDelayed](ruledummy=(Times@@DeleteDuplicates[{cca,ccb,ccc,ccd}])/(Times@@dummy);1);If[ruledummy===1,#,ruledummy=ruledummy/.aa_ bb_ Power[cc_,-1]Power[dd_,-1]\[RuleDelayed]{aa\[Rule]cc,bb\[Rule]dd}/.aa_ Power[bb_,-1]\[RuleDelayed]{aa\[Rule]bb};#/.ruledummy])&/@pole;*)

(* show as a table *)
pole=({(*(#/.{aa_vl->1,bb_vr->1})*)FCReplaceD[#/.{aa_vl->1,bb_vr->1},D->4],#/.aa_ vll_vl:>List@@vll,#/.aa_ vrr_vr:>List@@vrr}&/@pole)//QNormal;


(* after replace dummy indices, there are many terms can be combined to a single term. Combine the term with same vertices and same contraction pattern of color indices.*)
pole=Gather[pole,({Position[{#1[[3,1]], #1[[3,3]]}, #1[[2,1]]], Position[{#1[[3,1]], #1[[3,3]]}, #1[[2,3]]]}==={Position[{#2[[3,1]], #2[[3,3]]}, #2[[2,1]]], Position[{#2[[3,1]], #2[[3,3]]}, #2[[2,3]]]} && QSimplify[vl[#1[[2,2]] ]vr[#1[[3,2]] ] - vl[#2[[2,2]] ] vr[#2[[3,2]] ]]===0 )&];


(* combine the same term *)
pole=If[Length[#]==1,
	#[[1]]
,
	{Simplify[Plus@@(#[[1]]&/@#)],#[[1,2]],#[[1,3]]}
]&/@pole;


pole=ChangeDimension[pole,D];(* show the gamma-matrices as D-dimensional *)

If[ToLowerCase[ToString[astb]]=="true",
	pole
,
	pole=DeleteCases[pole,{0,0,0}];
	Plus@@(#[[1]]FourquarkCurrent[Current[f11,#[[2,1]],#[[2,2]],f12,#[[2,3]]],Current[-1,f21,#[[3,1]],#[[3,2]],-1,f22,#[[3,3]]]]&/@pole)
]


]/;Length[cur[[2]]]>=6


End[]
(*EndPackage[]*)
