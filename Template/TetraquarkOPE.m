(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


TetraquarkOPE::usage="TetraquarkLeading[q_,j1_,j2_] gives OPE of <j1 j2>, with j1, j2 are tetraquark currents. "
TetraquarkOPE::inderr="Dummy indices conflict!"
TetraquarkOPE::ipterr="Unknow expression!"
TetraquarkOPE::curerr="Unknow current structure!"
TetraquarkOPE::cerr="Can't recongize the assigned condensates"
TetraquarkOPE::mwarn="Momentum not match, momentum `1` in the files will be replaced to `2`"
TetraquarkOPE::merr="Fail to specify the momentum"


Begin["`Private`TetraquarkOPE`"]

Options[TetraquarkOPE]={
	NLO->True,
	Renormalization->True,
	Parallelized->True,
	Strategy->"Fourier",
	AutoNDR->True,
	Condensates->All,
	HoldFlavor->False,
	EpsOrder->0,
	TypeTags->"None"
}



TetraquarkOPE[qq_:"auto",factor1_ current1_FourquarkCurrent,factor2_ current2_FourquarkCurrent,ops:OptionsPattern[]]:=factor1 factor2 TetraquarkOPE[qq,current1,current2,ops]/;(FreeQ[factor1,Current]&&FreeQ[factor2,Current])



(*-----------------------------------------------*)
TetraquarkOPE[qq_:"auto",cur1_FourquarkCurrent,cur2_FourquarkCurrent,ops:OptionsPattern[]]:=Block[{null,tmpc1,tmpc2,tmp1,tmp2,tmp,trs,dia,x,q,v11=1,v12=1,c11,c12,c13,c14,f11,f12,f13,f14,v21=1,v22=1,c21,c22,c23,c24,f21,f22,f23,f24,hv21,hv22,ndr0,ndr1,ndr=OptionValue[AutoNDR],
tmp0,nlo=True,lofact=1,nlofact=1,renorm=True,condensates=OptionValue[Condensates],con1=0,con1a=0,conmqq=0,congg=0,conqq2=0,conmqgq=0,conggg=0,cond81=0,cond82=0,cond83=0,cond84=0,mqqgg=0,cond100=0,cond101=0,cond102=0,cond103=0,condmqq3=0,opslist,parallelized,pall,fdir,holdf=OptionValue[HoldFlavor],tmpq,conds,types,d01,d02,d03,tags},

If[qq==="auto"&&!DirectoryQ[OptionValue[Parallelized]//ToString],Message[TetraquarkOPE::merr];Abort[]];
(* momentum must be specified, except when load the date from a directory *)


(*--- get the falovrs, vertices, and color indices ---*)
tmpc1=cur1/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f11=aa;c11=bb;f12=cc;c12=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f11=aa;c11=bb;v11=vv;f12=cc;c12=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f13=aa;c13=bb;f14=cc;c14=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f13=aa;c13=bb;v12=vv;f14=cc;c14=dd;1)};


tmpc2=cur2/.{Current[{aa_,bb_},{cc_,dd_}]|Current[aa_,bb_,cc_,dd_]:>(f21=aa;c21=bb;f22=cc;c22=dd;1),Current[{aa_,bb_},vv_,{cc_,dd_}]|Current[aa_,bb_,vv_,cc_,dd_]:>(f21=aa;c21=bb;v21=vv;f22=cc;c22=dd;1),
					Current[{-1,aa_,bb_},{-1,cc_,dd_}]|Current[-1,aa_,bb_,-1,cc_,dd_]:>(f23=aa;c23=bb;f24=cc;c24=dd;1),Current[{-1,aa_,bb_},vv_,{-1,cc_,dd_}]|Current[-1,aa_,bb_,vv_,-1,cc_,dd_]:>(f23=aa;c23=bb;v22=vv;f24=cc;c24=dd;1)};


If[!FreeQ[{tmpc1,tmpc2},Current],
	Message[TetraquarkLeading::curerr];
	Abort[]
];


(*--- check dummy indices in input ---*)
tmp2=Length[#]&/@Gather[({c11,c12,c13,c14,c21,c22,c23,c24}//FCI)/.SUNFIndex[aa_]:>aa];

If[DummyCheck[v11,v12,v21,v22]==0||!FreeQ[tmp2,nn_Integer/;nn>2],
	Message[TetraquarkOPE::inderr];
	Abort[]
];


(*---------------- options & setting ------------------*)

Which[
	OptionValue[NLO]===False, nlofact=0,
	OptionValue[NLO]==="NLOonly", lofact=0
	];



If[OptionValue[Renormalization]=!=True,renorm=False];


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
		Message[TetraquarkOPE::cerr];
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

(* to do: allow add tags for each type of condensate *)
If[OptionValue[TypeTags],
	tags=True
,
	tags=False
];

(*----------------------------------------------------------*)

tmp={con1 lofact Hold[TetraquarkLeading[qq,cur1,cur2,AutoNDR->ndr0,Parallelized->parallelized,HoldFlavor->holdf,EpsOrder->OptionValue[EpsOrder]]],
	con1 nlofact Hold[TetraquarkAlphas[qq,cur1,cur2,AutoNDR->ndr1,Renormalization->renorm,Parallelized->parallelized,HoldFlavor->holdf,AlphasTypes->types]],
	congg Hold[Tetraquarkgg[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	conmqq Hold[Tetraquarkmqq[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	conggg Hold[Tetraquarkggg[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf,TypeTags->tags]],
	conqq2 Hold[Tetraquarkqq2[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond81 Hold[Tetraquarkd8type1[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond82 Hold[Tetraquarkd8type2[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond83 Hold[Tetraquarkd8type3[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond84 Hold[Tetraquarkd8type4[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	conmqgq Hold[Tetraquarkmqgq[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	mqqgg Hold[Tetraquarkmqqgg[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond100 Hold[Tetraquarkd10type0[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond101 Hold[Tetraquarkd10type1[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond102 Hold[Tetraquarkd10type2[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	cond103 Hold[Tetraquarkd10type3[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]],
	condmqq3 Hold[Tetraquarkmqq3[qq,cur1,cur2,Parallelized->parallelized,AutoNDR->ndr,HoldFlavor->holdf]]};
	

(*----------------------------------------------------------*)

If[pall===True,
		
	tmp=Plus@@#&/@WaitAll[ReleaseHold[tmp/.parallelized->"External"]];
	Plus@@((QGather[(#/.{CA->SUNN,CF->(SUNN^2-1)/(2CA)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0})&/@tmp)
	
				
,		
	

	If[pall==="External",
	(* !!! do not touch tmp before WaitAll[tmp], any function may modifiy the things in EvaluationObject[...] should not be applied here. !!! *)
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
	
			
	,(* If not parallelized *)
		

		tmp=ReleaseHold[tmp/.parallelized->False];
		Plus@@((QGather[(#/.{CA->SUNN,CF->(SUNN^2-1)/(2CA)})//Expand//SUNSimplify,qq,ShowasTable->False]/.CA-2CF->1/CA/.{-1+CA^2-2CA CF->0,1-CA^2+2CA CF->0})&/@tmp)
		
	]

]
]



End[]
(*EndPackage[]*)
