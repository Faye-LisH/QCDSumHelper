(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



QCDSumDiagram::usage="Generate the OPE diagrams for the correlator"
QCDSumDiagram::ierr="Dummy indices conflict"
QCDSumDiagram::ferr="The number of Fermions must equal to the number of Anti-Fermions"
QCDSumDiagram::gerr="The number of Ghosts must equal to the number of Anti-Ghosts"
QCDSumDiagram::cerr="Unrecogized condensates detected, each condensate should specified as {nq_,ng_,nm_}, Condensate[{nq_,ng_,nm_}], or product of basic Condensate[{nq_,ng_,nm_}]"
QCDSumDiagram::xerr="The coordinates of the operators conflict"
QCDSumDiagram::olisterr="The input operators should be a list of operators or a list of two sublists of operators."
QCDSumDiagram::owan="The lowest-order connected diagram has order higher than the specified order."
QCDSumDiagram::dcwan="Duplicated condensates will be ignored!"
(*QCDSumDiagram::derwan="The covariant derivatives acting on the fields in the operator will be symmetrized, for the case they are antisymmetrized, written them as gluonfieldstrength."*)
QCDSumDiagram::cderwan="The covariant derivatives acting on the fields will not be expanded for the case where the gluonfieldstrength is involved in the condensate and the diagram gives an nonlocal-condensate."
QCDSumDiagram::loopwan="(Sub)diagram that cannot integrated to a propagator involved, the automatic evaluation stoped."


Begin["`Private`QCDSumDiagram`"]	
Options[QCDSumDiagram]={
	Condensates->"d0",
	Draw->False,
	Evaluate->False,
	DiscardTadPole->"MassLess",
	ActiveFlavors->{"u","d","s"},
	Factorization->"Auto",
	MassOrder->1,
	Massive->False,
	AutomaticOrder->False,
	NonLocalCondensate->False,
	Order->0
}
(* Condensates -> {cond_1, cond_2, ...}, for the case that cond_i = Condensate[{...}] Condensate[{...}]... and each Condensate[{...}] is a basic condensate, the NonLocalCondensate will always be treated as False for the digram related to cond_i. For cond_i = Condensate[{...}], the order is given by Optrion[Order], for cond_i = {n_, Condensate[{...}]}, the diagrams corresponding to the order through 0 to n are generated, for cond_i = {{n}, Condensate[{...}]}, only the diagrams with order n are generated *)


(* for the input fields {o1, o2, ..., on}, identify it as correlator on-> o1+o2+... ;
for input fields {{o1_, o2, ...}, {oa_, ob_, ...}}, identify at as correlator oa+ob+... -> o1+o2+... *)
QCDSumDiagram[{op0_,op__},ops:OptionsPattern[]]:=Block[{tmp,tmpa,tmpb,null,ncompositeo,factor,olist,order=OptionValue[Order],pertm=0,z,tmpf,vertexlist={},condlist,fiedslist,pert,pertord,nonlocalc=OptionValue[NonLocalCondensate],list,dialist,in,out},

(* check position conflict *)
tmp = {op0,op}//Flatten;
If[Or@@Flatten[Table[!FreeQ[Join[tmp[[;;i-1]],tmp[[i+1;;]]],#]&/@DeleteDuplicates[Cases[{tmp[[i]]},Operator[xx_,__]:>xx,Infinity]//Flatten],{i,1,Length[tmp]}]],
	Message[QCDSumDiagram::xerr];
	Abort[]
];
(* each element in tmp should has unique position *)
If[Or@@((Length[DeleteDuplicates[Cases[#,Operator[xx_,__]:>xx,Infinity]]]>1)&/@tmp),
	Message[QCDSumDiagram::xerr];
	Abort[]
];

(* expand composite operators *)
tmp=null {op0,op}//Expand;

If[Or@@(MatchQ[#,_Plus]&/@tmp),
(*  for the case on-> o1+o2+... *)
	ncompositeo=0;
	tmp=tmp//.{aa___,bb_Plus,cc___}/;!FreeQ[bb, dd_ _Operator]:>(ncompositeo=ncompositeo+1;{aa,#,cc}&/@(List@@bb));
	QCDSumDiagram[#,ops]&/@(Flatten[tmp,ncompositeo-1]/.null->1)

,

	Which[!(MatchQ[tmp,{_List,_List}]||MatchQ[tmp,_List]),
		Message[QCDSumDiagram::olisterr];
		Abort[]
	,
		(* for the case involcing Operator[...] + Operator[...] *)
		MatchQ[tmp,{_List,_List}]&&Or[Or@@(MatchQ[#,_Plus]&/@tmp[[1]]), Or@@(MatchQ[#,_Plus]&/@tmp[[2]])],
		(* for the case oa+ob+... -> o1+o2+... *)
		ncompositeo=0;
		tmpa = tmp[[1]]//.{aa___,bb_Plus,cc___}/;!FreeQ[bb, dd_ _Operator]:>(ncompositeo=ncompositeo+1;{aa,#,cc}&/@(List@@bb));
		If[ncompositeo>0,tmpa=Flatten[tmpa,ncompositeo-1]/.null->1,tmpa={tmpa}/.null->1];
		ncompositeo=0;
		tmpb = tmp[[2]]//.{aa___,bb_Plus,cc___}/;!FreeQ[bb, dd_ _Operator]:>(ncompositeo=ncompositeo+1;{aa,#,cc}&/@(List@@bb));
		If[ncompositeo>0,tmpb=Flatten[tmpb,ncompositeo-1]/.null->1,tmpb={tmpb}/.null->1];

		Flatten[Outer[QCDSumDiagram[{#1,#2},ops]&,tmpa,tmpb,1],1]

	,
		True,
		(* each term in the list is a single operator *)
		factor=tmp/._Operator->1;(* all prefactors *)
		(* list of operators *)
		olist=(tmp/factor)/.Operator[xx_,vv_,os_List]:>Operator[xx,vv//DiracChainJoin//SUNSimplify,os];(* join the \[Gamma]-matrices and sun-matrices if possible *)
		factor=Times@@Flatten[factor/.null->1];

		(*-----------------*)
		(* rename the dummy indices *)
		olist=dummyRename[olist];
		olist=olist/.{Fermion[{f_},other__]:>Fermion[f,other],AntiFermion[{f_},other__]:>AntiFermion[f,other]};(* unify the form of (Anti)Fermion *)
		If[Length[olist]<2,Message[QCDSumDiagram::olisterr];Abort[]];
		(*--------------------------------------------------*)

		(* the number of fermions must = the number of anti-fermions *)
		tmp={Cases[olist,_Fermion,Infinity],Cases[olist,_AntiFermion,Infinity]};
		If[Length[tmp[[1]]]=!=Length[tmp[[2]]],Message[QCDSumDiagram::ferr];Abort[]];


		(* the number of ghosts must = the number of anti-ghosts, if they are manually introduced in the operators *)
		tmp={Cases[olist,_Ghost,Infinity],Cases[olist,_AntiGhost,Infinity]};
		If[Length[tmp[[1]]]=!=Length[tmp[[2]]],Message[QCDSumDiagram::gerr];Abort[]];

		(*----------------------------------------------------------------------*)
		condlist = OptionValue[Condensates];
		If[!FreeQ[condlist,"d0"],
			pert=True;
			condlist=Replace[condlist,{"d0"->(pertord=order;1),{ords_,"d0"}:>(pertord=ords;1)},1];
			condlist=DeleteCases[condlist,1]
		,
			If[!FreeQ[condlist,"m"],
				pert=True;
				pertm=1;
				condlist=Replace[condlist,{"m"->(pertord=order;1),{ords_,"m"}:>(pertord=ords;1)},1];
				condlist=DeleteCases[condlist,1]
			,
				pert=False
			]
		];
		(* unify the form of condensates *)
		condlist=Replace[condlist,{str_String:>{str},{nor_,str_String}:>{{nor,str}},{nq_Integer,ng_Integer,nm___Integer}:>{Condensate[{nq,ng,nm}]},{nor_,{nq_Integer,ng_Integer,nm___Integer}}:>{{nor,Condensate[{nq,ng,nm}]}}},{0}];
		condlist=condlist/.{"gg"->Condensate[{0,2,0}],"qq"->Condensate[{2,0,0}],"mqq"->Condensate[{2,0,1}],"qq2"->Condensate[{4,0,0}],"qgq"->Condensate[{2,1,0}],"ggg"->Condensate[{0,3,0}],"qq3"->Condensate[{6,0,0}]};

		condlist=Replace[condlist,{{nq_,ng:Except[_List|_Condensate|_NonLocalCondensate],nm___}:>Condensate[{nq,ng,nm}],{no_,{nq_,ng:Except[_List|_Condensate|_NonLocalCondensate],nm___}}:>{no,Condensate[{nq,ng,nm}]}},{1}];

		condlist=condlist/.{Condensate[nn_Integer]:>Condensate[{nn,0,0}],Condensate[{nn1_,nn2_}]:>Condensate[{nn1,nn2,0}],NonLocalCondensate[clist_List,c2list___]:>NonLocalCondensate[0,clist,c2list]};
		(* combine condensates with duplicated spedified order *)
		condlist=condlist//.{{aa___,{n1_Integer,cond_Condensate},bb___,{n2_Integer,cond_Condensate},cc___}:>{aa,{Max[n1,n2],cond},bb,cc},{aa___,{n1_Integer,cond_Condensate},bb___,{n2_List,cond_Condensate},cc___}/;n2[[1]]<=n1:>{aa,{n1,cond},bb,cc},{aa___,{n1_List,cond_Condensate},bb___,{n2_Integer,cond_Condensate},cc___}/;n1[[1]]<=n2:>{aa,{n2,cond},bb,cc}};
		condlist=condlist//.{{aa___,{n1_Integer,cond_NonLocalCondensate},bb___,{n2_Integer,cond_NonLocalCondensate},cc___}:>{aa,{Max[n1,n2],cond},bb,cc},{aa___,{n1_Integer,cond_NonLocalCondensate},bb___,{n2_List,cond_NonLocalCondensate},cc___}/;n2[[1]]<=n1:>{aa,{n1,cond},bb,cc},{aa___,{n1_List,cond_NonLocalCondensate},bb___,{n2_Integer,cond_NonLocalCondensate},cc___}/;n1[[1]]<=n2:>{aa,{n2,cond},bb,cc}};

		(* detect duplicated condensates *)
		If[!DuplicateFreeQ[condlist],
			Message[QCDSumDiagram::dcwan];
			condlist=DeleteDuplicates[condlist]
		];

		If[Length[DeleteCases[Replace[Replace[condlist,{_Integer|{_Integer},_Condensate|_NonLocalCondensate}:>1,1],_Condensate|_NonLocalCondensate->1,1],1]]>0,
			Message[QCDSumDiagram::cerr];
			Abort[]
		];


		(* Interpret the condensates into background fields vertices; if the order hasn't been specified, set it to the OptionValue[Order] *)
		condlist=If[MatchQ[#,_Condensate|_NonLocalCondensate],{order,condInterpretor[#,nonlocalc]},{#[[1]],condInterpretor[#[[2]],nonlocalc]}]&/@condlist;
		(* prepare the vertices for each Condensate; the operators will also be expanded; each element contains operators, condensates, and vertices, and the operators and condensates cannot be expanded *)
		If[pert&&pertm==0,
			condlist = Join[{{pertord,{{0,{}}}} },condlist]
		,
			If[pertm==1,
				condlist = Join[{{pertord,{{1,{}}}} },condlist]
			]
		];

		fiedslist=Flatten[diagramFieldsPrepare[olist,#[[2]],#[[1]]]&/@condlist]/.list->List;

		fiedslist=DeleteCases[fiedslist,ocv_/;ocv[[2]]==={-1}];(* delete the case vertices_list = {-1}, the -1 label means the lowest-order connected diagram has order higher than specified order *)

		fiedslist=DeleteCases[fiedslist,ocv_/;FreeQ[ocv[[1,1]],_Gluon|_GluonStrength|_Ghost|_AntiGhost]&&FreeQ[ocv[[2]],_Fermion|_AntiFermion]&&!FreeQ[ocv[[2]],_Gluon|_GluonStrength]];
		(* if the operators are pure fermions while the vertices are pure gluons, they cannot connected to diagrams *)
		fiedslist=DeleteCases[fiedslist,ocv_/;FreeQ[ocv[[1,1]],_Gluon|_GluonStrength]&&!FreeQ[ocv[[1,1]],_Ghost|_AntiGhost]&&FreeQ[ocv[[2]],_Fermion|_AntiFermion|_Ghost|_AntiGhost]];
		fiedslist=DeleteCases[fiedslist,ocv_/;FreeQ[ocv[[1,1]],_Ghost|_AntiGhost]&&OddQ[Length[Cases[ocv[[2]],_Ghost,Infinity]]]];

		(* the positions for in states and out states *)
		If[MatchQ[olist,{_List,_List}],
			in = Cases[olist[[1]],Operator[xx_,__]:>xx,Infinity]//DeleteDuplicates;
			out = Cases[olist[[2]],Operator[xx_,__]:>xx,Infinity]//DeleteDuplicates
		,
			in ={ olist[[-1]][[1]]};
			out =Cases[olist[[;;-2]],Operator[xx_,__]:>xx,Infinity]//DeleteDuplicates
		];
		(*Print[">>>",fiedslist];*)
		(* generate diagram for each case *)
		(* dialist=diagramGenerator[#,in,out]&/@fiedslist; *)
		
		
		fiedslist
		(* each element looks like {factors, {fields at x, fields at y, ...}, vertices, fields in condensate, {condensate}} or {factors, {fields at x, fields at y, ...}, vertices, fields in nonlocal-condensate, {nonlocal-condensate, condensate}} *)


	]
]
](*/;And@@(MatchQ[#,_Operator|factor_ o_Operator]&/@Flatten[{op0,op}])*)


(* generate diagrams for each Condensate *)
diagramFieldsPrepare[os_List,conds_List,order_]:=Block[{tmp,in,out,oplist,tmpo,vertices,condnotation,list,listcon,gluonnot,fermionnot,tmp2},
(* conds = {{0,{}}} or {{1,{}}} for perturbative diagram *)
vertices=(#/.{{nm_,clist__List}/;Length[{clist}[[1]]]>0:>{nm,Table[CondtoVertex[#[[2]]],#[[1]]]&/@{clist}},{nm_,nc_NonLocalCondensate}:>{nm,CondtoVertex[nc]}})&/@conds;
(*condnotation=Replace[conds,{{nm_,clist__List}:>quarkMass[]^nm Times@@Flatten[#[[2]]^#[[1]]&/@{clist}],{nm_,nc_NonLocalCondensate}:>quarkMass[]^nm nc},1];*)
vertices=vertices/.{{ncv_nlocalconVertex}}:>ncv;(* unify the form for e.g. the nlocalconVertex generated for Condensate[{2,0,1}] *)
(* each term has the form: {order, {{n_mass in propagators, {n_cond, condensate_1}, ...}, {n_mass in propagators, {condensate_1 vertices}, {condensate_2 vertices}, ...}}} *)
tmp=listcon[order,#]&/@({conds,vertices}//Transpose);


oplist=Replace[Flatten[os],Times->List,{2},Heads->True]//Flatten;
oplist=Flatten[OperatorFactorizor[#]&/@oplist];(* flattened operator list *)

(* expand the covairant derivaives in the operators if no nonlocalcondensate involved, 
or, 1) covariant derivarive act on fermions but no fermions in the nonlocalcondensate
2) covariant derivative act on gluon/gluonstrength but no gluon/gluonstrength in the nonlocalcondensate *)
tmp=(gluonnot=!FreeQ[#,ncon_NonLocalCondensate/;!FreeQ[ncon,_Gluon|_GluonStrength]];
		fermionnot=!FreeQ[#,ncon_NonLocalCondensate/;!FreeQ[ncon,_Fermion|_Antifermion]];

(*tmpo=list@@(CDExpand[#,fermionnot,gluonnot]&/@oplist);*)
		tmpo=list@@(CDExpand[#,False,False]&/@oplist);(* for D^u\psi, since A^(0) = 0, expand covariant derivative will not causing problem *)
(* expand composite fields, do not commute any operators *)
		tmpo=tmpo//.list[aa___,bb_List,cc___]:>{list[aa,#,cc]&/@bb};
		tmpo=Flatten[tmpo]/.list->List;
(* to a list of {operator_List,{condensate_List, condensates' fields list}} *)
		Table[list[tmpo[[i]],#],{i,1,Length[tmpo]}]
	)&/@tmp;

(* recall that order = {n} means g^n, order = n means all g^i with i<=n *)
tmp=tmp/.listcon[n_Integer/;n>0,condver_List]:>Table[listcon[i,condver],{i,0,n}]/.listcon[{n_Integer},condver_List]:>listcon[n,condver];
tmp=tmp/.list[op_List,conlist_List]:>(list[op,#]&/@conlist);
(* order = n means g^n from now on *)
tmp=Flatten[tmp]/.listcon[n_,{con_List,ver_List}]:>{Join[{n},con],Join[{n},ver]};
(* each term has the form {list of expaned operators, {{n_order, n_mass, {n_condensate, Condensate}, ...}, {n_order, n_mass, list of condensate vertices}} } *)


(* The field involving covariant derivative should be carefully treated, for D^uD^v\psi, expand D^u = d^u -igA^u yields A^u d^v\psi, d^ud^v\psi, ..., but if the Condensate <\bar{\psi}D^uD^v\psi> is involved, the expanded d^u and A^u should be collected, alternatively, writting d^ud^b\psi as d^ud^v(1+ x^a D_a + 1/2 x^a x^b D_a D_b+ ...)\psi(0)|x=0 and d^uA^v = d^u(1/2 x_aG^av + ...) |x=0 D^uD^v\psi yields D{uD^v}\psi - 1/2 ig G^uv\psi, which just (anti)symmetrize the original D^uD^v\psi, but the later term equivalent to using equation of motion before factorization. one can write gluon field as A^u+B^u where A^u is dynamic/perturbative while B^u is background gluon, so A^u only contributes to wilson coefficient while B^u only contributes to condensates, similarly for fermion field. Add tags in (Anti)Fermion, GLuon(Strength) to label whether the fields are dynamic or background, and the covariant derivative can be safely expandee as (d^u -igB^u) - igA^u, the d^uA^v is dynamic so it not affect thecondensate, while the B^u still contain in D^u so no G^uv will appear in the expansion. *)
(* to do: 
1) for each field, add option Contract -> "Dynamic"/"Background", the former can only connect to propagator while the later can only connect to Condensate
2) such expansion of covariant derivative may yield (D^u)^ab (D^v\psi)^c (T^n)^bc A^n, store longer color indices and derivative structure are required for the fields
3) for CDExpand, derAAExpand, derexpand, the expansion becomes complicate *)


tmp=Flatten[tmp]/.list->List;
(* each term has the form: {operator list, {{order, n_mass, {number of condensates, cond}, ...}, {order, n_mass, {number of condensates vertices, Condensate vertices}, ...}} }*)
(* preparing the vertices needed for each diagram *)

tmp2=tmp/.{Power[gStrong,nn_]:>(list@@Table[gStrong,nn]),NonLocalCondensate[fields1_List,cond_List]:>NonLocalCondensate[0,fields1,cond]};(* prepare the g^n and condensates for preparing the vertices *)
(*Print[#[[1]],"--->>>",#[[2]]]&/@Transpose[{tmp,tmp2}];*)
If[conds==={{0,{}}}||conds==={{1,{}}},(* perturbative *)
	tmp2=((*Print[#,"<><><><>",Total[Cases[#[[1]],gStrong->1,Infinity]],">>>>",Total[Cases[#[[1]],GluonStrength[_List,Order->1]->1,Infinity]],"--->>",Total[Cases[#[[1]],GluonStrength[_List,Order->0]->1,Infinity]],"--->>",Total[Cases[#[[1]],GluonStrength[_List,Order->1]->2,Infinity]],"---->>",Total[Cases[#[[1]],_Gluon->1,Infinity]]];*)VerticesPrepare@@Join[{Total[Cases[#[[1]],gStrong->1,Infinity]]+Total[Cases[#[[1]],GluonStrength[_,_,1]->1,Infinity]],(* all gStrong in operators *)
			Total[Cases[#[[1]],_Gluon->1,Infinity]]+Total[Cases[#[[1]],GluonStrength[_,_,0]->1,Infinity]]+Total[Cases[#[[1]],GluonStrength[_,_,1]->2,Infinity]],(* all gluons in operators *)
			#[[2,1,1]]},
			{{},{}}
			]
		)&/@tmp2
,
	tmp2=((*Print[#,"<><><><>",Total[Cases[#[[1]],gStrong->1,Infinity]],">>>>",Total[Cases[#[[1]],GluonStrength[_List,Order->1]->1,Infinity]],"--->>",Total[Cases[#[[1]],GluonStrength[_List,Order->0]->1,Infinity]],"--->>",Total[Cases[#[[1]],GluonStrength[_List,Order->1]->2,Infinity]],"---->>",Total[Cases[#[[1]],_Gluon->1,Infinity]]];*)VerticesPrepare@@Join[{Total[Cases[#[[1]],gStrong->1,Infinity]]+Total[Cases[#[[1]],GluonStrength[_,_,1]->1,Infinity]],
			Total[Cases[#[[1]],_Gluon->1,Infinity]]+Total[Cases[#[[1]],GluonStrength[_,_,0]->1,Infinity]]+Total[Cases[#[[1]],GluonStrength[_,_,1]->2,Infinity]],
			#[[2,1,1]]},
			If[MatchQ[#[[2,1,3]],_List],{#[[2,1,3;;]],#[[2,1,3;;]]},{#[[2,1,3,2]],#[[2,1,3,3]]}]
			]
		)&/@tmp2
];

(* for each case in tmp, there are multipule choice of verteices combination *)
(* restore list -> List at QCDSumDiagram[] *)

tmp=Table[list[tmp[[i]],#]&/@tmp2[[i]],{i,1,Length[tmp]}]//Flatten

]





dummyRename[os_List]:=Block[{tmp,olist,inds,inds1,dummy,len1=0},

(* dummy indices in each operators *)
olist=os/.Operator[xx_,vv_,oo_List,lors___]:>(inds=Flatten[oo/.{Fermion[fla_,others__]:>{others},AntiFermion[fla_,others__]:>{others},Gluon->List,GluonStrength->List,Ghost->List,AntiGhost->List}];
inds1=Flatten[Cases[vv,#,Infinity]&/@DeleteDuplicates[inds]];
inds1=Gather[Join[inds,inds1]];

If[!FreeQ[inds1,is_List/;Length[is]>2,{1}],(* dummy indices appears more than twice *)Message[QCDSumDiagram::ierr];Abort[]];

(* the dummy indices *)
inds=Cases[inds1,{aa_,aa_}:>aa,{1}];
inds=(#->dummy[Unique[#]])&/@inds;(* rename dummy indices, add dummy[] head for later distinguishment *)
If[Or@@(!FreeQ[{lors},#[[1]]]&/@inds),
(* the (anti)symmetrized indices should has no intersection with the dummy indices *)
Message[QCDSumDiagram::ierr];
Abort[];
,
If[Length[inds]>0,Operator[xx,FCI[vv/.inds],oo/.inds,lors],Operator[xx,vv//FCI,oo,lors]]
]

);

(* unify the form of the cases {operators...} and {{operators...}, {operators...}} *)
If[MatchQ[olist,{_List,_List}],
len1= Length[olist[[1]]];
olist=Flatten[olist]
];

(* the same indices appear in the operators at different positions are treated as different indices *)
Do[
(* dummy indices for different operators in same position *)
inds=Cases[olist[[i]],#,Infinity]&/@{LorentzIndex[lor_,___]:>lor,SUNIndex[su_]:>su,SUNFIndex[su_]:>su,Fermion[fla_,others__]:>{others},AntiFermion[fla_,others__]:>{others},gg_Gluon:>List@@gg,gg_GluonStrength:>List@@gg,gh_Ghost:>List@@gh,gh_AntiGhost:>List@@gh};
inds=DeleteCases[Flatten[inds],_dummy];(* ignore the indices that have already been renamed *)

inds=Gather[inds];
If[!FreeQ[inds,is_List/;Length[is]>2,{1}],(* dummy indices appears more than twice *)Message[QCDSumDiagram::ierr];Abort[]];

(* the dummy indices *)
inds=Cases[inds,{aa_,aa_}:>aa,{1}];
inds=(#->Unique[#])&/@inds;(* rename dummy indices *)
olist[[i]]=olist[[i]]/.inds;
(*------------------------------*)
(* remove the dummy[] head *)
olist[[i]]=olist[[i]]/.dummy[ids_]:>ids
,
{i,1,Length[olist]}
];

If[len1>0,olist={olist[[;;len1]],olist[[len1+1;;]]}];

olist

]


(* e.g, meson-meson operator -> {meson operator, meson operator} *)
(* Operator[xx_,vv_,os_List]:>Operator[xx,vv//DiracChainJoin//SUNSimplify,os] has already applied before *)
(* for fermions, the dirac indices is prior, so the tetraquark operator -> {diaquark operator, anti-diaquark operator} *)
OperatorFactorizor[op_Operator]:=Module[{tmp,null,olist,vlist,opslist={},fmatchlist,flist,gmatch,gmatchlist,glist,i,tmpvertex,dindx1,dindx2,sindx1,sindx2,gindx,remainv,remaino,tmpo,sign},
vlist=Expand[op[[2]]];
olist=op[[3]];

If[MatchQ[vlist,_Plus],
(* for composited vertex *)
Flatten[OperatorFactorizor[Operator[op[[1]],#,op[[3]]]]&/@(List@@vlist)]/.null->1
,

vlist=DeleteCases[List@@FCI[null vlist//Expand],null];

For[i=1,i<Length[vlist]+1,i++,
(* 1) find the (anti)fermions with dirac indices match DCHN[GAD[...],i,j], 
2) obtain the su(n) indices of these two (anti)fermions, e.g. call it a, b 
3) if find SUNTF[n,a,b], and find G^n_uv or A^n, encompass them as { DCHN[GAD[...],i,j]SUNTF[n,a,b], a list of (anti)fermions and G^n_uv or A^n}  *)
(* to do (maybe): 1) factorize pure gluons operators 2) dealing with the case involving ghosts *)
tmpvertex =1;
flist={};
glist={};
If[MatchQ[vlist[[i]],_DiracChain],
(* obtain dirac indices *)

dindx1=vlist[[i,-2,1]];
dindx2=vlist[[i,-1,1]];

fmatchlist= Boole[!FreeQ[#,dindx1|dindx2]]&/@olist;
If[Total[fmatchlist]==2,
(* selected fermions *)
flist= DeleteCases[fmatchlist olist,0];

(* obtain the sign generated by commute fermions *)
tmpo=Position[fmatchlist,1];
tmpo=olist[[tmpo[[1,1]];;tmpo[[2,1]] ]];
sign=(-1)^Length[Cases[tmpo,_AntiFermion|_Fermion,Infinity]];
(* update olist *)
olist=(1-fmatchlist)olist;

(* get the vertex; update vlist *)
tmpvertex = sign vlist[[i]];
vlist[[i]]=1;


(* obtain su(n) fundamental indices *)
sindx1=flist[[1,3]];
sindx2=flist[[2,3]];

gmatch=Position[vlist,sun_/;!FreeQ[sun,sindx1]&&!FreeQ[sun,sindx2],1];
If[Length[gmatch]==1,
(* su(n) adjoint indices *)
gindx=Cases[vlist[[gmatch[[1,1]]]],SUNIndex[ii_]:>ii,Infinity];
tmpvertex=tmpvertex vlist[[gmatch[[1,1]] ]];
vlist[[gmatch[[1,1]]]]=1;

(* find G^n_uv or A^n with matching adjoint indices *)
gmatchlist=Total[#]&/@Outer[Boole[!FreeQ[#1,#2]]&,olist,gindx];
If[Total[gmatchlist]>0,
glist=DeleteCases[gmatchlist olist,0];

olist=(1-gmatchlist)olist
]
]
]
];

If[tmpvertex=!=1,
opslist=Append[opslist,{tmpvertex,Join[flist,glist]}]
]
];


(* if the remaining olist contains fields, construct it to an operator with the remaining of vlist, otherwise append the vertex into the last operator in the opslist *)
remainv=Times@@vlist;
remaino=DeleteCases[olist,0];
If[Length[remaino]>0,
opslist=Append[opslist,{remainv,remaino}]
,
If[remainv=!=1&&Length[remaino]==0,
opslist[[-1]]={opslist[[-1,1]]remainv,opslist[[-1,2]]}
]
];


Operator[op[[1]],#[[1]],#[[2]]]&/@opslist
]
]


(* expand covariant derivative, also G^nv -> d^uA^v - d^vA^u  + g f A^u A^v *)
CDExpand[op_Operator,fermionnot_:False,gluonnot_:False]:=Block[{tmp,derlist,antiFermion,fermion,antiGhost,ghost,gluon,gLuon,gLuonS,gLuonSL,gLuonSD},
tmp={op};

(* keep covariant derivative act on (anti)fermions not expanded if fermionnot = True *)
If[!fermionnot,
tmp=tmp//.Operator[xx_,vv_,{os1___,AntiFermion[fla_,diracindx_,sunindx_,der_List],os2___}]:>(derlist=derexpand[der,True,False,sunindx];If[Length[#]==2,Operator[xx,vv,{os1,antiFermion[fla,diracindx,sunindx,#[[1]]],os2}],Operator[xx,vv #[[2,2]], Join[{os1,antiFermion[fla,diracindx,#[[2,1]],#[[1]]]},#[[3]],{os2}]]]&/@derlist);

tmp=tmp//.Operator[xx_,vv_,{os1___,Fermion[fla_,diracindx_,sunindx_,der_List],os2___}]:>(derlist=derexpand[der,False,False,sunindx];If[Length[#]==2,Operator[xx,vv,{os1,fermion[fla,diracindx,sunindx,#[[2]]],os2}],Operator[xx,vv #[[2,1]], Join[{os1,fermion[fla,diracindx,#[[2,2]],#[[3]]]},#[[1]],{os2}]]]&/@derlist)
];

tmp=tmp//.Operator[xx_,vv_,{os1___,AntiGhost[sunindx_,der_List],os2___}]:>(derlist=derexpand[der,True,True,sunindx];If[Length[#]==2,Operator[xx,vv,{os1,antiGhost[sunindx,#[[1]]],os2}],Operator[xx,vv #[[2,2]], Join[{os1,antiGhost[#[[2,1]],#[[1]]]},#[[3]],{os2}]]]&/@derlist);

tmp=tmp//.Operator[xx_,vv_,{os1___,Ghost[sunindx_,der_List],os2___}]:>(derlist=derexpand[der,False,True,sunindx];If[Length[#]==2,Operator[xx,vv,{os1,ghost[sunindx,#[[2]]],os2}],Operator[xx,vv #[[2,1]], Join[{os1,ghost[#[[2,2]],#[[3]]]},#[[1]],{os2}]]]&/@derlist);


(* keep covariant derivative act on gluon(strength) not expanded if gluonnot=True *)
If[!gluonnot,

tmp=tmp//.Operator[xx_,vv_,{os1___,Gluon[lorindx_,sunindx_,der_List],os2___}]:>(derlist=derexpand[der,False,True,sunindx];If[Length[#]==2,Operator[xx,vv,{os1,gLuon[lorindx,sunindx,#[[2]]],os2}],Operator[xx,vv #[[2,1]], Join[{os1,gLuon[lorindx,#[[2,2]],#[[3]]]},#[[1]],{os2}]]]&/@derlist);


(* for G^uv without derivative, set d^a A^b - d^b A^a as GluonStrength[..., Order->0] and gf A^b A^c as GluonStrength[..., Order->1] to reduce independent operators *)

tmp=tmp/.Operator[xx_,vv_,{os1___,GluonStrength[{lor1_,lor2_,sunindx_},{},nn_Integer],os2___}]:>Operator[xx,vv,{os1,gLuonS[{lor1,lor2,sunindx},nn],os2}];
tmp=tmp//.Operator[xx_,vv_,{os1___,GluonStrength[{lor1_,lor2_,sunindx_},{},opts_],os2___}]:>{Operator[xx,vv,{os1,gLuonS[{lor1,lor2,sunindx},0],os2}],Operator[xx,vv,{os1,gLuonS[{lor1,lor2,sunindx},1],os2}]};

(* for DG^uv, first Expand the D in DG^uv *)
tmp=tmp//.Operator[xx_,vv_,{os1___,GluonStrength[{lor1_,lor2_,sunindx_},der_List,opts_],os2___}]:>(derlist=derexpand[der,False,True,sunindx];If[Length[#]==2,Operator[xx,vv,{os1,gLuonSD[{lor1,lor2,sunindx},#[[2]],opts],os2}],Operator[xx,vv #[[2,1]], Join[{os1,gLuonSD[{lor1,lor2,#[[2,2]]},#[[3]] ,opts]},#[[1]],{os2}]]]&/@derlist);

(* then expand d^u... (d^a A^b - d^b A^a + g f A^b A^c) *)
tmp=tmp//.Operator[xx_,vv_,{os1___,gLuonSD[{lor1_,lor2_,sunindx_},der_List,0],os2___}]:>Operator[xx,vv,{os1,gLuonSL[{lor1,lor2,sunindx},der,0],os2}];

tmp=tmp//.Operator[xx_,vv_,{os1___,gLuonSD[{lor1_,lor2_,sunindx_},der_List,1],os2___}]:>(derlist=derAAExpand[{lor1,lor2,sunindx},der];Operator[xx,vv gStrong SUNF[sunindx,#[[1,2]],#[[2,2]]], Join[{os1},#,{os2}]]&/@derlist);

tmp=tmp//.Operator[xx_,vv_,{os1___,gLuonSD[{lor1_,lor2_,sunindx_},der_List,opts_],os2___}]:>(derlist=derAAExpand[{lor1,lor2,sunindx},der];{Operator[xx,vv,{os1,gLuonSL[{lor1,lor2,sunindx},der,0],os2}],Operator[xx,vv gStrong SUNF[sunindx,#[[1,2]],#[[2,2]]], Join[{os1},#,{os2}]]&/@derlist})
];


(* restore the form *)(*/.gLuonSL[lsu_List,{},Order->0]:>gLuonSL[lsu,Order->0]*)
Flatten[tmp/.{gLuonS[ind_List,nn_]:>GluonStrength[ind,Order->nn],gLuonSL[ind_List,der_List,nn_]:>GluonStrength[ind,der,Order->nn],gluon->Gluon,gLuon->Gluon,antiGhost->AntiGhost,ghost->Ghost,fermion->Fermion,antiFermion->AntiFermion}]

]


(* expand d^ud^v ... (g f^abc A^b_lor1 A^c_lor2) *)
derAAExpand[{lor1_,lor2_,sunindx_},der_List]:=Block[{tmp,gluon,list,col1=Unique["$Col"],col2=Unique["$Col"]},
tmp={list@@Join[Flatten[der],{gluon[lor1,col1,{}],gluon[lor2,col2,{}]}]};
tmp=tmp//.list[ders___,lora_,gl1_,gl2_]:>{list[ders,gluon[gl1[[1]],gl1[[2]],Append[gl1[[3]],{lora}]],gl2],list[ders,gl1,gluon[gl2[[1]],gl2[[2]],Append[gl2[[3]],{lora}]]]};
Flatten[tmp]/.gluon[inds_List,{}]:>gluon[inds]/.gluon->Gluon/.list->List
]


(* the covariant derivative should not be expandend if the gluon in covraiant derivative is not dynamic, 
i.e., for example, for <\bar{d}DDd> -> <\bar{s}Gs>, the D^u should not expanded to d^u -ig A^u, where d^u denote partial derivative,
D^uD^v s = (d^ud^v d + d^uA^v + ...) s(0), identifying d^u s as d^u s(x)|x->0, then (d^ud^v s(x)|x->0) -> D^{uD^v} s, 
the term d^uA^v s  gives 1/2 G^uv s, the expansion process equivalent to symmetrize and antisymmetrize the indices, but write D^uD^v to different form before or after factorization yields different recults, the factorization ambiguity involved *)
(* left=True if the derivative is act to the left, the output looks like { { {a},{b},...}, {T^nT^m...}, {A1, A2, ...}}, ...}; for left=False,  the output looks like { {{A1, A2, ...},, {T^nT^m...}, {a},{b},...}}, ...}; for the case no gluons, the {T^nT^m...}, {A1, A2, ...} become {} *)
derexpand[der_List,left_,adjoint_,sunindx_]:=Block[{tmp,tmp2,list,gluon,lindx,col,sunf,ig},
tmp=Replace[der,{{lor_}:>lindx[lor],lor:Except[{_}]:>{lindx[lor],gluon[lor,Unique["$Col"],{}]}},{1}];

If[left==True,
tmp=tmp//Reverse
];

tmp={list@@tmp};

(* expand covariante derivative *)
tmp=tmp//.{list[aa___,{lor1_,gl1_gluon},bb___]:>{list[aa,lor1,bb],list[aa,gl1,bb]}};

(* commute derivative and A  *)
tmp=tmp//.list[aa___,lor1_lindx,gg_gluon,bb___]:>{list[aa,gg,lor1,bb],list[aa,gluon[gg[[1]],gg[[2]],Append[gg[[3]],{lor1[[1]]}]],bb]};
tmp=Flatten[tmp]/.{gluon[lor1_,col1_,{}]:>gluon[lor1,col1]}/.list->List;


If[left==True,
tmp=Reverse[#]&/@tmp;
tmp=Replace[tmp,{ds___lindx,gs___gluon}:>{{#}&/@{ds},{gs}},1]
,
tmp=Replace[tmp,{gs___gluon,ds___lindx}:>{{gs},{#}&/@{ds}},1]
];

tmp=tmp/.lindx->Identity;

(* A_u A_v   -> A^n_u A^m_v T^n.T^m *)
Which[!adjoint&&!left,
col=Unique["$col"];
tmp=Replace[#,{glu_List,lor_List}/;Length[glu]>0:>{glu,{(-I gStrong)^Length[glu]SUNTF[glu/.gluon[_,suindx_,___]:>suindx,sunindx,col],col},lor},{0}]&/@tmp
,
!adjoint&&left,
col=Unique["$col"];
tmp=Replace[#,{lor_List,glu_List}/;Length[glu]>0:>{lor,{col,(I gStrong)^Length[glu]SUNTF[glu/.gluon[_,suindx_,___]:>suindx,col,sunindx]},glu},{0}]&/@tmp
,
adjoint&&!left,
tmp=Replace[#,{glu_List,lor_List}/;Length[glu]>0:>{glu,(sunf[{#[[2]]}]&/@glu)//.{
aa___,sunf[co1___,co2_List],sunf[co3_List,co4___],bb___}:>(col=Unique["$Col"];{aa,sunf[co1,co2,col],sunf[col,co3,co4],bb}),lor},
{0}]&/@tmp;
tmp=tmp/.sunf[co1_List,co2_]:>sunf[sunindx,co1,co2];
tmp=tmp/.{suns___,sunf[co1_,co2_List]}:>(col=Unique["$Col"];{Times@@{suns, sunf[co1,co2,col]},col});
tmp=tmp/.{sunf[co_List]}:>(col=Unique["$Col"];{sunf[sunindx,co,col],col});
tmp=tmp/.sus_sunf:>gStrong SUNF@@Flatten[List@@sus]
,
adjoint&&left,
tmp=Replace[#,{lor_List,glu_List}/;Length[glu]>0:>{lor,(sunf[{#[[2]]}]&/@glu)//.{
aa___,sunf[co1___,co2_List],sunf[co3_List,co4___],bb___}:>(col=Unique["$Col"];{aa,sunf[co1,co2,col],sunf[col,co3,co4],bb}),glu},
{0}]&/@tmp;
tmp=tmp/.sunf[co1_,co2_List]:>sunf[co1,co2,sunindx];
tmp=tmp/.{sunf[co1_List,co2_],suns___}:>(col=Unique["$Col"];{col,Times@@{sunf[co1,co2,col],suns}});
tmp=tmp/.{sunf[co_List]}:>(col=Unique["$Col"];{col,sunf[col,co,sunindx]});
tmp=tmp/.sus_sunf:>-gStrong SUNF@@Flatten[List@@sus]
];

tmp/.gluon->Gluon
]


(* generate the vertices in feynman rules *)
BasicVertex[cc_]:=Block[{tmp,z,tmploru,loru,tmplorv,lorv,tmpi,di,tmpj,dj,tmpa,suna,tmpb,sunb,tmpn,sunn,tmpm,sunm,tmpfa,tmpfb, tmpfc,tmpfd,fa,fb, fc,fd,tmpf,factor,flavor},
(* unique dummy indices *)
tmploru=Unique[loru];
tmplorv=Unique[lorv];
tmpi=Unique[di];
tmpj=Unique[dj];
tmpn=Unique[sunn];
tmpm=Unique[sunm];
tmpa=Unique[suna];
tmpb=Unique[sunb];
tmpfa=Unique[fa];
tmpfb=Unique[fb];
tmpfc=Unique[fc];
tmpfd=Unique[fd];
flavor=Unique["psi$"];

(*---------------------------*)
Which[cc==="v_qgq",
tmp=vertex[Unique[z],DCHN[DiracGamma[LorentzIndex[tmploru,D],D],tmpi,tmpj]SUNTF[tmpn,tmpa,tmpb],{AntiFermion[flavor,tmpi,tmpa],Fermion[flavor,tmpj,tmpb],Gluon[tmploru,tmpn]}];
factor=I gStrong
,
cc==="v_ggg",
tmp=vertex[Unique[z],SUNF[tmpfa,tmpfb, tmpfc],{Gluon[tmplorv,tmpfa,{{tmploru}}],Gluon[tmploru,tmpfb],Gluon[tmplorv,tmpfc]}];
factor=-I gStrong
,
cc==="v_gggg",
tmp=vertex[Unique[z],SUNF[tmpn,tmpfa, tmpfb]SUNF[tmpn,tmpfc, tmpfd],{Gluon[tmploru,tmpfa],Gluon[tmplorv,tmpfb],Gluon[tmploru,tmpfc],Gluon[tmplorv,tmpfd]}];(* ommit the 1/4; the duplicated contractions should be avoid in diagramGenerator *)
factor=-I gStrong^2
,
cc==="v_cgc",
tmp=vertex[Unique[z],SUNF[tmpfa,tmpfb, tmpfc],{Gluon[tmploru,tmpfa],AntiGhost[tmpfb,{{tmploru}}],Ghost[tmpfc]}];
factor=gStrong
];

{factor,tmp}
]



(* generate the vertices needed for the digram with specified order *)
VerticesPrepare[ngv_,ngluons_,order_,condlist_List,cond_]:=Block[{tmp1,tmp2,oddv,ngo1,ngo2,nglu,nq1,nq2,nv,k,l,i,vlist={{}}},
(* ngv = the number of g in operator list, condlist = the condensates and/or nonlocalcondensate involved in the diagram, cond = the required condensate *)
tmp1=(cond/.Condensate[n:Except[_List]]:>Condensate[{n,0,0}])/.{n_Integer/;n>0,con_Condensate}:>Table[con,n]//Flatten;
tmp2=(condlist/.Condensate[n:Except[_List]]:>Condensate[{n,0,0}])/.{n_Integer/;n>0,con_Condensate}:>Table[con,n]//Flatten;

ngo1=Total[Cases[tmp1,Condensate[{_,gg_,_}]:>gg,Infinity]];
ngo2=Total[Cases[tmp2,Condensate[{_,gg_,_}]:>gg,Infinity]]+Total[Cases[tmp2,GluonStrength[]->1,Infinity]];
nq1=Total[Cases[tmp1,Condensate[{nq_,_,_}]:>nq,Infinity]];
nq2=Total[Cases[tmp2,Condensate[{nq_,_,_}]:>nq,Infinity]]+Total[Cases[tmp2,Fermion[]|AntiFermion[]->1,Infinity]];
(*Print[order," ?????",ngo1,"---",ngo2,"---",nq1,"---",nq2];*)
(* nglu= Total[Cases[tmp2,Condensate[{_,ngl_,_}]:>ngl,Infinity]]+Total[Cases[tmp2,GluonStrength[]->1,Infinity]]; *)(* the number of gluon in condensates *)
nglu=ngluons-ngo2;
(*Print[ngv,">>>>",ngluons,"----",nglu,"-----",condlist,"-----",cond];*)


(* the g accompanied with G^nv is accounted into the Condensate, not accounted in the order of the diagram *)
nv=order+ngo1-ngv;
(*oddv= order-ngv-ngo1;*)
(* odd number of perturbative vertices is allowed, e.g. q^bar Gq and q^bar q give <q^bar q> Condensate at one loop level *)
(* for the case that the condensate is comes from the expansion of nonlocal-condensate *)
If[ngo2<ngo1||nq2<nq1,
(* extra g generated from DG -> g \bar{\psi}\psi, DD -> g G;
 e.g. GGG-> \bar{\psi}\psi GG generate one g, the number of G decrease, but the change of g should is not relevant in this case *)
nv=nv- Max[0,(nq1-nq2)/2]-Max[0,(ngo1-ngo2)](*;
oddv= oddv  -  Max[0,(nq1-nq2)/2]-Max[0,(ngo1-ngo2)]*)
];
(*Print[oddv,">=====>>>",nv];*)
If[nv<0,
(* lowest order connected diagram has order > than specified order *)
(* Message[QCDSumDiagram::owan]*)(* mute this warning since it's too common, e.g. expand G^nv, D^u ... will generate higher g^n in the operator *)
(* label this case as -1 *) 
vlist={{-1}}
,

For[k=0,k<Floor[nv/2]+1,k++,(* corresponding to g^2 AAAA vertex *)
For[l=0,l<nv-2k+1,l++,
For[i=0,i<nv-2k-l+1,i++,
vlist=Append[vlist,{Table[BasicVertex["v_qgq"],i],Table[BasicVertex["v_ggg"],nv-2k-l-i],Table[BasicVertex["v_gggg"],k],Table[BasicVertex["v_cgc"],l]}]
]
]
];
(* Print[OddQ[nglu+Total[Cases[#,_Gluon->1,Infinity]]]&/@vlist,"------",vlist]; *)


If[!FreeQ[vlist,_vertex],(* delete the cases involving odd number of perturbative gluons for nonempty vlist *)
vlist = DeleteCases[vlist,vs_/;OddQ[nglu+Total[Cases[vs,_Gluon->1,Infinity]]]];
If[FreeQ[vlist,_vertex],
vlist={{-1}}
]
];
(* cannot Contract to a diagram if odd number of perturbarive gluons involved for all cases *)
If[FreeQ[vlist,_vertex]&&OddQ[nglu],
vlist={{-1}}
];

vlist = DeleteCases[DeleteCases[#,{}]&/@vlist,{}];
If[vlist==={},vlist={{}}]
];

vlist
]


(* turn basic condensates <qq>, <GG>, <qgq>, <GGG> into conVertex[fields], and nonlocalcondensate <q...G...> into nlocalondvertex[fields] *)
(* the input comes from the output of condInterpretor *)
CondtoVertex[con_Condensate|con_NonLocalCondensate]:=Block[{tmp},

(* Condensate to vertex *)
Which[MatchQ[con,Condensate[2|{2,0}|{2,0,0}]],
conVertex[AntiFermion[],Fermion[]]
,
MatchQ[con,Condensate[{0,2}|{0,2,0}]],
conVertex[GluonStrength[],GluonStrength[]]
,
MatchQ[con,Condensate[{0,3}|{0,3,0}]],
conVertex[GluonStrength[],GluonStrength[],GluonStrength[]]
,
MatchQ[con,Condensate[{2,1}|{2,1,0}]],
conVertex[AntiFermion[],GluonStrength[],Fermion[]]
,
MatchQ[con,_NonLocalCondensate],
nlocalconVertex@@con[[2]]
]

]


(* for the product of basic condensates *)
condInterpretor[conds_Times,nonlocal_:False]:=Block[{tmp,nulllist},
tmp=conds//Expand;
If[MatchQ[tmp,_Plus]||!FreeQ[tmp,Power[Condensate[{_,_,Except[0]}],nn_]],
Message[QCDSumDiagram::cerr];
Abort[]
];

tmp=List@@tmp;
nulllist=tmp/.Power[con_Condensate,nn_Integer]/;nn>0:>con/.Condensate[nn_List]:>nn/.{{2,1,0}->0,{2,0,0}->0,{2,1,1}->1,{2,0,1}->1,{0,2,0}->0,{0,3,0}->0,{0,2,1}->1,{0,3,1}->1};
nulllist=DeleteCases[nulllist,0];
If[!MatchQ[nulllist,{1}|{}],
Message[QCDSumDiagram::cerr];
Abort[]
];

tmp=Replace[tmp,{con_Condensate:>{1,con},Power[con_Condensate,ni_]:>{ni,con}},1];
If[!FreeQ[tmp,Condensate[{_,_,1}]],
{{0,tmp},{1,tmp/.Condensate[{nq_,ng_,1}]:>Condensate[{nq,ng,0}]}}
,
{{0,tmp}}
]
]

condInterpretor[ncond_NonLocalCondensate,nonlocal_:False]:={{0,ncond}}

(* turn Condensates into condensate "vertices", e.g., "qgq" -> Cond[integer label, DCHN[DiracSigma[GAD[\[Mu],\[Nu]]],i,j] SUNTF[n,a,b], {AntiFermion[f,i,a], Fermion[f,j,b], GluonStrength[{\[Mu],\[Nu],n}] } ] *)
condInterpretor[conds_Condensate,nonlocal_:False]:=Block[{tmp,condlist={},condensate,nq=0,ngg,ng=0,nm=0,null,i,j,k,m,mcondlist={},tmpmcondlist={},nonlocallist={}},

If[FreeQ[conds,_Integer],
Message[QCDSumDiagram::cerr];
Abort[]
];

conds/.Condensate->condensate/.{condensate[n:Except[_List]]:>(nq=n;null),condensate[{n1_,n2_}]:>(nq=n1;ng=n2;null),condensate[{n1_,n2_,n3_}]:>(nq=n1;ng=n2;nm=n3;null)};
If[OddQ[nq],Message[Condensate::err];Abort[]];

(* factorization *)
For[i=0,i<=nq/2,i++,(* how many <\bar{q}G q > allowed *)
For[j=0,j<Ceiling[(ng-i)/2]+1,j++,(* how many <GG> allowed *)
If[IntegerQ[(ng-i-2j)/3]&&(ng-i-2j)/3>=0,
condlist=Append[condlist,{nq/2-i, i, j, (ng-i-2j)/3}](* {n1_qq, n2_qgq, n3_gg, n4_ggg} *)
]
]
];

condlist=Replace[DeleteCases[condlist,{0,0,0,0}|{}],{n1_,n2_,n3_,n4_}:>{nm(* mass^n *),{n1,Condensate[{2,0,0}]},{n2,Condensate[{2,1,0}]},{n3,Condensate[{0,2,0}]},{n4,Condensate[{0,3,0}]}},{1}];
condlist=DeleteCases[#,{0,_}]&/@condlist;

(*-------------------------------------------------*)
(* expansion of nonlocal condensate with lower dimension *)
If[nonlocal,

(* DG -> g \bar{q}q, the number of \bar{q}q is always increasing by expanding fields not at 0 *)
For[i=1,i<nq/2,i++,(* for the cases 0 < number of \bar{q}q < nq/2 *)
For[j=0,j<=Floor[(nq-2i) 3/4+ng],j++,(* the number of gluons *) 
nonlocallist=Append[nonlocallist,{nq/2-i,j}]
]
];
(* recall the equation of motion nabla_slash \[Psi] = -im\[Psi], so the number of alowed gluons is determined by nq3/2+2ng, e.g., m^3<\bar{q}Gq> cannot comes form the expansion of <\bar{q}GGq> *)
(* for i=0 and i=nq/2 *)
If[nq==0,(* no quarks involved at first *)
(* the allowed G less than ng *)
nonlocallist=Join[nonlocallist,Table[{0,j},{j,2,ng-1}]];
,
(* i=nq/2, the allowed G must less than ng *)
nonlocallist=Join[nonlocallist,Table[{nq/2,j},{j,0,ng-1}]];
(* i=0, the number of G must >=2, since no <G> *)
nonlocallist=Join[nonlocallist,Table[{0,j},{j,2,Floor[nq 3/4+ng]}]];
];


(* all possible m^(nm-n) <m^n qG...>, the m^(nm-n) comes from the mass in the propagator; <m^n qG...> means the m^n is generated form the expansion of \[Psi](x) and G(x) in nonlocal-condensate *)
nonlocallist=Flatten[Table[{mm,#},{mm,0,nm}]&/@nonlocallist,1];
(* to nonlocalCondensate symbols *)
nonlocallist={#[[1]],NonLocalCondensate[nm-#[[1]],Join[Table[{AntiFermion[],Fermion[]},#[[2,1]]]//Flatten,Table[GluonStrength[],#[[2,2]]]], {conds}]}&/@nonlocallist;
(* the element of condlist = {number of mass, number of \bar{q}q, number of G } *)
(* the element of nonlocallist = {number of mass comes from propagator, nonlocal-condensate with mass from the expansion of fields, required(the input) condensate} *)
];

(*----------------------------------------------*)
(* for the case involve quark mass and quark fields, the quark mass can also come from expansion of fermion fields, except the masses, the condensate structure and the diagram are same as the cases for condlist; e.g. m <\bar{q}Gq> for the case that m comes from propagators and the case that comes from expansion of \bar{q} and q *)
(* the difference between the condensates in nonlocallist and mcondlist is that in nonlocallist, the Condensate sturcture are different with the condensates in condlist, e.g. for m <\bar{q}Gq>, it can comes from the expansion of <\bar{q}(x)q(0)> *)
If[nm>0,
mcondlist={condlist}/.Condensate[cc_List]:>nlocalcon[0,cc];
For[i=1,i<nm+1,i++,
tmpmcondlist={};
For[j=1,j<Length[mcondlist[[-1]]]+1,j++,
For[k=2,k<Length[mcondlist[[-1,j]]]+1,k++,

If[MatchQ[mcondlist[[-1,j,k]],{1,nlocalcon[nmc_,{nnq_/;nnq>0,__}]}],

tmpmcondlist=Append[tmpmcondlist,Join[{mcondlist[[-1,j,1]]-1},mcondlist[[-1,j,2;;k-1]],{{1,nlocalcon[mcondlist[[-1,j,k,2,1]]+1,mcondlist[[-1,j,k,2,2]]]}},mcondlist[[-1,j,k+1;;]]]]
,
If[MatchQ[mcondlist[[-1,j,k]],{nnc_/;nnc>1,nlocalcon[nmc_,{nnq_/;nnq>0,__}]}],

(*Print[mcondlist[[-1,j,k]]];*)
tmpmcondlist=Append[tmpmcondlist,Join[{mcondlist[[-1,j,1]]-1},mcondlist[[-1,j,2;;k-1]],{{mcondlist[[-1,j,k,1]]-1,nlocalcon[mcondlist[[-1,j,k,2,1]],mcondlist[[-1,j,k,2,2]]]},{1,nlocalcon[mcondlist[[-1,j,k,2,1]]+1,mcondlist[[-1,j,k,2,2]]]}},mcondlist[[-1,j,k+1;;]]]]
]
]
]
];

mcondlist=Append[mcondlist,tmpmcondlist]
];
mcondlist=Flatten[mcondlist[[2;;]],1];
(* combine same nonlocal condensates *)
mcondlist=(#//.{aa___,{n1_,nc_nlocalcon},bb___,{n2_,nc_nlocalcon},cc___}:>{aa,{n1+n2,nc},bb,cc})&/@mcondlist
];

mcondlist=mcondlist/.{nlocalcon[0,con_]:>Condensate[con],nlocalcon[mmm_/;mmm>0,con_]:>NonLocalCondensate[mmm,Join[Table[{AntiFermion[],Fermion[]},con[[1]]/2]//Flatten,Table[GluonStrength[],con[[2]]]]]};
(* duplicates could be exist, expand n-th term and m-th term in different order yeilds same structure *)
(* to do: avoid duplicates in above code *)
mcondlist=DeleteDuplicates[mcondlist];


(*----------------------------------------------*)
If[nonlocal,
Flatten[DeleteCases[{condlist,mcondlist,nonlocallist},{}],1]
,
Flatten[DeleteCases[{condlist,mcondlist},{}],1]
]
]










(* generate diagrams for each Condensate *)
diagramGenerator[ovc_List,in_,out_,holdf_,tadpole_:"all"]:=Block[{tmp,gslist,suna,sunb,gtlist,fieldlist,factor,olist,tmpv,tmpc,operator,vertex,v3,v4,cv2,cv3,ncvqb,ncvq,ncvg,field,tmpfields,match,dia,fields,dialist,ppgator,list},

tmp={ovc[[1,1]],ovc[[2]]};
(* for the case that other A^u involved in the oeprator, restore the GluonStrength[{u,v,a},_,1] as g f^abc A^au A^bv *)
(* to do : combine all A^u in operator as a gluon vertex, to avoid generate diagrams merely differ by the indices of gluons *)
tmp=tmp/.Operator[xx_,vv_,fs_List]/;!FreeQ[fs,_Gluon]&&!FreeQ[fs,GluonStrength[_,_,1]]:>(gslist=Cases[fs,GluonStrength[{lor1_,lor2_,sunn_},{},1]:>(suna=Unique["$Col"];sunb=Unique["$Col"];{gStrong SUNF[sunn,suna,sunb],{Gluon[lor1,suna],Gluon[lor2,sunb]}}),Infinity];
gslist=Transpose[gslist];Operator[xx,vv(Times@@gslist[[1]]),Join[DeleteCases[fs,GluonStrength[_,_,1]],Flatten[gslist[[2]]]] ]);


factor=Cases[tmp,{gn_,_vertex}:>gn,Infinity];(* g factor in the vertex *)
(* for gluon vertices f^abc d^uA^av A^b_u A^c_v and f^nab f^ncd A^au A^bv A^c_u A^d_v, label them as v3 and v4,
the possible contractions for gluoms and vertices are A-v3 A-v3 A-v3,  A-v3 A-v3 B-v3,  A-v4 A-v4 A-v4 A-v4,  A-v4 A-v4 A-v4 B-v4,  A-v4 A-v4 B-v4 B-v4
where B is background gluon fields,
 insert the explict vertices only after all fields in the vertices are contracted *)
(* extract the gamma and su(N) f and t matrices, expect the f^abc in gluon vertices *)
tmp=tmp/.{vertex[zz_,sunt_,{g1_Gluon,g2_Gluon,g3_Gluon}]:>v3[zz,3],vertex[zz_,sunt_,{g1_Gluon,g2_Gluon,g3_Gluon,g4_Gluon}]:>v4[zz,4]};(* the m in vn[zz,m] labels how many gluons remain in the vertex *)
gtlist=Join[Cases[tmp,Operator[_,vv_,_List]:>vv,Infinity], Cases[tmp,vertex[_,sunt_,fi_List]:>sunt,Infinity]];
tmp=tmp/.{gs_, vv_vertex|vv_v3|vv_v4}:>vv;

(* label the operators by integers, since operators at same position occur after applying OperatorFactorizor; add integers for the vertices and condensates to unify the form *)
fieldlist={Table[{i,tmp[[1,i]]},{i,1,Length[tmp[[1]]]}],
Table[{j,#[[j]]},{j,1,Length[#]}]&/@tmp[[2]],
If[MatchQ[ovc[[1,2,2,3]],_nlocalconVertex],

{
DeleteCases[{
{field[0,1,"nc",ncvqb[Length[Cases[ovc[[1,2,2,3]],_AntiFermion,Infinity]]]]},
{field[0,2,"nc", ncvq[Length[Cases[ovc[[1,2,2,3]],_Fermion,Infinity]]]]},
 {field[0,1,"nc",ncvg[Length[Cases[ovc[[1,2,2,3]],_GluonStrength,Infinity]] ]]} 
},
 {field[_,_,_,nvv_/;nvv[[1]]==0]}]
} (* unify the structure; avoid duplicate and introducing symmetry factor; ncvqb = nonlocalcondensate-vertex_\bar{quark}, ncvqb = nonlocalcondensate-vertex_quark, ncvg = nonlocalcondensate-vertex_gluon, *)
,
Table[{j,#[[j]]},{j,1,Length[#]}]&/@ovc[[1,2,2,3]] 
]
};

fieldlist= {Flatten[fieldlist[[1]]/.{nn_Integer, Operator[xx_,vv_,fs_List]}:>(field[xx,nn,"o",#]&/@fs)/.GluonStrength[{lor1_,lor2_,cola_},{},1]:>ovg[{lor1,lor2,cola},2] ]
,
fieldlist[[2]]/.{{nn_Integer, vertex[zz_,suntf_, fieldss_List]}/;!FreeQ[fieldss,_Fermion]:>(field[zz,nn,"vq",#]&/@fieldss),{nn_Integer, vertex[zz_,suntf_, fieldss_List]}/;!FreeQ[fieldss,_Ghost]:>(field[zz,nn,"vc",#]&/@fieldss),{nn_Integer, vvv_v3}:>{field[vvv[[1]],nn,"v3",Drop[vvv,1]]},{nn_Integer, vvv_v4}:>{field[vvv[[1]],nn,"v4",Drop[vvv,1]]}}
,
fieldlist[[3]]/.{{nn_,conVertex[_AntiFermion,_Fermion]}:>{field[0,nn,"c2q",AntiFermion[]],field[0,nn,"c2q",Fermion[]]},{nn_,conVertex[_GluonStrength,_GluonStrength]}:>{field[0,nn,"c2g",cv2[2]]},{nn_,conVertex[_GluonStrength,_GluonStrength,_GluonStrength]}:>{field[0,nn,"c3g",cv3[3]]},{nn_,conVertex[_AntiFermion,_GluonStrength,_Fermion]}:>{field[0,nn,"cqgq",AntiFermion[]],field[0,nn,"cqgq",GluonStrength[]],field[0,nn,"cqgq",Fermion[]]}}(* the last case corresponding to nonlocalcondensate *)
};

(* the structure of fieldlist
fieldlist[[1]]: {field[xx_,nn_,"o",field_1], ... }
fieldlist[[2]]: { { {fields in the vertex}, second list of fields for same type of vertex, ... }, another type of vertices, ...}
fieldlist[[3]]: { { {fields in the condensate}, second list of fields for same type of condensate, ... }, another type of condensates, ...} 
*)

(*-------------------------------------------------------------*)
(* Wick contraction *)
(* for each field in fieldlist, contract it with the later fields, anticommute the fermions *)
(* to avoid introducing symmetry factors, for same type of vertices, contract the later vetex only when the former vertices have been contracted *)
(* for first term in fieldlist[[1]], find all contractable fields in {fiellist[[1, 2;;]], fieldlist[[2]], fieldlist[[3]]} *)
(* the positionfs for contractable terms in fieldlist[[2]] and fieldlist[[2]] looks like: {{indx in level_1, indx in level_2, indx in level_3}, {}, ... } *)
(* for each term in level_1, first) pick the first list of fields in level_2, second) drop the contracted field, append remain fields into fieldlist[[1]], for v3[n], v4[n], cv2[n], cv3[n], this means n -> n-1, discard n==0 *)
(* when contracting two fields, only the case that the positions of them are different are allowed *)

(*dialist={list[dia[1,{}],fieldlist]}//.list[di_dia,fid_List]:>fieldcontract[di,fid]*)
dialist=FixedPoint[Flatten[#/.list[di_dia,fid_List]:>fieldcontract[di,fid,holdf,tadpole]]&,{list[dia[1,{}],fieldlist]}];
(* delete the diagrams that are not single connected *)

dialist=DeleteCases[dialist,ddia_/;!validQ[ddia]];

dialist=duplicatremove[dialist];

(* ppgator[f1_,f2_, ...] -> ppgator[{f1,f2}, +-1, ...] after propagatorCombine, the contraction between field and Condensate gives condline, dia becomes dia[+-1, {propagators list}, {condline list}] *)
dialist=propagatorCombine[#]&/@dialist;
dialist=DeleteCases[dialist,aa_/;!FreeQ[aa,propagator[_,{xx_,xx_},{vv_,vv_},{nn_,nn_}]|propagator[_,{xx_,xx_},{"o","o"},{_,_}]|"tadpole"]];
dialist=dia[#[[1]],#[[2]],#[[3]]]&/@dialist;
dialist=DeleteCases[dialist,ddia_/;TadPoleQ[ddia]]
(*fieldlist*)

]


fieldcontract[di_,fieldlist_List,holdf_,tadpole_:"all"]:=Block[{field1,tmp1,tmp2,tmp3,match1,match2,match3,sign =1,fidlist1,fidlist2,fidlist3,tmpfields={},fhead,fhead2,fsign=1,fidsremain,gluremove,list},
If[fieldlist[[1]]==={}&&!FreeQ[fieldlist,_field],
{}(* all fields in fieldlist[[1]] are contracted, remaining fields in fieldlist[[2]] or fieldlist[[3]] can only conntract to a disconnected diagram *)
,
field1=fieldlist[[1,1]];

(* remove first field in fieldlist[[1]] *)
If[MatchQ[fieldlist[[1,1,4]],v3[n_/;n>1]|v4[n_/;n>1]|cv2[n_/;n>1]|cv3[n_/;n>1]|ncvqb[n_/;n>1]|ncvq[n_/;n>1]|ncvg[n_/;n>1]|ovg[_,n_/;n>1]],
tmp1=Join[{fieldlist[[1,1]]/.{v3[n_]:>v3[n-1],v4[n_]:>v4[n-1],cv2[n_]:>cv2[n-1],cv3[n_]:>cv3[n-1],ncvqb[n_]:>ncvqb[n-1],ncvq[n_]:>ncvq[n-1],ncvg[n_]:>ncvg[n-1],ovg[inds_,n_]:>ovg[inds,n-1]}},fieldlist[[1,2;;]]]
,
tmp1=fieldlist[[1,2;;]]
];

tmp2=fieldlist[[2]];
tmp3=fieldlist[[3]];


Which[
(* fermions and ghosts lines; ncvqb, ncvq means the \bar{\psi} \psi in nonlocalcondensate *)
MatchQ[field1[[4]],_Fermion|_AntiFermion|_Ghost|_AntiGhost|_ncvqb|_ncvq],
Which[MatchQ[field1[[4]],_Fermion|_ncvq],
fhead=AntiFermion;
fhead2=ncvqb
,
MatchQ[field1[[4]],_AntiFermion|_ncvqb],
fhead=Fermion;
fhead2=ncvq;
fsign=-1 (* extra -1 by anti commuting \bar{\psi}\psi -> - \psi\bar{\psi} *)
,
MatchQ[field1[[4]],_Ghost],
fhead=AntiGhost
,
MatchQ[field1[[4]],_AntiGhost],
fhead=Ghost;
fsign=-1
];

(* ------------------------------------------------------------ *)
If[FreeQ[{tmp1,tmp2,tmp3},_fhead|_fhead2],

tmpfields={}(* no contractable fields exist *)
,
(* the number of (anti)fermions or (anti)ghosts in fieldlist[[2]] and fieldlist[[3]] is always even *)

(* take the fields from fieldlist[[1]] *)
match1=Position[tmp1,field[_,_,_,_fhead|_fhead2]]//Flatten;

If[Length[match1]>0,
fidlist1=Table[{(* the -1 sign generated by anti commuting (anti)fermions and (anti)ghost *)
(-1)^(Length[Cases[tmp1[[;;ii-1]],_Fermion|_AntiFermion|_Ghost|_AntiGhost,Infinity]]+Total[Cases[tmp1[[;;ii-1]],ncvqb[nn_]:>nn,Infinity]]+Total[Cases[tmp1[[;;ii-1]],ncvq[nn_]:>nn,Infinity]])
,
If[fsign==1,
ppgator[field1[[4]],tmp1[[ii,4]],{field1[[1]],tmp1[[ii,1]]},{field1[[3]],tmp1[[ii,3]]},{field1[[2]],tmp1[[ii,2]]}]
,
ppgator[tmp1[[ii,4]],field1[[4]],{tmp1[[ii,1]],field1[[1]]},{tmp1[[ii,3]],field1[[3]]},{tmp1[[ii,2]],field1[[2]]}]
]
,
(* remain fields in tmp1 *)
Join[
tmp1[[;;ii-1]]
,
If[MatchQ[tmp1[[ii,4]],_Fermion|_AntiFermion|_Ghost|_AntiGhost],
{}
,(* remove the contracted (anti)fermions or (anti)ghost, for ncvqb[n], ncvq[n], this means n -> n-1, discard n==0 *)
{tmp1[[ii]]/.{ncvqb[n_]:>ncvqb[n-1],ncvq[n_]:>ncvq[n-1]}}/.{field[_,_,_,ncvqb[0]|ncvq[0]]}->{}
]
,
tmp1[[ii+1;;]] 
]
}
,{ii,match1}];
tmpfields=Join[tmpfields,list[dia[fsign di[[1]] #[[1]],Append[di[[2]],#[[2]]]],{#[[3]],tmp2,tmp3}]&/@fidlist1]
];


(* take the fields from fieldlist[[2]] *)
match2=Position[tmp2,field[_,_,_,_fhead]]//.{indxa___,{level1_,levels___},{level1_,___}..,indxb___}/;FreeQ[{indxa},{level1,___}]:>{indxa,{level1,levels},indxb};
(* keeping only the first vertices in the same type of vertices *)
If[Length[match2]>0,
fidlist2=Table[{(* the -1 sign generated by anti commuting (anti)fermions and (anti)ghost *)
(-1)^(Length[Cases[tmp1,_Fermion|_AntiFermion|_Ghost|_AntiGhost,Infinity]]+Total[Cases[tmp1,ncvqb[nn_]:>nn,Infinity]]+Total[Cases[tmp1,ncvq[nn_]:>nn,Infinity]])
,
If[fsign==1,
ppgator[field1[[4]],tmp2[[##,4]],{field1[[1]],tmp2[[##,1]]},{field1[[3]],tmp2[[##,3]]},{field1[[2]],tmp2[[##,2]]}]&@@inds
,
ppgator[tmp2[[##,4]],field1[[4]],{tmp2[[##,1]],field1[[1]]},{tmp2[[##,3]],field1[[3]]},{tmp2[[##,2]],field1[[2]]}]&@@inds
]
,
{
Join[tmp1, tmp2[[inds[[1]],inds[[2]],;;inds[[3]]-1 ]] , tmp2[[inds[[1]],inds[[2]],inds[[3]]+1;;-1 ]] ] (* append the remain fields in vertex into tmp1 *)
,
DeleteCases[Delete[tmp2,inds[[;;2]]],{}(* if all fields in this type are removed *)]
,
tmp3
}
}
,{inds,match2}
];
tmpfields=Join[tmpfields,list[dia[fsign di[[1]] #[[1]],Append[di[[2]],#[[2]]]],#[[3]]]&/@fidlist2]
];


(* take the fields from fieldlist[[3]] *)
match3=Position[tmp3,field[_,_,_,_fhead|_fhead2]]//.{indxa___,{level1_,levels___},{level1_,___}..,indxb___}/;FreeQ[{indxa},{level1,___}]:>{indxa,{level1,levels},indxb};
(* keeping only the first vertices in the same type of vertices *)
If[Length[match3]>0,
fidlist3=Table[{(* the -1 sign generated by anti commuting (anti)fermions and (anti)ghost *)
(-1)^(Length[Cases[tmp1,_Fermion|_AntiFermion|_Ghost|_AntiGhost,Infinity]]+Total[Cases[tmp1,ncvqb[nn_]:>nn,Infinity]]+Total[Cases[tmp1,ncvq[nn_]:>nn,Infinity]])
,
If[fsign==1,
ppgator[field1[[4]],tmp3[[##,4]],{field1[[1]],tmp3[[##,1]]},{field1[[3]],tmp3[[##,3]]},{field1[[2]],tmp3[[##,2]]}]&@@inds
,
ppgator[tmp3[[##,4]],field1[[4]],{tmp3[[##,1]],field1[[1]]},{tmp3[[##,3]],field1[[3]]},{tmp3[[##,2]],field1[[2]]}]&@@inds
]
,
{
(* append the remain fields into tmp1 *)
If[FreeQ[tmp3[[inds[[1]],inds[[2]] ]], _ncvqb|_ncvq],
Join[tmp1,  tmp3[[inds[[1]],inds[[2]],;;inds[[3]]-1 ]] ,  tmp3[[inds[[1]],inds[[2]],inds[[3]]+1;;-1 ]] ] 
,(* for ncvqb[n], ncvq[n], this means n -> n-1, discard n==0 *)
Join[tmp1,  tmp3[[inds[[1]],inds[[2]] ]]/.{ncvqb[n_]:>ncvqb[n-1],ncvq[n_]:>ncvq[n-1]}/.{field[_,_,_,ncvqb[0]|ncvq[0]]}->{} ]
]
,
tmp2
,
DeleteCases[Delete[tmp3,inds[[;;2]]],{}(* if all fields in this type are removed *)]
}

}
,{inds,match3}
];
tmpfields=Join[tmpfields,list[dia[fsign di[[1]] #[[1]],Append[di[[2]],#[[2]]]],#[[3]]]&/@fidlist3]
]
]
,
(*------------------------------------------------------------------------*)
(*------------------------------------------------------------------------*)
(* gluon lines *)

MatchQ[field1[[4]],_Gluon|_GluonStrength|_v3|_v4|_cv2|_cv3|_ncvg|_ovg],

If[FreeQ[{tmp1,tmp2,tmp3},_Gluon|_GluonStrength|_v3|_v4|_cv2|_cv3|_ncvg|_ovg],
tmpfields={}(* no contractable fields *)
,
(* take the fields from fieldlist[[1]] *)
match1=Position[tmp1,field[_,_,_,_Gluon|_GluonStrength|_v3|_v4|_cv2|_cv3|_ncvg|_ovg]]//Flatten;
If[Length[match1]>0,
fidlist1=Table[
{
ppgator[field1[[4]],tmp1[[ii,4]],{field1[[1]],tmp1[[ii,1]]},{field1[[3]],tmp1[[ii,3]]},{field1[[2]],tmp1[[ii,2]]}]
,
Join[
tmp1[[;;ii-1]]
,
If[MatchQ[tmp1[[ii,4]],_Gluon|_GluonStrength],
{}(* remove the contracted gluon(strength) *)
,(* for v3[n], v4[n], cv2[n], cv3[n], ncvg[n], this means n -> n-1; discard n==0 *)
{tmp1[[ii]]/.{v3[n_]:>v3[n-1],v4[n_]:>v4[n-1],cv2[n_]:>cv2[n-1],cv3[n_]:>cv3[n-1],ncvg[n_]:>ncvg[n-1],ovg[inds_,n_]:>ovg[inds,n-1]}}/.{field[_,_,_,v3[0]|v4[0]|cv2[0]|cv3[0]|ncvg[0]|ovg[_,0]]}->{}
]
,
tmp1[[ii+1;;]] 
]
}
,{ii,match1}];

tmpfields=Join[tmpfields,list[dia[di[[1]] ,Append[di[[2]],#[[1]]]],{#[[2]],tmp2,tmp3}]&/@fidlist1]
];



(* take the fields from fieldlist[[2]] *)
match2=Position[tmp2,field[_,_,_,_Gluon|_GluonStrength|_v3|_v4|_cv2|_cv3]]//.{indxa___,{level1_,levels___},{level1_,___}..,indxb___}/;FreeQ[{indxa},{level1,___}]:>{indxa,{level1,levels},indxb};
(* keeping only the first vertices in the same type of vertices *)
If[Length[match2]>0,
fidlist2=Table[
{
ppgator[field1[[4]],tmp2[[##,4]],{field1[[1]],tmp2[[##,1]]},{field1[[3]],tmp2[[##,3]]},{field1[[2]],tmp2[[##,2]]}]&@@inds
,
{
(* append the remain fields, (anti)fermions or (anti)ghosts in the vertex, into tmp1 *)
If[!FreeQ[tmp2[[inds[[1]],inds[[2]] ]], _Fermion|_Ghost],
Join[tmp1,  tmp2[[inds[[1]],inds[[2]],;;inds[[3]]-1 ]] ,  tmp2[[inds[[1]],inds[[2]],inds[[3]]+1;;-1 ]] ] 
,(* for v3[n], v4[n], cv2[n], cv3[n], this means n -> n-1, discard n==0 *)
Join[tmp1,  tmp2[[inds[[1]],inds[[2]] ]] /.{v3[n_]:>v3[n-1],v4[n_]:>v4[n-1],cv2[n_]:>cv2[n-1],cv3[n_]:>cv3[n-1],ncvg[n_]:>ncvg[n-1],ovg[inds_,n_]:>ovg[inds,n-1]}/.{field[_,_,_,v3[0]|v4[0]|cv2[0]|cv3[0]|ncvg[0]|ovg[_,0]]}->{} ]
]
,
DeleteCases[Delete[tmp2,inds[[;;2]]],{}(* if all fields in this type are removed *)]
,
tmp3
}
}
,{inds,match2}
];

tmpfields=Join[tmpfields,list[dia[di[[1]] ,Append[di[[2]],#[[1]]]],#[[2]]]&/@fidlist2]
];



(* take the fields from fieldlist[[3]] *)
match3=Position[tmp3,field[_,_,_,_Gluon|_GluonStrength|_v3|_v4|_cv2|_cv3|_ncvg]]//.{indxa___,{level1_,levels___},{level1_,___}..,indxb___}/;FreeQ[{indxa},{level1,___}]:>{indxa,{level1,levels},indxb};
(* keeping only the first vertices in the same type of vertices *)
If[Length[match3]>0,
fidlist3=Table[
{
ppgator[field1[[4]],tmp3[[##,4]],{field1[[1]],tmp3[[##,1]]},{field1[[3]],tmp3[[##,3]]},{field1[[2]],tmp3[[##,2]]}]&@@inds
,
{
(* append the remain fields into tmp1; ong qgq Condensate involving (anti)fermions *)
If[!FreeQ[tmp3[[inds[[1]],inds[[2]] ]], _Fermion],
Join[tmp1,  tmp3[[inds[[1]],inds[[2]],;;inds[[3]]-1 ]] ,  tmp3[[inds[[1]],inds[[2]],inds[[3]]+1;;-1 ]] ] 
,(* for v3[n], v4[n], cv2[n], cv3[n], ncvg[n], this means n -> n-1, discard n==0 *)
Join[tmp1,  tmp3[[inds[[1]],inds[[2]] ]] /.{v3[n_]:>v3[n-1],v4[n_]:>v4[n-1],cv2[n_]:>cv2[n-1],cv3[n_]:>cv3[n-1],ncvg[n_]:>ncvg[n-1],ovg[inds_,n_]:>ovg[inds,n-1]}/.{field[_,_,_,v3[0]|v4[0]|cv2[0]|cv3[0]|ncvg[0]|ovg[_,0]]}->{} ]
]
,
tmp2
,
DeleteCases[Delete[tmp3,inds[[;;2]]],{}(* if all fields in this type are removed *)]
}
}
,{inds,match3}
];

tmpfields=Join[tmpfields,list[dia[di[[1]],Append[di[[2]],#[[1]]]],#[[2]]]&/@fidlist3]
]
]

];
(*tmppppp=Append[tmppppp,Length[tmpfields]];*)
(* keeping only dia[...] if all fields have been contracted *)
tmpfields=Flatten[tmpfields/.list[dias_dia,fds_/;FreeQ[fds,_field]]:>dias];

(* delete samplest tadpole: one propagator loop *)
tmpfields=DeleteCases[tmpfields,ddia_/;!FreeQ[ddia,ppgator[_,_,_,{"c2q","c2q"}|{"c2g","c2g"}|{"c3g","c3g"}|{"cqgq","cqgq"}|{"nc","nc"},_]|ppgator[_,_,{zz_,zz_},{"v","v"},_]]];
If[tadpole==="all",
tmpfields=DeleteCases[tmpfields,ddia_/;!FreeQ[ddia,ppgator[_,_,{zz_,zz_},{"o","o"},_]|ppgator[_,_,{zz_,zz_},{"vq","vq"}|{"v3","v3"}|{"vc","vc"},_]]]
];
(* the connection between Condensate and Condensate is not alowed, the composited Condensate is accounted in the nonlocalcondensate, and the local condensates are always factorized *)
tmpfields=DeleteCases[tmpfields,ddia_/;!FreeQ[ddia,ppgator[_,_,{zz_,zz_},{v1_,v2_},_]/;FreeQ[{v1,v2},"o"|"vq"|"vc"|"v3"|"v4"]]];

If[holdf===False,(* delete the simplest case that flavors not match *)
tmpfields=DeleteCases[tmpfields,ddia_/;!FreeQ[ddia,ppgator[Fermion[fla_,___],AntiFermion[flb_,___],{zz_,yy_},{"o","o"},_]/;fla=!=flb]]
];

tmpfields
]
]


validQ[di_dia]:=momentumconnectedQ[di]&&connectedQ[di]&&momentumtransferQ[di]&&!condtadpoleQ[di]&&condnonvainshQ[di]&&FreeQ[di,ppgator[_,_,{_,_},{v1_,v2_},{_,_}]/;!FreeQ[{v1,v2},"vc"]&&!FreeQ[{v1,v2},"c2g"|"c3g"|"cqgq"|"nc"]]
(* the ghost vertex must connecting with perturbative gluon *)

(* whether a diagram is single connected by propagator, exclude the case that two loops connected by condensate lines *) 
momentumconnectedQ[di_dia]:=Block[{dias,tmp,tmp2,tmpdia,vs,edges,pindx,remains},

pindx=Position[di[[2]],ppgator[_,_,{_,_},{v1_,v2_},{_,_}]/;FreeQ[{v1,v2},"c2q"|"c2g"|"c3g"|"cqgq"|"nc"]];
If[Length[pindx]==0,
False(* if no propagator *)
,
(* discard the operator labels, all operators at same position form a single vertex *)
dias=di[[2]]/.ppgator[f1_,f2_,{xx_,yy_},{"o","o"},{n1_,n2_}]:>ppgator[f1,f2,{xx,yy},{"o","o"},{xx,yy}]/.{ppgator[f1_,f2_,{xx_,yy_},{"o",v2_},{n1_,n2_}]:>ppgator[f1,f2,{xx,yy},{"o",v2},{xx,n2}],ppgator[f1_,f2_,{xx_,yy_},{v1_,"o"},{n1_,n2_}]:>ppgator[f1,f2,{xx,yy},{v1,"o"},{n1,yy}]};

tmp=dias[[pindx[[1,1]]]];(* the first propagator *)
tmp={vs[tmp[[3,1]],tmp[[4,1]],tmp[[5,1]]], vs[tmp[[3,2]],tmp[[4,2]],tmp[[5,2]]]};

remains=Delete[dias,pindx[[1,1]]];
edges=FixedPoint[(tmp2={};
tmpdia=#//.{(* remove connected propagators, update the vertices list; remove the vertices that have already verified *)
ppgator[_,_,{x1_,x2_},{v1_,v2_},{n1_,n2_}]/;!FreeQ[tmp,vs[x1,v1,n1]]&&MatchQ[v2,"o"|"vq"|"vc"|"v3"|"v4"]:>(tmp2=Append[tmp2,vs[x2,v2,n2]];{})
,
ppgator[_,_,{x1_,x2_},{v1_,v2_},{n1_,n2_}]/;!FreeQ[tmp,vs[x2,v2,n2]]&&MatchQ[v1,"o"|"vq"|"vc"|"v3"|"v4"]:>(tmp2=Append[tmp2,vs[x1,v1,n1]];{})
 };tmp=tmp2;tmpdia
)&,remains(*di[[2,2;;]]*)];
(*tmp=DeleteCases[tmp2,vs[_,"c2q"|"c2g"|"c3g"|"cqgq"|"nc",_]](* remove connected Condensate vertices *)*)
(* if there still ramain propagators *)
If[FreeQ[edges,ppgator[_,_,{_,_},{v1_,v2_},{_,_}]/;FreeQ[{v1,v2},"c2q"|"c2g"|"c3g"|"cqgq"|"nc"]],
True
,
False
]
]
]


(* whether a diagram is single connected *)
connectedQ[di_dia]:=Block[{dias,tmp,tmp2,tmpdia,vs,edges,pindx,remains},

pindx=Position[di[[2]],ppgator[_,_,{_,_},{v1_,v2_},{_,_}]/;FreeQ[{v1,v2},"c2q"|"c2g"|"c3g"|"cqgq"|"nc"]];
If[Length[pindx]==0,
False(* if no propagator *)
,
(* discard the operator labels, all operators at same position form a single vertex *)
dias=di[[2]]/.ppgator[f1_,f2_,{xx_,yy_},{"o","o"},{n1_,n2_}]:>ppgator[f1,f2,{xx,yy},{"o","o"},{xx,yy}]/.{ppgator[f1_,f2_,{xx_,yy_},{"o",v2_},{n1_,n2_}]:>ppgator[f1,f2,{xx,yy},{"o",v2},{xx,n2}],ppgator[f1_,f2_,{xx_,yy_},{v1_,"o"},{n1_,n2_}]:>ppgator[f1,f2,{xx,yy},{v1,"o"},{n1,yy}]};
tmp=dias[[pindx[[1,1]]]];(* the first propagator *)
tmp={vs[tmp[[3,1]],tmp[[4,1]],tmp[[5,1]]], vs[tmp[[3,2]],tmp[[4,2]],tmp[[5,2]]]};

remains=Delete[dias,pindx[[1,1]]];
edges=FixedPoint[(tmp2={};
tmpdia=#//.{(* remove connected propagators, update the vertices list; remove the vertices that have already verified *)
ppgator[_,_,{x1_,x2_},{v1_,v2_},{n1_,n2_}]/;!FreeQ[tmp,vs[x1,v1,n1]]:>(tmp2=Append[tmp2,vs[x2,v2,n2]];{})
,
ppgator[_,_,{x1_,x2_},{v1_,v2_},{n1_,n2_}]/;!FreeQ[tmp,vs[x2,v2,n2]]:>(tmp2=Append[tmp2,vs[x1,v1,n1]];{})
 };tmp=tmp2;tmpdia
)&,remains(*di[[2,2;;]]*)];

(* if there still ramain propagators *)
If[FreeQ[edges,ppgator[_,_,{_,_},{v1_,v2_},{_,_}]],
True
,
False
]
]
]

(* whether an operator is connected by at least on propagator *)
momentumtransferQ[di_dia]:=Block[{tmp,tmp2,olabels},
olabels =Cases[di,ppgator[_,_,_,{"o",vv_},{nn1_,nn2_}]|ppgator[_,_,_,{vv_,"o"},{nn2_,nn1_}]:>{"o",nn1,vv},Infinity]//DeleteDuplicates;

olabels=olabels//.{aa___,{"o",nn_,vv1__},bb___,{"o",nn_,vv2__},cc___}:>{{"o",nn,vv1,vv2},aa,bb,cc};
If[!FreeQ[olabels,{"o",_,vvs__}/;FreeQ[{vvs},"o"|"v3"|"v4"|"vq"|"vc"]],
False
,
True
]
]

(* for vertex has n legs, the maxmum number of connected condensate lines is n-2 to guarantee momentum transfer *)
(* the G in DD\psi gives a qgq vertex has only one propagator, which is not allows since such situation is included in nonlocalcondensate diagram *)
condtadpoleQ[di_dia]:=Block[{tmp,tmp2,olabels},
olabels =Cases[di,ppgator[_,_,{xx1_,xx2_},{v_,vv_},{nn1_,nn2_}]|ppgator[_,_,{xx2_,xx1_},{vv_,v_},{nn2_,nn1_}]/;MatchQ[v,"vq"|"vc"|"v3"|"v4"]&&MatchQ[vv,"c2q"|"c2g"|"c3g"|"cqgq"|"nc"]:>{v,xx1,vv},Infinity];

olabels=olabels//.{aa___,{vv_,xx_,vv1__},bb___,{vv_,xx_,vv2__},cc___}:>{{vv,xx,vv1,vv2},aa,bb,cc};

If[FreeQ[olabels,({"vq"|"vc"|"v3",_,vvs__}/;Length[{vvs}]>1)|({"v4",_,vvs__}/;Length[{vvs}]>2)],
False
,
True
]
]


(* condline[G^nv, A^, {z, 0},...] vanish since A^(0) cannot give G^uv *)
condnonvainshQ[di_dia]:=Block[{tmp},
If[FreeQ[di,ppgator[_,_Gluon|_ovg,{_,0},{"c2g"|"c3g"|"cqgq"|"nc","o"},{_,_}]]&&FreeQ[di,ppgator[_Gluon|_ovg,_,{0,_},{"o","c2g"|"c3g"|"cqgq"|"nc"},{_,_}]],
True
,
False
]
]


(* if exist a propagator connected with two parts, one contains all operators while another doesn't contain operators *)
(* to do: can be more efficient *)
TadPoleQ[di_dia]:=Block[{tmp,tmp2,i,vtxlist,vtx1,vtx2,part1,part2,tf=False},
tmp=Flatten[di[[2]]/.propagator[aa_,__]:>aa];
For[i=1,i<Length[tmp]+1,i++,
{vtx1,vtx2}=Transpose[List@@tmp[[i,3;;]]];
vtx1={vtx1};
vtx2={vtx2};
tmp2=Delete[tmp,i];(* remove the selected propagator *)

part2=tmp2//.ppgator[_,_,{x1_,xx_},{v1_,vv_},{n1_,nn_}]|ppgator[_,_,{xx_,x1_},{vv_,v1_},{nn_,n1_}]/;!FreeQ[vtx1,{xx,vv,nn}]:>(vtx1=Append[vtx1,{x1,v1,n1}];{});
part2=part2//Flatten;

If[Length[part2]>0,(* if contain disjoint part *)
If[(!FreeQ[vtx1,{_,"o",_}]&&FreeQ[part2,ppgator[_,_,_,{_,"o"},_]|ppgator[_,_,_,{"o",_},_]]&&FreeQ[vtx2,{_,"o",_}])||(FreeQ[vtx1,{_,"o",_}]&&!FreeQ[part2,ppgator[_,_,_,{_,"o"},_]|ppgator[_,_,_,{"o",_},_]]&&!FreeQ[vtx2,{_,"o",_}]),
tf=True;
Break[]
]
]
];

tf
]


(* despite that every time Contract with vertices, it always pick up only the first uncontracted vertex, duplicate diagram could still generate, for example, the G^uv in hybrid Current Contract with v3(3), the remaining v3(2) can Contract with another v3(3) and vqgq, but different order of contraction generate duplicated diagram *)
(* to avoid this, one strategy is that instead of Contract the field once each time, for vn(m), list all contractable fields and enurmate all combination of m fields, contact them *)
(* a lazy way is delete the duplicate diagrams by suppress the integer label i in vertex vn(i), then such duplicated diagrams will be identical; the label n in ppgator[...,{n1_,n2_}] should not be suppressed, since for each type of vertex, the fieldcontract always pick the first uncontracted one *)
duplicatremove[dialist_List]:=Block[{tmp,operatorq},
tmp=dialist;
tmp={#/.{v3[n_]:>v3[0],v4[n_]:>v4[0],cv2[n_]:>cv2[0],cv3[n_]:>cv3[0],ncvqb[n_]:>ncvqb[0],ncvq[n_]:>ncvq[0],ncvg[n_]:>ncvg[0],ovg[inds_,n_]:>ovg[inds,0]},#}&/@tmp;
(*tmp={#[[1]]/.{ppgator[f1_,f2_,{x1_,x2_},{v1_,v2_},{n1_,n2_}]:>(operatorq=Boole[MatchQ[#,"o"]]&/@{v1,v2};ppgator[f1,f2,operatorq{x1,x2},{v1,v2},operatorq{n1,n2}])}, #[[2]]}&/@tmp;*)

(* for gluon propagators, the order of gluons is irrelevant *)
tmp={#[[1]]/.ppgator[f1_,f2_,{x1_,x2_},{v1_,v2_},{n1_,n2_}]/;FreeQ[{f1,f2},_Fermion|_AntiFermion|_Ghost|_AntiGhost|_ncvq|_ncvqb]&&!OrderedQ[{f1,f2}]:>ppgator[f2,f1,{x2,x1},{v2,v1},{n2,n1}],#[[2]]}&/@tmp;
tmp={#[[1,2]]//Sort,#[[2]]}&/@tmp;

(* gather identical diagram *)
tmp=Gather[tmp,#1[[1]]===#2[[1]]& ];
tmp=#[[1,2]]&/@tmp
]


(* for the vertex connected with two propagators and all other legs are connected with condensate, combine the two propagators to a single one *)
propagatorCombine[di_dia]:=Block[{tmp,sign, pglist,combined={},single={},condlist},
(* distinguish propagator and Condensate line *)
sign=di[[1]];
tmp=di[[2]]/.ppgator[f1_,f2_,{xx1_,xx2_},{v1_,v2_},{n1_,n2_}]/;!FreeQ[{v1,v2},"c2q"|"c2g"|"c3g"|"cqgq"|"nc"]:>condline[f1,f2,{xx1,xx2},{v1,v2},{n1,n2}];

pglist=DeleteCases[tmp,_condline]/.pg_ppgator:>propagator[{pg}/.ppgator[f1_,f2_,xx_List,vv_List,nn_List]:>ppgator[{f1,f2},{1(* the direction: f2->f1: +1, f1->f2: -1 *)},xx,vv,nn],pg[[3]],pg[[4]],pg[[5]]];
condlist=DeleteCases[tmp,_ppgator];
 (* pglist must be nonempty for diagram that connectedQ = True *)
{combined,pglist,single}=FixedPoint[combine[#[[1]],#[[2]],#[[3]],condlist]&,{combined,pglist,single}];

{sign,{combined,pglist,single},condlist}
]


combine[combined_,pgs_,single_,condlist_]:=Block[{tmp,tmplist,labels,idx,tadpole,vertexgrade},
tmp=pgs[[1]];
tmplist=pgs[[2;;]];
(* call a propagator without connecting with background fields as single propagator *)

If[And@@(MatchQ[#,{_,"o",_}]&/@Transpose[List@@tmp[[2;;]]]),
(* single propagator connected with two operator *)
{combined, tmplist, Append[single,tmp]}
,

labels=List@@tmp[[2;;]];
idx=Position[tmplist,propagator[_,labels[[1]],labels[[2]],labels[[3]]]|propagator[_,labels[[1]]//Reverse,labels[[2]]//Reverse,labels[[3]]//Reverse]];
If[Length[idx]>0,(* another propagator(s) connect with same two vertices *)

vertexgrade=Table[If[MatchQ[labels[[2,i]],"o"],0,(* the number of other propagators connected with the same internal vertex *)Length[Cases[{tmplist,single},propagator[_,{labels[[1,i]],_},{labels[[2,i]],_},{labels[[3,i]],_}]|propagator[_,{_,labels[[1,i]]},{_,labels[[2,i]]},{_,labels[[3,i]]}],Infinity]]],{i,1,2}];
(*Print[vertexgrade,Length[idx]];*)
If[MatchQ[vertexgrade,{_,Length[idx]}|{Length[idx],_}],
(* exist another single propagator(s) connect with same two vertices, and one of the two vertices connected with only these propagators, it is a tadpole *)
{Append[combined,"tadpole"],Delete[tmplist,idx],Join[single,{tmp},Part[tmplist,Flatten[idx]]]}
,
(* exist another single propagator(s) and they are connected to a loop(s), both vertices connected with other propagator(s) *)
{combined,Delete[tmplist,idx],Join[single,{tmp},Part[tmplist,Flatten[idx]]]}
]
,
{tmp,tmplist}=FixedPoint[combineonce[#[[1]],#[[2]],condlist]&,{tmp,tmplist}];

If[Length[tmp[[1]]]>1,
{Append[combined,tmp],tmplist,single}
,
{combined,tmplist,Append[single,tmp]}
]
]
]
]/;Length[pgs]>0

combine[combined_,{},single_,condlist_]:={combined,{},single}


combineonce[pg1_,pgslist_,condlist_]:=Block[{tmp1,tmp2,lr,matchp,indx,indx1,indx2,i,sign},
tmp1=pg1;
tmp2=pgslist;

lr=If[MatchQ[#,{_,"o",_}],"o",#]&/@Transpose[List@@tmp1[[2;;]]];
(* i=1: left index; i=2: right index *)
Do[(* find the propagator(s) connected to the same vertices; the Condensate must also connected with the corresponding vertex *)
If[MatchQ[lr[[i]],_List]&&!FreeQ[condlist,condline[_,_,{lr[[i,1]],_},{lr[[i,2]],_},{lr[[i,3]],_}]|condline[_,_,{_,lr[[i,1]]},{_,lr[[i,2]]},{_,lr[[i,3]]}]],
indx1=Position[tmp2,propagator[_,{lr[[i,1]],_},{lr[[i,2]],_},{lr[[i,3]],_}]];
indx2=Position[tmp2,propagator[_,{_,lr[[i,1]]},{_,lr[[i,2]]},{_,lr[[i,3]]}]];

Which[Length[indx1]==1&&Length[indx2]==0,
If[i==1,sign=-1,sign=1];
indx1=indx1[[1,1]];
matchp=momentuM[tmp2[[indx1]],sign];
tmp2=Delete[tmp2,indx1];
indx=If[i==1,
{{matchp[[2,2]],tmp1[[2,3-i]]},{ matchp[[3,2]],tmp1[[3,3-i]]},{ matchp[[4,2]],tmp1[[4,3-i]]}}
,(* the ramaining index, append it to the left or right if the matched index in tmp1 appears at left or right *)
{{tmp1[[2,3-i]], matchp[[2,2]]},{tmp1[[3,3-i]], matchp[[3,2]]},{tmp1[[4,3-i]], matchp[[4,2]]}}
];

tmp1=propagator[Join[tmp1[[1]],  matchp[[1]]],##]&@@indx
,
Length[indx1]==0&&Length[indx2]==1,
If[i==1,sign=1,sign=-1];
indx2=indx2[[1,1]];
matchp=momentuM[tmp2[[indx2]],sign];
tmp2=Delete[tmp2,indx2];
indx=If[i==1,
{{matchp[[2,1]],tmp1[[2,3-i]]},{ matchp[[3,1]],tmp1[[3,3-i]]},{ matchp[[4,1]],tmp1[[4,3-i]]}}
,(* the ramaining index, append it to the left or right if the matched index in tmp1 appears at left or right *)
{{tmp1[[2,3-i]], matchp[[2,1]]},{tmp1[[3,3-i]], matchp[[3,1]]},{tmp1[[4,3-i]], matchp[[4,1]]}}
];

tmp1=propagator[Join[tmp1[[1]],  matchp[[1]]],##]&@@indx
]
]
,{i,1,2}];


{tmp1,tmp2}
]


(* multiple propagators in propagator[] *)

propagator/:momentuM[propagator[pgs_List,xx_List,vv_List,nn_List]]:={pgs[[-1]]}/;(And@@(MatchQ[#,_ppgator|_propagator]&/@pgs[[;;-2]]))&&!MatchQ[pgs[[-1]],_ppgator|_propagator]

propagator/:momentuM[propagator[pgs_List,xx_List,vv_List,nn_List]]:={}/;(And@@(MatchQ[#,_ppgator|_propagator]&/@pgs))

propagator/:momentuM[propagator[pgs_List,xx_List,vv_List,nn_List],sign_Integer]:=propagator[Join[pgs[[;;-2]]/.{pgs[[-1]]->sign pgs[[-1]],-pgs[[-1]]->-sign pgs[[-1]]},{sign pgs[[-1]]}],xx,vv,nn]/;(And@@(MatchQ[#,_ppgator|_propagator]&/@pgs[[;;-2]]))&&!MatchQ[pgs[[-1]],_ppgator|_propagator]

propagator/:momentuM[propagator[pgs_List,xx_List,vv_List,nn_List],replace_Rule]:=propagator[pgs/.replace,xx,vv,nn](*/;(And@@(MatchQ[#,_ppgator|_propagator]&/@pgs[[;;-2]]))&&!MatchQ[pgs[[-1]],_ppgator|_propagator]*)

(*propagator/:momentuM[propagator[{pgs__ppgator|pgs__propagator},xx_List,vv_List,nn_List],k:Except[_Integer|_Rule]]:=propagator[Append[{pgs},k],xx,vv,nn]/;Length[{pgs}]>1*)

propagator/:momentuM[propagator[pgs_List,xx_List,vv_List,nn_List],k:Except[_Integer|_Rule]]:=If[MatchQ[pgs[[-1]],_Integer],
propagator[Join[pgs[[;;-2]],{k pgs[[-1]]}],xx,vv,nn]
,(* pgs[[-1]] can be composite; the replacement a+b-c+d-e/.a-c+d -> g works *)
If[MatchQ[pgs[[-1]],_ppgator|_propagator],
propagator[Join[ppgator[#[[1]],#[[2]]k,#[[3]],#[[4]],#[[5]]]&/@pgs,{k}],xx,vv,nn]
,
propagator[Join[pgs[[;;-2]]/.{pgs[[-1]]->k,-pgs[[-1]]->-k},{k}],xx,vv,nn]
]
]/;Length[pgs]>1&&(And@@(MatchQ[#,_ppgator|_propagator]&/@pgs[[;;-2]]))(*&&!MatchQ[pgs[[-1]],_ppgator|_propagator]*)

(* automatic assign momentums for the propagators in composite propagator, if the momentum haven't speificed *)
propagator[{pgs__ppgator,kk:Except[_ppgator]},xx_,vv_,nn_]:=propagator[Append[momentuM[#,kk]&/@{pgs},kk],xx,vv,nn]/;(And@@(IntegerQ[#[[2,1]]]&/@{pgs}))


(* single propagator in propagator[] *)
propagator/:momentuM[propagator[{pgs_ppgator},xx_List,vv_List,nn_List],sign_Integer]:=propagator[{ppgator[pgs[[1]],sign pgs[[2]], pgs[[3]],pgs[[4]],pgs[[5]]]},xx,vv,nn]
propagator/:momentuM[propagator[{pgs_ppgator},xx_List,vv_List,nn_List],replace_Rule]:=propagator[{ppgator[pgs[[1]],pgs[[2]]/.replace, pgs[[3]],pgs[[4]],pgs[[5]]]},xx,vv,nn]

propagator/:momentuM[propagator[{pgs_ppgator},xx_List,vv_List,nn_List],k:Except[_Integer|_Rule]]:=If[IntegerQ[pgs[[2,1]]],
propagator[{ppgator[pgs[[1]],k pgs[[2]], pgs[[3]],pgs[[4]],pgs[[5]]]},xx,vv,nn],
propagator[{ppgator[pgs[[1]],{k}, pgs[[3]],pgs[[4]],pgs[[5]]]},xx,vv,nn]
]

ppgator/:momentuM[pg_ppgator,sign_Integer]:=ppgator[pg[[1]],sign pg[[2]],pg[[3]],pg[[4]],pg[[5]]]
ppgator/:momentuM[pg_ppgator,replace_Rule]:=ppgator[pg[[1]],pg[[2]]/.replace,pg[[3]],pg[[4]],pg[[5]]]
ppgator/:momentuM[pg_ppgator,k:Except[_Integer|_Rule]]:=If[IntegerQ[pg[[2,1]]],
ppgator[pg[[1]],k pg[[2]],pg[[3]],pg[[4]],pg[[5]]],
ppgator[pg[[1]],{k},pg[[3]],pg[[4]],pg[[5]]]
]

(*-----------------------------*)
propagator/:subpropagatorS[propagator[pgs_List,xx_List,vv_List,nn_List]]:=Cases[pgs,_ppgator|_looppgator]


(* asign the momentums when two propagators are joined *)
(* sign = 1/-1 if pg2 parallel/antiparallel with pg1 *)
momentumcombine[pg1_,pg2_,labels_,sign_,reverse_:False]:=Block[{tmp,pgs1,pgs2,pk1,pk2},
If[MatchQ[pg1[[-1]],_ppgator|_propagator],
pk1={};
pgs1=pg1
,
pk1={pg1[[-1]]};
pgs1=pg1[[;;-2]]
];

If[MatchQ[pg2[[-1]],_ppgator|_propagator],
pk2={};
pgs2=pg2
,
pk2={pg2[[-1]]};
pgs2=pg2[[;;-2]]
];

(*Print[pk1,pk2];*)
Which[Length[pk1]==1&&Length[pk2]==0,
propagator[Append[Join[pgs1,momentuM[#,sign pk1[[1]] ]&/@pgs2],pk1[[1]]],##]&@@labels
,
Length[pk1]==0&&Length[pk2]==1,
If[sign==1,
propagator[Append[Join[momentuM[#,pk2[[1]] ]&/@pgs1,pgs2],pk2[[1]]],##]&@@labels
,(* flip the sign of pk2[[1]] *)
propagator[Append[Join[momentuM[#,-pk2[[1]] ]&/@pgs1,momentuM[#,-1 ]&/@pgs2],-pk2[[1]]],##]&@@labels
]
,
Length[pk1]==1&&Length[pk2]==1,
If[reverse,
propagator[Append[Join[momentuM[#,pk1[[1]]->sign pk2[[1]] ]&/@pgs1,pgs2],pk1[[1]]],##]&@@labels
,
propagator[Append[Join[pgs1,momentuM[#,pk2[[1]]->sign pk1[[1]] ]&/@pgs2],pk1[[1]]],##]&@@labels
]
]
]
(* combine the joined propagators into one composite progator; the external momentum may be replaced after join two propagators, but the internal loop momentums are not touched in this process *)
propconnect[pgs_List]:=Block[{tmp,tmppgator,proplist={},sublist,labels,prop1,prop2,pk1,pk2,pk3},
tmp=pgs;

(* combine two propagators connected at a same vertex, if they are connected at operator-vertex, do not combine them since external momentum flow into the operator-vertex *)
tmp=FixedPoint[Replace[#,{{aa___,propagator[ps1_,{xx1_,xx_},{vv1_,vv_},{nn1_,nn_}],bb___,propagator[ps2_,{xx_,xx2_},{vv_,vv2_},{nn_,nn2_}],cc___}/;FreeQ[{aa,bb,cc},propagator[_,{xx,_},{vv,_},{nn,_}]|propagator[_,{_,xx},{_,vv},{_,nn}]]&&vv=!="o":>((*Print["1-->",ps1," >> ",ps2," >> ",momentumcombine[ps1,ps2,{{xx1,xx2},{vv1,vv2},{nn1,nn2}},1]];*){aa,momentumcombine[ps1,ps2,{{xx1,xx2},{vv1,vv2},{nn1,nn2}},1],bb,cc})
,
{aa___,propagator[ps2_,{xx_,xx2_},{vv_,vv2_},{nn_,nn2_}],bb___,propagator[ps1_,{xx1_,xx_},{vv1_,vv_},{nn1_,nn_}],cc___}/;FreeQ[{aa,bb,cc},propagator[_,{xx,_},{vv,_},{nn,_}]|propagator[_,{_,xx},{_,vv},{_,nn}]]&&vv=!="o":>((*Print["2-->",ps1," >> ",ps2," >> ",momentumcombine[ps1,ps2,{{xx1,xx2},{vv1,vv2},{nn1,nn2}},1]];*){aa,momentumcombine[ps1,ps2,{{xx1,xx2},{vv1,vv2},{nn1,nn2}},1],bb,cc})
,
{aa___,propagator[ps1_,{xx1_,xx_},{vv1_,vv_},{nn1_,nn_}],bb___,propagator[ps2_,{xx2_,xx_},{vv2_,vv_},{nn2_,nn_}],cc___}/;FreeQ[{aa,bb,cc},propagator[_,{xx,_},{vv,_},{nn,_}]|propagator[_,{_,xx},{_,vv},{_,nn}]]&&vv=!="o":>((*Print["3-->",ps1," >> ",ps2," >> ",momentumcombine[ps1,ps2,{{xx1,xx2},{vv1,vv2},{nn1,nn2}},1]];*){aa,momentumcombine[ps1,ps2,{{xx1,xx2},{vv1,vv2},{nn1,nn2}},-1],bb,cc})
,
{aa___,propagator[ps1_,{xx_,xx1_},{vv_,vv1_},{nn_,nn1_}],bb___,propagator[ps2_,{xx_,xx2_},{vv_,vv2_},{nn_,nn2_}],cc___}/;FreeQ[{aa,bb,cc},propagator[_,{xx,_},{vv,_},{nn,_}]|propagator[_,{_,xx},{_,vv},{_,nn}]]&&vv=!="o":>((*Print["4-->",ps1," >> ",ps2," >> ",momentumcombine[ps1,ps2,{{xx1,xx2},{vv1,vv2},{nn1,nn2}},1]];*){aa,momentumcombine[ps1,ps2,{{xx1,xx2},{vv1,vv2},{nn1,nn2}},-1,True],bb,cc})
},{0}]&,
tmp]
]


noloop[pgs_List]:=Block[{tmp,lk},
tmp=pgs/.propagator[fds_List,labels__]:>propagator[Append[fds,Unique["k$"]],labels];(* set arbitrary momentum initially *)

If[Length[tmp]==1,
tmp
,
tmp=propconnect[tmp];
If[Length[tmp]==1,(* if the propagators can combine into one composite propagator *)
tmp
,
False
]
]
]


(* specify the momentums for each propagator, recursively find the loops, after finding a length-n loop, some of them connected with other propagators becomes a composite proapgator, *)
MomentumSpecify[di_dia]:=Block[{tmp,tree,sign,condlist,vlist,momentumlist={},loops,remain},
sign=di[[1]];
tmp=di[[2]]//Flatten;
condlist=di[[3]];
tree=noloop[tmp];

If[tree=!=False,
{sign,tree[[1,2]],tree,condlist}

,(* contain loops *)
tmp=loopfinder[tmp,3];
momentumlist=Join[momentumlist,tmp[[2]]];
tmp=propconnect[Join[tmp[[1]],tmp[[3]]]];

tmp=loopfinder[tmp,2];
momentumlist=Join[momentumlist,tmp[[2]]];
tmp=propconnect[Join[tmp[[1]],tmp[[3]]]];

(* give up more complicate structure, merely write delta[k1+ k2+ ...] for each vertex and k1, k2, ... to each propagator *)
(*If[!FreeQ[tmp,ppgator[{},{_Integer},__]],


];*)
If[!FreeQ[tmp,propagator[_,_,{xx_,xx_},{vv_,vv_},{nn_,nn_}]],(* detect tadpole *)
{sign,{"tadpole"},tmp,condlist}
,
{sign,momentumlist,tmp,condlist}
]
]
]
(* find all loops with n propagators; for propagator[_,{x,y},_,_], the momentum along y->x is positive *)
loopfinder[pgs_List,n_]:=Block[{tmp,indx1,indx2,len3loop,adjacentlen3loop={},labels,labellist,nprop,loopprop,looplist={},proplist={},loopmomentum,exmomentum,momentumlist,integralmomentum={},twoloop},
tmp=pgs;

If[n==2,
While[Length[tmp]>1,
labels=List@@tmp[[1,2;;]];
indx1=Position[tmp[[2;;]],propagator[_,labels[[1]],labels[[2]],labels[[3]]],1];(* parallel propagators *)
indx2=Position[tmp[[2;;]],propagator[_,labels[[1]]//Reverse,labels[[2]]//Reverse,labels[[3]]//Reverse],1];(* anti parallel propagators *)
nprop=1+Length[indx1]+Length[indx2];
If[nprop>1,
loopmomentum=Table[Unique["k$"],nprop];
momentumlist=loopmomentum;(* set loopmomentum[[1]] as external momentum, momentumlist[[2;;]] as loop momentums *)
integralmomentum=Append[integralmomentum,momentumlist[[2;;]]//Reverse];
loopmomentum=Flatten[{loopmomentum[[1]]+loopmomentum[[2]],loopmomentum[[3;;]]-loopmomentum[[2;;-2]],-loopmomentum[[-1]]}];
loopmomentum=Table[If[i<=1+Length[indx1],1,-1],{i,nprop}]loopmomentum;(* times -1 for anti parallel propagatots *)
(* Print[momentumlist]; *)
loopprop=Prepend[Extract[tmp[[2;;]],Join[indx1,indx2]],tmp[[1]]](* /.propagator[{pg_ppgator},xx_,vv_,nn_]:>pg *);
(*Print["<><><>",loopprop];*)
loopprop=Table[momentuM[loopprop[[i]],loopmomentum[[i]]],{i,1,nprop}];
(* combine to a composite propagator, clean the unnecessary nest structure *)
(* Print[loopprop,"--->",momentumlist];*)
looplist=Append[looplist,propagator[Append[loopprop,momentumlist[[1]]],labels[[1]],labels[[2]],labels[[3]]]];
(*Print[">:>>>",looplist];*)
tmp=Delete[tmp[[2;;]],Join[indx1,indx2]]
,
proplist=Append[proplist,tmp[[1]]];
tmp=tmp[[2;;]]
]
]

,
If[n==3,
While[Length[tmp]>2,
len3loop=loopn[tmp,3];

If[Length[len3loop]==0,(* tmp[[1]] is not contained in a length-3 loop *)
proplist=Append[proplist,tmp[[1]]];
tmp=tmp[[2;;]]
,

If[Length[len3loop]==1,
(* find the length-3 loops for second and third propagators in the loop of len3loop; discard other two propagators for each case to avoid obtain duplicated length-3 loop *)
adjacentlen3loop=Join[loopn[Prepend[len3loop[[1,-1]],len3loop[[1,2,2]]],3],loopn[Prepend[len3loop[[1,-1]],len3loop[[1,3,2]]],3] ]
];
Print[len3loop,"____",adjacentlen3loop];
(* the two loops sharing a same propagator *)
If[Length[len3loop]==2||Length[adjacentlen3loop]==1,
twoloop=twoloopprogator[len3loop,adjacentlen3loop];
If[Length[twoloop[[1]]]==0,
(* the two loops cannot be comibined as a propagator *)
Message[QCDSumDiagram::loopwan];
proplist=Join[proplist,twoloop[[3]]];
tmp={}
,
proplist=Append[proplist,twoloop[[1]]];
integralmomentum=Append[integralmomentum,twoloop[[2]]];
tmp=twoloop[[3]]
]
,
(* give up assigning loop momentum for complicat cases *)
Message[QCDSumDiagram::loopwan];
proplist=Join[proplist,tmp];
tmp={}
];
]
]
,
(* for multiple propagators connected to a loop, it gives a multiple legs vertex, generate a new vertex labels, for the propagators connected with this loop, change the vertex labels ...{xx_,yy_},{v1_,v2_},{n1_,n2_}] to connect the new vertex *)
(* to do: assign momentums for multi-propagators loop *)
Message[QCDSumDiagram::loopwan];
proplist=Join[proplist,tmp];
tmp={}
]
];

{looplist,integralmomentum,Join[tmp,proplist]}
]


(* find all length-n loops that the first propagator appears in *)
loopn[pgs_List,n_,inner_:True]:=Block[{list,start,tmp,labels,pathlist,connected,remain,i},
pathlist={{{1,pgs[[1]]},pgs[[2;;]]}};
(* start form the right vertex of the first propagator *)
(* for large n, it's better to store the paths as a tree: {{propagator_1, {{propagator_a, {propagators ...}}, propagator_b, the proagators connected with right end of propagator_1}}, ...}   *)

(* store each propagatos as {1, propagator} or {-1, propagator}, depending on the direction *)
Do[For[i=1,i<Length[pathlist]+1,i++,
If[pathlist[[i,-2,1]]==1,
labels={pathlist[[i,-2,2,2,2]],pathlist[[i,-2,2,3,2]],pathlist[[i,-2,2,4,2]]}
,
labels={pathlist[[i,-2,2,2,1]],pathlist[[i,-2,2,3,1]],pathlist[[i,-2,2,4,1]]}
];(* the vertex hasn't been connected for last propagator *)

(* the propagators connected with the last propagator *)
connected = Join[Cases[pathlist[[i,-1]],propagator[pg_,{labels[[1]],xx_},{labels[[2]],vv_},{labels[[3]],nn_}]:>{1,propagator[pg,{labels[[1]],xx},{labels[[2]],vv},{labels[[3]],nn}]}]
,
Cases[pathlist[[i,-1]],propagator[pg_,{xx_,labels[[1]]},{vv_,labels[[2]]},{nn_,labels[[3]]}]:>{-1,propagator[pg,{xx,labels[[1]]},{vv,labels[[2]]},{nn,labels[[3]]}]}]];

connected = {#,DeleteCases[pathlist[[i,-1]],#[[2]]]}&/@connected ;
(* update the path list *)
pathlist[[i]]=Join[pathlist[[i,;;-2]],#]&/@connected;

];
(* for each sublist in pathlist, the first n elemetns are {+-1, propagator}, which are the propagators of the length=n loop, and the last element is the ramining propagators *)
pathlist=Flatten[pathlist,1];
,
n-1];
(* ignore the propagator connecting two operators directely *)
If[inner===True,
pathlist=DeleteCases[#,{_Integer,propagator[_,{_,_},{"o","o"},{_,_}]}]&/@pathlist;
];

(* delete the path has length < n *)
pathlist=DeleteCases[pathlist,aa_/;Length[Cases[aa,{_Integer,_propagator}]]<n];
(* for the path connecting to a single loop, each verteix appears twice *)
DeleteCases[pathlist,aa_/;!(And@@(Length[#]==2&/@Gather[Flatten[Cases[aa[[1;;n]],{_,propagator[_,{x1_,x2_},{v1_,v2_},{n1_,n2_}]}:>{{x1,v1,n1},{x2,v2,n2}}],1]]))
]
]


(* for two length-3 loops sharing a same propagator, combine them as a propagator if possible *)
(* if loopn[..., 3] gives a length=2 list, the two length=3 loop sharing a same propagator, the duplicate part in the remaining fields in two sublist gives the propagators not in these two length=3 loops *)
(* if loopn[..., 3] gives a length=1 list, and adjacentlen3loop gives length=1 list, the remaining fields in adjacentlen3loop gives the propagators not in these two length=3 loops *)
twoloopprogator[pgs1_List,pgs2_:{}]:=Block[{tmp,loop1,loop2,remain,commonprogator,labels,ver1,ver2,interv,exmomentum=Unique["k$"],loopk1=Unique["k$"],loopk2=Unique["k$"]},
(* discard the +-1 label for each propagator *)
If[Length[pgs2]==0,
loop1=pgs1[[1,;;-2]]/.{_Integer,pg_propagator}:> pg;
loop2=pgs1[[2,;;-2]]/.{_Integer,pg_propagator}:> pg;
commonprogator=pgs1[[1,1,2]];
remain=DeleteCases[pgs1[[1,-1]], aa_/;!FreeQ[loop2,aa]]
,
loop1=pgs1[[1,;;-2]]/.{_Integer,pg_propagator}:> pg;
loop2=pgs2[[1,;;-2]]/.{_Integer,pg_propagator}:> pg;
commonprogator=loop2[[1]];
remain=pgs2[[1,-1]]
];

labels=Transpose[List@@commonprogator[[2;;]]];

(* the vertices of the commonpropagator should not exist in the remaining propagators *)
If[FreeQ[remain,propagator[_,{labels[[1,1]],_},{labels[[1,2]],_},{labels[[1,3]],_}]|propagator[_,{labels[[2,1]],_},{labels[[2,2]],_},{labels[[2,3]],_}]|propagator[_,{_,labels[[1,1]]},{_,labels[[1,2]]},{_,labels[[1,3]]}]|propagator[_,{_,labels[[2,1]]},{_,labels[[2,2]]},{_,labels[[2,3]]}]],
(* the two vertices that can be view as the two vertices of the two-loop propagator *)

ver1=Flatten[Transpose[List@@#[[2;;]]]&/@loop1,1];
ver1=DeleteDuplicates[DeleteCases[DeleteCases[ver1,labels[[1]]],labels[[2]]]][[1]];

ver2=Flatten[Transpose[List@@#[[2;;]]]&/@loop2,1];
ver2=DeleteDuplicates[DeleteCases[DeleteCases[ver2,labels[[1]]],labels[[2]]]][[1]];
loop1=DeleteCases[loop1,commonprogator];
loop2=DeleteCases[loop2,commonprogator];

(*loop1=loop1/.{propagator[pg_,{ver1[[1]],labels[[1,1]]},{ver1[[2]],labels[[1,2]]},{ver1[[3]],labels[[1,3]]}]:>looppgator[Append[pg,exmomentum+loopk1],{ver1[[1]],labels[[1,1]]},{ver1[[2]],labels[[1,2]]},{ver1[[3]],labels[[1,3]]}],
propagator[pg_,{labels[[1,1]],ver1[[1]]},{labels[[1,2]],ver1[[2]]},{labels[[1,3]],ver1[[3]]}]:>looppgator[Append[pg,-exmomentum-loopk1],{labels[[1,1]],ver1[[1]]},{labels[[1,2]],ver1[[2]]},{labels[[1,3]],ver1[[3]]}]
,
propagator[pg_,{labels[[2,1]],ver1[[1]]},{labels[[2,2]],ver1[[2]]},{labels[[2,3]],ver1[[3]]}]:>looppgator[Append[pg,loopk1],{labels[[2,1]],ver1[[1]]},{labels[[2,2]],ver1[[2]]},{labels[[2,3]],ver1[[3]]}],
propagator[pg_,{ver1[[1]],labels[[2,1]]},{ver1[[2]],labels[[2,2]]},{ver1[[3]],labels[[2,3]]}]:>looppgator[Append[pg,-loopk1],{ver1[[1]],labels[[2,1]]},{ver1[[2]],labels[[2,2]]},{ver1[[3]],labels[[2,3]]}]
};

loop2=loop2/.{
propagator[pg_,{labels[[1,1]],ver2[[1]]},{labels[[1,2]],ver2[[2]]},{labels[[1,3]],ver2[[3]]}]:>looppgator[Append[pg,exmomentum+loopk2],{labels[[1,1]],ver2[[1]]},{labels[[1,2]],ver2[[2]]},{labels[[1,3]],ver2[[3]]}],propagator[pg_,{ver2[[1]],labels[[1,1]]},{ver2[[2]],labels[[1,2]]},{ver2[[3]],labels[[1,3]]}]:>looppgator[Append[pg,-exmomentum-loopk2],{ver2[[1]],labels[[1,1]]},{ver2[[2]],labels[[1,2]]},{ver2[[3]],labels[[1,3]]}]
,
propagator[pg_,{ver2[[1]],labels[[2,1]]},{ver2[[2]],labels[[2,2]]},{ver2[[3]],labels[[2,3]]}]:>looppgator[Append[pg,loopk2],{ver2[[1]],labels[[2,1]]},{ver2[[2]],labels[[2,2]]},{ver2[[3]],labels[[2,3]]}],
propagator[pg_,{labels[[2,1]],ver2[[1]]},{labels[[2,2]],ver2[[2]]},{labels[[2,3]],ver2[[3]]}]:>looppgator[Append[pg,-loopk2],{labels[[2,1]],ver2[[1]]},{labels[[2,2]],ver2[[2]]},{labels[[2,3]],ver2[[3]]}]
};
(* combine to a composite propagator, clean the unnecessary nest structure *)
commonprogator=looppgator[{commonprogator[[1,1]],loopk1-loopk2},commonprogator[[2]],commonprogator[[3]],commonprogator[[4]]];
Print[loop1,"    ---   ",loop2,"   --    ",commonprogator];
{propagator[Append[Join[loop1,loop2,{commonprogator}],exmomentum],{ver1[[1]],ver2[[1]]},{ver1[[2]],ver2[[2]]},{ver1[[3]],ver2[[3]]}]/.looppgator[{pg__ppgator,kk:Except[_ppgator]},xlist_List,vlist_List,nlist_List]/;(And@@(IntegerQ[#[[2,1]]]&/@(pg))):>(momentuM[#,kk]&/@{pg}),{{loopk1,loopk2}},
remain
}*)

loop1=loop1/.{propagator[pg_,{ver1[[1]],labels[[1,1]]},{ver1[[2]],labels[[1,2]]},{ver1[[3]],labels[[1,3]]}]:>propagator[Append[pg,exmomentum+loopk1],{ver1[[1]],labels[[1,1]]},{ver1[[2]],labels[[1,2]]},{ver1[[3]],labels[[1,3]]}],
propagator[pg_,{labels[[1,1]],ver1[[1]]},{labels[[1,2]],ver1[[2]]},{labels[[1,3]],ver1[[3]]}]:>propagator[Append[pg,-exmomentum-loopk1],{labels[[1,1]],ver1[[1]]},{labels[[1,2]],ver1[[2]]},{labels[[1,3]],ver1[[3]]}]
,
propagator[pg_,{labels[[2,1]],ver1[[1]]},{labels[[2,2]],ver1[[2]]},{labels[[2,3]],ver1[[3]]}]:>propagator[Append[pg,loopk1],{labels[[2,1]],ver1[[1]]},{labels[[2,2]],ver1[[2]]},{labels[[2,3]],ver1[[3]]}],
propagator[pg_,{ver1[[1]],labels[[2,1]]},{ver1[[2]],labels[[2,2]]},{ver1[[3]],labels[[2,3]]}]:>propagator[Append[pg,-loopk1],{ver1[[1]],labels[[2,1]]},{ver1[[2]],labels[[2,2]]},{ver1[[3]],labels[[2,3]]}]
};

loop2=loop2/.{
propagator[pg_,{labels[[1,1]],ver2[[1]]},{labels[[1,2]],ver2[[2]]},{labels[[1,3]],ver2[[3]]}]:>propagator[Append[pg,exmomentum+loopk2],{labels[[1,1]],ver2[[1]]},{labels[[1,2]],ver2[[2]]},{labels[[1,3]],ver2[[3]]}],propagator[pg_,{ver2[[1]],labels[[1,1]]},{ver2[[2]],labels[[1,2]]},{ver2[[3]],labels[[1,3]]}]:>propagator[Append[pg,-exmomentum-loopk2],{ver2[[1]],labels[[1,1]]},{ver2[[2]],labels[[1,2]]},{ver2[[3]],labels[[1,3]]}]
,
propagator[pg_,{ver2[[1]],labels[[2,1]]},{ver2[[2]],labels[[2,2]]},{ver2[[3]],labels[[2,3]]}]:>propagator[Append[pg,loopk2],{ver2[[1]],labels[[2,1]]},{ver2[[2]],labels[[2,2]]},{ver2[[3]],labels[[2,3]]}],
propagator[pg_,{labels[[2,1]],ver2[[1]]},{labels[[2,2]],ver2[[2]]},{labels[[2,3]],ver2[[3]]}]:>propagator[Append[pg,-loopk2],{labels[[2,1]],ver2[[1]]},{labels[[2,2]],ver2[[2]]},{labels[[2,3]],ver2[[3]]}]
};
(* combine to a composite propagator *)
commonprogator=propagator[{commonprogator[[1,1]],loopk1-loopk2},commonprogator[[2]],commonprogator[[3]],commonprogator[[4]]];
(*Print[loop1,"    ---   ",loop2,"   --    ",commonprogator];*)
{propagator[Append[Join[loop1,loop2,{commonprogator}],exmomentum],{ver1[[1]],ver2[[1]]},{ver1[[2]],ver2[[2]]},{ver1[[3]],ver2[[3]]}],{{loopk1,loopk2}},
remain
}
,
Message[QCDSumDiagram::loopwan];
{{},{},Join[loop1,loop2,remain]}
]
]


VertexCreator[pgs_List]:=Block[{tmp,vertex,pglist,oplist,condlist,ncvtx,sunflist={},vlist,positionidx,bkgluon,lorindx,lor1,lor2,diracindx,diraci,diracj,sunindx,suna,sunb,sunc,sunn,i,flavor1,flavor2,pgx=False},
tmp=pgs;
(* assign the momentum if the if hasn't been specified, label it differently so that the momentum will be added into the vertex and give a delta[k1+k2+...] later *)
tmp=tmp/.ppgator[{f1_Fermion,f2_AntiFermion},{k_Integer},{x1_,x2_},{v1_,v2_},{n1_,n2_}]:>ppgator[{f1,f2},{{k Unique["k$"]}},{x1,x2},{v1,v2},{n1,n2}];
(*tmp=tmp/.ppgator[fields_List,labels___]:>ppgator[fields/.{v3[_Integer]|v4[_Integer]:>Gluon[Unique["lor$"],Unique["col$"]]},labels];*)
(* recover the gluon fields since the contraction has been done *)
tmp=tmp/.{v3[_Integer]|v4[_Integer]:>Gluon[Unique["lor$"],Unique["col$"]]};

(* extract the ovg *)
sunindx=Cases[tmp,ovg[{lo1_,lo2_,sun_},_]:>ovg[{lo1,lo2,sun}],Infinity]//DeleteDuplicates;
(* the rule that recover gf^nab A^a_u A^b_v *)
sunindx=sunindx/.ovg[{lo1_,lo2_,sun_}]:>(suna=Unique["Col$"];sunb=Unique["Col$"];lor1=Unique["lor$"];lor2=Unique["lor$"];sunflist=Append[sunflist,gStrong SUNF[sun,suna,sunb]];{ovg[{lo1,lo2,sun},1]->ovg[{lo1,lo2,sun},lor1,suna],ovg[{lo1,lo2,sun},2]->ovg[{lo1,lo2,sun},lor2,sunb]});
sunindx=Flatten[sunindx];
(* for ovg in propagators, recover the gluon via ovg[{lo1,lo2,sun},lor1,suna] -> Gluon[lor1,suna], for ovg in vertex, recover the vertex via: {..., ovg[{lo1,lo2,n},lor1,a], ovg[{lo1,lo2,n},lor2,b], ...} -> g f^nab (g_lo1lor1 g_lo2lor2 - g_lo1lor2 g_lo2lor1) *)

(* propagator list *)
pglist=tmp/.sunindx/.ovg[{__},lora_,sun_]:>Gluon[lora,sun];

(* recover the ovg later *)
tmp=tmp/.sunindx;

(* assign the indices for the fields in Condensate, add a -1 sign for <q^bar q> <q^bar G q>, since in both ppgator and conline the order of fields are q q^bar *)
tmp=tmp/.condline[f1_,f2_,{x1_,x2_},{v1_,v2_},{n1_,n2_}]:>If[Length[f1]==0,condline[f2,f2,{x1,x2},{v1,v2},{n1,n2}],condline[f1,f1 ,{x1,x2},{v1,v2},{n1,n2}]];
(* the overall sign is recorded in the dia, no other signs needed *)


vlist=Join[
Cases[tmp,ppgator[{f1_,f2_},k_,{x1_,x2_},{v1_,v2_},{n1_,n2_}]:>{vertex[v1,x1,n1,{f1,k}],vertex[v2,x2,n2,{f2,-k}]}
(* the momentum flow into the vertex is positive *)
(* the momentum needed for f^abc d_uA_v A^u A^v vertex, d^uA^v -d^vA^u term in G^uv, and the case that the momentum in propagators hasn't been assigned *)
,
Infinity]
,
(* for nonlocal Condensate, reacord each field's position *)
Cases[tmp,condline[f1_,f2_,{x1_,x2_},{v1_,v2_},{n1_,n2_}]:>{vertex[v1,x1,n1,{If[MatchQ[v1,"nc"],FieldsPut[f1,x2],f1],0}],vertex[v2,x2,n2,{If[MatchQ[v2,"nc"],FieldsPut[f2,x1],f2],0}]},Infinity]
]//Flatten;


(* gather the vertex; extract the indices from the fields *)
vlist=Gather[vlist,#1[[;;3]]===#2[[;;3]]&];
vlist=vlist/.{vertex[vv:Except["nc"],xx_,nn_,fk_],vers__}:>vertex[vv,xx,nn,Sort[#[[4]]&/@{vertex[vv,xx,nn,fk],vers},Order[#1[[1]],#2[[1]]]&]];


(* recover the ovg *)
vlist=vlist/.vertex["o",xx_,nn_,fk_List]:>vertex["o",xx,nn,fk//.{faa___,{ovg[{lo1_,lo2_,sun_},lr1_,sua_],k1_},fbb___,{ovg[{lo1_,lo2_,sun_},lr2_,sub_],k2_},fcc___}:>{gStrong SUNF[sun,sua,sub](MTD[lo1,lr1]MTD[lo2,lr2]-MTD[lo1,lr2]MTD[lo2,lr1]),faa,fbb,fcc}];(* discard the momentum for A^u A^v in gf^nab A^a_u A^b_v *)

(* all momentum are incoming; the connecting gluons are perturbative, when connecting with background gluons, the gluon vertices will be different *)
vlist=vlist/.{vertex["vq",xx_,nn_,{{AntiFermion[_,ii_,aa_,___],k1_},{Fermion[_,jj_,bb_,___],k2_},{Gluon[lr1_,sn_,___],k3_}}]/;k3=!=0:>vertex["vq",xx,nn,{I gStrong DCHN[GAD[lr1],ii,jj]SUNTF[sn,aa,bb]}]
,
vertex["vc",xx_,nn_,{{AntiGhost[cola_,{{lor_}}],k1_},{Ghost[colb_],k3_},{Gluon[lor_,sn2_,___],k2_}}]:>vertex["vc",xx,nn,{-gStrong  SUNF[sn2,cola,colb]FVD[k1,lor]}]
,
vertex["v3",xx_,nn_,{{Gluon[lr1_,sn1_,___],k1_},{Gluon[lr2_,sn2_,___],k2_},{Gluon[lr3_,sn3_,___],k3_}}]/;FreeQ[{k1,k2,k3},0]:>vertex["v3",xx,nn,{-gStrong SUNF[sn1,sn2,sn3](FVD[k2-k3,lr1]MTD[lr2,lr3]+FVD[k3-k1,lr2]MTD[lr3,lr1]+FVD[k1-k2,lr3]MTD[lr2,lr1])}]
,
vertex["v4",xx_,nn_,{{Gluon[lr1_,sn1_,___],k1_},{Gluon[lr2_,sn2_,___],k2_},{Gluon[lr3_,sn3_,___],k3_},{Gluon[lr4_,sn4_,___],k4_}}]/;FreeQ[{k1,k2,k3,k4},0]:>(sunn=Unique["Col$"];vertex["v4",xx,nn,{-I gStrong^2(SUNF[sn1,sn2,sunn]SUNF[sn3,sn4,sunn](MTD[lr1,lr3]MTD[lr2,lr4]-MTD[lr1,lr4]MTD[lr2,lr3])+SUNF[sn1,sn3,sunn]SUNF[sn2,sn4,sunn](MTD[lr1,lr2]MTD[lr3,lr4]-MTD[lr1,lr4]MTD[lr2,lr3])+SUNF[sn1,sn4,sunn]SUNF[sn3,sn2,sunn](MTD[lr1,lr3]MTD[lr2,lr4]-MTD[lr1,lr2]MTD[lr3,lr4]))}])
};

(* the case that connected with condensate *)
(* condline can be discard since the information about Condensate is already recorded in vertex; times -1 for q^bar Gq and q^bar q because all fermion fields are connected as \psi \bar{\psi} previously *)
oplist=Cases[vlist,vertex["o",__]];
condlist=Cases[vlist,vertex["c2q"|"c2g"|"c3g"|"cqgq"|"nc",__]];
vlist=DeleteElements[vlist,Join[oplist,condlist]];

(*-------------------------------------------------------------*)
(* below is temporarily  *)
(* gluon propagator with one background gluons *)
positionidx=Position[vlist,vertex["v3",_,_,gAlist_/;!FreeQ[gAlist,{_Gluon,0}]]];

If[Length[positionidx]==1,
bkgluon=vlist[[positionidx[[1,1]]]];
bkgluon=BackGroundGluon[bkgluon];
vlist[[positionidx[[1,1]]]]=bkgluon[[1]];
condlist=condlist/.bkgluon[[2]]
];

(* quark propagator with backgluond gluons *)
If[!FreeQ[condlist,vertex["nc",__]],
pgx=True
];

positionidx=Position[vlist,vertex["vq",_,_,{_,_,{_Gluon,0}}]];
If[Length[positionidx]>0,
For[i=1,i<Length[positionidx]+1,i++,
bkgluon=vlist[[positionidx[[i,1]]]];

lorindx=bkgluon[[4,3,1,1]];
sunn=bkgluon[[4,3,1,2]];

vlist[[positionidx[[i,1]]]]=bkgluon/.vertex["vq",xx_,nn_,{{AntiFermion[_,ii_,aa_,___],k1_},{Fermion[_,jj_,bb_,___],k2_},{Gluon[lr1_,sn_,___],k3_}}]:>vertex["vq",xx,nn,{1/2I gStrong DCHN[GAD[lr1],ii,jj]SUNTF[sn,aa,bb]}];(* interpret vertex; A(x) -> 1/2 x G(0) *)

lor1=Unique["lor$"];
pglist=pglist/.{pg1___,ppgator[{bkgluon[[4,2,1]],af_AntiFermion},kk_,{bkgluon[[2]],x2_},{"vq",vv_},{bkgluon[[3]],nn_}],pg2___}:>
If[pgx,
{pg1,FVD[bkgluon[[2]],lor1],ppgatorX[{bkgluon[[4,2,1]],af},kk,{bkgluon[[2]],x2},{"vq",vv},{bkgluon[[3]],nn}],
pg2}(* write the propagator in coordinate space *)
,
{pg1,partialD[lor1],ppgator[{bkgluon[[4,2,1]],af},kk,{bkgluon[[2]],x2},{"vq",vv},{bkgluon[[3]],nn}],
pg2}
]
;
condlist=condlist/.bkgluon[[4,3,1]]->GluonStrength[{lor1,lorindx,sunn}]
]
];
pglist=pglist/.propagator[pg_List,{xx_,0},{vv_,"o"},{n1_,n2_}]/;!FreeQ[pg,_partialD]:>propagator[pg,{xx,0},{vv,"o"},{n1,n2}];


{pglist,oplist,vlist,condlist}
]





End[]
