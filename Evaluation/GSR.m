(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


GSR::usage = 
	"GSR[expr,momentum,{parameters__}] apply Gaussian transformation to expr"


GSR::argerr="Momentum `1` is not involved in the expression!"
GSR::aerr="Analytic intergation failed, for numerical integration, the `1`, `2`, and `3` must be numeric."


Begin["`Private`GSR`"]



Options[GSR]={
	Tol->10^(-12),
	Method->"Circle",
	Normalize->True,
	Sort->False,
	Numerical->"Auto"
}


GSR[expr_,q_,{ss_,tau_,s0_},OptionsPattern[]]:=Block[{tmp,log,null,q2,qlog,qlist,qdelta,qdeltalog,qlogn,qlistn,qdeltan,qdeltalogn,sort=OptionValue[Sort],gsrnorm=1,norm=OptionValue[Normalize],tol=OptionValue[Tol],num=OptionValue[Numerical]},
If[FreeQ[expr,Momentum[q,___]],
	Message[GSR::argerr,q];
	expr

,

	tmp=Expand[null expr+log[null] null//FCI];
	tmp=FeynAmpDenominatorSplit[tmp]/.{FeynAmpDenominator[PropagatorDenominator[pp_Momentum,0]]:>(1/Pair[pp,pp])};
	tmp=List@@(tmp)/.Pair[Momentum[q,___],Momentum[q,___]]->q2; (* a list of (q^2)^n Log[-q^2/u^2]^m <O> *)
	tmp=DeleteCases[tmp,aa_/;FreeQ[aa,_Log|Power[q2,_?Negative]]]/.Log->log;
	
	
	(* gather the terms with same q2^n log[-q2/u^2]^m *)
	tmp=Map[qlogsort,tmp];
	tmp=Gather[tmp,#1[[2;;]]==#2[[2;;]]&]/.null->1;
	tmp={#[[1]]&/@#,#[[1,2;;]]}&/@tmp;(* each term in tmp looks like {{fator1 <O1>, factor2 <O2>, ...}, {u2, nq, nlog}} *)
	
	(* qlog: q^(2n)Log[-q^2/u^2]^m ;
	qlist: q^(2n)Log[-q^2/u^2];
	qdelta: 1/q^(2n);
	qdeltalog: Log[-q^2/u^2]^m/q^(2n) 
	*)
	(* select the corresponding terms; unify the form *)
	qlog=Replace[Select[tmp,#[[2,2]]>=0&&#[[2,3]]>1&],{}->{{{},{}}},{0}];
	qlist=Replace[Select[tmp,#[[2,2]]>=0&&#[[2,3]]==1&],{}->{{{},{}}},{0}];
	qdelta=Replace[Select[tmp,#[[2,2]]<0&&#[[2,3]]==0&],{}->{{{},{}}},{0}];
	qdeltalog=Replace[Select[tmp,#[[2,2]]<0&&#[[2,3]]>0&],{}->{{{},{}}},{0}];

	(* whether use numerical integration *)
	If[(!((MatchQ[qlog,{{{},{}}}])&&(MatchQ[qdeltalog,{{{},{}}}])))||num===True,
	
		If[!(NumberQ[ss]&&NumberQ[tau]&&NumberQ[s0]),
			Message[GSR::aerr,tau,ss,s0];
			Abort[]
		,
			(*---------------------*)
			(* numerical integration *)
			If[norm=!=False,
			
				qlogn={Plus@@#[[1]],ngsrn[#[[2]],{tau,s0}]}&/@qlog;
				qlistn={Plus@@#[[1]],ngsrn[#[[2]],{tau,s0}]}&/@qlist;
				qdeltan={Plus@@#[[1]],gsrn[#[[2]],{tau,s0}]//N}&/@qdelta;(* analytic expression for ~ 1/q^2n *)
				qdeltalogn={Plus@@#[[1]],ngsrn[#[[2]],{tau,s0},tol]}&/@qdeltalog;
				
				(* sum up each terms *)
				qlogn=#[[1]]#[[2]]&/@qlogn;
				qlistn=#[[1]]#[[2]]&/@qlistn;
				qdeltan=#[[1]]#[[2]]&/@qdeltan;
				qdeltalogn=#[[1]]#[[2]]&/@qdeltalogn;
				
				(* total *)
				gsrnorm=Total[qlogn]+Total[qlistn]+Total[qdeltan]+Total[qdeltalogn]
			];
			
			If[norm=!="NormOnly",
				qlog={Plus@@#[[1]],ngsr[#[[2]],{ss,tau,s0}]}&/@qlog;
				qlist={Plus@@#[[1]],ngsr[#[[2]],{ss,tau,s0}]}&/@qlist;
				qdelta={Plus@@#[[1]],gsr[#[[2]],{ss,tau,s0}]//N}&/@qdelta;(* analytic expression for ~ 1/q^2n *)
				qdeltalog={Plus@@#[[1]],ngsr[#[[2]],{ss,tau,s0},tol]}&/@qdeltalog
			]
		]
	,
	(*----------------------*)
	(* analytical integration *)
		
		If[norm=!=False,
			qlogn={Plus@@#[[1]],gsrn[#[[2]],{tau,s0}]}&/@qlog;
			qlistn={Plus@@#[[1]],gsrn[#[[2]],{tau,s0}]}&/@qlist;
			qdeltan={Plus@@#[[1]],gsrn[#[[2]],{tau,s0}]}&/@qdelta;
			qdeltalogn={Plus@@#[[1]],gsrn[#[[2]],{tau,s0},tol]}&/@qdeltalog;
			
			(* combine the factor and integration *)
			qlogn=Apply[Times,qlogn,{1}];
			qlistn=Apply[Times,qlistn,{1}];
			qdeltan=Apply[Times,qdeltan,{1}];
			qdeltalogn=Apply[Times,qdeltalogn,{1}];
			
			(* total *)
			gsrnorm=Total[qlogn]+Total[qlistn]+Total[qdeltan]+Total[qdeltalogn]
		];
		
		If[norm=!="NormOnly",
			qlog={Plus@@#[[1]],gsr[#[[2]],{ss,tau,s0}]}&/@qlog;
			qlist={Plus@@#[[1]],gsr[#[[2]],{ss,tau,s0}]}&/@qlist;
			qdelta={Plus@@#[[1]],gsr[#[[2]],{ss,tau,s0}]}&/@qdelta;
			qdeltalog={Plus@@#[[1]],gsr[#[[2]],{ss,tau,s0},tol]}&/@qdeltalog
		]
	];

	If[norm==="NormOnly",
		gsrnorm
	,
		If[sort,
			(* keep the integral for each cases not combined *)
			{{qlog,qlist,qdelta,qdeltalog},gsrnorm}
		,
		
			(* combine the factor and integration *)
			qlog=Apply[Times,qlog,{1}];
			qlist=Apply[Times,qlist,{1}];
			qdelta=Apply[Times,qdelta,{1}];
			qdeltalog=Apply[Times,qdeltalog,{1}];
			
			(* total *)
			(Total[qlog]+Total[qlist]+Total[qdelta]+Total[qdeltalog])/gsrnorm
		
		]
	]
]
]




(* (q^2)^n Log[-q^2/u^2]^m c <O> -> {c <O>, u^2, n , m} ; Log[-q^2/(4Pi u^2)] may involved *)
qlogsort[expr_]:=Block[{tmp,factor,q2term,logterm,u2=0,nq=0,nlog=0},
factor=expr/._log->1/.q2->1;
tmp=(expr/factor);
q2term=tmp/._log->1;
logterm=tmp/q2term;


q2term/.Power[q2,nn_]:>(nq=nn;1)/.q2:>(nq=1;1);
logterm/.Power[_log,nn_]:>(nlog=nn;1)/._log:>(nlog=1;1);
logterm/.log[q2 uterm_]/;!FreeQ[uterm,_?Negative,{1}]:>(u2=-1/uterm;1);

{factor,u2,nq,nlog}
]


ngsr[{},{ss_,tau_,s0_}]=0;
ngsrn[{},{tau_,s0_}]=1;

(* numerical integration 1/Pi Im\[CapitalPi](s) exp[-(s-t)^2/(4 tau)] *)
(* qlog and qlist *)
ngsr[{u2_,nq_,nlog_},{ss_,tau_,s0_}]:=Block[{tmp,null,s,im,imu2,ims},
im=Expand[ (log[s/u2]-Pi I)^nlog];
(* seperate log[u^2], avoid it enter into NIntegrate *)
im=List@@(null+Expand[ (log[s]-log[u2]-Pi I)^nlog]);

(*im=im-(im/.Complex[_,_]->0);*)
im=Select[im,MatchQ[#,_ Complex[0,_]]&]/.null->1;
im=im/.Complex[0,ii_]:>ii;

imu2=im/.log[s]->1/.log->Log;(* log[u2] list *)
ims=(im/imu2)/.log->Log;(* log[s] list *)

tmp=NIntegrate[# s^nq Exp[-(s-ss)^2/(4 tau)],{s,0,s0}]&/@ims;
tmp=tmp imu2//Total;

1/Sqrt[4 Pi tau] 1/Pi tmp
]/;(nq>=0&&nlog>0) 

(* Integrate[ngsr, {ss, -Infinity, Infinity}] *)
ngsrn[{u2_,nq_,nlog_},{tau_,s0_}]:=Block[{tmp,null,s,im},
im=Expand[ (log[s/u2]-Pi I)^nlog];
(* seperate log[u^2], avoid it enter into NIntegrate *)
im=List@@(null+Expand[ (log[s]-log[u2]-Pi I)^nlog]);

(*im=im-(im/.Complex[_,_]->0);*)
im=Select[im,MatchQ[#,_ Complex[0,_]]&]/.null->1;
im=im/.Complex[0,ii_]:>ii;

imu2=im/.log[s]->1/.log->Log;(* log[u2] list *)
ims=(im/imu2)/.log->Log;(* log[s] list *)

tmp=NIntegrate[# s^nq ,{s,0,s0}]&/@ims;
tmp=tmp imu2//Total;

1/Pi tmp
]/;(nq>=0&&nlog>0) 



(* qdeltalog *)
(* 
(* compare with integral along the |z|=\[Epsilon] and along s\[PlusMinus]0i *)
(* log[q^2/u^2]/q2 *)
ss=1;
s00=2;
u22=2;
Table[{n,NIntegrate[-Exp[-(s-ss)^2/(4)]/s,{s,10^(-n),s00}]-Exp[-(ss-10^(-n))^2/4]Log[10^(-n)/u22]},{n,0,10}];
ListLinePlot[%]
%%[[-1,2]]

nqq=-1;
tmp=-1/(2Pi)NIntegrate[s00^(1+nqq)Exp[I(1+nqq)theta](Log[s00/u22]+I (theta-Pi))^1 Exp[-(s00 Exp[I theta]-ss)^2/(4 )],{theta,0,2Pi}]
%/.Complex[aa_,bb_/;bb<10^(-12)]:>aa 
*)
(* Im(log^m/s^n) is singular at s=0 ; integrate along |z|=s0 instead *)
ngsr[{u2_,nq_,nlog_},{ss_,tau_,s0_},tol_]:=Block[{tmp,theta},
tmp=-1/(2Pi)NIntegrate[s0^(1+nq)Exp[I(1+nq)theta](Log[s0/u2]+I (theta-Pi))^nlog Exp[-(s0 Exp[I theta]-ss)^2/(4 tau)],{theta,0,2Pi}];

Expand[tmp/Sqrt[4 Pi tau] ]/.Complex[aa_,bb_/;bb<tol]:>aa
]/;(nq<0&&nlog>0)


(* Integrate[ngsr, {ss, -Infinity, Infinity}] *)
ngsrn[{u2_,nq_,nlog_},{tau_,s0_},tol_]:=Block[{tmp,theta},
tmp=-1/(2Pi)NIntegrate[s0^(1+nq)Exp[I(1+nq)theta](Log[s0/u2]+I (theta-Pi))^nlog ,{theta,0,2Pi}];

Expand[tmp]/.Complex[aa_,bb_/;bb<tol]:>aa
]/;(nq<0&&nlog>0)



gsr[{},{ss_,tau_,s0_}]=0;
gsrn[{},{tau_,s0_}]=1;
(* integration 1/Pi Im\[CapitalPi](s) exp[-(s-t)^2/(4 tau)] *)
(* qlist *)
gsr[{u2_,nq_,1},{ss_,tau_,s0_}]:=Block[{tmp,s,sss,tauu,s00},
tmp=Integrate[-s^nq Exp[-(s-sss)^2/(4 tauu)],{s,0,s00},Assumptions->{sss>0,tauu>0,s00>0}];

1/Sqrt[4 Pi tau]  tmp/.{sss->ss,tauu->tau,s00->s0}
]/;nq>=0

(* Integrate[ngsr, {ss, -Infinity, Infinity}] *)
gsrn[{u2_,nq_,1},{tau_,s0_}]:=Block[{tmp,s,sss,tauu,s00},
-s0^(nq+1)/(nq+1)
]/;nq>=0

(* qdelta*)
gsr[{u2_,nq_,0},{ss_,tau_,s0_}]:=Block[{tmp,s,im,nnq,sss,tauu},
nnq=-nq;
(* \[Delta]^n(s) -> \[Delta](s)\[PartialD]_s^(n-1) *)
tmp=-1/((nnq-1)!)D[Exp[-(s-sss)^2/(4 tauu)],{s,nnq-1}]/.s->0;

1/Sqrt[4 Pi tau]  tmp/.{sss->ss,tauu->tau}
]/;nq<0

(* Integrate[ngsr, {ss, -Infinity, Infinity}] for \[Delta]^n(s)exp[-(s-ss)^2/(4 tau)] *)
gsrn[{u2_,-1,0},{tau_,s0_}]=-1;
gsrn[{u2_,nq_/;nq<-1,0},{tau_,s0_}]=0;



End[]
