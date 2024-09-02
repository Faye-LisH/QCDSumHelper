(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)


ImExport::usage="ImExport[dir_,names_,{val_List,arg_List,types_List},{renval_,renarg_List,rens_List}] is an internal function QSSRHelper for Import and Export the results in Parallelozed evaluation."
ImExport::direrr="The directory does not exist!"

(*ToFile::usage=="Give a directory to Export and Import the results."*)


Begin["`Private`ImExport`"]



(*--- import & export ---*)
ImExport[dir_,names_,{val_List,arg_List,expr_List}]:=Block[{queue,tmpnames,files,results,filt,tmp,dirs,vals},	


If[!DirectoryQ[dir],Message[ImExport::direrr];Abort[]];
	
(*--------------------------------*)	
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
	{{Plus@@(1-filt),Length[filt]},tmp}
,

(*--- the experssions should be evaluated ---*)


(* the directory and file names for the expersssion which should be evaulated *)

	dirs=FileNameJoin[{dir,#<>".txt"}]&/@names;
	dirs=DeleteCases[filt dirs,0];
	
(* the experssion which should be evaluated *)	
	queue=DeleteCases[filt expr,0];


(* export the nonexist results(ParallelSubmit the experssions) *)
	DistributeDefinitions[#]&/@Join[val,arg];

	{{Plus@@(1-filt),Length[filt]},Table[ParallelSubmit[{queue,dirs,i},Export[dirs[[i]],FCI[queue[[i]]@@arg]]],{i,1,Length[queue]}]}

]
]


(* (*backup*)
(*--- import & export ---*)

ImExport[dir_,names_,{val_List,arg_List,types_List},{renval_,renarg_List,rens_List}]:=Block[{queue1,queue2,tmpnames,files,results,filt,tmp,fdir,vals,renvals,tmp2,tmp3},	

	
files=FileNames["*.txt",dir];
tmpnames=StringSplit[StringSplit[#,"\\"|"/"][[-1]],"."][[1]]&/@files;(* get the file names *)
tmpnames={tmpnames,files}//Transpose;

(* import the already exist results *)
tmp=(Plus@@Cases[tmpnames,{#,file_}\[RuleDelayed]ToExpression[Import[file]]])&/@names;
filt=Boole[#===0]&/@tmp;


(*----------------------------*)
If[(Plus@@filt)\[Equal]0,
(* if all results are already exist *)
       tmp
,

(*--- the experssions should be evaluated ---*)


(* the directory and file names for the expersssion which should be evaulated *)
(* names ~ Join[types,rens] *)
	queue1=FileNameJoin[{dir,#<>".txt"}]&/@names[[;;Length[types]]];
	queue2=FileNameJoin[{dir,#<>".txt"}]&/@names[[Length[types]+1;;]];


(* The ParallelSubmit is trick; here defines a directory function so that the directory can be distributed before submiting *)
	Do[fdir[types[[i]]]=queue1[[i]],{i,1,Length[queue1]}];
	Do[fdir[rens[[i]]]=queue2[[i]],{i,1,Length[queue2]}];
	vals=Append[val,fdir];
	renvals=Append[renval,fdir];


(* the experssion should be evaluated *)
	queue1=DeleteCases[filt[[;;Length[queue1]]]types,0];
	queue2=DeleteCases[filt[[Length[types]+1;;]]rens,0];


(* export the nonexist results(ParallelSubmit the experssions) *)
	DistributeDefinitions[#]&/@Append[Join[vals,renvals(*,queue1,queue2*)],fdir];



	tmp2=ParallelSubmit[Export[fdir[#],FCI[#@@arg]]]&/@queue1;
	tmp3=ParallelSubmit[Export[fdir[#],FCI[#@@renarg]]]&/@queue2;

(* Join[ParallelSubmit[...]&/@...,ParallelSubmit[...]&/@... doesn't work, it maps ParallelSubmit[] but doesn't use the definition of queue *)
	Join[tmp2,tmp3]

]
]*)


End[]
(*EndPackage[]*)
