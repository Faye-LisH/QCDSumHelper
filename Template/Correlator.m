(* ::Package:: *)

(* ::Code::Initialization::Plain:: *)
(* Wolfram Language package *)

(* Author: ShungHong Li *)


(* ::Code::Initialization::Plain:: *)
Correlator::filerr="Not a proper file name or there is no file match the name."
Correlator::direrr="Can not read the results from this directory."
Correlator::gerr="Fail to get the required data from the file(s)."
Correlator::mwan="Momentum not match, momentum `1` in files will be replaced to `2`"


Begin["`Private`Correlator`"]


(* ::Code::Initialization::Plain:: *)
Options[Correlator]={
	HoldFlavor->True,
	NLO->True,
	Renormalization->True,
	Condensates->All,
	ShowasTable->"Auto",
	Subtract->"None"
}


(* ::Code::Initialization::Plain:: *)
Correlator[q_:"detect",{j1_,j1list_String},{j2_,j2list_String},dir_String,OptionsPattern[]]:=Block[{tmp,ja,jafactor,jalist,jb,jbfactor,jblist,filedir,files,filexist,la,lb,jablist,folders,type,momentum,tmpq},

(*{jalist,jblist}={j1,j2};*)

(* check the format of the file names *)

{jalist,jblist}=If[StringMatchQ[#,"*.txt"],
					#
				,
					If[StringFreeQ[#,"."],
						#<>".txt"
					,
						Message[Correlator::filerr];
						Abort[]
					]
				]&/@{j1list,j2list};

(*-------------------------------------*)
(* read the files in the directory *)
If[DirectoryQ[dir],

(* unify the format of directory *)
	filedir=StringReplace[dir,"\\"->"/"];
	If[StringTake[filedir,-1]=!="/",filedir=filedir<>"/"];

(* find the files match the given names in the directory *)
	{jalist,jblist}=(filedir<>#)&/@{jalist,jblist};

	filexist=FileExistsQ[#]&/@{jalist,jblist};
	filexist=And@@filexist;

	If[!filexist,
		Message[Correlator::direrr];
		Abort[]
	]
,
	Message[Correlator::direrr];
	Abort[]
];

(* read the txt files to a list of currents *)
{jalist,jblist}=(StringInsert[StringInsert[StringReplace[Import[#],"\n"->","],"{",1],"}",-1]//ToExpression)&/@{jalist,jblist};

(*------------------------------------------*)
(* get the directories of the folders store the results *)
files=FileNames["*",filedir];
files=files (Boole[DirectoryQ[#]]&/@files);
files=DeleteCases[files,0];
files=StringReplace[StringReplace[#,"\\"->"/"],"//"->"/"]&/@files;
(* unify the format of directory, each one is a directory of calculated OPE result. *)

(* check the folder sturcture *)
{la,lb}={Length[jalist],Length[jblist]};
jablist=Flatten[Table[ToString[i]<>"_"<>ToString[j],{i,la},{j,lb}]];
files=Sort[files];

If[Length[files]=!=(la lb),
	Message[Correlator::direrr];
	Abort[]
,
	jablist=Transpose[{files,jablist}];
	folders=StringMatchQ[StringTake[#[[1]],-StringLength[#[[2]]]],#[[2]]]&/@jablist;
	folders=And@@folders;
	If[!folders,
		Message[Correlator::direrr];
		Abort[]
	]
];

(*---------------------------------------------*)
(* compare with j1 & j2, get the factors of each current in j1 & j2 *)
{ja,jalist,jafactor}=findmatch[j1,jalist];
{jb,jblist,jbfactor}=findmatch[j2,jblist];

(*---------------------------------------------*)
(* load all results *)
(* it's impossible to check every files in each folder before load the data; stop when encounter problem *)

Check[
(* get the type of currents. the Current[...] has different length in FourquarkCurrent for fourquark-current and tetraquark-current. *)
	type=(Length[#[[1]]]<=3)&/@Join[jalist,jblist];
	type=And@@type;

	If[type,
		type=FourquarkOPE
		,
		type=TetraquarkOPE
	];

(* load the data *)
	tmp=Table[type[jalist[[i]],jblist[[j]],Parallelized->filedir<>ToString[i]<>"_"<>ToString[j],HoldFlavor->OptionValue[HoldFlavor],NLO->OptionValue[NLO],Renormalization->OptionValue[Renormalization],Condensates->OptionValue[Condensates]],{i,Length[ja]},{j,Length[jb]}]

,
	Message[Correlator::gerr];
	Abort[]
];

(*Print[tmp//FullForm];*)

(*---------------------------------------------*)
(* Gather the results *)

(* tmp is a Length[ja]xLength[jb] table, each term in tmp is the output of FourquarkOPE or TetraquarkOPE, which has the form: { {n, m}, { {list of results for a certain type condensate}, {....} } }  } *)
tmp=((Plus@@((SUNSimplify[Plus@@#]/.CF->(CA^2-1)/(2CA))&/@#[[2]]))&/@#)&/@tmp;

(* times the corresponding factors *)
tmp=tmp TensorProduct[jafactor,ComplexConjugate[jbfactor]];
tmp=Plus@@Flatten[tmp];

(* check the momentum in files and input *)
momentum=DeleteDuplicates[Cases[tmp,Pair[Momentum[mm_,___],Momentum[mm_,___]]:>mm,Infinity]];
If[Length[momentum]>1,
	Message[Correlator::gerr];
	Abort[]
,
	tmpq=momentum[[1]]
];

(*-------------------------*)
tmp=QGather[tmp,tmpq,ShowasTable->OptionValue[ShowasTable],Subtract->OptionValue[Subtract]];
(* replace the momentum in data if it's not match with the input momentum *)
If[q=!="detect"&&tmpq=!=q,
	Message[Correlator::mwan,tmpq,q];
	tmp/.tmpq->q
,
	tmp
]

]


findmatch[ji_,jf_]:=Block[{ja,jafactor,jalist,exist,jftmp,tmp},
ja=DotSimplify[Expand[ji//FCI]];
jftmp=DotSimplify[jf];(* seprate dot product of SU(N) and gamma-matirces *)

(* turn input current to a list *)
If[Head[ja]===Plus,
	ja=List@@ja
,
	ja={ja}
];

jafactor=ja/._FourquarkCurrent->1;
ja=ja/jafactor;


(* for each current in ji, find the matched one in jf, label it to 1 *)
ja=(Boole[QSimplify[#]===0]&/@(jftmp-#))&/@ja;



(* Gather the same currents *)
tmp=Transpose[{jafactor,ja}];

tmp=tmp//.{aa___,{fac1_,labels_List},bb___,{fac2_,labels_List},cc___}:>{aa,{fac1+fac2,labels},bb,cc};
tmp=Transpose[tmp];

jafactor=tmp[[1]];
ja=tmp[[2]];


(* if there exist current can't be find in jf *)
exist=!FreeQ[#,1]&/@ja;
If[!(And@@exist),
	Message[Correlator::gerr];
	Abort[]
];

(* turn the {0,...,1,...,0} label to integer label *)
ja=Position[#,1][[1,1]]&/@ja;
jalist=jf[[#]]&/@ja;

(* {label, the current, the factor} *)
{ja,jalist,jafactor}

]


(* ::Code::Initialization::Plain:: *)
End[]
(*EndPackage[]*)
