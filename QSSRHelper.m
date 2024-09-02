(* ::Package:: *)

(* Wolfram Language Package *)


(*
	This software is covered by the GNU General Public License 3.
	Copyright (C) 2021 Shuang-Hong Li
*)


(* Created by the Wolfram Workbench 2021/12/14 *)

BeginPackage["QSSRHelper`",{"FeynCalc`"}]

Print[Style["QSSRHelper","Text",Bold],Style["  is a package used for QCD Sum Rules calculation","Text"]];



(*quarkMass::usage = "quarkMass[a_] give the symbol of quark mass with flavor a."*)


$AC::usage = "Head of dummy Color index in adjoint representation which is automatically generated."
$ac::usage = "Head of dummy Color index in fundmental representation which is automatically generated."
(*-------------------------------------------------------------------------------------------*)


Begin["`Private`"]
(* Implementation of the package *)


(*-------------------------------------------------------------------------------------------*)


End[]


(*-------------------------------------------------------------------------------------------*)


(* Quote all symbols before load their definitions, so that the symbols in different files "know each other". *)
If[Global`$LoadTemplate===True,
	files=Flatten[{FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Evaluation"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Fourier"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Loop"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Style"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Lorentz"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Auxiliary"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Template"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Symbols"}]}];
,			
	files=Flatten[{FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Evaluation"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Fourier"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Loop"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Style"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Lorentz"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Auxiliary"}],
				FileNames[ "*.m",FileNameJoin@{DirectoryName[$InputFileName],"Symbols"}]}];			
];


symbols=StringSplit[StringSplit[#,"\\"|"/"][[-1]],"."][[1]]&/@files;
symbols=StringJoin["QSSRHelper`",#]&/@symbols;


options=Import[DirectoryName[$InputFileName]<>"OptionsList.m","Text"]//ToExpression;


ToExpression[#]&/@symbols;


(* Load the functions *)
Get/@files



Remove[files];
Remove[symbols];
Remove[options];
Remove[Global`$LoadTemplate];


EndPackage[]
