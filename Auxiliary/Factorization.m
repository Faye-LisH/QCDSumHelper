(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



Factorization::usage = "Factorize the dimension-6, -7, and -8 condensate"


Begin["`Private`Factorization`"]
Options[Factorization] = {
	Method->"EMFirst"
}


(* factorization rules *)
rule={
 Condensate[{"Q81", f_}] :> (-2*Condensate["gg"]*Condensate[{f, f}]*quarkMass[f])/(CA*D), 
 Condensate[{"Q82", f_}] :> 0,
 Condensate[{"Q83", f_}] :> -1/2*(Condensate[{f, f}]*Condensate[{f, "G", f}]), 
 Condensate[{"Q85", f_}] :> (Condensate[{f, f}]*Condensate[{f, "G", f}])/4, 
 Condensate[{"Q86", f_}] :> -1/8*((-2 + CA^2)*(-2 + D)*Condensate[{f, f}]*Condensate[{f, "G", f}])/CA^2, 
 Condensate[{"A", f_}] -> 0, 
 Condensate[{"Q71", f_}] :> (Condensate["gg"]*Condensate[{f, f}])/(2*CA), 
 Condensate[{"Q72", f_}] -> 0,
 Condensate[{"Q73", f_}] -> 0, 
 Condensate[{"Q74", f_}] :> -((-1 + CA^2)*(-1 + D)*Condensate[{f, f}]^2*quarkMass[f])/(4*CA^2), 
 Condensate[{"Q6", f_}] :> -1/4*(CF*D*Condensate[{f, f}]^2)/CA, 
 Condensate[{"Q3", f_}] :> Condensate[{f, f}], 
 Condensate[{"Q5", f_}] :> Condensate[{f, "G", f}]
 };
 
 
 Q84EMrule=Condensate[{"Q84", f_}] :> (Condensate[{f, f}]*(-((4 + (-2 + CA^2)*D)*Condensate[{f, "G", f}]) + 
     2*(-1 + CA^2)*(-1 + D)*Condensate[{f, f}]*quarkMass[f]^2))/(8*CA^2);
 Q84VSrule=Condensate[{"Q84", f_}] :> ((-1 + CA^2)*Condensate[{f, f}]*(2*(-1 + D)*quarkMass[f]^2*Condensate[{f, f}] - D*Condensate[{f, "G", f}]))/(8*CA^2); 



Factorization[expr_,OptionsPattern[]]:=Block[{tmp},
tmp=expr/.rule;

If[OptionValue[Method]=="VSFirst",
	tmp/. Q84VSrule
,
	tmp/. Q84EMrule
]
]


End[]
