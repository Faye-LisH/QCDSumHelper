(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



D8Condensate::usage = "D8Condensate[lor1_,lor2_,lor3_,lor4_,lor5_,psi_] = \[LeftAngleBracket]\!\(\*SubsuperscriptBox[OverscriptBox[\(psi\), \(_\)], \(i\), \(a\)]\) (\!\(\*SuperscriptBox[\(\[Del]\), \(lor1\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor2\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor3\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor4\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor5\)]\)\!\(\*SubscriptBox[\(psi\), \(j\)]\)\!\(\*SuperscriptBox[\()\), \(b\)]\)\[RightAngleBracket]; "



Begin["`Private`D8Condensate`"]
Options[D8Condensate] = {
	Explicit->False,
	Massless->True,
	Factorization->"EMFirst"
}


(* <q_i DDDDDq_j > *)
D8Condensate[lor1:Except[_List],lor2:Except[_List],lor3:Except[_List],lor4_,lor5:Except[_Rule],f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=(((-1/4*I)*AGammaD[lor1, lor2, lor3, lor4, lor5]*Condensate[{"A", f}])/(CA*(-4 + D)*(-3 + D)*(-2 + D)*(-1 + D)*D) + 
 ((I/32)*(AGammaD[lor1, lor4, lor5]*MTD[lor2, lor3] + AGammaD[lor1, lor2, lor5]*MTD[lor3, lor4])*(-8*Condensate[{"A", f}] + 4*D*Condensate[{"A", f}] - 2*(-3 + D)*Condensate[{"Q81", f}] + 4*Condensate[{"Q82", f}] - 2*D*Condensate[{"Q82", f}] + 
    5*Condensate[{"Q83", f}] - D*Condensate[{"Q83", f}] + 10*Condensate[{"Q85", f}] - 2*D*Condensate[{"Q85", f}] - 12*Condensate[{"Q86", f}] + 4*D*Condensate[{"Q86", f}] - 
    2*(2*(-3 + D)*Condensate[{"Q71", f}] + (8 - 3*D)*Condensate[{"Q72", f}] - 4*(-4 + D)*Condensate[{"Q73", f}])*quarkMass[f] - 4*(-3 + D)*Condensate[{"Q6", f}]*quarkMass[f]^2 + 4*(-3 + D)*Condensate[{"Q5", f}]*quarkMass[f]^3))/
  (CA*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)) + ((I/32)*(AGammaD[lor2, lor4, lor5]*MTD[lor1, lor3] + AGammaD[lor1, lor2, lor4]*MTD[lor3, lor5])*(8*Condensate[{"A", f}] - 2*(-3 + D)*Condensate[{"Q81", f}] + 3*Condensate[{"Q83", f}] - 
    D*Condensate[{"Q83", f}] + 6*Condensate[{"Q85", f}] - 2*D*Condensate[{"Q85", f}] + 4*(-((-3 + D)*Condensate[{"Q71", f}]) + (-2 + D)*Condensate[{"Q72", f}] + 2*(-3 + D)*Condensate[{"Q73", f}])*quarkMass[f] - 
    4*(-3 + D)*Condensate[{"Q6", f}]*quarkMass[f]^2 + 4*(-3 + D)*Condensate[{"Q5", f}]*quarkMass[f]^3))/(CA*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)) + 
 ((I/32)*(AGammaD[lor3, lor4, lor5]*MTD[lor1, lor2] + AGammaD[lor1, lor2, lor3]*MTD[lor4, lor5])*(-4*D*Condensate[{"A", f}] - 2*(-3 + D)*Condensate[{"Q81", f}] - 4*Condensate[{"Q82", f}] + 2*D*Condensate[{"Q82", f}] + Condensate[{"Q83", f}] - 
    D*Condensate[{"Q83", f}] + 2*Condensate[{"Q85", f}] - 2*D*Condensate[{"Q85", f}] + 2*(-2*(-3 + D)*Condensate[{"Q71", f}] + (-6 + D)*Condensate[{"Q72", f}] + 4*(-2 + D)*Condensate[{"Q73", f}])*quarkMass[f] - 
    4*(-3 + D)*Condensate[{"Q6", f}]*quarkMass[f]^2 + 4*(-3 + D)*Condensate[{"Q5", f}]*quarkMass[f]^3))/(CA*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)) + 
 ((I/32)*AGammaD[lor2, lor3, lor4]*MTD[lor1, lor5]*(4*Condensate[{"A", f}] + 4*D*Condensate[{"A", f}] - 2*(-3 + D)*Condensate[{"Q81", f}] - 6*Condensate[{"Q82", f}] + 2*D*Condensate[{"Q82", f}] - 3*Condensate[{"Q83", f}] + D*Condensate[{"Q83", f}] + 
    6*Condensate[{"Q84", f}] - 2*D*Condensate[{"Q84", f}] - 6*Condensate[{"Q86", f}] + 2*D*Condensate[{"Q86", f}] - 2*(2*(-3 + D)*Condensate[{"Q71", f}] + (5 - 3*D)*Condensate[{"Q72", f}] + 2*(-3 + D)*Condensate[{"Q74", f}])*quarkMass[f] - 
    4*(-3 + D)*Condensate[{"Q6", f}]*quarkMass[f]^2 + 4*(-3 + D)*Condensate[{"Q5", f}]*quarkMass[f]^3))/(CA*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)) + 
 ((I/32)*(AGammaD[lor2, lor3, lor5]*MTD[lor1, lor4] + AGammaD[lor1, lor3, lor4]*MTD[lor2, lor5])*(-4*Condensate[{"A", f}] - 2*(-3 + D)*Condensate[{"Q81", f}] - 2*Condensate[{"Q82", f}] - Condensate[{"Q83", f}] + D*Condensate[{"Q83", f}] + 
    6*Condensate[{"Q84", f}] - 2*D*Condensate[{"Q84", f}] + 4*Condensate[{"Q85", f}] - 6*Condensate[{"Q86", f}] + 2*D*Condensate[{"Q86", f}] - 
    2*(2*(-3 + D)*Condensate[{"Q71", f}] + 7*Condensate[{"Q72", f}] - 2*D*Condensate[{"Q72", f}] + 4*Condensate[{"Q73", f}] - 6*Condensate[{"Q74", f}] + 2*D*Condensate[{"Q74", f}])*quarkMass[f] - 4*(-3 + D)*Condensate[{"Q6", f}]*quarkMass[f]^2 + 
    4*(-3 + D)*Condensate[{"Q5", f}]*quarkMass[f]^3))/(CA*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)) + 
 ((I/32)*AGammaD[lor1, lor3, lor5]*MTD[lor2, lor4]*(-20*Condensate[{"A", f}] + 4*D*Condensate[{"A", f}] - 2*(-3 + D)*Condensate[{"Q81", f}] + 10*Condensate[{"Q82", f}] - 2*D*Condensate[{"Q82", f}] - Condensate[{"Q83", f}] - D*Condensate[{"Q83", f}] + 
    6*Condensate[{"Q84", f}] - 2*D*Condensate[{"Q84", f}] + 4*Condensate[{"Q85", f}] - 4*D*Condensate[{"Q85", f}] - 18*Condensate[{"Q86", f}] + 6*D*Condensate[{"Q86", f}] - 
    2*(2*(-3 + D)*Condensate[{"Q71", f}] + 11*Condensate[{"Q72", f}] - 3*D*Condensate[{"Q72", f}] + 4*Condensate[{"Q73", f}] - 4*D*Condensate[{"Q73", f}] - 6*Condensate[{"Q74", f}] + 2*D*Condensate[{"Q74", f}])*quarkMass[f] - 
    4*(-3 + D)*Condensate[{"Q6", f}]*quarkMass[f]^2 + 4*(-3 + D)*Condensate[{"Q5", f}]*quarkMass[f]^3))/(CA*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)) - 
 ((I/16)*(GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5])*(8*Condensate[{"A", f}] - 8*D*Condensate[{"A", f}] - 2*(4 + D)*Condensate[{"Q81", f}] - 4*Condensate[{"Q82", f}] + 4*D*Condensate[{"Q82", f}] - 
    2*Condensate[{"Q83", f}] - D*Condensate[{"Q83", f}] - 4*Condensate[{"Q84", f}] - 16*Condensate[{"Q85", f}] + 2*D*Condensate[{"Q85", f}] + 16*Condensate[{"Q86", f}] - 6*D*Condensate[{"Q86", f}] + 
    2*(4 - 10*D + D^2)*Condensate[{"Q71", f}]*quarkMass[f] - 2*(-1 + D)*(-4*Condensate[{"Q73", f}] + D*(Condensate[{"Q72", f}] + 2*Condensate[{"Q73", f}]) + 4*Condensate[{"Q74", f}])*quarkMass[f] + 20*(-2 + D)*Condensate[{"Q6", f}]*quarkMass[f]^2 - 
    4*(2 - 3*D + D^2)*Condensate[{"Q5", f}]*quarkMass[f]^3 + 4*(2 - 3*D + D^2)*Condensate[{"Q3", f}]*quarkMass[f]^5))/(CA*(-2 + D)*(-1 + D)*D*(2 + D)*(4 + D)) - 
 ((I/16)*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*(-8*Condensate[{"A", f}] + 4*D*Condensate[{"A", f}] + 4*D^2*Condensate[{"A", f}] + 2*(-4 + 3*D + D^2)*Condensate[{"Q81", f}] + 4*Condensate[{"Q82", f}] - 2*D*Condensate[{"Q82", f}] - 
    2*D^2*Condensate[{"Q82", f}] - 6*Condensate[{"Q83", f}] + 2*D*Condensate[{"Q83", f}] + D^2*Condensate[{"Q83", f}] + 4*Condensate[{"Q84", f}] + 2*D*Condensate[{"Q84", f}] + 6*D*Condensate[{"Q85", f}] - 8*Condensate[{"Q86", f}] - 
    4*D*Condensate[{"Q86", f}] + 2*D^2*Condensate[{"Q86", f}] + 2*(20 - 6*D + D^2)*Condensate[{"Q71", f}]*quarkMass[f] - 4*(-2*(-1 + D)*Condensate[{"Q72", f}] + (2 - 3*D + D^2)*Condensate[{"Q73", f}] + (-6 + D)*Condensate[{"Q74", f}])*quarkMass[f] + 
    20*(-2 + D)*Condensate[{"Q6", f}]*quarkMass[f]^2 - 4*(2 - 3*D + D^2)*Condensate[{"Q5", f}]*quarkMass[f]^3 + 4*(2 - 3*D + D^2)*Condensate[{"Q3", f}]*quarkMass[f]^5))/(CA*(-2 + D)*(-1 + D)*D*(2 + D)*(4 + D)) - 
 ((I/16)*(GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5])*
   (D*(-2*Condensate[{"A", f}] + Condensate[{"Q82", f}] - Condensate[{"Q83", f}] + Condensate[{"Q84", f}] - Condensate[{"Q85", f}] - 3*Condensate[{"Q86", f}]) - 
    D^2*(-2*Condensate[{"A", f}] + Condensate[{"Q82", f}] + 2*Condensate[{"Q85", f}] - 2*Condensate[{"Q86", f}]) - 2*(Condensate[{"Q83", f}] - 2*Condensate[{"Q85", f}] + 2*Condensate[{"Q86", f}]) + 
    (-1 + D)*(4*(-2 + D)*Condensate[{"Q71", f}] + 4*Condensate[{"Q72", f}] - D*Condensate[{"Q72", f}] + 8*Condensate[{"Q73", f}] - 4*D*Condensate[{"Q73", f}] + 2*D*Condensate[{"Q74", f}])*quarkMass[f] + 
    2*(4 + 4*D - 3*D^2)*Condensate[{"Q6", f}]*quarkMass[f]^2 - 4*(2 - 3*D + D^2)*Condensate[{"Q5", f}]*quarkMass[f]^3 + 4*(2 - 3*D + D^2)*Condensate[{"Q3", f}]*quarkMass[f]^5))/(CA*(-2 + D)*(-1 + D)*D*(2 + D)*(4 + D)) - 
 ((I/16)*MTD[lor2, lor4]*(GAD[lor5]*MTD[lor1, lor3] + GAD[lor1]*MTD[lor3, lor5])*(16*Condensate[{"A", f}] - 6*D*Condensate[{"A", f}] + (-8 + 3*D)*Condensate[{"Q82", f}] + 2*Condensate[{"Q83", f}] + D*Condensate[{"Q84", f}] - 4*Condensate[{"Q85", f}] + 
    5*D*Condensate[{"Q85", f}] + 12*Condensate[{"Q86", f}] - 7*D*Condensate[{"Q86", f}] + (4*(2 - 3*D + D^2)*Condensate[{"Q71", f}] + (4 + 3*D - 2*D^2)*Condensate[{"Q72", f}] - 2*(-1 + D)*(2*(-2 + D)*Condensate[{"Q73", f}] - D*Condensate[{"Q74", f}]))*
     quarkMass[f] - 4*(2 - 3*D + D^2)*Condensate[{"Q6", f}]*quarkMass[f]^2 - 4*(2 - 3*D + D^2)*Condensate[{"Q5", f}]*quarkMass[f]^3 + 4*(2 - 3*D + D^2)*Condensate[{"Q3", f}]*quarkMass[f]^5))/(CA*(-2 + D)*(-1 + D)*D*(2 + D)*(4 + D)) - 
 ((I/16)*(GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4])*(16*Condensate[{"A", f}] - 6*D*Condensate[{"A", f}] + 2*(4 + D)*Condensate[{"Q81", f}] - 8*Condensate[{"Q82", f}] + 3*D*Condensate[{"Q82", f}] + 
    2*Condensate[{"Q83", f}] + D*Condensate[{"Q84", f}] - 4*Condensate[{"Q85", f}] + 5*D*Condensate[{"Q85", f}] + 12*Condensate[{"Q86", f}] - 7*D*Condensate[{"Q86", f}] + 
    ((8 - 4*D + 6*D^2)*Condensate[{"Q71", f}] + (4 + 3*D - 2*D^2)*Condensate[{"Q72", f}] - 2*(-1 + D)*(2*(-2 + D)*Condensate[{"Q73", f}] - D*Condensate[{"Q74", f}]))*quarkMass[f] - 4*(2 - 3*D + D^2)*Condensate[{"Q6", f}]*quarkMass[f]^2 - 
    4*(2 - 3*D + D^2)*Condensate[{"Q5", f}]*quarkMass[f]^3 + 4*(2 - 3*D + D^2)*Condensate[{"Q3", f}]*quarkMass[f]^5))/(CA*(-2 + D)*(-1 + D)*D*(2 + D)*(4 + D)) - 
 ((I/32)*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*(D*(-8*Condensate[{"A", f}] - 8*Condensate[{"Q81", f}] + 4*Condensate[{"Q82", f}] + 2*Condensate[{"Q83", f}] - 4*Condensate[{"Q84", f}]) + 
    4*(-4*Condensate[{"A", f}] + 2*Condensate[{"Q82", f}] - Condensate[{"Q83", f}] + 2*Condensate[{"Q84", f}] + 4*Condensate[{"Q85", f}] - 4*Condensate[{"Q86", f}]) + 
    D^2*(4*Condensate[{"A", f}] - 2*Condensate[{"Q81", f}] - 2*Condensate[{"Q82", f}] + Condensate[{"Q83", f}] - 2*Condensate[{"Q84", f}] - 4*Condensate[{"Q85", f}] + 6*Condensate[{"Q86", f}]) + 
    2*(4*(-2 - 4*D + D^2)*Condensate[{"Q71", f}] - (8 - 4*D + D^2)*Condensate[{"Q72", f}] - 4*(2 - 3*D + D^2)*Condensate[{"Q73", f}] + 2*(-4 - 2*D + D^2)*Condensate[{"Q74", f}])*quarkMass[f] + 8*(6 + D - 2*D^2)*Condensate[{"Q6", f}]*quarkMass[f]^2 - 
    8*(2 - 3*D + D^2)*Condensate[{"Q5", f}]*quarkMass[f]^3 + 8*(2 - 3*D + D^2)*Condensate[{"Q3", f}]*quarkMass[f]^5))/(CA*(-2 + D)*(-1 + D)*D*(2 + D)*(4 + D)) - 
 ((I/32)*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*(-16*Condensate[{"A", f}] - 8*D*Condensate[{"A", f}] + 4*D^2*Condensate[{"A", f}] - 2*(-8 + 2*D + D^2)*Condensate[{"Q81", f}] + 8*Condensate[{"Q82", f}] + 4*D*Condensate[{"Q82", f}] - 
    2*D^2*Condensate[{"Q82", f}] + 4*Condensate[{"Q83", f}] - 4*D*Condensate[{"Q83", f}] - D^2*Condensate[{"Q83", f}] + 8*Condensate[{"Q84", f}] - 4*D*Condensate[{"Q84", f}] - 2*D^2*Condensate[{"Q84", f}] - 4*D*Condensate[{"Q85", f}] - 
    4*D^2*Condensate[{"Q85", f}] - 16*Condensate[{"Q86", f}] + 6*D^2*Condensate[{"Q86", f}] + 2*((-8 - 8*D + 6*D^2)*Condensate[{"Q71", f}] - (8 - 4*D + D^2)*Condensate[{"Q72", f}] - 4*(2 - 3*D + D^2)*Condensate[{"Q73", f}] + 
      2*(-4 - 2*D + D^2)*Condensate[{"Q74", f}])*quarkMass[f] + 8*(6 + D - 2*D^2)*Condensate[{"Q6", f}]*quarkMass[f]^2 - 8*(2 - 3*D + D^2)*Condensate[{"Q5", f}]*quarkMass[f]^3 + 8*(2 - 3*D + D^2)*Condensate[{"Q3", f}]*quarkMass[f]^5))/
  (CA*(-2 + D)*(-1 + D)*D*(2 + D)*(4 + D)) - ((I/32)*(GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5])*(-16*Condensate[{"A", f}] - 8*D*Condensate[{"A", f}] + 4*D^2*Condensate[{"A", f}] + 
    2*(-8 + 2*D + D^2)*Condensate[{"Q81", f}] + 8*Condensate[{"Q82", f}] + 4*D*Condensate[{"Q82", f}] - 2*D^2*Condensate[{"Q82", f}] - 4*Condensate[{"Q83", f}] + 2*D*Condensate[{"Q83", f}] + D^2*Condensate[{"Q83", f}] - 8*Condensate[{"Q84", f}] + 
    4*D*Condensate[{"Q85", f}] - 2*D^2*Condensate[{"Q85", f}] - 4*D*Condensate[{"Q86", f}] + 4*D^2*Condensate[{"Q86", f}] + 4*(12 - 8*D + D^2)*Condensate[{"Q71", f}]*quarkMass[f] - 2*(8 - 4*D + D^2)*Condensate[{"Q72", f}]*quarkMass[f] - 
    8*(-1 + D)*((-2 + D)*Condensate[{"Q73", f}] + 2*Condensate[{"Q74", f}])*quarkMass[f] - 4*(12 - 8*D + D^2)*Condensate[{"Q6", f}]*quarkMass[f]^2 - 8*(2 - 3*D + D^2)*Condensate[{"Q5", f}]*quarkMass[f]^3 + 
    8*(2 - 3*D + D^2)*Condensate[{"Q3", f}]*quarkMass[f]^5))/(CA*(-2 + D)*(-1 + D)*D*(2 + D)*(4 + D)) - 
 ((I/32)*MTD[lor1, lor5]*(GAD[lor4]*MTD[lor2, lor3] + GAD[lor2]*MTD[lor3, lor4])*(-4*D*Condensate[{"A", f}] + 4*D^2*Condensate[{"A", f}] - 2*(-8 + 2*D + D^2)*Condensate[{"Q81", f}] + 2*D*Condensate[{"Q82", f}] - 2*D^2*Condensate[{"Q82", f}] + 
    12*Condensate[{"Q83", f}] - 2*D*Condensate[{"Q83", f}] - D^2*Condensate[{"Q83", f}] + 2*D*Condensate[{"Q84", f}] + 8*Condensate[{"Q85", f}] + 6*D*Condensate[{"Q85", f}] - 2*D^2*Condensate[{"Q85", f}] - 8*Condensate[{"Q86", f}] - 
    6*D*Condensate[{"Q86", f}] + 4*D^2*Condensate[{"Q86", f}] + 4*(-4 - 4*D + 3*D^2)*Condensate[{"Q71", f}]*quarkMass[f] - 2*(-1 + D)*((-4 + D)*Condensate[{"Q72", f}] + 4*(-2 + D)*Condensate[{"Q73", f}] - 2*D*Condensate[{"Q74", f}])*quarkMass[f] - 
    4*(-4 - 4*D + 3*D^2)*Condensate[{"Q6", f}]*quarkMass[f]^2 - 8*(2 - 3*D + D^2)*Condensate[{"Q5", f}]*quarkMass[f]^3 + 8*(2 - 3*D + D^2)*Condensate[{"Q3", f}]*quarkMass[f]^5))/(CA*(-2 + D)*(-1 + D)*D*(2 + D)*(4 + D)));
    
    
If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* <q_i (DG)^m G^n q_j> *)
D8Condensate[{colm_,{lor3_,{lor4_,lor5_}}},{coln_,{lor1_,lor2_}},f:Except[_Rule],ops___Rule]:=D8Condensate[{coln,{lor1,lor2}},{colm,{lor3,{lor4,lor5}}},f,ops]

(* <q_i G^n (DG)^m q_j> *)
D8Condensate[{coln_,{lor1_,lor2_}},{colm_,{lor3_,{lor4_,lor5_}}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=((I/4)*(4*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] - 4*CA^2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*
    MTD[lor1, lor5] + 6*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q86", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] - 
   2*CA^2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q86", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] + 4*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*
    MTD[lor1, lor5] - 4*CA^2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5] + 
   6*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q86", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5] - 2*CA^2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q86", f}]*SUNT[coln] . SUNT[colm]*
    MTD[lor1, lor5] + 8*CA^2*AGammaD[lor1, lor4, lor5]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor3] - 
   8*CA^2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor3] + 12*CA^2*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*SUNT[colm] . SUNT[coln]*
    MTD[lor2, lor3] - 4*CA^2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor3] + 
   8*CA^2*AGammaD[lor1, lor4, lor5]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor3] - 8*CA^2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*
    MTD[lor2, lor3] + 12*CA^2*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor3] - 
   4*CA^2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor3] - 6*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor4]*
    MTD[lor2, lor3] + 3*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
   8*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 4*CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor4]*
    MTD[lor2, lor3] - 2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
   CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 12*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor4]*
    MTD[lor2, lor3] - 6*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
   4*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor4]*
    MTD[lor2, lor3] + 6*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 3*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*
    MTD[lor1, lor4]*MTD[lor2, lor3] - 8*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
   4*CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor4]*
    MTD[lor2, lor3] - CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
   12*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 6*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor4]*
    MTD[lor2, lor3] + 4*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
   2*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 6*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor5]*
    MTD[lor2, lor3] - 3*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
   8*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 4*CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor5]*
    MTD[lor2, lor3] + 2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
   CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 12*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor5]*
    MTD[lor2, lor3] + 6*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
   4*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor5]*
    MTD[lor2, lor3] - 6*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 3*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*
    MTD[lor1, lor5]*MTD[lor2, lor3] + 8*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
   4*CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor5]*
    MTD[lor2, lor3] + CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
   12*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 6*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor5]*
    MTD[lor2, lor3] - 4*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
   2*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 4*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*
    MTD[lor2, lor4] - 4*CA^2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4] + 
   6*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4] - 2*CA^2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*SUNT[colm] . SUNT[coln]*
    MTD[lor2, lor4] + 4*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] - 
   4*CA^2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] + 6*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*SUNT[coln] . SUNT[colm]*
    MTD[lor2, lor4] - 2*CA^2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] + 
   6*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 3*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*
    MTD[lor2, lor4] - 8*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
   4*CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*
    MTD[lor2, lor4] - CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
   12*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 6*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*
    MTD[lor2, lor4] + 4*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
   2*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 6*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*
    MTD[lor2, lor4] + 3*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
   8*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 4*CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*
    MTD[lor2, lor4] - 2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
   CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 12*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*
    MTD[lor2, lor4] - 6*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
   4*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*
    MTD[lor2, lor4] + 12*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 6*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*
    MTD[lor1, lor5]*MTD[lor2, lor4] - 16*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
   8*CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 4*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*
    MTD[lor2, lor4] - 2*CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
   24*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 12*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*
    MTD[lor2, lor4] + 8*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
   4*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 12*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*
    MTD[lor2, lor4] + 6*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
   16*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 8*CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*
    MTD[lor2, lor4] - 4*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
   2*CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 24*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*
    MTD[lor2, lor4] - 12*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
   8*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 4*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*
    MTD[lor2, lor4] - 4*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] + 
   4*CA^2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] - 6*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*SUNT[colm] . SUNT[coln]*
    MTD[lor2, lor5] + 2*CA^2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] - 
   4*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] + 4*CA^2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*
    MTD[lor2, lor5] - 6*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] + 
   2*CA^2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] - 6*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*
    MTD[lor2, lor5] + 3*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
   8*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 4*CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*
    MTD[lor2, lor5] - 2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
   CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 12*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*
    MTD[lor2, lor5] - 6*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
   4*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*
    MTD[lor2, lor5] + 6*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 3*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*
    MTD[lor1, lor3]*MTD[lor2, lor5] - 8*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
   4*CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*
    MTD[lor2, lor5] - CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
   12*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 6*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*
    MTD[lor2, lor5] + 4*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
   2*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 12*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*
    MTD[lor2, lor5] + 6*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
   16*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 8*CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*
    MTD[lor2, lor5] - 4*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
   2*CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 24*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*
    MTD[lor2, lor5] - 12*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
   8*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 4*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*
    MTD[lor2, lor5] + 12*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 6*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*
    MTD[lor1, lor4]*MTD[lor2, lor5] - 16*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
   8*CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 4*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*
    MTD[lor2, lor5] - 2*CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
   24*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 12*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*
    MTD[lor2, lor5] + 8*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
   4*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 24*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*
    MTD[lor3, lor4] + 12*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q86", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4] - 
   4*CA^2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q86", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4] + 24*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*
    MTD[lor3, lor4] + 12*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q86", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4] - 
   4*CA^2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q86", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4] + 18*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*
    MTD[lor3, lor4] - 9*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
   6*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 3*CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*
    MTD[lor3, lor4] + 12*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 6*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*
    MTD[lor1, lor5]*MTD[lor3, lor4] + 8*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
   4*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 4*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*
    MTD[lor3, lor4] + 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
   18*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 9*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*
    MTD[lor3, lor4] + 6*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
   3*CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 12*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*
    MTD[lor3, lor4] + 6*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
   8*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 4*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*
    MTD[lor3, lor4] + 4*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
   2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 18*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*
    MTD[lor3, lor4] + 9*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
   6*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 3*CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*
    MTD[lor3, lor4] - 12*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 6*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*
    MTD[lor2, lor5]*MTD[lor3, lor4] - 8*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
   4*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 4*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*
    MTD[lor3, lor4] - 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
   18*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 9*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*
    MTD[lor3, lor4] - 6*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
   3*CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 12*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*
    MTD[lor3, lor4] - 6*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
   8*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 4*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*
    MTD[lor3, lor4] - 4*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
   2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 24*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*
    MTD[lor3, lor5] - 12*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5] + 
   4*CA^2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5] - 24*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*
    MTD[lor3, lor5] - 12*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5] + 
   4*CA^2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5] - 18*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*
    MTD[lor3, lor5] + 9*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
   6*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 3*CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*
    MTD[lor3, lor5] - 12*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 6*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*
    MTD[lor1, lor4]*MTD[lor3, lor5] - 8*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
   4*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 4*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*
    MTD[lor3, lor5] - 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
   18*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 9*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*
    MTD[lor3, lor5] - 6*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
   3*CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 12*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*
    MTD[lor3, lor5] - 6*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
   8*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 4*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*
    MTD[lor3, lor5] - 4*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
   2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 18*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*
    MTD[lor3, lor5] - 9*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
   6*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 3*CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*
    MTD[lor3, lor5] + 12*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 6*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*
    MTD[lor2, lor4]*MTD[lor3, lor5] + 8*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
   4*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 4*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*
    MTD[lor3, lor5] + 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
   18*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 9*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*
    MTD[lor3, lor5] + 6*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
   3*CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 12*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*
    MTD[lor3, lor5] + 6*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
   8*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 4*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*
    MTD[lor3, lor5] + 4*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
   2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 2*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*
    MTD[lor1, lor5]*quarkMass[f] - 2*CA^2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5]*quarkMass[f] + 
   2*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5]*quarkMass[f] - 2*CA^2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*
    SUNT[coln] . SUNT[colm]*MTD[lor1, lor5]*quarkMass[f] + 4*CA^2*AGammaD[lor1, lor4, lor5]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor3]*quarkMass[f] - 
   4*CA^2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor3]*quarkMass[f] + 
   4*CA^2*AGammaD[lor1, lor4, lor5]*Condensate[{"Q72", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor3]*quarkMass[f] - 4*CA^2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q72", f}]*
    SUNT[coln] . SUNT[colm]*MTD[lor2, lor3]*quarkMass[f] + 2*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4]*quarkMass[f] - 
   2*CA^2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4]*quarkMass[f] + 
   2*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4]*quarkMass[f] - 2*CA^2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*
    SUNT[coln] . SUNT[colm]*MTD[lor2, lor4]*quarkMass[f] - 2*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5]*quarkMass[f] + 
   2*CA^2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5]*quarkMass[f] - 
   2*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5]*quarkMass[f] + 2*CA^2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*
    SUNT[coln] . SUNT[colm]*MTD[lor2, lor5]*quarkMass[f] + 12*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4]*quarkMass[f] + 
   12*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q72", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4]*quarkMass[f] - 12*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q72", f}]*
    SUNT[colm] . SUNT[coln]*MTD[lor3, lor5]*quarkMass[f] - 12*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q72", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5]*quarkMass[f] + 
   4*CA^2*AGammaD[lor2, lor4, lor5]*(SUNT[colm] . SUNT[coln] + SUNT[coln] . SUNT[colm])*MTD[lor1, lor3]*(2*(-1 + D)*Condensate[{"A", f}] + (-3 + D)*Condensate[{"Q86", f}] + 
     (-1 + D)*Condensate[{"Q72", f}]*quarkMass[f]) + 2*CA^2*AGammaD[lor2, lor3, lor5]*(SUNT[colm] . SUNT[coln] + SUNT[coln] . SUNT[colm])*MTD[lor1, lor4]*
    (2*(-1 + D)*Condensate[{"A", f}] + (-3 + D)*Condensate[{"Q86", f}] + (-1 + D)*Condensate[{"Q72", f}]*quarkMass[f])))/(CA*(2 - 3*CA^2 + CA^4)*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D));
    
    
If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* <q_i G^n G^m Dq_j> *)
(*D8Condensate[{coln_,{lor1_,lor2:Except[_List]}},{colm_,{lor3_,lor4:Except[_List]}},lor5_,f:Except[_Rule],ops___Rule]:=D8Condensate[lor5,{coln,{lor3,lor4}},{colm,{lor1,lor2}},f,ops]/.{ga_GAD:>-ga,AGammaD[lors__]/;Length[{lors}]==5:>-AGammaD[lors]}
*)
D8Condensate[{coln_,{lor1_,lor2:Except[_List]}},{colm_,{lor3_,lor4:Except[_List]}},lor5_,f:Except[_Rule],OptionsPattern[]]:=Block[{cond},
cond=D8Condensate[lor5,{coln,{lor3,lor4}},{colm,{lor1,lor2}},f,Explicit->False,Massless->OptionValue[Massless],Factorization->OptionValue[Factorization]]/.{ga_GAD:>-ga,AGammaD[lors__]/;Length[{lors}]==5:>-AGammaD[lors]};

(*If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];*)
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]




(* <q_i(<-D) G^n G^m q_j> *)
D8Condensate[lor1_,{coln_,{lor2_,lor3:Except[_List]}},{colm_,{lor4_,lor5:Except[_List]}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=(((I/8)*(-32*CA^2*(2 + D)*AGammaD[lor1, lor2, lor3, lor4, lor5]*Condensate[{"A", f}]*(SUNT[colm] . SUNT[coln] + SUNT[coln] . SUNT[colm]) + 
   (-4 + D)*(-16*CA^2*AGammaD[lor2, lor3, lor5]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4] - 8*CA^2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4] - 
     16*AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4] + 8*CA^2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4] + 
     8*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4] - 4*CA^2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4] - 
     8*AGammaD[lor2, lor3, lor5]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4] + 4*CA^2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4] - 
     16*AGammaD[lor2, lor3, lor5]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4] + 8*CA^2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4] - 
     16*CA^2*AGammaD[lor2, lor3, lor5]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4] - 8*CA^2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4] + 
     16*AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4] - 8*CA^2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4] - 
     8*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4] + 4*CA^2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4] + 
     8*AGammaD[lor2, lor3, lor5]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4] - 4*CA^2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4] + 
     16*AGammaD[lor2, lor3, lor5]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4] - 8*CA^2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4] + 
     16*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] + 8*CA^2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] + 
     16*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] - 8*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] - 
     8*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] + 4*CA^2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] + 
     8*AGammaD[lor2, lor3, lor4]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] - 4*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] + 
     16*AGammaD[lor2, lor3, lor4]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] - 8*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5] + 
     16*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5] + 8*CA^2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5] - 
     16*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5] + 8*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5] + 
     8*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5] - 4*CA^2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5] - 
     8*AGammaD[lor2, lor3, lor4]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5] + 4*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5] - 
     16*AGammaD[lor2, lor3, lor4]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5] + 8*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5] - 
     16*AGammaD[lor1, lor3, lor5]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4] + 8*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4] + 
     4*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4] - 2*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4] + 
     4*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4] - 2*CA^2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4] + 
     8*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4] - 4*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4] + 
     8*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4] - 4*CA^2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4] + 
     16*AGammaD[lor1, lor3, lor5]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] - 8*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] - 
     4*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] + 2*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] - 
     4*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] + 2*CA^2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] - 
     8*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] + 4*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] - 
     8*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] + 4*CA^2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4] + 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
     12*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 6*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
     2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
     2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
     24*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 12*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
     4*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
     4*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
     12*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 6*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
     2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
     2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
     24*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 12*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
     4*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
     4*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
     12*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 6*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
     2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
     2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
     24*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 12*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
     4*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
     4*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
     12*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 6*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
     2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
     2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
     24*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 12*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
     4*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
     4*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
     16*AGammaD[lor1, lor3, lor4]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] - 8*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] - 
     4*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] + 2*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] - 
     4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] + 2*CA^2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] - 
     8*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] + 4*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] - 
     8*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] + 4*CA^2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5] - 
     16*AGammaD[lor1, lor3, lor4]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] + 8*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] + 
     4*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] - 2*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] + 
     4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] - 2*CA^2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] + 
     8*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] - 4*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] + 
     8*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] - 4*CA^2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5] - 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
     12*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 6*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
     2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
     2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
     24*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 12*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
     4*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
     4*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
     12*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 6*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
     2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
     2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
     24*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 12*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
     4*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
     4*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
     12*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 6*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
     2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
     2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
     24*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 12*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
     4*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
     4*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
     12*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 6*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
     2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
     2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
     24*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 12*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
     4*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
     4*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
     16*AGammaD[lor1, lor2, lor5]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4] - 8*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4] - 
     4*AGammaD[lor1, lor2, lor5]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4] + 2*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4] - 
     4*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4] + 2*CA^2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4] - 
     8*AGammaD[lor1, lor2, lor5]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4] + 4*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4] - 
     8*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4] + 4*CA^2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4] - 
     16*AGammaD[lor1, lor2, lor5]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4] + 8*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4] + 
     4*AGammaD[lor1, lor2, lor5]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4] - 2*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4] + 
     4*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4] - 2*CA^2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4] + 
     8*AGammaD[lor1, lor2, lor5]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4] - 4*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4] + 
     8*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4] - 4*CA^2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4] - 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
     12*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 6*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
     2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
     2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
     24*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 12*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
     4*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
     4*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
     12*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 6*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
     2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
     2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
     24*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 12*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
     4*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
     4*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
     12*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 6*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
     2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
     2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
     24*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 12*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
     4*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
     4*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
     12*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 6*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
     2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
     2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
     24*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 12*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
     4*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
     4*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
     24*CA^2*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 8*CA^2*D*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
     24*CA^2*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 8*CA^2*D*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 
     16*AGammaD[lor1, lor2, lor4]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5] + 8*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q82", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5] + 
     4*AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5] - 2*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5] + 
     4*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5] - 2*CA^2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5] + 
     8*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5] - 4*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5] + 
     8*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5] - 4*CA^2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5] + 
     16*AGammaD[lor1, lor2, lor4]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5] - 8*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q82", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5] - 
     4*AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5] + 2*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5] - 
     4*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5] + 2*CA^2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5] - 
     8*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5] + 4*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5] - 
     8*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5] + 4*CA^2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5] + 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
     12*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 6*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
     2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
     2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
     24*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 12*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
     4*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
     4*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
     12*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 6*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
     2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
     2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
     24*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 12*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
     4*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
     4*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
     12*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 6*CA^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
     2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - CA^2*D*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
     2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + CA^2*D^2*Condensate[{"Q83", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
     24*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 12*CA^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
     4*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
     4*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
     6*CA^2*D*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 2*CA^2*D^2*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
     12*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 6*CA^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
     2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + CA^2*D*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
     2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - CA^2*D^2*Condensate[{"Q83", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
     24*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 12*CA^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
     4*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 2*CA^2*D*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
     4*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 2*CA^2*D^2*Condensate[{"Q85", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
     24*CA^2*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 8*CA^2*D*Condensate[{"Q81", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
     24*CA^2*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 8*CA^2*D*Condensate[{"Q81", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
     8*CA^2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4]*quarkMass[f] - 4*CA^2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4]*quarkMass[f] + 
     32*AGammaD[lor2, lor3, lor5]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4]*quarkMass[f] - 16*CA^2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor4]*quarkMass[f] - 
     8*CA^2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q72", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4]*quarkMass[f] - 4*CA^2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q72", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4]*quarkMass[f] - 
     32*AGammaD[lor2, lor3, lor5]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4]*quarkMass[f] + 16*CA^2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor4]*quarkMass[f] + 
     8*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5]*quarkMass[f] + 4*CA^2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5]*quarkMass[f] - 
     32*AGammaD[lor2, lor3, lor4]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5]*quarkMass[f] + 16*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor1, lor5]*quarkMass[f] + 
     8*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5]*quarkMass[f] + 4*CA^2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5]*quarkMass[f] + 
     32*AGammaD[lor2, lor3, lor4]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5]*quarkMass[f] - 16*CA^2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor1, lor5]*quarkMass[f] - 
     16*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4]*quarkMass[f] + 8*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4]*quarkMass[f] - 
     16*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4]*quarkMass[f] + 8*CA^2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor4]*quarkMass[f] + 
     16*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4]*quarkMass[f] - 8*CA^2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4]*quarkMass[f] + 
     16*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4]*quarkMass[f] - 8*CA^2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor4]*quarkMass[f] + 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] - 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] + 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] - 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] + 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] - 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] + 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] - 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] + 
     16*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5]*quarkMass[f] - 8*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5]*quarkMass[f] + 
     16*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5]*quarkMass[f] - 8*CA^2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor2, lor5]*quarkMass[f] - 
     16*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5]*quarkMass[f] + 8*CA^2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5]*quarkMass[f] - 
     16*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5]*quarkMass[f] + 8*CA^2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor2, lor5]*quarkMass[f] - 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] + 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] + 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] + 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] - 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] + 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] + 
     16*AGammaD[lor1, lor2, lor5]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4]*quarkMass[f] - 8*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4]*quarkMass[f] + 
     16*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4]*quarkMass[f] - 8*CA^2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor4]*quarkMass[f] - 
     16*AGammaD[lor1, lor2, lor5]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4]*quarkMass[f] + 8*CA^2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4]*quarkMass[f] - 
     16*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4]*quarkMass[f] + 8*CA^2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor4]*quarkMass[f] - 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] + 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] - 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] + 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] - 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*quarkMass[f] + 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*quarkMass[f] - 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*quarkMass[f] + 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*quarkMass[f] + 
     24*CA^2*D*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4]*quarkMass[f] - 8*CA^2*D^2*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4]*quarkMass[f] + 
     24*CA^2*D*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4]*quarkMass[f] - 8*CA^2*D^2*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4]*quarkMass[f] - 
     16*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5]*quarkMass[f] + 8*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5]*quarkMass[f] - 
     16*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5]*quarkMass[f] + 8*CA^2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*SUNT[colm] . SUNT[coln]*MTD[lor3, lor5]*quarkMass[f] + 
     16*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5]*quarkMass[f] - 8*CA^2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5]*quarkMass[f] + 
     16*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5]*quarkMass[f] - 8*CA^2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*SUNT[coln] . SUNT[colm]*MTD[lor3, lor5]*quarkMass[f] + 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] + 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] + 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] - 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
     24*CA^2*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] - 8*CA^2*D*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
     24*CA^2*D*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] + 8*CA^2*D^2*Condensate[{"Q71", f}]*SUNT[colm] . SUNT[coln]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
     24*CA^2*D*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] + 8*CA^2*D^2*Condensate[{"Q71", f}]*SUNT[coln] . SUNT[colm]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
     4*AGammaD[lor3, lor4, lor5]*MTD[lor1, lor2]*(SUNT[colm] . SUNT[coln]*(-4*Condensate[{"Q82", f}] + 2*D*Condensate[{"Q82", f}] - 2*Condensate[{"Q83", f}] - 4*Condensate[{"Q85", f}] + 8*Condensate[{"Q73", f}]*quarkMass[f] + 
         CA^2*(4*Condensate[{"A", f}] + 2*D*Condensate[{"A", f}] - (-2 + D)*Condensate[{"Q82", f}] + Condensate[{"Q83", f}] + 2*Condensate[{"Q85", f}] + (2 + D)*Condensate[{"Q72", f}]*quarkMass[f] - 4*Condensate[{"Q73", f}]*quarkMass[f])) + 
       SUNT[coln] . SUNT[colm]*(4*Condensate[{"Q82", f}] - 2*D*Condensate[{"Q82", f}] + 2*Condensate[{"Q83", f}] + 4*Condensate[{"Q85", f}] - 8*Condensate[{"Q73", f}]*quarkMass[f] + 
         CA^2*(4*Condensate[{"A", f}] + 2*D*Condensate[{"A", f}] + (-2 + D)*Condensate[{"Q82", f}] - Condensate[{"Q83", f}] - 2*Condensate[{"Q85", f}] + (2 + D)*Condensate[{"Q72", f}]*quarkMass[f] + 4*Condensate[{"Q73", f}]*quarkMass[f]))) + 
     4*AGammaD[lor2, lor4, lor5]*MTD[lor1, lor3]*(SUNT[colm] . SUNT[coln]*(-4*Condensate[{"Q82", f}] + 2*D*Condensate[{"Q82", f}] - 2*Condensate[{"Q83", f}] - 4*Condensate[{"Q85", f}] + 8*Condensate[{"Q73", f}]*quarkMass[f] + 
         CA^2*(4*Condensate[{"A", f}] + 2*D*Condensate[{"A", f}] - (-2 + D)*Condensate[{"Q82", f}] + Condensate[{"Q83", f}] + 2*Condensate[{"Q85", f}] + (2 + D)*Condensate[{"Q72", f}]*quarkMass[f] - 4*Condensate[{"Q73", f}]*quarkMass[f])) + 
       SUNT[coln] . SUNT[colm]*(4*Condensate[{"Q82", f}] - 2*D*Condensate[{"Q82", f}] + 2*Condensate[{"Q83", f}] + 4*Condensate[{"Q85", f}] - 8*Condensate[{"Q73", f}]*quarkMass[f] + 
         CA^2*(4*Condensate[{"A", f}] + 2*D*Condensate[{"A", f}] + (-2 + D)*Condensate[{"Q82", f}] - Condensate[{"Q83", f}] - 2*Condensate[{"Q85", f}] + (2 + D)*Condensate[{"Q72", f}]*quarkMass[f] + 4*Condensate[{"Q73", f}]*quarkMass[f]))))))/
 (CA*(2 - 3*CA^2 + CA^4)*(-4 + D)*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)));
    
    
If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* <q_iG(DG)q_j> *)
D8Condensate[{lor1_,lor2:Except[_List]},{lor3_,{lor4_,lor5:Except[_List]}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=(((I/16)*(4*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*MTD[lor1, lor5] - 4*D*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*MTD[lor1, lor5] + 
   6*AGammaD[lor2, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor1, lor5] - 2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor1, lor5] + 
   8*AGammaD[lor1, lor4, lor5]*Condensate[{"A", f}]*MTD[lor2, lor3] - 8*D*AGammaD[lor1, lor4, lor5]*Condensate[{"A", f}]*MTD[lor2, lor3] + 
   12*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor3] - 4*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor3] + 
   3*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 4*D*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
   D^2*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 6*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
   2*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 3*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
   4*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - D^2*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
   6*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 2*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
   4*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*MTD[lor2, lor4] - 4*D*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*MTD[lor2, lor4] + 
   6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor4] - 2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor4] - 
   3*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 4*D*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
   D^2*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 6*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
   2*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 6*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
   8*D*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 2*D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
   12*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 4*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
   4*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*MTD[lor2, lor5] + 4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*MTD[lor2, lor5] - 
   6*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor2, lor5] + 2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor2, lor5] + 
   3*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 4*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
   D^2*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 6*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
   2*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 6*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
   8*D*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 2*D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
   12*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 4*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
   24*AGammaD[lor1, lor2, lor5]*Condensate[{"A", f}]*MTD[lor3, lor4] + 12*AGammaD[lor1, lor2, lor5]*Condensate[{"Q86", f}]*MTD[lor3, lor4] - 
   4*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q86", f}]*MTD[lor3, lor4] - 9*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
   3*D*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 6*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
   4*D*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 2*D^2*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
   9*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 3*D*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
   6*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 4*D*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 
   2*D^2*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 24*AGammaD[lor1, lor2, lor4]*Condensate[{"A", f}]*MTD[lor3, lor5] - 
   12*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*MTD[lor3, lor5] + 4*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*MTD[lor3, lor5] + 
   9*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 3*D*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
   6*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 4*D*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
   2*D^2*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 9*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
   3*D*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 6*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
   4*D*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 2*D^2*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
   2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor1, lor5]*quarkMass[f] - 2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor1, lor5]*quarkMass[f] + 
   4*AGammaD[lor1, lor4, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor3]*quarkMass[f] - 4*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor3]*quarkMass[f] + 
   2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor4]*quarkMass[f] - 2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor4]*quarkMass[f] - 
   2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor2, lor5]*quarkMass[f] + 2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor2, lor5]*quarkMass[f] + 
   12*AGammaD[lor1, lor2, lor5]*Condensate[{"Q72", f}]*MTD[lor3, lor4]*quarkMass[f] - 12*AGammaD[lor1, lor2, lor4]*Condensate[{"Q72", f}]*MTD[lor3, lor5]*quarkMass[f] + 
   4*AGammaD[lor2, lor4, lor5]*MTD[lor1, lor3]*(2*(-1 + D)*Condensate[{"A", f}] + (-3 + D)*Condensate[{"Q86", f}] + (-1 + D)*Condensate[{"Q72", f}]*quarkMass[f]) + 
   2*AGammaD[lor2, lor3, lor5]*MTD[lor1, lor4]*(2*(-1 + D)*Condensate[{"A", f}] + (-3 + D)*Condensate[{"Q86", f}] + (-1 + D)*Condensate[{"Q72", f}]*quarkMass[f])))/
 ((-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)));
    
    
If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]



(* <q_i(DG)Gq_j> *)
D8Condensate[{lor1_,{lor2_,lor3:Except[_List]}},{lor4_,lor5:Except[_List]},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=(((I/16)*(-8*AGammaD[lor2, lor3, lor5]*Condensate[{"A", f}]*MTD[lor1, lor4] + 8*D*AGammaD[lor2, lor3, lor5]*Condensate[{"A", f}]*MTD[lor1, lor4] - 
   12*AGammaD[lor2, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor1, lor4] + 4*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor1, lor4] + 
   8*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*MTD[lor1, lor5] - 8*D*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*MTD[lor1, lor5] + 
   12*AGammaD[lor2, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor1, lor5] - 4*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor1, lor5] - 
   4*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*MTD[lor2, lor4] + 4*D*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*MTD[lor2, lor4] - 
   6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor4] + 2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor4] - 
   9*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 3*D*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
   6*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 4*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
   2*D^2*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 3*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
   4*D*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
   6*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 2*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
   4*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*MTD[lor2, lor5] - 4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*MTD[lor2, lor5] + 
   6*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor2, lor5] - 2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor2, lor5] + 
   9*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 3*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
   6*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 4*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
   2*D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 3*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
   4*D*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
   6*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 2*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
   4*AGammaD[lor1, lor2, lor5]*Condensate[{"A", f}]*MTD[lor3, lor4] - 4*D*AGammaD[lor1, lor2, lor5]*Condensate[{"A", f}]*MTD[lor3, lor4] + 
   6*AGammaD[lor1, lor2, lor5]*Condensate[{"Q86", f}]*MTD[lor3, lor4] - 2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q86", f}]*MTD[lor3, lor4] + 
   9*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 3*D*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
   6*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 4*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
   2*D^2*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 3*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
   4*D*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + D^2*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
   6*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 2*D*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
   6*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 8*D*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
   2*D^2*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 12*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
   4*D*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 4*AGammaD[lor1, lor2, lor4]*Condensate[{"A", f}]*MTD[lor3, lor5] + 
   4*D*AGammaD[lor1, lor2, lor4]*Condensate[{"A", f}]*MTD[lor3, lor5] - 6*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*MTD[lor3, lor5] + 
   2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*MTD[lor3, lor5] - 9*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
   3*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 6*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
   4*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 2*D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
   3*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 4*D*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
   D^2*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 6*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
   2*D*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 6*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
   8*D*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 2*D^2*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
   12*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 4*D*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
   4*AGammaD[lor2, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor1, lor4]*quarkMass[f] + 4*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor1, lor4]*quarkMass[f] + 
   4*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor1, lor5]*quarkMass[f] - 4*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor1, lor5]*quarkMass[f] - 
   2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor4]*quarkMass[f] + 2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor4]*quarkMass[f] + 
   2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor2, lor5]*quarkMass[f] - 2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor2, lor5]*quarkMass[f] + 
   2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q72", f}]*MTD[lor3, lor4]*quarkMass[f] - 2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q72", f}]*MTD[lor3, lor4]*quarkMass[f] - 
   2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q72", f}]*MTD[lor3, lor5]*quarkMass[f] + 2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q72", f}]*MTD[lor3, lor5]*quarkMass[f] + 
   4*AGammaD[lor3, lor4, lor5]*MTD[lor1, lor2]*(6*Condensate[{"A", f}] - (-3 + D)*Condensate[{"Q86", f}] + 3*Condensate[{"Q72", f}]*quarkMass[f]) - 
   4*AGammaD[lor2, lor4, lor5]*MTD[lor1, lor3]*(6*Condensate[{"A", f}] - (-3 + D)*Condensate[{"Q86", f}] + 3*Condensate[{"Q72", f}]*quarkMass[f])))/
 ((-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)));
    
    
If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* \[LeftAngleBracket]q_i (DDDG)^n q_j *)
D8Condensate[{coln_,{lor1_,{lor2_,{lor3_,{lor4_,lor5_}}}}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=SUNT[coln]/(CA CF)((Condensate[{"Q83", f}]*(GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - D*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
    GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + D*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 3*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
    GAD[lor5]*(2*(-1 + D)*MTD[lor1, lor4]*MTD[lor2, lor3] + (-1 + D)*MTD[lor1, lor3]*MTD[lor2, lor4] - 3*MTD[lor1, lor2]*MTD[lor3, lor4]) - 
    3*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + GAD[lor4]*(-2*(-1 + D)*MTD[lor1, lor5]*MTD[lor2, lor3] - (-1 + D)*MTD[lor1, lor3]*MTD[lor2, lor5] + 
      3*MTD[lor1, lor2]*MTD[lor3, lor5])) - 
  2*(Condensate[{"Q85", f}]*(GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
      GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - D*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
      GAD[lor5]*(-(D*MTD[lor1, lor4]*MTD[lor2, lor3]) - (-1 + D)*MTD[lor1, lor3]*MTD[lor2, lor4] + 3*MTD[lor1, lor2]*MTD[lor3, lor4]) + 
      GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + D*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
      GAD[lor4]*(D*MTD[lor1, lor5]*MTD[lor2, lor3] + (-1 + D)*MTD[lor1, lor3]*MTD[lor2, lor5] - 3*MTD[lor1, lor2]*MTD[lor3, lor5])) + 
    (-2 + D)*Condensate[{"Q84", f}]*(GAD[lor5]*(MTD[lor1, lor4]*MTD[lor2, lor3] + MTD[lor1, lor3]*MTD[lor2, lor4] + MTD[lor1, lor2]*MTD[lor3, lor4]) - 
      GAD[lor4]*(MTD[lor1, lor5]*MTD[lor2, lor3] + MTD[lor1, lor3]*MTD[lor2, lor5] + MTD[lor1, lor2]*MTD[lor3, lor5]))))/(8*(-2 + D)*(-1 + D)*D*(2 + D)));
    
    
If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]



(* <q_i (DDG)^n Dq_j> *)
(*D8Condensate[{coln_,{lor1_,{lor2_,{lor3_,lor4_}}}},lor5_,f:Except[_Rule],ops___Rule]:=-D8Condensate[lor1,{coln,{lor2,{lor3,{lor4,lor5}}}},f,ops]/.{lor2->lor1,lor3->lor2,lor4->lor3,lor5->lor4,lor1->lor5}/.{ga_GAD:>-ga,AGammaD[lors__]/;Length[{lors}]==5:>-AGammaD[lors]}
*)

D8Condensate[{coln_,{lor1_,{lor2_,{lor3_,lor4_}}}},lor5_,f:Except[_Rule],OptionsPattern[]]:=Block[{cond},
cond=-D8Condensate[lor5,{coln,{lor1,{lor2,{lor3,lor4}}}},f,Explicit->False,Massless->OptionValue[Massless],Factorization->OptionValue[Factorization]]/.{ga_GAD:>-ga,AGammaD[lors__]/;Length[{lors}]==5:>-AGammaD[lors]};

(*If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];*)
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]



(* <q_i (<-D) (DDG)^n q_j> *)
D8Condensate[lor1_,{coln_,{lor2_,{lor3_,{lor4_,lor5_}}}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=SUNT[coln]/(CA CF)((-4*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor1, lor5] + 2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor1, lor5] - 
  2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor1, lor5] - 4*AGammaD[lor2, lor3, lor4]*Condensate[{"Q85", f}]*MTD[lor1, lor5] + 
  4*AGammaD[lor1, lor4, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor3] + 8*AGammaD[lor1, lor4, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor3] - 
  4*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor3] - 12*AGammaD[lor1, lor4, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor3] + 
  4*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor3] + 4*AGammaD[lor1, lor4, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor3] - 
  4*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor3] - 6*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  8*D*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 2*D^2*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  12*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 10*D*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  2*D^2*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 6*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
  2*D^2*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 6*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  8*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 2*D^2*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  12*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 10*D*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  2*D^2*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 6*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
  2*D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor4] + 
  4*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor4] - 2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor4] - 
  6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor4] + 2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor4] + 
  2*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor4] - 2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor4] - 
  3*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 4*D*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  D^2*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 12*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  10*D*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 2*D^2*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  6*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 8*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  2*D^2*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 3*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
  4*D*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
  6*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 2*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
  2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor2, lor5] - 4*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor2, lor5] + 
  2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor2, lor5] + 6*AGammaD[lor1, lor3, lor4]*Condensate[{"Q84", f}]*MTD[lor2, lor5] - 
  2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q84", f}]*MTD[lor2, lor5] - 2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*MTD[lor2, lor5] + 
  2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*MTD[lor2, lor5] + 3*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  4*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + D^2*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  12*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 10*D*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  2*D^2*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 6*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  8*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 2*D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  3*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 4*D*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 6*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  2*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 6*AGammaD[lor1, lor2, lor5]*Condensate[{"Q82", f}]*MTD[lor3, lor4] + 
  6*AGammaD[lor1, lor2, lor5]*Condensate[{"Q83", f}]*MTD[lor3, lor4] - 6*AGammaD[lor1, lor2, lor5]*Condensate[{"Q84", f}]*MTD[lor3, lor4] + 
  2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q84", f}]*MTD[lor3, lor4] + 6*AGammaD[lor1, lor2, lor5]*Condensate[{"Q85", f}]*MTD[lor3, lor4] + 
  2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q85", f}]*MTD[lor3, lor4] - 9*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
  3*D*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 12*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
  10*D*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 2*D^2*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
  18*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 6*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
  9*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 3*D*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
  6*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 4*D*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
  2*D^2*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 6*AGammaD[lor1, lor2, lor4]*Condensate[{"Q82", f}]*MTD[lor3, lor5] - 
  6*AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*MTD[lor3, lor5] + 6*AGammaD[lor1, lor2, lor4]*Condensate[{"Q84", f}]*MTD[lor3, lor5] - 
  2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q84", f}]*MTD[lor3, lor5] - 6*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*MTD[lor3, lor5] - 
  2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*MTD[lor3, lor5] + 9*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  3*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 12*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  10*D*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 2*D^2*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  18*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 6*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  9*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 3*D*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
  6*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 4*D*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
  2*D^2*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 8*AGammaD[lor2, lor3, lor4]*Condensate[{"Q73", f}]*MTD[lor1, lor5]*
   quarkMass[f] - 32*AGammaD[lor1, lor4, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor3]*quarkMass[f] + 
  16*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor3]*quarkMass[f] - 24*AGammaD[lor1, lor4, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor3]*
   quarkMass[f] + 8*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor3]*quarkMass[f] - 
  16*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor4]*quarkMass[f] + 8*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor4]*
   quarkMass[f] - 12*AGammaD[lor1, lor3, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor4]*quarkMass[f] + 
  4*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor4]*quarkMass[f] + 16*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*MTD[lor2, lor5]*
   quarkMass[f] - 8*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*MTD[lor2, lor5]*quarkMass[f] + 
  12*AGammaD[lor1, lor3, lor4]*Condensate[{"Q74", f}]*MTD[lor2, lor5]*quarkMass[f] - 4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q74", f}]*MTD[lor2, lor5]*
   quarkMass[f] - 24*AGammaD[lor1, lor2, lor5]*Condensate[{"Q73", f}]*MTD[lor3, lor4]*quarkMass[f] - 
  12*AGammaD[lor1, lor2, lor5]*Condensate[{"Q74", f}]*MTD[lor3, lor4]*quarkMass[f] + 4*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q74", f}]*MTD[lor3, lor4]*
   quarkMass[f] + 24*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*MTD[lor3, lor5]*quarkMass[f] + 
  12*AGammaD[lor1, lor2, lor4]*Condensate[{"Q74", f}]*MTD[lor3, lor5]*quarkMass[f] - 4*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q74", f}]*MTD[lor3, lor5]*
   quarkMass[f] + 4*AGammaD[lor2, lor4, lor5]*MTD[lor1, lor3]*(-((-2 + D)*Condensate[{"Q82", f}]) + Condensate[{"Q83", f}] + 2*Condensate[{"Q85", f}] - 
    4*Condensate[{"Q73", f}]*quarkMass[f]) + 2*AGammaD[lor2, lor3, lor5]*MTD[lor1, lor4]*(-((-2 + D)*Condensate[{"Q82", f}]) + Condensate[{"Q83", f}] + 
    2*Condensate[{"Q85", f}] - 4*Condensate[{"Q73", f}]*quarkMass[f]))/(16*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)));
    
    
If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* <q_i (DG)^n DDq_j > *)
(*D8Condensate[{coln_,{lor1_,{lor2_,lor3_}}},lor4_,lor5_,f:Except[_Rule],ops___Rule]:= -D8Condensate[lor1,lor2,{coln,{lor3,{lor4,lor5}}},f,ops]/.{lor3->lor1,lor4->lor2,lor5->lor3,lor2->lor4,lor1->lor5}/.{ga_GAD:>-ga,AGammaD[lors__]/;Length[{lors}]==5:>-AGammaD[lors]}
*)(* the Option Explicit->True causing the AgammaD written explicitly, thus the charge conjugation by replacement rule {ga_GAD:>-ga,AGammaD[lors__]/;Length[{lors}]==5:>-AGammaD[lors]} will not work *)
D8Condensate[{coln_,{lor1_,{lor2_,lor3_}}},lor4_,lor5_,f:Except[_Rule],OptionsPattern[]]:=Block[{cond},
cond= -D8Condensate[lor5,lor4,{coln,{lor1,{lor2,lor3}}},f,Explicit->False,Massless->OptionValue[Massless],Factorization->OptionValue[Factorization]]/.{ga_GAD:>-ga,AGammaD[lors__]/;Length[{lors}]==5:>-AGammaD[lors]};

(*If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];*)
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* <q_i(<-D)(<-D) (DG)^n q_j > *)
D8Condensate[lor1_,lor2_,{coln_,{lor3_,{lor4_,lor5_}}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=((2*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*MTD[lor1, lor5] - 2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*MTD[lor1, lor5] + 
  3*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor1, lor5] - D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor1, lor5] + 
  3*AGammaD[lor2, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor1, lor5] - D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor1, lor5] - 
  3*AGammaD[lor2, lor3, lor4]*Condensate[{"Q84", f}]*MTD[lor1, lor5] + D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q84", f}]*MTD[lor1, lor5] + 
  3*AGammaD[lor2, lor3, lor4]*Condensate[{"Q85", f}]*MTD[lor1, lor5] - D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q85", f}]*MTD[lor1, lor5] + 
  3*AGammaD[lor2, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor1, lor5] - D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor1, lor5] + 
  4*AGammaD[lor1, lor4, lor5]*Condensate[{"A", f}]*MTD[lor2, lor3] - 4*D*AGammaD[lor1, lor4, lor5]*Condensate[{"A", f}]*MTD[lor2, lor3] - 
  6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor3] + 2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor3] - 
  6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor3] + 2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor3] + 
  6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor3] - 2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor3] - 
  6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor3] + 2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor3] + 
  6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor3] - 2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor3] - 
  6*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 8*D*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
  2*D^2*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 3*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
  4*D*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + D^2*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  3*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 4*D*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  D^2*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 3*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  4*D*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - D^2*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  3*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 4*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  D^2*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 9*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  6*D*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - D^2*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  6*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 8*D*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
  2*D^2*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 3*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
  4*D*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - D^2*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  3*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 4*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  D^2*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 3*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  4*D*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + D^2*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  3*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 4*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 9*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  6*D*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + D^2*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
  2*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*MTD[lor2, lor4] - 2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*MTD[lor2, lor4] - 
  3*AGammaD[lor1, lor3, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor4] + D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor4] - 
  3*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor4] + D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor4] + 
  3*AGammaD[lor1, lor3, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor4] - D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor4] - 
  3*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor4] + D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor4] + 
  3*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor4] - D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor4] - 
  6*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 8*D*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  2*D^2*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 3*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  4*D*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + D^2*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  3*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 4*D*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  D^2*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 9*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  6*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + D^2*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  9*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 6*D*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  D^2*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 3*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
  4*D*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
  6*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 2*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
  2*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*MTD[lor2, lor5] + 2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*MTD[lor2, lor5] + 
  3*AGammaD[lor1, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor2, lor5] - D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor2, lor5] + 
  3*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor2, lor5] - D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor2, lor5] - 
  3*AGammaD[lor1, lor3, lor4]*Condensate[{"Q84", f}]*MTD[lor2, lor5] + D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q84", f}]*MTD[lor2, lor5] + 
  3*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*MTD[lor2, lor5] - D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*MTD[lor2, lor5] - 
  3*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor2, lor5] + D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor2, lor5] + 
  6*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 8*D*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  2*D^2*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 3*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  4*D*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - D^2*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  3*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 4*D*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  D^2*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 9*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  6*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  9*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 6*D*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  D^2*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 3*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  4*D*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  6*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 2*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
  12*AGammaD[lor1, lor2, lor5]*Condensate[{"A", f}]*MTD[lor3, lor4] + 6*AGammaD[lor1, lor2, lor5]*Condensate[{"Q86", f}]*MTD[lor3, lor4] - 
  2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q86", f}]*MTD[lor3, lor4] - 12*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
  4*D*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 6*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
  2*D*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 3*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
  D*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 6*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
  2*D*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 12*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
  10*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 2*D^2*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
  6*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 8*D*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
  2*D^2*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 6*Condensate[{"A", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
  2*D*Condensate[{"A", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 3*Condensate[{"Q82", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
  D*Condensate[{"Q82", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 6*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
  2*D*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 3*Condensate[{"Q84", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
  D*Condensate[{"Q84", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 3*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
  7*D*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 2*D^2*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
  3*Condensate[{"Q86", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + D*Condensate[{"Q86", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
  6*Condensate[{"A", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 2*D*Condensate[{"A", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 
  3*Condensate[{"Q82", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + D*Condensate[{"Q82", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
  3*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - D*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
  3*Condensate[{"Q84", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - D*Condensate[{"Q84", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
  9*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 3*D*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 
  3*Condensate[{"Q86", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + D*Condensate[{"Q86", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 
  12*AGammaD[lor1, lor2, lor4]*Condensate[{"A", f}]*MTD[lor3, lor5] - 6*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*MTD[lor3, lor5] + 
  2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*MTD[lor3, lor5] + 12*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  4*D*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 6*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  2*D*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 3*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 6*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  2*D*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 12*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  10*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 2*D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  6*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 8*D*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  2*D^2*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 6*Condensate[{"A", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
  2*D*Condensate[{"A", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 3*Condensate[{"Q82", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
  D*Condensate[{"Q82", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 6*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
  2*D*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 3*Condensate[{"Q84", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
  D*Condensate[{"Q84", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 3*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
  7*D*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 2*D^2*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
  3*Condensate[{"Q86", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - D*Condensate[{"Q86", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
  6*Condensate[{"A", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 2*D*Condensate[{"A", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
  3*Condensate[{"Q82", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - D*Condensate[{"Q82", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
  3*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + D*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
  3*Condensate[{"Q84", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + D*Condensate[{"Q84", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
  9*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 3*D*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
  3*Condensate[{"Q86", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - D*Condensate[{"Q86", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
  AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor1, lor5]*quarkMass[f] - D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor1, lor5]*
   quarkMass[f] - 12*AGammaD[lor2, lor3, lor4]*Condensate[{"Q73", f}]*MTD[lor1, lor5]*quarkMass[f] + 
  4*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q73", f}]*MTD[lor1, lor5]*quarkMass[f] - 6*AGammaD[lor2, lor3, lor4]*Condensate[{"Q74", f}]*MTD[lor1, lor5]*
   quarkMass[f] + 2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q74", f}]*MTD[lor1, lor5]*quarkMass[f] + 
  2*AGammaD[lor1, lor4, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor3]*quarkMass[f] - 2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor3]*
   quarkMass[f] + 24*AGammaD[lor1, lor4, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor3]*quarkMass[f] - 
  8*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor3]*quarkMass[f] + 12*AGammaD[lor1, lor4, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor3]*
   quarkMass[f] - 4*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor3]*quarkMass[f] - 
  3*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] + 4*D*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*
   quarkMass[f] - D^2*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] + 
  6*Condensate[{"Q74", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] - 2*D*Condensate[{"Q74", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*
   quarkMass[f] + 3*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] - 
  4*D*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] + D^2*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*
   quarkMass[f] - 6*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] + 
  2*D*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] + AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor4]*
   quarkMass[f] - D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor4]*quarkMass[f] + 
  12*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor4]*quarkMass[f] - 4*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor4]*
   quarkMass[f] + 6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor4]*quarkMass[f] - 
  2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor4]*quarkMass[f] - 3*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*
   quarkMass[f] + 4*D*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] - 
  D^2*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] + 6*Condensate[{"Q74", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*
   quarkMass[f] - 2*D*Condensate[{"Q74", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] - 
  AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor2, lor5]*quarkMass[f] + D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor2, lor5]*
   quarkMass[f] - 12*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*MTD[lor2, lor5]*quarkMass[f] + 
  4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*MTD[lor2, lor5]*quarkMass[f] - 6*AGammaD[lor1, lor3, lor4]*Condensate[{"Q74", f}]*MTD[lor2, lor5]*
   quarkMass[f] + 2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q74", f}]*MTD[lor2, lor5]*quarkMass[f] + 
  3*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 4*D*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*
   quarkMass[f] + D^2*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 
  6*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] + 2*D*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*
   quarkMass[f] + 6*AGammaD[lor1, lor2, lor5]*Condensate[{"Q72", f}]*MTD[lor3, lor4]*quarkMass[f] - 
  6*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] + 2*D*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*
   quarkMass[f] - 12*Condensate[{"Q74", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] + 
  4*D*Condensate[{"Q74", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] + 3*Condensate[{"Q72", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*
   quarkMass[f] - D*Condensate[{"Q72", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*quarkMass[f] - 
  6*Condensate[{"Q74", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*quarkMass[f] + 8*D*Condensate[{"Q74", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*
   quarkMass[f] - 2*D^2*Condensate[{"Q74", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*quarkMass[f] + 
  3*Condensate[{"Q72", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4]*quarkMass[f] - D*Condensate[{"Q72", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4]*
   quarkMass[f] - 6*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4]*quarkMass[f] + 
  8*D*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4]*quarkMass[f] - 2*D^2*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor5]*
   MTD[lor3, lor4]*quarkMass[f] - 6*AGammaD[lor1, lor2, lor4]*Condensate[{"Q72", f}]*MTD[lor3, lor5]*quarkMass[f] + 
  6*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 2*D*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*
   quarkMass[f] + 12*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 
  4*D*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 3*Condensate[{"Q72", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*
   quarkMass[f] + D*Condensate[{"Q72", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
  6*Condensate[{"Q74", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] - 8*D*Condensate[{"Q74", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*
   quarkMass[f] + 2*D^2*Condensate[{"Q74", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
  3*Condensate[{"Q72", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] + D*Condensate[{"Q72", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*
   quarkMass[f] + 6*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
  8*D*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] + 2*D^2*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor4]*
   MTD[lor3, lor5]*quarkMass[f] + 24*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f]^2 - 
  20*D*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f]^2 + 4*D^2*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor2]*
   MTD[lor3, lor4]*quarkMass[f]^2 + 24*Condensate[{"Q6", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*quarkMass[f]^2 - 
  20*D*Condensate[{"Q6", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*quarkMass[f]^2 + 4*D^2*Condensate[{"Q6", f}]*GAD[lor2]*MTD[lor1, lor5]*
   MTD[lor3, lor4]*quarkMass[f]^2 + 24*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4]*quarkMass[f]^2 - 
  20*D*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4]*quarkMass[f]^2 + 4*D^2*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor5]*
   MTD[lor3, lor4]*quarkMass[f]^2 - 24*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f]^2 + 
  20*D*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f]^2 - 4*D^2*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor2]*
   MTD[lor3, lor5]*quarkMass[f]^2 - 24*Condensate[{"Q6", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 + 
  20*D*Condensate[{"Q6", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 - 4*D^2*Condensate[{"Q6", f}]*GAD[lor2]*MTD[lor1, lor4]*
   MTD[lor3, lor5]*quarkMass[f]^2 - 24*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 + 
  20*D*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 - 4*D^2*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor4]*
   MTD[lor3, lor5]*quarkMass[f]^2 + 2*AGammaD[lor2, lor4, lor5]*MTD[lor1, lor3]*(-2*Condensate[{"A", f}] + 2*D*Condensate[{"A", f}] + 
    (-3 + D)*Condensate[{"Q82", f}] - 3*Condensate[{"Q83", f}] + D*Condensate[{"Q83", f}] + 3*Condensate[{"Q84", f}] - D*Condensate[{"Q84", f}] - 
    3*Condensate[{"Q85", f}] + D*Condensate[{"Q85", f}] - 3*Condensate[{"Q86", f}] + D*Condensate[{"Q86", f}] + 
    (-1 + D)*Condensate[{"Q72", f}]*quarkMass[f] - 2*(-3 + D)*(2*Condensate[{"Q73", f}] + Condensate[{"Q74", f}])*quarkMass[f]) + 
  AGammaD[lor2, lor3, lor5]*MTD[lor1, lor4]*(-2*Condensate[{"A", f}] + 2*D*Condensate[{"A", f}] + (-3 + D)*Condensate[{"Q82", f}] - 
    3*Condensate[{"Q83", f}] + D*Condensate[{"Q83", f}] + 3*Condensate[{"Q84", f}] - D*Condensate[{"Q84", f}] - 3*Condensate[{"Q85", f}] + 
    D*Condensate[{"Q85", f}] - 3*Condensate[{"Q86", f}] + D*Condensate[{"Q86", f}] + (-1 + D)*Condensate[{"Q72", f}]*quarkMass[f] - 
    2*(-3 + D)*(2*Condensate[{"Q73", f}] + Condensate[{"Q74", f}])*quarkMass[f])) SUNT[coln]/(16*CA*CF*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)));


    
    
If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]




(* <q_i(<-D) (DG)^n Dq_j > *)
D8Condensate[lor1_,{coln_,{lor2_,{lor3_,lor4_}}},lor5_,f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=((-2*AGammaD[lor2, lor3, lor5]*Condensate[{"A", f}]*MTD[lor1, lor4] + 
  2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"A", f}]*MTD[lor1, lor4] + 
  AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*MTD[lor1, lor4] - 
  D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*MTD[lor1, lor4] - 
  AGammaD[lor2, lor3, lor5]*Condensate[{"Q83", f}]*MTD[lor1, lor4] + 
  D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q83", f}]*MTD[lor1, lor4] + 
  3*AGammaD[lor2, lor3, lor5]*Condensate[{"Q84", f}]*MTD[lor1, lor4] - 
  D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q84", f}]*MTD[lor1, lor4] + 
  AGammaD[lor2, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor1, lor4] + 
  D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor1, lor4] - 
  3*AGammaD[lor2, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor1, lor4] + 
  D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor1, lor4] + 
  12*AGammaD[lor1, lor4, lor5]*Condensate[{"A", f}]*MTD[lor2, lor3] - 
  6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor3] + 
  6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor3] - 
  6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor3] + 
  2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor3] + 
  6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor3] + 
  2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor3] + 
  6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor3] - 
  2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor3] - 
  6*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  2*D*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  3*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
  D*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
  3*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  D*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
  3*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  D*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
  9*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  3*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  3*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
  D*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
  12*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  4*D*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  6*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
  2*D*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
  6*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  2*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  6*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
  8*D*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  2*D^2*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
  6*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
  4*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  2*D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
  6*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  8*D*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
  2*D^2*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
  12*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*MTD[lor2, lor4] + 
  6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor4] - 
  6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor4] + 
  6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor4] - 
  2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor4] - 
  6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor4] - 
  2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor4] - 
  6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor4] + 
  2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor4] + 
  6*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  2*D*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  3*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
  D*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
  3*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  D*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
  3*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  D*Condensate[{"Q84", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
  9*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  3*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  3*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
  D*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
  12*Condensate[{"A", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
  4*D*Condensate[{"A", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
  6*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
  2*D*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
  6*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
  2*D*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
  6*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
  8*D*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
  2*D^2*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
  6*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
  4*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
  2*D^2*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
  6*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
  8*D*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
  2*D^2*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
  4*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*MTD[lor2, lor5] + 
  4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*MTD[lor2, lor5] + 
  2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor2, lor5] - 
  2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor2, lor5] - 
  2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor2, lor5] + 
  2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor2, lor5] + 
  6*AGammaD[lor1, lor3, lor4]*Condensate[{"Q84", f}]*MTD[lor2, lor5] - 
  2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q84", f}]*MTD[lor2, lor5] + 
  2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*MTD[lor2, lor5] + 
  2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*MTD[lor2, lor5] - 
  6*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor2, lor5] + 
  2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor2, lor5] + 
  6*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  8*D*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  2*D^2*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  3*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  4*D*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  D^2*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  3*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  4*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  D^2*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  9*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  6*D*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  D^2*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  3*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  2*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  9*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  6*D*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
  D^2*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
  6*Condensate[{"A", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
  8*D*Condensate[{"A", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  2*D^2*Condensate[{"A", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
  3*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  4*D*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
  D^2*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  3*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
  4*D*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
  9*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  6*D*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
  D^2*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
  3*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
  2*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  D^2*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  9*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
  6*D*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  D^2*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
  2*AGammaD[lor1, lor2, lor4]*Condensate[{"A", f}]*MTD[lor3, lor5] + 
  2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"A", f}]*MTD[lor3, lor5] + 
  AGammaD[lor1, lor2, lor4]*Condensate[{"Q82", f}]*MTD[lor3, lor5] - 
  D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q82", f}]*MTD[lor3, lor5] - 
  AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*MTD[lor3, lor5] + 
  D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*MTD[lor3, lor5] + 
  3*AGammaD[lor1, lor2, lor4]*Condensate[{"Q84", f}]*MTD[lor3, lor5] - 
  D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q84", f}]*MTD[lor3, lor5] + 
  AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*MTD[lor3, lor5] + 
  D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*MTD[lor3, lor5] - 
  3*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*MTD[lor3, lor5] + 
  D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*MTD[lor3, lor5] + 
  6*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  8*D*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  2*D^2*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  3*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  4*D*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  D^2*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  3*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  4*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  D^2*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  9*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  6*D*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  D^2*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  3*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  2*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  9*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
  6*D*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  D^2*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
  6*Condensate[{"A", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
  2*D*Condensate[{"A", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
  3*Condensate[{"Q82", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
  D*Condensate[{"Q82", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
  3*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
  D*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
  3*Condensate[{"Q84", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
  D*Condensate[{"Q84", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
  9*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
  3*D*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
  3*Condensate[{"Q86", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
  D*Condensate[{"Q86", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
  2*AGammaD[lor1, lor2, lor3]*Condensate[{"A", f}]*MTD[lor4, lor5] - 
  2*D*AGammaD[lor1, lor2, lor3]*Condensate[{"A", f}]*MTD[lor4, lor5] - 
  AGammaD[lor1, lor2, lor3]*Condensate[{"Q82", f}]*MTD[lor4, lor5] + 
  D*AGammaD[lor1, lor2, lor3]*Condensate[{"Q82", f}]*MTD[lor4, lor5] + 
  AGammaD[lor1, lor2, lor3]*Condensate[{"Q83", f}]*MTD[lor4, lor5] - 
  D*AGammaD[lor1, lor2, lor3]*Condensate[{"Q83", f}]*MTD[lor4, lor5] - 
  3*AGammaD[lor1, lor2, lor3]*Condensate[{"Q84", f}]*MTD[lor4, lor5] + 
  D*AGammaD[lor1, lor2, lor3]*Condensate[{"Q84", f}]*MTD[lor4, lor5] - 
  AGammaD[lor1, lor2, lor3]*Condensate[{"Q85", f}]*MTD[lor4, lor5] - 
  D*AGammaD[lor1, lor2, lor3]*Condensate[{"Q85", f}]*MTD[lor4, lor5] + 
  3*AGammaD[lor1, lor2, lor3]*Condensate[{"Q86", f}]*MTD[lor4, lor5] - 
  D*AGammaD[lor1, lor2, lor3]*Condensate[{"Q86", f}]*MTD[lor4, lor5] - 
  6*Condensate[{"A", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
  8*D*Condensate[{"A", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
  2*D^2*Condensate[{"A", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
  3*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
  4*D*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
  D^2*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
  3*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
  4*D*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
  D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
  9*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
  6*D*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
  D^2*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
  3*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
  2*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
  D^2*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
  9*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
  6*D*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
  D^2*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
  6*Condensate[{"A", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 
  2*D*Condensate[{"A", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 
  3*Condensate[{"Q82", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] - 
  D*Condensate[{"Q82", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] - 
  3*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 
  D*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] - 
  3*Condensate[{"Q84", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 
  D*Condensate[{"Q84", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] - 
  9*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 
  3*D*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 
  3*Condensate[{"Q86", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] - 
  D*Condensate[{"Q86", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] - 
  AGammaD[lor2, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor1, lor4]*quarkMass[f] + 
  D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor1, lor4]*quarkMass[f] + 
  4*AGammaD[lor2, lor3, lor5]*Condensate[{"Q73", f}]*MTD[lor1, lor4]*quarkMass[f] - 
  4*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q73", f}]*MTD[lor1, lor4]*quarkMass[f] + 
  6*AGammaD[lor2, lor3, lor5]*Condensate[{"Q74", f}]*MTD[lor1, lor4]*quarkMass[f] - 
  2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q74", f}]*MTD[lor1, lor4]*quarkMass[f] + 
  6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor3]*quarkMass[f] - 
  24*AGammaD[lor1, lor4, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor3]*quarkMass[f] - 
  12*AGammaD[lor1, lor4, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor3]*quarkMass[f] + 
  4*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor3]*quarkMass[f] - 
  3*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] + 
  D*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] + 
  6*Condensate[{"Q74", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] - 
  8*D*Condensate[{"Q74", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] + 
  2*D^2*Condensate[{"Q74", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] + 
  6*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] - 
  2*D*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] + 
  12*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] - 
  4*D*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] - 
  6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor4]*quarkMass[f] + 
  24*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor4]*quarkMass[f] + 
  12*AGammaD[lor1, lor3, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor4]*quarkMass[f] - 
  4*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor4]*quarkMass[f] + 
  3*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] - 
  D*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] - 
  6*Condensate[{"Q74", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] + 
  8*D*Condensate[{"Q74", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] - 
  2*D^2*Condensate[{"Q74", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] - 
  6*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] + 
  2*D*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] - 
  12*Condensate[{"Q74", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] + 
  4*D*Condensate[{"Q74", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] - 
  2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor2, lor5]*quarkMass[f] + 
  2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor2, lor5]*quarkMass[f] + 
  8*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*MTD[lor2, lor5]*quarkMass[f] - 
  8*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*MTD[lor2, lor5]*quarkMass[f] + 
  12*AGammaD[lor1, lor3, lor4]*Condensate[{"Q74", f}]*MTD[lor2, lor5]*quarkMass[f] - 
  4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q74", f}]*MTD[lor2, lor5]*quarkMass[f] + 
  3*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 
  4*D*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] + 
  D^2*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 
  6*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] + 
  2*D*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 
  3*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] + 
  4*D*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] - 
  D^2*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] + 
  6*Condensate[{"Q74", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] - 
  2*D*Condensate[{"Q74", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] - 
  AGammaD[lor1, lor2, lor4]*Condensate[{"Q72", f}]*MTD[lor3, lor5]*quarkMass[f] + 
  D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q72", f}]*MTD[lor3, lor5]*quarkMass[f] + 
  4*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*MTD[lor3, lor5]*quarkMass[f] - 
  4*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*MTD[lor3, lor5]*quarkMass[f] + 
  6*AGammaD[lor1, lor2, lor4]*Condensate[{"Q74", f}]*MTD[lor3, lor5]*quarkMass[f] - 
  2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q74", f}]*MTD[lor3, lor5]*quarkMass[f] + 
  3*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 
  4*D*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] + 
  D^2*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 
  6*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] + 
  2*D*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] + 
  3*Condensate[{"Q72", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
  D*Condensate[{"Q72", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
  6*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
  8*D*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
  2*D^2*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
  AGammaD[lor1, lor2, lor3]*Condensate[{"Q72", f}]*MTD[lor4, lor5]*quarkMass[f] - 
  D*AGammaD[lor1, lor2, lor3]*Condensate[{"Q72", f}]*MTD[lor4, lor5]*quarkMass[f] - 
  4*AGammaD[lor1, lor2, lor3]*Condensate[{"Q73", f}]*MTD[lor4, lor5]*quarkMass[f] + 
  4*D*AGammaD[lor1, lor2, lor3]*Condensate[{"Q73", f}]*MTD[lor4, lor5]*quarkMass[f] - 
  6*AGammaD[lor1, lor2, lor3]*Condensate[{"Q74", f}]*MTD[lor4, lor5]*quarkMass[f] + 
  2*D*AGammaD[lor1, lor2, lor3]*Condensate[{"Q74", f}]*MTD[lor4, lor5]*quarkMass[f] - 
  3*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f] + 
  4*D*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f] - 
  D^2*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f] + 
  6*Condensate[{"Q74", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f] - 
  2*D*Condensate[{"Q74", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f] - 
  3*Condensate[{"Q72", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f] + 
  D*Condensate[{"Q72", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f] + 
  6*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f] - 
  8*D*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f] + 
  2*D^2*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f] - 
  24*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f]^2 + 
  20*D*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f]^2 - 
  4*D^2*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f]^2 - 
  24*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f]^2 + 
  20*D*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f]^2 - 
  4*D^2*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f]^2 + 
  24*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f]^2 - 
  20*D*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f]^2 + 
  4*D^2*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f]^2 + 
  24*Condensate[{"Q6", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f]^2 - 
  20*D*Condensate[{"Q6", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f]^2 + 
  4*D^2*Condensate[{"Q6", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f]^2 + 
  24*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 - 
  20*D*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 + 
  4*D^2*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 - 
  24*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f]^2 + 
  20*D*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f]^2 - 
  4*D^2*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f]^2 - 
  2*AGammaD[lor3, lor4, lor5]*MTD[lor1, lor2]*(-2*Condensate[{"A", f}] + 2*D*Condensate[{"A", f}] + 
    Condensate[{"Q82", f}] - D*Condensate[{"Q82", f}] - Condensate[{"Q83", f}] + D*Condensate[{"Q83", f}] + 
    3*Condensate[{"Q84", f}] - D*Condensate[{"Q84", f}] + Condensate[{"Q85", f}] + D*Condensate[{"Q85", f}] - 
    3*Condensate[{"Q86", f}] + D*Condensate[{"Q86", f}] + 
    ((-1 + D)*Condensate[{"Q72", f}] - 4*(-1 + D)*Condensate[{"Q73", f}] - 2*(-3 + D)*Condensate[{"Q74", f}])*
     quarkMass[f]) - AGammaD[lor2, lor4, lor5]*MTD[lor1, lor3]*(-2*Condensate[{"A", f}] + 
    2*D*Condensate[{"A", f}] + Condensate[{"Q82", f}] - D*Condensate[{"Q82", f}] - Condensate[{"Q83", f}] + 
    D*Condensate[{"Q83", f}] + 3*Condensate[{"Q84", f}] - D*Condensate[{"Q84", f}] + Condensate[{"Q85", f}] + 
    D*Condensate[{"Q85", f}] - 3*Condensate[{"Q86", f}] + D*Condensate[{"Q86", f}] + 
    ((-1 + D)*Condensate[{"Q72", f}] - 4*(-1 + D)*Condensate[{"Q73", f}] - 2*(-3 + D)*Condensate[{"Q74", f}])*
     quarkMass[f])) SUNT[coln]/(16*CA*CF*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)));


    
    
If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* <q_i G^n DDDq_j > *)
(*D8Condensate[{coln_,{lor1_,lor2_}},lor3_,lor4_,lor5_,f:Except[_Rule],ops___Rule]:=-D8Condensate[lor1,lor2,lor3,{coln,{lor4,lor5}},f,ops]/.{lor4->lor1,lor5->lor2,lor2->lor4,lor1->lor5}/.{ga_GAD:>-ga,AGammaD[lors__]/;Length[{lors}]==5:>-AGammaD[lors]}
*)

D8Condensate[{coln_,{lor1_,lor2_}},lor3_,lor4_,lor5_,f:Except[_Rule],OptionsPattern[]]:=Block[{cond},
cond= -D8Condensate[lor5,lor4,lor3,{coln,{lor1,lor2}},f,Explicit->False,Massless->OptionValue[Massless],Factorization->OptionValue[Factorization]]/.{ga_GAD:>-ga,AGammaD[lors__]/;Length[{lors}]==5:>-AGammaD[lors]};

(*If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];*)
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]




(* <q_i(<-D)(<-D)(<-D) G^n q_j > *)
D8Condensate[lor1_,lor2_,lor3_,{coln_,{lor4_,lor5_}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=((-16*(2 + D)*AGammaD[lor1, lor2, lor3, lor4, lor5]*Condensate[{"A", f}] + 
  (-4 + D)*(-8*AGammaD[lor2, lor3, lor5]*Condensate[{"A", f}]*MTD[lor1, lor4] - 
    4*D*AGammaD[lor2, lor3, lor5]*Condensate[{"A", f}]*MTD[lor1, lor4] + 
    4*AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*MTD[lor1, lor4] - 
    2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*MTD[lor1, lor4] + 
    2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q83", f}]*MTD[lor1, lor4] + 
    4*AGammaD[lor2, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor1, lor4] + 
    8*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*MTD[lor1, lor5] + 4*D*AGammaD[lor2, lor3, lor4]*
     Condensate[{"A", f}]*MTD[lor1, lor5] - 4*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor1, lor5] + 
    2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor1, lor5] - 
    2*AGammaD[lor2, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor1, lor5] - 
    4*AGammaD[lor2, lor3, lor4]*Condensate[{"Q85", f}]*MTD[lor1, lor5] - 
    16*AGammaD[lor1, lor4, lor5]*Condensate[{"A", f}]*MTD[lor2, lor3] + 
    8*D*AGammaD[lor1, lor4, lor5]*Condensate[{"A", f}]*MTD[lor2, lor3] + 
    12*AGammaD[lor1, lor4, lor5]*Condensate[{"Q81", f}]*MTD[lor2, lor3] - 
    4*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q81", f}]*MTD[lor2, lor3] + 
    8*AGammaD[lor1, lor4, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor3] - 
    4*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor3] + 
    10*AGammaD[lor1, lor4, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor3] - 
    2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor3] + 
    20*AGammaD[lor1, lor4, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor3] - 
    4*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor3] - 
    24*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor3] + 
    8*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor3] + 
    24*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 20*D*Condensate[{"A", f}]*GAD[lor5]*
     MTD[lor1, lor4]*MTD[lor2, lor3] + 4*D^2*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
    6*D*Condensate[{"Q81", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
    2*D^2*Condensate[{"Q81", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
    12*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
    10*D*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
    2*D^2*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
    6*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 5*D*Condensate[{"Q83", f}]*GAD[lor5]*
     MTD[lor1, lor4]*MTD[lor2, lor3] - D^2*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
    12*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
    10*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
    2*D^2*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
    24*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
    20*D*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
    4*D^2*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
    24*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 20*D*Condensate[{"A", f}]*GAD[lor4]*
     MTD[lor1, lor5]*MTD[lor2, lor3] - 4*D^2*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
    6*D*Condensate[{"Q81", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
    2*D^2*Condensate[{"Q81", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
    12*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
    10*D*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
    2*D^2*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
    6*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 5*D*Condensate[{"Q83", f}]*GAD[lor4]*
     MTD[lor1, lor5]*MTD[lor2, lor3] + D^2*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
    12*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
    10*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
    2*D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
    24*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
    20*D*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
    4*D^2*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
    16*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*MTD[lor2, lor4] + 
    4*D*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*MTD[lor2, lor4] + 
    12*AGammaD[lor1, lor3, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor4] - 
    2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor4] - 
    2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor4] - 
    4*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor4] - 
    12*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor4] + 
    4*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor4] + 
    24*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 20*D*Condensate[{"A", f}]*GAD[lor5]*
     MTD[lor1, lor3]*MTD[lor2, lor4] + 4*D^2*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
    12*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
    10*D*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
    2*D^2*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
    6*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 2*D*Condensate[{"Q83", f}]*GAD[lor5]*
     MTD[lor1, lor3]*MTD[lor2, lor4] - 12*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
    16*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
    4*D^2*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
    24*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
    20*D*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
    4*D^2*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
    12*Condensate[{"Q81", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
    4*D*Condensate[{"Q81", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
    6*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 8*D*Condensate[{"Q83", f}]*GAD[lor3]*
     MTD[lor1, lor5]*MTD[lor2, lor4] + 2*D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
    12*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
    4*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
    16*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*MTD[lor2, lor5] - 
    4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*MTD[lor2, lor5] - 
    12*AGammaD[lor1, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor2, lor5] + 
    2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor2, lor5] + 
    2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor2, lor5] + 
    4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*MTD[lor2, lor5] + 
    12*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor2, lor5] - 
    4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor2, lor5] - 
    24*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 20*D*Condensate[{"A", f}]*GAD[lor4]*
     MTD[lor1, lor3]*MTD[lor2, lor5] - 4*D^2*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
    12*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
    10*D*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
    2*D^2*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
    6*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 2*D*Condensate[{"Q83", f}]*GAD[lor4]*
     MTD[lor1, lor3]*MTD[lor2, lor5] + 12*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
    16*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
    4*D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
    24*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
    20*D*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
    4*D^2*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
    12*Condensate[{"Q81", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
    4*D*Condensate[{"Q81", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
    6*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 8*D*Condensate[{"Q83", f}]*GAD[lor3]*
     MTD[lor1, lor4]*MTD[lor2, lor5] - 2*D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
    12*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
    4*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
    16*AGammaD[lor1, lor2, lor5]*Condensate[{"A", f}]*MTD[lor3, lor4] + 
    4*D*AGammaD[lor1, lor2, lor5]*Condensate[{"A", f}]*MTD[lor3, lor4] + 
    4*AGammaD[lor1, lor2, lor5]*Condensate[{"Q82", f}]*MTD[lor3, lor4] - 
    2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q82", f}]*MTD[lor3, lor4] + 
    2*AGammaD[lor1, lor2, lor5]*Condensate[{"Q83", f}]*MTD[lor3, lor4] + 
    4*AGammaD[lor1, lor2, lor5]*Condensate[{"Q85", f}]*MTD[lor3, lor4] - 
    12*AGammaD[lor1, lor2, lor5]*Condensate[{"Q86", f}]*MTD[lor3, lor4] + 
    4*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q86", f}]*MTD[lor3, lor4] + 
    24*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 20*D*Condensate[{"A", f}]*GAD[lor5]*
     MTD[lor1, lor2]*MTD[lor3, lor4] + 4*D^2*Condensate[{"A", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
    6*D*Condensate[{"Q81", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
    2*D^2*Condensate[{"Q81", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
    12*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
    10*D*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
    2*D^2*Condensate[{"Q82", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
    3*D*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
    D^2*Condensate[{"Q83", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
    24*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
    14*D*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
    2*D^2*Condensate[{"Q85", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
    24*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
    20*D*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
    4*D^2*Condensate[{"Q86", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4] + 
    12*Condensate[{"Q81", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
    10*D*Condensate[{"Q81", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
    2*D^2*Condensate[{"Q81", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
    12*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
    7*D*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
    D^2*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
    6*D*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] - 
    2*D^2*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4] + 
    12*Condensate[{"Q81", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] - 
    4*D*Condensate[{"Q81", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4] + 
    16*AGammaD[lor1, lor2, lor4]*Condensate[{"A", f}]*MTD[lor3, lor5] - 
    4*D*AGammaD[lor1, lor2, lor4]*Condensate[{"A", f}]*MTD[lor3, lor5] - 
    4*AGammaD[lor1, lor2, lor4]*Condensate[{"Q82", f}]*MTD[lor3, lor5] + 
    2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q82", f}]*MTD[lor3, lor5] - 
    2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*MTD[lor3, lor5] - 
    4*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*MTD[lor3, lor5] + 
    12*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*MTD[lor3, lor5] - 
    4*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q86", f}]*MTD[lor3, lor5] - 
    24*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 20*D*Condensate[{"A", f}]*GAD[lor4]*
     MTD[lor1, lor2]*MTD[lor3, lor5] - 4*D^2*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    6*D*Condensate[{"Q81", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
    2*D^2*Condensate[{"Q81", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    12*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
    10*D*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    2*D^2*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    3*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
    D^2*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    24*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
    14*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    2*D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
    24*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    20*D*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
    4*D^2*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
    12*Condensate[{"Q81", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
    10*D*Condensate[{"Q81", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
    2*D^2*Condensate[{"Q81", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
    12*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
    7*D*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
    D^2*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
    6*D*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
    2*D^2*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
    12*Condensate[{"Q81", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
    4*D*Condensate[{"Q81", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
    4*AGammaD[lor2, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor1, lor4]*quarkMass[f] - 
    2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor1, lor4]*quarkMass[f] - 
    8*AGammaD[lor2, lor3, lor5]*Condensate[{"Q73", f}]*MTD[lor1, lor4]*quarkMass[f] + 
    4*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor1, lor5]*quarkMass[f] + 
    2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor1, lor5]*quarkMass[f] + 
    8*AGammaD[lor2, lor3, lor4]*Condensate[{"Q73", f}]*MTD[lor1, lor5]*quarkMass[f] + 
    24*AGammaD[lor1, lor4, lor5]*Condensate[{"Q71", f}]*MTD[lor2, lor3]*quarkMass[f] - 
    8*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q71", f}]*MTD[lor2, lor3]*quarkMass[f] - 
    32*AGammaD[lor1, lor4, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor3]*quarkMass[f] + 
    12*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor3]*quarkMass[f] - 
    64*AGammaD[lor1, lor4, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor3]*quarkMass[f] + 
    16*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor3]*quarkMass[f] + 
    24*Condensate[{"Q71", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] - 
    8*D*Condensate[{"Q71", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] + 
    12*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] - 
    10*D*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] + 
    2*D^2*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] - 
    24*Condensate[{"Q71", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] + 
    8*D*Condensate[{"Q71", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] - 
    12*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] + 
    10*D*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] - 
    2*D^2*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] - 
    8*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor4]*quarkMass[f] + 
    2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor4]*quarkMass[f] + 
    8*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor4]*quarkMass[f] + 
    12*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] - 
    10*D*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] + 
    2*D^2*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] + 
    12*D*Condensate[{"Q71", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] - 
    4*D^2*Condensate[{"Q71", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] + 
    8*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor2, lor5]*quarkMass[f] - 
    2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor2, lor5]*quarkMass[f] - 
    8*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*MTD[lor2, lor5]*quarkMass[f] - 
    12*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] + 
    10*D*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 
    2*D^2*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 
    12*D*Condensate[{"Q71", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] + 
    4*D^2*Condensate[{"Q71", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] - 
    8*AGammaD[lor1, lor2, lor5]*Condensate[{"Q72", f}]*MTD[lor3, lor4]*quarkMass[f] + 
    2*D*AGammaD[lor1, lor2, lor5]*Condensate[{"Q72", f}]*MTD[lor3, lor4]*quarkMass[f] - 
    8*AGammaD[lor1, lor2, lor5]*Condensate[{"Q73", f}]*MTD[lor3, lor4]*quarkMass[f] - 
    24*Condensate[{"Q71", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] + 
    8*D*Condensate[{"Q71", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] + 
    12*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] - 
    10*D*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] + 
    2*D^2*Condensate[{"Q72", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] - 
    24*Condensate[{"Q71", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*quarkMass[f] + 
    20*D*Condensate[{"Q71", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*quarkMass[f] - 
    4*D^2*Condensate[{"Q71", f}]*GAD[lor2]*MTD[lor1, lor5]*MTD[lor3, lor4]*quarkMass[f] + 
    12*D*Condensate[{"Q71", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4]*quarkMass[f] - 
    4*D^2*Condensate[{"Q71", f}]*GAD[lor1]*MTD[lor2, lor5]*MTD[lor3, lor4]*quarkMass[f] + 
    8*AGammaD[lor1, lor2, lor4]*Condensate[{"Q72", f}]*MTD[lor3, lor5]*quarkMass[f] - 
    2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q72", f}]*MTD[lor3, lor5]*quarkMass[f] + 
    8*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*MTD[lor3, lor5]*quarkMass[f] + 
    24*Condensate[{"Q71", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 
    8*D*Condensate[{"Q71", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 
    12*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] + 
    10*D*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 
    2*D^2*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] + 
    24*Condensate[{"Q71", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
    20*D*Condensate[{"Q71", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
    4*D^2*Condensate[{"Q71", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
    12*D*Condensate[{"Q71", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
    4*D^2*Condensate[{"Q71", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
    24*AGammaD[lor1, lor4, lor5]*Condensate[{"Q6", f}]*MTD[lor2, lor3]*quarkMass[f]^2 - 
    8*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q6", f}]*MTD[lor2, lor3]*quarkMass[f]^2 - 
    24*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f]^2 + 
    20*D*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f]^2 - 
    4*D^2*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f]^2 + 
    24*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f]^2 - 
    20*D*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f]^2 + 
    4*D^2*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f]^2 - 
    24*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f]^2 + 
    20*D*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f]^2 - 
    4*D^2*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f]^2 + 
    24*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f]^2 - 
    20*D*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f]^2 + 
    4*D^2*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f]^2 - 
    24*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f]^2 + 
    20*D*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f]^2 - 
    4*D^2*Condensate[{"Q6", f}]*GAD[lor5]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f]^2 + 
    24*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f]^2 - 
    20*D*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f]^2 + 
    4*D^2*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f]^2 - 
    24*AGammaD[lor1, lor4, lor5]*Condensate[{"Q5", f}]*MTD[lor2, lor3]*quarkMass[f]^3 + 
    8*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q5", f}]*MTD[lor2, lor3]*quarkMass[f]^3 + 
    2*AGammaD[lor2, lor4, lor5]*MTD[lor1, lor3]*(8*Condensate[{"A", f}] + 6*Condensate[{"Q81", f}] - 
      2*D*Condensate[{"Q81", f}] + 3*Condensate[{"Q83", f}] - D*Condensate[{"Q83", f}] + 
      6*Condensate[{"Q85", f}] - 2*D*Condensate[{"Q85", f}] + 
      4*(-((-3 + D)*Condensate[{"Q71", f}]) + (-2 + D)*Condensate[{"Q72", f}] + 
        2*(-3 + D)*Condensate[{"Q73", f}])*quarkMass[f] - 4*(-3 + D)*Condensate[{"Q6", f}]*quarkMass[f]^2 + 
      4*(-3 + D)*Condensate[{"Q5", f}]*quarkMass[f]^3) + 2*AGammaD[lor3, lor4, lor5]*MTD[lor1, lor2]*
     (-4*D*Condensate[{"A", f}] + 6*Condensate[{"Q81", f}] - 2*D*Condensate[{"Q81", f}] - 
      4*Condensate[{"Q82", f}] + 2*D*Condensate[{"Q82", f}] + Condensate[{"Q83", f}] - 
      D*Condensate[{"Q83", f}] + 2*Condensate[{"Q85", f}] - 2*D*Condensate[{"Q85", f}] + 
      2*(-2*(-3 + D)*Condensate[{"Q71", f}] + (-6 + D)*Condensate[{"Q72", f}] + 
        4*(-2 + D)*Condensate[{"Q73", f}])*quarkMass[f] - 4*(-3 + D)*Condensate[{"Q6", f}]*quarkMass[f]^2 + 
      4*(-3 + D)*Condensate[{"Q5", f}]*quarkMass[f]^3))) SUNT[coln]/(32*CA*CF*(-4 + D)*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)));


    
    
If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* <q_i(<-D) G^n DDq_j > *)
(*D8Condensate[lor1_,{coln_,{lor2_,lor3_}},lor4:Except[_List],lor5:Except[_List],f:Except[_Rule],ops___Rule]:=-D8Condensate[lor1,lor2,{coln,{lor3,lor4}},lor5,f,ops]/.{lor3->lor2,lor4->lor3,lor5->lor1,lor2->lor4,lor1->lor5}/.{ga_GAD:>-ga,AGammaD[lors__]/;Length[{lors}]==5:>-AGammaD[lors]}
*)

D8Condensate[lor1_,{coln_,{lor2_,lor3_}},lor4:Except[_List],lor5:Except[_List],f:Except[_Rule],OptionsPattern[]]:=Block[{cond},
cond= -D8Condensate[lor5,lor4,{coln,{lor2,lor3}},lor1,f,Explicit->False,Massless->OptionValue[Massless],Factorization->OptionValue[Factorization]]/.{ga_GAD:>-ga,AGammaD[lors__]/;Length[{lors}]==5:>-AGammaD[lors]};

(*If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];*)
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]




(* <q_i(<-D)(<-D) G^n Dq_j > *)
D8Condensate[lor1_,lor2_,{coln_,{lor3_,lor4_}},lor5:Except[_List],f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=((16*(2 + D)*AGammaD[lor1, lor2, lor3, lor4, lor5]*Condensate[{"A", f}] - 
  (-4 + D)*(-12*AGammaD[lor2, lor3, lor5]*Condensate[{"A", f}]*MTD[lor1, lor4] - 
    2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q82", f}]*MTD[lor1, lor4] - 
    4*AGammaD[lor2, lor3, lor5]*Condensate[{"Q83", f}]*MTD[lor1, lor4] + 
    2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q83", f}]*MTD[lor1, lor4] + 
    6*AGammaD[lor2, lor3, lor5]*Condensate[{"Q84", f}]*MTD[lor1, lor4] - 
    2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q84", f}]*MTD[lor1, lor4] - 
    2*AGammaD[lor2, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor1, lor4] + 
    2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor1, lor4] - 
    6*AGammaD[lor2, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor1, lor4] + 
    2*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor1, lor4] + 
    8*AGammaD[lor2, lor3, lor4]*Condensate[{"A", f}]*MTD[lor1, lor5] + 8*D*AGammaD[lor2, lor3, lor4]*
     Condensate[{"A", f}]*MTD[lor1, lor5] + 12*AGammaD[lor2, lor3, lor4]*Condensate[{"Q81", f}]*
     MTD[lor1, lor5] - 4*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q81", f}]*MTD[lor1, lor5] - 
    12*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor1, lor5] + 
    4*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor1, lor5] - 
    6*AGammaD[lor2, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor1, lor5] + 
    2*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor1, lor5] + 
    12*AGammaD[lor2, lor3, lor4]*Condensate[{"Q84", f}]*MTD[lor1, lor5] - 
    4*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q84", f}]*MTD[lor1, lor5] - 
    12*AGammaD[lor2, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor1, lor5] + 
    4*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor1, lor5] + 
    12*AGammaD[lor1, lor4, lor5]*Condensate[{"A", f}]*MTD[lor2, lor3] - 
    6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor3] + 
    6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor3] - 
    6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor3] + 
    2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor3] + 
    6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor3] + 
    2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor3] + 
    6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor3] - 
    2*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor3] + 
    12*Condensate[{"Q81", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] - 
    4*D*Condensate[{"Q81", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3] + 
    12*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 4*D*Condensate[{"A", f}]*GAD[lor4]*
     MTD[lor1, lor5]*MTD[lor2, lor3] - 6*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
    2*D*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
    6*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 2*D*Condensate[{"Q83", f}]*GAD[lor4]*
     MTD[lor1, lor5]*MTD[lor2, lor3] - 6*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
    8*D*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
    2*D^2*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
    6*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 4*D*Condensate[{"Q85", f}]*GAD[lor4]*
     MTD[lor1, lor5]*MTD[lor2, lor3] - 2*D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] + 
    6*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 8*D*Condensate[{"Q86", f}]*GAD[lor4]*
     MTD[lor1, lor5]*MTD[lor2, lor3] + 2*D^2*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3] - 
    12*AGammaD[lor1, lor3, lor5]*Condensate[{"A", f}]*MTD[lor2, lor4] + 
    6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q82", f}]*MTD[lor2, lor4] - 
    6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q83", f}]*MTD[lor2, lor4] + 
    6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor4] - 
    2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q84", f}]*MTD[lor2, lor4] - 
    6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor4] - 
    2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q85", f}]*MTD[lor2, lor4] - 
    6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor4] + 
    2*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q86", f}]*MTD[lor2, lor4] - 
    12*Condensate[{"Q81", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
    4*D*Condensate[{"Q81", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4] - 
    12*Condensate[{"A", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 4*D*Condensate[{"A", f}]*GAD[lor3]*
     MTD[lor1, lor5]*MTD[lor2, lor4] + 6*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
    2*D*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
    6*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 2*D*Condensate[{"Q83", f}]*GAD[lor3]*
     MTD[lor1, lor5]*MTD[lor2, lor4] + 6*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
    8*D*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 
    2*D^2*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
    6*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 4*D*Condensate[{"Q85", f}]*GAD[lor3]*
     MTD[lor1, lor5]*MTD[lor2, lor4] + 2*D^2*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
    6*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] + 8*D*Condensate[{"Q86", f}]*GAD[lor3]*
     MTD[lor1, lor5]*MTD[lor2, lor4] - 2*D^2*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4] - 
    8*AGammaD[lor1, lor3, lor4]*Condensate[{"A", f}]*MTD[lor2, lor5] + 12*AGammaD[lor1, lor3, lor4]*
     Condensate[{"Q81", f}]*MTD[lor2, lor5] - 4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q81", f}]*
     MTD[lor2, lor5] - 4*AGammaD[lor1, lor3, lor4]*Condensate[{"Q82", f}]*MTD[lor2, lor5] - 
    2*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor2, lor5] + 
    2*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q83", f}]*MTD[lor2, lor5] + 
    12*AGammaD[lor1, lor3, lor4]*Condensate[{"Q84", f}]*MTD[lor2, lor5] - 
    4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q84", f}]*MTD[lor2, lor5] + 
    8*AGammaD[lor1, lor3, lor4]*Condensate[{"Q85", f}]*MTD[lor2, lor5] - 
    12*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor2, lor5] + 
    4*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q86", f}]*MTD[lor2, lor5] + 
    12*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 4*D*Condensate[{"A", f}]*GAD[lor4]*
     MTD[lor1, lor3]*MTD[lor2, lor5] + 6*D*Condensate[{"Q81", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
    2*D^2*Condensate[{"Q81", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
    6*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 2*D*Condensate[{"Q82", f}]*GAD[lor4]*
     MTD[lor1, lor3]*MTD[lor2, lor5] - 3*D*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
    D^2*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
    6*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 8*D*Condensate[{"Q84", f}]*GAD[lor4]*
     MTD[lor1, lor3]*MTD[lor2, lor5] - 2*D^2*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
    6*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 2*D*Condensate[{"Q85", f}]*GAD[lor4]*
     MTD[lor1, lor3]*MTD[lor2, lor5] + 6*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
    8*D*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] + 
    2*D^2*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5] - 
    12*Condensate[{"A", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 4*D*Condensate[{"A", f}]*GAD[lor3]*
     MTD[lor1, lor4]*MTD[lor2, lor5] - 6*D*Condensate[{"Q81", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
    2*D^2*Condensate[{"Q81", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
    6*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 2*D*Condensate[{"Q82", f}]*GAD[lor3]*
     MTD[lor1, lor4]*MTD[lor2, lor5] + 3*D*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
    D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
    6*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 8*D*Condensate[{"Q84", f}]*GAD[lor3]*
     MTD[lor1, lor4]*MTD[lor2, lor5] + 2*D^2*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
    6*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 2*D*Condensate[{"Q85", f}]*GAD[lor3]*
     MTD[lor1, lor4]*MTD[lor2, lor5] - 6*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
    8*D*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] - 
    2*D^2*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5] + 
    8*AGammaD[lor1, lor2, lor4]*Condensate[{"A", f}]*MTD[lor3, lor5] + 4*D*AGammaD[lor1, lor2, lor4]*
     Condensate[{"A", f}]*MTD[lor3, lor5] + 4*AGammaD[lor1, lor2, lor4]*Condensate[{"Q82", f}]*MTD[lor3, lor5] - 
    2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q82", f}]*MTD[lor3, lor5] + 
    2*AGammaD[lor1, lor2, lor4]*Condensate[{"Q83", f}]*MTD[lor3, lor5] + 
    4*AGammaD[lor1, lor2, lor4]*Condensate[{"Q85", f}]*MTD[lor3, lor5] - 
    12*D*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    4*D^2*Condensate[{"A", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
    6*D*Condensate[{"Q81", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    2*D^2*Condensate[{"Q81", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    6*D*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
    2*D^2*Condensate[{"Q82", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    6*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 5*D*Condensate[{"Q83", f}]*GAD[lor4]*
     MTD[lor1, lor2]*MTD[lor3, lor5] + D^2*Condensate[{"Q83", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
    12*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    4*D*Condensate[{"Q84", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
    6*D*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    2*D^2*Condensate[{"Q85", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    12*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] - 
    4*D*Condensate[{"Q86", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5] + 
    12*Condensate[{"A", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 4*D*Condensate[{"A", f}]*GAD[lor2]*
     MTD[lor1, lor4]*MTD[lor3, lor5] + 12*Condensate[{"Q81", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
    10*D*Condensate[{"Q81", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
    2*D^2*Condensate[{"Q81", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
    6*Condensate[{"Q82", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 2*D*Condensate[{"Q82", f}]*GAD[lor2]*
     MTD[lor1, lor4]*MTD[lor3, lor5] - 3*D*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
    D^2*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
    6*Condensate[{"Q84", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 2*D*Condensate[{"Q84", f}]*GAD[lor2]*
     MTD[lor1, lor4]*MTD[lor3, lor5] + 6*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
    8*D*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 
    2*D^2*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] - 
    6*Condensate[{"Q86", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5] + 2*D*Condensate[{"Q86", f}]*GAD[lor2]*
     MTD[lor1, lor4]*MTD[lor3, lor5] + 12*Condensate[{"A", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
    4*D*Condensate[{"A", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
    12*Condensate[{"Q81", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
    4*D*Condensate[{"Q81", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
    6*Condensate[{"Q82", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 2*D*Condensate[{"Q82", f}]*GAD[lor1]*
     MTD[lor2, lor4]*MTD[lor3, lor5] + 6*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
    2*D*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 
    6*Condensate[{"Q84", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 2*D*Condensate[{"Q84", f}]*GAD[lor1]*
     MTD[lor2, lor4]*MTD[lor3, lor5] + 18*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
    6*D*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] - 
    6*Condensate[{"Q86", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5] + 2*D*Condensate[{"Q86", f}]*GAD[lor1]*
     MTD[lor2, lor4]*MTD[lor3, lor5] - 8*AGammaD[lor1, lor2, lor3]*Condensate[{"A", f}]*MTD[lor4, lor5] - 
    4*D*AGammaD[lor1, lor2, lor3]*Condensate[{"A", f}]*MTD[lor4, lor5] - 
    4*AGammaD[lor1, lor2, lor3]*Condensate[{"Q82", f}]*MTD[lor4, lor5] + 
    2*D*AGammaD[lor1, lor2, lor3]*Condensate[{"Q82", f}]*MTD[lor4, lor5] - 
    2*AGammaD[lor1, lor2, lor3]*Condensate[{"Q83", f}]*MTD[lor4, lor5] - 
    4*AGammaD[lor1, lor2, lor3]*Condensate[{"Q85", f}]*MTD[lor4, lor5] + 
    12*D*Condensate[{"A", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
    4*D^2*Condensate[{"A", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
    6*D*Condensate[{"Q81", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
    2*D^2*Condensate[{"Q81", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
    6*D*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
    2*D^2*Condensate[{"Q82", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
    6*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 5*D*Condensate[{"Q83", f}]*GAD[lor3]*
     MTD[lor1, lor2]*MTD[lor4, lor5] - D^2*Condensate[{"Q83", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
    12*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
    4*D*Condensate[{"Q84", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
    6*D*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
    2*D^2*Condensate[{"Q85", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
    12*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] + 
    4*D*Condensate[{"Q86", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5] - 
    12*Condensate[{"A", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5] + 4*D*Condensate[{"A", f}]*GAD[lor2]*
     MTD[lor1, lor3]*MTD[lor4, lor5] - 12*Condensate[{"Q81", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5] + 
    10*D*Condensate[{"Q81", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5] - 
    2*D^2*Condensate[{"Q81", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5] + 
    6*Condensate[{"Q82", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5] - 2*D*Condensate[{"Q82", f}]*GAD[lor2]*
     MTD[lor1, lor3]*MTD[lor4, lor5] + 3*D*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5] - 
    D^2*Condensate[{"Q83", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5] - 
    6*Condensate[{"Q84", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5] + 2*D*Condensate[{"Q84", f}]*GAD[lor2]*
     MTD[lor1, lor3]*MTD[lor4, lor5] - 6*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5] + 
    8*D*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5] - 
    2*D^2*Condensate[{"Q85", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5] + 
    6*Condensate[{"Q86", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5] - 2*D*Condensate[{"Q86", f}]*GAD[lor2]*
     MTD[lor1, lor3]*MTD[lor4, lor5] - 12*Condensate[{"A", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 
    4*D*Condensate[{"A", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] - 
    12*Condensate[{"Q81", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 
    4*D*Condensate[{"Q81", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 
    6*Condensate[{"Q82", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] - 2*D*Condensate[{"Q82", f}]*GAD[lor1]*
     MTD[lor2, lor3]*MTD[lor4, lor5] - 6*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 
    2*D*Condensate[{"Q83", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] - 
    6*Condensate[{"Q84", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 2*D*Condensate[{"Q84", f}]*GAD[lor1]*
     MTD[lor2, lor3]*MTD[lor4, lor5] - 18*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 
    6*D*Condensate[{"Q85", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] + 
    6*Condensate[{"Q86", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5] - 2*D*Condensate[{"Q86", f}]*GAD[lor1]*
     MTD[lor2, lor3]*MTD[lor4, lor5] - 6*AGammaD[lor2, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor1, lor4]*
     quarkMass[f] + 16*AGammaD[lor2, lor3, lor5]*Condensate[{"Q73", f}]*MTD[lor1, lor4]*quarkMass[f] - 
    8*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q73", f}]*MTD[lor1, lor4]*quarkMass[f] + 
    12*AGammaD[lor2, lor3, lor5]*Condensate[{"Q74", f}]*MTD[lor1, lor4]*quarkMass[f] - 
    4*D*AGammaD[lor2, lor3, lor5]*Condensate[{"Q74", f}]*MTD[lor1, lor4]*quarkMass[f] + 
    24*AGammaD[lor2, lor3, lor4]*Condensate[{"Q71", f}]*MTD[lor1, lor5]*quarkMass[f] - 
    8*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q71", f}]*MTD[lor1, lor5]*quarkMass[f] - 
    20*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor1, lor5]*quarkMass[f] + 
    12*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor1, lor5]*quarkMass[f] + 
    24*AGammaD[lor2, lor3, lor4]*Condensate[{"Q74", f}]*MTD[lor1, lor5]*quarkMass[f] - 
    8*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q74", f}]*MTD[lor1, lor5]*quarkMass[f] + 
    6*AGammaD[lor1, lor4, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor3]*quarkMass[f] - 
    24*AGammaD[lor1, lor4, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor3]*quarkMass[f] - 
    12*AGammaD[lor1, lor4, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor3]*quarkMass[f] + 
    4*D*AGammaD[lor1, lor4, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor3]*quarkMass[f] + 
    12*D*Condensate[{"Q71", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] - 
    4*D^2*Condensate[{"Q71", f}]*GAD[lor5]*MTD[lor1, lor4]*MTD[lor2, lor3]*quarkMass[f] + 
    6*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] - 
    2*D*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] + 
    12*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] - 
    4*D*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f] - 
    6*AGammaD[lor1, lor3, lor5]*Condensate[{"Q72", f}]*MTD[lor2, lor4]*quarkMass[f] + 
    24*AGammaD[lor1, lor3, lor5]*Condensate[{"Q73", f}]*MTD[lor2, lor4]*quarkMass[f] + 
    12*AGammaD[lor1, lor3, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor4]*quarkMass[f] - 
    4*D*AGammaD[lor1, lor3, lor5]*Condensate[{"Q74", f}]*MTD[lor2, lor4]*quarkMass[f] - 
    12*D*Condensate[{"Q71", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] + 
    4*D^2*Condensate[{"Q71", f}]*GAD[lor5]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] - 
    6*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] + 
    2*D*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] - 
    12*Condensate[{"Q74", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] + 
    4*D*Condensate[{"Q74", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f] + 
    24*AGammaD[lor1, lor3, lor4]*Condensate[{"Q71", f}]*MTD[lor2, lor5]*quarkMass[f] - 
    8*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q71", f}]*MTD[lor2, lor5]*quarkMass[f] - 
    28*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor2, lor5]*quarkMass[f] + 
    8*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q72", f}]*MTD[lor2, lor5]*quarkMass[f] - 
    16*AGammaD[lor1, lor3, lor4]*Condensate[{"Q73", f}]*MTD[lor2, lor5]*quarkMass[f] + 
    24*AGammaD[lor1, lor3, lor4]*Condensate[{"Q74", f}]*MTD[lor2, lor5]*quarkMass[f] - 
    8*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q74", f}]*MTD[lor2, lor5]*quarkMass[f] + 
    24*Condensate[{"Q71", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 
    8*D*Condensate[{"Q71", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] + 
    6*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 
    2*D*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] + 
    12*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 
    4*D*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f] - 
    24*Condensate[{"Q71", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] + 
    8*D*Condensate[{"Q71", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] - 
    6*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] + 
    2*D*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] - 
    12*Condensate[{"Q74", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] + 
    4*D*Condensate[{"Q74", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f] + 
    4*AGammaD[lor1, lor2, lor4]*Condensate[{"Q72", f}]*MTD[lor3, lor5]*quarkMass[f] + 
    2*D*AGammaD[lor1, lor2, lor4]*Condensate[{"Q72", f}]*MTD[lor3, lor5]*quarkMass[f] - 
    8*AGammaD[lor1, lor2, lor4]*Condensate[{"Q73", f}]*MTD[lor3, lor5]*quarkMass[f] - 
    24*Condensate[{"Q71", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] + 
    8*D*Condensate[{"Q71", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 
    6*D*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] + 
    2*D^2*Condensate[{"Q72", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 
    24*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] + 
    8*D*Condensate[{"Q74", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f] - 
    24*Condensate[{"Q71", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
    20*D*Condensate[{"Q71", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
    4*D^2*Condensate[{"Q71", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
    6*Condensate[{"Q72", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
    2*D*Condensate[{"Q72", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
    12*Condensate[{"Q74", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
    16*D*Condensate[{"Q74", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
    4*D^2*Condensate[{"Q74", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
    12*D*Condensate[{"Q71", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
    4*D^2*Condensate[{"Q71", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
    6*Condensate[{"Q72", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
    2*D*Condensate[{"Q72", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
    12*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] + 
    16*D*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
    4*D^2*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f] - 
    4*AGammaD[lor1, lor2, lor3]*Condensate[{"Q72", f}]*MTD[lor4, lor5]*quarkMass[f] - 
    2*D*AGammaD[lor1, lor2, lor3]*Condensate[{"Q72", f}]*MTD[lor4, lor5]*quarkMass[f] + 
    8*AGammaD[lor1, lor2, lor3]*Condensate[{"Q73", f}]*MTD[lor4, lor5]*quarkMass[f] + 
    24*Condensate[{"Q71", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f] - 
    8*D*Condensate[{"Q71", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f] + 
    6*D*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f] - 
    2*D^2*Condensate[{"Q72", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f] + 
    24*Condensate[{"Q74", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f] - 
    8*D*Condensate[{"Q74", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f] + 
    24*Condensate[{"Q71", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5]*quarkMass[f] - 
    20*D*Condensate[{"Q71", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5]*quarkMass[f] + 
    4*D^2*Condensate[{"Q71", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5]*quarkMass[f] - 
    6*Condensate[{"Q72", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5]*quarkMass[f] + 
    2*D*Condensate[{"Q72", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5]*quarkMass[f] + 
    12*Condensate[{"Q74", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5]*quarkMass[f] - 
    16*D*Condensate[{"Q74", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5]*quarkMass[f] + 
    4*D^2*Condensate[{"Q74", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5]*quarkMass[f] - 
    12*D*Condensate[{"Q71", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f] + 
    4*D^2*Condensate[{"Q71", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f] - 
    6*Condensate[{"Q72", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f] + 
    2*D*Condensate[{"Q72", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f] + 
    12*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f] - 
    16*D*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f] + 
    4*D^2*Condensate[{"Q74", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f] + 
    24*AGammaD[lor2, lor3, lor4]*Condensate[{"Q6", f}]*MTD[lor1, lor5]*quarkMass[f]^2 - 
    8*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q6", f}]*MTD[lor1, lor5]*quarkMass[f]^2 - 
    24*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f]^2 + 
    20*D*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f]^2 - 
    4*D^2*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor5]*MTD[lor2, lor3]*quarkMass[f]^2 + 
    24*Condensate[{"Q6", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f]^2 - 
    20*D*Condensate[{"Q6", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f]^2 + 
    4*D^2*Condensate[{"Q6", f}]*GAD[lor3]*MTD[lor1, lor5]*MTD[lor2, lor4]*quarkMass[f]^2 + 
    24*AGammaD[lor1, lor3, lor4]*Condensate[{"Q6", f}]*MTD[lor2, lor5]*quarkMass[f]^2 - 
    8*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q6", f}]*MTD[lor2, lor5]*quarkMass[f]^2 - 
    24*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f]^2 + 
    20*D*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f]^2 - 
    4*D^2*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor3]*MTD[lor2, lor5]*quarkMass[f]^2 + 
    24*Condensate[{"Q6", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f]^2 - 
    20*D*Condensate[{"Q6", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f]^2 + 
    4*D^2*Condensate[{"Q6", f}]*GAD[lor3]*MTD[lor1, lor4]*MTD[lor2, lor5]*quarkMass[f]^2 + 
    24*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f]^2 - 
    20*D*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f]^2 + 
    4*D^2*Condensate[{"Q6", f}]*GAD[lor4]*MTD[lor1, lor2]*MTD[lor3, lor5]*quarkMass[f]^2 + 
    48*Condensate[{"Q6", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 - 
    40*D*Condensate[{"Q6", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 + 
    8*D^2*Condensate[{"Q6", f}]*GAD[lor2]*MTD[lor1, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 + 
    48*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 - 
    40*D*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 + 
    8*D^2*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor4]*MTD[lor3, lor5]*quarkMass[f]^2 - 
    24*Condensate[{"Q6", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f]^2 + 
    20*D*Condensate[{"Q6", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f]^2 - 
    4*D^2*Condensate[{"Q6", f}]*GAD[lor3]*MTD[lor1, lor2]*MTD[lor4, lor5]*quarkMass[f]^2 - 
    48*Condensate[{"Q6", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5]*quarkMass[f]^2 + 
    40*D*Condensate[{"Q6", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5]*quarkMass[f]^2 - 
    8*D^2*Condensate[{"Q6", f}]*GAD[lor2]*MTD[lor1, lor3]*MTD[lor4, lor5]*quarkMass[f]^2 - 
    48*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f]^2 + 
    40*D*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f]^2 - 
    8*D^2*Condensate[{"Q6", f}]*GAD[lor1]*MTD[lor2, lor3]*MTD[lor4, lor5]*quarkMass[f]^2 - 
    24*AGammaD[lor2, lor3, lor4]*Condensate[{"Q5", f}]*MTD[lor1, lor5]*quarkMass[f]^3 + 
    8*D*AGammaD[lor2, lor3, lor4]*Condensate[{"Q5", f}]*MTD[lor1, lor5]*quarkMass[f]^3 - 
    24*AGammaD[lor1, lor3, lor4]*Condensate[{"Q5", f}]*MTD[lor2, lor5]*quarkMass[f]^3 + 
    8*D*AGammaD[lor1, lor3, lor4]*Condensate[{"Q5", f}]*MTD[lor2, lor5]*quarkMass[f]^3 + 
    2*AGammaD[lor2, lor4, lor5]*MTD[lor1, lor3]*(6*Condensate[{"A", f}] + Condensate[{"Q82", f}] + 
      2*Condensate[{"Q83", f}] - D*Condensate[{"Q83", f}] - 3*Condensate[{"Q84", f}] + 
      D*Condensate[{"Q84", f}] + Condensate[{"Q85", f}] - D*Condensate[{"Q85", f}] + 3*Condensate[{"Q86", f}] - 
      D*Condensate[{"Q86", f}] + (3*Condensate[{"Q72", f}] + 4*(-2 + D)*Condensate[{"Q73", f}] + 
        2*(-3 + D)*Condensate[{"Q74", f}])*quarkMass[f]) + 2*AGammaD[lor3, lor4, lor5]*MTD[lor1, lor2]*
     (-4*D*Condensate[{"A", f}] + 6*Condensate[{"Q81", f}] - 2*D*Condensate[{"Q81", f}] - 
      4*Condensate[{"Q82", f}] + 2*D*Condensate[{"Q82", f}] + Condensate[{"Q83", f}] - 
      D*Condensate[{"Q83", f}] + 2*Condensate[{"Q85", f}] - 2*D*Condensate[{"Q85", f}] + 
      2*(-2*(-3 + D)*Condensate[{"Q71", f}] + (-6 + D)*Condensate[{"Q72", f}] + 
        4*(-2 + D)*Condensate[{"Q73", f}])*quarkMass[f] - 4*(-3 + D)*Condensate[{"Q6", f}]*quarkMass[f]^2 + 
      4*(-3 + D)*Condensate[{"Q5", f}]*quarkMass[f]^3))) SUNT[coln]/(32*CA*CF*(-4 + D)*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D)));


    
    
If[OptionValue[Factorization]==="EMFirst",
   cond=Factorization[cond,Method->"EMFirst"]//Simplify
   ,
   If[OptionValue[Factorization]==="VSFirst",
        cond=Factorization[cond,Method->"VSFirst"]//Simplify
   ]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


End[]
