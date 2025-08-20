(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



D7Condensate::usage = "D7Condensate[lor1_,lor2_,lor3_,lor4_,psi_] = \[LeftAngleBracket]\!\(\*SubsuperscriptBox[OverscriptBox[\(psi\), \(_\)], \(i\), \(a\)]\) (\!\(\*SuperscriptBox[\(\[Del]\), \(lor1\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor2\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor3\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor4\)]\)\!\(\*SubscriptBox[\(psi\), \(j\)]\)\!\(\*SuperscriptBox[\()\), \(b\)]\)\[RightAngleBracket]; "


Begin["`Private`D7Condensate`"]
Options[D7Condensate] = {
	Explicit->False,
	Massless->True,
	Factorization->False
}


(* \[LeftAngleBracket]q^a_i DDDDq^b_j\[RightAngleBracket]*)
D7Condensate[lor1_,lor2_,lor3_,lor4:Except[_Rule],f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=(-1/8*(AGammaD[lor1, lor2, lor3, lor4]*Condensate[{"Q72", f}])/(D*(-6 + 11*D - 6*D^2 + D^3)) + 
 ((I/8)*DiracSigma[GAD[lor1], GAD[lor4]]*MTD[lor2, lor3]*(-3*(-2 + D)*Condensate[{"Q71", f}] - 4*Condensate[{"Q72", f}] + D*Condensate[{"Q72", f}] - 8*Condensate[{"Q73", f}] + 
    2*D*Condensate[{"Q73", f}] + 2*Condensate[{"Q74", f}] - D*Condensate[{"Q74", f}] + 2*(-1 + D)*Condensate[{"Q6", f}]*quarkMass[f] + 
    (-4 + D)*Condensate[{"Q5", f}]*quarkMass[f]^2))/((-2 + D)*(-1 + D)*D*(2 + D)) + 
 ((I/8)*(DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] + DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4])*(-2*(-2 + D)*Condensate[{"Q71", f}] - 2*Condensate[{"Q72", f}] + 
    D*Condensate[{"Q72", f}] - 2*Condensate[{"Q73", f}] + 2*D*Condensate[{"Q73", f}] + 2*Condensate[{"Q74", f}] - D*Condensate[{"Q74", f}] + 
    (-2 + D)*Condensate[{"Q6", f}]*quarkMass[f] + (-2 + D)*Condensate[{"Q5", f}]*quarkMass[f]^2))/((-2 + D)*(-1 + D)*D*(2 + D)) + 
 (MTD[lor1, lor2]*MTD[lor3, lor4]*((-4 + D)*Condensate[{"Q71", f}] + Condensate[{"Q72", f}] - D*Condensate[{"Q72", f}] + 2*Condensate[{"Q73", f}] - 2*D*Condensate[{"Q73", f}] - 
    2*Condensate[{"Q74", f}] + 4*Condensate[{"Q6", f}]*quarkMass[f] - 2*(-1 + D)*Condensate[{"Q5", f}]*quarkMass[f]^2 + 2*(-1 + D)*Condensate[{"Q3", f}]*quarkMass[f]^4))/
  (8*(-1 + D)*D*(2 + D)) + (MTD[lor1, lor3]*MTD[lor2, lor4]*(2*(-1 + D)*Condensate[{"Q71", f}] + Condensate[{"Q72", f}] - D*Condensate[{"Q72", f}] + 2*Condensate[{"Q73", f}] - 
    2*D*Condensate[{"Q73", f}] + D*Condensate[{"Q74", f}] - 2*D*Condensate[{"Q6", f}]*quarkMass[f] - 2*(-1 + D)*Condensate[{"Q5", f}]*quarkMass[f]^2 + 
    2*(-1 + D)*Condensate[{"Q3", f}]*quarkMass[f]^4))/(8*(-1 + D)*D*(2 + D)) + 
 (MTD[lor1, lor4]*MTD[lor2, lor3]*(Condensate[{"Q72", f}] + 2*Condensate[{"Q73", f}] + D*(3*Condensate[{"Q71", f}] - Condensate[{"Q72", f}] - 2*Condensate[{"Q73", f}] + 
      Condensate[{"Q74", f}]) - 2*D*Condensate[{"Q6", f}]*quarkMass[f] - 2*(-1 + D)*Condensate[{"Q5", f}]*quarkMass[f]^2 + 2*(-1 + D)*Condensate[{"Q3", f}]*quarkMass[f]^4))/
  (8*(-1 + D)*D*(2 + D)) + ((I/8)*(DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2] + DiracSigma[GAD[lor1], GAD[lor2]]*MTD[lor3, lor4])*
   (-2*(-Condensate[{"Q71", f}] + Condensate[{"Q73", f}] + Condensate[{"Q6", f}]*quarkMass[f]) + 
    D*(-Condensate[{"Q71", f}] + Condensate[{"Q72", f}] + 2*Condensate[{"Q73", f}] + Condensate[{"Q5", f}]*quarkMass[f]^2)))/((-2 + D)*(-1 + D)*D*(2 + D)) + 
 ((I/8)*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4]*(2*(Condensate[{"Q71", f}] + Condensate[{"Q74", f}] - Condensate[{"Q6", f}]*quarkMass[f]) + 
    D*(-Condensate[{"Q71", f}] + Condensate[{"Q72", f}] - Condensate[{"Q74", f}] + Condensate[{"Q5", f}]*quarkMass[f]^2)))/((-2 + D)*(-1 + D)*D*(2 + D)))/CA;

        
If[OptionValue[Factorization],
   cond=Factorization[cond]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* \[LeftAngleBracket]q^a_i G^n G^m q^b_j\[RightAngleBracket] *)
D7Condensate[{coln_,{lor1_,lor2_}},{colm_,{lor3_,lor4_}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=(((SUNT[colm] . SUNT[coln] + SUNT[coln] . SUNT[colm])*(2*AGammaD[lor1, lor2, lor3, lor4]*Condensate[{"Q72", f}] + 
   (-3 + D)*((-I)*Condensate[{"Q73", f}]*(DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] - DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] - 
       DiracSigma[GAD[lor1], GAD[lor4]]*MTD[lor2, lor3] + DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4]) - 
     (-2 + D)*Condensate[{"Q71", f}]*(MTD[lor1, lor4]*MTD[lor2, lor3] - MTD[lor1, lor3]*MTD[lor2, lor4]))))/(2*(-2 + CA^2)*CF*(-3 + D)*(-2 + D)*(-1 + D)*D));

        
If[OptionValue[Factorization],
   cond=Factorization[cond]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* \[LeftAngleBracket]q_i GG q_j\[RightAngleBracket] *)
D7Condensate[{lor1_,lor2:Except[_List]},{lor3_,lor4:Except[_List]},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=((2*AGammaD[lor1, lor2, lor3, lor4]*Condensate[{"Q72", f}] + 
  (-3 + D)*((-I)*Condensate[{"Q73", f}]*(DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] - DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] - 
      DiracSigma[GAD[lor1], GAD[lor4]]*MTD[lor2, lor3] + DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4]) - 
    (-2 + D)*Condensate[{"Q71", f}]*(MTD[lor1, lor4]*MTD[lor2, lor3] - MTD[lor1, lor3]*MTD[lor2, lor4])))/(4*(-3 + D)*(-2 + D)*(-1 + D)*D));

        
If[OptionValue[Factorization],
   cond=Factorization[cond]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* \[LeftAngleBracket]q^a_j (DDG)^n q^b_j\[RightAngleBracket] *)
D7Condensate[{coln_,{lor1_,{lor2_,{lor3_,lor4_}}}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=-1/4*(((-2 + D)*Condensate[{"Q74", f}]*(2*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2] + DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] - 
      DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] + DiracSigma[GAD[lor1], GAD[lor4]]*MTD[lor2, lor3] - DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4]) + 
    2*Condensate[{"Q73", f}]*(2*(-1 + D)*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2] + (-1 + D)*DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] + 
      DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] - D*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] - 3*DiracSigma[GAD[lor1], GAD[lor4]]*MTD[lor2, lor3] + 
      3*DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4]))*SUNT[coln])/(CA*CF*(-2 + D)*(-1 + D)*D*(2 + D));

        
If[OptionValue[Factorization],
   cond=Factorization[cond]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* \[LeftAngleBracket]q^a_j (DG)^n (Dq_j)^b\[RightAngleBracket] *)
D7Condensate[{coln_,{lor1_,{lor2_,lor3_}}},lor4_,f:Except[_Rule],ops___Rule]:=-D7Condensate[lor1,{coln,{lor2,{lor3,lor4}}},f,ops]/.{lor2->lor1,lor3->lor2,lor4->lor3,lor1->lor4}/.{ds_DiracSigma:>-ds}

(* \[LeftAngleBracket]q_j(<-D)^a (DG)^n q_j^b\[RightAngleBracket] *)
D7Condensate[lor1_,{coln_,{lor2_,{lor3_,lor4_}}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=((2*Condensate[{"Q73", f}]*(2*(-1 + D)*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2] + (-1 + D)*DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] + 
     DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] - D*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] - 3*DiracSigma[GAD[lor1], GAD[lor4]]*MTD[lor2, lor3] + 
     3*DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4]) + 
   (-2 + D)*(Condensate[{"Q74", f}]*(2*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2] + DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] - 
       DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] + DiracSigma[GAD[lor1], GAD[lor4]]*MTD[lor2, lor3] + (2*I)*MTD[lor1, lor4]*MTD[lor2, lor3] + 
       I*D*MTD[lor1, lor4]*MTD[lor2, lor3] - DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4] - (2*I)*MTD[lor1, lor3]*MTD[lor2, lor4] - I*D*MTD[lor1, lor3]*MTD[lor2, lor4]) - 
     (2*I)*(2 + D)*Condensate[{"Q6", f}]*(MTD[lor1, lor4]*MTD[lor2, lor3] - MTD[lor1, lor3]*MTD[lor2, lor4])*quarkMass[f]))*SUNT[coln])/(8*CA*CF*(-2 + D)*(-1 + D)*D*(2 + D));

        
If[OptionValue[Factorization],
   cond=Factorization[cond]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* \[LeftAngleBracket]q^a_j (G)^n (DDq_j)^b\[RightAngleBracket] *)
D7Condensate[{coln_,{lor1_,lor2_}},lor3_,lor4_,f_,ops___Rule]:=-D7Condensate[lor1,lor2,{coln,{lor3,lor4}},f,ops]/.{lor3->lor1,lor4->lor2,lor1->lor4,lor2->lor3}/.{ds_DiracSigma:>-ds}

(* \[LeftAngleBracket]q_j(<-D)(<-D)^a (DG)^n q_j^b\[RightAngleBracket] *)
D7Condensate[lor1_,lor2_,{coln_,{lor3_,lor4_}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=(((-2*I)*(2 + D)*AGammaD[lor1, lor2, lor3, lor4]*Condensate[{"Q72", f}] + (-3 + D)*(4*Condensate[{"Q73", f}]*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2] - 
     4*D*Condensate[{"Q73", f}]*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2] + 2*Condensate[{"Q73", f}]*DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] - 
     2*D*Condensate[{"Q73", f}]*DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] - 2*Condensate[{"Q73", f}]*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] + 
     2*D*Condensate[{"Q73", f}]*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] + 6*Condensate[{"Q73", f}]*DiracSigma[GAD[lor1], GAD[lor4]]*MTD[lor2, lor3] - 
     6*Condensate[{"Q73", f}]*DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4] - 2*Condensate[{"Q72", f}]*(D*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2] - 
       DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] + DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] - DiracSigma[GAD[lor1], GAD[lor4]]*MTD[lor2, lor3] + 
       DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4]) + (-2 + D)*Condensate[{"Q71", f}]*(2*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2] + 
       DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] - DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] + DiracSigma[GAD[lor1], GAD[lor4]]*MTD[lor2, lor3] + 
       (2*I)*MTD[lor1, lor4]*MTD[lor2, lor3] + I*D*MTD[lor1, lor4]*MTD[lor2, lor3] - DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4] - (2*I)*MTD[lor1, lor3]*MTD[lor2, lor4] - 
       I*D*MTD[lor1, lor3]*MTD[lor2, lor4]) + 4*Condensate[{"Q6", f}]*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2]*quarkMass[f] - 
     D*Condensate[{"Q6", f}]*DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3]*quarkMass[f] + D*Condensate[{"Q6", f}]*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4]*
      quarkMass[f] - D*Condensate[{"Q6", f}]*DiracSigma[GAD[lor1], GAD[lor4]]*MTD[lor2, lor3]*quarkMass[f] + D*Condensate[{"Q6", f}]*DiracSigma[GAD[lor1], GAD[lor3]]*
      MTD[lor2, lor4]*quarkMass[f] - 2*D*Condensate[{"Q5", f}]*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2]*quarkMass[f]^2 + 
     2*Condensate[{"Q5", f}]*DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3]*quarkMass[f]^2 - 2*Condensate[{"Q5", f}]*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4]*
      quarkMass[f]^2 + 2*Condensate[{"Q5", f}]*DiracSigma[GAD[lor1], GAD[lor4]]*MTD[lor2, lor3]*quarkMass[f]^2 - 2*Condensate[{"Q5", f}]*DiracSigma[GAD[lor1], GAD[lor3]]*
      MTD[lor2, lor4]*quarkMass[f]^2))*SUNT[coln])/(8*CA*CF*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D));

        
If[OptionValue[Factorization],
   cond=Factorization[cond]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
If[OptionValue[Explicit]===True,
   cond=cond/.{AGammaD[lors__]:>AGammaD[lors,Explicit->True],AGamma[lors__]:>AGamma[lors,Explicit->True]}
];

cond
]


(* \[LeftAngleBracket]q_j(<-D)^a (DG)^n (Dq_j)^b\[RightAngleBracket] *)
D7Condensate[lor1_,{coln_,{lor2_,lor3_}},lor4_,f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=(((-2*I)*(2 + D)*AGammaD[lor1, lor2, lor3, lor4]*Condensate[{"Q72", f}] + (-3 + D)*(2*Condensate[{"Q74", f}]*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2] - 
     D*Condensate[{"Q74", f}]*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2] - 2*Condensate[{"Q74", f}]*DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] + 
     D*Condensate[{"Q74", f}]*DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] - 4*Condensate[{"Q74", f}]*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] + 
     2*D*Condensate[{"Q74", f}]*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] - 2*Condensate[{"Q74", f}]*DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4] + 
     D*Condensate[{"Q74", f}]*DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4] - (4*I)*Condensate[{"Q74", f}]*MTD[lor1, lor3]*MTD[lor2, lor4] + 
     I*D^2*Condensate[{"Q74", f}]*MTD[lor1, lor3]*MTD[lor2, lor4] + 2*Condensate[{"Q74", f}]*DiracSigma[GAD[lor1], GAD[lor2]]*MTD[lor3, lor4] - 
     D*Condensate[{"Q74", f}]*DiracSigma[GAD[lor1], GAD[lor2]]*MTD[lor3, lor4] + (4*I)*Condensate[{"Q74", f}]*MTD[lor1, lor2]*MTD[lor3, lor4] - 
     I*D^2*Condensate[{"Q74", f}]*MTD[lor1, lor2]*MTD[lor3, lor4] - 2*Condensate[{"Q72", f}]*(DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2] - 
       DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] + D*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] - DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4] + 
       DiracSigma[GAD[lor1], GAD[lor2]]*MTD[lor3, lor4]) + (-2 + D)*Condensate[{"Q71", f}]*(-(DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2]) + 
       DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3] + 2*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4] + DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4] + 
       (2*I)*MTD[lor1, lor3]*MTD[lor2, lor4] + I*D*MTD[lor1, lor3]*MTD[lor2, lor4] - DiracSigma[GAD[lor1], GAD[lor2]]*MTD[lor3, lor4] - (2*I)*MTD[lor1, lor2]*MTD[lor3, lor4] - 
       I*D*MTD[lor1, lor2]*MTD[lor3, lor4]) + D*Condensate[{"Q6", f}]*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2]*quarkMass[f] - 
     D*Condensate[{"Q6", f}]*DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3]*quarkMass[f] + 4*Condensate[{"Q6", f}]*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4]*
      quarkMass[f] - D*Condensate[{"Q6", f}]*DiracSigma[GAD[lor1], GAD[lor3]]*MTD[lor2, lor4]*quarkMass[f] + (8*I)*Condensate[{"Q6", f}]*MTD[lor1, lor3]*MTD[lor2, lor4]*
      quarkMass[f] - (2*I)*D^2*Condensate[{"Q6", f}]*MTD[lor1, lor3]*MTD[lor2, lor4]*quarkMass[f] + D*Condensate[{"Q6", f}]*DiracSigma[GAD[lor1], GAD[lor2]]*MTD[lor3, lor4]*
      quarkMass[f] - (8*I)*Condensate[{"Q6", f}]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] + (2*I)*D^2*Condensate[{"Q6", f}]*MTD[lor1, lor2]*MTD[lor3, lor4]*quarkMass[f] - 
     2*Condensate[{"Q5", f}]*DiracSigma[GAD[lor3], GAD[lor4]]*MTD[lor1, lor2]*quarkMass[f]^2 + 2*Condensate[{"Q5", f}]*DiracSigma[GAD[lor2], GAD[lor4]]*MTD[lor1, lor3]*
      quarkMass[f]^2 - 2*D*Condensate[{"Q5", f}]*DiracSigma[GAD[lor2], GAD[lor3]]*MTD[lor1, lor4]*quarkMass[f]^2 + 2*Condensate[{"Q5", f}]*DiracSigma[GAD[lor1], GAD[lor3]]*
      MTD[lor2, lor4]*quarkMass[f]^2 - 2*Condensate[{"Q5", f}]*DiracSigma[GAD[lor1], GAD[lor2]]*MTD[lor3, lor4]*quarkMass[f]^2))*SUNT[coln])/
 (8*CA*CF*(-3 + D)*(-2 + D)*(-1 + D)*D*(2 + D));

        
If[OptionValue[Factorization],
   cond=Factorization[cond]
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
