(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



D6Condensate::usage = "D6Condensate[lor1_,lor2_,lor3_,psi_] = \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor1\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor2\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor3\)]\)q\[RightAngleBracket]"


Begin["`Private`D6Condensate`"]
Options[D6Condensate] = {
	Explicit->False,
	Massless->True,
	Factorization->False
}


(* \[LeftAngleBracket]q^a_i (DDDq_j)^b\[RightAngleBracket]*)
D6Condensate[lor1:Except[_List],lor2:Except[_List],lor3:Except[_Rule],f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=(((I/8)*Condensate[{"Q6", f}]*((2 + D)*AGammaD[lor1, lor2, lor3] + (-2 + D)*(GAD[lor3]*MTD[lor1, lor2] - (1 + D)*GAD[lor2]*MTD[lor1, lor3] + GAD[lor1]*MTD[lor2, lor3])))/(CA*(-2 + D)*(-1 + D)*D*(2 + D)) - 
 ((I/8)*Condensate[{"Q5", f}]*((2 + D)*AGammaD[lor1, lor2, lor3] + (2 - 3*D + D^2)*(GAD[lor3]*MTD[lor1, lor2] + GAD[lor2]*MTD[lor1, lor3] + GAD[lor1]*MTD[lor2, lor3]))*quarkMass[f])/
  (CA*(-2 + D)*(-1 + D)*D*(2 + D)) + ((I/4)*Condensate[{"Q3", f}]*(GAD[lor3]*MTD[lor1, lor2] + GAD[lor2]*MTD[lor1, lor3] + GAD[lor1]*MTD[lor2, lor3])*quarkMass[f]^3)/(CA*D*(2 + D)));

        
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


(* \[LeftAngleBracket]q^a_i G^n (Dq_j)^b\[RightAngleBracket] *)
D6Condensate[{coln_,{lor1_,lor2_}},lor3:Except[_Rule],f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=(I*(((I/8)*Condensate[{"Q6", f}]*(2*AGammaD[lor1, lor2, lor3] - (-2 + D)*(GAD[lor2]*MTD[lor1, lor3] - GAD[lor1]*MTD[lor2, lor3])))/((-2 + D)*(-1 + D)*D) - 
   ((I/4)*AGammaD[lor1, lor2, lor3]*Condensate[{"Q5", f}]*quarkMass[f])/(D*(2 - 3*D + D^2)))*SUNT[coln])/(CA*CF);
        
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


(* \[LeftAngleBracket]q_i(<-D)^a G^n q^b_j\[RightAngleBracket] *)
D6Condensate[lor1:Except[_List],{coln_,{lor2_,lor3_}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=((-I)*(((I/8)*Condensate[{"Q6", f}]*(2*AGammaD[lor1, lor2, lor3] + (-2 + D)*(GAD[lor3]*MTD[lor1, lor2] - GAD[lor2]*MTD[lor1, lor3])))/((-2 + D)*(-1 + D)*D) - 
   ((I/4)*AGammaD[lor1, lor2, lor3]*Condensate[{"Q5", f}]*quarkMass[f])/(D*(2 - 3*D + D^2)))*SUNT[coln])/(CA*CF);
        
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


(* FCE[FCI[D6Condensate[lor1,{coln,{lor2,lor3}},f]]/.gg_DiracGamma:>-gg]/.{lor1->lor3,lor3->lor1}
D6Condensate[{coln,{lor1,lor2}},lor3,f]-%//Simplify*)


(* \[LeftAngleBracket]q^a_i (DG)^n q^b_j\[RightAngleBracket] *)
D6Condensate[{coln_,{lor1_,{lor2_,lor3_}}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=-(Condensate[{"Q6",f}]*(GAD[lor3]*MTD[lor1,lor2]-GAD[lor2]*MTD[lor1,lor3])*SUNT[coln])/(4*CA*CF*(-1+D)*D);

        
If[OptionValue[Factorization],
   cond=Factorization[cond]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   
(* no AGammaD involved *)

cond
]


End[]
