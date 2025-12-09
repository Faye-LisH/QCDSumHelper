(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



D5Condensate::usage = "D5Condensate[lor1_,lor2_,psi_] = \[LeftAngleBracket]\!\(\*OverscriptBox[\(q\), \(_\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor1\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor2\)]\)\!\(\*SuperscriptBox[\(\[Del]\), \(lor3\)]\)q\[RightAngleBracket]"


Begin["`Private`D5Condensate`"]
Options[D5Condensate] = {
	Explicit->False,
	Massless->True,
	Factorization->False
}


(* \[LeftAngleBracket]q^a_i (DDq_j)^b\[RightAngleBracket]*)
D5Condensate[lor1:Except[_List],lor2:Except[_List],f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=(1/(8 D CA)(MTD[lor1,lor2]+I/(D-1) DiracSigma[GAD[lor1,lor2]])Condensate[{"Q5",f}]
	-quarkMass[f]^2 Condensate[{"Q3",f}]MTD[lor1,lor2]/(4 D CA));
  
              
If[OptionValue[Factorization],
   cond=Factorization[cond]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
 

cond
]


(* \[LeftAngleBracket]q^a_i G^n q_j^b\[RightAngleBracket] *)
D5Condensate[{coln_,{lor1_,lor2_}},f:Except[_Rule],OptionsPattern[]]:=Block[{cond},

cond=1/(4(D-1)D CA CF)SUNT[coln]Condensate[{f,"G",f}]DiracSigma[GAD[lor1,lor2]];
        
If[OptionValue[Factorization],
   cond=Factorization[cond]
];
    
If[OptionValue[Massless]===True,
   cond=cond/._quarkMass->0
];
   

cond
]


End[]
