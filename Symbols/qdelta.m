(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



qdelta::usage = "Delta symbol used in FourierXP, FourierPX, QNormal, and QEvaluate."
Begin["`Private`qdelta`"]


qdelta/: MakeBoxes[qdelta,TraditionalForm]="\[Delta]"


End[]
