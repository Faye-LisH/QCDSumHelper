(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



ColorDelta::usage =
	"ColorDelta[i_,j_] is a Kroncker-delta for SU(n) with the indices i, j in basic represention."	

Begin["`Private`ColorDelta`"]	
(*Options[ColorDelta] = {}*)



SetAttributes[{ColorDelta}, Orderless];

ColorDelta/: MakeBoxes[ColorDelta[i_,j_],TraditionalForm]:=SuperscriptBox[ToBoxes["\[Delta]",TraditionalForm],RowBox[{ToBoxes[i,TraditionalForm],ToBoxes[j,TraditionalForm]}]]




ColorDelta/:ColorDelta[i_,j_]ColorDelta[j_,k_]:=ColorDelta[i,k]
ColorDelta/:Power[ColorDelta[i_,j_],2]=SUNN
(*ColorDelta/:Power[others_ ColorDelta[i_,j_],2]:=SUNN Power[others,2]*)

ColorDelta/:SUNTF[ca_,i_,j_]ColorDelta[j_,k_]:=SUNTF[ca,i,k]
ColorDelta/:SUNTF[ca_,i_,j_]ColorDelta[i_,k_]:=SUNTF[ca,k,j]


ColorDelta/:SUNTF[ca_List,SUNFIndex[i_],SUNFIndex[j_]]ColorDelta[j_,k_]:=SUNTF[ca,i,k]
ColorDelta/:SUNTF[ca_List,SUNFIndex[i_],SUNFIndex[j_]]ColorDelta[i_,k_]:=SUNTF[ca,k,j]

ColorDelta/:SUNTF[ca_,i_,j_]ColorDelta[i_,j_]:=0
ColorDelta/:SUNTF[ca_,j_,i_]ColorDelta[i_,j_]:=0


ColorDelta/:SUNTF[ca_List,SUNFIndex[i_],SUNFIndex[j_]]ColorDelta[i_,j_]:=SUNTrace[Dot@@(SUNT[#]&/@ca)]
ColorDelta/:SUNTF[ca_List,SUNFIndex[j_],SUNFIndex[i_]]ColorDelta[i_,j_]:=SUNTrace[Dot@@(SUNT[#]&/@ca)]


ColorDelta/:SUNFDelta[SUNFIndex[i_],SUNFIndex[j_]]ColorDelta[j_,k_]:=ColorDelta[i,k]
ColorDelta/:SUNFDelta[SUNFIndex[j_],SUNFIndex[i_]]ColorDelta[j_,k_]:=ColorDelta[i,k]


ColorDelta[i_,i_]=CA



End[]
