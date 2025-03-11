(* ::Package:: *)

(* Wolfram Language package *)

(* Author: ShungHong Li *)



NonLocalCondensate::usage="Combine several fields into a nonlocal condensate"


Begin["`Private`NonLocalCondensateCondensate`"]	


NonLocalCondensate/:MakeBoxes[NonLocalCondensate[fields_List],TraditionalForm]:=RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",UnderscriptBox[RowBox[{##}],"_"],"\[ThinSpace]","\[RightAngleBracket]"}]&@@(ToBoxes[#,TraditionalForm]&/@fields)
NonLocalCondensate/:MakeBoxes[NonLocalCondensate[0,fields_List],TraditionalForm]:=RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",UnderscriptBox[RowBox[{##}],"_"],"\[ThinSpace]","\[RightAngleBracket]"}]&@@(ToBoxes[#,TraditionalForm]&/@fields)
NonLocalCondensate/:MakeBoxes[NonLocalCondensate[1,fields_List],TraditionalForm]:=RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",UnderscriptBox[RowBox[{##}],"_"],SubscriptBox["|","m"],"\[RightAngleBracket]"}]&@@(ToBoxes[#,TraditionalForm]&/@fields)
NonLocalCondensate/:MakeBoxes[NonLocalCondensate[m_Integer,fields_List],TraditionalForm]:=RowBox[{"\[LeftAngleBracket]","\[ThinSpace]",UnderscriptBox[RowBox[{##}],"_"],SubscriptBox["|",SuperscriptBox["m",ToBoxes[m]]],"\[RightAngleBracket]"}]&@@(ToBoxes[#,TraditionalForm]&/@fields)/;m>1


NonLocalCondensate/:MakeBoxes[NonLocalCondensate[fields_List,cond_List],TraditionalForm]:=RowBox[Join[{"\[LeftAngleBracket]","\[ThinSpace]",UnderscriptBox[RowBox[ToBoxes[#,TraditionalForm]&/@fields],"_"],"\[ThinSpace]","\[RightAngleBracket]","\[Rule]"},ToBoxes[#,TraditionalForm]&/@cond]]
NonLocalCondensate/:MakeBoxes[NonLocalCondensate[0,fields_List,cond_List],TraditionalForm]:=RowBox[Join[{"\[LeftAngleBracket]","\[ThinSpace]",UnderscriptBox[RowBox[ToBoxes[#,TraditionalForm]&/@fields],"_"],"\[ThinSpace]","\[RightAngleBracket]","\[Rule]"},ToBoxes[#,TraditionalForm]&/@cond]]
NonLocalCondensate/:MakeBoxes[NonLocalCondensate[1,fields_List,cond_List],TraditionalForm]:=RowBox[Join[{"\[LeftAngleBracket]","\[ThinSpace]",UnderscriptBox[RowBox[ToBoxes[#,TraditionalForm]&/@fields],"_"],SubscriptBox["|","m"],"\[RightAngleBracket]","\[Rule]"},ToBoxes[#,TraditionalForm]&/@cond]]
NonLocalCondensate/:MakeBoxes[NonLocalCondensate[m_Integer,fields_List,cond_List],TraditionalForm]:=RowBox[Join[{"\[LeftAngleBracket]","\[ThinSpace]",UnderscriptBox[RowBox[ToBoxes[#,TraditionalForm]&/@fields],"_"],SubscriptBox["|",SuperscriptBox["m",ToBoxes[m]]],"\[RightAngleBracket]","\[Rule]"},ToBoxes[#,TraditionalForm]&/@cond]]/;m>1


End[]
