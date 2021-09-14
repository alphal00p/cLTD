(* ::Package:: *)

(* ::Section:: *)
(*Evaluate cLTD Expression*)


(* ::Subsection::Closed:: *)
(*Path*)


$CLTDPATH = NotebookDirectory[];


(* ::Subsection::Closed:: *)
(*Functions*)


(* ::Subsubsection::Closed:: *)
(*Numerator*)


ClearAll[NumVarUnfold]
NumVarUnfold[NumFun_,y_List,yvar_]:=Module[{},
	If[Length[y] == 1,
		NumFun/.yvar->y[[1]],
		-(NumVarUnfold[NumFun,Drop[y,{-1}],yvar]-NumVarUnfold[NumFun,Drop[y,{-2}],yvar])/Apply[Subtract,y[[-2;;]]]]//Cancel
]

ClearAll[num]
num[ys__List,Num_,vars_List]:=Module[{xs,x,vi,numtmp},
	(*Check inputs*)
	If[Length[{ys}]!=Length[vars],Print["ys list must be of same length as vars"];Abort[]];
	
	(*Top numerator N0 initialized with the provided Num function*)
	numtmp=Num@@vars;
	(*Apply steps*)
	Table[
		xs = vars[[i+1;;]]/.List->Sequence;
		numtmp=NumVarUnfold[numtmp,{ys}[[i]],vars[[i]]];
		(*num[i][xs]*),
	{i,Length[{ys}]}];
	(*Return Final expression*)
	numtmp
]


(*MyFun[x_,y_]:=x^2
MyFun[x_,y_]:=f[x,y]
num[{x1+y,x2,x3},{1},MyFun,{x,y}]
num[{x1,x2},{1},MyFun,{x,y}]*)


(* ::Subsubsection::Closed:: *)
(*Evaluate Element*)


ClearAll[evaluate]
evaluate[{factor_,dens_List,zs_List},NumFunction_,vars_]:=Module[{
		denominator = Times@@dens,
		numerator=num @@Join[zs,{NumFunction,vars}]
	},
	prefactor[dens,zs]factor *numerator/denominator
]

ClearAll[prefactor]
prefactor[dens_List,zs_List]:=1/Product[2 ToExpression["E"<>ToString[n]],{n,0,Length[Select[dens,#=!=1&]]+Apply[Plus,Map[Length[#]&,zs]]-1}];


(* ::Subsection:: *)
(*Import and Evaluate*)


(* ::Subsubsection::Closed:: *)
(*Box*)


ClearAll[pf,MyNum]
(*Import instructions*)
pf= Import[$CLTDPATH <> "cLTD_1l_box.m"]/.a_Real:>Rationalize[a];
(*Define Numerator and Evaluate for different powers*)
MyNum[x_]:=x^0; pf0=evaluate[#,MyNum,{k0}]&/@pf//Simplify//Cancel;
MyNum[x_]:=x^1; pf1=evaluate[#,MyNum,{k0}]&/@pf//Simplify//Cancel;
MyNum[x_]:=x^2; pf2=evaluate[#,MyNum,{k0}]&/@pf//Simplify//Cancel;
MyNum[x_]:=x^3; pf3=evaluate[#,MyNum,{k0}]&/@pf//Simplify//Cancel;


(*Compute Residue Explicitly*)
ConvertNotation =  {e[n_]:>ToExpression["E"<>ToString[n]],p[n_]:>ToExpression["p"<>ToString[n]]};
rExp=1/Product[k+e[i]+p[i],{i,0,3}]/Product[k-e[i]+p[i],{i,0,3}]/.ConvertNotation;
resPositions = Table[-e[i]-p[i],{i,0,3}]/.ConvertNotation;
res0=Table[Residue[k^0 rExp,{k,r}],{r,resPositions}];
res1=Table[Residue[k^1 rExp,{k,r}],{r,resPositions}];
res2=Table[Residue[k^2 rExp,{k,r}],{r,resPositions}];
res3=Table[Residue[k^3 rExp,{k,r}],{r,resPositions}];


(*Compare with PF expression*)
Cancel[Apply[Plus,res0]/Apply[Plus,pf0]]
Cancel[Apply[Plus,res1]/Apply[Plus,pf1]]
Cancel[Apply[Plus,res2]/Apply[Plus,pf2]]
Cancel[Apply[Plus,res3]/Apply[Plus,pf3]]


(*Compare with degenerate edges*)
MyNum[x_]:= Array[c,4,0].Array[x^#&,4,0]
vacuumBox=Residue[MyNum[k]/(k-E0+p0)^4/(k+E0+p0)^4,{k,-E0-p0}];
PFvacuumBox=Apply[Plus,evaluate[#,MyNum,{k0}]&/@pf/.{E0|E1|E2|E3->E0,p0|p1|p2|p3->p0}]//Simplify;
Cancel[vacuumBox/PFvacuumBox]


(* ::Subsubsection::Closed:: *)
(*Sunrise*)


ClearAll[pf,MyNum]
(*Import instructions*)
pf= Import[$CLTDPATH <> "cLTD_2l_sunrise.m"]/.a_Real:>Rationalize[a];
(*Define Numerator*)
MyNum[x_,y_]:=c[1]x+c[2]y^2x^2
(*MyNum[x_,y_]:=1*)
(*Evaluate*)
Monitor[pf2=Table[evaluate[pf[[ipf]], MyNum, {k0, k1}], {ipf, Length[pf]}], PercentForm[N[ipf/Length[pf]]]]


(* ::Input:: *)
(**)


(* ::Subsubsection::Closed:: *)
(*PentaBox*)


ClearAll[pf, MyNum]
pf = Import[$CLTDPATH <> "cLTD_2l_pentabox.m"]/.a_Real:>Rationalize[a]; 
MyNum[x_, y_] := x^2y^2
(*MyNum[x_, y_] := 1*)
(*Evaluate*)
Monitor[pf2=Table[evaluate[pf[[ipf]], MyNum, {k0, k1}], {ipf, Length[pf]}], PercentForm[N[ipf/Length[pf]]]]


blocks=Select[Chop[pf2]/.a_Real:>Round[a],#=!=0&]//Cancel
dens=Union[Flatten[Chop[pf[[All,2]]]/.a_Real:>Round[a]]][[2;;]]//Simplify;
subDen=Table[dens[[i]]->Den[i],{i,Length[dens]}];


Plus@@(blocks/.subDen)


(* ::Subsubsection::Closed:: *)
(*2x2 Fishnet*)


ClearAll[pf,MyNum]
(*Import instructions*)
pf= Import[$CLTDPATH <> "cLTD_4l_2x2_fishnet.m"]/.a_Real:>Rationalize[a];
(*Define Numerator*)
MyNum[x1_,x2_,x3_,x4_]:=1
(*Evaluate*)
Print["Evaluate:"]
Monitor[ipf=0;pf2=Table[evaluate[pf[[ipf]],MyNum,{k0,k1,k2,k3}],{ipf,Length[pf]}],PercentForm[N[ipf/Length[pf]]]]
