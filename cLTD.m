(* ::Package:: *)

(* ::Chapter:: *)
(*cLTD package*)


BeginPackage["cLTD`"];


(* ::Section:: *)
(*Global *)


(* ::Subsection::Closed:: *)
(*Usage*)


cLTD::usage="\
Compute the corresponding cLTD expression using FORM.
"<>Style["Input Format",Bold]~ToString~StandardForm<>":
The scalar propagators of the expression must be expressed as: 
	-> prop[4momentum,mass] (e.g. prop[k1+p1, 0])

Scalar products among 4-vectors in the numerator can be expressed as:
    -> SP4[p1,p2]

The numerator can also be a polynomial in the loop momenta's energy 
components
    -> e.g.: k[0]*p1[0] - p1[0]*p2[0] (k: loop momentum) 

"<>Style["Options",Bold]~ToString~StandardForm<>": 
	To change change the default options to different values use SetOptions[cLTD,...]:
	e.g.:
		Use `form` in the alphaloop directory "<>Style["SetOptions[cLTD,\"FORMpath\"\[Rule]\"/home/dude/dir1/dir2/form\"]",Bold]~ToString~StandardForm<>"

"<>Style["Example",Bold]~ToString~StandardForm<>":
 expr = SP4[k1,p1] prop[k1-p1,0]prop[k1-p2,0] - k1[0] p1[0] prop[k1-p3,m]prop[k1-p4,m]
=> cLTD[expr]

WARNING: when summing multiple diagrams, they must contain the same 
         number of loops, if this is not possible, then call cLTD[] 
         for each of them individually.

"<>Style["Output",Bold]~ToString~StandardForm<>":
  { cLTD expression, Energies definition} 
  The expression contains the functions den[] and cLTDnorm[] which are to be evaluated as 
		den[a_] :> 1/a, cLTDnorm[a_] :> 1/a
";

Print[" :::::::::::::::::::::::: cLTD ::::::::::::::::::::::::"]
Print["Authors: Z. Capatti, V. Hirschi, D. Kermanschah, A. Pelloni, B. Ruijl"];
Print["\tA Mathematica front end for cLTD [arxiv:2009.05509].\n"];
(*Print[cLTD::usage];*)


(* ::Section:: *)
(*Private*)


Begin["cLTDPrivate`"];


(* ::Subsection::Closed:: *)
(*FORM file*)


FORMhead = "Format 100;
Auto S n,m,y,E,p,k;
S iComplex;
CF ncmd, num, Eres, Echain, Diff;
CF den, prop, prop0, norm, error;
CF x, xbar;

set energies:E0,...,E1000;
set shifts:p0,...,p1000;
";


(* ::Subsubsection::Closed:: *)
(*cLTD*)


FORMpf[filename_,OptimizationLVL_] := "multiply norm(1)*num();
repeat id norm(E?)*prop(y0?,y1?) = norm(2*E*y1)*(den(y0-y1)-den(y0+y1)); 
.sort

#ifdef `NoNumerator'
	#define NumeratorSwitch \"0\"
#else
	#define NumeratorSwitch \"1\"
#endif
* NOTE:
*    This procedure follows what presented in the paper 2009.05509
*    The only difference is that the numerator function is that the
*    divided difference is definded with a minus sing.
*    As a result we have that each time a ncmd gain an extra term 
*    we need to compensate with an additional minus sing.
#do i=0,{`LOOPS'-1}
* Normalise by the contour orientation 
    id norm(n?) = norm(-n);
  
* Activate numerator
    id num(?y) = num(`NumeratorSwitch',?y,ncmd());

* Identify poles in the loop momentum energy 
    splitarg (k`i') den;
    id den(?a,p?) = den(p,?a);
    factarg (-1) den 1;
    id den(p?,y?,p1?) = x(p1/y,0,0)/y;
    id den(p?,y?) = den(y*p);
    .sort:get-poles;

* Split the poles by their position in the complex plane
    splitarg x; 
    .sort
    id x(?y1,+E?energies,?y2,E1?,E2?) = x(?y1,?y2,E1+E,E2);
    id x(?y1,-E?energies,?y2,E1?,E2?) = x(?y1,?y2,E1,E2+E);

* check for errors
    id x(?y,0,0) = error(x(?y,0,0));
    id x(?y,E1?!{,0},E2?!{,0}) = error(x(?y, E1,E2));
    if (count(error,1));
         print \"[k`i'] Cannot find energy signature %t\";
         exit \"Critical ERROR\";
    endif;

* re-absorb 
    id x(?y,0,E?!{,0}) = x(?y,-E);
    id x(?y,E?!{,0},0) = xbar(?y,E);
    transform, x,xbar , addargs(1,last);
    id x(y?) = x(-y);
    id xbar(y?) = xbar(-y);
*    print +s;

* Remove elements with zero residue
    if (count(x,1) == 0) discard;
    if ( count(xbar,1) == 0 ) 
        repeat id num(?n,ncmd(?n1))*x(y1?) = -num(?n,ncmd(?n1,y1));
    .sort:x-xbar;

    
* Start unfolding the expression 
    repeat;
        if ( count(x,1) == 1 );
            repeat id x(y1?)*xbar(y2?) = x(y1)*den(y1-y2);
            id num(0,?n,ncmd(?n1))*x(y1?) = num(?n,ncmd(?n1));
            id num(1,?n,ncmd(?n1))*x(y1?) = -num(?n,ncmd(?n1,y1));
        endif;


        id once num(0,?n,ncmd(?n1))*x(y1?)*xbar(y2?) = 
                     -num(0,?n,ncmd(?n1))*Eres(y1,y2)*xbar(y2);
    
        id once num(1,?n,ncmd(?n1))*x(y1?)*xbar(y2?) = 
                     +num(0,?n,ncmd(?n1,y1))*Eres(y1,y2)*xbar(y2)
                     -num(1,?n,ncmd(?n1,y1))*xbar(y2);

        repeat id Eres(y1?,y2?)*xbar(y2?)*xbar(y3?)= 
                      den(y1-y2)*xbar(y3)*(xbar(y2)+Eres(y1,y3));
        id once Eres(y1?,y2?)*xbar(y2?)= den(y1-y2)*xbar(y2);
    endrepeat;
    id num(n0?{0,1},?n) = num(?n);
    .sort:recursive-expansion;

* Sanity check
    if (count(xbar,1) || count(x,1));
         print \"[k`i'] Some roots were not absorbed: %t\";
         exit \"Critical ERROR\";
    endif;
#enddo

* If there is no numerator call it a day
#ifdef `NoNumerator'
    id num(?x1, ncmd(n1?,n2?,?n3),?x2) = 0;
	id num(?x1) = 1;
	.sort
	Format mathematica;
	Format O"<>ToString[OptimizationLVL]<>", stats=on, method=greedy;
    ExtraSymbols, array, FormEvalStep;
	#Optimize F;
	#write<"<>filename<>".out> \"%O\";
	#write<"<>filename<>".out> \"%E\" F;
	.end
#endif

";


FORMnum[nLoops_] := "*off statistics;
CF ncmd, kval, num, a;


Auto S invd;
Auto S f,z,r,x,s,c;

#$LOOPS="<>ToString[nLoops]<>";

** Additional numerator function in LMB
*#$numF = 1;
** Multiply numerator functions
*Multiply $numF;
*print;
.sort:start;

* Start unfoding the numerator instructions loop by loop
#do i = 0,{`$LOOPS'-1}
    id once num(ncmd(?z),?x) =ncmd(?z,0)*num(?x);
    B+ ncmd, k`i';
    .sort:collect-ncmd;
    keep brackets;
    
    id k`i'^r?*ncmd(?z,z1?,0) = sum_(s,0, r - nargs_(?z) + 0, a(r,s,?z)*z1^s);
    B+ a;
    .sort:collect-a;
    keep brackets;

    repeat id a(r?,s?,?z,z1?) = -sum_(k, s+1,r-nargs_(?z) + 0, a(r,k,?z)*z1^(k-s-1));
    id a(r?,s?) = delta_(r,s);
    .sort:energy-k`i';
#enddo
id num=1;

* check that the substitution is complete
if (count(ncmd, 1,num,1,a, 1));
    Print \"Unsubstituted ncmd: %t\";
    exit \"Critical error\";
endif;
.sort:rm-num;
";


(* ::Subsubsection:: *)
(*LTD*)


FORMltd[filename_] := "multiply norm(1)*num();

* NOTE:
*   Compute the LTD expression for arbitrary topologies 
id prop(?a) = prop(?a,1);
#do i=0,{`LOOPS'-1}
    B+ prop k`i' norm num;
    .sort:loop-{`i'+1};
    keep brackets;
* Deal with powers 
    repeat id prop(?a,n1?)*prop(?a,n2?) = prop(?a,n1+n2);

* Identify poles in the loop momentum energy 
    splitarg (k`i') prop;
    id prop(k`i',E?,n?) = prop0(k`i',0,E,n); 
    id prop(p?,k?,E?,n?) = prop0(k,p,E,n); 

    factarg (-1) prop0 1;
    id prop0(k`i',y?{-1,1},E?,n?) = prop0(k`i',y,0,E,n); 
    id all prop0(k?,y?{-1,1},p?,E?,n?)*norm(n0?) = 
        sum_(m,0,n-1,
            (-1)^(n-m-1)*fac_(2*n-m-2)/fac_(n-m-1)/fac_(m)*
            Diff(k`i',m)*norm(n0*(2*E)^(2*n-m-2))
            )/fac_(n-1)*
        x(y*p,0,E)*x(y*p,E,0); 

* Apply differential
    repeat;
         repeat id Diff(k`i',n1?{>0})*prop0(k`i',y?{-1,1},p?,E?,n2?) =  
             +Diff(k`i',n1)*prop(k`i',y,p,E,n2)
             -Diff(k`i',n1-1,0)*(2*n2*(k`i'+y*p))*prop0(k`i',y,p,E,n2+1) 
         ; 
         id Diff(k`i',n1?{>0})*k`i'^n3? =
             +n3*Diff(k`i',n1-1,0)*k`i'^(n3-1)
             ; 
         id Diff(k`i',n?{>0})= 0;
         id Diff(k`i',n?,0)=Diff(k`i',n);
         id prop(k`i',y?{-1,1},p?,E?,n2?) = prop0(k`i',y,p,E,n2);
    endrepeat;
    id Diff(k`i',0)=1;
    id prop0(k?,y?{-1,1},?a) = prop(y*k,?a); 

* check for errors
    if (count(Diff,1));
         print \"Remaining differenctial operators: %t\";
         exit \"Critical ERROR\";
    endif;

* Split the poles by their position in the complex plane
    splitarg x;
    repeat id x(?y1,+E?energies,?y2,E1?,E2?) = x(?y1,?y2,E1+E,E2);
    repeat id x(?y1,-E?energies,?y2,E1?,E2?) = x(?y1,?y2,E1,E2+E);
    id x(E1?,E2?) = x(0,E1,E2);

* check for errors
    id x(?y,0,0) = error(x(?y,0,0));
    if (count(error,1));
         print \"[k`i'] Cannot find energy signature %t\";
         exit \"Critical ERROR\";
    endif;
    transform, x, prop, addargs(1,last-2);

* take residue
    id x(p?,0,E?)*x(p?,E1?,E2?)*norm(n?)*num(?y) = replace_(k`i',E-p)*norm(n*(E+E1-E2))*num(?y,E-p);
    if ( count(x,1) ) discard;
#enddo
id prop(?a,n?) = prop(?a)^n;
id prop(p?,E?) = den(p^2-E^2);
id num(?a) = 1;

";


(* ::Subsubsection::Closed:: *)
(*Generate FORM file*)


FORMfile[filename_,expr_,nLoops_,FORMvars_,OptimizationLVL_,stdLTD_:False]:=Module[{
file="#--\n"},
file = file <> FORMhead;
file = file <> "Auto S sp,"<>StringRiffle[ToString/@FORMvars,","]<>";\n";
file = file <> "#define LOOPS \""<>ToString[nLoops]<>"\"\n";
file = file <> "#define topoID \""<>filename<>"\"\n";
file = file <> "L F = "<>expr<>";\n";
If[stdLTD, 
	file = file <> FORMltd[filename] <> "\n.sort\n";	
,
	file = file <> FORMpf[filename,OptimizationLVL] <> "\n.sort\n";
	file = file <> "L firstF = firstterm_(F);\n.sort\n#if ( termsin(firstF) != 0 )\n";
	file = file <> "#redefine oldextrasymbols \"`extrasymbols_'\"
ExtraSymbols, underscore, den;
hide firstF;
argtoextrasymbol den;
id den(y?) = y;
.sort:den;\n#redefine newextrasymbols \"`extrasymbols_'\"\noff statistics;\n";
	file = file <> FORMnum[nLoops]; 
	file = file <> "#if `oldextrasymbols' != `newextrasymbols'\n";
	file = file <> "multiply replace_(<den{`oldextrasymbols'+1}_,den(extrasymbol_({`oldextrasymbols'+1}))>\\
                  ,...,<den`extrasymbols_'_,den(extrasymbol_(`extrasymbols_'))>);\n";
	file = file <> "#endif\n";
	file = file <> "#endif\n.sort";
];
	file = file <> "
Format mathematica;
Format O"<>ToString[OptimizationLVL]<>", stats=on, method=greedy;
ExtraSymbols, array, FormEvalStep;
#Optimize F;
#write<"<>filename<>".out> \"%O\";
#write<"<>filename<>".out> \"%E\" F;\n.end";
file
]


(* ::Subsection::Closed:: *)
(*Functions*)


toLTDprop[loopmom_]:=Module[{SP},
	SetAttributes[SP,Orderless];
	SP[p1_+p2_,p3_]:=SP[p1,p3]+SP[p2,p3];
	SP[-p1_,p2_]:=-SP[p1,p2];
	Global`prop[mom_?(!FreeQ[#,Alternatives@@loopmom]&),m_]:>LTDprop[
		mom/.Table[v->v[0],{v,Complement[Variables[mom],loopmom]}],
		Sqrt[SP[mom,mom]+m^2]/.SP->Dot]
	];


SetAttributes[SP4,Orderless];
SP4[p1_+p2_,q_][loopmom_]:=SP4[p1,q][loopmom]+SP4[p2,q][loopmom]
SP4[a_Integer p_,q_][loopmom_]:=a SP4[p,q][loopmom]
SP4[p_,q_][loopmom_]:=If[FreeQ[loopmom,p],p[0],p]*If[FreeQ[loopmom,q],q[0],q] - p . q 


Options[cLTD]={
	"loopmom"->{Global`k0,Global`k1,Global`k2,Global`k3},
	"FORMpath"->"form",
	"tFORMpath"->"tform",
	"WorkingDirectory"->Check[NotebookDirectory[],Print["Cannot find Notebook Directory!\nMust define the WorkingDirectory option before running."]; None] // Quiet,
	"FORM_ID"->None,
	"FORMcores"-> 1, 
	"OptimizationLVL"->0,
	"keep_FORM_script"->False,
	"EvalAll"->False,
	"NoNumerator"->False,
	"FORMsubs"->False,
	"stdLTD"->False
};
cLTD[expression_,OptionsPattern[]]:=Module[{
expr=If[Head[expression]===Plus,List@@expression, {expression}]/.Global`SP4[p1_,p2_]:> SP4[p1,p2][OptionValue["loopmom"]]/.toLTDprop[OptionValue["loopmom"]],
energies,i=0,loop0subs, p0subs,spsubs,funsubs,FORMinput, 
filenameID,
runfilename,cLTDfilename,
FORMpath,
exec,result, return,cleanKs,cleanProps,FormEvalStep,
FORMvars={}},
(*Mathematica has a problem with ~*)
SetOptions[cLTD,"WorkingDirectory"->StringReplace["~"->$HomeDirectory][OptionValue["WorkingDirectory"]]];

(*Check that the reserved variables are not being used in the input*)
If[!FreeQ[expression,Global`x|Global`xbar],Print["x and xbar are used internally by FORM. Pease use different variables in your expression"];Abort[]];

PrintTemporary["Mapping variables for FORM expression"];
(*Sanitisation*)
cleanKs=Table[k->ToExpression["nonLoopk"<>ToString[i++]],{k, Complement[Table[ToExpression["k"<>ToString[n]],{n,0,9}],OptionValue["loopmom"]]}];i=0;
cleanProps=Table[noLoopProp->ToExpression["noLoopProp"<>ToString[i++]],{noLoopProp,Union[Cases[expr,Global`prop[__],Infinity]]}];i=0;
expr = expr/.cleanKs/.cleanProps;

(*map FORM symbols to input values*)
energies = Table[(While[!FreeQ[expr,ToExpression["E"<>ToString[i++]]],Null];ToExpression["E"<>ToString[i-1]])->e,{e,Last/@Union[Cases[expr,LTDprop[__],Infinity]]}];i=0;
loop0subs = Table[l->ToExpression["LTDk"<>ToString[i++]],{l,Select[OptionValue["loopmom"],!FreeQ[expr,#]&]}];i=0;
If[Length[loop0subs] == 0, Print["No loop momenta"]; Abort[]];
p0subs = Table[p0->ToExpression["p0Component"<>ToString[i++]],{p0, Union[Cases[expr,p_?(FreeQ[#,Alternatives@@OptionValue["loopmom"]]&)[0],Infinity]]}];i=0;
(*Print[cleanProps];*)
expr = expr/.ReplaceAll[energies,Rule[a_,b_]:>Rule[b,a]]/.p0subs/.Table[k[0]->k,{k,OptionValue["loopmom"]}];

(*Check for scalar products in the numerator*)
spsubs = Table[x->(While[!FreeQ[expr,ToExpression["sp"<>ToString[i++]]],Null];ToExpression["sp"<>ToString[i-1]]),{x,Union[Cases[expr,_Dot,Infinity]]}];i=0;
expr = expr/.spsubs/.loop0subs;
(*Remove Mathematica Functions*)
funsubs=Table[fun->ToExpression["MMfunction"<>ToString[i++]],{fun,Select[Variables[expr],MatchQ[#,_?(#=!=cLTDPrivate`LTDprop&)[___]]&]}];i=0;
expr = expr/.funsubs;
funsubs=funsubs/.Reverse/@loop0subs;

(*Extract FORM variable*)
FORMvars = Join[FORMvars, Union[Flatten[Variables[#[[1]]]&/@Union[Cases[expr,_LTDprop,Infinity]]]]/.{
	Alternatives@@Last/@loop0subs:>Nothing,
	Alternatives@@Last/@cleanKs:>Nothing
	}];
FORMvars = Join[FORMvars, Variables[expr]/.{
	_LTDprop:>Nothing,
	Alternatives@@Last/@spsubs:>Nothing
	}];
FORMvars = Union[FORMvars];

If[OptionValue["FORMsubs"],Return[{Reverse/@Join[funsubs,spsubs,p0subs,cleanProps,cleanKs], energies/.Reverse/@cleanKs}]];

(*Create strings to send to FORM*)
PrintTemporary["Create FORM expression"];
FORMinput=StringReplace[{"cLTDPrivate`"->"","LTD"->"","[":>"(","]"->")","I"->"iComplex"}][ToString[Plus@@expr,InputForm]];
(*Print[FORMinput];*)

PrintTemporary["Execute FORM"];
(*FORM Executable*)
If[Head[OptionValue["FORMcores"]]=!=Integer, Print["FORMcores must be an integer"]; Abort[]];
If[OptionValue["FORMcores"]<1, Print["FORMcores must be a positive integer"]; Abort[]];
If[Head[OptionValue["OptimizationLVL"]]=!=Integer, Print["OptimizationLVL must be an integer"]; Abort[]];
If[!(0<=OptionValue["OptimizationLVL"]<=4), Print["OptimizationLVL must be an integer between 0 and 4"]; Abort[]];
FORMpath=If[OptionValue["FORMcores"]>1,
 If[OptionValue["tFORMpath"]=="tform","tform",
    FileInformation[OptionValue["tFORMpath"],"AbsoluteFileName"]],
 If[OptionValue["FORMpath"]=="form","form",
    FileInformation[OptionValue["FORMpath"],"AbsoluteFileName"]]
];
If[FORMpath=!="form"&& FORMpath=!="tform"&&
FileInformation[OptionValue["FORMpath"],"AbsoluteFileName"]==={},Print["Cannot find ",OptionValue["FORMpath"]];Abort[]];

(*Call FORM*)
filenameID=If[OptionValue["FORM_ID"]===None,
	StringJoin@@RandomSample[Join[Alphabet[], ToUpperCase /@ Alphabet[]], 10],
	OptionValue["FORM_ID"]
];
runfilename=OptionValue["WorkingDirectory"]<>"/cLTD_"<>ToString[filenameID]<>".frm";
cLTDfilename=OptionValue["WorkingDirectory"]<>"/cLTD_out_"<>ToString[filenameID]<>".out";
(*Print[{FORMpath,runfilename}];*)
(*Print[FORMfile["cLTD_out_"<>ToString[filenameID],FORMinput,Length[loop0subs],alPATH]];*)
Export[runfilename,FORMfile["cLTD_out_"<>ToString[filenameID],FORMinput,Length[loop0subs],FORMvars,OptionValue["OptimizationLVL"],OptionValue["stdLTD"]],"Text"];

exec=If[OptionValue["NoNumerator"],
	{FORMpath,If[OptionValue["FORMcores"]>1,"-w"<>ToString[OptionValue["FORMcores"]],Nothing],"-D", "NoNumerator",runfilename},
	{FORMpath,If[OptionValue["FORMcores"]>1,"-w"<>ToString[OptionValue["FORMcores"]],Nothing],runfilename}];
return = RunProcess[exec,ProcessDirectory->OptionValue["WorkingDirectory"]];
If[return["ExitCode"] != 0, Print[return]];

PrintTemporary["Retrieve result"];
result=ToExpression[StringReplace[{" "|"\\"|"\n"->"","FormEvalStep("~~x:DigitCharacter..~~")":>ToString[FormEvalStep]~~"["~~x~~"]"}][Import[cLTDfilename,"Text"]]];
If[!OptionValue["keep_FORM_script"],DeleteFile[{runfilename,cLTDfilename}]];
PrintTemporary["Map back to mathematica variables"];
result = result/.Reverse/@funsubs\
				/.Reverse/@spsubs\
				/.Reverse/@p0subs\
				/.Reverse/@cleanProps\
				/.Reverse/@cleanKs\
				/.Global`norm->cLTD`cLTDnorm\
				/.cLTD`cLTDnorm[a_]:>cLTD`cLTDnorm[a *Power[-I,Length[loop0subs]]]\
				/.Global`iComplex->I;

{If[OptionValue["OptimizationLVL"]==0,Collect[result,_cLTDnorm],result],
 energies/.Reverse/@cleanKs}//.If[OptionValue["EvalAll"],{Global`den[a_]:>1/a,cLTD`cLTDnorm[a_]:>1/a}, {}]
]
SetAttributes[cLTD,Protected];


End[];
EndPackage[];


(* ::Subsection::Closed:: *)
(*Cross - free family LTD expression generator*)


(* ::Subsubsection::Closed:: *)
(*Graph builders*)


(* ::Input::Initialization:: *)
cFFAssignSignaturesStartingFromEdge[g_Graph,startingVertex_,InputSignatures_]:=Module[
{allIncomingEdges, allOutgoingEdges,otherVertexAndDirection,edgesWithoutSignature,signatures},
signatures=Association[Table[
k->InputSignatures[k],{k,Keys[InputSignatures]}]];
allIncomingEdges=EdgeList[g,DirectedEdge[_,startingVertex,_]];
allOutgoingEdges=EdgeList[g,DirectedEdge[startingVertex,_,_]];
edgesWithoutSignature=Select[Join[allIncomingEdges,allOutgoingEdges],Not[MemberQ[Keys[signatures],#]]&];

If[Length[edgesWithoutSignature]==1,
Do[
otherVertexAndDirection=e/.{DirectedEdge[startingVertex,a_,_]:>{a,1},DirectedEdge[a_,startingVertex,_]:>{a,-1}};
AssociateTo[signatures,e->(otherVertexAndDirection[[2]]*(
Total[Table[signatures[eIN],{eIN,Select[allIncomingEdges,MemberQ[Keys[signatures],#]&]}]]
-Total[Table[signatures[eOUT],{eOUT,Select[allOutgoingEdges,MemberQ[Keys[signatures],#]&]}]]
))];
signatures=cFFAssignSignaturesStartingFromEdge[g,otherVertexAndDirection[[1]],signatures];
,{e,edgesWithoutSignature}];
];
signatures
]


(* ::Input::Initialization:: *)
cFFGetEdgeSignatures[g_Graph,OptionsPattern[{lmb->Automatic}]]:=Module[
{spanningTree,selectedLMB, leaves,startEdge, signatures,externalVertices,externalEdges, tmpA,leafEdge},
externalVertices=Sort[Select[VertexList[g],VertexDegree[g,#]===1&]];
If[
OptionValue[lmb]===Automatic,
spanningTree=FindSpanningTree[g];
selectedLMB=SortBy[Complement[
EdgeList[g],
EdgeList[spanningTree]
],(#/.{DirectedEdge[a_,b_,name_]:>name})&];
,
spanningTree=DirectedGraph[Select[EdgeList[g],Not[MemberQ[OptionValue[lmb],(#/.{DirectedEdge[a_,b_,id_]:>id})]]&]];
selectedLMB=SortBy[
Select[EdgeList[g],MemberQ[OptionValue[lmb],(#/.{DirectedEdge[a_,b_,id_]:>id})]&],
Position[OptionValue[lmb],(#/.{DirectedEdge[a_,b_,id_]:>id})][[1]][[1]]&
];
];

(* Assign signatures to externals now *)
signatures=Association[Table[
tmpA=EdgeList[g,DirectedEdge[externalVertices[[iev]],_,_]];
(* It is much simpler if we force all externals incoming so check for this here *)
If[ Length[tmpA]!=1,Print["ERROR: please define all externals as incoming in your topology."];Return[Null];];
tmpA = If[Length[tmpA]==1,tmpA[[1]],EdgeList[g,DirectedEdge[_,externalVertices[[iev]],_]][[1]]];
tmpA->{Table[0,{i,Length[selectedLMB]}],Table[If[iev===j,1,0],{j,Length[externalVertices]}]}
,{iev,Length[externalVertices]}
]];

(* And LMB edges now *)
signatures=AssociateTo[signatures,
Table[
selectedLMB[[ie]]->{Table[If[ie===j,1,0],{j,Length[selectedLMB]}],Table[0,{i,Length[externalVertices]}]},{ie,Length[selectedLMB]}
]
];

(* And all others now *)
leaves=Select[VertexList[spanningTree],VertexDegree[spanningTree,#]===1&];
leaves=Table[
leafEdge=EdgeList[spanningTree,DirectedEdge[leaf,_,_]];
If[Length[leafEdge]==1,
leafEdge=leafEdge[[1]];
If[MemberQ[Keys[signatures],leafEdge],
leafEdge/.{DirectedEdge[_,x_,_]:>x},leaf],
leafEdge=EdgeList[spanningTree,DirectedEdge[_,leaf,_]][[1]];
If[MemberQ[Keys[signatures],leafEdge],
leafEdge/.{DirectedEdge[x_,_,_]:>x},leaf]
]
,{leaf,leaves}];
Do[signatures=cFFAssignSignaturesStartingFromEdge[g,leaf,signatures];,{leaf,leaves}];

signatures

]


(* ::Input::Initialization:: *)
AssignSignatures[INg_Graph,OptionsPattern[{masses->Null,lmb->Automatic}]]:=Module[{g,signatures},
g=DirectedGraph[Table[EdgeList[INg][[ie]]/.{DirectedEdge[a_,b_]:>DirectedEdge[a,b,ie]},{ie,Length[EdgeList[INg]]}]];
signatures= cFFGetEdgeSignatures[g,lmb->OptionValue[lmb]];
DirectedGraph[Table[e/.{DirectedEdge[a_,b_,c_]:>DirectedEdge[a,b,
<|"id"->c,"sig"->signatures[e],"mass"->If[OptionValue[masses]===Null,0,If[MemberQ[Keys[OptionValue[masses]],c],OptionValue[masses][c],0]]|>
]},{e,EdgeList[g]}],VertexLabels->Automatic,EdgeLabels->None]
]


(* ::Input::Initialization:: *)
cFFGenerateMomentaLabels[g_Graph]:=Module[
{oneSig},
oneSig=(EdgeList[g][[1]]/.{DirectedEdge[_,_,props_]:>props})["sig"];
{
Table[ToExpression["k"<>ToString[i]],{i,Length[oneSig[[1]]]}],
Table[ToExpression["p"<>ToString[i]],{i,Length[oneSig[[2]]]}]
}
]


(* ::Input::Initialization:: *)
cFFGetLMBEdges[g_Graph]:=Module[
{lmbEdges},
lmbEdges=Select[EdgeList[g],(Function[s,And[AllTrue[s[[2]],Function[si,si===0]],Count[s[[1]],0]===Length[s[[1]]]-1,Count[s[[1]],1]===1]]@((#/.{DirectedEdge[_,_,props_]:>props})["sig"]))&];
lmbEdges=SortBy[lmbEdges,Position[(#/.{DirectedEdge[_,_,props_]:>props})["sig"][[1]],1][[1]][[1]]&];
DeleteDuplicatesBy[lmbEdges,(#/.{DirectedEdge[_,_,props_]:>props})["sig"]&]
]


(* ::Input::Initialization:: *)
cFFGetExternalEdges[g_Graph]:=Module[
{externalEdges},
externalEdges=Select[EdgeList[g],(Function[s,And[AllTrue[s[[1]],Function[si,si===0]],Count[s[[2]],0]===Length[s[[2]]]-1,Count[s[[2]],1]===1]]@((#/.{DirectedEdge[_,_,props_]:>props})["sig"]))&];
externalEdges=SortBy[externalEdges,Position[(#/.{DirectedEdge[_,_,props_]:>props})["sig"][[2]],1][[1]][[1]]&]
]


(* ::Subsubsection::Closed:: *)
(*cFF Generator*)


(* ::Input::Initialization:: *)
cFFBugFixedEdgeContract[g_Graph,edgeToContract_]:=Module[{nodeToMap,newEdges},
nodeToMap=(edgeToContract/.{DirectedEdge[a_,b_,_]:>{a->b}});
newEdges=Table[e/.{DirectedEdge[a_,b_,props_]:>DirectedEdge[a/.nodeToMap,b/.nodeToMap,props]},{e,DeleteCases[EdgeList[g],edgeToContract]}];
(* Remove self-loops too here *)
newEdges=Select[newEdges,(#/.DirectedEdge[a_,b_,_]:>a!=b)&];
DirectedGraph[newEdges]
]


(* ::Input::Initialization:: *)
cFFGraphContributeToCFFQ[g_Graph]:=AcyclicGraphQ[g]


(* ::Input::Initialization:: *)
cFFEnumerateAcyclicOrientations[g_Graph,OptionsPattern[{Orientations->Automatic}]]:=Module[
{el,allOrientations},
el=EdgeList[g];
allOrientations=Select[
Table[
{o,DirectedGraph[Table[If[o[[iE]],el[[iE]]/.{DirectedEdge[a_,b_,name_]:>DirectedEdge[b,a,name]},el[[iE]]],{iE,Length[el]}],VertexLabels->Automatic,EdgeLabels->Automatic]},{o,Tuples[{True,False},Length[el]]}]
,cFFGraphContributeToCFFQ[#[[2]]]&
];
If[OptionValue[Orientations]===Automatic,allOrientations,
Select[allOrientations,
MemberQ[OptionValue[Orientations],#[[1]]]&
]
]
]


(* ::Input::Initialization:: *)
cFFCrossFreeFamilyReduce[g_Graph,OptionsPattern[{ConvertToNormalisedFormat->False,DEBUG->False}]]:=Module[
{EsurfaceSign,sources, sinks,nodeSelected,allEdgesToDelete,allTerms,Esurf,allContractedGraphs},

If[OptionValue[DEBUG],
Print[StringJoin[Table["-",{i,Range[60]}]]];
Print["Now considering the following graph:"];
Print[DirectedGraph[g,VertexLabels->Automatic,
EdgeLabels->Table[e->(e/.{DirectedEdge[_,_,props_]:>props})["id"],{e,EdgeList[g]}]
]];
];

If[Length[VertexList[g]]==2,
Esurf=Table[e/.{DirectedEdge[a_,b_,name_]:>name["id"]},{e,EdgeList[g]}];
If[OptionValue[DEBUG],Print["Reached termination condition, returning Esurf:",Esurf];];
EsurfaceSign=1;
If[OptionValue[ConvertToNormalisedFormat],
Return[{{{EsurfaceSign,Esurf}}}];,Return[EsurfaceSign(1/Total[Table[OSE[name],{name,Esurf}]])];
];
];

sources=GraphComputation`SourceVertexList[g];
sinks=GraphComputation`SinkVertexList[g];
If[OptionValue[DEBUG],Print["Found the following sources and sinks(first eligible one will be considered):\nsource ",sources,"\nsinks ",sinks];];
nodeSelected=None;
Do[
If[
Length[WeaklyConnectedComponents[Graph[VertexList[g],Select[EdgeList[g],Not[MatchQ[#,DirectedEdge[s,_,_]]]&]]]]<=2,
nodeSelected=s;
allEdgesToDelete=EdgeList[g,DirectedEdge[nodeSelected,_,_]];
EsurfaceSign=1;
Break[];
];
,{s,sources}
];
Do[
If[
Length[WeaklyConnectedComponents[Graph[VertexList[g],Select[EdgeList[g],Not[MatchQ[#,DirectedEdge[_,s,_]]]&]]]]<=2,
nodeSelected=s;
allEdgesToDelete=EdgeList[g,DirectedEdge[_,nodeSelected,_]];
EsurfaceSign=1;
Break[];
];
,{s,sinks}
];

If[nodeSelected===None,Print["Could not find an eligible source or sink for the following graph:"];Print[g];];

allContractedGraphs=Select[
Table[{edgeToDelete,cFFBugFixedEdgeContract[g,edgeToDelete]},{edgeToDelete,allEdgesToDelete}]
,cFFGraphContributeToCFFQ[#[[2]]]&];
allContractedGraphs=DeleteDuplicatesBy[allContractedGraphs,#[[2]]&];

If[OptionValue[DEBUG],
Print["Edges to be contracted (possibly won't contribute all): ",Table[(e/.{(DirectedEdge[a_,b_,props_]:>props)})["id"],{e,allEdgesToDelete}]];Print["Contributing edge contractions: ",
Table[eInfo[[1]],{eInfo,allContractedGraphs}
]];
];

Esurf={EsurfaceSign,Table[e/.{DirectedEdge[a_,b_,name_]:>name["id"]},{e,allEdgesToDelete}]};
If[OptionValue[DEBUG],Print["Esurf for this term: ",Esurf];];
If[Not[OptionValue[ConvertToNormalisedFormat]],
Esurf=Esurf[[1]]*Total[Table[OSE[name],{name,Esurf[[2]]}]];
];
allTerms=((cFFCrossFreeFamilyReduce[#[[2]],ConvertToNormalisedFormat->OptionValue[ConvertToNormalisedFormat],DEBUG->OptionValue[DEBUG]])&)/@allContractedGraphs;

If[Not[OptionValue[ConvertToNormalisedFormat]],
If[OptionValue[DEBUG],
Print["Complete term returned at this level: ",1/Esurf Total[allTerms]];
Print[StringJoin[Table["-",{i,Range[60]}]]];
];
1/Esurf Total[allTerms],
Flatten[Table[Table[Prepend[ti,Esurf],{ti,t}],{t,allTerms}],1]
]

]


(* ::Input::Initialization:: *)
cFFGenerateEnergyShiftSignatureForEsurf[g_Graph,edgeCuts_,flipConfiguration_]:=Module[
{edgesCut},
edgesCut=Select[EdgeList[g],MemberQ[edgeCuts,(#/.DirectedEdge[_,_,props_]:>props)["id"]]&];
-Total[Table[flipConfiguration[e]*((e/.DirectedEdge[_,_,props_]:>props)["sig"][[2]]),{e,edgesCut}]]
]


(* ::Input::Initialization:: *)
CrossFreeFamilyLTD[g_Graph,OptionsPattern[{ConvertToNormalisedFormat->False,DEBUG->False,Orientations->Automatic}]]:=Module[
{allTerms,reducedGraph, flipConfiguration,allOrientations},

(* act on the reduced graph stripped off external for now *)
reducedGraph=GraphDifference[g,DirectedGraph[cFFGetExternalEdges[g]]];
allOrientations=cFFEnumerateAcyclicOrientations[reducedGraph,Orientations->OptionValue[Orientations]];
Monitor[
allTerms=Table[{allOrientations[[iacg]][[1]],
If[OptionValue[DEBUG],Print["Now considering orientation: "<>ToString[Table[If[o,-1,1],{o,allOrientations[[iacg]][[1]]}]]]];
cFFCrossFreeFamilyReduce[allOrientations[[iacg]][[2]],
ConvertToNormalisedFormat->OptionValue[ConvertToNormalisedFormat],DEBUG->OptionValue[DEBUG]
]},{iacg,Length[allOrientations]}];,
ToString["Evaluating orientation #"<>ToString[iacg]<>"/"<>ToString[Length[allOrientations]]<>" ("<>ToString[Round[100.(iacg/Length[allOrientations]),0.01]]<>"%)"]
];
If[Not[OptionValue[ConvertToNormalisedFormat]],
allTerms,
Table[
flipConfiguration=Association[(Table[#[[iE]]->If[o[[1]][[iE]],-1,1],{iE,Length[#]}]&)@EdgeList[reducedGraph]];
<|
"Orientation"->Association[Table[(EdgeList[reducedGraph][[iOrientation]]/.{DirectedEdge[_,_,props_]:>props})["id"]->If[o[[1]][[iOrientation]],-1,1],{iOrientation,Length[o[[1]]]}]],
"Terms"->Table[<|
"Num"->(1),
"Orderings"->Null,
"Esurfs"->Table[
<|
"overall_sign"->ti[[1]],
"OSE"->ti[[2]],
"shift_signature"->cFFGenerateEnergyShiftSignatureForEsurf[g,ti[[2]],flipConfiguration]
|>
,{ti,t}]
|>,{t,o[[2]]}]
|>
,{o,allTerms}
]
]
]


(* ::Subsubsection::Closed:: *)
(*Evaluators*)


(* ::Input::Initialization:: *)
cFFCountOPs[expr_]:=Module[{leaves},
leaves=Level[expr,{-1},Heads->True];
Select[Tally[leaves],MemberQ[{Times,Power,Plus},#[[1]]]&]
]


(* ::Input::Initialization:: *)
cFFCompareOPCount[cLTDexpr_,cFFFactorisedExpr_]:=Module[
{cLTDOPcount,cFFOPcount,opTypes},
cLTDOPcount=Association[Table[count[[1]]->count[[2]],{count,cFFCountOPs[cLTDexpr[[1]]]}]];
cFFOPcount=Association[Table[count[[1]]->count[[2]],{count,cFFCountOPs[Total[Table[t[[2]],{t,cFFFactorisedExpr}]]]}]];
Dataset[
Table[
<|"op type"->opType,
"cLTD"->cLTDOPcount[opType],"cFF"->cFFOPcount[opType],
"cFF-cLTD"->cFFOPcount[opType]-cLTDOPcount[opType],
"\[CapitalDelta](cFF-cLTD) %"->Round[200. (cFFOPcount[opType]-cLTDOPcount[opType])/(cFFOPcount[opType]+cLTDOPcount[opType]),0.1]
|>
,{opType,Keys[cLTDOPcount]}]
]
]


(* ::Input::Initialization:: *)
GenerateRandomSample[g_Graph,OptionsPattern[{Seed->Automatic}]]:=Module[{momLabels,RNSeed,RN,NumericalSample},
momLabels=cFFGenerateMomentaLabels[g];
RNSeed=If[OptionValue[Seed]===Automatic,RandomInteger[{1,100}],OptionValue[Seed]];
RN = Table[Prime[i]/Prime[i+1],{i,RNSeed,(3*(Length[momLabels[[1]]]+Max[Length[momLabels[[2]]]-1,0])+Max[Length[momLabels[[2]]]-1,0])+RNSeed-1}];
NumericalSample=Join[
Table[momLabels[[1]][[ki+1]]->Table[RN[[ki*3+j]],{j,3}],{ki,0,Length[momLabels[[1]]]-1}],
If[Length[momLabels[[2]]]>1,
Join[
Table[momLabels[[2]][[pi+1]]->Table[RN[[(Length[momLabels[[1]]]+pi)*3+j]],{j,3}],{pi,0,Length[momLabels[[2]]]-2}],
Table[ToExpression[ToString[momLabels[[2]][[pi]]]<>"E"]->RN[[(Length[momLabels[[1]]]+Length[momLabels[[2]]]-1)*3+pi]],{pi,1,Length[momLabels[[2]]]-1}]
],
{}
]
];
If[Length[momLabels[[2]]]>1,
NumericalSample=Join[NumericalSample,
{
momLabels[[2]][[-1]]->Simplify[-Total[momLabels[[2]][[;;-2]]]/.NumericalSample],
ToExpression[ToString[momLabels[[2]][[-1]]]<>"E"]->Simplify[-Total[Table[ToExpression[ToString[pi]<>"E"],{pi,momLabels[[2]][[;;-2]]}]]/.NumericalSample]
}
];,
If[Length[momLabels[[2]]]==1,
NumericalSample=Join[NumericalSample,
{
momLabels[[2]][[-1]]->{0,0,0},
ToExpression[ToString[momLabels[[2]][[-1]]]<>"E"]->0
}
];
];
];
NumericalSample
]


Options[GeneratecLTDExpression]={
Num->None,
	"FORMpath"->"form",
	"tFORMpath"->"tform","WorkingDirectory"->Check[NotebookDirectory[],Print["Cannot find Notebook Directory!\nMust define the WorkingDirectory option before running."]; None] // Quiet,
	"FORM_ID"->None,
	"FORMcores"-> 1, 
	"OptimizationLVL"->0,
	"keep_FORM_script"->False,
	"stdLTD"->False
};


(* ::Input::Initialization:: *)
GeneratecLTDExpression[g_Graph,OptionsPattern[]]:=Module[
{props,momLabels,momEnergyLabels,num,cLTDExpr,reducedGraph,signatures},
momLabels=cFFGenerateMomentaLabels[g];
momEnergyLabels={Table[ml[0],{ml,momLabels[[1]]}],Table[ml[0],{ml,momLabels[[2]]}]};
num=If[OptionValue[Num]===None,1,
signatures=Association[Table[((#["id"]->#["sig"])&)[e/.{DirectedEdge[_,_,props_]:>props}],{e,EdgeList[g]}]];
Expand[(OptionValue[Num]
/.{q[i_]:>(signatures[i][[1]] . momLabels[[1]]+signatures[i][[2]] . momLabels[[2]])}
/.{qE[i_]:>(signatures[i][[1]] . momEnergyLabels[[1]]+signatures[i][[2]] . momEnergyLabels[[2]])}
)]
];
reducedGraph=GraphDifference[g,DirectedGraph[cFFGetExternalEdges[g]]];
props=Table[(
 prop[#["sig"][[1]] . momLabels[[1]]+#["sig"][[2]] . momLabels[[2]],#["mass"]]&)@(e/.{DirectedEdge[_,_,props_]:>props}),{e,EdgeList[reducedGraph]}
];

cLTDExpr = cLTD[num( Times@@ props), NoNumerator -> OptionValue[Num]===None, loopmom -> momLabels[[1]],EvalAll -> True,
"FORMpath"->OptionValue["FORMpath"],
"tFORMpath"->OptionValue["tFORMpath"],
"WorkingDirectory"->OptionValue["WorkingDirectory"],
"FORM_ID"->OptionValue["FORM_ID"],
"FORMcores"->OptionValue["FORMcores"],
"OptimizationLVL"->OptionValue["OptimizationLVL"],
"keep_FORM_script"->OptionValue["keep_FORM_script"],
"stdLTD"->OptionValue["stdLTD"]
];
cLTDExpr=cLTDExpr/.Table[pi[0]->ToExpression[ToString[pi]<>"E"],{pi,momLabels[[2]]}];

(* Correct for an overall (-1)^Nloops missing when NoNumerator is set to False *)
cLTDExpr={If[Not[OptionValue[Num]===None],(-1)^Length[momLabels[[1]]],1]cLTDExpr[[1]],cLTDExpr[[2]]};

cLTDExpr
]


(* ::Input::Initialization:: *)
EvalcLTD[cLTDexpr_,numerics_]:=cLTDexpr[[1]]/.(cLTDexpr[[2]]/.numerics)/.numerics


(* ::Input::Initialization:: *)
ComputeOSEReplacements[g_Graph,numerics_]:=Module[{momLabels,ks,ps,edgeProperties,OSRepl},
momLabels=cFFGenerateMomentaLabels[g];
ks=momLabels[[1]];
ps=momLabels[[2]];
OSRepl=Table[
edgeProperties=(e/.DirectedEdge[_,_,props_]:>props);
OSE[edgeProperties["id"]]:>Evaluate[Sqrt[
(edgeProperties["sig"][[1]] . ks+edgeProperties["sig"][[2]] . ps) .
(edgeProperties["sig"][[1]] . ks+edgeProperties["sig"][[2]] . ps)+edgeProperties["mass"]^2
]/.numerics],{e,EdgeList[g]}
];
OSRepl
]


(* ::Input::Initialization:: *)
EvalcFF[g_Graph,cFFexpr_,numerics_,OptionsPattern[{DEBUG->False, Num->None}]]:=Module[
{ResPerOrientations,OSEReplacement,ShiftsReplacement,lmbEdges,externalEdges,
signatures,momLabels,reducedGraph,processedNum, numForThisOrientation,internalOSEs},
lmbEdges=cFFGetLMBEdges[g];
externalEdges=cFFGetExternalEdges[g];
reducedGraph=GraphDifference[g,DirectedGraph[externalEdges]];
OSEReplacement = ComputeOSEReplacements[g,numerics];
If[OptionValue[DEBUG],Print["Onshell energy values: ",OSEReplacement];];
ShiftsReplacement = Table[pE[iExt]:>Evaluate[ToExpression["p"<>ToString[iExt]<>"E"]],{iExt,Length[externalEdges]}]/.numerics;
If[Not[OptionValue[Num]===None],
momLabels=cFFGenerateMomentaLabels[g];
signatures=Association[Table[((#["id"]->#["sig"])&)[e/.{DirectedEdge[_,_,props_]:>props}],{e,EdgeList[g]}]];
processedNum=TensorExpand[OptionValue[Num]/.{SP4:>Cross}]/.{Cross:>SP4};
processedNum=((((processedNum
/.{SP4[a_,b_]:>(a[0]*b[0]-a . b)})/.{q[i_][0]:>qE[i]})
/.{q[i_]:>(signatures[i][[1]] . momLabels[[1]]+signatures[i][[2]] . momLabels[[2]])})
/.Table[qE[(externalEdges[[iExt]]/.{DirectedEdge[_,_,props_]:>props})["id"]]:>Evaluate[pE[iExt]],{iExt,Length[externalEdges]}]);
If[OptionValue[DEBUG],Print["Numerator before numerics substition: ",processedNum];];
processedNum=processedNum/.ShiftsReplacement/.numerics;
If[OptionValue[DEBUG],Print["Numerator after numerics substition: ",processedNum];];
];
ResPerOrientations=Table[
If[OptionValue[Num]===None,
numForThisOrientation=1;
,
numForThisOrientation=processedNum/.Table[qE[edgeID]->(t["Orientation"][edgeID]*OSE[edgeID]),{edgeID,Keys[t["Orientation"]]}];
numForThisOrientation=numForThisOrientation/.OSEReplacement;
];
If[OptionValue[DEBUG],Print["Numerator with signed on-shell energies substituted for orienations "<>ToString[t["Orientation"]]<>" : ",numForThisOrientation];];
Total[Table[
numForThisOrientation ((st["Num"]/.OSEReplacement)/.numerics)/
(Times@@Table[
(eta["overall_sign"]*((Total[Table[OSE[e],{e,eta["OSE"]}]]/.OSEReplacement)+
((Table[pE[s],{s,Length[externalEdges]}] . eta["shift_signature"])/.ShiftsReplacement)))
,{eta, st["Esurfs"]}])
,{st,t["Terms"]}
]]
,{t,cFFexpr}
];
( (I)^Length[lmbEdges] / ((Times@@Table[-2*OSE[(e/.{DirectedEdge[_,_,props_]:>props})["id"]],{e,EdgeList[reducedGraph]}])/.OSEReplacement) )*Total[ResPerOrientations]
];
