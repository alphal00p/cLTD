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

"<>Style["Options",Bold]~ToString~StandardForm<>": {
          \"loopmom\"->{loopmomsymbol1,loopmomsymbol2,...} (default: {k1,k2,k3,k4})
          \"FORMpath\"\[Rule]\"form\",
          \"WorkingDirectory\"\[Rule]Directory[],
          \"FORM_ID\"\[Rule]None,
          \"keep_FORM_script\"\[Rule]False,
          \"EvalAll\"\[Rule]False
}
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
  The expression contains the den[] and cLTDnorm[] function which are to be evaluated as 
		den[a_] :> 1/a, cLTDnorm[a_] :> 1/a
";

Print[" :::::::::::::::::::::::: cLTD ::::::::::::::::::::::::"]
Print["Authors: Z. Capatti, V. Hirschi, D. Kermanschah, A. Pelloni, B. Ruijl\n"];
Print[cLTD::usage];


(* ::Section:: *)
(*Private*)


Begin["cLTDPrivate`"];


(* ::Subsection::Closed:: *)
(*FORM file*)


FORMhead = "Format 100;
Auto S n,y,E,p,k;
CF ncmd, num, Eres, Echain;
CF den, prop, norm, error;
CF x, xbar;

set energies:E0,...,E1000;
set shifts:p0,...,p1000;
";


FORMpf = "multiply norm(1)*num();
repeat id norm(E?)*prop(y0?,y1?) = norm(2*E*y1)*(den(y0-y1)-den(y0+y1)); 
.sort

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
    id num(?y) = num(1,?y,ncmd());

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


FORMfile[filename_,expr_,nLoops_,FORMvars_]:=Module[{
file="#--\n"},
file = file <> FORMhead;
file = file <> "Auto S sp,"<>StringRiffle[ToString/@FORMvars,","]<>";\n";
file = file <> "#define LOOPS \""<>ToString[nLoops]<>"\"\n";
file = file <> "#define topoID \""<>filename<>"\"\n";
file = file <> "L F = "<>expr<>";\n";
file = file <> FORMpf <> "\n.sort\n";
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
file = file <> "#endif\n.sort\nFormat mathematica;\n#write<"<>filename<>".out> \"%E\" F;\n.end";
file
]


(* ::Subsection:: *)
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
	"WorkingDirectory"->NotebookDirectory[],
	"FORM_ID"->None, 
	"keep_FORM_script"->False,
	"EvalAll"->False};
cLTD[expression_,OptionsPattern[]]:=Module[{
expr=If[Head[expression]===Plus,List@@expression, {expression}]/.Global`SP4[p1_,p2_]:> SP4[p1,p2][OptionValue["loopmom"]]/.toLTDprop[OptionValue["loopmom"]],
energies,i=0,loop0subs, p0subs,spsubs,FORMinput, 
filenameID,
runfilename,cLTDfilename,
FORMpath,
result, return,cleanKs,cleanProps,
FORMvars={}},

(*Check that the reserved variables are not being used in the input*)
If[!FreeQ[expression,Global`x|Global`xbar],Print["x and xbar are used internally by FORM. Pease use different variables in your expression"];Abort[]];

PrintTemporary["Mapping variables for FORM expression"];
(*Sanitisation*)
cleanKs=Table[k->ToExpression["nonLoopk"<>ToString[i++]],{k, Complement[Table[ToExpression["k"<>ToString[n]],{n,0,9}],OptionValue["loopmom"]]}];i=0;
cleanProps=Table[noLoopProp->ToExpression["noLoopProp"<>ToString[i++]],{noLoopProp,Union[Cases[expr,Global`prop[__],Infinity]]}];i=0;

(*map FORM symbols to input values*)
energies = Table[(While[!FreeQ[expr,ToExpression["E"<>ToString[i++]]],Null];ToExpression["E"<>ToString[i-1]])->e,{e,Last/@Union[Cases[expr,LTDprop[__],Infinity]]}];i=0;
loop0subs = Table[l->ToExpression["LTDk"<>ToString[i++]],{l,Select[OptionValue["loopmom"],!FreeQ[expr,#]&]}];i=0;
If[Length[loop0subs] == 0, Print["No loop momenta"]; Abort[]];
p0subs = Table[p0->ToExpression["p0Component"<>ToString[i++]],{p0, Union[Cases[expr,p_?(FreeQ[#,Alternatives@@OptionValue["loopmom"]]&)[0],Infinity]]}];i=0;
(*Print[cleanProps];*)
expr = expr/.ReplaceAll[energies,Rule[a_,b_]:>Rule[b,a]]/.cleanProps/.cleanKs/.p0subs/.Table[k[0]->k,{k,OptionValue["loopmom"]}];

(*Check for scalar products in the numerator*)
spsubs = Table[x->(While[!FreeQ[expr,ToExpression["sp"<>ToString[i++]]],Null];ToExpression["sp"<>ToString[i-1]]),{x,Union[Cases[expr,_Dot,Infinity]]}];i=0;
expr = expr/.spsubs/.loop0subs;

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

(*Create strings to send to FORM*)
PrintTemporary["Create FORM expression"];
FORMinput=StringReplace[{"cLTDPrivate`"->"","LTD"->"","[":>"(","]"->")"}][ToString[Plus@@expr,InputForm]];
(*Print[FORMinput];*)

PrintTemporary["Execute FORM"];
(*FORM Executable*)
FORMpath=If[
OptionValue["FORMpath"]=="form", "form",
FileInformation[OptionValue["FORMpath"],"AbsoluteFileName"]
];

(*Call FORM*)
filenameID=If[OptionValue["FORM_ID"]===None,
	StringJoin@@RandomSample[Join[Alphabet[], ToUpperCase /@ Alphabet[]], 10],
	OptionValue["FORM_ID"]
];
runfilename=OptionValue["WorkingDirectory"]<>"/cLTD_"<>ToString[filenameID]<>".frm";
cLTDfilename=OptionValue["WorkingDirectory"]<>"/cLTD_out_"<>ToString[filenameID]<>".out";
(*Print[{FORMpath,runfilename}];*)
(*Print[FORMfile["cLTD_out_"<>ToString[filenameID],FORMinput,Length[loop0subs],alPATH]];*)
Export[runfilename,FORMfile["cLTD_out_"<>ToString[filenameID],FORMinput,Length[loop0subs],FORMvars],"Text"];
return = RunProcess[{FORMpath,runfilename},ProcessDirectory->OptionValue["WorkingDirectory"]];
If[return["ExitCode"] != 0, Print[return]];

PrintTemporary["Retrieve result"];
result=ToExpression[StringReplace[" "|"\\"|"\n"->""][Import[cLTDfilename,"Text"]]];
If[!OptionValue["keep_FORM_script"],DeleteFile[{runfilename,cLTDfilename}]];
PrintTemporary["Map back to mathematica variables"];
result = result/.Reverse/@spsubs\
				/.Reverse/@p0subs\
				/.Reverse/@cleanKs\
				/.Reverse/@cleanProps\
				/.Global`norm->cLTD`cLTDnorm\
				/.cLTD`cLTDnorm[a_]:>cLTD`cLTDnorm[a *Power[-I,Length[loop0subs]]];
				
{Collect[result,_cLTDnorm], energies}//.If[OptionValue["EvalAll"],{Global`den[a_]:>1/a,cLTD`cLTDnorm[a_]:>1/a}, {}]
]


End[];
EndPackage[];
