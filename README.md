[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7229437.svg)](https://doi.org/10.5281/zenodo.7229437)
#
# cLTD - Causal Loop-Tree Duality
**Authors**: Z. Capatti, V. Hirschi, D. Kermanschah, A. Pelloni, B. Ruijl
![double box diagrams](./img/diags.svg)

# Table of Contents
1. [Introduction](#introduction)
2. [Dependencies](#dependencies)
4. [Usage](#usage)
    1. [Import](#import)
    2. [Options](#options)
    2. [Input format](#input)
    2. [Output format](#output)
4. [Example](#example)
    1. [Scalar integral: Double Box](#exampleDoubleBox)
    2. [Raised power propagators: Triangle](#exampleTriangle)
    3. [**Numerators**](#numerator)
    4. [Cross-free Family generation and LTD from graphs](#LTDfromGraphs)
5. [Reference](#reference)


## Introduction
The integration of a `n`-loop Feynman diagram requires the integrating over `4*n` dimensions.
By means of the Loop-Tree duality, it's possible to bring down this integration to only `3*n` spatial dimensions. 

This new representation consists of the sum over all possible spanning tree of the original diagram. One of the major problems that one can encounter is that these terms will introduce a number of spurious singularities that could lead to a poor numerical integration convergence.

Here we present a **Mathematica Package** with the implementation of the algorithm presented in the paper [Manifestly Causal Loop-Tree Duality](https://arxiv.org/abs/2009.05509), where it was already implemented into a python script.

The master formula reads:

![master-formula](./img/master.svg)

## Dependencies
* **Mathematica**: tested on Mathematica 11.0 and higher.
* **FORM**: can be downloaded from [here](https://github.com/vermaseren/form/releases).

## Usage

### Import
This package represents a Mathematica frontend and requires a working installation of [FORM](https://github.com/vermaseren/form/releases) in order to run successfully.
If the cLTD folder is added to the libraries in Mathematica, one can import it with:
```Mathematica
<< cLTD`
```
or by defining the location of the package:
```Mathematica
Get["<path-to-package-folder>/cLTD.m"]
```
or alternatively directly from github:
```mathematica
Import["https://github.com/apelloni/cLTD/raw/main/cLTD.m"]
SetOptions[cLTD, WorkingDirectory -> "./"]
```


### Options

There are several options that can be set when calling the `cLTD` function. The default values are:
```Mathematica 
In[1]:= Options[cLTD]
Out[1]= {"loopmom" -> {k0, k1, k2, k3}, 
         "FORMpath" -> "form", 
         "tFORMpath" -> "tform", 
         "FORMcores" -> 1 
         "WorkingDirectory" -> NotebookDirectory[], 
         "FORM_ID" -> None,
         "keep_FORM_script" -> False, 
         "EvalAll" -> False,
         "NoNumerator" -> False}
         "stdLTD" -> False}
```

Their usage is as follows:

* `loopmom`: defines the variables that need to be considered loop momenta in the evaluated expression

* `FORMpath`: By default, it assumes that the path to `form` is included in `$PATH`. If this is not the case, one can set this option to the required path.
* `tFORMpath`: By default, it assumes that the path to `tform` is included in `$PATH`. If this is not the case, one can set this option to the required path.

```Mathematica
SetOptions[cLTD, 
    "FORMpath" -> "/home/dude/dir1/dir2/bin/form",
    "tFORMpath" -> "/home/dude/dir1/dir2/bin/tform"]
```

* `FORMcores`: Is positive integer defining with how many cores to run the FORM script.

* `WorkingDirectory`: FORM will need to create some temporary files in order to process the expression. By default, this happens in the same folder as the current Mathematica notebook.

* `"FORM_ID"` and `"keep_FORM_script"`: the temporary file generate by the program have the name `cLTD_<FORM_ID>.frm` and `cLTD_<FORM_ID>.out`.
By default, the `FORM_ID` is filled with a series of random characters in order to avoid possible conflicts with pre-existing files.
As soon as the `cLTD` call retrieves the output from FORM it deletes these files.
It may happen that you want to keep the script and also give it a more precise name.
One way to do it is:
```Mathematica
SetOptions[cLTD, "FORM_ID" -> 1]
SetOptions[cLTD, "keep_FORM_script" -> True]
```

* `EvalAll`: It substitute the `LTDnorm` and `den` functions.
* `NoNumerator`: When computing the manifestly causal expression with a **constant numerator** one can se this value to `True`.  
* `stdLTD`: Produce the standard LTD representation.  

### Input

The input must contain the **propagator** in the form:
 - `prop[momentum, mass]`

One can multiply several such propagators to build its own topology. One such use can be found in the [double-box example](#exampleDoubleBox), or even sum several topologies together. In the latter case, they all need to contain the same amount of loops.

The **numerator** can be expressed by means of one helper function:
 - `Dot[p1,p2]` or `p1.p2`: Euclidean scalar product of the spatial component of two momenta
 - `SP4[p1,p2]`: representation of the 4-dimensional scalar product between two momenta `p1` and `p2` using the metric `(1,-1,-1,-1)`.
 

In case it's not possible to express the numerator in terms of scalar products, one can always write it as an energy polynomial.

Some examples of valid numerator inputs (`k,l`: loop momenta, `p1,p2`: external momenta)
```mathematica
(*Unity*)
num0 = 1;
(*Using scalar products*)
num1 = SP4[p1,k]+SP4[p2,k]+SP4[k,l];
(*Energy polynomial*)
num2 = k[0]*c1+l[0]*c2+k[0]*l[0]*c3+p1[0]*p2[0];
(*Mixed*)
num3 = k[0]*SP4[p1,l] + p2.l;
num4 = k[0] p[0] - k.p (*Equivalent to SP4[k,p]*)
```

Note that only the energy components of the loop momenta are actively involved during the evaluation of the residue.

### Output
The output contains:
 - `cLTDnorm[expr]`: normalization factor, it corresponds to `1/expr`.
 - `den[expr]`: denominators, can be evaluated as `1/expr`.
 - `p.q`: Euclidean scalar product of the spatial component of two momenta
 - `p[0]`: energy component of momentum `p`

The output consists of two values:
`{<cLTD expression>, <Energy substitutions>}`
```mathematica 
(*Bubble output example*)
<< cLTD`
SetOptions[cLTD, "FORMpath" -> "/home/dude/dir1/dir2/bin/form"];

In[1]:= cLTD[prop[k, 0]*prop[k + p1, 0], loopmom -> {k}, EvalAll -> True]
Out[1]= {
        (1/(E0 + E1 - p1[0]) + 1/(E0 + E1 + p1[0]))/(I 4 E0 E1)), 
        {E0 -> Sqrt[k.k], E1 -> Sqrt[k.k + 2 k.p1 + p1.p1]}
    }
```

## Example

### Scalar integral: Double Box <a name="exampleDoubleBox"></a>
For scalar propagators, one should use the helper function `prop[mom, mass]` in order to define the propagators that are present in the diagram:
```Mathematica
(* Get cLTD and define path to FORM *)
<< cLTD`
SetOptions[cLTD, "FORMpath" -> "/home/dude/dir1/dir2/bin/form"];

(* Double box *)
expr = prop[k, 0] * prop[k-p1, 0] * prop[k+p2, 0] * prop[k-l,0]\
     * prop[l+p2, 0] * prop[l+p2+p3, 0] * prop[l-p1, 0];

(* Get the cLTD expression *)
cDoubleBox = cLTD[expr, loopmom -> {k,l}];      
```
In order to numerically integrate the expression above in the remaining 6 dimensions, one will need to regulate the IR divergences.

It is always possible to consider a massive box:
```Mathematica
(* Massive Double box *)
expr = prop[k, m] * prop[k-p1, m] * prop[k+p2, m] * prop[k-l,0]\
     * prop[l+p2, m] * prop[l+p2+p3, m] * prop[l-p1, m];
```

### Raised power propagators: Triangle <a name="exampleTriangle"></a>
The code also allows for raised powers in the propagators. 
For the case of a triangle with off-shell externals, one will need to regulate the UV divergences in order to obtain a finite integral:

```Mathematica
(* Get cLTD and define path to FORM *)
<< cLTD`
SetOptions[cLTD, "FORMpath" -> "/home/dude/dir1/dir2/bin/form"];

(* TriangleBox - UV CT *)
expr = prop[k, 0] * prop[k-p1, 0] * prop[k+p2, 0] - prop[k, mUV]^3;

(* Get the cLTD expression *)
cTriangleUV = cLTD[expr, loopmom -> {k}];      
```

### Numerators <a name="numerator"></a>
The code also supports numerators. 

One can include **Minkowski scalar products** by using `SP4[mom1,mom2]` in the expression.

For the one-loop photon self energy, the numerator takes the form:

![](./img/bubble.svg)

Which can be expressed as: 
```mathematica
(* One-loop photon self energy *)
numerator = - 4 * SP4[k,p1] + 4 * m^2 \
            + 8 * SP4[pol,k] * SP4[cpol,k]\
            + 4 * SP4[pol,p1] * SP4[cpol,k]\
            + 4 * SP4[pol,k] * SP4[cpol,p1];
props = prop[k,m] * prop[k+p1,m];

(* Get the cLTD expression *)
cPhotonSelf = cLTD[numerator * props, loopmom -> {k}];      

```

It could be that the numerator cannot be written in terms of scalar products; in this case, one should write the numerator as a polynomial in the loop momenta's energy components.

### Cross-free Family generation and LTD from graphs <a name="LTDfromGraphs"></a>

The package also offers to generate the cLTD expression directly from a Mathematica directed graph, which can be generated as follow:

```mathematica
DoubleBox2L = DirectedGraph[{
    DirectedEdge[1, 2, 1], DirectedEdge[2, 3, 2],
    DirectedEdge[4, 5, 3], DirectedEdge[5, 6, 4],
    DirectedEdge[1, 4, 5], DirectedEdge[2, 5, 6], 
    DirectedEdge[3, 6, 7],
    DirectedEdge[101, 1, 8], DirectedEdge[102, 4, 9], 
    DirectedEdge[103, 3, 10], DirectedEdge[104, 6, 11]
    }, VertexLabels -> Automatic, EdgeLabels -> Automatic];
```

After which the kinematics and masses can be automatically generated by dressing the edges with extra metadata, including momenta signatures:

```mathematica
DoubleBox2L = DirectedGraph[AssignSignatures[DoubleBox2L,
    masses -> Association[Table[i -> i, {i, Range[7]}]],
    lmb -> {5, 6}]];
```

which can then conveniently be displayed directly in Mathematica:

```mathematica
DirectedGraph[DoubleBox2L, VertexLabels -> Automatic, 
 ImageSize -> 275, GraphLayout -> "HighDimensionalEmbedding",
 EdgeLabels -> 
  Table[e -> 
    Style[("q" <> 
       ToString[((e /. {DirectedEdge[_, _, props_] :> props})[
          "id"])]), FontSize -> 20], {e, EdgeList[DoubleBox2L]}]]
```
![double box diagram from Mathematica](./img/doublebox.png)

The graph object `DoubleBox2L` can then be used for various operations, for instance for generating a random kinematic sample as follows:

```mathematica
DoubleBox2Lnumerics = GenerateRandomSample[DoubleBox2L, Seed -> 1]
```

and then also generate its (c)LTD representation, and then evaluate it for the kinematics above with:

```mathematica
DoubleBox2LcLTDexpr = GeneratecLTDExpression[DoubleBox2L, Num -> DoubleBox2LNumerator];
EvalcLTD[DoubleBox2LcLTDexpr, DoubleBox2Lnumerics] // N // FullForm
```

but also a different LTD representation whose E-surface denominator all correspond to cross-free families of subsets of the graph vertices (see ref. entitled `Threshold singularity structure of Feynman diagrams from triangulations of convex cones`):

```mathematica
DoubleBox2LcFFexpr = CrossFreeFamilyLTD[DoubleBox2L, ConvertToNormalisedFormat -> True];
EvalcFF[DoubleBox2L, DoubleBox2LcFFexpr, DoubleBox2Lnumerics, 
   Num -> DoubleBox2LNumerator] // N // FullForm
```

And a direct mathematica expression for this representation can be obtained with:

```mathematica
DoubleBox2LcFFexprFactorised = 
  CrossFreeFamilyLTD[DoubleBox2L, ConvertToNormalisedFormat -> False];
```

Note that the option `Orientations` can be specified so as to compute only a subset of the contributing acyclic graphs, e.g.:

 `Orientations -> {{True, True, True, True, True, True, True}}`

 Finally, one can also provide numerators with this representation.
 However, because it is TOPT-based, the numerator must be supplied in the edge-momentum representation, that is using the momenta `q[i]` for each edge `i` of the graph.
 Note however that the representation only supports numerator linear at most in `q[i]`. 
 Here is an example for a one-loop box:

 ```mathematica
 Box1L = DirectedGraph[{
    DirectedEdge[1, 2, 1], DirectedEdge[2, 3, 2], 
    DirectedEdge[3, 4, 3], DirectedEdge[4, 1, 4],
    DirectedEdge[101, 1, 5], DirectedEdge[102, 2, 6], 
    DirectedEdge[103, 3, 7], DirectedEdge[104, 4, 8]
}];
Box1L = AssignSignatures[Box1L, masses -> <|1 -> 1, 2 -> 2, 3 -> 3, 4 -> 4|>, lmb -> {1}];
Box1Lnumerics = GenerateRandomSample[Box1L, Seed -> 1]

Print["Note that q[1] = q[4] + q[5]"];
Print[""];
Box1Lnumerator = SP4[q[1], q[1]] ;
Print["Numerator considered: ", Box1Lnumerator]
Box1LcLTDexpr = GeneratecLTDExpression[Box1L, Num -> Box1Lnumerator];
Print["cLTD result: ", EvalcLTD[Box1LcLTDexpr, Box1Lnumerics] // N // FullForm];
Box1LcFFexpr = CrossFreeFamilyLTD[Box1L, ConvertToNormalisedFormat -> True];
Print["cFF result : ", EvalcFF[Box1L, Box1LcFFexpr, Box1Lnumerics, Num -> Box1Lnumerator] // N // FullForm];
Print[""];
Box1Lnumerator = SP4[q[1], q[4] + q[5]];
(* 
Note that it can be given like this alternatively:
Box1Lnumerator = qE[1]*qE[4] + q[1][0]*q[5][0] - q[1] . q[4] - q[1] . q[5];
*)
Print["Numerator considered: ", Box1Lnumerator]
Box1LcLTDexpr = GeneratecLTDExpression[Box1L, Num -> Box1Lnumerator];
Print["cLTD result: ", EvalcLTD[Box1LcLTDexpr, Box1Lnumerics] // N // FullForm];
Box1LcFFexpr = CrossFreeFamilyLTD[Box1L, ConvertToNormalisedFormat -> True];
Print["cFF result : ", EvalcFF[Box1L, Box1LcFFexpr, Box1Lnumerics, Num -> Box1Lnumerator] // N // FullForm];
Print[""];
```

giving:

```
Note that q[1] = q[4] + q[5]

Numerator considered: SP4[q[1],q[1]]
cLTD result: Complex[0.`,-0.00014091625878931386`]
cFF result : 0.`

Numerator considered: SP4[q[1],q[4]+q[5]]
cLTD result: Complex[0.`,-0.00014091625878931386`]
cFF result : Complex[0.`,-0.00014091625878931378`]
```

## Reference 
If you use this program, please cite the corresponding work on which the procedure is based on:

> Z. Capatti, V. Hirschi, D. Kermanschah, A. Pelloni, B. Ruijl, Manifestly Causal Loop-Tree Duality. 
arXiv: 2009.05509 [hep-ph] (Sept. 2020)

> Z. Capatti, Threshold singularity structure of Feynman diagrams from triangulations of convex cones.
arXiv: TO APPEAR [hep-ph]
