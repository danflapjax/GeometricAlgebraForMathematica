(* ::Package:: *)

BeginPackage["GeometricAlgebra`"]

(*core functions*)
e::usage =
"e[i] represents the basis vector with the identifier i.
e[i,j,...] automatically expands to a basis k-vector."
Metric::usage =
"You may assign a value to Metric[i] to change the value of e[i]^2. E.g. Metric[1]=-1 results in e[1]**e[1] being evaluated as -1."
ExpandNCM::usage =
"ExpandNCM[expr] expands geometric products in expr."

(*decomposition of multivectors*)
KVectorPart::usage =
"KVectorPart[expr,k] gives the k-vector part of expr."
ScalarPart::usage =
"ScalarPart[expr] gives the scalar part of expr."
MaxGrade::usage =
"MaxGrade[expr] returns the largest grade in an expression. For example MaxGrade[e[1]+e[2]**e[3]] gives 2."
KVectorDecomposition::usage =
"KVectorDecomposition[m] decomposes a multivector m into a list of its k-vector parts."

(*products*)
GeometricProduct::usage =
"GeometricProduct[a,b,...] gives the geometric product of an arbitrary number of multivectors.
It is equivalent to NonCommutativeMultiply with an automatic expansion."
OuterProduct::usage =
"OuterProduct[a,b] gives the standard outer product of two multivectors."
CommutatorProduct::usage =
"CommutatorProduct[a,b] gives the commutator product of two multivectors, defined as (a**b-b**a)/2."
LeftContraction::usage =
"LeftContraction[a,b] gives the left contraction of two multivectors."
RightContraction::usage =
"RightContraction[a,b] gives the right contraction of two multivectors."
FatDotProduct::usage =
"FatDotProduct gives the 'fat dot' product of two multivectors."
ScalarProduct::usage =
"ScalarProduct[a,b,...] gives the scalar product of an arbitrary number of multivectors, defined as the scalar part of the geometric product."

(*involutions*)
Reversion::usage =
"Reversion[m] gives the reversion of a multivector m, defined as multiplying each k-vector part by (-1)^(k(k-1)/2)"
Involution::usage =
"Involution[m] gives the involution of a multivector m, defined as multiplying each k-vector part by (-1)^k"
Conjugation::usage =
"Reversion[m] gives the reversion of a multivector m, defined as multiplying each k-vector part by (-1)^(k(k+1)/2)"

(*norms*)
MultivectorNorm::usage =
"MultivectorNorm[m] gives the norm of a multivector m."
MultivectorNormSquared::usage =
"MultivectorNormSquared[m] gives the squared norm of a multivector m."

(*geometric operations*)
Angle::usage =
"Angle[a,b] gives the angle between two blades a and b. It permits general multivectors as arguments but may give meaningless answers for non-blades."

(*sets, not yet implemented*)
Multivectors::usage =
"Multivectors[] symbolically represents the set of all multivectors.
Multivectors[k] represents the set of k-vectors.
Multivectors[{grades}] represents the set of set of multivectors with a set of given grades.
Multivectors[grade(s), dim] represents a set of multivectors in 'dim' dimensions.
Multivectors[grade(s), {basis}] represents a set of multivectors with listed basis vectors."

Begin["`Private`"]

e/: pre___**e[i_]**e[j_]**post___/;Order[i,j]==-1:= -NonCommutativeMultiply[pre,e[j],e[i],post]
e/: pre__**e[i_]**e[i_]**post___:= Metric[i]*(pre**post)
e/: pre___**e[i_]**e[i_]**post__:= Metric[i]*(pre**post)
e/: e[i_]**e[i_]:= Metric[i]

e[n__]:= NonCommutativeMultiply@@e/@{n}/;Length[{n}]>1

SetAttributes[e,{HoldAll,NHoldAll}]

e::times = "Basis vectors multiplied using Times rather than NonCommutativeMultiply. Order may not have been preserved."
e/: e[i_]*e[j_] := (Message[e::times]; e[i]**e[j])

Metric[_]=1

number=_Integer|_Real|_Rational
NcmExpansionRules={
	a_Plus**b_:>(#**b&/@a),(*distributive rule*)
	a_**b_Plus:>(a**#&/@b),
	x___**(a:number)**y___:>a x**y,(*move scalars out of products*)
	x___**((a:number)y_)**z___:>a x**y**z,
	x___**a_?FreeQ[e|OverVector]**y___:>a x**y,(*here until I add $Assumptions support*)
	x___**(a_?FreeQ[e|OverVector] y_)**z___:>a x**y**z,
	NonCommutativeMultiply[a_]:>a,
	x_^n_Integer/;!FreeQ[x,e]&&n>1:>NonCommutativeMultiply@@Table[x,n]
};
ExpandNCM[expr_]:= expr//.NcmExpansionRules

KVectorPart[a_?(FreeQ[e|OverVector]),k_Integer]:= If[k===0,a,0]
KVectorPart[(a_:1)*(b:_e|_OverVector),k_Integer]:= If[k===1,a*b,0]/;FreeQ[a,e|OverVector]
KVectorPart[(a_:1)*(b:NonCommutativeMultiply[__e]),k_Integer]:= If[Length[b]===k,a*b,0]/;FreeQ[a,e|OverVector]
KVectorPart[(a_:1)*(b:NonCommutativeMultiply[(_e|_OverVector)..]),k_Integer]:= 0/;FreeQ[{a},e|OverVector]&&(Mod[Length[b]+1,2]==Mod[k,2]||k>Length[b])
KVectorPart[a_+b_,k_Integer]:= KVectorPart[a,k]+KVectorPart[b,k]

ScalarPart[a_]:= KVectorPart[a,0]

SetAttributes[KVectorPart,Listable]

MaxGrade[mv_]:= Length@DeleteDuplicates@Cases[mv,_e|_OverVector,{0,\[Infinity]}]

KVectorDecomposition[mv_]:= Table[KVectorPart[mv,n],{n,0,MaxGrade[mv]}]

Reversion[mv_]:= Sum[(-1)^k KVectorPart[mv,k],{k,0,MaxGrade[mv]}]
Involution[mv_]:= Sum[(-1)^(k (k-1)/2) KVectorPart[mv,k],{k,0,MaxGrade[mv]}]
Conjugation[mv_]:= Sum[(-1)^(k (k+1)/2) KVectorPart[mv,k],{k,0,MaxGrade[mv]}]

MultivectorNorm[mv_]:= Sqrt[KVectorPart[ExpandNCM[Reversion[mv]**mv],0]]
MultivectorNormSquared[mv_]:= KVectorPart[ExpandNCM[Reversion[mv]**mv],0]

GeometricProduct=ExpandNCM@*NonCommutativeMultiply

OuterProduct[a_,b_]:= ExpandNCM@Sum[KVectorPart[KVectorPart[a,j]**KVectorPart[b,k],j+k],{j,0,MaxGrade[a]},{k,0,MaxGrade[b]}]
CommutatorProduct[a_,b_]:= ExpandNCM[(a**b-b**a)/2]

LeftContraction[a_,b_]:= ExpandNCM@Sum[KVectorPart[KVectorPart[a,j]**KVectorPart[b,k],k-j],{j,0,MaxGrade[a]},{k,0,MaxGrade[b]}]
RightContraction[a_,b_]:= ExpandNCM@Sum[KVectorPart[KVectorPart[a,j]**KVectorPart[b,k],j-k],{j,0,MaxGrade[a]},{k,0,MaxGrade[b]}]
FatDotProduct[a_,b_]:= ExpandNCM@Sum[KVectorPart[KVectorPart[a,j]**KVectorPart[b,k],Abs[k-j]],{j,0,MaxGrade[a]},{k,0,MaxGrade[b]}]
ScalarProduct[a__]:= ExpandNCM@ScalarPart[GeometricProduct[a]]

Angle[a_,b_]:= ArcCos[MultivectorNorm[LeftContraction[a,b]]/(MultivectorNorm[a]MultivectorNorm[b])]

Multivectors/:{a__}\[Element]Multivectors[args___]:= Alternatives[a]\[Element]Multivectors[args]

(*for internal use only*)
(*assumedMultivectorQ[a_]:= Switch[Head@$Assumptions, And, , Element, MemberQ[$Assumptions[[1]], a]*)

End[]

EndPackage[]
