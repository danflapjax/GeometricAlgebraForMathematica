# GeometricAlgebraForMathematica
A Wolfram Language/Mathematica package for working with Geometric Algebra

To use this package, put it in your working directory and import it using:
```
<<GeometricAlgebra`
```
## The basics

Basis vectors are represented using `e[i]` where i is an index/label. In typical usage, i is a non-negative integer, e.g. `e[3]`, `e[0]`, but arbitrary symbols or expressions are permitted: `e[x]`, `e[∞]`, and `e[time]` would all be valid. Additionally, basis k-vectors may, for convenience, be given as sequences of arguments: `e[1,2]` (a basis bivector), `e[1,3,5]` (a basis trivector). These are expanded into geometric products of basis vectors.

Geometric products are represented using the built-in Wolfram Language function [NonCommutativeMultiply](https://reference.wolfram.com/language/ref/NonCommutativeMultiply.html), which has operator form `**`. These may then be expanded using `ExpandNCM`. Basis vectors automatically order themselves when multiplied: `e[2]**e[1]` becomes `-e[1]**e[2]`.

This library presently assumes orthogonality between all basis vectors. What each one squares to may be changed, however, by assigning values to `Metric[i]`, which may be arbitrary expressions. This can be useful for Spacetime Algebra, Projective Geometric Algebra, and working in alternative charts, manifolds, or spaces. For example:

```Mathematica
In[2]:= {e[0]**e[0], e[1]**e[1], e[2]**e[2]}
Out[2]= {1, 1, 1}

In[3]:= Metric[0]=0; Metric[1]=-1; Metric[2]=r;
In[4]:= {e[0]**e[0], e[1]**e[1], e[2]**e[2]}
Out[4]= {0, -1, r}
```

## Multivector decomposition

`KVectorPart[mv,k]` will extract the k-vector part of a multivector `mv`. `ScalarPart[mv]` will give the scalar part of a multivector `mv` (equivalent to `KVectorPart[mv,0]`). KVectorDecomposition is also provided which splits a multivector into a list of its parts:

```Mathematica
In[2]:= KVectorDecomposition[1+e[1]-e[2,3]+4e[1,2,3,4]]                                                                
Out[2]= {1, e[1], -e[2] ** e[3], 0, 4 e[1] ** e[2] ** e[3] ** e[4]}
```

## Products

This package includes the following products: GeometricProduct, OuterProduct, CommutatorProduct, LeftContraction, RightContraction, FatDotProduct, and ScalarProduct.

## Other operations

The package includes various involutions: Reversion, Involution, and Conjugation. MultivectorNorm and MultivectorNormSquared are provided for calculating norms. Angle computes the angle between two blades.

## Symbolic manipulation

Symbolic vectors may be represented using [OverVector](https://reference.wolfram.com/language/ref/OverVector.html), which renders an arrow over a symbol. These will behave properly as though they are vectors in functions like KVectorPart and ExpandNCM. I also plan on expanding symbolic capabilities further as I continue to work on this package.
