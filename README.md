# Efficient First-order Methods for Convex Minimization: a Constructive Approach

This code can be used to reproduce the results from the work:

> [1] Yoel Drori, Adrien B Taylor, "Efficient First-order Methods for Convex Minimization: a Constructive Approach", 2018.

## Getting started

To use the code, download the repository and add it to the path in Matlab.

**Notes:** 
- This code requires [YALMIP](https://yalmip.github.io/) along with a suitable SDP solver (e.g., Sedumi, SDPT3, Mosek).
- The files in the PESTO_files folder requires the installation of the Performance Estimation Toolbox ([PESTO](https://github.com/AdrienTaylor/Performance-Estimation-Toolbox)).

## List of files

- [`GreedyFirstOrderMethod_NonSmooth`](PESTO_files/pesto_GreedyFirstOrderMethod_NonSmooth.m) Numerically computes the worst-case performance of the GFOM for non-smooth convex minimization.
- [`SSEP_SubgradientMethod`](PESTO_files/pesto_SSEP_SubgradientMethod.m) Numerically computes the worst-case performance of the SSEP-based subgradient method for non-smooth convex minimization.
- [`GreedyFirstOrderMethod_Smooth`](PESTO_files/pesto_GreedyFirstOrderMethod_Smooth.m) Numerically computes the worst-case performance of the GFOM for smooth (possibly strongly) convex minimization.
- [`OptimizedGradientMethod`](PESTO_files/pesto_OptimizedGradientMethod.m) Numerically computes the worst-case performance of OGM for smooth (possibly strongly) convex minimization.
- [`StronglyConvexFastGradientMethod`](PESTO_files/pesto_StronglyConvexFastGradientMethod.m) Numerically computes the worst-case performance of FGM for smooth strongly convex minimization.
- [`StronglyConvexTripleMomentumMethod`](PESTO_files/pesto_StronglyConvexTripleMomentumMethod.m) Numerically computes the worst-case performance of TMM for smooth strongly convex minimization.

- [`demo_script`](SSEP_files/demo_script.m) This script show how to reproduce the results from Tables 1-2 in the paper (numerical computations of efficient formulations).
- [`GFOM_SmoothStronglyConvex`](SSEP_files/GFOM_SmoothStronglyConvex.m) This function generates worst-case certificates for the GFOM in the smooth (possibly strongly) convex case, along with the SSEP-base fixed-step method.
- [`Efficient_form_SmoothStronglyConvex`](SSEP_files/Efficient_form_SmoothStronglyConvex.m) This function generates an efficiently-implementable (efficient form 1) method from the SSEP-based method for smooth (possibly strongly) convex minimization.
- [`Efficient_form2_SmoothStronglyConvex`](SSEP_files/Efficient_form2_SmoothStronglyConvex.m) This function generates an efficiently-implementable (efficient form 2) method from the SSEP-based method for smooth (possibly strongly) convex minimization.

## Example

The following code is contained in [`demo_script.m`](SSEP_files/demo_script.m)

```Matlab
% This script allows reproducing results from Tables 1-2.
% (Illustration of the SSEP procedure for generating a
% fixed-step method for smooth strongly convex minimization)


%%  Parameters

N     = 10; 
L     = 1;
kappa = 100;
R     = 1;
verb  = 1;
mu    = L/kappa;


%% Efficient form 1

[Algo, wc, err] = Efficient_form_SmoothStronglyConvex(R,mu,L,N,verb);

Algo.B.b1.', Algo.B.b2.', Algo.B.b3.'
Algo.C.c1.', Algo.C.c2.', Algo.C.c3.', Algo.C.c4.'

%% Efficient form 2

[Algo, wc, err] = Efficient_form2_SmoothStronglyConvex(R,mu,L,N,verb);
Algo.zeta.'
Algo.eta.'
```

## Authors
- **Yoel Drori**
- [**Adrien Taylor**](http://www.di.ens.fr/~ataylor/)

