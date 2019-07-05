# Efficient First-order Methods for Convex Minimization: a Constructive Approach

This code can be used to reproduce the results from the work (on [arXiv](https://arxiv.org/abs/1803.05676), [Mathematical Programming](https://link.springer.com/article/10.1007/s10107-019-01410-2)):

> [1] Yoel Drori, Adrien B Taylor, "Efficient First-order Methods for Convex Minimization: a Constructive Approach",  	Mathematical Programming, 2019.

## Authors

- **Yoel Drori**
- [**Adrien Taylor**](http://www.di.ens.fr/~ataylor/)

## Getting started

To use the code, download the repository and add it to the path in Matlab.

**Notes:** 
- This code requires [YALMIP](https://yalmip.github.io/) along with a suitable SDP solver (e.g., Sedumi, SDPT3, Mosek).
- The files in the PESTO_files folder requires the installation of the Performance Estimation Toolbox ([PESTO](https://github.com/AdrienTaylor/Performance-Estimation-Toolbox)).

## List of files

(1) Subspace-search elimination procedure (SSEP):

- [`demo_script`](SSEP_files/demo_script.m) This script shows how to reproduce the results from Table 1 and Figure 1 in the paper (numerical computations of efficient formulations). This script also allows generating new Data/* files (saved SSEP step-sizes and worst-case guarantees), that can be used for validating the results with the pesto files.
- [`GFOM_SmoothStronglyConvex`](SSEP_files/GFOM_SmoothStronglyConvex.m) This function generates worst-case certificates for the GFOM in the smooth (possibly strongly) convex case, along with the SSEP-based fixed-step method.
- [`FactoredSSEP_SmoothStronglyConvex`](SSEP_files/FactoredSSEP_SmoothStronglyConvex.m) This function generates an efficiently-implementable (factored) method from the SSEP-based method for smooth (possibly strongly) convex minimization.
- [`FixedSteps_SmoothStronglyConvex`](SSEP_files/FixedSteps_SmoothStronglyConvex.m) This function generates worst-case certificates for fixed-step methods for smooth (possibly strongly) convex minimization. 

(2) PESTO_* files:

- [`pesto_GreedyFirstOrderMethod_NonSmooth`](PESTO_files/pesto_GreedyFirstOrderMethod_NonSmooth.m) Numerically computes the worst-case performance of the GFOM for non-smooth convex minimization.
- [`pesto_GreedyFirstOrderMethod_Smooth`](PESTO_files/pesto_GreedyFirstOrderMethod_Smooth.m) Numerically computes the worst-case performance of the GFOM for smooth (possibly strongly) convex minimization.
- [`pesto_GreedyFirstOrderMethod_SmoothStrCvx`](PESTO_files/pesto_GreedyFirstOrderMethod_SmoothStrCvx.m) Numerically computes the worst-case performance of the GFOM for smooth strongly convex minimization.
- [`pesto_SSEP_SubgradientMethod`](PESTO_files/pesto_SSEP_SubgradientMethod.m) Numerically computes the worst-case performance of the SSEP-based subgradient method for non-smooth convex minimization.
- [`pesto_SSEP_SmoothStronglyConvex.m`](PESTO_files/pesto_SSEP_SmoothStronglyConvex.m.m) Numerically computes the worst-case performance of the SSEP-based gradient method for smooth strongly convex minimization (using Data/* files, see below).
- [`pesto_FactoredSSEP_SmoothStronglyConvex.m`](PESTO_files/pesto_FactoredSSEP_SmoothStronglyConvex.m.m) Numerically computes the worst-case performance of the (factored) SSEP-based gradient method for smooth strongly convex minimization (using Data/* files, see below).
- [`pesto_OptimizedGradientMethod`](PESTO_files/pesto_OptimizedGradientMethod.m) Numerically computes the worst-case performance of OGM for smooth (possibly strongly) convex minimization.
- [`pesto_StronglyConvexFastGradientMethod`](PESTO_files/pesto_StronglyConvexFastGradientMethod.m) Numerically computes the worst-case performance of FGM for smooth strongly convex minimization.
- [`pesto_StronglyConvexTripleMomentumMethod`](PESTO_files/pesto_StronglyConvexTripleMomentumMethod.m) Numerically computes the worst-case performance of TMM for smooth strongly convex minimization.

(3) Data/* files: contains step-size and worst-case guarantees for a few values of N and kappa.

## Example

The following code is contained in [`demo_script.m`](SSEP_files/demo_script.m)

```Matlab
% This script allows reproducing results from Table 1 and Figure 1.
% (Illustration of the SSEP procedure for generating a
% fixed-step method for smooth strongly convex minimization)
clc; clear all;

N     = 10;     % Number of iterations
L     = 1;      % Lipschitz constant
kappa = 100;    % Condition number 1 < kappa <= Inf
R     = 1;      % Initial condition; i.e., ||x0-x*|| <= R
verb  = 0;      % Verbose solver ? [0/1]
saved = 0;      % Save data ? [0/1]

mu    = L/kappa;% Strong convexity (defined from condition number)
fprintf('Computing ... (a few seconds for N=10, possibly much longer when increasing N)\n')

[Algo, wc_GFOM, err, h] = FactoredSSEP_SmoothStronglyConvex(R,mu,L,N,verb);
fprintf('Computation for GFOM (1/3): DONE\n')
[wc_SSEP]          = FixedSteps_SmoothStronglyConvex(h,N,L,mu,R,verb);
fprintf('Computation for SSEP (2/3): DONE\n')
[wc_FactoredSSEP]  = FixedSteps_SmoothStronglyConvex(Algo.h,N,L,mu,R,verb);
fprintf('Computation for Factored SSEP (3/3): DONE\n')

fprintf(['Worst case guarantee for GFOM: L||x0-x*||^2/%5.2f \t for SSEP:'...
    'L||x0-x*||^2/%5.2f \t for factored SSEP: L||x0-x*||^2/%5.2f\n'],...
    1/wc_GFOM, 1/wc_SSEP, 1/wc_FactoredSSEP)
fprintf('Error max|h_ij - h_ij''|=%5.4e\n',max(max(abs(Algo.h-h))))
fprintf('Values for zeta_i (for 1<=i<=N in increasing i ordering):\n')
Algo.zeta.'
fprintf('Values for eta (for 1<=i<=N in increasing i ordering):\n')
Algo.eta.'
fprintf('Saving data\n')

if saved
    fileName = sprintf('../Data/Stepsizes_GFOM_N%d_kappa%d',N,round(kappa));
    eta  = Algo.eta;
    zeta = Algo.zeta;
    save(fileName,'L','mu','zeta','eta','h','wc_GFOM','wc_SSEP','wc_FactoredSSEP');
end





```

