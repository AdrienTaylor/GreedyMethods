% This script allows reproducing results from Table 1 and Figure 1.
% (Illustration of the SSEP procedure for generating a
% fixed-step method for smooth strongly convex minimization)
clc; clear all;

N     = 20;     % Number of iterations
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




