% This script allows reproducing results from Tables 1-2.
% (Illustration of the SSEP procedure for generating a
% fixed-step method for smooth strongly convex minimization)


%%  Parameters

N     = 10;       % Number of iterations
L     = 1;        % Lipschitz constant
kappa = 100;      % Condition ratio (can be set to 'kappa = Inf')
R     = 1;        % Initial distance
mu    = L/kappa;  % Strong convexity constant

verb  = 1;        % [0/1] Let the solver talk ?

%% Efficient form 1

[Algo, wc, err] = Efficient_form_SmoothStronglyConvex(R,mu,L,N,verb);

Algo.B.b1.', Algo.B.b2.', Algo.B.b3.'
Algo.C.c1.', Algo.C.c2.', Algo.C.c3.', Algo.C.c4.'

%% Efficient form 2

[Algo, wc, err] = Efficient_form2_SmoothStronglyConvex(R,mu,L,N,verb);
Algo.zeta.'
Algo.eta.'
