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
