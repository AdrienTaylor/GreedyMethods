clear all; clc;
% In this example, we use the factored SSEP-based gradient method for
% solving the smooth strongly convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x);
% where F(x) is smooth and strongly convex 
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of a subgradient method starting with an initial
% iterate satisfying ||x0-xs||<=1.

N     = 10; % we do not advise trying N>40 with pesto, but rather use the specialized code instead.
% for kappa = 100, N with saved values: N=1,...,10,12,14,16,18,21,25,28,33,38,43
% for other values, the data can be generated using the SSEP_files
kappa = 100; %
fileName = sprintf('../Data/Stepsizes_GFOM_N%d_kappa%d.mat',N,round(kappa));
load(fileName)


% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
param.L  = L;	% 'radius'-type constraint on the subgradient norms: ||g||<=M
param.mu = mu;
% F is the objective function
F = P.DeclareFunction('SmoothStronglyConvex',param); 


% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();            % x0 is some starting point
[xs,fs] = F.OptimalPoint();        % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1);% Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm and (4) performance measure
x = cell(N+1); y = cell(N+1); g = cell(N+1);
x{1} = x0; y{1} = x0;
for i=1:N
    [g{i},~] = F.oracle(x{i});
    y{i+1} = x{i} - 1/L * g{i};
    x{i+1} = y{i+1} + zeta(i) * (y{i+1} - y{i}) + eta(i) * (y{i+1} - x{i});
end
[g{N+1},fN]=F.oracle(x{N+1});
P.PerformanceMetric(fN-fs);

% (5) Solve the PEP
P.solve();

fprintf('Performance of GFOM L||x0-x*||^2/%5.3f \t performance of factored SSEP (saved): L||x0-x*||^2/%5.3f \t performance of factored SSEP (reproduced in PESTO): L||x0-x*||^2/%5.3f\n',1/wc_GFOM,1/wc_FactoredSSEP,1/double(fN-fs));
fprintf('Difference between PESTO and saved result: %5.3e\n',abs(double(fN-fs)-wc_SSEP));

