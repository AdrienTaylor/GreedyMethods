clear all; clc;
% In this example, we use the SSEP-based subgradient method for
% solving the non-smooth convex minimization problem
%   min_x F(x); for notational convenience we denote xs=argmin_x F(x);
% where F(x) satisfies a Lipschitz condition; i.e., it has a bounded
% gradient ||g||<=M for all g being a subgradient of F at some point.
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of a subgradient method starting with an initial
% iterate satisfying ||x0-xs||<=1.

% (0) Initialize an empty PEP
P=pep();

% (1) Set up the objective function
M=1;
param.R=M;	% 'radius'-type constraint on the subgradient norms: ||g||<=M

% F is the objective function
F=P.DeclareFunction('ConvexBoundedGradient',param); 

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();            % x0 is some starting point
[xs,fs]=F.OptimalPoint();        % xs is an optimal point, and fs=F(xs)
P.InitialCondition((x0-xs)^2<=1);% Add an initial condition ||x0-xs||^2<= 1

% (3) Algorithm and (4) performance measure
N=10; % number of iterations
x=x0;
y=x0;
d=Point('Point',0); %Initialize a vector of zeros

for i=1:N
    y=i/(i+1)*x+1/(i+1)*x0;
    [g,~]=F.oracle(x);
    d=(i*d+g)/(i+1);
    x=y-1/M/sqrt(N+1)*d;
end
[g,f]=F.oracle(x);
P.PerformanceMetric(f-fs);

% (5) Solve the PEP
P.solve();

% The result should be (and is) 1/sqrt(N+1).
[double(f-fs); 1/sqrt(N+1)]
