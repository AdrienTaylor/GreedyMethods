function [wc]=FixedSteps_SmoothStronglyConvex(h,N,L,mu,R,verb)

% Input:
%   - N:        Number of iterations N>0
%   - L:        Smoothness constant L>=0
%   - mu:       Strong convexity constant in the interval [0,L]
%   - R:        Initial condition: ||x_0-x_*||<=R with R>=0
%   - verb:     Let the solver talk? (Verbose mode) [0/1]
%
% Output:
%   - wc:       Worst-case guarantee on f(x_N)-f(x_*)
%
%   - h:        Step sizes. Their form is given such that
%               iterate x_i (i=0,...,N) is described by
%               0=P*(h(i,:)^T) with P=[x0 g0 ... gN ].
%
%   -beta:     dual variables in the text; corresponding constraint:
%               <g_i,g_j>==0 :  beta(i,j) (1<=j<=i=1,...N)
%
%   -gamma:     dual variables in the text; corresponding constraint:
%               <g_i,x_j-x_0>==0 :  gamma(i+1,j) (0<=j<i=1,...N)
%
% Usage:
%   N = 1; L = 1; mu = .1; R = 1; verb = 0;
%   [wc, h, beta, gamma]=GFOM_SmoothStronglyConvex(N,L,mu,R,verb)
%
% ** To reproduce the numerics from the paper, one needs to use mosek and
% to increase the default accuracy (see options for the solver below).
% ** If you don't have mosek installed, you change comment line #90 and
% comment out line #91

%% Initialize coordinates

dim              = N + 2; %[x0 g0 ... gN]
x                = [h; zeros(1,dim)];
g                = zeros(N+2,dim);
g(1:N+1,2:dim) = eye(N+1);

%% Write the PEP for GFOM

G = sdpvar(dim);
f = sdpvar(N+2,1);

cons       = (G>=0);
cons       = cons+(G(1,1)-R^2<=0);  % initial condition
cons_count = 1;               % count the number of constraints

% Interpolation constraints:
for i = 1:N+2
    for j = 1:N+2
        if (i ~= j)
            cons = cons+(f(j)-f(i)+g(j,:)*G*(x(i,:)-x(j,:)).'+...
                1/(2*(1-mu/L))*(1/L*(g(i,:)-g(j,:))*G*(g(i,:)-g(j,:)).'+...
                mu*(x(i,:)-x(j,:))*G*(x(i,:)-x(j,:)).'-...
                2*mu/L*(g(i,:)-g(j,:))*G*(x(i,:)-x(j,:)).')<=0);
            
            cons_count = cons_count+1;
        end
        
    end
end


% Set the objective function and solve the PEP
obj=f(N+1)-f(N+2); %f(xN)-f(x*)

% ** To reproduce the numerics from the paper:
ops = sdpsettings('verbose',verb,'solver','mosek','mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-10);
% ops = sdpsettings('verbose',verb);

optimize(cons,-obj,ops);

wc = double(obj);






















