function [Algo, wc, err, h]=FactoredSSEP_SmoothStronglyConvex(R,mu,L,N,verb)

% Input:
%   - N:        Number of iterations N > 0
%   - mu:       Strong convexity constant in the interval [0,L]
%   - L:        Smoothness constant L >= 0
%   - R:        Initial condition: ||x_0-x_*|| <= R with R >= 0
%   - verb:     Let the solver talk? (Verbose mode) [0/1]

% Output: 
%   
%   - Algo.zeta: values of zeta for the factored SSEP method
%
%   - Algo.eta:  values of eta for the previous algorithm
%            y_i = x_{i-1}-1/L f'(x_{i-1})
%            x_i = y_i+Algo.zeta(i)*(y_i-y_{i-1})+Algo.eta(i)*(y_i-x_{i-1})
%
%   - Algo.h:    Step sizes for the canonical form of the fixed-step method
%                in "efficient form 2".              
%                They are such that :
%                  0 = P*(Algo.h(i,:)^T) with P=[x0 g0 ... gN ].
%
%   - wc:       Worst-case guarantee on f(x_N)-f(x_*) (inherited from the
%               analysis of GFOM)
%
%   - err:      Absolute error between the canonical representation of the
%               "efficient form 1" algorithm and the one in "efficient 
%               form 2":
%                 err = max_{i,j} [h_{i,j}(efficient form 1)-
%                                                h_{i,j}(efficient form 2)]
%


% Usage:
%   N = 10; L = 1; mu = .1; R = 1; verb = 0; 
%   [Algo, wc, err] = Efficient_form2_SmoothStronglyConvex(R,mu,L,N,verb)


%% Efficient form 2

% (1) Obtain a worst-case certificate along with the corresponding SSEP
% method
[wc, h, ~, ~] = GFOM_SmoothStronglyConvex(N, L, mu, R, verb);


% (2) Obtain values for zeta's and eta's using necessary conditions
zeta = zeros(N,1);
eta  = zeros(N,1);

hh = - h(2:end,2:end);
eta(1) = L*hh(1,1)-1;
for i = 1:N-1
    zeta(i+1) = (hh(i+1,i)-hh(i,i))/(hh(i,i)-1/L);
    eta(i+1)  = L*hh(i+1,i+1)-1-zeta(i+1);
end

% (3) Validate the method via the canonical form
x_reconstructed = zeros(N+1,N+2);
x_reconstructed(1,1) = 1; y_reconstructed(1,:) = x_reconstructed(1,:);

% Initialize coordinates
grads   = [zeros(N+1,1) eye(N+1)];

for i=1:N
    y_reconstructed(i+1,:) = x_reconstructed(i,:)-grads(i,:)/L;
    x_reconstructed(i+1,:) = y_reconstructed(i+1,:)+...
        zeta(i)*(y_reconstructed(i+1,:)-y_reconstructed(i,:))+...
        eta(i)*(y_reconstructed(i+1,:)-x_reconstructed(i,:));
end

err   = max(max(abs(h-x_reconstructed)));

% Output the canonical step sizes corresponding to the efficient form
Algo.h      = x_reconstructed;
Algo.eta    = eta;
Algo.zeta   = zeta;













