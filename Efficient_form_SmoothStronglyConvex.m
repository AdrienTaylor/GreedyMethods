function [Algo, wc, err]=Efficient_form_SmoothStronglyConvex(R,mu,L,N,verb)

% Input:
%   - R:        Initial condition: ||x_0-x_*|| <= R with R >= 0
%   - mu:       Strong convexity constant in the interval [0,L]
%   - L:        Smoothness constant L >= 0
%   - N:        Number of iterations N > 0
%   - verb:     Let the solver talk? (Verbose mode) [0/1]
%
% Output:
%   
%   - Algo.B    Algo.B.bj(i) provides the values of b^{(j)}_i for the 
%               algorithm in the "efficient form 1".
%
%   - Algo.C    Algo.C.cj(i) provides the values of c^{(j)}_i for the 
%               algorithm in the "efficient form 1".
%
%   - Algo.h:   Step sizes for the canonical form of the fixed-step method
%               in "efficient form 1".              
%               They are such that :
%                   0=P*(Algo.h(i,:)^T) with P=[x0 g0 ... gN ].
%
%   - wc:       Worst-case guarantee on f(x_N)-f(x_*) (inherited from the
%               analysis of GFOM)
%
%   - err:      Absolute error between the canonical representation of the
%               SSEP-based algorithm and the one in "efficient form 1":
%               err = max_{i,j} [h_{i,j}(SSEP)- h_{i,j}(efficient form 1)]
%
% Usage:
%   N = 10; L = 1; mu = .1; R = 1; verb = 0; 
%   [Algo, wc, err] = Efficient_form_SmoothStronglyConvex(R,mu,L,N,verb)


%% Efficient form 1

% (1) Obtain worst-case certificate from GFOM using PEP
[wc, h, beta, gamma] = GFOM_SmoothStronglyConvex(N, L, mu, R, verb);

% (2) Construction of b1, b2, b3
beta = beta(2:end, :);
gamma = gamma(2:end, 2:end);
b1(:,1) = [0; beta(2:N,1)];
b2(:,1) = [(beta(N,1:N-1)./beta(N,1)).'; 0];

beta_rec = tril(b1*b2.', 0);
b3 = diag(beta - beta_rec);


% (3) Construction of c1, c2, c3, c4
c1 = gamma(1:N, 1);
c2 = [(gamma(N, 1:N-1)./gamma(N, 1)) 0].';
c1 = c1*c2(end-1, 1);  % Normalization (we fix the last nonzero c2 to 1)
c2 = c2/c2(end-1, 1);  % Normalization (we fix the last nonzero c2 to 1)
gamma_rec = c1*c2.';
gamma_rec = tril(gamma_rec, -1);

c3 = diag(gamma - gamma_rec, -1);
c3 =[ 0; c3];
c4 = diag(gamma);

% (4) Validate using canonical forms
h_reconstructed = zeros(N+1,N+2);
h_reconstructed(1,1) = 1;
g = zeros(N+1,N+2);
g(:,2:N+2)=eye(N+1);

for i=1:N
    if i>1
        yi = c2(1:i-1).'*(h_reconstructed(2:i,:)-...
            ones(i-1,1)*h_reconstructed(1,:));
    else
        yi = zeros(1,N+2);
    end
    di                      = b2(1:i).'*g(1:i,:);
    h_reconstructed(i+1,:)  = h_reconstructed(1,:)-...
        1/c4(i)*(c1(i)*yi+c3(i)*(h_reconstructed(i,:)-h_reconstructed(1,:))+...
        b1(i)*di+b3(i)*g(i,:));
end

err = max(max(abs(h_reconstructed-h)));
% Output the canonical step sizes corresponding to the efficient form
Algo.h = h_reconstructed;

Algo.B.b1 = b1;
Algo.B.b2 = b2;
Algo.B.b3 = b3;

Algo.C.c1 = c1;
Algo.C.c2 = c2;
Algo.C.c3 = c3;
Algo.C.c4 = c4;
