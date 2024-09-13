function [A,A_s,X,X_s] = generate_time_series_beta(n,N,T,seed)
%%% function to generate beta graphs with same seed for latent dyn.
% n<=N
rng(seed);
ALPHA = zeros(N,1);%normrnd(0,0.5,[N,1]);
LAMBDA = normrnd(0,0.4,[N,1]);
a = 0.9;
Q = 1-a^2;
eta = sqrt(Q)*randn(1,T);
lambda = LAMBDA(1:n);
alpha = ALPHA(1:n);
alpha_sparse = alpha - log(n/10);
Ft = nan(1,T);
X = nan(n,T);
X_s = nan(n,T);
for t = 1:T
    if t == 1
        Ft(t) = eta(1);
        X(:,t) = alpha + lambda*Ft(t);
        X_s(:,t) = alpha_sparse + lambda*Ft(t);
    else
        Ft(t) = a*Ft(t-1) + eta(t);
        X(:,t) = alpha + lambda*Ft(t);
        X_s(:,t) = alpha_sparse + lambda*Ft(t);
    end  
end

%%% dense
rng('shuffle');
A = zeros(n,n,T);
for t = 1:T
    x = exp(X(:,t));
    p = (x*x')./(1+x*x');
    p = p -diag(diag(p));
    nets = binornd(1,p);
    nets = triu(nets) + triu(nets)';
    A(:,:,t) = nets;
end

%%% sparse
rng('shuffle');
A_s = zeros(n,n,T);
for t = 1:T
    x = exp(X_s(:,t));
    p = (x*x')./(1+x*x');
    p = p -diag(diag(p));
    nets = binornd(1,p);
    nets = triu(nets) + triu(nets)';
    A_s(:,:,t) = nets;
end
end

