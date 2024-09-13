function [Yt, Xt, Ft, Lambda0, P0] = simulateBeta(N,r,T,alpha,A,F0)
%%% factor model with Markov order 1 for the factor dynamics
% A needs to be diagonal
A = diag(diag(A));
Q = zeros(r);
for i = 1:r
    Q(i,i) = 1-A(i,i)^2;
end
% Allocate memory
Ft = nan(r,T+1);
Xt = nan(N,T);
Ft(:,1) = F0;
Yt = zeros(N,N,T);

% Factor loadings

%LambdaStar = 2+randn(N,r);
%LambdaStar = lognrnd(0,0.25,[N,r]);
LambdaStar = normrnd(0,0.4,[N,r]);
auxLambda = (LambdaStar')*LambdaStar;
[auxQ, auxD] = eig(auxLambda);
[~,indexSort] = sort(diag(auxD),'descend');
d0 =diag(auxD);
d0 = d0(indexSort);
D0 = diag(d0);
Lambda0 = LambdaStar*auxQ(:,indexSort);
P0 = Lambda0*diag(1./(diag(sqrt(D0))));
% Fundamental innovations
eta = sqrt(Q)*randn(r,T);

    for t = 1:T
        Xt(:,t) = alpha + Lambda0*Ft(:,t);
        Ft(:,t+1) = A*Ft(:,t) + eta(:,t);
        %%% simulating the network
        x = exp(Xt(:,t));
        p = (x*x')./(1+x*x');
        p = p -diag(diag(p));
        nets = binornd(1,p);
        nets = triu(nets) + triu(nets)';
        Yt(:,:,t) = nets;
        %%%
    end
    Ft(:,end) = [];
end