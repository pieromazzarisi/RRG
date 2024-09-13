function [X] = betadMLE(A)
%MLE_beta_model Directed simple graphs
%   A is the adjacency matrix (symmetric), k the degree sequence
%%% control parameters
i_s = 3;
precision0 = 10^(-10);
max_ite = 10000;
%%% obs
[N,~] = size(A);
kout = sum(A,2);
%%%
kout(i_s) = 1/2;
%%%
kin = sum(A,1)';
k = [kout;kin];
[n,~] = size(k);
%%% output
X = nan(n,1);
%%% excluding 0-degree or N-1-degree (maximum) 
%%% because MLE does not exists in that case
k0 = k == 0;
kN = k == N-1;
%kS = zeros(n,1); kS(i_s) = 1; kS = logical(kS);
%ki = ~(k0 | kN | kS);
ki = ~(k0 | kN);
%%% sistemare che p = 0 o p = 1 a seconda di k0 o kN
%%%
x0=0.5 +rand([n,1]);
theta0 = log(x0);
theta0(k0) = -1e3;
theta0(kN) =  1e2;
x0 = exp(theta0);
theta = theta0;
temp_pre=1;
temp_ite=1;
while (temp_pre>precision0) && (temp_ite<max_ite)
    matrix_gout=(ones(N,1)*x0(N+1:n)')./(1+x0(1:N)*x0(N+1:n)');
    matrix_gout=matrix_gout-diag(diag(matrix_gout));
    matrix_gin=(ones(N,1)*x0(1:N)')./(1+x0(N+1:n)*x0(1:N)');
    matrix_gin=matrix_gin-diag(diag(matrix_gin));
    matrix_g = [matrix_gout;matrix_gin];
    temp_g = log(sum(matrix_g,2));
    theta(ki) = log(k(ki)) - temp_g(ki);
    x = exp(theta);
    g = (theta-theta0)./theta0;
    temp_pre=max(abs(g(ki)));
    temp_ite=temp_ite+1;
    x0=x;
    theta0 = theta;
    clear g matrix_g;
end
% if temp_ite == max_ite
%     disp('maximum number of iterations reached!');
% end
%%% define a global shift
% x = nan(n,1);
% x(ki) = theta0(ki);
% nn = sum(~isnan(x));
% xout = sum(x(1:N),'omitnan');
% xin = sum(x(N+1:n),'omitnan');
% a = (1/nn)*(xin-xout);
% x(1:N) = x(1:N) + a;
% x(N+1:n) = x(N+1:n) - a;
%%%
theta0(1:N) = theta0(1:N) - theta(i_s);
theta0(N+1:end) = theta0(N+1:end) + theta(i_s);
%%%
X(ki) = theta0(ki);
X(i_s) = [];
end