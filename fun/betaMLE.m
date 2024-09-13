function [X] = betaMLE(A)
%MLE_beta_model Undirected simple graphs
%   A is the adjacency matrix (symmetric), k the degree sequence
%%% control parameters
precision0 = 10^(-10);
max_ite = 10000;
%%% obs
[N,~] = size(A);
k = sum(A,2);
%%% output
X = nan(N,1);
%%% excluding 0-degree or N-1-degree (maximum) 
%%% because MLE does not exists in that case
k0 = k == 0;
kN = k == N-1;
ki = ~(k0 | kN);
%%% sistemare che p = 0 o p = 1 a seconda di k0 o kN
%%%
x0=0.5 +rand([N,1]);
theta0 = log(x0);
theta0(k0) = -1e3;
theta0(kN) =  1e2;
x0 = exp(theta0);
theta = theta0;
temp_pre=1;
temp_ite=1;
while (temp_pre>precision0) && (temp_ite<max_ite)
    matrix_g=(ones(N,1)*x0')./(1+x0*x0');
    matrix_g=matrix_g-diag(diag(matrix_g));
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
if temp_ite == max_ite
    disp('maximum number of iterations reached!');
    disp('check for non converging solution!');
end
X(ki) = theta0(ki);
end