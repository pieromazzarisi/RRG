function [X] = p1MLE(A)
%MLE_p1_model Directed simple graphs
%   A is the adjacency matrix (asymmetric)
%   theta(1) is put equal to zero for normalization
%%% control parameters
precision0 = 10^(-10);
max_ite = 10000;
%%% obs
[N,~] = size(A);
indA = (A == A') & (A == 1);
As = A;
As(indA) = 0;
kout = sum(As,2);
kin = sum(As,1)';
krec = sum(indA,2);
k = [kout;kin;krec];
%%% output
X = nan(3*N,1);
%%% excluding 0-degree or N-1-degree (maximum) 
%%% because MLE does not exists in that case
k0out = kout == 0;
k0in = kin == 0;
k0rec = krec == 0;
kNout = kout == N-1;
kNin = kin ==N-1;
kNrec = krec == N-1;
kiout = ~(k0out | kNout);
kiin = ~(k0in | kNin);
kirec = ~(k0rec | kNrec);
k0 = [k0out;k0in;k0rec];
kN = [kNout;kNin;kNrec];
ki = [kiout;kiin;kirec];
%%% sistemare che p = 0 o p = 1 a seconda di k0 o kN
%%%
x0=0.5 +rand([3*N,1]);
theta0 = log(x0);
theta0(k0) = -1e3;
theta0(kN) =  1e2;
x0 = exp(theta0);
temp_pre=1;
temp_ite=1;
while (temp_pre>precision0) && (temp_ite<max_ite)
    gout = (ones(N,1)*x0(N+1:2*N)')./(1+x0(1:N)*x0(N+1:2*N)'+...
        x0(N+1:2*N)*x0(1:N)'+x0(2*N+1:end)*x0(2*N+1:end)');
    gout = gout - diag(diag(gout));
    gin = (ones(N,1)*x0(1:N)')./(1+x0(1:N)*x0(N+1:2*N)'+...
        x0(N+1:2*N)*x0(1:N)'+x0(2*N+1:end)*x0(2*N+1:end)');
    gin = gin - diag(diag(gin));
    grec = (ones(N,1)*x0(2*N+1:end)')./(1+x0(1:N)*x0(N+1:2*N)'+...
        x0(N+1:2*N)*x0(1:N)'+x0(2*N+1:end)*x0(2*N+1:end)');
    grec = grec - diag(diag(grec));
    g = [gout;gin;grec];
    x = k./sum(g,2);
    theta = log(x);
    prec = (theta-theta0)./theta0;
    temp_pre=max(abs(prec(ki)));
    temp_ite=temp_ite+1;
    x0=x;
    theta0 = theta;
    clear gout gin grec g prec;
end
% if temp_ite == max_ite
%     disp('maximum number of iterations reached!');
% end
%%% NORMALIZATION
theta0(1:N) = theta0(1:N) - theta(1);
theta0(N+1:2*N) = theta0(N+1:2*N) + theta(1);
X(ki) = theta0(ki);
end