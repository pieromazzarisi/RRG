function [Ff, Fs, Sf, Ss, logL] = kfs_includenan(Xt,N,r,T,Lambda,Psi,A,Q,Fpca,computeLogL)
% TO DO
Ff = nan(r,T+1);
Fs = nan(r,T);
v = zeros(N,T);
V = zeros(N,N,T);
K = zeros(r,N,T);
Sf = nan(r,r,T+1);
Ss = nan(r,r,T);
q = zeros(r,1);
M = zeros(r);

%%% approaching missing values
Lambda0 = Lambda;
Psi0 = Psi;
Wt = cell(T,1);
It = cell(T,1);
In = eye(N);
for t = 1:T
    temp_nan = isnan(Xt(:,t));
    It{t} = ~temp_nan;
    Wt{t} = In(~temp_nan,:);
end
Xt(isnan(Xt)) = 0; % in order to avoid sum with 1*Xt(i,j)+ 0*nan = nan
%%%

y = Wt{1}*Xt(:,1);
Lambda = Wt{1}*Lambda0;
Psi = Wt{1}*Psi0*Wt{1}';

Ff(:,1) = Fpca(:,1);
v(It{1},1) = y-Lambda*Ff(:,1);
Sf(:,:,1) = cov(Fpca');
V(It{1},It{1},1) = Lambda*Sf(:,:,1)*Lambda' + Psi;
K(:,It{1},1) = A*Sf(:,:,1)*Lambda'/V(It{1},It{1},1);


Ff(:,2) = A*Ff(:,1) + K(:,It{1},1)*v(It{1},1);
Sf(:,:,2) = A*Sf(:,:,1)*(A-K(:,It{1},1)*Lambda)' + Q;

% filtering

for t = 2:T
    y = Wt{t}*Xt(:,t);
    Lambda = Wt{t}*Lambda0;
    Psi = Wt{t}*Psi0*Wt{t}';
    
    v(It{t},t) = y-Lambda*Ff(:,t);
    V(It{t},It{t},t) = Lambda*Sf(:,:,t)*Lambda' + Psi;
    K(:,It{t},t) = A*Sf(:,:,t)*Lambda'/V(It{t},It{t},t);

    Ff(:,t+1) = A*Ff(:,t) + K(:,It{t},t)*v(It{t},t);
    Sf(:,:,t+1) = A*Sf(:,:,t)*(A-K(:,It{t},t)*Lambda)' + Q;
    
end

% smoothing

for t = T:-1:1
    
    Lambda = Wt{t}*Lambda0;
    
    Lt = A - K(:,It{t},t)*Lambda;
    q = (Lambda')*(V(It{t},It{t},t)\v(It{t},t)) + Lt'*q;
    Fs(:,t) = Ff(:,t) + Sf(:,:,t)*q;
    
    M = (Lambda')*(V(It{t},It{t},t)\Lambda) + Lt'*M*Lt;
    Ss(:,:,t) = Sf(:,:,t) - Sf(:,:,t)*M*Sf(:,:,t);
end

% log-likelihood

if computeLogL
    logL = 0;
    for t = 1:T
        Ft = squeeze(V(:,:,t));
        vt = v(:,t);
        logL = logL - log(det(Ft)) - vt'*(Ft\vt);
    end
    logL = .5*logL-.5*T*N*log(2*pi);
else
    logL = nan;
end
end
        
    
    