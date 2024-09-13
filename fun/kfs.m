function [Ff, Fs, Sf, Ss, logL] = kfs(Xt,N,r,T,Lambda,Psi,A,Q,Fpca,computeLogL)

Ff = nan(r,T+1); % +1 is the (out-of-sample) prediction
Fs = nan(r,T);
v = nan(N,T);
V = nan(N,N,T);
K = nan(r,N,T);
Sf = nan(r,r,T+1);
Ss = nan(r,r,T);
q = zeros(r,1);
M = zeros(r);

Ff(:,1) = Fpca(:,1);
v(:,1) = Xt(:,1)-Lambda*Ff(:,1);
Sf(:,:,1) = cov(Fpca');
V(:,:,1) = Lambda*Sf(:,:,1)*Lambda' + Psi;
K(:,:,1) = A*Sf(:,:,1)*Lambda'/V(:,:,1);


Ff(:,2) = A*Ff(:,1) + K(:,:,1)*v(:,1);
Sf(:,:,2) = A*Sf(:,:,1)*(A-K(:,:,1)*Lambda)' + Q;

% filtering

for t = 2:T
    
    v(:,t) = Xt(:,t)-Lambda*Ff(:,t);
    V(:,:,t) = Lambda*Sf(:,:,t)*Lambda' + Psi;
    K(:,:,t) = A*Sf(:,:,t)*Lambda'/V(:,:,t);

    Ff(:,t+1) = A*Ff(:,t) + K(:,:,t)*v(:,t);
    Sf(:,:,t+1) = A*Sf(:,:,t)*(A-K(:,:,t)*Lambda)' + Q;
    
end

% smoothing

for t = T:-1:1
    
    Lt = A - K(:,:,t)*Lambda;
    q = (Lambda')*(V(:,:,t)\v(:,t)) + Lt'*q;
    Fs(:,t) = Ff(:,t) + Sf(:,:,t)*q;
    
    M = (Lambda')*(V(:,:,t)\Lambda) + Lt'*M*Lt;
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
        
    
    