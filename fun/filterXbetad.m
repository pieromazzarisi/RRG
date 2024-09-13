function [Xt] = filterXbetad(Yt)

% Allocate memory
[N,~,T] = size(Yt);
Xt = nan(2*N-1,T);

for t = 1:T
    Y = Yt(:,:,t);
    X = betadMLE(Y);
    Xt(:,t) = X;
end
end