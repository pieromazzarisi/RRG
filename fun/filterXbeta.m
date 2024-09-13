function [Xt] = filterXbeta(Yt)

% Allocate memory
[N,~,T] = size(Yt);
Xt = nan(N,T);

for t = 1:T
    Y = Yt(:,:,t);
    X = betaMLE(Y);
    Xt(:,t) = X;
end
end