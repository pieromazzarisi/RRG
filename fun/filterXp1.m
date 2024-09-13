function [Xt] = filterXp1(Yt)

% Allocate memory
[N,~,T] = size(Yt);
Xt = nan(3*N,T);

for t = 1:T
    Y = Yt(:,:,t);
    X = p1MLE(Y);
    Xt(:,t) = X;
end
end