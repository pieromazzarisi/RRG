function [AL,Q] = olsVAR(F,orderVAR)
[T,~] = size(F);
if T == 1
    F = F';
end
y = F;
X = zeros(T,orderVAR);
for i = 1:orderVAR
    X(:,i) = circshift(F,i);
end
y(1:orderVAR) = [];
X(1:orderVAR,:) = [];
% X = F(1:end-1);
% y = F(2:end);

[AL,~,r] = regress(y,X);
Q = var(r);
end