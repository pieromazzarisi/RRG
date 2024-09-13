function [output] = estBeta(Yt,r,input)
%ESTIMATION_BETA_MODEL 
%%% OUTPUT
output = struct('Xt',[],'alpha',[],'Fpca',[],'Fks',[],...
    'Lambda',[],'A',[],'Q',[],'Psi',[],...
    'P',[],'D',[],'flag',[]);
%%% INPUT
if nargin == 3
    do1step = input.do1step;
    doComparisonSim = input.doComparisonSim;
    computeLogL = input.computeLogL;
else
    do1step = false;
    doComparisonSim = false;
    computeLogL = false;
end
%%%
[N,~,T] = size(Yt);
Xt = filterXbeta(Yt);
alpha = mean(Xt,2,'omitnan');
Xt = Xt-alpha;
output.Xt = Xt;
output.alpha = alpha;
%Ir = eye(r);

%% 1-step: PCA
try
% C = cov(Xt');
C = cov(Xt','partialrows');
[P,D] = eig(C);
[~,indexSort] = sort(diag(D),'descend');
auxD = diag(D);
D = diag(auxD(indexSort(1:r)));
P = P(:,indexSort(1:r));
sqrtD = sqrt(D);
if doComparisonSim
    tempSign = diag(P'*input.P0);
    for i = 1:length(tempSign)
        P(:,i) = sign(tempSign(i)).*P(:,i);
    end
end

Xtf = Xt;
Xtf(isnan(Xt(:,1)),1) = 0;
temp_nan = isnan(Xtf);
while any(temp_nan(:))
    Xtf(temp_nan) = Xtf(circshift(temp_nan,-1,2));
    temp_nan = isnan(Xtf);
end
Fpca = diag(1./(diag(sqrtD)))*P'*Xtf;
Lambda = P*sqrtD;

output.Fpca = Fpca;
output.Lambda = Lambda;
output.P = P;
output.D = D;

% Estimation of VAR(1) coefficients
% we consider pVAR = 1 for all implementation at the moment
A = zeros(r);
Q = zeros(r);
for j=1:r
    [auxA, auxQ] = olsVAR(Fpca(j,:)',1);
    A(j,j) = auxA;
    Q(j,j) = auxQ;
end
output.A = A;
output.Q = Q;
%%% cov matrix of observation noise
Psi = C - Lambda*Lambda'; % called R before
output.Psi = Psi;

flag = 'First-step PCA estimation';
output.flag = flag;
catch
    disp('Error in PCA (1-step)');
    flag = 'Error in the PCA step';
    output.flag = flag;
end
%% 2-step: Kalman Smoothing (population results)
if do1step
    return
end
try
Psi = frobproj(Psi,1e-10);
psiHat = diag(diag(Psi));

if any(isnan(Xt(:)))
    [~,Fks,~,~,~] = kfs_includenan(Xt,N,r,T,Lambda,psiHat,A,Q,Fpca,computeLogL);
else
    [~,Fks,~,~,~] = kfs(Xt,N,r,T,Lambda,psiHat,A,Q,Fpca,computeLogL);
end
output.Fks = Fks;
flag = 'Two-step PCA + KS estimation';
output.flag = flag;
catch
    disp('Error in KS (2-step)');
    flag = 'Error in the KS step';
    output.flag = flag;
end
end