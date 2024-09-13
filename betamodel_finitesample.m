%%% Script to obtain the results for the beta model in finite sample
%%% in the case of dense graphs
%% inputs
clear all;
clc;
close all;
addpath(genpath(pwd));
%%%%
%file = 'MLE_beta_finite_sample.mat';
%%%
%% simulations
vec_n = ceil(logspace(1,2.3,10));
vec_n(6) = 50;
vec_n(8) = 100;
N = max(vec_n);
S = 100;
T = 400;
Nlags = 20;
%
errori_cell = cell(length(vec_n),S);
Mphi = nan(length(vec_n),2,S);

%
for s = 1:S
    %disp(s);
    %tic
    seed = s;
    for i = 1:length(vec_n)
        n = vec_n(i);   
        [A,~,X,~] = generate_time_series_beta(n,N,T,seed);
        %%% estimation of the beta model
        % errors
        [Y] = filterXbeta(A);
        errors = (Y-X);
        errori_cell{i,s} = errors(:);
        % scaling variance
        try
            Ve = nan(n,1);
            for j = 1:n
                Ve(j) = var(errors(j,:),'omitnan');
            end
            Mphi(i,1,s) = max(Ve);
            Mphi(i,2,s) = min(Ve);
        catch
            disp('Error variance noise');
        end
    end
    %toc
    %save(file);
end
%save(file);
%% plotting errors

figure
for i = [1,3,6,8,10]
    [M,~] = size(errori_cell{i,1});
    pts = nan(M,S);
    for s = 1:S
        pts(:,s) = errori_cell{i,s};
    end
    pts = pts(:);
    [f,x] = ksdensity(pts,'Bandwidth',0.25);
    plot(x,f,'LineWidth',1);
    if i == 1
        hold on
    end
end
xlabel('estimation error','interpreter','latex')
set(gca,'Fontsize',15,'yscale','log','xscale','lin');
xlim([-4,4]);
ylim([0.001,5]);
legend('$n=10$','$n=20$','$n=50$','$n=100$','$n=200$',...
'interpreter','latex','FontSize',15);
title('Dense regime, finite number of nodes $n$','interpreter','latex');
hold off

%% scaling of the variance of errors
figure
mD = mean(Mphi,3,'omitnan');

squareRootEq = 'a*x^(-1)';
fDmax = fit(vec_n(4:end)',mD(4:end,1),squareRootEq);
fDmin = fit(vec_n(4:end)',mD(4:end,2),squareRootEq);

plot(vec_n,mD(:,1),'ko-')
hold on
plot(vec_n,mD(:,2),'bo-')
set(gca,'Fontsize',15,'yscale','log','xscale','log');
plot(vec_n,fDmax(vec_n),'r--')
plot(vec_n,fDmin(vec_n),'r--')
xlabel('number of nodes','interpreter','latex')
title('Scaling of $\lambda_{max}[\Phi_\epsilon(0)]$ and $\lambda_{min}[\Phi_\epsilon(0)]$ for the $p_1$-model','interpreter','latex');
legend('$\lambda_{max}$ (dense)','$\lambda_{min}$ (dense)','$\lambda \propto 1/n$','interpreter','latex','FontSize',15)


