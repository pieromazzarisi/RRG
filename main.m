%%% main Matlab implementation of RRG (published version)
%% this is an example to show how the code works
% Here, we show the case of the \beta-model as an example
% (the method is general, it is enough to change the function on ...
% the MLE estimation of the \beta-model to generalize to any ERG model.
%
% Y is the input of the algorithm as time series of network snapshots,
% i.e. [m,m,T] = size(Y) with m the number of nodes, T the sample size.
% 
% simulateBeta: function to simulate Y according to the \beta-model
% [Yt, Xt, Ft, Lambda0, P0] = simulateBeta(m,r,T,alpha,A,F0)
%   inputs:
%       m: # nodes (n # time varying parameters, m = n)
%       r: # factors
%       T: sample size
%       alpha: m x 1 vector of constants
%       A: diagonal r x r matrix of autoregressive coefficients
%       F0: r x 1 vector of initial points for the factor VAR(1) dynamics
%   outputs:
%       Y: m x m x T time series of network snapshots
%       Xt: n x T time series of latent time-varying parameters
%       Ft: r x T time series of unobserved factors
%       Lambda0: n x r matrix of factor loadings
%       P0: n x r matrix (needs for comparison in simulation framework)
%
% estBeta: function to estimate the dynamic factor \beta-model
%   inputs:
%       Y: m x m x T matrix, time series of network snapshots
%       r: # factors
%       input: structure of inputs (if any, see below)
%   outputs:
%       output: structure with the following fields:
%           Xt: n x T, single-snapshot MLE estimates
%           alpha: n x 1 vector, estimates of the constants
%           Fpca: r x T matrix, time series of latent factors estimated
%               with PCA
%           Fks: r x T matrix, time series of latent factors re-estimated
%               with Kalman smoother
%           Lambda: n x r matrix, estimated factor loadings
%           A: r x r matrix, estimated autoregressive coefficients
%           Q, Psi, P, D, matrices used in PCA
%           flag: output flag
%
%   WARNING:
%   The algorithm works with any other ERG model as long as the function
%   betaMLE.m is replaced in filterXbeta.m and estBeta.m 
%   by the corresponding single-snapshot MLE estimation method 
%   for the considered ERG model.
%      
%   e.g. betaMLE.m has input Y(:,:,t), i.e. m x m adjacency matrix, and
%   returns a n x 1 vector X of single-snapshot MLE estimate
%   - for the \beta model n = m
%   - but the method works in general for n \neq m
%%
clear all;
clc;
close all;
addpath(genpath(pwd));
%% simulate Y according to the \beta-model
m = 20;
r = 1;
T = 1000;
A = 0.9.*eye(r);
F0 = rand([r,1]);
%alpha = -1 + 2.*rand([m,1]);
alpha = zeros(m,1);
[Y, Xt, Ft, Lambda0, P0] = simulateBeta(m,r,T,alpha,A,F0);

%% input structure for the two-step approach (if any)
% such an input structure is only to define the same rotation for the
% factors in the simulation setting, for a direct comparison between the
% time series, otherwise leave empty input = [];

input = struct('Xt',Xt,'Ft',Ft,'alpha0',alpha,...
    'Lambda0',Lambda0,'P0',P0,'A',A,...
    'do1step',false,'doComparisonSim',true,'computeLogL',false);

%% Two-step approach for estimation: PCA + Kalman smoother

[output] = estBeta(Y,r,input);
