% This code performs MaxVol-Dual on Jasper Ridge database for hyperspectral unmixing
clc
clear all
close all

%% Load data
addpath(genpath('..'));
load('Urban');
r = 4;
%% Run MVDual
tic;
options.num_workers = 1;
[v, W2, theta, iter] = maxvoldual(X,r,0.0001);
% show results
W_est = max(W2,0);
% sprintf("MRSA is %.2f", mrsa(W, matchCol(W_est,W)))
% W_est = matchCol(W_est,W);
% H_est = FGMfcnls(M,W_est);
time = toc
% Vaff = affichage(H_est',r,nRow,nCol,1);
% sprintf("Reconstruction err is %.2f", norm(M-W_est*H_est,'fro')/norm(M,'fro')*100)

% show spectral signatures
% figure;plot(W_est(:,1))
% hold on;
% plot(W_est(:,2),'--')
% plot(W_est(:,3),'-.')
% plot(W_est(:,4)',':')
% legend('#1 Tree','#2 Water','#3 Soil','#4 Road');
% xlabel('# bands');
% ylabel('Reflectance');
% grid on;
% show 2d projection of data and convex hull
in.plotall = 1;
in.plotspec = {'k.','rx','r--'};figure;plot2d(X,W_est,in)