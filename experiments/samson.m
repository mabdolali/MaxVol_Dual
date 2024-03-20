% This code performs MaxVol_Dual on Samson database for hyperspectral unmixing
clc
clear all
close all
%% Setting
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultLineLineWidth', 2);
addpath(genpath('..'));
%% Load data
load('samson_1');
load('end3');
W = M; %endmembers
M = V; %input data
r = 3;
%% Run MV_Dual
X = M;
tic;
[v, W2, theta, iter] = maxvoldual(X,r,0.002);
%show results
W_est = max(W2,0);%*MAX;
H_est = FGMfcnls(M,W_est);
time = toc
Vaff = affichage(H_est',r,nRow,nCol,1);
sprintf("MRSA is %.2f", mrsa(W, matchCol(W_est,W)))
sprintf("Reconstruction err is %.2f", norm(M-W_est*H_est,'fro')/norm(M,'fro')*100)
% show spectral signatures
figure;plot(W_est(:,1)/max(W_est(:,1)),'--')
hold on;
plot(W_est(:,2)/max(W_est(:,2)),'-.')
plot(W_est(:,3)/max(W_est(:,3)))
legend('#1 Tree','#2 Soil','#3 Water');
xlabel('# bands');
ylabel('Reflectance');
grid on;
% show 2d projections
in.plotall = 1;
in.plotspec = {'k.','rx','r--'};figure;plot2d(M,W_est,in)