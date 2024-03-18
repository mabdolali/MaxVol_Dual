% This code performs GFPI on Jasper Ridge database for hyperspectral unmixing
clc
clear all
close all

%% Load data
addpath(genpath('library'));

load('jasperRidge2_R198');
load('end4');
W = M;
M = Y;
r = 4;
%% Run Greedy Facet-based Polytope Identification
% set parameters
%tic;
%MAX =  max(max(M));
%X = M / MAX ; %scale to [0,1]
X = M;
tic;
[v, W2, theta, iter] = maxvoldual(X,r,0.0015,5); %ghabl = 0.0015
% [Z,theta,Delta] = maxvol_dual_sparse_noisy(X,r,0.0005);
% show results
W_est = max(W2,0);
H_est = FGMfcnls(M,W_est);
time = toc
Vaff = affichage(H_est',r,nRow,nCol,1);
sprintf("MRSA is %.2f", mrsa(W, matchCol(W_est,W)))
sprintf("Reconstruction err is %.2f", norm(M-W_est*H_est,'fro')/norm(M,'fro')*100)
W_est = matchCol(W_est,W);
% show spectral signatures
figure;plot(W_est(:,1))
hold on;
plot(W_est(:,2),'--')
plot(W_est(:,3),'-.')
plot(W_est(:,4)',':')
legend('#1 Water','#2 Soil','#3 Road','#4 Tree');
xlabel('# bands');
ylabel('Reflectance');
grid on;
% show 2d projection of data and convex hull
in.plotall = 1;
in.plotspec = {'k.','rx','r--'};figure;plot2d(M,W_est,in)