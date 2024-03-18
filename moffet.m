% This code performs GFPI on Samson database for hyperspectral unmixing
clc
clear all
close all
%% Setting
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultLineLineWidth', 2);
addpath(genpath('library'));
%% Load data
load('Moffet');
% load('end3');

% W = M; %endmembers
M = X; %input data

r = 3;
%% Run Greedy Facet-based Polytope Identification
tic;

MAX =  max(max(M));
X = M / MAX ; %scale to [0,1]
[v, W2, theta, iter] = maxvoldual(X,r,0.01);
time = toc
%show results
W_est = max(W2,0)*MAX;
H_est = FGMfcnls(M,W_est);
nRow = 50;
nCol = 50;
Vaff = affichage(H_est',r,nRow,nCol,1);
% sprintf("MRSA is %.2f", mrsa(W, matchCol(W_est,W)))
% sprintf("Reconstruction err is %.2f", norm(M-W_est*H_est,'fro')/norm(M,'fro')*100)

% show spectral signatures
figure;plot(W_est(:,1)/max(W_est(:,1)),'--')
hold on;
plot(W_est(:,2)/max(W_est(:,2)),'-.')
plot(W_est(:,3)/max(W_est(:,3)))
legend('#1 Water','#2 Soil','#3 Vegetation');
xlabel('# bands');
ylabel('Reflectance');
grid on;

% show 2d projections
in.plotall = 1;
in.plotspec = {'k.','rx','r--'};figure;plot2d(M,W_est,in)
% % %% Other Approaches
% figure;
% [J,~]=SNPA(M,r);
% W_est = M(:,J);
% H_est = FGMfcnls(M,W_est);
% sprintf("MRSA of SNPA is %.2f", mrsa(W, matchCol(W_est,W)))
% sprintf("Reconstruction err of SNPA is %.2f", norm(M-W_est*H_est,'fro')/norm(M,'fro')*100)
% in.plotall = 1;
% in.plotspec = {'k.','rx','r--'};figure;plot2d(M,W_est,in)
% 
% [W_est, ~, ~] = HyperCSI(M,r);
% H_est = FGMfcnls(M,W_est);
% sprintf("Reconstruction err of HyperCSI is %.2f", norm(M-W_est*H_est,'fro')/norm(M,'fro')*100)
% sprintf("MRSA of HyperCSI is %.2f", mrsa(W, matchCol(W_est,W)))