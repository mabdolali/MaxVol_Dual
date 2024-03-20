% This code compares MV-Dual with GFPI in noisy cases where facet-based
% approaches might not succeed
clc
clear all
close all
addpath(genpath('..'));
%% Setting
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultLineLineWidth', 2);
%% generate data
SNR = 30;
m = 3; % dimension
r = m; % # of endmembers
rng(59); 
purity = 0.8;
Ni1 = 50*ones(r,1); % # of points on each facets
Ni2 = 50; % # of points within polytope
while(true)
    [M, W] = gendata_rnd(m,r,purity,Ni1,Ni2); %generating the data points
    if cond(W) <r*10 %limiting the condition number
        break;
    end
end
disp("data generation finished");
%adding noise
[m,N]=size(M);
varianc = sum(M(:).^2)/10^(SNR/10) /m/N ;
n = sqrt(varianc)*randn([m N]);
M = M + n;
Wg = W;
%% Max vol dual
tic;
if SNR ==40
    lambda = 1;
elseif SNR == 50
    lambda = 5;
elseif SNR ==60
    lambda = 10;
elseif SNR ==30
    lambda = 0.5;
else
    lambda = 1e2;
end

[v, West, theta, iter] = maxvoldual(M,r,lambda);

%%
tic;
if SNR == 40
    vals = [10,0.5,0.1];
elseif SNR ==50
    vals = [10,0.5,0.05];
elseif SNR ==60
    vals = [100,0.5,0.01];
elseif SNR == 30
    vals = [1,0.5,0.1];
end
gfpi_options.lambda=vals(1);
gfpi_options.eta = vals(2); %margin
gfpi_options.gamma=vals(3); %safety gap
gfpi_options.no_show = true; % do not show intermediate results
gfpi_options.timelimit = 100; % timelimit of cplex
gfpi_options.centerstrategy = 'mean'; % center selection strategy
gfpi_options.outlier = false; % no consideration of outliers
% gfpi_options.MIPsolver = 'matlab';
W1 = GFPI(M,r,gfpi_options);


in.plotall = 1;
in.plotspec = {'k.','rx','r--'};figure;hold on;plot2d(M,W,in);
in.plotspec = {'k.','gd','g-'};plot2d(M,West,in)
in.plotspec = {'k.','bO','b:'};plot2d(M,W1,in)
legend('GT','MV-Dual','GFPI');

