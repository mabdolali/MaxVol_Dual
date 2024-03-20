% This code compares facet-based criterion used in GFPI with the
% volume-based MV-Dual for a simple example in 2 dimensions. When the
% purity is low, GFPI inds the correct endmembers, whereas the volume-based
% MV-Dual selects the other enclosing simplex with smaller volume 
clc
clear all
close all
addpath(genpath('..'));
%% Setting
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultLineLineWidth', 2);
%% Parameters
%% generate data

m = 3; % dimension
r = m; % # of endmembers

num_experiments=1; % # of trials
startp = (1/(r-1)+0.01); % starting purity value
endp = 1; % ending purity value
step = (endp-startp)/6;

ind = 0;
purity = 0.6;
Ni1 = 50*ones(r,1); % # of points on each facets
Ni2 = 20; % # of points within polytope
while(true)
    [M, W] = gendata_rnd(m,r,purity,Ni1,Ni2); %generating the data points
    if cond(W) <r*10 %limiting the condition number
        break;
    end
end
disp("data generation finished");
%adding noise
[m,N]=size(M);
Wg = W;
%% Max vol dual
[v, West, theta, iter] = maxvoldual(M,r,1e2);

%% GFPI
vals = [100,0.5,0.01];
gfpi_options.lambda=vals(1);
gfpi_options.eta = vals(2); %margin
gfpi_options.gamma=vals(3); %safety gap
gfpi_options.no_show = true; % do not show intermediate results
gfpi_options.timelimit = 100; % timelimit of cplex
gfpi_options.centerstrategy = 'mean'; % center selection strategy
gfpi_options.outlier = false; % no consideration of outliers
W1 = GFPI(M,r,gfpi_options);

%% Compare results
in.plotall = 1;
in.plotspec = {'k.','rx','r--'};figure;hold on;plot2d(M,W,in);
in.plotspec = {'k.','gd','g-'};plot2d(M,West,in)
in.plotspec = {'k.','bO','b:'};plot2d(M,W1,in)
legend('GT','MV-Dual','GFPI');
axis equal

