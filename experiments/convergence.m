% This code studies the convergence of MV-Dual. In particular, the evolution of 
% $\frac{\|v_{k}-v_{k-1}\|_2}{\|v_{k-1}\|_2}$ and $\frac{\|W_{k}-W_{t}\|_F}{\|W_{t}\|_F}$
% where $W_t$ is the ground truth, for different iterations.

clc
clear all
close all
addpath(genpath('..'));
%% Setting
set(0, 'DefaultAxesFontSize', 22);
set(0, 'DefaultLineLineWidth', 2);

m = 3; % dimension
r = m; % # of vertices
num_iter = 20; % # number of considered iterations
num_experiments=10; % # of trials
purity = 0.76;
num_workers = 10;
for no= 1 : num_experiments
    Ni1 = 30*ones(r,1);% use [500;100;30] for imbalanced case; # of data points on each facet
    Ni2 = 10; % # of data points within polytope
    while(true)
        [M, W, ~] = gendata_rnd(m,r,purity,Ni1,Ni2); %generating the data points
        if cond(W) <r*10 %limiting the condition number
            break;
        end
    end
    disp('finished generating data');
    Wg = W;
    Mg = M;
    z_init = randn(r,r);
    [values1,values2] = maxvoldual_iter(M,r,1e2,Wg,num_workers);
    result1(1,no,:)=values1;
    result11(1,no,:)=values2;

    [m,N]=size(M);
    SNR = 60;
    varianc = sum(M(:).^2)/10^(SNR/10) /m/N ;
    n = sqrt(varianc)*randn([m N]);
    M = Mg + n;
    [values1,values2] = maxvoldual_iter(M,r,10,Wg,num_workers);
    result1(2,no,:)=values1;
    result11(2,no,:)=values2;

    
    SNR = 40;
    [m,N]=size(M);
    varianc = sum(M(:).^2)/10^(SNR/10) /m/N ;
    n = sqrt(varianc)*randn([m N]);
    M = Mg + n;
    [values1,values2] = maxvoldual_iter(M,r,1,Wg,num_workers);
    result(3,no,:)=values1;
    result11(3,no,:)=values2;
    
    
    SNR = 30;
    [m,N]=size(M);
    varianc = sum(M(:).^2)/10^(SNR/10) /m/N ;
    n = sqrt(varianc)*randn([m N]);
    M = Mg + n;
    [values1,values2] = maxvoldual_iter(M,r,0.5,Wg,num_workers);
    result(4,no,:)=values1;
    result11(4,no,:)=values2;

end
figure;
result5=mean(result,2);
result5(result5==0)=eps;
semilogy(1:num_iter,result5(1,:),'rO-',1:num_iter,result5(2,:),'b*--',1:num_iter,result5(3,:),'kd:',1:num_iter,result5(4,:),'ms-.');
legend('SNR=\infty','SNR=60','SNR=40','SNR=30');
xlabel('iteration');
ylabel('$\frac{||W_{k}-W_{t}||_F}{||W_{t}||_F}$','interpreter','latex');

result2=mean(result11,2);
figure;
semilogy(1:num_iter,result2(1,:),'rO-',1:num_iter,result2(2,:),'b*--',1:num_iter,result2(3,:),'kd:',1:num_iter,result2(4,:),'ms-.');
legend('SNR=\infty','SNR=60','SNR=40','SNR=30');
xlabel('iteration');
ylabel('$\frac{||v_{k}-v_{k-1}||_F}{||v_{k-1}||_F}$','interpreter','latex');
