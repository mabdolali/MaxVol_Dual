clc
clear all
close all
addpath(genpath('library'));
%% Setting
set(0, 'DefaultAxesFontSize', 22);
set(0, 'DefaultLineLineWidth', 2);

m = 3; % dimension
r = m; % # of vertices
num_iter = 20;
num_experiments=50; % # of trials
purity = 0.76;
num_workers = 10;
for no= 1 : num_experiments
    Ni1 = 30*ones(r,1);%[500;100;30];%30*ones(r,1);%# of data points on each facet
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
    [v, West, theta, iter, Y, C, values1,values2] = maxvoldual_iter(M,r,1e2,Wg,num_workers);
    result1(1,no,:)=values1;
    result11(1,no,:)=values2;

    [m,N]=size(M);
    SNR = 60;
    varianc = sum(M(:).^2)/10^(SNR/10) /m/N ;
    n = sqrt(varianc)*randn([m N]);
    M = Mg + n;
    [v, West, theta, iter, Y, C, values1,values2] = maxvoldual_iter(M,r,10,Wg,num_workers);
    result1(2,no,:)=values1;
    result11(2,no,:)=values2;

    
    SNR = 40;
    [m,N]=size(M);
    varianc = sum(M(:).^2)/10^(SNR/10) /m/N ;
    n = sqrt(varianc)*randn([m N]);
    M = Mg + n;
    [v, West, theta, iter, Y, C, values1,values2] = maxvoldual_iter(M,r,1,Wg,num_workers);
    result(3,no,:)=values1;
    result11(3,no,:)=values2;
    
    
    SNR = 30;
    [m,N]=size(M);
    varianc = sum(M(:).^2)/10^(SNR/10) /m/N ;
    n = sqrt(varianc)*randn([m N]);
    M = Mg + n;
    [v, West, theta, iter, Y, C, values1,values2] = maxvoldual_iter(M,r,0.5,Wg,num_workers);
    result(4,no,:)=values1;
    result11(4,no,:)=values2;

end
result5=(mean(log2(result),2));
result5(isinf(result5))=0;

figure;
plot(1:num_iter,result5(1,:),'rO-',1:num_iter,result5(2,:),'b*--',1:num_iter,result5(3,:),'kd:',1:num_iter,result5(4,:),'ms-.');
legend('SNR=\infty','SNR=60','SNR=40','SNR=30');
xlabel('iteration');
ylabel('$log2(\frac{||W_{k}-W_{t}||_F}{||W_{t}||_F})$','interpreter','latex');
%ylabel('log_2(MRSA(W_{k},W_{t}))');

result2=(mean(log2(result11),2));

figure;
plot(1:num_iter,result2(1,:),'rO-',1:num_iter,result2(2,:),'b*--',1:num_iter,result2(3,:),'kd:',1:num_iter,result2(4,:),'ms-.');
legend('SNR=\infty','SNR=60','SNR=40','SNR=30');
xlabel('iteration');
ylabel('$log2(\frac{||v_{k}-v_{k-1}||_F}{||v_{k-1}||_F})$','interpreter','latex');
%ylabel('log_2(MRSA(W_{k},W_{t}))');