% This code compares the performance of MV-Dual with GFPI, Min-Vol, HyperCSI, SNPA and MVIE in the noiseless case for the full rank
% matrices. 
clc
clear all
close all
addpath(genpath('..'));
%% Setting
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultLineLineWidth', 2);
%% generate data & initialization
m = 3; % dimension
r = m; % # of vertices
num_experiments=10; % # of trials
startp = (1/(r-1))+0.01; % starting purity
endp = 1; % ending purity
step = (endp-startp)/6;
range_purity = startp:step:endp; %define purity range
time = zeros(num_experiments,8,length(range_purity));
result = zeros(num_experiments,8, length(range_purity));
ind = 0;
for purity = range_purity
    ind = ind +1;
    for no= 1 : num_experiments
        Ni1 = 30*ones(r,1); % # of data points on each facet
        Ni2 = 10; % # of data points within polytope
        while(true)
            [M, W, ~] = gendata_rnd(m,r,purity,Ni1,Ni2); %generating the data points
            if cond(W) <r*10 %limiting the condition number
                break;
            end
        end
        disp('finished generating data');
        Wg = W;
        %% max vol dual
        tic;
        [v, West, theta, iter] = maxvoldual(M,r,1e2);
        result(no,1,ind)=compareWs(Wg, West);
        time(no,1,ind) = toc;
        %% GFPI
        tic;
        gfpi_options.lambda=1000;
        gfpi_options.eta = 0.5; % margin
        gfpi_options.gamma=0.001; %saftey gap
        gfpi_options.no_show = true; % do not display intermediate results
        gfpi_options.timelimit = 1000; % cplex timelimit
        gfpi_options.outlier = false; % no consideration of outliers
        % gfpi_options.MIPsolver = 'matlab';
        West = GFPI(M,r,gfpi_options);
        result(no,2,ind)=compareWs(Wg, West);
        time(no,2,ind) = toc;
        %%
        tic;
        options.lambda=0.1;
        [W2,H,e,er1,er2,lambda] = minvolNMF(M,r,options);
        result(no,3,ind)=compareWs( Wg, W2);
        time(no,3,ind) = toc;
        %%
        tic;
        options.lambda=1;
        [W2,H,e,er1,er2,lambda] = minvolNMF(M,r,options);
        result(no,4,ind)=compareWs( Wg, W2);
        time(no,4,ind) = toc;
        %%
        tic;
        options.lambda=5;
        [W2,H,e,er1,er2,lambda] = minvolNMF(M,r,options);
        result(no,5,ind)=compareWs( Wg, W2);
        time(no,5,ind) = toc;
        %% SNPA
        tic;
        K = SNPA(M,r);
        W3 = M(:,K);
        result(no,6,ind)=compareWs( Wg, W3);
        time(no,6,ind) = toc;
        %% MVIE
        tic;
        W4 = MVIE(M,r);
        result(no,7,ind)=compareWs( Wg, W4);
        time(no,7,ind) = toc;
        %% HyperCSI
        tic;
        [W5, ~, ~] = HyperCSI(M,r);
        result(no,8,ind)=compareWs( Wg, W5);
        time(no,8,ind) = toc;
        %% MVES
        tic;
        [W6,S_est,iter_cnt] = MVES(M,r,1);
        result(no,9,ind)=compareWs( Wg, W6);
        time(no,9,ind) = toc;
    end
end
result2=zeros(size(result,2),length(range_purity));
result2(:,:) = mean(result,1);
figure;
plot(range_purity,result2(1,:),'--*',range_purity,result2(2,:),'--O',range_purity,result2(3,:),'--+',range_purity,result2(4,:),'--<',range_purity,result2(5,:),'--d',range_purity,result2(6,:),'-->',range_purity,result2(7,:),'--s',range_purity,result2(8,:),'--x',range_purity,result2(9,:),'--^','markersize',10);
xlabel('purity');
ylabel('|| W - W_t || / || W_t ||');
% title(strcat('r=',num2str(r),',m=',num2str(m)));
legend('MV-Dual','GFPI','min vol \lambda = 0.1','min vol \lambda = 1','min vol \lambda = 5','SNPA','MVIE','HyperCSI','MVES');
xlim([startp,endp]);
xticks(range_purity);
xticklabels(range_purity);
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
%%
result2=zeros(size(time,2),length(range_purity));
result2(:,:) = mean(log2(time),1);
figure;
plot(range_purity,result2(1,:),'--*',range_purity,result2(2,:),'--O',range_purity,result2(3,:),'--+',range_purity,result2(4,:),'--<',range_purity,result2(5,:),'--d',range_purity,result2(6,:),'-->',range_purity,result2(7,:),'--s',range_purity,result2(8,:),'--x',range_purity,result2(9,:),'--^','markersize',10);
xlabel('purity');
ylabel('log time(s)');
% title(strcat('r=',num2str(r),',m=',num2str(m)));
legend('MV-Dual','GFPI','min vol \lambda = 0.1','min vol \lambda = 1','min vol \lambda = 5','SNPA','MVIE','HyperCSI','MVES');
xlim([startp,endp]);
xticks(range_purity);
xticklabels(range_purity);
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))