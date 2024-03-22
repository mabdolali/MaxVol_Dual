% This code compares the performance of MV-Dual with GFPI, Min-Vol, HyperCSI, SNPA and MVIE in the noisy case for the full rank
% matrices. 
clc
clear all
close all
addpath(genpath('..'));
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultLineLineWidth', 2);


%% generate data & initialization
m = 4; % dimension
r = m; % # of endmembers
SNR = 60;
num_experiments=10; % # of trials
startp = (1/(r-1)+0.01); % starting purity value
endp = 1; % ending purity value
step = (endp-startp)/6;
range_purity = startp:step:endp; % define range of purity values
time = zeros(num_experiments,length(range_purity));
result = zeros(num_experiments,9, length(range_purity));
ind = 0;
mvdual_options.num_workers = 10;

for purity = range_purity
    ind = ind +1;
    for no= 1 : num_experiments
        Ni1 = 30*ones(r,1); % # of points on each facets
        Ni2 = 10; % # of points within polytope
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
        end 
        [v, West, theta, iter] = maxvoldual(M,r,lambda,mvdual_options);
        result(no,1,ind)=compareWs(Wg, West);
        time(no,1,ind) = toc;
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
        West = min(W1,1);
        result(no,2,ind)=compareWs(Wg, West);
        time(no,2,ind) = toc;
        %%
        tic;
        options.lambda=0.1;
        [W2,~,~,~,~,~] = minvolNMF(M,r,options);
        result(no,3,ind)=compareWs( Wg, W2 );%norm(W2 - Wg,'fro')/norm(Wg,'fro');
        time(no,3,ind) = toc;
        %%
        tic;
        options.lambda=1;
        [W2,~,~,~,~,~] = minvolNMF(M,r,options);
        result(no,4,ind)=compareWs( Wg, W2 );
        time(no,4,ind) = toc;

        %%
        tic;
        options.lambda=5;
        [W2,~,~,~,~,~] = minvolNMF(M,r,options);
        result(no,5,ind)=compareWs( Wg, W2 );
        time(no,5,ind) = toc;

        %% SNPA
        tic;
        K = SNPA(M,r);
        W3 = M(:,K);
        result(no,6,ind)=compareWs( Wg, W3 );
        time(no,6,ind) = toc;

        %% MVIE
        tic;
        W4 = MVIE(M,r);
        result(no,7,ind)=compareWs( Wg, W4 );
        time(no,7,ind) = toc;

        %% HyperCSI
        tic;
        [W5, ~, ~] = HyperCSI(M,r);
        result(no,8,ind)=compareWs( Wg, W5 );
        time(no,8,ind) = toc;
        %% MVES
        tic;
        [W6,S_est,iter_cnt] = MVES(M,r,1);
        result(no,9,ind)=compareWs( Wg, W6);
        time(no,9,ind) = toc;
    end
end
% show results
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
result2(:,:) = mean(log(time),1);
figure;
plot(range_purity,result2(1,:),'--*',range_purity,result2(2,:),'--O',range_purity,result2(3,:),'--+',range_purity,result2(4,:),'--<',range_purity,result2(5,:),'--d',range_purity,result2(6,:),'-->',range_purity,result2(7,:),'--s',range_purity,result2(8,:),'--x',range_purity,result2(9,:),'--^','markersize',10);
xlabel('purity');
ylabel('log time(s)');
legend('MV-Dual','GFPI','min vol \lambda = 0.1','min vol \lambda = 1','min vol \lambda = 5','SNPA','MVIE','HyperCSI','MVES');
xlim([startp,endp]);
xticks(range_purity);
xticklabels(range_purity);
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))