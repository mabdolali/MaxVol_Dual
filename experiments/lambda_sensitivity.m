% Sensetivity analysis of MV-Dual vs parameter lambda
% This code analyzes the performance of MV-Dual with respect to different
% lambda values
clc
clear all
close all
addpath(genpath('library'));
%% Setting
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultLineLineWidth', 2);
%% Parameters
SNR = 30;
%% generate data

m = 3; % dimension
r = m; % # of endmembers

num_experiments=10; % # of trials
startp = (1/(r-1)+0.01); % starting purity value
endp = 1; % ending purity value
step = (endp-startp)/6;
range_purity = startp:step:endp; % define range of purity values
time = zeros(num_experiments,length(range_purity));
result = zeros(num_experiments,10, length(range_purity));
ind = 0;
lambda_range = [0.01, 0.05, 0.1, 0.4, 0.7, 1, 1.5, 2, 5, 10];
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
        [v, West, theta, iter] = maxvoldual(M,r,lambda_range(1));
        result(no,1,ind)=compareWs(Wg, West);
        time(no,1,ind) = toc;
        %%
        tic;              
        [v, West, theta, iter] = maxvoldual(M,r,lambda_range(2));
        result(no,2,ind)=compareWs(Wg, West);
        time(no,2,ind) = toc;
        %%
        tic;      
        [v, West, theta, iter] = maxvoldual(M,r,lambda_range(3));
        result(no,3,ind)=compareWs(Wg, West);
        time(no,3,ind) = toc;
        %%
        tic;      
        [v, West, theta, iter] = maxvoldual(M,r,lambda_range(4));
        result(no,4,ind)=compareWs(Wg, West);
        time(no,4,ind) = toc;

        %%
        tic;      
        [v, West, theta, iter] = maxvoldual(M,r,lambda_range(5));
        result(no,5,ind)=compareWs(Wg, West);
        time(no,5,ind) = toc;

        %% 
        tic;      
        [v, West, theta, iter] = maxvoldual(M,r,lambda_range(6));
        result(no,6,ind)=compareWs(Wg, West);
        time(no,6,ind) = toc;

        %% 
        tic;      
        [v, West, theta, iter] = maxvoldual(M,r,lambda_range(7));
        result(no,7,ind)=compareWs(Wg, West);
        time(no,7,ind) = toc;
        %% 
        tic;      
        [v, West, theta, iter] = maxvoldual(M,r,lambda_range(8));
        result(no,8,ind)=compareWs(Wg, West);
        time(no,8,ind) = toc;
        %% 
        tic;      
        [v, West, theta, iter] = maxvoldual(M,r,lambda_range(9));
        result(no,9,ind)=compareWs(Wg, West);
        time(no,9,ind) = toc;
        %% 
        tic;      
        [v, West, theta, iter] = maxvoldual(M,r,lambda_range(10));
        result(no,10,ind)=compareWs(Wg, West);
        time(no,10,ind) = toc;
    end
end
result2=zeros(size(result,2),length(range_purity));
result2(:,:) = mean(result,1);
figure;
plot(range_purity,result2(1,:),'--*',range_purity,result2(2,:),'--O',range_purity,result2(3,:),'--+',range_purity,result2(4,:),'--<',range_purity,result2(5,:),'--d',range_purity,result2(6,:),'-->',range_purity,result2(7,:),'--s',range_purity,result2(8,:),'--x',range_purity,result2(9,:),'--^',range_purity,result2(10,:),'--v','markersize',10);
xlabel('purity');
ylabel('|| W - W_t || / || W_t ||');
% title(strcat('r=',num2str(r),',m=',num2str(m)));
a = mat2cell(string(lambda_range),1);
legend(cellstr(a{1}))
xlim([startp,endp]);
xticks(range_purity);
xticklabels(range_purity);
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))
%%
result2=zeros(size(time,2),length(range_purity));
result2(:,:) = mean(log(time),1);
figure;
plot(range_purity,result2(1,:),'--*',range_purity,result2(2,:),'--O',range_purity,result2(3,:),'--+',range_purity,result2(4,:),'--<',range_purity,result2(5,:),'--d',range_purity,result2(6,:),'-->',range_purity,result2(7,:),'--s',range_purity,result2(8,:),'--x',range_purity,result2(9,:),'--^',range_purity,result2(10,:),'--v','markersize',10);
xlabel('purity');
ylabel('log time(s)');
a = mat2cell(string(lambda_range),1);
legend(cellstr(a{1}));
xlim([startp,endp]);
xticks(range_purity);
xticklabels(range_purity);
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'))