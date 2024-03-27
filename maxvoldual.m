%   function [v1, West, best_theta, iter, Y, C] = maxvoldual(X,r,lambda,options)
% 
%    This code solves simplex-structured matrix factorization (SSMF) 
%    via a dual approach. Given the input matrix, X, and a factorization
%    rank r, it looks for W and H such that WH approximates X and H is column 
%    stochastic. 
%    To do so, it first reduces the dimension of the problem: 
%        Y = C' (X - v e') where e is the vector of all ones, v is in conv(X)
%                                C' contains the first r singular vectors of 
%                                X - v e'. 
% 
%    Then it solves the maximum-volume dual problems: 
%
%        max_{Z,Theta,Delta}  det(Z)^2 - lambda ||Delta||^2
%                  such that  Z = [Theta; e'] and Y' Theta <= 1 + Delta. 
%
%    where Theta represents the polar of conv(W) whose volumns is maximized. 
% 
% See the paper "Dual Simplex Volume Maximization for Simplex-Structured 
% Matrix Factorization", by M. Abdolali, G. Barbarino and N. Gillis, 2024.
% 
% ****** Input ******
% X      :  the input matrix
% r      :  the rank of the sought approximation
% lambda :  the regularization parameter
% ---Options---
% .maxiter    : the maximum number of iterations performed by the algorithm
%             -default = 20.
% .epsilon    : the tolerance level for convergence
%             -default = 1e-2.
% .num_workers: number of parallelized solutions used
%             -default = 5.
% .timelimit  : maximum alloted time for the outer loops
%             -default = 300.
% 
% ****** Output ******
% v1          :    estimated center vector
% West        :    estimated W
% best_theta  :    final Theta at the convergence with maximum dual volume
% iter        :    number of iterations till convergence
% Y           :    projected points in r-1 dimensions
% C           :    projection matrix

function [v1, West, best_theta, iter, Y, C] = maxvoldual(X,r,lambda,options)
tic; 
if nargin <= 3
    options = [];
end
if ~isfield(options,'maxiter')
    options.maxiter = 20;
end
if ~isfield(options,'epsilon')
    options.epsilon = 1e-2;
end
if ~isfield(options,'num_workers')
    options.num_workers = 5;
end
if ~isfield(options,'timelimit')
    options.timelimit = 300;
end
% pre-processing
MAX =  max(max(X));
X = X / MAX ; %scaling for numerical stability
v = mean(X,2);
if issparse(X)
    [m,n]=size(X); %if X is sparse, compute svds implicitly
    [C,~,~]=svds(@afun,[m,n],r-1);
else
    Y = X - v; %otherwise compute svds normally
    [C,~,~]=svds(Y,r-1);
end
% initialization
iter = 0; %initializing iteration counter
v1 = 0; % initializing v_{k-1} (keeps center vector of previous iteration)
nn = 0; % a counter for # of endmembers in SNPA initalization if average initialization goes wrong
theta = cell(options.num_workers,1); % Theta_i for i=1,...num_workers 
ignore = cell(options.num_workers,1); % ignore_i for i=1,...num_workers (ignoring i-th solution or not)
delta = cell(options.num_workers,1); % delta_i for i=1,...num_workers (estimated noise matrix)
flag = cell(options.num_workers,1); % flag_i for i=1,...num_workers (optimization was successful or not)
z = cell(options.num_workers,1); % z_i for i=1,...num_workers
for i = 1: options.num_workers
    z{i}=[randn(r-1,r);ones(1,r)]; %initialize z_i
end
CtX = C'*X;
O = ones(r-1,1);
cro = nchoosek(1:r,r-1); % for calculating intersection later
outeriter = 1;
West = []; 
% main loop
while norm(v-v1,'fro')/norm(v1,'fro') > options.epsilon ...
        && iter < options.maxiter ...
        && toc <= options.timelimit  
    fprintf('Outer iteration %1.0d started. \n',outeriter); 
    % projection
    v1 = v;
    Y = CtX - C'*v;
    iter = iter + 1;
    % parallel computation of Z_i for i = 1: num_workers
    parfor i=1:options.num_workers
        [z{i},theta{i},delta{i},flag{i}] = algorithm2_update_theta(Y,r,lambda,z{i}); %solves a QP (refer to the paper for more info)
        ignore{i} = ~flag{i} || any(is_correct(normc(theta{i}))<0.01); % check if optimization was successful
    end
    % select the best candidate till now (maximum dual volume)
    best_theta = []; 
    best = -Inf; 
    for i = 1: options.num_workers
        if ~ignore{i} %if optimization for i-th solution was successful
            vol = (det(z{i}))^2 - lambda*sum(delta{i}(:).^2); % evaluate the objective function for i-th candidate
            if vol > best
                best_theta = theta{i}; %update the best candidate
                best = vol;
            end
        else
            z{i}=[randn(r-1,r);ones(1,r)]; % if optimization was unseccessful, change initialization for next iteration
        end
    end
    % if all candidates failed, use alternative initialization
    if isempty(best_theta)
        nn = nn + 1;
        [J,~] = SNPA(X,nn*r); %estimate endmembers with SNPA
        v = mean(X(:,J),2);
        if issparse(X)
            [m,n]=size(X); %if X is sparse, compute svds implicitly
            [C,~,~]=svds(@afun,[m,n],r-1);
        else
            Y = X - v;
            [C,~,~]=svds(Y,r-1); %if not sparse, compute svds normally
        end
        CtX = C'*X; %update CtX for next iterations
        v1 = 0;
    else
        % find intersections (W)
        W_e = [];
        for i=1:size(cro,1) %for each r-1 facets compute intersection
            c = cro(i,:);
            U = best_theta(:,c);
            coef = U' \ O;
            W_e = [W_e coef];
        end
        % update mean vector
        West = C * W_e + v;
        v = mean(West,2);
    end
    outeriter = outeriter + 1;
end
West = West * MAX; %re-scale to the original space

function Ax = afun(x, cond)
    % function for implicitly computing svds of X-v*e^T where X is a sparse
    % matrix. The function afun satisfies these required conditions:
    % Afun(x,'notransp') accepts a vector x and returns the product A*x.
    % Afun(x,'transp') accepts a vector x and returns the product A'*x.

    if strcmp(cond,'notransp')
        Ax = X * x - v * sum(x);
    else
        Ax = X' * x - v' * x;
    end
end
end
