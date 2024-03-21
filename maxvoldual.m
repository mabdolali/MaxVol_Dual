%   This code solves the maxvol-dual problem
%
%        max_{Z,Theta,Delta} det(Z)^2 - lambda ||Delta||^2 
%                      such that Z = [Theta; e'] and Y' Theta <= 1 + Delta 
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
% ****** Output ******
% 
% v1          :    estimated center vector
% West        :    estimated W
% best_theta  :    final Theta at the convergence with maximum dual volume
% iter        :    number of iterations till convergence
% Y           :    projected points in r-1 dimensions
% C           :    projection matrix

function [v1, West, best_theta, iter, Y, C] = maxvoldual(X,r,lambda,options)
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
% pre-processing
MAX =  max(max(X));
X = X / MAX ;
MEAN = mean(X,2);
v = MEAN;
Y = X - v;
[C,~,~]=svds(Y,r-1);
% initialization
iter = 0;
v1 = ones(size(MEAN));
v1 = v1 / norm(v1,'fro');
nn = 0;
theta = cell(options.num_workers,1);
ignore = cell(options.num_workers,1);
delta = cell(options.num_workers,1);
flag = cell(options.num_workers,1);
z = cell(options.num_workers,1);
for i = 1: options.num_workers
    z{i}=[randn(r-1,r);ones(1,r)];
end
% main loop
while norm(v-v1,'fro')/norm(v1,'fro') > options.epsilon && iter < options.maxiter 
    % projection
    v1 = v;
    Y = X - v;
    %U = Y*Y';
    %[C,~] = eigs(U,r-1);
    %[C,~,~]=svds(Y,r-1);
    Y = C'*Y;
    iter = iter + 1;
    % parallel computation of Z_i for i = 1: num_workers
    parfor i=1:options.num_workers
        [z{i},theta{i},delta{i},flag{i}] = algorithm2_update_theta(Y,r,lambda,z{i});
        ignore{i} = ~flag{i} || any(is_correct(normc(theta{i}))<0.01);
    end
    % select the best candidate till now (maximum dual volume)
    best_theta = [];
    best = 0;
        for i = 1: options.num_workers
            if ~ignore{i}
                vol = (det(z{i}))^2;
                if vol > best
                    best_theta = theta{i};
                    best = vol;
                end
            else
                z{i}=[randn(r-1,r);ones(1,r)];
            end
        end
     % if all candidates failed, use alternative initialization
        if isempty(best_theta)
            nn = nn + 1;
            [J,~] = SNPA(X,nn*r);
            v = mean(X(:,J),2);
            Y = X - v;
%             U = Y*Y';
%             [C,~] = eigs(U,r-1);
            [C,~,~]=svds(Y,r-1);
            Y = C'*Y;
            v1=randn(r,1);
        else
      % find intersections (W)
            W_e = [];
            cro=nchoosek(1:r,r-1); % for each r-1 facets from r facets
            for i=1:size(cro,1)
                c=cro(i,:);
                U = best_theta(:,c);
                O = ones(r-1,1);
                coef=U' \ O;
                W_e=[W_e coef];
            end
       % update mean vector
            W2 = W_e;
            West = C * W_e + v;
            v = mean(West,2);
        end
end
West = West * MAX;

