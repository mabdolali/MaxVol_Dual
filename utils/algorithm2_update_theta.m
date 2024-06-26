% This code is for the optimization of the following problem for a fixed
% center vector v, after projecting samples X on the r-1 dimensional space
% and obtaining matrix Y:
%        max_{Z,Theta,Delta} det(Z)^2 - lambda ||Delta||^2 
%                      such that Z = [Theta; e'] and Y' Theta <= 1 + Delta 
% ****** Input ******
% Y      :  the projected input matrix
% r      :  the rank of the sought approximation
% lambda0:  the regularization parameter
% Z_tilde: initial Z matrix
% options.maxiter 
% 
% ****** Output ******
% 
% Z_tilde     :    estimated Z matrix
% Theta       :    estimated Theta matrix in the dual space (vertices of
% the simplex in the dual space)
% Delta       :    estimated noise matrix
% flag        : a boolean flag which indicates whether the optimization was
% successful or not

function [Z_tilde,Theta,Delta,flag] = algorithm2_update_theta(Y,r,lambda0,Z_tilde,options)

if nargin <= 4
    options = [];
end
if ~isfield(options,'maxiter')
    options.maxiter = 5; 
end
if ~isfield(options,'timelimit')
    options.timelimit = 60; 
end
% initialization
n = size(Y,2);
p = r;
m = r-1;
Theta =zeros(r-1,r);
% initialize Z (and theta)
[Z_tilde,flag] = initial_theta(Y,r,lambda0,Z_tilde,options);
% checking whether initial optimization problem was successful
if ~flag
    Theta = randn(r-1,r);
    Delta = randn(r,n);
    return 
end
Z_0 = Z_tilde;
iter = 1;
Z_pre = randn(p,r);
Delta = [];
H = sparse(p+m+n+m,p+m+n+m);
H(p+m+1:end-m,p+m+1:end-m)=speye(n);
while iter < options.maxiter ... 
        && norm(Z_tilde - Z_pre,'fro')/norm(Z_pre,'fro')>1e-3
    Z_pre = Z_tilde;
    iter = iter + 1;
    for i = 1 : r %updating each column of Z sequentially
        lambda = lambda0 * det(Z_tilde)^2/det(Z_0)^2; %update lambda
        inv_Z = pinv(Z_tilde); % compute Z^-1
        zz = Z_tilde;
        zz(:,i)=[]; % taking the i-th (current) column out (for second equality constraint on alpha)
        f = inv_Z(i,:); %Z^-1(i,:), refer to the paper for the details
        f = [f zeros(1,m) zeros(1,n) zeros(1,m)]; %optimizaton parameters: [Z, Theta, Delta, alpha]
        A = [zeros(n,p) Y' -speye(n) zeros(n,m)]; %inequality constraint: Y' * Theta <= 1+ Delta --> Y' * Theta - Delta <= 1
        b = ones(n,1);
        Aeq1 = [eye(p) -eye(p,r-1) zeros(p,n) zeros(p,m)]; %equality constraint 1: Z = [Theta; e']
        beq1 = [zeros(m,1);1];
        Aeq2 = [zeros(m,p) eye(m) zeros(m,n), zz(1:m,:)]; %equality constraint 2: 0 to be in the interior of estimated Thetas
        beq2 = zeros(m,1);
        opts = optimoptions('quadprog','Display','none'); 
        [Z_tilde_i,~,e,~] = quadprog(lambda*H,-f',A,b,[Aeq1;Aeq2],[beq1;beq2],[-inf*ones(p+m+n,1); 0.01*ones(m,1)],[],[],opts);
        if e<0 && e~=-4 %check whether optimization was successful
            Z_tilde = randn(p,r);
            flag = 0;
            return
        else %update i-th column of Z, Theta, Delta and i-th element of alpha
            Z_tilde(:,i) = Z_tilde_i(1:p);
            theta = Z_tilde_i(p+1:p+m);
            Delta(:,i) = Z_tilde_i(p+m+1:end-m);
            Theta(:,i) = theta;
            alpha(:,i) = Z_tilde_i(end-m+1:end);
        end
    end
end
flag = 1;
end
