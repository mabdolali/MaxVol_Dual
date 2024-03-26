% This code performs an initial estimate of Z without any updating for
% lambda. The code is very similar to algorithm2_update_theta.m. 
function [Z_tilde,flag] = initial_theta(Y,r,lambda0,Z_tilde,options)

%cputime0 = cputime; 
if nargin <= 4
    options = [];
end
if ~isfield(options,'maxiter')
    options.maxiter = 100; 
end
% if ~isfield(options,'timelimit')
%     options.timelimit = 60; 
% end
%initialization
n = size(Y,2);
p = r;
m = r-1;
Theta =zeros(r-1,r);
iter = 1;
Z_pre = randn(p,r);
Delta = [];
% initializing H matrix for the QP optimization problem
H = sparse(p+m+n+m,p+m+n+m);
H(p+m+1:end-m,p+m+1:end-m)=speye(n);
lambda = lambda0;
%main loop
while iter < options.maxiter ... 
        && norm(Z_tilde - Z_pre,'fro')/norm(Z_pre,'fro')>1e-2 %...
        %&& cputime-cputime0 <= options.timelimit 
    Z_pre = Z_tilde;
    iter = iter + 1;
    for i = 1 : r %updating each column of Z sequentially
        inv_Z = pinv(Z_tilde); %compute z^-1
        zz = Z_tilde;
        zz(:,i)=[]; % taking the i-th (current) column out
        f = inv_Z(i,:); %Z^-1(i,:), refer to the paper for the details
        f = [f zeros(1,m) zeros(1,n) zeros(1,m)]; %optimizaton parameters: [Z, Theta, Delta, alpha]
        A = [zeros(n,p) Y' -speye(n) zeros(n,m)]; %inequality constraint: Y' * Theta <= 1+ Delta --> Y' * Theta - Delta <= 1
        b = ones(n,1); 
        Aeq1 = [eye(p) -eye(p,r-1) zeros(p,n) zeros(p,m)]; %equality constraint 1: Z = [Theta; e']
        beq1 = [zeros(m,1);1];
        Aeq2 = [zeros(m,p) eye(m) zeros(m,n), zz(1:m,:)]; %equality constraint 2: 0 to be in the interior of estimated Thetas
        beq2 = zeros(m,1);
        opts = optimoptions('quadprog','Display','none'); %Turn off display information of the optimizer
        [Z_tilde_i,~,e,~] = quadprog(lambda*H,-f',A,b,[Aeq1;Aeq2],[beq1;beq2],[-inf*ones(p+m+n,1); 0.01*ones(m,1)],[],[],opts);
        if e<0 && e~=-4 % check if optimization has failed
            Z_tilde = zeros(p,r);
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
