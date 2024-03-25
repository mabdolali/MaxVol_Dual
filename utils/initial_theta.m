function [Z_tilde,flag] = initial_theta(Y,r,lambda0,Z_tilde,options)

cputime0 = cputime; 
if nargin <= 4
    options = [];
end
if ~isfield(options,'maxiter')
    options.maxiter = 100; 
end
if ~isfield(options,'timelimit')
    options.timelimit = 60; 
end
n = size(Y,2);
p = r;
m = r-1;
Theta =zeros(r-1,r);
Z_0 = Z_tilde;
iter = 1;
Z_pre = randn(p,r);
Delta = [];
vals = [];
H = sparse(p+m+n+m,p+m+n+m);
H(p+m+1:end-m,p+m+1:end-m)=speye(n);
lambda = lambda0;

while iter < options.maxiter ... 
        && norm(Z_tilde - Z_pre,'fro')/norm(Z_pre,'fro')>1e-2 ...
        && cputime-cputime0 <= options.timelimit 
    Z_pre = Z_tilde;
    iter = iter + 1;
    for i = 1 : r
        inv_Z = pinv(Z_tilde);
        zz = Z_tilde;
        zz(:,i)=[];
        f = inv_Z(i,:);
        f = [f zeros(1,m) zeros(1,n) zeros(1,m)];
        A = [zeros(n,p) Y' -speye(n) zeros(n,m)];
        b = ones(n,1);
        Aeq1 = [eye(p) -eye(p,r-1) zeros(p,n) zeros(p,m)];
        beq1 = [zeros(m,1);1];
        Aeq2 = [zeros(m,p) eye(m) zeros(m,n), zz(1:m,:)];
        beq2 = zeros(m,1);
        opts = optimoptions('quadprog','Display','none'); 
        [Z_tilde_i,~,e,~] = quadprog(lambda*H,-f',A,b,[Aeq1;Aeq2],[beq1;beq2],[-inf*ones(p+m+n,1); 0.01*ones(m,1)],[],[],opts);
        if e<0 && e~=-4
            Z_tilde = zeros(p,r);
            i=1;
            flag = 0;
            return
        else
            Z_tilde(:,i) = Z_tilde_i(1:p);
            theta = Z_tilde_i(p+1:p+m);
            Delta(:,i) = Z_tilde_i(p+m+1:end-m);
            Theta(:,i) = theta;
            alpha(:,i) = Z_tilde_i(end-m+1:end);
        end
        
    end
end
flag = 1;
Z = Z_tilde;
end
