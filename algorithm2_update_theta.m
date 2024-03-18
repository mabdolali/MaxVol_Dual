function [Z_tilde,Theta,Delta,flag] = algorithm2_update_theta(Y,r,lambda0,Z_tilde)

n = size(Y,2);
p = r;
m = r-1;
Theta =zeros(r-1,r);
MAX_ITER = 100;
[Z_tilde,flag] = initial_theta(Y,r,lambda0,Z_tilde);
if ~flag
    Theta = randn(r-1,r);
    Delta = randn(r,n);
    return 
end
%Z_tilde = normc(Z_tilde)/abs(det(Z_tilde)^(1/r));
Z_0 = Z_tilde;
iter = 1;
Z_pre = randn(p,r);
Delta = [];
vals = [];
H = sparse(p+m+n+m,p+m+n+m);
H(p+m+1:end-m,p+m+1:end-m)=speye(n);
lambda = lambda0;

while(iter < MAX_ITER && norm(Z_tilde - Z_pre)/norm(Z_pre)>1e-3)
    Z_pre = Z_tilde;
    iter = iter + 1;
    for i = 1 : r
        lambda = lambda0 * det(Z_tilde)^2/det(Z_0)^2;
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
        [Z_tilde_i,~,e,~] = quadprog(lambda*H,-f',A,b,[Aeq1;Aeq2],[beq1;beq2],[-inf*ones(p+m+n,1); 0.01*ones(m,1)],[]);
        if e<0 && e~=-4
            Z_tilde = randn(p,r);
            i=1;
            flag = 0;
            return
        else
            Z_0 = Z_tilde;
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