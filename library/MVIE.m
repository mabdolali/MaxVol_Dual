function W = MVIE(X,r)
%% Dimension Reduction
b=mean(X,2);
X = X - b;
% MM = M;
[Xp, Phi] = lindimred(X,r-1,1);
%% Facet Enumaration
[G,h] = vert2con(Xp');
G = G';
K = size(G,2);
%% Convex Optimization
n = r-1;
cvx_begin
variable F(n,n) symmetric;
variable c(n);
minimize(-log_det(F));
subject to
for k =1 : K
    norm(F*G(:,k))+G(:,k)'*c <= h(k);
end
F == semidefinite(n);
cvx_end
Qp=[];
for k = 1: K
    if norm(norm(F*G(:,k))+G(:,k)'*c - h(k))<1e-3
        Qp = [Qp F*(F*G(:,k)/norm(F*G(:,k)))+c];
    end
end
Q = Phi*Qp + b;

if size(Q,2)> r
    [idx,~]=kmeans(Q',r);
    Qs=zeros(r,r);
    for k=1:r
        Qs(:,k)=mean(Q(:,find(idx==k)),2);
    end
    Q = Qs;
end
abar = sum(Q,2);
W=[];
for k = 1: size(Q,2)
    W = [W abar - (r-1)*Q(:,k)];
end

% figure(1);
% plot3(W(1,:),W(2,:),W(3,:),'rO');
end
