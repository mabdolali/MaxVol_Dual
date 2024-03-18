% Find vertices V of the set P = { x in R^{k-1} | Cx+d >= 0 } with brute force
function V = vertices(W)
[m,r] = size(W);
cro=nchoosek(1:r,r-1);
V=[];
for i=1:size(cro,1)
    c=cro(i,:);
    A = W(:,c);
    b = ones(m,1);
    theta = A' \ b;
    V = [V theta];
    
end