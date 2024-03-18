function [v1, West, best_theta, iter, Y, C] = maxvoldual(X,r,lambda,num_workers)
MAX =  max(max(X));
X = X / MAX ;
MEAN = mean(X,2);
v = MEAN;
%Y = X - v;
%U = Y*Y';
%[C,~] = eigs(U,r-1);
%Y = C'*Y;
iter = 0;
v1 = ones(size(MEAN));
v1 = v1 / norm(v1,'fro');
nn = 0;
theta = cell(num_workers,1);
ignore = cell(num_workers,1);
delta = cell(num_workers,1);
flag = cell(num_workers,1);
z = cell(num_workers,1);
for i = 1: num_workers
    z{i}=[randn(r-1,r);ones(1,r)];
end
while norm(v-v1,'fro')/norm(v1,'fro') > 1e-2 && iter < 20 %norm(W1 - W2, 'fro')/norm(W1)
    v1 = v;
    Y = X - v;
    U = Y*Y';
    [C,~] = eigs(U,r-1);
    Y = C'*Y;
    iter = iter + 1;
    % find theta


    parfor i=1:num_workers
        [z{i},theta{i},delta{i},flag{i}] = algorithm2_update_theta(Y,r,lambda,z{i});
        ignore{i} = ~flag{i} || any(is_correct(normc(theta{i}))<0.01);
    end
    best_theta = [];
    best = 0;
        for i = 1: num_workers
            if ~ignore{i}
                vol = (det(z{i}))^2-lambda*norm(delta{i},'fro');
                if vol > best
                    best_theta = theta{i};
                    best = vol;
                end
            else
                z{i}=[randn(r-1,r);ones(1,r)];
            end
        end
        if isempty(best_theta)
            nn = nn + 1;
            [J,~] = SNPA(X,nn*r);
            v = mean(X(:,J),2);
            Y = X - v;
            U = Y*Y';
            [C,~] = eigs(U,r-1);
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
            W2 = W_e;
            West = C * W_e + v;
            v = mean(West,2);
        end
end
West = West * MAX;

