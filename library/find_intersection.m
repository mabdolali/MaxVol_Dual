function W = find_intersection(Indices,M,r)
d = size(M,1);
U=[];O=[];
W = [];
if d == r-1 % if full-rank
    cro=nchoosek(1:r,r-1); % for each r-1 facets from r facets
    for i=1:size(cro,1)
        c=cro(i,:);
        U=[];
        O=[];
        for j = 1:length(c)
            data = M(:,Indices{c(j)}); % find data points on each facet
            mean_data = mean(data,2);  % subtract mean
            [u,~,~]=svd(data-mean_data); % calculate left singular vector as normal vector
            U=[U u(:,end)];
            O =[O ; mean(u(:,end)'*data)]; % calculate offset
        end
        coef=U' \ O; % solve a linear system of equations
        W=[W coef]; % add the founded vertex in primal space
    end
    
else % if rank deficient
    
    for j = 1:r % find the H representation
        data = M(:,Indices{j});
        mean_data = mean(data,2);
        [u,~,~]=svd(data-mean_data);
        if mean(u(:,end)'*data) > 0
            u(:,end) = - u(:,end);
        end
        U=[U u(:,end)];
        O =[O ; mean(u(:,end)'*data)];
        
    end
    W=con2vert(U',O); % convert the H representation to V representation
    W = W';
end
end

