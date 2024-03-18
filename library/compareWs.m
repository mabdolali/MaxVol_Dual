% Compute
%
% errW = min_perm ||W - West(:,perm)||_F / ||W||_F

function errW = compareWs( W, West );
West(isnan(West))=0;
r = size(W,2); 

% Compute the MRSA: V vs H 
for i = 1 : r
    for j = 1 : r
        Dist(i,j) = norm( W(:,i) - West(:,j) , 'fro' )^2; 
    end
end
 
[assignment,cost] = munkres(Dist); 

errW = norm( W - West(:,assignment) , 'fro' ) / norm(W,'fro'); 