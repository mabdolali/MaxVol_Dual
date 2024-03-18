% Linear dimensionality reduction 
% 
% *** Description ***
% Given an m-by-n matrix and a dimension r, perform a linear dimensionality
% reduction to dimension r. 
% 
% [Mr, u] = lindimred(M,r,dimred) 
%
% ****** Input ******
% M          : m-by-n matrix
% r          : dimension of the space to which the data is projected. 
% dimred     : parameter allowing to choose the dimensionality reduction
%              technique used (for m > r): 
%              dimred=1: use the truncated SVD, that is, svds(M,r). 
%              dimred=2: use the truncated SVD, but computation based on MM'.
%              dimred~=1,2: use a random projection (u = randn(m,r))
%              (default = 1.)
%
% ****** Output ******
% Mr          : r-by-n reduced data matrix M. 
% u           : corresponding projection, that is, Mr = u'*M. 

function [Mr, u] = lindimred(M,r,dimred) 

[m,n] = size(M); 
if r > min(m,n)
    error('r has to be smaller than m and n'); 
end
if nargin <= 2, 
    dimred = 1; 
    % If m << n, use SVD of MM^T. It is computationally much cheaper although
    % it is numerically much less stable (condition number squared). 
    if m^2 <= n, dimred = 2; end
    % If m,n,r are too large, use a random projection
    if m*n*r*min(m,n) >= 1e12, dimred = 3; end
end

if m == r % dimensionality reduction not applicable
    u = eye(r);
else   
    if dimred == 1 % Using truncated SVD
        [u,s,v] = svds(M,r); 
    elseif dimred == 2 % Using truncated SVD of MM'
        MMt = M*M'; 
        [u,s,v] = svds(MMt,r); 
    else  % Using random projection
        u = randn(m,r);  
    end
end
Mr = u'*M;  