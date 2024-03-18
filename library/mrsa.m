function [mean_score score]= mrsa( A,B )
% Compute MRSA between matrix A and B (B is colum-matched with A)
na = size(A,2);
nb = size(B,2);
if na ~= nb
   error('Input matrices have unmatched number of columns'); 
end

for i = 1 : na
 a = A(:,i);
 b = B(:,i);
 a = a-mean(a);
 b = b-mean(b);
 a = a/norm(a,2);
 b = b/norm(b,2);
 score(i) = abs(acos(a'*b))*100/pi;
end

 mean_score = mean(score);
end