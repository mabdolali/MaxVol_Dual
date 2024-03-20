% checks whether 0 is in the interior of the convex hull of columns of
% Theta_c
function [c] = is_correct(Theta_c)
[m,num_c] = size(Theta_c);
f = [zeros(num_c,1)];
Aeq1 = [Theta_c];
beq1 = zeros(m,1);
Aeq2 = [ones(1,num_c)];
beq2 =1;
Aeq = [Aeq1;Aeq2];
beq = [beq1;beq2];
lb = zeros(num_c,1);
ub = [ones(num_c,1)];
[sol,fval] = linprog(-f,[],[],Aeq,beq,lb,ub);
try
    c = sol(1:num_c);
catch
    c = nan;
end
end

