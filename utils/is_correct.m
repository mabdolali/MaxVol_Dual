% checks whether 0 is in the interior of the convex hull of columns of
% Theta_c
function [c] = is_correct(Theta_c)
[m,num_c] = size(Theta_c);
f = [zeros(num_c,1)]; %no objective function, this is only a feasibility checking, parameter: c \in R^num_c
Aeq1 = [Theta_c]; % equality constraint 1: Theta_c * c = 0
beq1 = zeros(m,1);
Aeq2 = [ones(1,num_c)]; % equality constraint 2: c'*e = 1
beq2 =1;
Aeq = [Aeq1;Aeq2];
beq = [beq1;beq2];
lb = zeros(num_c,1); % c>=0 : convexity constraint
ub = [ones(num_c,1)]; % c<=1 : convexity constraint
[sol,fval] = linprog(-f,[],[],Aeq,beq,lb,ub); %LP with constraints
try
    c = sol(1:num_c); %if successful return coefficients
catch
    c = nan; %if not successful return nan
end
end

