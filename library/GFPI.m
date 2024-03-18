function [W,H,thetas,num,facets_mean,outflag] ... 
    = GFPI(X,r,options)
% This function implements the greedy facet-based identification (GFPI) 
% algorithm for simplex-structured matrix factorization (SSMF). 
% 
% Assumption: The input matrix, X = WH, satisfies the facet-based
% condition (FBC), that is, 
% a) The r columns of W are vertices of conv(W). 
% b) The columns of H belong to the unit simplex. 
% c) H is sufficiently sparse so that each facet of conv(W) contains 
%      s >= d=rank(W) columns of X. 
% d) Every facet of conv(X) which is not a facet of conv(W) contains
%      strictly less than s columns of X. 
%
% Under this assumption, GFPI recovers the columns of W by sequentially
% identifying the facets of conv(W) (approximately in the presence of
% noise). The matrix H is recovered solving a convex optimization problem. 
%
% See the paper 'Simplex-Structured Matrix Factorization: Sparsity-based  
% Identifiability and Provably Correct Algorithms', by Maryam Abdolali and
% Nicolas Gillis, arXiv, July 2020. 
% 
% ****** Input ****** 
% X : m-by-n matrix to factorize 
% r : # of vertices
% options:
%   - d : intrinsic dimension of X
%   - T : number of facets to be extracted
%   - gamma : saftey gap
%   - lambda : trade off between noise and # points on facets
%   - eta : margin parameter which controls how far the facets are expected
%   to be from each other
%   - no_show = (default: true) not displaying the intermediate figures for the
%   selected facets in each step
%   - timelimit= (default: 1000) maximum allowed time (in seconds) for
%   finding each facet. This is a parameter to CPLEX.
%   - centerstrategy= (default: mean) the strategy used for finding the
%   center of a facet, can be average of points ('mean') or average of
%   points selected by SNPA ('snpa')
%   - outlier= (default: false) is data occluded with outliers or not
% 
% ****** Output ****** 
% (W,H) : low-rank approximation W*H of X, and columns of H
% belong to the unit simplex.
% thetas : computed thetas
% num : total number of points on all computed facets
% facets_mean : the center points on each facet
% outflag : indicating whether the timelimit is exceeded or not

%setting paramters and variables
if nargin < 2
    error('Not enough inputs.');
elseif nargin < 3
    options = [];
end
if ~isfield(options,'d')
    options.d= r;
end
if ~isfield(options,'T')
    options.T= r;
end
if ~isfield(options,'gamma')
    options.gamma= 0.2;
end
if ~isfield(options,'lambda')
    options.lambda= 10;
end
if ~isfield(options,'eta')
    options.eta= 0.5;
end
if ~isfield(options,'no_show')
    options.no_show= true;
end
if ~isfield(options,'timelimit')
    options.timelimit= 1000;
end
if ~isfield(options,'centerstrategy')
    options.centerstrategy= 'mean';
end
if ~isfield(options,'outlier')
    options.outlier= false;
end
if ~isfield(options,'BIGM')
    options.BIGM=10;
end
d = options.d; T = options.T; gamma = options.gamma; lambda = options.lambda; eta = options.eta; BIG_M = options.BIGM;
% preprocessing : centering data points & dimension reduction
orig_data = X;
MEAN=mean(X,2);
X = X - MEAN;
U = X*X';
[C,~]=eigs(U,d-1);
X = pinv(C)*X;

% initialization 
[m,n] = size(X);
Theta=[];
inner_iter=0;
Indices={};
facets_mean=[];
saved_Points=[];
thetas=[];
outflag = false;

%visualize the dimension reduced data points
if ~options.no_show
    points=X;
    if size(points,1) >= 3
        plot3(points(1,:),points(2,:),points(3,:),'bO'); hold on;
    else
        plot(points(1,:),points(2,:),'bO');hold on;
    end
    axis equal
    hold on;
end
% check if it is rank deficient
rank_def = (T>(m+1)); 

while inner_iter < T
    inner_iter=inner_iter+1;
    if inner_iter < r || rank_def % if not the last facet or if it is rank_def problem
        % parameters of the optimizations are: y (n x 1), delta (n x 1),
        % \theta (m x 1)
        f=[ones(n,1); lambda .* ones(n,1); sparse(m,1)]; %sum (y_i + \lambda delta_i)
        A1=[-BIG_M*eye(n) sparse(n,n) -X']; % X^T * \theta >= 1- \gamma - A * y (BIG_M is A)
        b1=-(1-gamma)*ones(n,1); 
        A2=[sparse(n,n) -eye(n) X' ]; % X^T * \theta < 1 + delta
        b2=(1)*ones(n,1); 
        if ~isempty(facets_mean) % if not the first facet
            A3=[sparse(size(facets_mean,2),2*n) facets_mean']; % M^T \theta <= 1 - gamma - eta
            b3=(1-gamma-eta)*ones(size(facets_mean,2),1);
        else
            A3=[];b3=[];
        end
        if (options.outlier) % if outliers are considered
            A4 = [-BIG_M*eye(n) eye(n) sparse(n,m)]; % delta <= A * y + gamma
            b4 = gamma*ones(n,1);
        else
            A4=[];b4=[];
        end
        A=[A1;A2;A3;A4]; %concat inequality constraints
        b=[b1;b2;b3;b4];
        lb=[0*ones(2*n,1); -Inf*ones(m,1)]; %y>=0, delta>=0
        ub=[1*ones(n,1);Inf*ones(n+m,1)]; %y <1
        Aeq=[];
        beq=[];
        t = [repmat('B',[1 n]),repmat('C',[1 n+m])]; % define y as binary variable, delta and \theta as continuous variables
    else  % if the final facet is going to be selected in the full rank case   
        % parameters of the optimizations are: y (n x 1), delta (n x 1),
        % \theta (m x 1), \mu (r-1 x 1)
        f=[ones(n,1); lambda .* ones(n,1); sparse(m,1);sparse(r-1,1)]; %sum (y_i + \lambda delta_i)
        A1=[-BIG_M*eye(n) sparse(n,n) -X' sparse(n,r-1)]; % X^T * \theta >= 1- \gamma - A * y (BIG_M is A)
        b1=-(1-gamma)*ones(n,1);
        A2=[sparse(n,n) -eye(n) X' sparse(n,r-1)]; % X^T * \theta < 1 + delta
        b2=(1)*ones(n,1);
        A3=[sparse(size(facets_mean,2),2*n) facets_mean' sparse(size(facets_mean,2),r-1)]; % M^T \theta <= 1 - gamma - eta
        b3=(1-gamma-eta)*ones(size(facets_mean,2),1);
        if (options.outlier) % if outliers are considered
            A4 = [-BIG_M*eye(n) eye(n) sparse(n,m+r-1)]; % delta <= A * y + gamma
            b4 = gamma*ones(n,1);
        else
            A4=[];b4=[];
        end
        A=[A1;A2;A3;A4]; % concat inequality constraints
        b=[b1;b2;b3;b4];
        lb=[0*ones(2*n,1); -Inf*ones(m,1);0.1*ones(r-1,1)]; %define lower bounds: y>=0, delta>=0, mu >=0.1
        ub=[1*ones(n,1);Inf*ones(n+m,1);Inf*ones(r-1,1)]; % y <1
        Aeq=[zeros(m,n) zeros(m,n) eye(m) thetas(1:m,:)]; % bounded condition: \theta = - sum \mu * \theta ^(i)
        beq=zeros(m,1);
        t = [repmat('B',[1 n]),repmat('C',[1 n+m+r-1])]; % define y as binary variable, delta, \theta and \mu as continuous variables
    end
    cplex_options = cplexoptimset('cplex'); %get the cplex options default values
    cplex_options.timelimit = options.timelimit; %set the timelimit of cplex
    [x_opt,fval,~,output]=cplexmilp(f,A,b,Aeq,beq,[],[],[],lb,ub,t,[],cplex_options); %optimize using cplex
    if ~outflag && isequal(output.cplexstatusstring,'time limit exceeded') %is timelimit exceeded
        outflag = true;
    end
    display(fval)
    if isempty(fval) % if no feasible solution
        inner_iter = inner_iter - 1; 
        eta = eta - 0.2; %reduce the margin (eta)
    else
        theta=x_opt(2*n+1:2*n+m); 
        ntetha=theta./norm(theta);
        if (size(Theta,1)==0 || (min(rad2deg(acos(Theta(1:m,:)'*ntetha)))>5) ) % if the solution is enough different from previous solutions
            ind=find(abs(1-X'*theta)<=gamma+0.001); %find the points on the facet (added 0.001 for numerical stability)
            if strcmp(options.centerstrategy,'mean') % if the center selection strategy is mean
                current_facet = X(:,ind); %calculate mean of the points on current facet
                current_mean = mean(current_facet,2);
            else % if the center selection strategy is snpa
                current_facet = X(:,ind);
                [J,~] = SNPA(current_facet,r-1);
                current_facet_ = current_facet(:,J);
                current_mean = mean(current_facet_,2); %calculate mean of the points selected by snpa
            end
            if ~options.no_show % if showing intermediate results is requested
                figure(1); % show the points on the selected facet
                points = current_facet;
                try
                plot3(points(1,:),points(2,:),points(3,:),'.','MarkerSize',10,'HandleVisibility','off')
                catch
                plot(points(1,:),points(2,:),'.','MarkerSize',10,'HandleVisibility','off')
                end
                title(strcat('after iteration ',num2str(inner_iter)));
                drawnow;
                pause(1)
            end
            % update 
            facets_mean = [facets_mean current_mean];
            Indices{inner_iter}=ind; %indices of the points on each facet
            Theta=[Theta [ntetha;mean(X(:,ind)'*ntetha)]]; %normalized Thetas
            points_on_current_facet = zeros(n,1);
            points_on_current_facet(abs(1-X'*theta)<gamma)=1;
            saved_Points = [saved_Points points_on_current_facet]; % save the points that are indicated on current facet as a binary matrix
            thetas = [thetas theta]; 
        else % if not very different solution from previous solutions (facets)
            inner_iter = inner_iter - 1;
            eta = eta + 0.2; % increase the margin (eta)
        end
    end
end
%%
num = sum(min(sum(saved_Points,2),1)); % sum the total points that are on all the facets
%%
W_tilde = find_intersection(Indices,X,r); % find intersection 
if ~options.no_show % if displaying the results is requested 
    figure(1); % show the Ws
    if size(W_tilde,1)>=3
        plot3(W_tilde(1,:),W_tilde(2,:),W_tilde(3,:),'d','MarkerSize',20);
        hold on;
    else
        figure(1);
        plot(W_tilde(1,:),W_tilde(2,:),'d','MarkerSize',20);
        hold on;
    end
    drawnow;
end
W = max(C*W_tilde+MEAN,0); % bring back to original space
H = FGMfcnls(orig_data,W); % find H
return;

end