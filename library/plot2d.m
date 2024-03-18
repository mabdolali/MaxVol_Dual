function plot2d(X,V,in)
% X : data point, m-by-n matrix
% V : vertices,   m-by-r matrix
% if X,V are alread 2D (or projected onto 2D) then plot them directly
% otherwise X,V are projected onto 2D by PCA
% in: plotting parameters
%  in.plotall 
%   in.plotall = 1 means plot all data point, vertices and hull
%   in.plotall = 0 means only plot vertices and hul
%   in.plotall default is 0
%  in.plotspec 
%  in.plotspec = {'data point color and plot specification','vertex color and plot specification','hull line specification'}
%  e.g. in.plotspec = {'k.','rx','r--'} means data points are black dots,
%  vertices are red cross, hull are red broken lines
%% Input handling
m = size(X,1);
if m == 1 
   error('X need to have at least 2 rows.');  
end
%% Check need project or not
if m > 2 % project to 2D   
 [X V] = proj2d(X, V);
end
%% Plot
plotThe2d(X,V,in);
end % EOF

%% Projectiop to 2D
function [X_proj V_proj] = proj2d(X, V)
% 2D projection by PCA
% X : data matrix
% V : vertices, a matrix

% *** eigs use random starts, solution may not be unique 
[m n]= size(X);
r = size(V,2);
OPTS.disp = 0; 
OPTS.v0 = ones(m,1); 
d = mean(X,2);
U = X-d*ones(1,n);
[C D] = eigs(U*U',r-1,'LM',OPTS);
X_proj = C'*(X-d*ones(1,n));
V_proj = C'*(V-d*ones(1,r));    
end %EOF

%% Plot 2d
function plotThe2d(X,V,in)
% X : 2d projected data matrix
% V : 2d projected vertices, a matrix
% in: plotting parameters
%  in.pltall 
%   in.plotall = 1 means plot all data point, vertices and hull
%   in.plotall = 0 means only plot vertices and hul
%   in.plotall default is 0
%  in.plotspec 
%  in.plotspec = {'data point color and plot specification','vertex color and plot specification','hull line specification'}
%  e.g. in.plotspec = {'k.','rx','r--'} means data points are black dots,
%  vertices are red cross, hull are red broken lines
%% Input handling
if nargin <= 2
   in = [];
end

if ~isfield(in,'plotall')
    in.plotall = 0 ;
end

if ~isfield(in,'plotspec')
    in.plotspec = {'k.','k*','k-.'};
end

if ~isfield(in,'showlabel')
    in.showlabel = 0;
end
%% 
 % colorcode : line specification, a string {'data point color','vertex color'}
   c_data   = in.plotspec{1};
   c_vertex = in.plotspec{2};
   c_hull   = in.plotspec{3};
%% plot 

% data points
 if in.plotall == 1
   % plot(X(1,:),X(2,:),c_data,'Color',[0.55 0.55 0.55],'markersize',2), hold on, 
   plot(X(1,:),X(2,:),c_data,'markersize',5,'HandleVisibility','off'), hold on,  
 end

% plot vertex
 plot(V(1,:),V(2,:),c_vertex,'markersize',13), hold on

% plot hull
 hullid = convhull(V(1,:),V(2,:));
 plot( V(1,hullid), V(2,hullid), c_hull,'linewidth',3,'HandleVisibility','off')

% some code to show column order
if in.showlabel == 1
hold on
for i  = 1 : size(V,2)
 text(V(1,i),V(2,i),num2str(i),'fontsize',20);
end
end

grid on
end % EOF