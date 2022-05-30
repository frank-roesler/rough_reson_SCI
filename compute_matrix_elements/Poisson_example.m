clear
addpath('matplotlib')

% Load and plot mesh:
load('julia0.02.mat')
n4e = n4e_filled_left_right;
TR = triangulation(n4e,c4n);
figure
triplot(TR);
drawnow

% Construct Dirichlet and Neumann boundaries:
boundary = freeBoundary(TR);
b1 = boundary(:,1);
b2 = boundary(:,2);
db1 = vecnorm(c4n(b1,:),2,2)<0.9*r_ball;
db2 = vecnorm(c4n(b2,:),2,2)<0.9*r_ball;
Db = [b1(db1),b2(db2)];
Nb = [b1(~db1),b2(~db2)];

% % Plot Dirichlet boundary:
% figure
% patch('vertices',c4n,'faces',Db,'edgecol','k','LineWidth',0.1,'facecol',[.8,.9,1]);
% axis off
% drawnow

[nC,d]  = size(c4n);            % number of nodes
nE      = size(n4e,1);          % number of elements
dNodes  = unique(Db);           % Dirichlet boundary
fNodes  = setdiff(1:nC,dNodes); % free nodes


% Solve Dirichlet Poisson Problem:
[s,m,b,vol_T,mp_T] = fe_matrices(c4n,n4e,Nb);

max_u = [];
figure
for k=-10:0.1:10
%     Solve problem for M_inner:
    k = k - 0.1i
    u         = zeros(nC,1);
    S         = s - k^2*m;   % weak version of -∆+1
    u(fNodes) = S(fNodes,fNodes)\b(fNodes); 

    % Plot Results:
    subplot(1,2,1)
    trisurf(n4e,c4n(:,1),c4n(:,2),real(u),'LineWidth',0.01,'EdgeColor','none');
    title('Real Part')
    colormap('viridis')
    view([0,90])
    caxis([-0.5,0.5])
    subplot(1,2,2)
    trisurf(n4e,c4n(:,1),c4n(:,2),imag(u),'LineWidth',0.01,'EdgeColor','none');
    title('Imaginary Part')
    colormap('viridis')
    view([0,90])
    caxis([-0.5,0.5])
    drawnow
    
    max_u = [max_u,max(abs(u(:)))];
end

figure
plot(max_u)























