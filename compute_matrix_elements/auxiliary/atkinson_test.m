clear

upper_bound_eigs = 20;
k_0 = -1i;
k = 1-1i;
N = 10;

% Load and plot mesh:
load('meshes/helmholtz0.05.mat')
n4e = n4e_filled_left_right;
% Construct Dirichlet and Neumann boundaries:
TR = triangulation(n4e,c4n);clc

boundary = freeBoundary(TR);
b1 = boundary(:,1);
b2 = boundary(:,2);
db1 = vecnorm(c4n(b1,:),2,2)<0.9*r_ball;
db2 = vecnorm(c4n(b2,:),2,2)<0.9*r_ball;
Db = [b1(db1),b2(db2)];
Nb = [b1(~db1),b2(~db2)];
[nC,d]  = size(c4n);            % number of nodes
nE      = size(n4e,1);          % number of elements
dNodes  = unique(Db);           % Dirichlet boundary
fNodes  = setdiff(1:nC,dNodes); % free nodes
% Solve Dirichlet Poisson Problem:
[s,m,b,vol_T,mp_T] = fe_matrices(c4n,n4e,Nb);
S = s(fNodes,fNodes);
M = m(fNodes,fNodes);
[V, spectrum, iresult] = sptarn(S,M,0,upper_bound_eigs,1,'jmax',300);
W = zeros(nC,size(V,2));
W(fNodes,:) = V;
W = W./sqrt(diag(W'*m*W))';

M_in_atk = atkinson(k_0,k,c4n,Nb,N,W,spectrum,r_ball);
M_in_direct_0 = compute_M_in_fem(k_0,N,c4n,n4e,Nb,fNodes,r_ball);
M_in_direct_k = compute_M_in_fem(k,N,c4n,n4e,Nb,fNodes,r_ball);

%
M_in_appr = M_in_direct_0 + M_in_atk;

error = max(max(abs(M_in_appr - M_in_direct_k)))

%%
subplot(1,2,1)
surf(-N:N,-N:N,real(M_in_appr))
subplot(1,2,2)
surf(-N:N,-N:N,imag(M_in_appr))

figure
subplot(1,2,1)
surf(-N:N,-N:N,real(M_in_direct_k))
subplot(1,2,2)
surf(-N:N,-N:N,imag(M_in_direct_k))

%%

wa = zeros(nC,1);
wa(fNodes) = Wa(:,19);
trisurf(n4e,c4n(:,1),c4n(:,2),real(wa),'LineWidth',0.01,'EdgeColor','none');

figure
ws = zeros(nC,1);
ws(fNodes) = -Ws(:,19);
trisurf(n4e,c4n(:,1),c4n(:,2),real(ws),'LineWidth',0.01,'EdgeColor','none');

max(abs(wa-ws))















% %% Plot eigenfunctions:
% 
% for j=1:size(V,2)
%     j
%     e = zeros(nC,1);
%     e(fNodes) = V(:,j);
%     e = e/sqrt(e'*(m*e));
%     
%     subplot(1,2,1)
%     trisurf(n4e,c4n(:,1),c4n(:,2),real(e), 'LineStyle','none');
%     zlim([-2,2])
%     subplot(1,2,2)
%     patch('vertices',c4n,'faces',n4e,'FaceVertexCData',real(e),'FaceColor','flat','EdgeColor','none')
%     zlim([-2,2])
%     drawnow
%     pause(0.5)
% end
















