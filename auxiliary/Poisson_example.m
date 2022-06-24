clear

% Load and plot mesh:
load('meshes/julia0.02.mat')
TR = triangulation(n4e,c4n);
figure
triplot(TR);
drawnow
size(c4n,1)

%% Construct Dirichlet and Neumann boundaries:
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
[s,m] = fe_matrices(c4n,n4e); % old version; doesn't work anymore.
b = zeros(nC,1);
        for j = 1:size(Nb,1)
            vol_S = norm(c4n(Nb(j,1),:)-c4n(Nb(j,2),:));
            mp_S  = sum(c4n(Nb(j,:),:),1)/d;
            for i = 1:d
                b(Nb(j,i)) = b(Nb(j,i))+(1/d)*vol_S*e_alpha(5,r_ball,mp_S);
            end
        end
max_u = [];
figure
k_values = -10:0.1:10;
for k=k_values
%     Solve problem for M_inner:
    k = k - 0.1i
    u         = zeros(nC,1);
    S         = s - k^2*m;   % weak version of -∆+1
    u(fNodes) = S(fNodes,fNodes)\b(fNodes); 

    % Plot Results:
    subplot(1,2,1)
    trisurf(n4e,c4n(:,1),c4n(:,2),real(u),'LineWidth',0.01,'EdgeColor','none');
    title('Real Part')
    view([0,90])
    caxis([-1,1])
    subplot(1,2,2)
    trisurf(n4e,c4n(:,1),c4n(:,2),imag(u),'LineWidth',0.01,'EdgeColor','none');
    title('Imaginary Part')
    view([0,90])
    caxis([-1,1])
    drawnow
    
    max_u = [max_u,max(abs(u(:)))];
end

%%
figure
plot(k_values,max_u)























