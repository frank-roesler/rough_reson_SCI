function [W,spectrum,c4n,n4e,fNodes,Nb,r_ball] = compute_eigenfunctions(mesh_path,upper_bound_eigs)
    % Load and plot mesh:
    load(mesh_path)
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
    [s,m] = fe_matrices(c4n,n4e,Nb);
    S = s(fNodes,fNodes);
    M = m(fNodes,fNodes);
    [V, spectrum, iresult] = sptarn(S,M,0,upper_bound_eigs,1,'jmax',300);
    W = zeros(nC,size(V,2));
    W(fNodes,:) = V;
    W = W./sqrt(diag(W'*m*W))';
end