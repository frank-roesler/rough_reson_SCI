clear
addpath('/Users/frankrosler/Documents/MATLAB/matplotlib')
mesh_path = 'meshes/julia0.02.mat';
upper_bound_eigs = 60;
[W,spectrum,c4n,n4e,fNodes,Nb,r_ball] = compute_eigenfunctions(mesh_path,upper_bound_eigs);

%% Precompute M_in(k_0):
N = 33; % half size of matrix M_in^n(k)
k_0 = 2-1i;
M_in_direct_0 = compute_matrix_fem(k_0,N,c4n,n4e,Nb,fNodes,r_ball);

%% Compute det(M_in + M_out) in complex plane:
z1 = -2i;
z2 = 4;
hres = 1000;
L = build_lattice(z1, z2, hres);
dets = zeros(size(L));

cN = diag(abs(-N:N)); cN(N+1,N+1)=1; cN = sqrt(cN);

parfor i=1:length(L(:))
    k = L(i);
    corrector = atkinson(k_0,k,c4n,Nb,N,W,spectrum,r_ball);
    M_in = M_in_direct_0 + corrector;
%     M_in = compute_matrix_fem(k,N,c4n,n4e,Nb,fNodes,r_ball);
    M_out = compute_M_outer(k,N,r_ball);
    A = cN*(M_in + M_out)*cN/(2*r_ball);
    dets(i) = det(A);
end

%% Plot results:
figure
contour(real(L),imag(L),log(abs(dets)),100)
colormap viridis
colorbar
hold on
resonances = islocalmin(abs(dets),1) & islocalmin(abs(dets),2);
plot(L(resonances),'.','MarkerSize',20,'color','r')
hold off
L(resonances)

%%
figure
k_test = 2-0.5i;
M_in_direct_t = compute_matrix_fem(k_test,N,c4n,n4e,Nb,fNodes,r_ball);
M_out_t = compute_M_outer(k_test,N,r_ball);
A = cN*(M_in_direct_t + M_out_t)*cN/(2*r_ball);
plot(abs(diag(A-eye(2*N+1))))



























