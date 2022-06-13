% -------------------------------------------------------------------------
% VARIABLES:
% N:             half the size of the matrix M_in^n(k)
% k_0:           point in the complex plane where M_in is precomputed via FEM
% M_in_direct_0: precomputed matrix for M_in at k_0
% cN:            diagonal matrix eith entries max(1,sqrt(|n|))
% M_in:          matrix M_in^n(k), computed via Atkinson acceleration
% dets:          determinant of cN(M_in+M_out)cN/(2r_ball) 
%                (i.e. resonances are zeros of dets)
% -------------------------------------------------------------------------

clear
addpath('/Users/frankrosler/Documents/MATLAB/matplotlib')
mesh_path = 'meshes/helmholtz0.05.mat';
upper_bound_eigs = 80;
[W,spectrum,c4n,n4e,fNodes,Nb,r_ball] = compute_eigenfunctions(mesh_path,upper_bound_eigs);
disp(['Eigenfunctions done. Found ',num2str(size(W,2)),' eigenfunctions.'])

%% Precompute M_in(k_0):
N = 30; 
k_0 = 2-0.05i;
M_in_direct_0 = compute_M_in_fem(k_0,N,c4n,n4e,Nb,fNodes,r_ball);
cN = diag(abs(-N:N)); cN(N+1,N+1)=1; cN = sqrt(cN);

%% Compute det(M_in + M_out) in complex plane:
z1 = 1-0.1i;
z2 = 3;
hres = 500;
L = build_lattice(z1, z2, hres);
dets = zeros(size(L));

parfor i=1:length(L(:))
    k = L(i);
    corrector = atkinson_corrector(k_0,k,c4n,Nb,N,W,spectrum,r_ball);
    M_in = M_in_direct_0 + corrector;
    M_out = compute_M_outer(k,N,r_ball);
    A = cN*(M_in + M_out)*cN/(2*r_ball);
    dets(i) = det(A);
end

%% Plot results:
figure('Position',[500,600,1000,300])
contour(real(L),imag(L),log(abs(dets)),100)
colormap viridis
colorbar
hold on
resonances = islocalmin(abs(dets),1) & islocalmin(abs(dets),2);
plot(L(resonances),'.','MarkerSize',10,'color','r')
hold off
Resonances = L(resonances)
drawnow

%%
figure
k_test = -2-0.2i;
M_in_direct_t = compute_M_in_fem(k_test,N,c4n,n4e,Nb,fNodes,r_ball);
M_out_t = compute_M_outer(k_test,N,r_ball);
A = cN*(M_in_direct_t + M_out_t)*cN/(2*r_ball);
plot(abs(diag(A-eye(2*N+1))))
ylim([0,1])


