clear

mesh_path = 'meshes/julia0.02.mat';
upper_bound_eigs = 50;

[W,spectrum,c4n,n4e,fNodes,Nb,r_ball] = compute_eigenfunctions(mesh_path,upper_bound_eigs);


%% Precompute M_in(k_0):
N = 50;
k_0 = 0.5-0.2i;
M_in_direct_0 = compute_matrix_fem(k_0,N,c4n,n4e,Nb,fNodes,r_ball);

%% Compute det(M_in + M_out) in complex plane:
z1 = 0-0.4i;
z2 = 1;
hres = 1000;
L = build_lattice(z1, z2, hres);
dets = zeros(size(L));

cN = diag(abs(-N:N));
cN(N+1,N+1)=1;
cN = sqrt(cN);

parfor i=1:length(L(:))
    k = L(i);
    corrector = atkinson(k_0,k,c4n,Nb,N,W,spectrum,r_ball);
    M_in = M_in_direct_0 + corrector;
    M_out = compute_M_outer(k,N,r_ball);
    A = cN*(M_in + M_out)*cN/(2*r_ball);
    dets(i) = det(A);
end

%%
figure
contour(real(L),imag(L),log(abs(dets)),100)
colorbar

%%
figure
M_out_0 = compute_M_outer(k_0,N,r_ball);
A = cN*(M_in_direct_0 + M_out_0)*cN/(2*r_ball);
surf(abs(A-eye(2*N+1)))



























