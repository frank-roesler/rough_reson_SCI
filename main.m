%----------------------------------------------------------------------------
% VARIABLES:                                                                  
% N:             half the size of the matrix M_in^n(k)
% k_0:           point in the complex plane where M_in is precomputed via FEM
% M_in_direct_0: precomputed matrix for M_in at k_0
% cN:            diagonal matrix eith entries max(1,sqrt(|n|))
% M_in:          matrix M_in^n(k), computed via Atkinson acceleration
% dets:          determinant of cN(M_in+M_out)cN/(2r_ball) 
%                (i.e. resonances are zeros of dets)
%----------------------------------------------------------------------------

clear
addpath('/Users/frankrosler/Documents/MATLAB/matplotlib')
mesh_path = 'meshes/koch4_uniform.mat';
upper_bound_eigs = 1500;
[W,spectrum,c4n,n4e,fNodes,Nb,r_ball] = compute_eigenfunctions(mesh_path,upper_bound_eigs);
disp(['Eigenfunctions done. Found ',num2str(size(W,2)),' eigenfunctions.'])

% patch('vertices',c4n,'faces',Nb,'edgecol','k','LineWidth',0.001,'facecol',[.8,.9,1]);

%% Precompute M_in(k_0):
N = 15; 
k_0 = 0-1i;
M_in_direct_0 = compute_M_in_fem(k_0,N,c4n,n4e,Nb,fNodes,r_ball);
cN = diag(abs(-N:N)); cN(N+1,N+1)=1; cN = sqrt(cN);

%% Compute det(M_in+M_out) in complex plane:
z1 = -5-3i;
z2 = 5-0.02i;
hres = 500;
L = build_lattice(z1, z2, hres);
dets = zeros(size(L));
det_out = zeros(size(L));
parfor i=1:length(L(:))
    k = L(i);
    corrector = atkinson_corrector(k_0,k,c4n,Nb,N,W,spectrum,r_ball);
    M_in = M_in_direct_0 + corrector;
    M_out = compute_M_outer(k,N,r_ball);
    det_out(i) = det(cN*M_out*cN/r_ball); % for comparison and regularisation
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

% % Regularise determinant by dividing by 1+det_out:
% figure('Position',[500,600,1000,300])
% contour(real(L),imag(L),log(abs(det_out)),100)
% colorbar

%% Plot diagonal entries of K(z). Should go to 0.
figure
N = 20;
k_test = -1-0.1i;
cN = diag(abs(-N:N)); cN(N+1,N+1)=1; cN = sqrt(cN);
M_in_direct_t = compute_M_in_fem(k_test,N,c4n,n4e,Nb,fNodes,r_ball);
M_out_t = compute_M_outer(k_test,N,r_ball);
A = cN*(M_in_direct_t + M_out_t)*cN/(2*r_ball);
plot(abs(diag(A-eye(2*N+1))),'LineWidth',1)
ylim([0,1])


