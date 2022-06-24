clear
addpath /Users/frankrosler/Documents/MATLAB/distmesh-Original
% addpath /Users/frankrosler/Documents/MATLAB/Poisson_Finite_Elements

d = 2;
r_ball = 1.85;
h = 0.01;
npts = round(1/h)+1;
Rx = 1.46+h/2;
Ry = 0.95+h/2;
L = build_lattice_mesh(-Rx-Ry*1i, Rx+Ry*1i, h);
M = ~in_Julia(L, 1000, 2);

fd=@(p) ddiff(dcircle(p,0,0,r_ball), drectangle(p,-Rx,Rx,-Ry,Ry)); 

box=[-r_ball,-r_ball;r_ball,r_ball];
fix=[real(L(1,:))' , imag(L(1,:))'  ;  
    real(L(end,:))', imag(L(end,:))'; 
    real(L(:,1))   , imag(L(:,1))   ; 
    real(L(:,end)) , imag(L(:,end))];

disp('distmeshing...')
[c4n_outer,n4e_outer]=distmesh2d(fd,@huniform,h,box,fix, false); 
disp('distmesh done!')
close all
[c4n_outer,n4e_outer]=fixmesh(c4n_outer,n4e_outer);
nC = length(c4n_outer);
nE = length(n4e_outer);

%
[c4n_inner, n4e_inner] = build_mesh_reson(h,Rx-h/2,Ry-h/2);

[c4n_inner_new, n4e_inner_new] = remove_inner_boundary(c4n_inner,n4e_inner,c4n_outer,h);
c4n = [c4n_outer; c4n_inner_new];
n4e = [n4e_outer; n4e_inner_new + max(n4e_outer(:))];

n4e_filled_left_right = join_inner_outer_mesh(c4n,n4e, Rx,Ry,h);
n4e = n4e_filled_left_right;
%% Plot Mesh:
TR = triangulation(n4e_filled_left_right,c4n);
figure
triplot(TR);
xlim([-r_ball-.1,r_ball+.1])
ylim([-r_ball-.1,r_ball+.1])

%%
save(['../meshes/julia',num2str(h),'.mat'], 'L', 'h', 'M', 'c4n', 'n4e','r_ball');











