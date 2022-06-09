clear
addpath /Users/frankrosler/Documents/MATLAB/distmesh-Original
% addpath /Users/frankrosler/Documents/MATLAB/Poisson_Finite_Elements

d = 2;
r_ball = 2.5;
r = 0.2;
h = 0.1;
npts = round(1/h)+1;
R = 1.5+h/2;
L = build_lattice(-R*(1+1i), R*(1+1i), h);
M = ~in_Julia(L, 1000, 2);

fd=@(p) ddiff(dcircle(p,0,0,r_ball), drectangle(p,-R,R,-R,R)); 

box=[-r_ball,-r_ball;r_ball,r_ball];
fix=[real(L(1,:))' , imag(L(1,:))'  ;  
    real(L(end,:))', imag(L(end,:))'; 
    real(L(:,1))   , imag(L(:,1))   ; 
    real(L(:,end)) , imag(L(:,end))];

disp('distmeshing...')
[c4n_outer,n4e_outer]=distmesh2d(fd,@huniform,h,box,fix, true); 
disp('distmesh done!')
close all
[c4n_outer,n4e_outer]=fixmesh(c4n_outer,n4e_outer);
nC = length(c4n_outer);
nE = length(n4e_outer);

%
[c4n_inner, n4e_inner] = build_mesh_reson(h);

[c4n_inner_new, n4e_inner_new] = remove_inner_boundary(c4n_inner,n4e_inner,c4n_outer,h,R);
c4n = [c4n_outer; c4n_inner_new];
n4e = [n4e_outer; n4e_inner_new + max(n4e_outer(:))];

n4e_filled_left_right = join_inner_outer_mesh(c4n,n4e, R,h);
n4e = n4e_filled_left_right;
%% Plot Mesh:
TR = triangulation(n4e_filled_left_right,c4n);
figure
triplot(TR);
xlim([-r_ball-.1,r_ball+.1])
ylim([-r_ball-.1,r_ball+.1])

%%
save(['julia_test',num2str(h),'.mat'], 'L', 'h', 'M', 'c4n', 'n4e','r_ball');











