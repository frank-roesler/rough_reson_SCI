clear
clc
addpath /Users/frankroesler/Documents/MATLAB/distmesh-Original

Rx=1.5;
Ry=1;
r_ball = 2;
h=0.01;
box=[-r_ball,-r_ball;r_ball,r_ball];
L = build_lattice_mesh(-Rx-Ry*1i, Rx+Ry*1i, h);
M = in_Julia(L, 100, 2);

XX = real(L(M));
YY = imag(L(M));
bdry = boundary(XX,YY,1);
XY = [XX(bdry),YY(bdry)];
% XY = unique(XY,'rows','stable');
X = XY(:,1);
Y = XY(:,2);

[Xi,Yi] = interpolate_poly(X,Y);
% [Xi,Yi] = interpolate_poly(Xi,Yi);
XYi = [Xi,Yi];

XY = unique(XY,'rows','stable');
X = XY(:,1);
Y = XY(:,2);
XYi = unique(XYi,'rows','stable');
Xi = XYi(:,1);
Yi = XYi(:,2);

%%
fd=@(p) ddiff(dcircle(p,0,0,r_ball), distsjulia(p,Xi,Yi));
fh=@(p) h*(0.5 + 5*distsjulia(p,Xi,Yi).^2);

fix = XY;

disp('distmeshing...')
[c4n,n4e] = distmesh2doriginal(fd,fh,h,box,fix); 
disp('distmesh done!')
%%

TR = triangulation(n4e,c4n);
figure
triplot(TR);

%%

function dists = distsjulia(p,X,Y)
    dists = zeros(size(p,1),1);
    in = inpolygon(p(:,1),p(:,2),X,Y);
    dists(in) = -pdist2([X,Y],p(in,:),'euclidean','Smallest',1).';
    dists(~in) = pdist2([X,Y],p(~in,:),'euclidean','Smallest',1).';
end


function [X_interp,Y_interp] = interpolate_poly(X,Y)
    X_interp = zeros(2*size(X,1)-1,size(X,2));
    Y_interp = zeros(2*size(X,1)-1,size(X,2));
    X_midp = (X(1:end-1,:)+X(2:end,:))/2;
    Y_midp = (Y(1:end-1,:)+Y(2:end,:))/2;
    X_interp(1:2:end) = X;
    X_interp(2:2:end-1) = X_midp;
    Y_interp(1:2:end) = Y;
    Y_interp(2:2:end-1) = Y_midp;
end
