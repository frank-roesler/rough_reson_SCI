clear

XY = Koch(6);
XY = unique(XY,'rows','stable');
XY = [XY ; XY(1,:)];
X = XY(:,1);
Y = XY(:,2);
plot(X,Y)

% [X,Y] = interpolate_koch(X,Y);
% XY = [X,Y];

%%
% H = Koch_curve_GUI;
% %%
% Xf = start_data(:,1)-144;
% Yf = start_data(:,2)-137;
% Xf = Xf(~isnan(Yf))/40;
% Yf = Yf(~isnan(Yf))/40;
% XYf = [Xf,Yf];
% XYf = unique(XYf,'rows','stable');
% Xf = XYf(:,1);
% Yf = XYf(:,2);

%%
addpath /Users/frankrosler/Documents/MATLAB/distmesh-Original

p1 = XY(1:end-1,:);
p2 = XY(2:end,:);
deltas = vecnorm(p1-p2,2,2);
h = mean(deltas)
%
r_ball = 0.7;

dkoch =@(p) inkoch(p,X,Y);
fd=@(p)  ddiff(dcircle(p,0,0,r_ball), dkoch(p)); 

adap =@(p) 0.01*(distskoch(p,X,Y)./distskoch(r_ball*p./(vecnorm(p,2,2)),X,Y)).^2 + h/2;

fix = XY;
box=[-r_ball,-r_ball;r_ball,r_ball];
disp('distmeshing...')
tic
[c4n,n4e] = distmesh2d(fd,adap,h,box,fix, true); 
toc

hold on
plot(fix(:,1), fix(:,2),'.','color','r')
hold off

%
save(['/Users/frankrosler/Desktop/rough_resonance_matlab/meshes/koch6_uniform',num2str(h),'.mat']);
%%

function val = distskoch(p,X,Y)
%     [X,Y] = interpolate_koch(X,Y);
%     [X,Y] = interpolate_koch(X,Y);
    dists = sqrt(pdist2([X,Y],p,'euclidean','Smallest',1).');
    dists(inpolygon(p(:,1),p(:,2),X,Y)) = -dists(inpolygon(p(:,1),p(:,2),X,Y));
    val = dists;
end

function dists = inkoch(p,X,Y)
    dists = zeros(size(p,1),1);
    in = inpolygon(p(:,1),p(:,2),X,Y);
    dists(in) = -sqrt(pdist2([X,Y],p(in,:),'euclidean','Smallest',1).');
    dists(~in) = sqrt(pdist2([X,Y],p(~in,:),'euclidean','Smallest',1).');
end

function [X_interp,Y_interp] = interpolate_koch(X,Y)
    X_interp = zeros(2*size(X,1)-1,size(X,2));
    Y_interp = zeros(2*size(X,1)-1,size(X,2));
    X_midp = (X(1:end-1,:)+X(2:end,:))/2;
    Y_midp = (Y(1:end-1,:)+Y(2:end,:))/2;
    X_interp(1:2:end) = X;
    X_interp(2:2:end-1) = X_midp;
    Y_interp(1:2:end) = Y;
    Y_interp(2:2:end-1) = Y_midp;
end

