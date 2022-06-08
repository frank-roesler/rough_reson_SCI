function n4e_new = join_inner_outer_mesh(c4n,n4e, R,h)
[r_left1,c_left1] = find(abs(c4n(:,1)+(R-h))<h/1000);
[r_left2,c_left2] = find(abs(c4n(:,1)+R)<h/1000);
[c_left1,perm1] = sort(c4n(r_left1,:));
[c_left2,perm2] = sort(c4n(r_left2,:));
r_left1 = r_left1(perm1(:,2));
r_left2 = r_left2(perm2(:,2));
%
left_triangles = zeros(2*length(r_left1)-2, 3);
for i=1:length(r_left1)-1
    Tu = [r_left1(i+1), r_left2(i), r_left1(i)];
    Tl = [r_left1(i+1), r_left2(i+1), r_left2(i)];
    left_triangles(i,:) = Tu;
    left_triangles(i+length(r_left1)-1,:) = Tl;
end
%

[r_right1,c_right1] = find(abs(c4n(:,1)-(R-h))<h/1000);
[r_right2,c_right2] = find(abs(c4n(:,1)-R)<h/1000);
[c_right1,perm1] = sort(c4n(r_right1,:));
[c_right2,perm2] = sort(c4n(r_right2,:));
r_right1 = r_right1(perm1(:,2));
r_right2 = r_right2(perm2(:,2));
%
right_triangles = zeros(2*length(r_right1)-2, 3);
for i=1:length(r_right1)-1
    Tu = [r_right1(i), r_right2(i+1), r_right1(i+1)];
    Tl = [r_right1(i), r_right2(i), r_right2(i+1)];
    right_triangles(i,:) = Tu;
    right_triangles(i+length(r_right1)-1,:) = Tl;
end
%

[r_up1,c_up1] = find(abs(c4n(:,2)-(R-h))<h/1000 & abs(c4n(:,1))<R-h/2);
[r_up2,c_up2] = find(abs(c4n(:,2)-R)<h/1000 & abs(c4n(:,1))<R-h/2);
[c_up1,perm1] = sort(c4n(r_up1,:));
[c_up2,perm2] = sort(c4n(r_up2,:));
r_up1 = r_up1(perm1(:,2));
r_up2 = r_up2(perm2(:,2));
%
up_triangles = zeros(2*length(r_up1)-2, 3);
for i=1:length(r_up1)-1
    Tu = [r_up1(i), r_up1(i+1), r_up2(i+1)];
    Tl = [r_up1(i), r_up2(i+1), r_up2(i)];
    up_triangles(i,:) = Tu;
    up_triangles(i+length(r_up1)-1,:) = Tl;
end
%

[r_low1,c_low1] = find(abs(c4n(:,2)+(R-h))<h/1000 & abs(c4n(:,1))<R-h/2);
[r_low2,c_low2] = find(abs(c4n(:,2)+R)<h/1000 & abs(c4n(:,1))<R-h/2);
[c_low1,perm1] = sort(c4n(r_low1,:));
[c_low2,perm2] = sort(c4n(r_low2,:));
r_low1 = r_low1(perm1(:,2));
r_low2 = r_low2(perm2(:,2));
%
low_triangles = zeros(2*length(r_low1)-2, 3);
for i=1:length(r_low1)-1
    Tu = [r_low1(i), r_low2(i+1), r_low1(i+1)];
    Tl = [r_low1(i), r_low2(i), r_low2(i+1)];
    low_triangles(i,:) = Tu;
    low_triangles(i+length(r_low1)-1,:) = Tl;
end

%%
n4e_new = [n4e ; left_triangles; right_triangles; up_triangles; low_triangles];