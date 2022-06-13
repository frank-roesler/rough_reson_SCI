function [c4n_inner_new, n4e_inner_new] = remove_inner_boundary(c4n_inner,n4e_inner,c4n_outer,h)
join_pts = [];
for i=1:size(c4n_outer,1)
    eq = abs(c4n_inner-c4n_outer(i,:))<h/10;
    eq = eq(:,1) & eq(:,2);
    join_pt = find(eq==1);
    join_pts = [join_pts, join_pt];
end

c4n_inner_new = c4n_inner;
n4e_inner_new = n4e_inner;
for i=1:length(join_pts)
    p1 = join_pts(i);
    [c4n_inner_new, n4e_inner_new] = delete_row(c4n_inner_new, n4e_inner_new, p1);
    join_pts = join_pts-1;
end

% subplot(2,1,1)
% TR = triangulation(n4e_inner,c4n_inner);
% triplot(TR);
% xlim([-R,R])
% ylim([-R,R])
% subplot(2,1,2)
% TR = triangulation(n4e_inner_new,c4n_inner_new);
% triplot(TR);
% xlim([-R,R])
% ylim([-R,R])
% drawnow
end