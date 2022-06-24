function P = Koch(n)
%Koch draws Koch snowflake
%n is how many iterations of the creative process we want 
% initialize P to an equilateral triangle 
P = [ 0 0; 
      1 0; 
      cos(-pi/3), sin(-pi/3); 
      0 0  ];
for iteration=1:n
   newP = zeros( size(P,1)*4+1, 2);
   
   for i=1:size(P,1)-1
       newP(4*i+1,:) = P(i,:);
       newP(4*i+2,:) = (2*P(i,:) + P(i+1,:) )/3;
       link = P(i+1,:)-P(i,:);
       ang = atan2( link(2), link(1) );   
       linkLeng = sqrt( sum(link.^2) );   
       newP(4*i+3,:) = newP(4*i+2,:) + (linkLeng/3)*[ cos(ang+pi/3), sin(ang+pi/3) ];
       newP(4*i+4,:) = (P(i,:) + 2*P(i+1,:) )/3;
   end
   newP( 4*size(P,1)+1,:) = P(size(P,1),:);
   P = newP;
end
% % now join up the points in P
%    clf;        % clear the figure window   
%    plot( P(:,1), P(:,2) );   % plot P
%    axis equal; % make the x- and y-scale the same
   P = P + [-0.5, 1/(2*sqrt(3))];
end