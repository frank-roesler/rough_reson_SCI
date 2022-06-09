function M_in = compute_matrix_fem(k,N,c4n,n4e,Nb,fNodes,r_ball)
% Computes matrix approximation for M_inner at specified value of k
% directly from the FEM solution u.
    [nC,d] = size(c4n);
    M_in = zeros(2*N+1);
    
    x = c4n(Nb(:,1),:); % boundary of disk
    dx = vecnorm(x - c4n(Nb(:,2),:),2,2);
    e_bdry = exp(1i*(-N:N).*atan2(x(:,2),x(:,1)))/sqrt(2*pi*r_ball);
    
    [s,m,b,vol_T,mp_T] = fe_matrices(c4n,n4e,Nb);
    S = s - k^2*m;   % weak version of -∆-k^2
    for alpha=-N:N
        b = zeros(nC,1);
        for j = 1:size(Nb,1)
            vol_S = norm(c4n(Nb(j,1),:)-c4n(Nb(j,2),:));
            mp_S  = sum(c4n(Nb(j,:),:),1)/d;
            for i = 1:d
                b(Nb(j,i)) = b(Nb(j,i))+(1/d)*vol_S*e_alpha(alpha,r_ball,mp_S);
            end
        end
        u         = zeros(nC,1);
        u(fNodes) = S(fNodes,fNodes)\b(fNodes); 
        u_bdry    = (u(Nb(:,1)) + u(Nb(:,2)))/2;
        M_in(alpha+N+1,:) = (e_bdry.*dx)'*u_bdry;
    end
    M_in = M_in.';
end
