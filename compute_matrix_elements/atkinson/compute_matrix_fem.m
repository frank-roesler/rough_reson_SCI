function M = compute_matrix_fem(k,N,c4n,n4e,Nb,fNodes,r_ball)
% Computes matrix approximation for M_inner at specified value of k
% directly from the FEM solution u.
    [nC,d] = size(c4n);
    M = zeros(2*N+1);
    
    x = c4n(Nb(:,1),:);
    dx = vecnorm(x - c4n(Nb(:,2),:),2,2);
    e_bdry = exp(1i*(-N:N).*atan2(x(:,2),x(:,1)))/sqrt(2*pi*r_ball);
    
    [s,m,b,vol_T,mp_T] = fe_matrices(c4n,n4e,Nb);
    S = s - k^2*m;   % weak version of -∆-k^2
    for alpha=-N:N
        b = zeros(nC,1);
        for j = 1:size(Nb,1)
            vol_S = norm(c4n(Nb(j,1),:)-c4n(Nb(j,2),:));
            mp_S  = sum(c4n(Nb(j,:),:),1)/d;
            for k = 1:d
                b(Nb(j,k)) = b(Nb(j,k))+(1/d)*vol_S*e_alpha(alpha,r_ball,mp_S);
            end
        end
        u         = zeros(nC,1);
        u(fNodes) = S(fNodes,fNodes)\b(fNodes); 
        u_bdry    = u(Nb(:,1));
%         for beta=-N:N
%             e_beta = e_bdry(:,beta+N+1);
%             M(alpha+N+1,beta+N+1) = sum(conj(e_beta).*u_bdry.*dx);
%         end
        M(alpha+N+1,:) = (conj(e_bdry).*u_bdry).'*dx;
    end
end
