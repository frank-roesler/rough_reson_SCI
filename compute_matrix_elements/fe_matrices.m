function [s,mass,b,vol_T,mp_T] = fe_matrices(c4n,n4e,Nb) 

% Computes the mass and stiffness matrices [s,m,vol_T,mp_T] for
% triangulation given by [c4n,n4e]

[nC,d] = size(c4n);
nE = size(n4e,1);
nNb = size(Nb,1);
m_loc = (ones(d+1,d+1)+eye(d+1))/((d+1)*(d+2));
ctr = 0; 
ctr_max = (d+1)^2*nE;
I = zeros(ctr_max,1); 
J = zeros(ctr_max,1);
X_s = zeros(ctr_max,1); 
X_m = zeros(ctr_max,1); 
vol_T = zeros(nE,1); 
mp_T = zeros(nE,2);
b = zeros(nC,1);
for j = 1:nE
    X_T = [ones(1,d+1);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,d);eye(d)]; 
    vol_T(j) = det(X_T)/factorial(d); 
    mp_T(j,:)= sum(c4n(n4e(j,:),:),1)/(d+1);
    for m = 1:d+1
        b(n4e(j,m)) = b(n4e(j,m)) + (1/(d+1))*vol_T(j)*rhs(mp_T(j,:));
        for n = 1:d+1
            ctr = ctr+1; 
            I(ctr) = n4e(j,m); 
            J(ctr) = n4e(j,n); 
            X_s(ctr) = vol_T(j)*grads_T(m,:)*grads_T(n,:)'; 
            X_m(ctr) = vol_T(j)*m_loc(m,n);
        end
    end
end
s = sparse(I,J,X_s,nC,nC); 
mass = sparse(I,J,X_m,nC,nC); 
for j = 1:nNb
    vol_S = norm(c4n(Nb(j,1),:)-c4n(Nb(j,2),:));
    mp_S  = sum(c4n(Nb(j,:),:),1)/d;
    for k = 1:d
        b(Nb(j,k)) = b(Nb(j,k))+(1/d)*vol_S*N_data(mp_S);
    end
end
end