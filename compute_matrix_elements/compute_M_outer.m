function Mout = compute_M_outer(k,N,X)
%     Computes matrix elements of the outer NtD map
    alpha = -N:N;
    d = k*besselh(abs(alpha)-1,X*k)./besselh(abs(alpha),X*k) - abs(alpha)/X;
    Mout = diag(-1./d);
end