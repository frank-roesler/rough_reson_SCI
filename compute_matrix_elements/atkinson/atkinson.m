function Min = atkinson(k_0,k,c4n,Nb,N,W,spectrum,r_ball)
    x = c4n(Nb(:,1),:); % boundary of disk
    dx = vecnorm(x - c4n(Nb(:,2),:),2,2);
    W_bdry = W(Nb(:,1),:);
    Min = zeros(2*N+1);
    pref = (k_0-k)./((spectrum-k).*(spectrum-k_0));

    for alpha=-N:N
        ea = exp(1i*alpha.*atan2(x(:,2),x(:,1)))/sqrt(2*pi*r_ball);
        eaW = (ea.*dx)'*W_bdry;
        for beta=-N:N
            eb = exp(1i*beta.*atan2(x(:,2),x(:,1)))/sqrt(2*pi*r_ball);
            ebW = W_bdry'*(eb.*dx);
            Min(alpha+N+1,beta+N+1) = sum((pref.*eaW.').*ebW);
        end
    end
end