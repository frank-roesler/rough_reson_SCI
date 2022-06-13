function Min = atkinson_corrector(k_0,k,c4n,Nb,N,W,spectrum,r_ball)
% Computes M_in(k)-M_in(k_0) via Atkinson acceleration.
    x = c4n(Nb(:,1),:); % boundary of disk
    dx = vecnorm(x - c4n(Nb(:,2),:),2,2);
    W_bdry = (W(Nb(:,1),:) + W(Nb(:,2),:))/2;
    pref = (k^2-k_0^2)./((spectrum-k^2).*(spectrum-k_0^2));
	e_bdry = exp(1i*(-N:N).*atan2(x(:,2),x(:,1)))/sqrt(2*pi*r_ball);
    eaW_big = (e_bdry.*dx)'*W_bdry;
    Min = (pref.'.*eaW_big)*eaW_big';
end