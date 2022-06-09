function val = e_alpha(alpha,r_ball,x)
    val = exp(1i*alpha.*atan2(x(2),x(1)))/sqrt(2*pi*r_ball);
end