function tf = in_Julia(z, iters, escape_radius)
% Returns "True" if z is in Julia set.
    phi = (1+sqrt(5))/2;
    c = (1-phi);
    f = z.^2+c;
    for n=1:iters
        f = f.^2+c;
%         maxf(abs(f)>maxf) = abs(f(abs(f)>maxf))
    end
    tf = abs(f) < escape_radius;
end