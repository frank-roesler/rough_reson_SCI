function L = build_lattice(z1, z2, h)
%   Defines lattice L in the complex plane from corner points z1, z2
%   with horizontal_resolution grid points between re(z1) and re(z2)
    r1 = real(z1);
    r2 = real(z2);
    i1 = imag(z1);
    i2 = imag(z2);
    X = r1:h:r2;
    Y = i1:h:i2;
    [XX,YY] = meshgrid(X,Y);
    L = XX+1i*YY;
end