function t2 = hotellingWilliams(jk, jh, kh, n)

% compute determinant of 3x3 matrix being tested
R = (1 - jk^2 - jh^2 - kh^2) + (2*jk*jh*kh);

% definition of r
r = 0.5*(jk+jh);


t2 = (jk - jh) * sqrt(((n-1)*(1+kh))/(2*((n-1)/(n-3))*abs(R) + (r^2)*((1-kh)^3)));

end