function [d] = safaee(a, b, X0, Y0)

t = 0:0.0001:2*pi;
Xe = a*cos(t);
Ye = b*sin(t);

m = 1 - a^2/b^2;
a5 = b^2 * m^2;
a4 = - (a^6)*X0 /(b^2*a5);
a3 = - 2*a^4*m*X0/a5;
a2 = (a^2*Y0^2 - a^2*b^2*m^2 + a^4*X0^2/b^2)/a5;
a1 = 2*a^2*m*X0/a5;

[x1, x2, x3, x4] = solve4(1, a1, a2,a3, a4);
x = [x1, x2, x3, x4];
y = zeros(4, 1);
d = zeros(4, 1);

for i = 1:1:4
    y(i) = (b/a)*sqrt(a^2 - x(i)^2);
    d(i) = sqrt((X0-x(i))^2 + (Y0-y(i))^2);
    
    if isreal(x(i))
        x(i)^2/a^2 + y(i)^2/b^2 - 1
        Xn = x(i):0.0001:X0;
        kn = a^2*Y0/b^2/X0;
        Yn = kn * (Xn - x(i)) +y(i);
        
        X = x(i):0.0001:X0;
        k = (Y0 - y(i)) / (X0 - x(i));
        Y = k * (X - x(i)) + y(i);
        
        plot(Xe, Ye, Xn, Yn, X, Y);
    end
end

end