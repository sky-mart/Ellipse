% Решаем уравнение четвертой степени вида x^4 + a*x^3 + b*x^2 + c*x + d = 0

function [x1, x2, x3, x4] = solve4(a, b, c, d, e)

b = b/a; c = c/a; d = d/a; e = e/a;
a = b; b = c; c = d; d = e;

% Приводим к виду без третьей степени
% y^4 + p*y^2 + q*y + r = 0
p = b - 3*a^2/8;
q = a^3/8 - a*b/2 + c;
r = - 3*a^4/256 + a^2*b/16 - c*a/4 + d;

% Получаем кубическую резольвенту
% A*s^3 + B*s^2 + C*s + D = 0
A = 2;
B = -p;
C = -2*r;
D = r*p - q^2/4;

[s1, s2, s3] = solve3(A, B, C, D);

%prec = A*s1^3 + B*s1^2 + C*s1 + D
s = 0;
if s1 > p/2 
    s = s1;
elseif s2 > p/2
    s = s2;
elseif s3 > p/2
    s = s3;
end

c1 = [1, -sqrt(2*s-p), q/(2*sqrt(2*s-p)) + s];
r1 = roots(c1);

c2 = [1, sqrt(2*s-p), -q/(2*sqrt(2*s-p)) + s];
r2 = roots(c2);

x1 = r1(1) - a/4;
x2 = r1(2) - a/4;
x3 = r2(1) - a/4;
x4 = r2(2) - a/4;

% disp('Четвертая степень');
% e1 = x1^4 + a*x1^3 + b*x1^2 + c*x1 + d;
% e2 = x2^4 + a*x2^3 + b*x2^2 + c*x2 + d;
% e3 = x3^4 + a*x3^3 + b*x3^2 + c*x3 + d;
% e4 = x4^4 + a*x4^3 + b*x4^2 + c*x4 + d;

end





