% Тестирование алгоритма аппроксимации эллипса
pi=3.14159265358979323846;

Xc = -7;
Yc = 21;
a = 130;
b = 90;
alpha = pi/8;
C = cos(alpha);
S = sin(alpha);

step = pi / 1239;
t = 0:pi/1239:2*pi;

x = a*cos(t);
y = b*sin(t);

%rnd = rand();
X = C*x - S*y + Xc;
%rnd = rand();
Y = S*x + C*y + Yc;
R = [X; Y];
%N = wgn(2, size(R, 2), -5)*2 + (randn(size(R))-0.5)*4;
%R = R + N;

[Xcans aAns bAns alAns] = ellipseFitting(R)







