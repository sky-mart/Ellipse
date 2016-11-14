points_num = 30;
Xc = -1.0;
Yc = -1.0;
a = 1.0;
b = 6.0;
alpha = 0.8;
noise_level = 0;
prec = 1e-3; 

X = zeros(2, points_num);
step = 2 * pi / (points_num - 1);

t = 0;
for i = 1:1:points_num
    x = a*cos(t);
    y = b*sin(t);
    X(1,i) = (cos(alpha)*x - sin(alpha)*y + Xc);
    X(2,i) = (sin(alpha)*x + cos(alpha)*y + Yc);
    t = t + step;
end

[rest_Xc, rest_a, rest_b, rest_alpha] = ellipseFitting(X);