% Ёллипс дл€ графика

function [x, y] = ellipseCurve(a, b)

t = 0:0.0001:2*pi;
x = a*cos(t);
y = b*sin(t);

end