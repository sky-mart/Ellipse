% Прямая для графика через две точки

function [x, y] = lineCurve(M, N)

minX = min(M(1), N(1));
maxX = max(M(1), N(1));

x = minX:0.0001:maxX;
y = (

end