% Расстояние от точки до эллипса

function [d] = distToEllipse(a, b, X0, Y0)

if X0 == 0 && Y0 == 0
    d = min(a, b);
else
%     m = max(a, b);
%     a = a / m;
%     b = b / m;
%     X0 = X0 / m;
%     Y0 = Y0 / m;
    
    % Решаем уравнение a0*l^4 + a1*l^3 + a2*l^2 + a3*l + a4 = 0
%     a0 = - 1 / (a^2 * b^2);
%     a1 = - 2/a^2 - 2/b^2;
%     a2 = X0^2/b^2 + Y0^2/a^2 - b^2/a^2 - a^2/b^2 - 4;
%     a3 = 2*X0^2 + 2*Y0^2 - 2*b^2 - 2*a^2;
%     a4 = X0^2*b^2 + Y0^2*a^2 - a^2*b^2;
    a0 = 1;
    a1 = 2 * (a^2 + b^2);
    a2 = a^4 + b^4 + 4*a^2*b^2 - X0^2*a^2 - Y0^2*b^2;
    a3 = 2*a^2*b^2 * (a^2 + b^2 - X0^2 - Y0^2);
    a4 = a^2*b^2 * (a^2*b^2 - X0^2*b^2 - Y0^2*a^2);

    c = [a0, a1, a2, a3, a4];
    [l1, l2, l3, l4] = solve4(a0, a1, a2, a3, a4);
    l = [l1, l2, l3, l4];
     disp(l);
    r = roots(c)

    minX = - realmax;
    minY = - realmax;
    minVal = (minX - X0)^2 + (minY - Y0)^2;

    for k = 1:1:4
        root = l(k);
        if isreal(root) || imag(root) < 1e-20
            root = complex(real(root), 0);
            if ((2 + 2*root/a^2 > 0) && (2 + 2*root/b^2 >= 0)) || ((2 + 2*root/a^2 >= 0) && (2 + 2*root/b^2 > 0))
                if root == -(a^2)
                    tmpY = Y0 / (1 + root/b^2);
                    tmpX = sign(tmpY) * (b/a) * sqrt(a^2 - tmpY^2);
                elseif root == -(b^2)
                    tmpX = X0 / (1 + root/a^2);
                    tmpY = sign(tmpX) * (a/b) * sqrt(b^2 - tmpX^2);
                else
                    tmpX = X0 / (1 + root/a^2);
                    tmpY = Y0 / (1 + root/b^2);
                end
    %             tmpX^2/a^2 + tmpY^2/b^2
                if tmpX * X0 >= 0 && tmpY * Y0 >= 0
                    val = (tmpX - X0)^2 + (tmpY - Y0)^2;
                    if val < minVal
                        minX = tmpX;
                        minY = tmpY;
                        minVal = val;
                    end
                end
            end
        end
    end
%     minX = minX * m;
%     minY = minY * m;
%     a = a * m;
%     b = b * m;
%     X0 = X0 * m;
%     Y0 = Y0 * m;
    d = sqrt((minX-X0)^2 + (minY-Y0)^2);
%end
%end
%     d = [minX; minY];
%     Xc = minX
%     Yc = minY
    
%     X = Xc:0.0001:X0;
%     k = (Y0 - Yc) / (X0 - Xc);
%     Y = k * (X - Xc) + Yc;
%     %[X, Y] = lineCurve([Xc, Yc], [X0, Y0]);
%     
%     Xn = X;
%     kn = a^2*Yc/b^2/Xc;
%     Yn = kn * (Xn - Xc) + Yc;
%     
%     dk = k - kn
%     
%     Xr = 0:0.0001:Xc;
%     Yr = (Yc/Xc) * Xr;
% %     [Xr, Yr] = lineCurve([0, 0], [Xc, Yc]);
%     
%     t = 0:0.0001:2*pi;
%     Xe = a*cos(t);
%     Ye = b*sin(t);
%     
%     plot(Xe, Ye, X, Y, Xr, Yr, Xn, Yn);
% 
%     Xt = minX - 0.001;
%     Yt = (b/a)*sqrt(a^2 - Xt^2);
%     
%     val = (Xt - X0)^2 + (Yt - Y0)^2;
%     if val < minVal
%         val
%         minVal
%     end
%     
%     Xt = minX + 0.001;
%     Yt = (b/a)*sqrt(a^2 - Xt^2);
%     
%     val = (Xt - X0)^2 + (Yt - Y0)^2;
%     if val < minVal
%         val
%         minVal
    end
end
