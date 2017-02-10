% Аппроксимация набора точек эллипсом

function [Xc, a, b, alpha] = ellipseFitting(X)

    m = size(X, 2);
% Начальные параметры берем из аппроксимации кругом
    [Cs Rs] = massCenter(X);%circleFitting(X);%
    Xc = Cs;
    a = Rs;
    b = Rs;
    alpha = 0;
%     Xc = zeros(2, 1);
%     Xc(1) = -7.000022;
%     Xc(2) = 21.000046;
%     a = 129.819432;
%     b = 90.102916;
%     alpha = 0.394887;
    
    
    J = zeros(2*m, 5);
    e = zeros(2*m, 1);
    
    iter = 0;
    while 1
        fprintf('Итерация %d\n', iter + 1);
        for i = 1:1:m
            C = cos(alpha);
            S = sin(alpha);

            R = [C S
                -S C];
            x = conversion(X(:,i), R, Xc);

            xs = distToEllipse(a, b, x(1), x(2));
            Xs = invConversion(xs, R, Xc);

            Q = [b^2 * xs(1)                  a^2 * xs(2)
                (a^2-b^2)*xs(2) + b^2*x(2)   (a^2-b^2)*xs(1) - a^2*x(1)];

            M = B(a, b, C, S, xs, x);

            Ji = inv(R) * inv(Q) * M; 
            ei = X(:,i) - Xs;

            for j = 1:1:2
                J(2*(i-1)+j,:) = Ji(j,:);
                e(2*(i-1)+j,:) = ei(j,:);
            end
        end

        if abs(a - b) < 1e-3
            Js = J(1:2*m,1:4);
            da = linsolve(Js,e);
            for j = 1:1:4
                fprintf('da(%d) = %f\n', j, da(j));
            end
        else
            da = linsolve(J,e);
            alpha = alpha + da(5); 
            for j = 1:1:5
                fprintf('da(%d) = %f\n', j, da(j));
            end
        end
       
   
        Xc(1) = Xc(1) + da(1);
        Xc(2) = Xc(2) + da(2);
        a = a + da(3);
        b = b + da(4);
        
        
        
        ea = da(3) / a;
        eb = da(4) / b;
        s = sqrt(da(1)^2 + da(2)^2);
        
        iter = iter + 1;
        % Условие выхода из цикла
        if ea < 1e-3 && eb < 1e-3 && 2*s/(a+b) < 1e-5 %da'*da < 1e-20
            break
        end
    end
    
    fprintf('Эллипс: %d итераций\n', iter);  
    
%     t = 0:0.0001:2*pi;
%     x = a*cos(t);
%     y = b*sin(t);
%     r = [x; y];
%     el = zeros(size(r));
%     for i = 1:1:size(r, 2)
%         el(:,i) = inv(R)*r(:,i) + Xc;
%     end
%     
%     plot(el(1,:), el(2,:), X(1,:), X(2,:), 'o');
end