% Аппроксимация набора точек кругом

function [C R] = circleFitting(X)

% Начальные параметры берем из центра масс
    C = zeros(2, 1);
    m = size(X, 2);
    for i = 1:1:m
        C = C + X(:,i);       
    end
    C = C / m;
    
    R = 0;
    for i = 1:1:m
        R = R + distance(C, X(:,i))^2;
    end
    R = sqrt(R/m);
    
    dRda = [1 0 0];
    dCda = [0 1 0;
            0 0 1];
    I = eye(2);
    
    J = zeros(2*m, 3);
    e = zeros(2*m, 1);
    da = zeros(3, 1);
    
    iter = 0;
    while 1
        for i = 1:1:m
            d = distance(X(:,i), C);
            Ji = zeros(2,3);
            Ji = dCda + (X(:,i) - C)*dRda/d - R*(I - (X(:,i) - C)*(X(:,i) - C)'/d^2)*dCda/d;
            ei = zeros(2,1);
            ei = (d - R)*(X(:,i) - C)/d;
            for j = 1:1:2
                J(2*(i-1)+j,:) = Ji(j,:);
                e(2*(i-1)+j,:) = ei(j,:);
            end
        end

        da = linsolve(J,e);
        R = R + da(1);
        C(1) = C(1) + da(2);
        C(2) = C(2) + da(3);
        
        iter = iter + 1;
        if da'*da < 1e-20
            break
        end
    end
    
    fprintf('Круг: %d итераций\n', iter); 
end