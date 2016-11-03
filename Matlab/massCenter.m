% Центр масс набора точек

function [C R] = massCenter(X)
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
end