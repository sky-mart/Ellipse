% Расстояние между двумя точками

function d = distance(X, Y)
    d = sqrt((X(1) - Y(1))^2 + (X(2) - Y(2))^2);
end