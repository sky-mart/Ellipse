% ѕреобразование координат из системы координат, св€занной с эллипсом, в глобальную систему координат

function Rglob = invConversion(Rel, R, C)
    Rglob = inv(R)*Rel + C;
end