% �������������� ��������� �� ������� ���������, ��������� � ��������, � ���������� ������� ���������

function Rglob = invConversion(Rel, R, C)
    Rglob = inv(R)*Rel + C;
end