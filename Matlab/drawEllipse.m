function drawEllipse(Xc, Yc, a, b, alpha)
    
    C1 = cos(alpha(1));
    S1 = sin(alpha(1));
    
    C2 = cos(alpha(2));
    S2 = sin(alpha(2));

    t = 0:pi/100:2*pi;
    
    x1 = a(1)*cos(t);
    y1 = b(1)*sin(t);
    x2 = a(2)*cos(t);
    y2 = b(2)*sin(t);

    X1 = C1*x1 - S1*y1 + Xc(1);
    Y1 = S1*x1 + C1*y1 + Yc(1);
    
    X2 = C2*x2 - S2*y2 + Xc(2);
    Y2 = S2*x2 + C2*y2 + Yc(2);
    
    plot(X1, Y1, X2, Y2);
end