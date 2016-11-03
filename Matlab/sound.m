[d, Fs] = audioread('C:\Users\user\Desktop\y2.wav');
delete('C:\Users\user\Desktop\out.wav');
y = d(1:32768*2,1);
T = 1/Fs;                     % Sample time
L = size(y, 1);                     % Length of signal
t = (0:L-1)*T;                % Time vector

%NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,L);
f = Fs/2*linspace(0,1,L);

A = abs(Y);
Re = real(Y);
Im = imag(Y);
Re2 = Re;
Im2 = Im;
MaxA = max(A);
ph = Y;
for i = 1:1:L
    ph(i) = atan(Im(i)/Re(i)) + pi/2*sign(Im(i))*(1-sign(Re(i)));
   % ph(i) = rand()*pi;
    %ph(i) = 0;
    A(i) = rand()*MaxA*1000;
    Re2(i) = A(i)*cos(ph(i));
    Im2(i) = A(i)*sin(ph(i));
end
plot(f,A,f,ph);
%plot(f, ph);

%A = rand(size(A))*(max(A)/2);
%ph = rand(size(A))*pi;
Y2 = complex(Re2, Im2);

y2 = ifft(Y2);
y2=y2/max(abs(y2));
y2Int = y2 * 32767;
audiowrite('C:\Users\user\Desktop\out.wav', 10*y2, Fs);
audiowrite('C:\Users\user\Desktop\out2.wav', 10*y, Fs);
