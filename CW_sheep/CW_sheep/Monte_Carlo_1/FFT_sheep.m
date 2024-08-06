rawTable = readmatrix('D:\Code\TFO\CW_sheep\07252022_sheep\extract_07-25-22_R2');
period = 60; 
X=rawTable(:,2);
L=length(X);
Fs=80;
Y = fft(X);

P2 = (1/(Fs*L))*abs(Y).^2;
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

t=linspace(0,60,L);
figure
f = Fs*(0:(L/2))/L;
plot(f,10*log10(P1)-15)
title("Ch10, SD=10cm")
xlabel("Frequency (Hz)")
ylabel("Power (dBm)")
grid on
xlim([0 6])
ylim([-120 -20])

figure
plot(t,X)