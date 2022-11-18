clear all; close all; clc;

T0 = [12.5, 25];    %ps
for x = 1 : length(T0)
    L = 25;   %km
    D = 17;
    lambda = 1550;
    c = 3*10^8;
    B2 = -lambda*lambda*D/(2*pi*c*1e-3);
    T = 512;   %FFT window size (ps)
    N = 2048;   % number of sampling
    
    dt = T/N;   %time interval
    df = 1/T;   %frequency interval
    t = (-N/2:N/2-1)*dt; %time vector    
    f = (-N/2:N/2-1)*df; %frequency vector
    w = 2*pi*f;
    
    % cumpute dispersion effect
    A0T = exp(-t.*t./(2*T0(x)*T0(x)));  
    
    A0f = fftshift(ifft(A0T));
    
    Disp = exp(1i*B2*w.^2*L/2);
    
    ALf_1 = A0f .* Disp;   %propagation in frequency domain %output spectrum 1
    ALt_1 = fft(fftshift(ALf_1));   %output pulse 1 
    ALf_2 = ALf_1 .* Disp;  %propagation in frequency %output spectrum 2
    ALt_2 = fft(fftshift(ALf_2));   %output pulse 2
    
    %Plot in time domain
    figure('Position',[100 100 800 300])
    subplot(121)
    plot(t,abs(A0T).^2)
    grid on;
    hold on;
    plot(t,abs(ALt_1).^2, '--')
    hold on;
    plot(t,abs(ALt_2).^2, '-.')
    ylim([0 1.3])
    xlim([-200 200])
    xticks([-200:50:200])
    legend('0km','25km','50km')
    title("A(z,T) | T0 = "+ T0(x) )
    xlabel('T[Ps]')
    ylabel('P[mW]')
    
    %Plot PSD
    subplot(122)
    plot(f, 10*log10(abs(N*A0f).^2/N*dt.^2))
    grid on;
    hold on;
    plot(f, 10*log10(abs(N*ALf_1).^2/N*dt.^2), '--')
    hold on;
    plot(f, 10*log10(abs(N*ALf_2).^2/N*dt.^2), '-.')
    ylim([-40 20])
    xlim([-0.05 0.05])
    xticks([-0.05:0.01:0.05])
    legend('0km','25km','50km')
    title("A(z,w) | T0 = "+ T0(x) )
    xlabel('f[THz]')
    ylabel('PSD [dBm/THz]')
end