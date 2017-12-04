function f = plot_freqdomain(Rx_Signal,fs,TITLE,colour,UpsamplingFactor,POWER_PLOT_1MHZ,PA_Power_Measured)

if POWER_PLOT_1MHZ
    
    L = length(Rx_Signal);
    NFFT = 2^nextpow2(L); % Next power of 2 from length of Rx_Signal
    Power_Rx_Signal = 10*log10(mean(abs(Rx_Signal).^2));
    PwrScaleFactor = 10^((PA_Power_Measured - Power_Rx_Signal)*0.1);
    Rx_Signal = Rx_Signal*sqrt(PwrScaleFactor);
    Rx_Signal_FD = fftshift(fft(Rx_Signal,NFFT)/L);
    Rx_Power_FD = (abs(Rx_Signal_FD).^2).';
    SubcarrierSpacing = fs/NFFT;
    % Bins_1MHz = round(UpsamplingFactor*1000/15);   
    Bins_1MHz = round(1e6/SubcarrierSpacing); 
    Rx_1MHz_Power = Bins_1MHz*smooth(Rx_Power_FD,Bins_1MHz,'moving');
    
    f = (fs/2)*linspace(-1,1,NFFT);
    % Plot double-sided amplitude spectrum.
    plot(f/1e6,10*log10(Rx_1MHz_Power),colour) 
    title(TITLE)
    xlabel('Frequency (MHz)')
    ylabel('Power (dBm/MHz)')
    axis([-fs/(2*1e6) fs/(2*1e6) -Inf Inf])
    
else
    
    Nfft    = 2048;
    Window  = kaiser(2000,9); 
    Signal_PSD = 10*log10(fftshift(pwelch(Rx_Signal,Window)));
    plot((-1:2/Nfft:1-2/Nfft)*(fs/(2e6)),Signal_PSD,colour,'LineWidth',2); 
    xlabel('Frequency (MHz)')
    ylabel('PSD')
    axis([-fs/(2*1e6) fs/(2*1e6) -Inf Inf])
    
end