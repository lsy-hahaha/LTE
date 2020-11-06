function [ sc_snr ] = calculate_SC_SNR( h,snr )

sigma=10^(-snr/10);

H=fft(h,2048);
effeH=H(((2048-1200)/2+1):((2048-1200)/2+1200));
size(effeH);
sc_snr=abs(effeH).^2/sigma;

end






















