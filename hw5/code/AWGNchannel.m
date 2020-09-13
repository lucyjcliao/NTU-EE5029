function r = AWGNchannel(s,snr)
% r = s+wgn(1,length(s),1/snr);
noise_var= 1/(2*(1/2)*10^(snr/10));
N = sqrt(noise_var)*(randn(1,length(s)));                                   % Gaussian noise
r = N + double(s);
end