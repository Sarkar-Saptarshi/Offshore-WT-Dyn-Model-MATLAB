function  [freq,Amp] = FS(t,x)
% Inputs: 
% t = time vector
% x = response vector
% Outputs:
% freq = frequency vector in Hertz
% Amp  = Fourier amplitude of response
%  ====== X == X =====
if length(t) ~= length(x)
    error('The vectors must have same length');
end

%% Zero padding
Fs = 1/(t(2)-t(1));
lpad = 2^nextpow2(length(x));
xdft = fft(x,lpad);
xdft = xdft(1:lpad/2+1);
xdft = xdft/length(x);
xdft(2:end-1) = 2*xdft(2:end-1);
freq = 0:Fs/lpad:Fs/2;
Amp = abs(xdft);

