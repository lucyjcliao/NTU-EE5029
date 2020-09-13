function s = Modulation(c)
% [Moldulation]
% Useing BPSK: 0 -> 1(+sqrt(Eb)), 1 -> -1(-sqrt(Eb))
% Argument:
%   c: encode codeword
% Returen:
%   s: modulation results
c(c==1) = -1;
c(c==0) = 1;
s = c;
end