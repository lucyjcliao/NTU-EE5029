function bit_error = ErrorCalculate(u,u_decode)
% [Error Calculator]
% Argument:
%   - u: original signal
%   - u_decode: decode signal
% Return:
%   - bit_error: sum of errors
    bit_error = sum(xor(u,u_decode));
end