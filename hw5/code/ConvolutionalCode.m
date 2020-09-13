function [y,new_state] = ConvolutionalCode(input,state)
% [Convolutional code encoder]
% Argument:
%   - input: input bit.
%   - state: the state in array form
% Return:
%   - y: encoder codewords
%   - new_state: new state in array form
G = [1 0 0 1 1; 1 1 1 0 1];                                                 % generator matrix
y_1 = mod(sum(and([input state],G(1,:))),2);                                % codeword number 1
y_2 = mod(sum(and([input state],G(2,:))),2);                                % codeword number 2
y = [y_1 y_2];
new_state = [input state(1:end-1)];
end
