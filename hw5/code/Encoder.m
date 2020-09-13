function c = Encoder(u)
% (2,1,4) Convolution code
memory = zeros(1,4);
v = [];
for i = 1:length(u)
    [y,new_state] = ConvolutionalCode(u(i),memory);
    v = [v y];
    memory = new_state;
    
end

c = v;

end