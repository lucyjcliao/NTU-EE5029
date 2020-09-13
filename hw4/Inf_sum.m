function y = Inf_sum(a,b,m)
for i = 1:length(a)
    if isinf(a(i)) && isinf(b(i))
        y(i) = -Inf;
    elseif isinf(a(i))
        y(i) = mod(b(i),2^m-1);
    elseif isinf(b(i))
        y(i) = mod(a(i),2^m-1);
    else
        y(i) = GF_add([a(i);b(i)],m);
    end
end
        