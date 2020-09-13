function y = GF_mul(x,m)
    for i = 1:size(x,2)
        if any(isinf(x(:,i)))
            y(i) = -Inf;
        else
            y(i) = mod(sum(x(:,i)),2^m-1);
        end
    end
end
        