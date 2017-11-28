function x = assignall(len, val)
    x = zeros(len,1);
    for i=1:len
        x(i) = val;
    end
end
