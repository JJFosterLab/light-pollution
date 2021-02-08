%for troubleshooting only
function y = sub_feval(fun, x)
    y = polyval([fun.p1 fun.p2 fun.p3 fun.p4], x);
end