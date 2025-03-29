function acc = NRMSE(rk, yk)
    y_max = max(yk);
    y_min = min(yk);

    N = length(yk);

    acc = sqrt(1/N * sum((rk-yk)^2)/(y_max - y_min)) * 100;
end