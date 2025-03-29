function PW = InverseIRC(u)
    c1 = 4.14;
    c2 = 2655.8;
    PW = 2 / c2 * atan2(sqrt(u^2/c1^2));
end