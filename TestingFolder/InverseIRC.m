function [PW_f, PW_e] = InverseIRC(u_f, u_e)
    c1_flex = 4.14;
    c2_flex = 2655.88;
    c1_ext  = 4.75;
    c2_ext  = 913.2;

    PW_f = (2 / c2_flex) * atanh(sqrt(u_f^2/c1_flex^2));
    PW_e = (2 / c2_ext) * atanh(sqrt(u_e^2/c1_ext^2));
end