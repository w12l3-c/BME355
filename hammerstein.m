classdef hammerstein
    %% setting up constants for hammerstein model
    %% Setting up constants for the Hammerstein model
    properties (Constant)
        % Input vector constants
        c1_flex = 4.14
        c2_flex = 2655.88
        c1_ext  = 4.75
        c2_ext  = 913.2
       
        % State matrix
        Phi = [0.82 0.008 0 0;
               0 0.82 0 0;
               0 0 0.78 0.008;
               0 0 0 0.78]
        
        % Input matrix
        Gamma = [0 0;
                 0 0.009;
                 0 0;
                 0 0.099]
        
        % Output matrix
        C = [5436.56 0 -6795.7 0]
    end
    properties
        PW {mustBeNumeric}
    end

%% approximation of isometric muscle recruitment curve determines peak value of muscle response for a given stimulus
function u = IRC(c1, c2, PW)
    u = c1 * abs(tanh(c2 * PW / 2));
end

end

