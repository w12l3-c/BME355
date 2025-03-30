classdef hammersteinModel

    %% Setting up constants for the Hammerstein model
    properties (Constant)
        % Input vector constants
        c1_flex = 4.14;
        c2_flex = 2655.88;
        c1_ext  = 4.75;
        c2_ext  = 913.2;
       
        % State matrix
        Phi = [0.82 0.008 0 0;
               0 0.82 0 0;
               0 0 0.78 0.008;
               0 0 0 0.78];
        
        % Input matrix
        Gamma = [0 0;
                 0 0.009;
                 0 0;
                 0 0.099];
        
        % Output matrix
        C = [5436.56 0 -6795.7 0];
    end

    properties
        PW {mustBeNumeric} % pulsewidth
        xk_bar (4, 1) {mustBeNumeric}
    end

    methods

        function obj = hammersteinModel(PW, xk_bar)
            obj.PW = PW;
            obj.xk_bar = xk_bar;
        end

        % input vector
        % approximation of isometric muscle recruitment curve determines peak value of muscle response for a given stimulus
        function u_bar = u_bar(obj)
            u_bar = [obj.c1_flex * abs(tanh(obj.c2_flex * obj.PW / 2));
                 obj.c1_ext  * abs(tanh(obj.c2_ext  * obj.PW / 2))];
        end

        % state-space model of muscle activation dynamics
        function xk1_bar = xk1_bar(obj) % x k+1 bar
            u_bar = [obj.c1_flex * abs(tanh(obj.c2_flex * obj.PW / 2));
                     obj.c1_ext  * abs(tanh(obj.c2_ext  * obj.PW / 2))];
            xk1_bar = (obj.Phi * obj.xk_bar) + (obj.Gamma * u_bar);
        end

        function yk = yk(obj)
            yk = obj.C * obj.xk_bar;
        end

    end

end
