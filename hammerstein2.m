classdef hammerstein2
    %% Setting up constants for the Hammerstein model
    properties (Constant)
        % Recruitment curve constants
        c1_flex = 4.14;
        c2_flex = 2655.88;
        c1_ext  = 4.75;
        c2_ext  = 913.2;

        % State matrix
        Phi = [0.82  0.008  0     0;
               0     0.82  0     0;
               0     0     0.78  0.008;
               0     0     0     0.78];

        % Input matrix
        Gamma = [0    0;
                 0.009    0;
                 0    0;
                 0    0.009];

        % Output matrix
        C = [5436.56  0  -6795.7  0];
    end

    properties
        % Separate pulse widths for flexor and extensor
        PW_f (1,1) {mustBeNumeric} = 0;   % Flexor pulse width
        PW_e (1,1) {mustBeNumeric} = 0;   % Extensor pulse width

        % State vector
        xk_bar (4,1) {mustBeNumeric} = zeros(4,1);
    end

    methods
        %% Constructor
        function obj = hammerstein2(PW_f, PW_e, xk_bar)
            if nargin > 0
                obj.PW_f = PW_f;
            end
            if nargin > 1
                obj.PW_e = PW_e;
            end
            if nargin > 2
                obj.xk_bar = xk_bar;
            end
        end

        %% Recruitment curve -> input vector
        function u_bar = u_bar(obj)
            u_bar = [
                obj.c1_flex * abs(tanh(obj.c2_flex * obj.PW_f / 2));
                obj.c1_ext  * abs(tanh(obj.c2_ext  * obj.PW_e / 2))
            ];
        end

        %% State update
        function xk1_bar = xk1_bar(obj)
            u = obj.u_bar();  % Get the input vector
            xk1_bar = obj.Phi * obj.xk_bar + obj.Gamma * u;
        end

        %% Output (grip force)
        function yk = yk(obj)
            yk = obj.C * obj.xk_bar;
        end
    end
end
