function [PW_f, PW_e] = InverseIRC(u_f, u_e)
    % InverseIRC calculates the required pulse widths for stimulation 
    % (PW_f for the flexor and PW_e for the extensor) from the desired 
    % normalized force inputs (u_f and u_e) based on the inverse isometric 
    % recruitment curve (IRC). This function uses a hyperbolic tangent-based 
    % model where the maximum force the muscle can generate is capped by c1.
    %
    % Inputs:
    %   u_f - desired normalized force for the flexor muscle (scalar)
    %   u_e - desired normalized force for the extensor muscle (scalar)
    %
    % Outputs:
    %   PW_f - computed pulse width for the flexor muscle (in appropriate units)
    %   PW_e - computed pulse width for the extensor muscle (in appropriate units)
    
    % Define recruitment curve parameters for the flexor muscle.
    % c1_flex represents the maximum force and c2_flex controls the slope.
    c1_flex = 4.14;            % Maximum normalized force for flexor (N)
    c2_flex = 2655.88;         % Slope parameter (scaling factor) for flexor
    
    % c1_ext represents the maximum force and c2_ext controls the slope.
    c1_ext  = 4.75;            % Maximum normalized force for extensor (N)
    c2_ext  = 913.2;           % Slope parameter (scaling factor) for extensor

    % Saturate the input force commands so they don't exceed the maximum force
    % (multiplying by 0.999 to avoid exactly hitting the saturation limit).
    u_f = min([u_f, c1_flex * 0.999]);
    u_e = min([u_e, c1_ext  * 0.999]);

    % Compute the required pulse width for the muscles.
    PW_f = (2 / c2_flex) * atanh( sqrt((u_f^2) / (c1_flex^2)) );
    PW_e = (2 / c2_ext) * atanh( sqrt((u_e^2) / (c1_ext^2)) );

    fprintf("pw_f: %d\n", PW_f);
    fprintf("pw_e: %d\n", PW_e);

    % If no force is required for a muscle, explicitly set the corresponding pulse width to zero.
    if u_f == 0
        PW_f = 0;
    end
    if u_e == 0
        PW_e = 0;
    end
end
